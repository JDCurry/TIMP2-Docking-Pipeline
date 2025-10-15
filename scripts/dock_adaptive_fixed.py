#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Modified dock_chunk_adaptive.py that preserves original molecule names from SMI files.
Key changes:
1. PDBQT files are named using the original molecule ID (e.g., ZINC001176034359_1.pdbqt)
2. Summary CSV uses original names, making it possible to map back to SMILES
3. Added safety for filesystem-incompatible characters in molecule names
"""

import sys, re, subprocess, shutil, time, os
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import psutil
from datetime import datetime

# ----------- CONFIG (edit these paths/knobs) -----------
OBABEL   = r"D:\timp2\OpenBabel-3.1.1\obabel.exe"
VINA     = r"D:\timp2\Tools\vina\vina.exe"
RECEPTOR = r"D:\timp2\receptor.pdbqt"
WORK_ROOT = r"D:\timp2\work_strict"   # each chunk gets its own subdir

# Site boxes
SITE1 = dict(cx=27.991873, cy=18.851075, cz=15.488124, sx=24, sy=24, sz=24)
SITE2 = dict(cx=23.255162, cy=27.655522, cz=14.632348, sx=18, sy=18, sz=18)

# Vina
EXHAUSTIVENESS = 6
NUM_MODES      = 1          # fast triage
VINA_THREADS   = 2

# Adaptive CPU settings
TARGET_CPU_UTIL = 0.93      # 93% target CPU usage
DEFAULT_N_PROCESSES = 32    # fallback if benchmarking fails
BENCHMARK_LIGANDS = 50      # number of ligands for benchmarking
MIN_PROCESSES = 8           # minimum process count to test
MAX_PROCESSES = min(48, os.cpu_count() * 2)  # reasonable upper bound

# Pose saving policy
KEEP_THRESH        = -7.0    # save pose if score <= this
TOP_KEEP_PER_SITE  = 2000    # cap winners per site per chunk
# -------------------------------------------------------

CRASH_EXIT_CODES = {3221225477, -1073741819}  # 0xC0000005 variations

def sanitize_filename(name: str) -> str:
    """
    Make a molecule name safe for use as a filename.
    Replaces problematic characters but preserves the essence of the name.
    """
    # Replace filesystem-unsafe characters
    safe = re.sub(r'[<>:"/\\|?*]', '_', name)
    # Remove any control characters
    safe = re.sub(r'[\x00-\x1f\x7f]', '', safe)
    # Limit length to avoid filesystem issues (keep it under 200 chars)
    if len(safe) > 200:
        safe = safe[:200]
    return safe.strip() or "unknown"

def log_performance(log_path: Path, message: str, cpu_percent=None, memory_mb=None):
    """Log performance metrics with timestamp"""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    cpu_str = f"CPU:{cpu_percent:.1f}%" if cpu_percent is not None else ""
    mem_str = f"MEM:{memory_mb:.0f}MB" if memory_mb is not None else ""
    metrics = f"[{cpu_str} {mem_str}]" if cpu_str or mem_str else ""
    
    with open(log_path, "a", encoding="utf-8") as f:
        f.write(f"{timestamp} {metrics} {message}\n")
    print(f"[{timestamp}] {metrics} {message}")

def get_system_metrics():
    """Get current CPU and memory usage"""
    cpu_percent = psutil.cpu_percent(interval=1)
    memory_info = psutil.virtual_memory()
    memory_mb = (memory_info.total - memory_info.available) / (1024 * 1024)
    return cpu_percent, memory_mb

def benchmark_process_count(ligands_subset: list, work_dir: Path, perf_log: Path) -> int:
    """
    Benchmark different process counts to find optimal CPU utilization.
    Returns the best process count for target CPU usage.
    """
    log_performance(perf_log, f"Starting benchmark with {len(ligands_subset)} ligands")
    
    # Test range: start conservative, step up
    test_counts = []
    step = max(4, MAX_PROCESSES // 8)
    for n in range(MIN_PROCESSES, MAX_PROCESSES + 1, step):
        test_counts.append(n)
    
    # Always test the default as a baseline
    if DEFAULT_N_PROCESSES not in test_counts:
        test_counts.append(DEFAULT_N_PROCESSES)
    test_counts.sort()
    
    log_performance(perf_log, f"Testing process counts: {test_counts}")
    
    results = {}
    
    for n_proc in test_counts:
        log_performance(perf_log, f"Benchmarking {n_proc} processes...")
        
        # Create temp directories for this benchmark
        bench_out_dir = work_dir / "bench_out"
        bench_log_dir = work_dir / "bench_logs"
        bench_out_dir.mkdir(exist_ok=True)
        bench_log_dir.mkdir(exist_ok=True)
        
        start_time = time.time()
        completed_jobs = 0
        cpu_readings = []
        
        try:
            with ProcessPoolExecutor(max_workers=n_proc) as ex:
                # Submit benchmark jobs
                futures = []
                for lig_path, orig_name in ligands_subset:
                    futures.append(ex.submit(dock_one, lig_path, orig_name, "site1", SITE1, 42, bench_out_dir, bench_log_dir))
                
                # Monitor progress and CPU usage
                monitor_start = time.time()
                for f in as_completed(futures):
                    try:
                        f.result()  # Just to trigger any exceptions
                        completed_jobs += 1
                        
                        # Sample CPU usage every few completions
                        if completed_jobs % max(1, len(futures) // 10) == 0:
                            cpu_pct, _ = get_system_metrics()
                            cpu_readings.append(cpu_pct)
                            
                    except Exception as e:
                        log_performance(perf_log, f"Benchmark job failed: {e}")
                
        except Exception as e:
            log_performance(perf_log, f"Benchmark failed for {n_proc} processes: {e}")
            continue
        
        elapsed_time = time.time() - start_time
        
        # Calculate metrics
        if elapsed_time > 0 and completed_jobs > 0:
            throughput = completed_jobs / (elapsed_time / 60)  # ligands per minute
            avg_cpu = sum(cpu_readings) / len(cpu_readings) if cpu_readings else 0
            
            # Score: balance throughput and CPU efficiency
            # Penalize if we're way over/under target CPU
            cpu_efficiency = 1.0 - abs(avg_cpu/100 - TARGET_CPU_UTIL) * 2
            cpu_efficiency = max(0.1, cpu_efficiency)  # don't go negative
            
            score = throughput * cpu_efficiency
            
            results[n_proc] = {
                'throughput': throughput,
                'avg_cpu': avg_cpu,
                'efficiency': cpu_efficiency,
                'score': score,
                'completed': completed_jobs
            }
            
            log_performance(perf_log, 
                f"  {n_proc} processes: {throughput:.1f} lig/min, "
                f"CPU {avg_cpu:.1f}%, efficiency {cpu_efficiency:.3f}, score {score:.2f}")
        
        # Cleanup benchmark files
        shutil.rmtree(bench_out_dir, ignore_errors=True)
        shutil.rmtree(bench_log_dir, ignore_errors=True)
        
        # Brief pause between tests
        time.sleep(2)
    
    if not results:
        log_performance(perf_log, f"All benchmarks failed! Using default {DEFAULT_N_PROCESSES}")
        return DEFAULT_N_PROCESSES
    
    # Select best performer
    best_n = max(results.keys(), key=lambda k: results[k]['score'])
    best_result = results[best_n]
    
    log_performance(perf_log, 
        f"SELECTED: {best_n} processes (score {best_result['score']:.2f}, "
        f"{best_result['throughput']:.1f} lig/min, {best_result['avg_cpu']:.1f}% CPU)")
    
    # Log all results for analysis
    with open(work_dir / "benchmark_results.csv", "w", encoding="utf-8") as f:
        f.write("processes,throughput_lig_per_min,avg_cpu_percent,efficiency,score,completed_jobs\n")
        for n_proc, data in sorted(results.items()):
            f.write(f"{n_proc},{data['throughput']:.2f},{data['avg_cpu']:.1f},"
                   f"{data['efficiency']:.3f},{data['score']:.2f},{data['completed']}\n")
    
    return best_n

def read_smi_lines(smi_path: Path):
    """Yield (idx, smiles, title) 1-based. Title (column 2+) optional."""
    with open(smi_path, "r", encoding="utf-8", errors="ignore") as fh:
        for i, line in enumerate(fh, 1):
            line = line.strip()
            if not line: 
                continue
            parts = line.split()
            smi = parts[0]
            title = parts[1] if len(parts) > 1 else f"mol_{i:08d}"
            yield i, smi, title

def obabel_convert_single(obabel, smiles: str, out_path: Path, title: str, log: Path, timeout=600):
    """
    Convert a single SMILES to a specific PDBQT path deterministically.
    We pass the SMILES via -: to avoid file slicing and we set the title.
    """
    # Quick sanity parse first
    cmd = [obabel, "-:'{}'".format(smiles), "-osmi"]
    try:
        run_logged(cmd, log, timeout=60)  # will fail fast on parse errors
    except RuntimeError:
        return False  # bad SMILES, skip

    tmp_out = out_path.with_suffix(".tmp.pdbqt")
    cmd = [obabel, "-:'{}'".format(smiles),
           "-opdbqt", "-O", str(tmp_out),
           "--gen3d", "-h", "--partialcharge", "gasteiger", "--title", title]
    try:
        run_logged(cmd, log, timeout=timeout)
        tmp_out.replace(out_path)
        return True
    except RuntimeError:
        # native crash or 3D gen failure
        if tmp_out.exists():
            tmp_out.unlink(missing_ok=True)
        return False

def convert_smi_to_pdbqt_preserving_names(obabel, smi_file: Path, out_dir: Path, log: Path):
    """
    Convert SMILES to PDBQT files, preserving the original molecule names.
    Returns a mapping of safe_filename -> original_name for later use.
    """
    bad_log = out_dir.parent / "logs" / "bad_smiles.txt"
    name_mapping_file = out_dir.parent / "logs" / "name_mapping.tsv"
    
    # Build an index of lines
    indexed = list(read_smi_lines(smi_file))
    if not indexed:
        return {}, 0, 0
    
    name_mapping = {}
    total_ok = 0
    total_bad = 0
    
    # Write name mapping header
    with open(name_mapping_file, "w", encoding="utf-8") as f:
        f.write("safe_filename\toriginal_name\tsmiles\n")
    
    for idx, smi, title in indexed:
        # Create a safe filename from the original title
        safe_name = sanitize_filename(title)
        
        # Handle potential duplicates by appending index if needed
        if safe_name in name_mapping:
            safe_name = f"{safe_name}_{idx}"
        
        out_path = out_dir / f"{safe_name}.pdbqt"
        
        # Try to convert
        ok = obabel_convert_single(obabel, smi, out_path, title, log)
        
        if ok:
            name_mapping[safe_name] = title
            total_ok += 1
            
            # Record the mapping
            with open(name_mapping_file, "a", encoding="utf-8") as f:
                f.write(f"{safe_name}\t{title}\t{smi}\n")
        else:
            with open(bad_log, "a", encoding="utf-8") as bf:
                bf.write(f"{idx}\t{title}\t{smi}\tFAILED_conversion\n")
            total_bad += 1
    
    print(f"[{smi_file.stem}] PDBQT ok: {total_ok:,} | bad: {total_bad:,}")
    print(f"[{smi_file.stem}] Name mapping saved to: {name_mapping_file}")
    
    return name_mapping, total_ok, total_bad

def run_logged(cmd, log_path: Path, timeout=None):
    """Run a command streaming output to a logfile (avoid pipe deadlocks)."""
    with open(log_path, "a", encoding="utf-8") as lf:
        lf.write("\n>> " + " ".join([f'"{c}"' if " " in str(c) else str(c) for c in cmd]) + "\n")
        proc = subprocess.Popen(cmd, stdout=lf, stderr=lf)
        proc.wait(timeout=timeout)
        if proc.returncode != 0:
            raise RuntimeError(f"Command failed (exit {proc.returncode}): {cmd}")

def vina_score_from_pdbqt(p: Path):
    pat = re.compile(r"REMARK VINA RESULT:\s+(-?\d+(?:\.\d+)?)")
    try:
        with open(p, "r", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                m = pat.search(line)
                if m: return float(m.group(1))
    except FileNotFoundError:
        return None
    return None

def dock_one(lig_path: Path, orig_name: str, site_tag: str, box: dict, seed: int,
             out_dir: Path, log_dir: Path):
    """
    Dock one ligand at one site.
    orig_name: the original molecule name (for reporting)
    lig_path: the PDBQT file path (with safe filename)
    """
    safe_name = lig_path.stem
    out_p = out_dir / f"{site_tag}_{safe_name}.pdbqt"
    log_p = log_dir / f"{site_tag}_{safe_name}.log"
    
    if out_p.exists():
        return (orig_name, site_tag, vina_score_from_pdbqt(out_p), out_p)

    cmd = [
        VINA, "--receptor", RECEPTOR, "--ligand", str(lig_path),
        "--center_x", str(box["cx"]), "--center_y", str(box["cy"]), "--center_z", str(box["cz"]),
        "--size_x", str(box["sx"]), "--size_y", str(box["sy"]), "--size_z", str(box["sz"]),
        "--exhaustiveness", str(EXHAUSTIVENESS),
        "--num_modes",      str(NUM_MODES),
        "--seed",           str(seed),
        "--cpu",            str(VINA_THREADS),
        "--out",            str(out_p),
        "--log",            str(log_p),
    ]
    # let vina write its own log; we don't capture
    proc = subprocess.run(cmd)
    if proc.returncode != 0:
        return (orig_name, site_tag, None, out_p)
    return (orig_name, site_tag, vina_score_from_pdbqt(out_p), out_p)

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("chunk_path", help="Path to strict_XXXX.smi chunk file")
    parser.add_argument("--no-clean", action="store_true", help="Do not delete intermediate files (poses, lig_pdbqt)")
    parser.add_argument("--use-existing-only", action="store_true", help="Forbid OBabel conversion; require all PDBQTs to exist")
    parser.add_argument("--processes", type=int, help="Manual override for process count (skips benchmarking)")
    parser.add_argument("--skip-benchmark", action="store_true", help="Skip benchmarking, use default process count")
    args = parser.parse_args()

    chunk_path = Path(args.chunk_path)
    if not chunk_path.exists():
        print(f"[FATAL] Chunk not found: {chunk_path}")
        sys.exit(2)

    chunk_id = chunk_path.stem
    work_dir = Path(WORK_ROOT) / chunk_id
    lig_dir  = work_dir / "lig_pdbqt"
    out_dir  = work_dir / "vina_out"
    keep_dir = work_dir / "poses_keep"
    log_dir  = work_dir / "logs"
    for d in (lig_dir, out_dir, keep_dir, log_dir):
        d.mkdir(parents=True, exist_ok=True)
    
    obabel_log = log_dir / "obabel_convert.log"
    perf_log = log_dir / "performance.log"

    # Initialize performance logging
    cpu_pct, mem_mb = get_system_metrics()
    log_performance(perf_log, f"Starting chunk {chunk_id}", cpu_pct, mem_mb)

    # Sanity: tool paths
    for p in (Path(OBABEL), Path(VINA), Path(RECEPTOR)):
        if not p.exists():
            print(f"[WARN] Missing path: {p}")

    # 1) Convert SMILES -> ligand PDBQTs (resume-safe) with name preservation
    name_mapping = {}
    name_mapping_file = log_dir / "name_mapping.tsv"
    
    if not any(lig_dir.glob("*.pdbqt")):
        if args.use_existing_only:
            print(f"[{chunk_id}] --use-existing-only set, but no PDBQTs found. Aborting.")
            sys.exit(10)
        print(f"[{chunk_id}] Converting SMILES to PDBQT (preserving original names)...")
        log_performance(perf_log, "Starting SMILES conversion")
        name_mapping, total_ok, total_bad = convert_smi_to_pdbqt_preserving_names(
            OBABEL, chunk_path, lig_dir, obabel_log
        )
        cpu_pct, mem_mb = get_system_metrics()
        log_performance(perf_log, f"SMILES conversion completed: {total_ok} ok, {total_bad} failed", cpu_pct, mem_mb)
    else:
        print(f"[{chunk_id}] Found existing PDBQTs; loading name mapping...")
        log_performance(perf_log, "Loading existing name mapping")
        
        # Try to load existing name mapping
        if name_mapping_file.exists():
            with open(name_mapping_file, "r", encoding="utf-8") as f:
                next(f)  # skip header
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) >= 2:
                        name_mapping[parts[0]] = parts[1]
            print(f"[{chunk_id}] Loaded {len(name_mapping)} name mappings")
        else:
            print(f"[{chunk_id}] WARNING: No name mapping file found. Using filename as molecule name.")
            # Fallback: use the PDBQT filename as the molecule name
            for p in lig_dir.glob("*.pdbqt"):
                name_mapping[p.stem] = p.stem
    
    # Get list of ligands with their original names
    ligs_with_names = []
    for pdbqt_path in sorted(lig_dir.glob("*.pdbqt")):
        safe_name = pdbqt_path.stem
        orig_name = name_mapping.get(safe_name, safe_name)  # fallback to safe_name if not in mapping
        ligs_with_names.append((pdbqt_path, orig_name))
    
    print(f"[{chunk_id}] Ligands: {len(ligs_with_names):,}")
    log_performance(perf_log, f"Total ligands available: {len(ligs_with_names):,}")
    
    if not ligs_with_names:
        print(f"[{chunk_id}] No ligands found. Check OBABEL path / permissions.")
        sys.exit(3)

    # 2) Determine optimal process count
    if args.processes:
        n_processes = args.processes
        log_performance(perf_log, f"Using manual process count: {n_processes}")
    elif args.skip_benchmark:
        n_processes = DEFAULT_N_PROCESSES
        log_performance(perf_log, f"Skipping benchmark, using default: {n_processes}")
    else:
        # Select subset for benchmarking (random sample to avoid bias)
        import random
        benchmark_subset = random.sample(ligs_with_names, min(BENCHMARK_LIGANDS, len(ligs_with_names)))
        log_performance(perf_log, f"Selected {len(benchmark_subset)} ligands for benchmarking")
        
        n_processes = benchmark_process_count(benchmark_subset, work_dir, perf_log)

    # 3) Dock in parallel (both sites) with optimized process count
    log_performance(perf_log, f"Starting main docking with {n_processes} processes")
    keep_counts = {"site1": 0, "site2": 0}
    results = []
    completed_jobs = 0
    start_time = time.time()
    
    # Progress monitoring setup
    total_jobs = len(ligs_with_names) * 2  # both sites
    last_progress_time = start_time
    progress_interval = 300  # log progress every 5 minutes
    
    with ProcessPoolExecutor(max_workers=n_processes) as ex:
        futs = []
        for lig_path, orig_name in ligs_with_names:
            futs.append(ex.submit(dock_one, lig_path, orig_name, "site1", SITE1, 42, out_dir, log_dir))
            futs.append(ex.submit(dock_one, lig_path, orig_name, "site2", SITE2, 43, out_dir, log_dir))
        
        for f in as_completed(futs):
            try:
                orig_name, site, score, out_p = f.result()
                results.append((orig_name, site, score))
                completed_jobs += 1
                
                # Keep only winners
                if score is not None and score <= KEEP_THRESH and keep_counts[site] < TOP_KEEP_PER_SITE:
                    dest = keep_dir / out_p.name
                    if not dest.exists() and out_p.exists():
                        shutil.copy2(out_p, dest)
                        keep_counts[site] += 1
                
                # Clean up intermediate files
                if not args.no_clean and out_p.exists():
                    out_p.unlink(missing_ok=True)
                
                # Progress logging
                current_time = time.time()
                if current_time - last_progress_time >= progress_interval:
                    elapsed = (current_time - start_time) / 60
                    progress_pct = (completed_jobs / total_jobs) * 100
                    rate = completed_jobs / elapsed if elapsed > 0 else 0
                    eta_min = (total_jobs - completed_jobs) / rate if rate > 0 else 0
                    
                    cpu_pct, mem_mb = get_system_metrics()
                    log_performance(perf_log, 
                        f"Progress: {completed_jobs}/{total_jobs} ({progress_pct:.1f}%) "
                        f"Rate: {rate:.1f} jobs/min ETA: {eta_min:.0f}min", cpu_pct, mem_mb)
                    last_progress_time = current_time
                    
            except Exception as e:
                print("[ERROR]", e, file=sys.stderr)
                log_performance(perf_log, f"Docking job failed: {e}")

    # Final performance summary
    total_time = (time.time() - start_time) / 60
    final_rate = completed_jobs / total_time if total_time > 0 else 0
    cpu_pct, mem_mb = get_system_metrics()
    log_performance(perf_log, 
        f"Docking completed: {completed_jobs} jobs in {total_time:.1f}min "
        f"(avg {final_rate:.1f} jobs/min)", cpu_pct, mem_mb)

    # 4) Write per-chunk CSV with ORIGINAL molecule names
    csv_p = work_dir / f"{chunk_id}_summary.csv"
    with open(csv_p, "w", encoding="utf-8") as fh:
        fh.write("ligand,site,score_kcal_per_mol\n")
        for orig_name, site, score in sorted(results):
            fh.write(f"{orig_name},{site},{'' if score is None else score}\n")

    # 5) Clean up ligand PDBQTs
    if not args.no_clean:
        for p in lig_dir.glob("*.pdbqt"):
            p.unlink(missing_ok=True)

    # 6) Sentinel and final summary
    (work_dir / "DONE.txt").write_text("ok", encoding="utf-8")
    
    log_performance(perf_log, 
        f"Chunk completed. Kept poses: {keep_counts}. "
        f"Optimal processes: {n_processes}. CSV: {csv_p}")

    print(f"[{chunk_id}] Done. Kept poses: {keep_counts}. CSV -> {csv_p}")
    print(f"[{chunk_id}] Performance log: {perf_log}")

if __name__ == "__main__":
    main()