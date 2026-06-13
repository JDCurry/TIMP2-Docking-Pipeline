#!/usr/bin/env python3
"""
Stage 2 of 3: Ensemble docking driver for TIMP2 Site 1.
=======================================================

Docks a FOCUSED ligand set into every receptor conformation from Stage 1, using
the SAME Site 1 grid as your primary screen but with elevated sampling
(exhaustiveness, multiple modes, multiple seeds) because the set is small enough
to afford it. The crystal baseline is docked at the identical elevated settings
so the crystal-vs-ensemble comparison is apples-to-apples.

Focused set (recommended): the seven MD-characterized compounds, because you
already have an MD-derived ordering (persistent > rebinding > transient >
dissociating) to interpret the result against:
  ZINC000583127796   persistent (thiirane,    Phe103-anchored, pocket-opening)
  ZINC001228612316   persistent (cyclopropane, crystal-fit)
  ZINC000826139982   rebinding  (-SCF3)
  ZINC000870248439   transient  (best Vina score, dissociates ~20 ns)
  ZINC000560333501   control    (non-benzamide cyanobenzyl)
  ZINC000329563636   dissociating (piperazine-urea)
  ZINC000543478179   dissociating (piperazine-urea)

NOTE: dock from NEUTRAL gen3d conformers (not your saved docked poses) so the
search is not biased toward the crystal-docked pose. Provide a .smi with one
"SMILES ZINCID" per line via --ligand-smi, OR a directory of prepared ligand
PDBQTs via --ligand-dir.

Vina CLI mirrors dock_adaptive_fixed.py:dock_one exactly except RECEPTOR varies
and exhaustiveness/num_modes/seed are configurable.

Dependencies: OpenBabel + Vina executables. (RDKit not required.)

Usage
-----
python dock_ensemble.py \
  --receptor-dir ensemble_receptors/receptor_pdbqt \
  --ligand-smi   focused7.smi \
  --out-dir      ensemble_docking \
  --vina vina --obabel obabel \
  --exhaustiveness 32 --num-modes 9 --seeds 1 2 3 4 5 --processes 16
"""

from __future__ import annotations
import argparse
import re
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

# Site 1 grid (identical to dock_adaptive_fixed.py)
# Updated for full 100ns ensemble: center on average Phe103 CB position across all receptors
SITE1 = dict(cx=25.15, cy=20.86, cz=13.75, sx=40, sy=40, sz=40)

VINA_RESULT_RE = re.compile(r"REMARK VINA RESULT:\s+(-?\d+(?:\.\d+)?)")


def sanitize(name: str) -> str:
    safe = re.sub(r'[<>:"/\\|?*]', "_", name)
    return re.sub(r"[\x00-\x1f\x7f]", "", safe).strip() or "lig"


def prep_ligands_from_smi(obabel: str, smi_path: Path, out_dir: Path, log_dir: Path):
    """One PDBQT per SMILES line via gen3d, matching your original ligand prep."""
    out_dir.mkdir(parents=True, exist_ok=True)
    made = []
    with open(smi_path, encoding="utf-8", errors="ignore") as fh:
        for i, raw in enumerate(fh, 1):
            line = raw.strip()
            if not line:
                continue
            parts = line.split()
            smi = parts[0]
            name = sanitize(parts[1] if len(parts) > 1 else f"lig_{i:04d}")
            out_pdbqt = out_dir / f"{name}.pdbqt"
            if out_pdbqt.exists():
                made.append(out_pdbqt)
                continue
            log = log_dir / f"{name}.obabel.log"
            cmd = [obabel, f"-:{smi}", "-O", str(out_pdbqt),
                   "--gen3d", "-h", "--partialcharge", "gasteiger", "--title", name]
            with open(log, "w", encoding="utf-8") as lf:
                lf.write(">> " + " ".join(cmd) + "\n")
                proc = subprocess.run(cmd, stdout=lf, stderr=lf)
            if proc.returncode == 0 and out_pdbqt.exists():
                made.append(out_pdbqt)
            else:
                print(f"  [WARN] ligand prep failed: {name}")
    return made


def parse_all_scores(pose_path: Path):
    """Return list of mode scores from a multi-mode Vina output PDBQT."""
    scores = []
    try:
        with open(pose_path, encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                m = VINA_RESULT_RE.search(line)
                if m:
                    scores.append(float(m.group(1)))
    except FileNotFoundError:
        pass
    return scores


def dock_one(vina: str, receptor: Path, ligand: Path, seed: int,
             exhaustiveness: int, num_modes: int, vina_threads: int,
             out_dir: Path, log_dir: Path):
    rid, lid = receptor.stem, ligand.stem
    out_p = out_dir / f"{rid}__{lid}__seed{seed}.pdbqt"
    log_p = log_dir / f"{rid}__{lid}__seed{seed}.log"
    if out_p.exists():
        return (rid, lid, seed, parse_all_scores(out_p), str(out_p))
    cmd = [
        vina, "--receptor", str(receptor), "--ligand", str(ligand),
        "--center_x", str(SITE1["cx"]), "--center_y", str(SITE1["cy"]),
        "--center_z", str(SITE1["cz"]),
        "--size_x", str(SITE1["sx"]), "--size_y", str(SITE1["sy"]),
        "--size_z", str(SITE1["sz"]),
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(num_modes),
        "--seed", str(seed),
        "--cpu", str(vina_threads),
        "--out", str(out_p),
        "--log", str(log_p),
    ]
    proc = subprocess.run(cmd)
    if proc.returncode != 0:
        return (rid, lid, seed, [], str(out_p))
    return (rid, lid, seed, parse_all_scores(out_p), str(out_p))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--receptor-dir", required=True, type=Path)
    ap.add_argument("--ligand-smi", type=Path, help="focused set as SMILES (preferred)")
    ap.add_argument("--ligand-dir", type=Path, help="dir of prepared ligand PDBQTs")
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--vina", default="vina")
    ap.add_argument("--obabel", default="obabel")
    ap.add_argument("--exhaustiveness", type=int, default=32)
    ap.add_argument("--num-modes", type=int, default=9)
    ap.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3, 4, 5])
    ap.add_argument("--vina-threads", type=int, default=2)
    ap.add_argument("--processes", type=int, default=8)
    args = ap.parse_args()

    out_dir = args.out_dir.resolve()
    pose_dir = out_dir / "poses"
    log_dir = out_dir / "logs"
    lig_dir = out_dir / "ligand_pdbqt"
    for d in (pose_dir, log_dir, lig_dir):
        d.mkdir(parents=True, exist_ok=True)

    receptors = sorted(args.receptor_dir.glob("*.pdbqt"))
    if not receptors:
        sys.exit(f"No receptor PDBQTs in {args.receptor_dir}")

    if args.ligand_smi:
        print("[1] Preparing ligands from SMILES (neutral gen3d) ...")
        ligands = prep_ligands_from_smi(args.obabel, args.ligand_smi, lig_dir, log_dir)
    elif args.ligand_dir:
        ligands = sorted(args.ligand_dir.glob("*.pdbqt"))
    else:
        sys.exit("Provide --ligand-smi or --ligand-dir")
    if not ligands:
        sys.exit("No ligands to dock.")

    n_jobs = len(receptors) * len(ligands) * len(args.seeds)
    print(f"[2] Docking {len(ligands)} ligands x {len(receptors)} receptors "
          f"x {len(args.seeds)} seeds = {n_jobs} jobs "
          f"(exhaustiveness {args.exhaustiveness}, num_modes {args.num_modes})")

    results = []
    with ProcessPoolExecutor(max_workers=args.processes) as ex:
        futs = []
        for rec in receptors:
            for lig in ligands:
                for seed in args.seeds:
                    futs.append(ex.submit(
                        dock_one, args.vina, rec, lig, seed,
                        args.exhaustiveness, args.num_modes, args.vina_threads,
                        pose_dir, log_dir))
        done = 0
        for f in as_completed(futs):
            rid, lid, seed, scores, posep = f.result()
            best = min(scores) if scores else ""
            for rank, sc in enumerate(scores, 1):
                results.append((rid, lid, seed, rank, sc, posep))
            if not scores:
                results.append((rid, lid, seed, 0, "", posep))
            done += 1
            if done % 25 == 0 or done == n_jobs:
                print(f"    {done}/{n_jobs} done")

    csv_p = out_dir / "ensemble_scores.csv"
    with open(csv_p, "w", encoding="utf-8") as fh:
        fh.write("receptor,ligand,seed,mode_rank,score_kcal_per_mol,pose_path\n")
        for rid, lid, seed, rank, sc, posep in sorted(results):
            fh.write(f"{rid},{lid},{seed},{rank},{sc},{posep}\n")
    print(f"\n[3] Wrote scores: {csv_p}")
    print(f"    Poses in: {pose_dir}")
    print("    Next: analyze_ensemble.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
