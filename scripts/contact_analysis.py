#!/usr/bin/env python3
"""
Contact Occupancy & Minimum Distance Analysis for TIMP2 Ligand-Bound MD
========================================================================

Reads existing trajectory files and computes:
  1. Per-residue contact occupancy (% frames with any heavy atom < cutoff)
  2. Minimum heavy-atom distance between ligand and pocket over time
  3. Per-residue minimum distance over time
  4. Contact frequency heatmap

These metrics are PBC-immune for short-range contacts and definitively
answer whether the ligand is in the pocket regardless of COM artifacts.

Usage:
    python3 contact_analysis.py --dir ~/Desktop/timp2/site1_analog_100ns --ligand site1_analog
    python3 contact_analysis.py --dir ~/Desktop/timp2/best_binder_100ns --ligand best_binder
    
    # Run on all directories at once for comparison:
    python3 contact_analysis.py --compare \
        ~/Desktop/timp2/site1_analog_100ns \
        ~/Desktop/timp2/best_binder_100ns \
        ~/Desktop/timp2/site1_analog_10ns \
        ~/Desktop/timp2/best_binder_10ns

Output:
    contact_analysis.json    - All numerical results
    contact_analysis.png     - Publication figure (4 panels)
    contact_analysis.pdf     - Vector version
"""

import os
import sys
import argparse
import json
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import MDAnalysis as mda
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from openmm.app import PDBFile, PDBxFile

# ============================================================
# ARGUMENTS
# ============================================================
parser = argparse.ArgumentParser(description='Contact occupancy analysis')
parser.add_argument('--dir', type=str, help='Single trajectory directory')
parser.add_argument('--ligand', type=str, default=None, help='Ligand name for title')
parser.add_argument('--cutoff', type=float, default=4.0, help='Contact distance cutoff in Å (default: 4.0)')
parser.add_argument('--pdb', type=str, default='md_start.pdb', help='Reference structure filename')
parser.add_argument('--dcd', type=str, default='md_production.dcd', help='Trajectory filename')
parser.add_argument('--compare', nargs='+', help='Multiple directories for comparison')
args = parser.parse_args()

# ============================================================
# POCKET DEFINITION
# ============================================================
TRIAD_RESIDS = [6, 100, 103]
TRIAD_NAMES = {6: 'Val6', 100: 'Leu100', 103: 'Phe103'}
SITE1_RESIDS = [3, 4, 5, 6, 8, 11, 99, 100, 101, 102, 103, 104, 105]

CRYSTAL = {
    'd_6_100': 11.98, 'd_6_103': 8.22, 'd_100_103': 7.61
}


def kabsch_align(mobile, reference):
    """Kabsch alignment: returns R, t such that R @ mobile + t ≈ reference."""
    mob_c = mobile - mobile.mean(axis=0)
    ref_c = reference - reference.mean(axis=0)
    H = mob_c.T @ ref_c
    U, S, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    sign = np.diag([1, 1, np.sign(d)])
    R = Vt.T @ sign @ U.T
    t = reference.mean(axis=0) - R @ mobile.mean(axis=0)
    return R, t


def resolve_structure_paths(work_dir, pdb_name='md_start.pdb', dcd_name='md_production.dcd'):
    """Resolve topology/trajectory filenames with sensible fallbacks."""
    candidates = [
        os.path.join(work_dir, pdb_name),
        os.path.join(work_dir, pdb_name.replace('.pdb', '.cif')),
    ]

    if pdb_name == 'md_start.pdb':
        stem = os.path.basename(work_dir.rstrip('/'))
        candidates.extend([
            os.path.join(work_dir, f'{stem}_solvated.pdb'),
            os.path.join(work_dir, f'{stem}.pdb'),
            os.path.join(work_dir, f'{stem}.cif'),
        ])

    topology_path = next((path for path in candidates if os.path.exists(path)), None)
    trajectory_path = os.path.join(work_dir, dcd_name)
    return topology_path, trajectory_path


def analyze_contacts(work_dir, ligand_name=None, cutoff=4.0, pdb_name='md_start.pdb', dcd_name='md_production.dcd'):
    """Run full contact analysis on a single trajectory."""
    
    pdb_path, dcd_path = resolve_structure_paths(work_dir, pdb_name=pdb_name, dcd_name=dcd_name)
    
    if pdb_path is None:
        print(f"  ERROR: No topology file found in {work_dir}")
        return None
    if not os.path.exists(dcd_path):
        print(f"  ERROR: {dcd_path} not found")
        return None
    
    print(f"\n  Loading: {os.path.basename(pdb_path)} + {os.path.basename(dcd_path)}")
    topology_source = PDBxFile(pdb_path) if pdb_path.lower().endswith('.cif') else PDBFile(pdb_path)
    u = mda.Universe(topology_source, dcd_path)
    n_frames = len(u.trajectory)
    total_ns = u.trajectory.totaltime / 1000
    print(f"  Frames: {n_frames}, Duration: {total_ns:.1f} ns")
    
    # Selections - heavy atoms only (no hydrogens)
    ligand_heavy = u.select_atoms(
        'not protein and not resname HOH WAT TIP3 TIP3P SOL NA CL K and not name H*'
    )
    
    if len(ligand_heavy) == 0:
        print(f"  ERROR: No ligand heavy atoms found")
        return None
    
    print(f"  Ligand heavy atoms: {len(ligand_heavy)}")
    
    # Pocket residue heavy atoms
    pocket_sels = {}
    for resid in SITE1_RESIDS:
        sel = u.select_atoms(f'resid {resid} and not name H*')
        if len(sel) > 0:
            pocket_sels[resid] = sel
    
    # Protein CA for alignment
    protein_ca = u.select_atoms("protein and name CA")
    u.trajectory[0]
    ref_ca_pos = protein_ca.positions.copy()
    
    print(f"  Pocket residues found: {len(pocket_sels)}")
    print(f"  Analyzing contacts (cutoff = {cutoff} Å)...")
    
    # Storage
    times = []
    min_dist_per_frame = []           # global minimum ligand-pocket distance
    triad_min_dist = {r: [] for r in TRIAD_RESIDS}  # per-triad-residue min distance
    all_residue_min_dist = {r: [] for r in pocket_sels}  # all pocket residues
    contact_frames = {r: 0 for r in pocket_sels}     # count frames in contact
    
    for ts in u.trajectory:
        times.append(ts.time / 1000)  # ns
        
        # Align to remove PBC drift
        R, t = kabsch_align(protein_ca.positions, ref_ca_pos)
        aligned_lig = (R @ ligand_heavy.positions.T).T + t
        
        frame_min = np.inf
        
        for resid, sel in pocket_sels.items():
            aligned_res = (R @ sel.positions.T).T + t
            
            # Compute all pairwise distances between ligand and residue atoms
            # Shape: (n_lig_atoms, n_res_atoms)
            diff = aligned_lig[:, np.newaxis, :] - aligned_res[np.newaxis, :, :]
            dists = np.sqrt(np.sum(diff**2, axis=2))
            min_d = dists.min()
            
            all_residue_min_dist[resid].append(float(min_d))
            
            if resid in TRIAD_RESIDS:
                triad_min_dist[resid].append(float(min_d))
            
            if min_d < cutoff:
                contact_frames[resid] += 1
            
            if min_d < frame_min:
                frame_min = min_d
        
        min_dist_per_frame.append(float(frame_min))
    
    times = np.array(times)
    min_dist_per_frame = np.array(min_dist_per_frame)
    
    # Compute occupancies
    occupancy = {r: contact_frames[r] / n_frames * 100 for r in pocket_sels}
    
    # Summary stats
    triad_occ = {r: occupancy.get(r, 0) for r in TRIAD_RESIDS}
    any_contact = np.array(min_dist_per_frame) < cutoff
    overall_contact_pct = np.sum(any_contact) / n_frames * 100
    
    # Print results
    print(f"\n  {'='*60}")
    print(f"  CONTACT ANALYSIS RESULTS ({cutoff} Å cutoff)")
    print(f"  {'='*60}")
    print(f"\n  Overall pocket contact: {overall_contact_pct:.1f}% of frames")
    print(f"  Min distance mean: {np.mean(min_dist_per_frame):.2f} ± {np.std(min_dist_per_frame):.2f} Å")
    print(f"  Min distance median: {np.median(min_dist_per_frame):.2f} Å")
    print(f"  Min distance range: {np.min(min_dist_per_frame):.2f} - {np.max(min_dist_per_frame):.2f} Å")
    
    print(f"\n  Per-residue contact occupancy:")
    print(f"  {'Residue':>12} {'Occupancy':>10} {'Min Dist Mean':>14} {'Note':>10}")
    print(f"  {'-'*50}")
    for resid in sorted(pocket_sels.keys()):
        name = TRIAD_NAMES.get(resid, f'Res{resid}')
        occ = occupancy[resid]
        mean_d = np.mean(all_residue_min_dist[resid])
        triad_flag = " ← TRIAD" if resid in TRIAD_RESIDS else ""
        print(f"  {name:>12} {occ:>9.1f}% {mean_d:>13.2f} Å{triad_flag}")
    
    print(f"\n  Hydrophobic triad summary:")
    for resid in TRIAD_RESIDS:
        name = TRIAD_NAMES[resid]
        occ = triad_occ[resid]
        dists = np.array(triad_min_dist[resid])
        print(f"    {name}: {occ:.1f}% occupancy, min dist {np.mean(dists):.2f} ± {np.std(dists):.2f} Å")
    
    # Verdict
    mean_triad_occ = np.mean([triad_occ[r] for r in TRIAD_RESIDS])
    if mean_triad_occ > 50:
        verdict = "BOUND (strong triad contact)"
    elif mean_triad_occ > 25:
        verdict = "ASSOCIATED (partial triad contact)"
    elif overall_contact_pct > 50:
        verdict = "PERIPHERAL (pocket contact, weak triad)"
    else:
        verdict = "UNBOUND (minimal contact)"
    
    print(f"\n  VERDICT: {verdict}")
    print(f"  Mean triad occupancy: {mean_triad_occ:.1f}%")
    print(f"  {'='*60}")
    
    # ============================================================
    # BUILD RESULTS DICT
    # ============================================================
    results = {
        'ligand': ligand_name or os.path.basename(work_dir),
        'simulation_ns': float(total_ns),
        'n_frames': n_frames,
        'cutoff_A': cutoff,
        'overall_contact_pct': float(overall_contact_pct),
        'min_distance': {
            'mean': float(np.mean(min_dist_per_frame)),
            'std': float(np.std(min_dist_per_frame)),
            'median': float(np.median(min_dist_per_frame)),
            'min': float(np.min(min_dist_per_frame)),
            'max': float(np.max(min_dist_per_frame)),
        },
        'triad_occupancy': {TRIAD_NAMES[r]: float(triad_occ[r]) for r in TRIAD_RESIDS},
        'triad_mean_occupancy': float(mean_triad_occ),
        'triad_min_distance': {
            TRIAD_NAMES[r]: {
                'mean': float(np.mean(triad_min_dist[r])),
                'std': float(np.std(triad_min_dist[r])),
            } for r in TRIAD_RESIDS
        },
        'all_residue_occupancy': {
            f"Res{r}": float(occupancy[r]) for r in sorted(pocket_sels.keys())
        },
        'verdict': verdict,
    }
    
    # ============================================================
    # FIGURE
    # ============================================================
    print(f"\n  Generating figure...")
    
    fig = plt.figure(figsize=(16, 14))
    gs = GridSpec(3, 2, figure=fig, hspace=0.4, wspace=0.3,
                  height_ratios=[1, 1, 0.8])
    
    window = max(1, len(times) // 100)
    
    def rmean(x, w):
        return np.convolve(x, np.ones(w)/w, mode='valid') if len(x) >= w else x
    
    title = f"TIMP2 Contact Analysis — {ligand_name or os.path.basename(work_dir)}\n{total_ns:.0f} ns MD, {cutoff} Å cutoff"
    
    # --- Panel A: Minimum distance over time ---
    ax = fig.add_subplot(gs[0, :])
    ax.plot(times, min_dist_per_frame, alpha=0.2, color='steelblue', linewidth=0.3)
    ax.plot(times[window-1:], rmean(min_dist_per_frame, window), color='navy', linewidth=2,
            label='Running average')
    ax.axhline(y=cutoff, color='red', linestyle='--', linewidth=1.5, alpha=0.7,
               label=f'Contact cutoff ({cutoff} Å)')
    ax.fill_between(times, 0, cutoff, alpha=0.05, color='green')
    ax.set_ylabel('Min Heavy-Atom Distance (Å)', fontsize=12)
    ax.set_xlabel('Time (ns)', fontsize=12)
    ax.set_title('(A) Minimum Ligand–Pocket Heavy-Atom Distance', fontsize=13, fontweight='bold')
    ax.legend(fontsize=10, loc='upper right')
    ax.grid(alpha=0.2)
    ax.set_xlim(0, times[-1])
    ax.set_ylim(0, min(max(min_dist_per_frame) * 1.1, 30))
    
    # Annotate contact percentage
    ax.text(0.02, 0.95, f'Contact: {overall_contact_pct:.0f}% of frames',
            transform=ax.transAxes, fontsize=11, fontweight='bold',
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='lightgreen' if overall_contact_pct > 50 else 'lightyellow',
                      alpha=0.8))
    
    # --- Panel B: Per-triad-residue min distance over time ---
    ax = fig.add_subplot(gs[1, :])
    colors = {'Val6': '#e41a1c', 'Leu100': '#377eb8', 'Phe103': '#4daf4a'}
    for resid in TRIAD_RESIDS:
        name = TRIAD_NAMES[resid]
        data = np.array(triad_min_dist[resid])
        ax.plot(times, data, alpha=0.15, color=colors[name], linewidth=0.3)
        ax.plot(times[window-1:], rmean(data, window), color=colors[name], linewidth=2,
                label=f'{name} ({triad_occ[resid]:.0f}% contact)')
    ax.axhline(y=cutoff, color='red', linestyle='--', linewidth=1.5, alpha=0.5)
    ax.fill_between(times, 0, cutoff, alpha=0.05, color='green')
    ax.set_ylabel('Min Heavy-Atom Distance (Å)', fontsize=12)
    ax.set_xlabel('Time (ns)', fontsize=12)
    ax.set_title('(B) Triad Residue–Ligand Minimum Distance', fontsize=13, fontweight='bold')
    ax.legend(fontsize=10, loc='upper right')
    ax.grid(alpha=0.2)
    ax.set_xlim(0, times[-1])
    ax.set_ylim(0, min(max(max(triad_min_dist[r]) for r in TRIAD_RESIDS) * 1.1, 30))
    
    # --- Panel C: Contact occupancy bar chart ---
    ax = fig.add_subplot(gs[2, 0])
    resids_sorted = sorted(pocket_sels.keys())
    occ_values = [occupancy[r] for r in resids_sorted]
    bar_colors = ['#e41a1c' if r in TRIAD_RESIDS else '#4393c3' for r in resids_sorted]
    
    labels = []
    for r in resids_sorted:
        sel = u.select_atoms(f"name CA and resid {r}")
        if len(sel) > 0:
            labels.append(f"{sel.residues[0].resname}{r}")
        else:
            labels.append(f"R{r}")
    
    x_pos = np.arange(len(resids_sorted))
    ax.bar(x_pos, occ_values, color=bar_colors, alpha=0.8, edgecolor='white')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=9)
    ax.set_ylabel('Contact Occupancy (%)', fontsize=11)
    ax.set_title('(C) Per-Residue Contact Occupancy', fontsize=12, fontweight='bold')
    ax.axhline(y=50, color='gray', linestyle=':', alpha=0.5, label='50% threshold')
    
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#e41a1c', alpha=0.8, label='Hydrophobic triad'),
        Patch(facecolor='#4393c3', alpha=0.8, label='Other pocket residues'),
    ]
    ax.legend(handles=legend_elements, fontsize=9, loc='upper right')
    ax.grid(alpha=0.2, axis='y')
    ax.set_ylim(0, 105)
    
    # --- Panel D: Min distance distribution ---
    ax = fig.add_subplot(gs[2, 1])
    ax.hist(min_dist_per_frame, bins=80, color='steelblue', alpha=0.7,
            edgecolor='navy', linewidth=0.3, density=True)
    ax.axvline(x=cutoff, color='red', linestyle='--', linewidth=2, label=f'Cutoff ({cutoff} Å)')
    ax.axvline(x=np.median(min_dist_per_frame), color='orange', linestyle='--', linewidth=2,
               label=f'Median ({np.median(min_dist_per_frame):.1f} Å)')
    ax.set_xlabel('Min Heavy-Atom Distance (Å)', fontsize=11)
    ax.set_ylabel('Density', fontsize=11)
    ax.set_title('(D) Distance Distribution', fontsize=12, fontweight='bold')
    ax.legend(fontsize=9)
    ax.set_xlim(0, min(max(min_dist_per_frame) * 1.1, 25))
    
    fig.suptitle(title, fontsize=15, fontweight='bold', y=0.99)
    
    # Save
    png_path = os.path.join(work_dir, 'contact_analysis.png')
    pdf_path = os.path.join(work_dir, 'contact_analysis.pdf')
    fig.savefig(png_path, dpi=300, bbox_inches='tight')
    fig.savefig(pdf_path, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {png_path}")
    print(f"  Saved: {pdf_path}")
    
    # Save JSON
    json_path = os.path.join(work_dir, 'contact_analysis.json')
    with open(json_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"  Saved: {json_path}")
    
    return results


def compare_runs(directories):
    """Run contact analysis on multiple directories and generate comparison figure."""
    
    all_results = []
    for d in directories:
        name = os.path.basename(d.rstrip('/'))
        print(f"\n{'='*60}")
        print(f"Analyzing: {name}")
        print(f"{'='*60}")
        result = analyze_contacts(d, ligand_name=name)
        if result:
            all_results.append(result)
    
    if len(all_results) < 2:
        print("Need at least 2 results for comparison")
        return
    
    # Comparison figure
    print(f"\nGenerating comparison figure...")
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    names = [r['ligand'] for r in all_results]
    x = np.arange(len(names))
    
    # Panel 1: Overall contact %
    ax = axes[0]
    vals = [r['overall_contact_pct'] for r in all_results]
    colors = ['green' if v > 50 else 'orange' if v > 25 else 'red' for v in vals]
    ax.bar(x, vals, color=colors, alpha=0.8, edgecolor='black')
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=30, ha='right', fontsize=9)
    ax.set_ylabel('Contact Occupancy (%)')
    ax.set_title('Overall Pocket Contact', fontweight='bold')
    ax.axhline(y=50, color='gray', linestyle=':', alpha=0.5)
    ax.set_ylim(0, 105)
    
    # Panel 2: Triad occupancy
    ax = axes[1]
    width = 0.25
    for i, resid in enumerate(TRIAD_RESIDS):
        name_r = TRIAD_NAMES[resid]
        vals = [r['triad_occupancy'].get(name_r, 0) for r in all_results]
        ax.bar(x + i * width - width, vals, width, label=name_r, alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=30, ha='right', fontsize=9)
    ax.set_ylabel('Contact Occupancy (%)')
    ax.set_title('Triad Residue Contact', fontweight='bold')
    ax.legend(fontsize=9)
    ax.set_ylim(0, 105)
    
    # Panel 3: Median min distance
    ax = axes[2]
    vals = [r['min_distance']['median'] for r in all_results]
    colors = ['green' if v < 4 else 'orange' if v < 8 else 'red' for v in vals]
    ax.bar(x, vals, color=colors, alpha=0.8, edgecolor='black')
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=30, ha='right', fontsize=9)
    ax.set_ylabel('Median Min Distance (Å)')
    ax.set_title('Median Ligand–Pocket Distance', fontweight='bold')
    ax.axhline(y=4.0, color='red', linestyle=':', alpha=0.5, label='4 Å cutoff')
    ax.legend(fontsize=9)
    
    fig.suptitle('TIMP2 Ligand Contact Analysis — Compound Comparison', 
                 fontsize=14, fontweight='bold')
    fig.tight_layout()
    
    out_dir = os.path.dirname(directories[0].rstrip('/'))
    png_path = os.path.join(out_dir, 'contact_comparison.png')
    pdf_path = os.path.join(out_dir, 'contact_comparison.pdf')
    fig.savefig(png_path, dpi=300, bbox_inches='tight')
    fig.savefig(pdf_path, bbox_inches='tight')
    plt.close()
    print(f"\nSaved: {png_path}")
    print(f"Saved: {pdf_path}")
    
    # Summary table
    print(f"\n{'='*80}")
    print(f"COMPARISON SUMMARY")
    print(f"{'='*80}")
    print(f"{'Compound':<25} {'Contact%':>9} {'Val6':>7} {'Leu100':>7} {'Phe103':>7} {'Med Dist':>9} {'Verdict'}")
    print(f"{'-'*80}")
    for r in sorted(all_results, key=lambda x: -x['triad_mean_occupancy']):
        print(f"{r['ligand']:<25} {r['overall_contact_pct']:>8.1f}% "
              f"{r['triad_occupancy']['Val6']:>6.1f}% "
              f"{r['triad_occupancy']['Leu100']:>6.1f}% "
              f"{r['triad_occupancy']['Phe103']:>6.1f}% "
              f"{r['min_distance']['median']:>8.2f}  "
              f"{r['verdict']}")


# ============================================================
# MAIN
# ============================================================
if __name__ == '__main__':
    if args.compare:
        compare_runs(args.compare)
    elif args.dir:
        analyze_contacts(args.dir, ligand_name=args.ligand, cutoff=args.cutoff,
                         pdb_name=args.pdb, dcd_name=args.dcd)
    else:
        print("Provide --dir for single analysis or --compare for multiple directories")
        sys.exit(1)
