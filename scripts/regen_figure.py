#!/usr/bin/env python3
"""
Regenerate pocket stability figure with RMSF panel replacing RMSD.
================================================================

Reads existing trajectory files (md_start.pdb + md_production.dcd)
and regenerates the 5-panel publication figure with:
  A: Triad distances over time
  B: Pocket triangle area over time
  C: Per-residue RMSF for pocket region (NEW - replaces broken RMSD)
  D: Area distribution histogram
  E: Distance violin plots

Usage (apo):
    python3 regen_figure.py --dir /path/to/apo_100ns

Usage (ligand-bound):
    python3 regen_figure.py --dir /path/to/best_binder_100ns --ligand best_binder

The script auto-detects ligand atoms and adds Panels F+G if present.
"""

import os
import sys
import argparse
import json
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import MDAnalysis as mda
from MDAnalysis import transformations as tr
from MDAnalysis.lib.distances import minimize_vectors
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from openmm.app import PDBFile, PDBxFile

# ============================================================
# ARGUMENTS
# ============================================================
parser = argparse.ArgumentParser(description='Regenerate MD pocket figure')
parser.add_argument('--dir', type=str, required=True, help='Directory with md_start.pdb and md_production.dcd')
parser.add_argument('--ligand', type=str, default=None, help='Ligand name (if ligand-bound run)')
parser.add_argument('--pdb', type=str, default='md_start.pdb', help='Reference PDB filename')
parser.add_argument('--dcd', type=str, default='md_production.dcd', help='Trajectory DCD filename')
args = parser.parse_args()

WORK = args.dir
pdb_path = os.path.join(WORK, args.pdb)
dcd_path = os.path.join(WORK, args.dcd)

# Try CIF if PDB doesn't exist
if not os.path.exists(pdb_path) and os.path.exists(pdb_path.replace('.pdb', '.cif')):
    pdb_path = pdb_path.replace('.pdb', '.cif')
    print(f"  Using CIF: {pdb_path}")

for f in [pdb_path, dcd_path]:
    if not os.path.exists(f):
        print(f"ERROR: {f} not found")
        sys.exit(1)

# ============================================================
# LOAD TRAJECTORY
# ============================================================
print(f"Loading trajectory...")
topology_source = PDBxFile(pdb_path) if pdb_path.lower().endswith('.cif') else PDBFile(pdb_path)
u = mda.Universe(topology_source, dcd_path)
u.trajectory.add_transformations(tr.unwrap(u.atoms))
print(f"  {len(u.trajectory)} frames, {u.trajectory.totaltime/1000:.1f} ns")

# ============================================================
# SELECTIONS
# ============================================================
triad_resids = [6, 100, 103]
site1_resids = [3, 4, 5, 6, 8, 11, 99, 100, 101, 102, 103, 104, 105]

triad_ca = u.select_atoms(f"name CA and resid {' '.join(str(r) for r in triad_resids)}")
site1_ca = u.select_atoms(f"name CA and resid {' '.join(str(r) for r in site1_resids)}")
pocket_all = u.select_atoms(f"resid {' '.join(str(r) for r in site1_resids)}")

# Ligand detection
ligand_sel = u.select_atoms('not protein and not resname HOH WAT TIP3 TIP3P SOL NA CL K')
has_ligand = len(ligand_sel) > 0 and args.ligand is not None
if has_ligand:
    print(f"  Ligand atoms: {len(ligand_sel)}")
else:
    print(f"  Apo simulation (no ligand)")

print(f"  Triad CA: {len(triad_ca)}, Site1 CA: {len(site1_ca)}, Pocket all-atom: {len(pocket_all)}")

# Crystal references
crystal = {'d_6_100': 11.98, 'd_6_103': 8.22, 'd_100_103': 7.61}

def tri_area(a, b, c):
    s = (a + b + c) / 2
    return np.sqrt(np.maximum(s * (s - a) * (s - b) * (s - c), 0))

crystal_area = tri_area(crystal['d_6_100'], crystal['d_6_103'], crystal['d_100_103'])

# ============================================================
# COLLECT PER-FRAME DATA
# ============================================================
print("Analyzing frames...")
times = []
d_6_100, d_6_103, d_100_103 = [], [], []
ligand_rmsds, lig_pocket_dists = [], []

# For RMSF: collect all positions
pocket_positions_all = []  # list of (n_frames, n_atoms, 3)

# For alignment: use all protein CA atoms
protein_ca = u.select_atoms("protein and name CA")

# Reference frame
u.trajectory[0]
ref_protein_ca_pos = protein_ca.positions.copy()
if has_ligand:
    ref_lig_pos = ligand_sel.positions.copy()
    ref_lig_com = ligand_sel.center_of_mass()

def kabsch_align(mobile, reference):
    """Align mobile onto reference using Kabsch algorithm. Returns rotation matrix and translation."""
    mob_center = mobile.mean(axis=0)
    ref_center = reference.mean(axis=0)
    mob_c = mobile - mob_center
    ref_c = reference - ref_center
    H = mob_c.T @ ref_c
    U, S, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    sign = np.diag([1, 1, np.sign(d)])
    R = Vt.T @ sign @ U.T
    t = ref_center - R @ mob_center
    return R, t

for ts in u.trajectory:
    times.append(ts.time / 1000)  # ps -> ns
    
    # Triad distances
    pos = {}
    for atom in triad_ca:
        pos[atom.resid] = atom.position
    if len(pos) == 3:
        d_6_100.append(np.linalg.norm(pos[triad_resids[0]] - pos[triad_resids[1]]))
        d_6_103.append(np.linalg.norm(pos[triad_resids[0]] - pos[triad_resids[2]]))
        d_100_103.append(np.linalg.norm(pos[triad_resids[1]] - pos[triad_resids[2]]))
    
    # Align this frame's protein backbone onto reference (fixes PBC drift)
    R, t = kabsch_align(protein_ca.positions, ref_protein_ca_pos)
    
    # Apply transform to pocket CA positions
    aligned_pocket = (R @ site1_ca.positions.T).T + t
    pocket_positions_all.append(aligned_pocket.copy())
    
    # Ligand tracking (also align)
    if has_ligand:
        aligned_lig = (R @ ligand_sel.positions.T).T + t
        aligned_pocket_com = aligned_pocket.mean(axis=0)
        aligned_lig_com = aligned_lig.mean(axis=0)

        # Reimage the ligand to the nearest periodic image of the pocket COM.
        # This prevents artificial 100+ Å spikes when the ligand wraps across the box.
        if ts.dimensions is not None:
            lig_vec = aligned_lig_com - aligned_pocket_com
            lig_vec_min = minimize_vectors(lig_vec[None, :], ts.dimensions)[0]
            aligned_lig = aligned_lig - (lig_vec - lig_vec_min)
            aligned_lig_com = aligned_lig.mean(axis=0)

        lig_rmsd = np.sqrt(np.mean(np.sum((aligned_lig - ref_lig_pos)**2, axis=1)))
        ligand_rmsds.append(lig_rmsd)
        lig_pocket_dists.append(np.linalg.norm(aligned_lig_com - aligned_pocket_com))

times = np.array(times)
d_6_100, d_6_103, d_100_103 = np.array(d_6_100), np.array(d_6_103), np.array(d_100_103)
areas = tri_area(d_6_100, d_6_103, d_100_103)

# ============================================================
# COMPUTE PER-RESIDUE RMSF
# ============================================================
print("Computing per-residue RMSF...")
pocket_positions_all = np.array(pocket_positions_all)  # (n_frames, n_atoms, 3)
mean_positions = pocket_positions_all.mean(axis=0)  # (n_atoms, 3)
rmsf_per_atom = np.sqrt(np.mean(np.sum((pocket_positions_all - mean_positions)**2, axis=2), axis=0))

# Map atoms to residues
pocket_resids = site1_ca.resids
pocket_resnames = [r.resname for r in site1_ca.residues]
rmsf_per_residue = rmsf_per_atom  # Already per-CA since we selected CA only

# Identify triad residues for highlighting
triad_mask = np.isin(pocket_resids, triad_resids)

print(f"  RMSF range: {rmsf_per_residue.min():.2f} - {rmsf_per_residue.max():.2f} Å")
print(f"  Triad RMSF: ", end="")
for rid, rmsf in zip(pocket_resids[triad_mask], rmsf_per_residue[triad_mask]):
    print(f"Res{rid}={rmsf:.2f}Å  ", end="")
print()

# ============================================================
# STATISTICS
# ============================================================
area_pres = abs(np.mean(areas) - crystal_area) / crystal_area
cv = np.std(areas) / np.mean(areas)

print(f"\n  RESULTS:")
print(f"  {'Measurement':>25} {'Mean':>8} {'StdDev':>8} {'Crystal':>8}")
print(f"  {'-'*55}")
print(f"  {'Val6-Leu100 (Å)':>25} {np.mean(d_6_100):>8.2f} {np.std(d_6_100):>8.2f} {crystal['d_6_100']:>8.2f}")
print(f"  {'Val6-Phe103 (Å)':>25} {np.mean(d_6_103):>8.2f} {np.std(d_6_103):>8.2f} {crystal['d_6_103']:>8.2f}")
print(f"  {'Leu100-Phe103 (Å)':>25} {np.mean(d_100_103):>8.2f} {np.std(d_100_103):>8.2f} {crystal['d_100_103']:>8.2f}")
print(f"  {'Triangle area (Å²)':>25} {np.mean(areas):>8.2f} {np.std(areas):>8.2f} {crystal_area:>8.2f}")
print(f"  Area preservation: {(1-area_pres)*100:.0f}%, CV: {cv*100:.0f}%")

verdict = "STABLE" if area_pres < 0.20 and cv < 0.30 else ("DYNAMIC" if area_pres < 0.40 else "UNSTABLE")
print(f"  Verdict: {verdict}")

if has_ligand:
    ligand_rmsds = np.array(ligand_rmsds)
    lig_pocket_dists = np.array(lig_pocket_dists)
    print(f"\n  Ligand RMSD:     {np.mean(ligand_rmsds):.1f} ± {np.std(ligand_rmsds):.1f} Å")
    print(f"  Pocket-lig dist: {np.mean(lig_pocket_dists):.1f} ± {np.std(lig_pocket_dists):.1f} Å")
    bound = np.mean(lig_pocket_dists) < 12.0
    print(f"  Ligand status:   {'BOUND' if bound else 'UNBOUND'}")

# ============================================================
# FIGURE
# ============================================================
print(f"\nGenerating publication figure...")

n_rows = 4 if not has_ligand else 5
fig = plt.figure(figsize=(14, 4 * n_rows))
gs = GridSpec(n_rows, 2, figure=fig, hspace=0.4, wspace=0.3,
              height_ratios=[1, 1, 0.8] + [0.8] + ([1] if has_ligand else []))

window = max(1, len(times) // 100)

def rmean(x, w):
    return np.convolve(x, np.ones(w)/w, mode='valid') if len(x) >= w else x

# --- Title ---
sim_label = f"{times[-1]:.0f} ns"
if args.ligand:
    title = f"TIMP2 Cryptic Pocket + Ligand Stability\n{args.ligand}\n{sim_label} MD (AMBER14/TIP3P-FB, 300 K, 1 bar, 0.15 M NaCl)"
else:
    title = f"TIMP2 Cryptic Allosteric Pocket Stability\n{sim_label} MD Simulation (AMBER14/TIP3P-FB, 300 K, 1 bar, 0.15 M NaCl)"

# --- Panel A: Triad distances ---
ax = fig.add_subplot(gs[0, :])
for data, label, cry, color in [
    (d_6_100, 'Val6–Leu100', crystal['d_6_100'], '#1f77b4'),
    (d_6_103, 'Val6–Phe103', crystal['d_6_103'], '#ff7f0e'),
    (d_100_103, 'Leu100–Phe103', crystal['d_100_103'], '#2ca02c'),
]:
    ax.plot(times, data, alpha=0.15, color=color, linewidth=0.3)
    ax.plot(times[window-1:], rmean(data, window), color=color, linewidth=2, label=label)
    ax.axhline(y=cry, color=color, linestyle='--', alpha=0.4, linewidth=1)
ax.set_ylabel('Distance (Å)', fontsize=12)
ax.set_title('(A) Hydrophobic Triad Cα Distances', fontsize=13, fontweight='bold')
ax.legend(fontsize=10, loc='upper right')
ax.grid(alpha=0.2)
ax.set_xlim(0, times[-1])

# --- Panel B: Pocket area ---
ax = fig.add_subplot(gs[1, :])
ax.fill_between(times, areas, alpha=0.15, color='steelblue')
ax.plot(times, areas, alpha=0.2, color='steelblue', linewidth=0.3)
ax.plot(times[window-1:], rmean(areas, window), color='navy', linewidth=2, label='Running average')
ax.axhline(y=crystal_area, color='red', linestyle='--', linewidth=1.5,
           label=f'Crystal structure ({crystal_area:.1f} Å²)', alpha=0.7)
ax.axhline(y=np.mean(areas), color='orange', linestyle=':', linewidth=1.5,
           label=f'MD mean ({np.mean(areas):.1f} ± {np.std(areas):.1f} Å²)', alpha=0.7)
ax.set_ylabel('Triangle Area (Å²)', fontsize=12)
ax.set_title('(B) Pocket Openness (Triad Triangle Area)', fontsize=13, fontweight='bold')
ax.legend(fontsize=10, loc='upper right')
ax.grid(alpha=0.2)
ax.set_xlim(0, times[-1])

# --- Panel C: Per-residue RMSF (REPLACES broken RMSD) ---
ax = fig.add_subplot(gs[2, :])
x_pos = np.arange(len(pocket_resids))
colors = ['#d62728' if m else '#1f77b4' for m in triad_mask]
bars = ax.bar(x_pos, rmsf_per_residue, color=colors, alpha=0.8, edgecolor='white', linewidth=0.5)

# Labels
residue_labels = []
for rid in pocket_resids:
    # Map resid to 3-letter name
    sel = u.select_atoms(f"name CA and resid {rid}")
    if len(sel) > 0:
        residue_labels.append(f"{sel.residues[0].resname}{rid}")
    else:
        residue_labels.append(f"R{rid}")

ax.set_xticks(x_pos)
ax.set_xticklabels(residue_labels, rotation=45, ha='right', fontsize=9)
ax.set_ylabel('RMSF (Å)', fontsize=12)
ax.set_title('(C) Per-Residue Cα RMSF — Site 1 Pocket', fontsize=13, fontweight='bold')
ax.axhline(y=np.mean(rmsf_per_residue), color='gray', linestyle=':', alpha=0.5,
           label=f'Mean RMSF ({np.mean(rmsf_per_residue):.1f} Å)')

# Legend for triad highlighting
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='#d62728', alpha=0.8, label='Hydrophobic triad'),
    Patch(facecolor='#1f77b4', alpha=0.8, label='Other pocket residues'),
]
ax.legend(handles=legend_elements, fontsize=9, loc='upper right')
ax.grid(alpha=0.2, axis='y')

# --- Panel D: Area histogram ---
ax = fig.add_subplot(gs[3, 0])
ax.hist(areas, bins=60, color='steelblue', alpha=0.7, edgecolor='navy',
        linewidth=0.3, density=True)
ax.axvline(x=crystal_area, color='red', linestyle='--', linewidth=2, label='Crystal')
ax.axvline(x=np.mean(areas), color='orange', linestyle='--', linewidth=2, label='MD mean')
ax.set_xlabel('Triangle Area (Å²)', fontsize=11)
ax.set_ylabel('Density', fontsize=11)
ax.set_title('(D) Area Distribution', fontsize=12, fontweight='bold')
ax.legend(fontsize=9)

# --- Panel E: Distance violins ---
ax = fig.add_subplot(gs[3, 1])
parts = ax.violinplot([d_6_100, d_6_103, d_100_103], positions=[1, 2, 3],
                       showmeans=True, showmedians=True)
for pc, color in zip(parts['bodies'], ['#1f77b4', '#ff7f0e', '#2ca02c']):
    pc.set_facecolor(color)
    pc.set_alpha(0.5)
ax.scatter([1, 2, 3], [crystal['d_6_100'], crystal['d_6_103'], crystal['d_100_103']],
           color='red', s=80, zorder=5, marker='D', label='Crystal structure')
ax.set_xticks([1, 2, 3])
ax.set_xticklabels(['V6–L100', 'V6–F103', 'L100–F103'], fontsize=10)
ax.set_ylabel('Distance (Å)', fontsize=11)
ax.set_title('(E) Distance Distributions', fontsize=12, fontweight='bold')
ax.legend(fontsize=9, loc='upper right')

# --- Panels F+G: Ligand (if present) ---
if has_ligand:
    # Panel F: Ligand RMSD
    ax = fig.add_subplot(gs[4, 0])
    ax.plot(times, ligand_rmsds, color='crimson', alpha=0.3, linewidth=0.5)
    ax.plot(times[window-1:], rmean(ligand_rmsds, window), color='crimson', linewidth=2)
    ax.set_xlabel('Time (ns)', fontsize=11)
    ax.set_ylabel('RMSD (Å)', fontsize=11)
    ax.set_title('(F) Ligand RMSD', fontsize=12, fontweight='bold')
    ax.grid(alpha=0.2)
    ax.set_xlim(0, times[-1])
    
    # Panel G: Ligand-pocket COM distance
    ax = fig.add_subplot(gs[4, 1])
    ax.plot(times, lig_pocket_dists, color='darkgreen', alpha=0.3, linewidth=0.5)
    ax.plot(times[window-1:], rmean(lig_pocket_dists, window), color='darkgreen', linewidth=2)
    ax.axhline(y=np.mean(lig_pocket_dists), color='green', linestyle='--', alpha=0.5,
               label=f'Mean {np.mean(lig_pocket_dists):.1f} Å')
    ax.set_xlabel('Time (ns)', fontsize=11)
    ax.set_ylabel('Distance (Å)', fontsize=11)
    ax.set_title('(G) Ligand COM – Pocket COM Distance', fontsize=12, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(alpha=0.2)
    ax.set_xlim(0, times[-1])

fig.suptitle(title, fontsize=15, fontweight='bold', y=0.99)

# Save
png_path = os.path.join(WORK, 'pocket_stability.png')
pdf_path = os.path.join(WORK, 'pocket_stability.pdf')
fig.savefig(png_path, dpi=300, bbox_inches='tight')
fig.savefig(pdf_path, bbox_inches='tight')
plt.close()
print(f"  Saved: {png_path}")
print(f"  Saved: {pdf_path}")

# ============================================================
# SAVE UPDATED JSON
# ============================================================
results = {
    "simulation": {
        "length_ns": float(times[-1]),
        "n_frames": len(times),
        "force_field": "AMBER14-all + TIP3P-FB",
        "temperature_K": 300, "pressure_bar": 1.0,
        "ionic_strength_M": 0.15, "timestep_fs": 4.0, "hmr": True
    },
    "triad_distances": {
        "Val6_Leu100": {"mean": float(np.mean(d_6_100)), "std": float(np.std(d_6_100)), "crystal": 11.98},
        "Val6_Phe103": {"mean": float(np.mean(d_6_103)), "std": float(np.std(d_6_103)), "crystal": 8.22},
        "Leu100_Phe103": {"mean": float(np.mean(d_100_103)), "std": float(np.std(d_100_103)), "crystal": 7.61},
    },
    "pocket_area": {
        "mean": float(np.mean(areas)), "std": float(np.std(areas)),
        "crystal": float(crystal_area), "preservation": float(1 - area_pres), "cv": float(cv),
    },
    "rmsf": {
        "residues": [int(r) for r in pocket_resids],
        "values": [float(v) for v in rmsf_per_residue],
        "mean": float(np.mean(rmsf_per_residue)),
        "triad_mean": float(np.mean(rmsf_per_residue[triad_mask])),
    },
    "verdict": verdict,
}

if has_ligand:
    results["ligand_dynamics"] = {
        "rmsd_mean": float(np.mean(ligand_rmsds)),
        "rmsd_std": float(np.std(ligand_rmsds)),
        "pocket_distance_mean": float(np.mean(lig_pocket_dists)),
        "pocket_distance_std": float(np.std(lig_pocket_dists)),
        "bound": bool(np.mean(lig_pocket_dists) < 12.0),
    }

json_path = os.path.join(WORK, 'md_analysis.json')
with open(json_path, 'w') as f:
    json.dump(results, f, indent=2)
print(f"  Saved: {json_path}")

print(f"\nDone!")
