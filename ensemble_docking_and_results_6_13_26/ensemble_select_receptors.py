#!/usr/bin/env python3
"""
Stage 1 of 3: Receptor ensemble generation for TIMP2 Site 1 ensemble docking.
=============================================================================

Purpose
-------
Convert a single 100 ns apo trajectory into a small, reproducible ensemble of
receptor conformations for ensemble (relaxed-complex) docking, so the Site 1
screening prioritization can be tested for robustness to receptor flexibility
rather than dependence on one static crystal conformation.

What it does
------------
1. Loads (md_start.pdb, md_production.dcd) with MDAnalysis.
2. Clusters the trajectory on Site 1 backbone RMSD (Ward hierarchical) and picks
   one medoid frame per cluster. Also adds the single most-open and most-closed
   frames by triad triangle area so the ensemble brackets the sampled range.
3. Aligns each selected protein-only frame to a crystal-frame reference
   (1br9_fixed.pdb) by Kabsch superposition on shared-residue CA atoms, so the
   FIXED Site 1 grid box still lands on the pocket.
4. Re-preps the crystal reference itself with the identical OpenBabel command,
   giving an apples-to-apples crystal baseline receptor.
5. Converts each aligned protein PDB to a rigid receptor PDBQT via OpenBabel.

IMPORTANT alignment assumption
------------------------------
The Site 1 grid box (center 27.99/18.85/15.49) is defined in the coordinate
frame of your receptor.pdbqt, which derives from 1BR9 with no atom movement.
1br9_fixed.pdb (your MD input) shares that frame and your residue numbering, so
aligning to it puts every snapshot in the grid frame. After running, sanity
check by docking ONE frame and confirming the top pose sits near Phe103/the
triad. If your receptor.pdbqt was built from a differently-oriented structure,
pass that structure (converted to PDB) as --ref-pdb instead.

Dependencies: numpy, scipy, MDAnalysis. OpenBabel (obabel) on PATH or via --obabel.

Usage
-----
python ensemble_select_receptors.py \
  --md-pdb  md_start.pdb \
  --md-dcd  md_production.dcd \
  --ref-pdb 1br9_fixed.pdb \
  --out-dir ensemble_receptors \
  --n-clusters 6 \
  --obabel obabel
"""

from __future__ import annotations
import argparse
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import MDAnalysis as mda
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

# ----- Pocket definition (1BR9 numbering); superset of the two definitions in
# your codebase so clustering captures the whole Site 1 backbone -----
SITE1_RESIDS = [3, 4, 5, 6, 8, 11, 99, 100, 101, 102, 103, 104, 105]
TRIAD_RESIDS = [6, 100, 103]
BACKBONE_NAMES = ("N", "CA", "C", "O")

# Crystal reference triad distances (from your analysis) -> reference area
CRYSTAL = dict(d_6_100=11.98, d_6_103=8.22, d_100_103=7.61)


def kabsch(mobile: np.ndarray, reference: np.ndarray):
    """Return R, t with R @ mobile + t ~= reference. Matches your codebase."""
    mob_c = mobile - mobile.mean(axis=0)
    ref_c = reference - reference.mean(axis=0)
    H = mob_c.T @ ref_c
    U, S, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    sign = np.diag([1.0, 1.0, np.sign(d)])
    R = Vt.T @ sign @ U.T
    t = reference.mean(axis=0) - R @ mobile.mean(axis=0)
    return R, t


def triangle_area(a, b, c):
    s = (a + b + c) / 2.0
    return np.sqrt(np.maximum(s * (s - a) * (s - b) * (s - c), 0.0))


def site1_backbone_coords(u: mda.Universe) -> np.ndarray:
    """(n_frames, n_atoms, 3) superposed-later Site 1 backbone coordinates."""
    resid_sel = " ".join(str(r) for r in SITE1_RESIDS)
    name_sel = " ".join(BACKBONE_NAMES)
    bb = u.select_atoms(f"protein and resid {resid_sel} and name {name_sel}")
    if len(bb) == 0:
        sys.exit("ERROR: no Site 1 backbone atoms selected. Check residue numbering.")
    coords = np.empty((len(u.trajectory), len(bb), 3), dtype=np.float64)
    for i, _ in enumerate(u.trajectory):
        coords[i] = bb.positions
    return coords


def triad_areas(u: mda.Universe) -> np.ndarray:
    """Triad triangle area per frame (CA atoms of resid 6, 100, 103)."""
    sel = {r: u.select_atoms(f"protein and resid {r} and name CA") for r in TRIAD_RESIDS}
    for r, ag in sel.items():
        if len(ag) != 1:
            sys.exit(f"ERROR: triad residue {r} CA selection returned {len(ag)} atoms.")
    areas = np.empty(len(u.trajectory), dtype=np.float64)
    for i, _ in enumerate(u.trajectory):
        p6, p100, p103 = (sel[r].positions[0] for r in TRIAD_RESIDS)
        d_6_100 = np.linalg.norm(p6 - p100)
        d_6_103 = np.linalg.norm(p6 - p103)
        d_100_103 = np.linalg.norm(p100 - p103)
        areas[i] = triangle_area(d_6_100, d_6_103, d_100_103)
    return areas


def superpose_all(coords: np.ndarray, ref_idx: int = 0) -> np.ndarray:
    """Superpose every frame onto coords[ref_idx] on the same atom set."""
    ref = coords[ref_idx]
    out = np.empty_like(coords)
    for i in range(coords.shape[0]):
        R, t = kabsch(coords[i], ref)
        out[i] = (R @ coords[i].T).T + t
    return out


def pairwise_rmsd_matrix(superposed: np.ndarray) -> np.ndarray:
    """Condensed RMSD matrix from common-frame superposed coords.

    With all frames superposed to a common reference, RMSD between frames i,j
    equals the Euclidean distance of their flattened coordinate vectors divided
    by sqrt(n_atoms). This is the standard 'align then cluster on Cartesian'
    approach (cf. gmx cluster -fit).
    """
    n_frames, n_atoms, _ = superposed.shape
    flat = superposed.reshape(n_frames, n_atoms * 3)
    # pdist-style, chunked to keep memory modest
    from scipy.spatial.distance import pdist
    d = pdist(flat, metric="euclidean") / np.sqrt(n_atoms)
    return d  # condensed


def cluster_medoids(condensed: np.ndarray, n_frames: int, k: int):
    """Ward clustering -> k clusters; return medoid frame index per cluster."""
    Z = linkage(condensed, method="ward")
    labels = fcluster(Z, t=k, criterion="maxclust")
    full = squareform(condensed)
    medoids = []
    for c in sorted(set(labels)):
        members = np.where(labels == c)[0]
        # medoid = member minimizing mean distance to clustermates
        sub = full[np.ix_(members, members)]
        medoid_local = members[np.argmin(sub.mean(axis=1))]
        medoids.append((int(c), int(medoid_local), len(members)))
    return labels, medoids


def write_aligned_protein_pdb(u: mda.Universe, frame_idx: int,
                              ref_ca_pos: np.ndarray, ref_ca_resids: np.ndarray,
                              out_pdb: Path):
    """Go to frame_idx, align protein to crystal-frame ref, write protein PDB."""
    u.trajectory[frame_idx]
    protein = u.select_atoms("protein")
    prot_ca = u.select_atoms("protein and name CA")
    # match CA by shared resid
    traj_resids = prot_ca.resids
    shared = np.intersect1d(traj_resids, ref_ca_resids)
    if len(shared) < 3:
        sys.exit("ERROR: <3 shared CA residues with reference; cannot align.")
    mob = np.array([prot_ca.positions[np.where(traj_resids == r)[0][0]] for r in shared])
    ref = np.array([ref_ca_pos[np.where(ref_ca_resids == r)[0][0]] for r in shared])
    R, t = kabsch(mob, ref)
    aligned = (R @ protein.positions.T).T + t
    # write a temporary copy with aligned coords
    protein_copy = protein.universe.copy().select_atoms("protein")
    protein_copy.positions = aligned
    protein_copy.write(str(out_pdb))


def obabel_receptor(obabel: str, in_pdb: Path, out_pdbqt: Path, log_dir: Path):
    """Convert protein PDB -> rigid receptor PDBQT, matching your prep route."""
    log = log_dir / (out_pdbqt.stem + ".obabel.log")
    cmd = [obabel, str(in_pdb), "-O", str(out_pdbqt),
           "-xr", "--partialcharge", "gasteiger"]
    with open(log, "w", encoding="utf-8") as lf:
        lf.write(">> " + " ".join(cmd) + "\n")
        proc = subprocess.run(cmd, stdout=lf, stderr=lf)
    if proc.returncode != 0 or not out_pdbqt.exists():
        print(f"  [WARN] OpenBabel failed for {in_pdb.name} (see {log.name})")
        return False
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--md-pdb", required=True, type=Path, help="md_start.pdb (solvated topology)")
    ap.add_argument("--md-dcd", required=True, type=Path, help="md_production.dcd")
    ap.add_argument("--ref-pdb", required=True, type=Path,
                    help="Crystal-frame reference (1br9_fixed.pdb). Defines the grid frame.")
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--n-clusters", type=int, default=6)
    ap.add_argument("--obabel", default="obabel", help="OpenBabel executable")
    ap.add_argument("--no-pdbqt", action="store_true",
                    help="Stop after writing aligned PDBs (skip OpenBabel).")
    ap.add_argument("--save-ps", type=float, default=50.0,
                    help="Frame save interval in ps (for reporting trajectory length). Default 50.")
    ap.add_argument("--min-frames", type=int, default=0,
                    help="Abort if fewer than this many frames load (guards against a "
                         "silently truncated / wrong DCD). E.g. 1900 for a 100 ns / 2,000-frame run.")
    args = ap.parse_args()

    out_dir = args.out_dir.resolve()
    pdb_dir = out_dir / "aligned_pdb"
    pdbqt_dir = out_dir / "receptor_pdbqt"
    log_dir = out_dir / "logs"
    for d in (pdb_dir, pdbqt_dir, log_dir):
        d.mkdir(parents=True, exist_ok=True)

    print("[1] Loading trajectory ...")
    u = mda.Universe(str(args.md_pdb), str(args.md_dcd))
    n_frames = len(u.trajectory)
    est_ns = n_frames * args.save_ps / 1000.0
    print(f"    Trajectory: {args.md_dcd}")
    print(f"    Frames loaded: {n_frames}  (~{est_ns:.0f} ns at {args.save_ps:g} ps/frame)")
    if args.min_frames and n_frames < args.min_frames:
        sys.exit(f"ABORT: only {n_frames} frames loaded but --min-frames={args.min_frames}. "
                 f"This usually means a truncated or wrong DCD was supplied (e.g. a 50 ns copy "
                 f"in the working dir instead of the full 100 ns trajectory). Check the path and "
                 f"cross-check n_frames against the apo md_analysis.json.")

    print("[2] Computing triad areas and Site 1 backbone coordinates ...")
    areas = triad_areas(u)
    coords = site1_backbone_coords(u)
    crystal_area = float(triangle_area(CRYSTAL["d_6_100"], CRYSTAL["d_6_103"], CRYSTAL["d_100_103"]))
    cv = 100.0 * areas.std() / areas.mean()
    print(f"    Triad area: mean {areas.mean():.1f} +/- {areas.std():.1f} A^2 (CV {cv:.0f}%) "
          f"(crystal {crystal_area:.1f}); range {areas.min():.1f}-{areas.max():.1f}")
    print(f"    [sanity] full 100 ns apo run should give CV ~10%; CV ~18% indicates only "
          f"the first ~50 ns was read.")

    print("[3] Clustering on Site 1 backbone RMSD ...")
    superposed = superpose_all(coords, ref_idx=0)
    condensed = pairwise_rmsd_matrix(superposed)
    labels, medoids = cluster_medoids(condensed, n_frames, args.n_clusters)

    # Selected frames: cluster medoids + extremes by openness
    selected = {}  # frame_idx -> label string
    for c, medoid, size in medoids:
        selected[medoid] = f"cluster{c}_n{size}"
    open_idx = int(np.argmax(areas))
    closed_idx = int(np.argmin(areas))
    selected.setdefault(open_idx, "extreme_open")
    selected.setdefault(closed_idx, "extreme_closed")

    print(f"    {len(medoids)} medoids + open/closed extremes "
          f"-> {len(selected)} unique receptor conformations")
    for fi in sorted(selected):
        print(f"      frame {fi:5d}  area {areas[fi]:5.1f} A^2  [{selected[fi]}]")

    # Reference for alignment
    ref_u = mda.Universe(str(args.ref_pdb))
    ref_ca = ref_u.select_atoms("protein and name CA")
    ref_ca_pos = ref_ca.positions.copy()
    ref_ca_resids = ref_ca.resids.copy()

    print("[4] Writing aligned protein PDBs ...")
    manifest = out_dir / "ensemble_manifest.csv"
    rows = ["receptor_id,frame_index,triad_area_A2,kind"]

    # crystal baseline (re-prepped with identical command)
    crystal_pdb = pdb_dir / "receptor_crystal.pdb"
    ref_u.select_atoms("protein").write(str(crystal_pdb))
    rows.append(f"receptor_crystal,-1,{crystal_area:.2f},crystal_baseline")

    for fi in sorted(selected):
        rid = f"receptor_frame{fi:05d}"
        out_pdb = pdb_dir / f"{rid}.pdb"
        write_aligned_protein_pdb(u, fi, ref_ca_pos, ref_ca_resids, out_pdb)
        rows.append(f"{rid},{fi},{areas[fi]:.2f},{selected[fi]}")

    manifest.write_text("\n".join(rows) + "\n", encoding="utf-8")
    print(f"    Manifest: {manifest}")

    if args.no_pdbqt:
        print("[5] --no-pdbqt set; skipping OpenBabel conversion.")
        print(f"\nDone. Aligned PDBs in {pdb_dir}")
        return 0

    print("[5] Converting to rigid receptor PDBQT (OpenBabel) ...")
    ok = 0
    for pdb in sorted(pdb_dir.glob("*.pdb")):
        out_pdbqt = pdbqt_dir / (pdb.stem + ".pdbqt")
        if obabel_receptor(args.obabel, pdb, out_pdbqt, log_dir):
            ok += 1
    print(f"    Converted {ok} receptor PDBQTs -> {pdbqt_dir}")
    print("\nSanity check before docking the full set:")
    print("  dock ONE frame and confirm the top pose sits near Phe103/the triad.")
    print("  If it lands elsewhere, your --ref-pdb is not in the grid frame.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
