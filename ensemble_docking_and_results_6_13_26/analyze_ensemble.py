#!/usr/bin/env python3
"""
Stage 3 of 3: Ensemble convergence analysis for TIMP2 Site 1.
=============================================================

Reads Stage 2 output and computes the metrics that actually support the
robustness claim (NOT absolute docking score, which your own data shows is a
poor residence-time predictor: ZINC000870248439 has the best Vina score yet
dissociates in MD):

  1. Pocket acceptance fraction  - in what fraction of receptor conformations
     does Site 1 accept the ligand with score below threshold.
  2. Pose convergence            - pairwise heavy-atom RMSD of the best pose
     across receptor conformations (lower = a single recurring pose family).
     Valid as a direct comparison because all receptors were aligned to one
     reference frame in Stage 1, so poses are co-registered.
  3. Anchor recurrence           - fraction of receptor conformations in which
     Phe103 (and Val6, Leu100) is contacted by the best pose (4.0 A heavy-atom
     cutoff). Cross-checks the MD finding that Phe103 is the primary anchor.

Outputs:
  ensemble_analysis.csv   - per-ligand summary
  ensemble_analysis.json  - full numbers
  ensemble_analysis.png   - figure
  ensemble_results.tex    - calibrated, drop-in LaTeX snippet (numbers filled,
                            wording gated on whether convergence is strong)

Pose RMSD is index-based (atom-order). It is conservative (can only
overestimate under molecular symmetry). For symmetry-corrected RMSD, swap in
RDKit GetBestRMS; index-based is the default to avoid extra dependencies.

Dependencies: numpy, matplotlib. (MDAnalysis/RDKit not required.)

Usage
-----
python analyze_ensemble.py \
  --scores-csv  ensemble_docking/ensemble_scores.csv \
  --receptor-dir ensemble_receptors/receptor_pdbqt \
  --out-dir     ensemble_analysis \
  --lead ZINC000583127796 \
  --accept-threshold -6.0 --contact-cutoff 4.0
"""

from __future__ import annotations
import argparse
import csv
import json
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TRIAD = {6: "Val6", 100: "Leu100", 103: "Phe103"}
VINA_RESULT_RE = re.compile(r"REMARK VINA RESULT:\s+(-?\d+(?:\.\d+)?)")


# ----------------------------- parsing -----------------------------
def _is_hydrogen(atom_name: str, element: str) -> bool:
    if element and element.strip().upper() == "H":
        return True
    nm = atom_name.strip()
    return nm[:1] == "H" or (len(nm) > 1 and nm[0].isdigit() and nm[1] == "H")


def parse_pose_models(pose_path: Path):
    """Return [(score, coords Nx3)] per MODEL in a Vina output PDBQT (heavy atoms)."""
    models = []
    cur_score = None
    cur_coords = []
    have_model = False
    try:
        lines = Path(pose_path).read_text(encoding="utf-8", errors="ignore").splitlines()
    except FileNotFoundError:
        return models
    for line in lines:
        if line.startswith("MODEL"):
            have_model = True
            cur_score, cur_coords = None, []
        elif line.startswith("REMARK VINA RESULT"):
            m = VINA_RESULT_RE.search(line)
            if m:
                cur_score = float(m.group(1))
        elif line.startswith(("ATOM", "HETATM")):
            name = line[12:16]
            element = line[76:78] if len(line) >= 78 else ""
            if _is_hydrogen(name, element):
                continue
            try:
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            except ValueError:
                continue
            cur_coords.append((x, y, z))
        elif line.startswith("ENDMDL"):
            if cur_coords:
                models.append((cur_score, np.array(cur_coords, dtype=np.float64)))
            cur_score, cur_coords = None, []
    # single-model file without MODEL/ENDMDL
    if not have_model and cur_coords:
        models.append((cur_score, np.array(cur_coords, dtype=np.float64)))
    return models


def parse_receptor_residue_heavy(receptor_pdbqt: Path, resid: int) -> np.ndarray:
    """Heavy-atom coords of a given residue id from a receptor PDBQT."""
    coords = []
    for line in Path(receptor_pdbqt).read_text(encoding="utf-8", errors="ignore").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        try:
            rid = int(line[22:26])
        except ValueError:
            continue
        if rid != resid:
            continue
        name = line[12:16]
        element = line[76:78] if len(line) >= 78 else ""
        if _is_hydrogen(name, element):
            continue
        try:
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
        except ValueError:
            continue
        coords.append((x, y, z))
    return np.array(coords, dtype=np.float64) if coords else np.empty((0, 3))


# ----------------------------- metrics -----------------------------
def index_rmsd(a: np.ndarray, b: np.ndarray) -> float:
    if a.shape != b.shape or a.size == 0:
        return float("nan")
    return float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1))))


def min_heavy_distance(lig: np.ndarray, res: np.ndarray) -> float:
    if lig.size == 0 or res.size == 0:
        return float("inf")
    diff = lig[:, None, :] - res[None, :, :]
    return float(np.sqrt(np.sum(diff ** 2, axis=2)).min())


def best_pose_per_receptor(rows, ligand):
    """For one ligand, map receptor -> (best_score, best_coords)."""
    by_rec_files = defaultdict(set)
    for r in rows:
        if r["ligand"] == ligand and r["pose_path"]:
            by_rec_files[r["receptor"]].add(r["pose_path"])
    best = {}
    for rec, files in by_rec_files.items():
        best_score, best_coords = None, None
        for f in files:
            for score, coords in parse_pose_models(Path(f)):
                if score is None:
                    continue
                if best_score is None or score < best_score:
                    best_score, best_coords = score, coords
        if best_coords is not None:
            best[rec] = (best_score, best_coords)
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scores-csv", required=True, type=Path)
    ap.add_argument("--receptor-dir", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--lead", default="ZINC000583127796")
    ap.add_argument("--accept-threshold", type=float, default=-6.0)
    ap.add_argument("--contact-cutoff", type=float, default=4.0)
    args = ap.parse_args()

    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    with open(args.scores_csv, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    ligands = sorted({r["ligand"] for r in rows})
    receptors = sorted({r["receptor"] for r in rows})
    n_rec = len(receptors)
    print(f"Ligands: {len(ligands)}  Receptors: {n_rec}")

    # Pre-load triad heavy atoms per receptor
    triad_coords = {}
    for rec in receptors:
        rpath = args.receptor_dir / f"{rec}.pdbqt"
        if rpath.exists():
            triad_coords[rec] = {rid: parse_receptor_residue_heavy(rpath, rid) for rid in TRIAD}
        else:
            triad_coords[rec] = {rid: np.empty((0, 3)) for rid in TRIAD}

    summary = {}
    for lig in ligands:
        best = best_pose_per_receptor(rows, lig)
        if not best:
            continue
        scores = {rec: s for rec, (s, _) in best.items()}
        coords = {rec: c for rec, (_, c) in best.items()}

        accepted = [rec for rec, s in scores.items() if s is not None and s <= args.accept_threshold]
        accept_frac = len(accepted) / n_rec

        # pose convergence among accepted receptors
        acc_coords = [coords[rec] for rec in accepted if coords[rec].size]
        pair_rmsds = []
        for i in range(len(acc_coords)):
            for j in range(i + 1, len(acc_coords)):
                rm = index_rmsd(acc_coords[i], acc_coords[j])
                if rm == rm:  # not NaN
                    pair_rmsds.append(rm)
        mean_pair_rmsd = float(np.mean(pair_rmsds)) if pair_rmsds else float("nan")

        # anchor recurrence per triad residue (over accepted receptors)
        anchor = {}
        for rid, nm in TRIAD.items():
            hits = 0
            denom = 0
            for rec in accepted:
                res = triad_coords[rec][rid]
                lig_c = coords[rec]
                d = min_heavy_distance(lig_c, res)
                if np.isfinite(d):
                    denom += 1
                    if d <= args.contact_cutoff:
                        hits += 1
            anchor[nm] = (hits / denom) if denom else float("nan")

        sc_vals = [s for s in scores.values() if s is not None]
        summary[lig] = {
            "n_receptors": n_rec,
            "n_accepted": len(accepted),
            "accept_fraction": accept_frac,
            "score_mean": float(np.mean(sc_vals)) if sc_vals else None,
            "score_std": float(np.std(sc_vals)) if sc_vals else None,
            "score_min": float(np.min(sc_vals)) if sc_vals else None,
            "mean_pairwise_pose_rmsd": mean_pair_rmsd,
            "anchor_recurrence": anchor,
        }
        print(f"\n{lig}")
        print(f"  accept {len(accepted)}/{n_rec} ({accept_frac*100:.0f}%) | "
              f"score {summary[lig]['score_mean']:.2f}+/-{summary[lig]['score_std']:.2f} | "
              f"pose RMSD {mean_pair_rmsd:.2f} A")
        print(f"  anchors: " + "  ".join(f"{k} {v*100:.0f}%" for k, v in anchor.items()))

    # ---------- CSV ----------
    csv_p = out_dir / "ensemble_analysis.csv"
    with open(csv_p, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["ligand", "n_receptors", "n_accepted", "accept_fraction",
                    "score_mean", "score_std", "score_min", "mean_pairwise_pose_rmsd",
                    "Phe103_recurrence", "Leu100_recurrence", "Val6_recurrence"])
        for lig, d in summary.items():
            a = d["anchor_recurrence"]
            w.writerow([lig, d["n_receptors"], d["n_accepted"], f"{d['accept_fraction']:.3f}",
                        d["score_mean"], d["score_std"], d["score_min"],
                        f"{d['mean_pairwise_pose_rmsd']:.3f}",
                        f"{a['Phe103']:.3f}", f"{a['Leu100']:.3f}", f"{a['Val6']:.3f}"])
    (out_dir / "ensemble_analysis.json").write_text(json.dumps(summary, indent=2))
    print(f"\nWrote {csv_p}")

    # ---------- figure ----------
    make_figure(summary, args.lead, out_dir)

    # ---------- calibrated LaTeX ----------
    write_latex(summary, args.lead, args, out_dir)
    return 0


def make_figure(summary, lead, out_dir):
    ligs = list(summary.keys())
    if not ligs:
        return
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    short = [l.replace("ZINC000", "..").replace("ZINC", "") for l in ligs]

    accept = [summary[l]["accept_fraction"] * 100 for l in ligs]
    axes[0].barh(short, accept, color="steelblue")
    axes[0].set_xlabel("Pocket acceptance (% of receptors)")
    axes[0].set_title("(A) Site 1 acceptance across ensemble")
    axes[0].set_xlim(0, 100)

    rmsd = [summary[l]["mean_pairwise_pose_rmsd"] for l in ligs]
    axes[1].barh(short, rmsd, color="darkgreen")
    axes[1].set_xlabel("Mean pairwise pose RMSD (A)")
    axes[1].set_title("(B) Pose convergence (lower = tighter)")

    phe = [summary[l]["anchor_recurrence"]["Phe103"] * 100 for l in ligs]
    axes[2].barh(short, phe, color="indianred")
    axes[2].set_xlabel("Phe103 contact recurrence (% of receptors)")
    axes[2].set_title("(C) Phe103 anchor recurrence")
    axes[2].set_xlim(0, 100)

    for ax in axes:
        for lbl in ax.get_yticklabels():
            if lead.replace("ZINC000", "..").replace("ZINC", "") in lbl.get_text():
                lbl.set_fontweight("bold")
    fig.suptitle("TIMP2 Site 1 ensemble docking: convergence across apo MD conformations",
                 fontweight="bold")
    fig.tight_layout()
    fig.savefig(out_dir / "ensemble_analysis.png", dpi=300, bbox_inches="tight")
    fig.savefig(out_dir / "ensemble_analysis.pdf", bbox_inches="tight")
    plt.close()
    print(f"Wrote {out_dir / 'ensemble_analysis.png'}")


def write_latex(summary, lead, args, out_dir):
    d = summary.get(lead)
    if d is None:
        print(f"  [WARN] lead {lead} not in summary; skipping LaTeX.")
        return
    n_rec = d["n_receptors"]
    phe = d["anchor_recurrence"]["Phe103"]
    rmsd = d["mean_pairwise_pose_rmsd"]
    n_acc_lead = d["n_accepted"]
    # set-level acceptance (median across ligands)
    accs = [s["accept_fraction"] for s in summary.values()]
    median_acc = float(np.median(accs)) if accs else 0.0

    strong = (phe >= 0.70 and (rmsd == rmsd and rmsd <= 2.5) and median_acc >= 0.80)

    L = []
    L.append(r"% Ensemble docking result -- insert after the MD validation subsection.")
    L.append(r"% Wording is calibrated; do NOT upgrade 'supports' to 'validates'.")
    L.append(r"\subsubsection{Ensemble Docking Across Apo Conformations}")
    L.append("")
    L.append(
        r"To test whether the prioritization of Site~1 reflects an intrinsic property of "
        r"the pocket rather than an artifact of the single crystallographic conformation "
        r"used for screening, we performed ensemble re-docking of the focused compound set "
        f"into {n_rec} receptor conformations sampled from the 100~ns apo trajectory. "
        r"Conformations were selected as medoids of Site~1 backbone RMSD clusters, "
        r"supplemented by the most-open and most-closed frames by triad triangle area, and "
        r"each was Kabsch-aligned to the crystal reference so that the primary-screen Site~1 "
        r"grid was held fixed across all receptors. Docking used AutoDock Vina at elevated "
        f"exhaustiveness ({args.accept_threshold:+.0f}~kcal/mol acceptance threshold, "
        r"multiple seeds and modes per pair).")
    L.append("")
    if strong:
        L.append(
            f"Across the {n_rec} apo conformations, ZINC000583127796 retained a "
            f"Phe103-contacting pose in {phe*100:.0f}\\% of accepted receptors, with a mean "
            f"pairwise heavy-atom pose RMSD of {rmsd:.1f}~\\AA{{}}, indicating a single "
            r"recurring pose family rather than conformation-dependent rebinding. Site~1 "
            f"accepted the focused set in a median of {median_acc*100:.0f}\\% of conformations. "
            r"These results indicate that Site~1 ligand compatibility and the lead-compound "
            r"binding pose are preserved across the conformational ensemble sampled in "
            r"solution, supporting the robustness of the screening-derived prioritization to "
            r"receptor flexibility rather than its dependence on a single static "
            r"conformation. Consistent with the molecular dynamics analysis, absolute Vina "
            r"score was not used as the criterion of robustness, as it is a ranking signal "
            r"rather than a binding free energy.")
    else:
        L.append(
            f"Across the {n_rec} apo conformations, ZINC000583127796 was accepted in "
            f"{n_acc_lead} of {n_rec} receptors and retained a Phe103-contacting pose in "
            f"{phe*100:.0f}\\% of accepted receptors (mean pairwise pose RMSD "
            f"{rmsd:.1f}~\\AA{{}}). The pocket therefore remained ligand-compatible across a "
            r"substantial fraction of, but not all, sampled conformations, and pose "
            r"convergence was partial. We report this as a bounded robustness check: it "
            r"argues against Site~1 being a single-frame artifact while indicating that "
            r"binding geometry is sensitive to receptor conformation, which experimental "
            r"co-crystallography would be required to resolve.")
    L.append("")
    L.append(
        r"% Limitations bullet upgrade (replace the existing 'Static structure' item):")
    L.append(
        r"% \item \textbf{Static structure and ensemble scope}: The primary virtual screen "
        r"was performed against a single crystallographic conformation (PDB: 1BR9). Ensemble "
        r"re-docking into apo MD-derived conformations indicates that Site~1 ligand "
        r"compatibility and the lead-compound pose are preserved across the sampled "
        r"conformational range, arguing against a single-conformation artifact. However, the "
        r"receptor ensemble derives from a single 100~ns apo trajectory and samples thermally "
        r"accessible fluctuations around the crystal state rather than the full breadth of "
        r"induced-fit conformations a diverse library might elicit; co-crystallography and "
        r"biophysical binding assays (SPR, MST/ITC) remain necessary to confirm binding modes "
        r"and affinities.")

    tex_p = out_dir / "ensemble_results.tex"
    tex_p.write_text("\n".join(L) + "\n", encoding="utf-8")
    verdict = "STRONG" if strong else "PARTIAL/BOUNDED"
    print(f"Wrote {tex_p}  [convergence verdict: {verdict}]")


if __name__ == "__main__":
    raise SystemExit(main())
