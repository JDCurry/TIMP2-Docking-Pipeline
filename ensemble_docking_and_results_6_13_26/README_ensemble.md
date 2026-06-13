# TIMP2 Site 1 Ensemble Docking

Three scripts that convert your existing 100 ns apo trajectory + Vina pipeline
into an ensemble-docking robustness check, addressing the "single static
crystal structure" limitation (item 4) without new wet-lab work or new MD.

**What this experiment is for:** testing whether the *screening prioritization*
of Site 1 is robust to receptor conformation. It is NOT re-validation of the
lead — your 100 ns ligand-bound MD already does that better. The readout is
pose/anchor convergence, NOT docking score (your own data shows score is a poor
residence predictor: ZINC000870248439 has the best Vina score yet dissociates).



## Prerequisites

- Files you already have: `md_start.pdb`, `md_production.dcd`, `1br9_fixed.pdb`,
  `receptor.pdbqt`, and Vina + OpenBabel.
- A `focused7.smi` with the seven MD compounds, one `SMILES ZINCID` per line.
  The lead is: `NC(=O)c1ccc(C(=O)N2CCC[C@@]3(CCSC3)C2)cc1 ZINC000583127796`
  (pull the other six SMILES from your Tier 1A `File_S3` CSV).
- Python env with: numpy, scipy, MDAnalysis (Stage 1); numpy, matplotlib
  (Stage 3). Stage 2 needs only the Vina/OpenBabel executables.

---

## Stage 1 — Build the receptor ensemble

```bash
python ensemble_select_receptors.py \
  --md-pdb  md_start.pdb \
  --md-dcd  md_production.dcd \
  --ref-pdb 1br9_fixed.pdb \
  --out-dir ensemble_receptors \
  --n-clusters 6 \
  --obabel obabel
```

Produces ~8 receptor PDBQTs in `ensemble_receptors/receptor_pdbqt/`
(6 cluster medoids + most-open + most-closed frames) plus a re-prepped
`receptor_crystal.pdbqt` baseline, all aligned to the crystal grid frame.

**CRITICAL alignment assumption:** the Site 1 grid lives in your
`receptor.pdbqt` frame; `1br9_fixed.pdb` shares that frame and numbering, so
aligning to it lands the fixed grid on the pocket. **Sanity check before
proceeding:** dock ONE frame and confirm the top pose sits near Phe103/the
triad. If it lands elsewhere, your `--ref-pdb` is not in the grid frame.

Also note: this uses ONLY apo frames. Ligand-bound frames are deliberately
excluded — docking a ligand back into the conformation it induced is circular.
(Cross-docking into ligand-induced frames is a separate optional analysis.)

---

## Stage 2 — Dock the focused set into every conformation

```bash
python dock_ensemble.py \
  --receptor-dir ensemble_receptors/receptor_pdbqt \
  --ligand-smi   focused7.smi \
  --out-dir      ensemble_docking \
  --vina vina --obabel obabel \
  --exhaustiveness 32 --num-modes 9 --seeds 1 2 3 4 5 \
  --processes 16
```

Same Site 1 grid as your screen, but elevated sampling (exhaustiveness 32, 9
modes, 5 seeds) because the set is small. The crystal baseline is docked at the
identical settings, so crystal-vs-ensemble is apples-to-apples. ~8 receptors x
7 ligands x 5 seeds ≈ 280 dockings; a few hours on your 36-core box.

Ligands are prepped fresh from SMILES (neutral gen3d), not from your saved
docked poses, so the search isn't biased toward the crystal-docked pose.

---

## Stage 3 — Convergence analysis

```bash
python analyze_ensemble.py \
  --scores-csv   ensemble_docking/ensemble_scores.csv \
  --receptor-dir ensemble_receptors/receptor_pdbqt \
  --out-dir      ensemble_analysis \
  --lead ZINC000583127796 \
  --accept-threshold -6.0 --contact-cutoff 4.0
```

Produces:
- `ensemble_analysis.csv` / `.json` — per-ligand acceptance fraction, pose
  convergence (mean pairwise RMSD), and Phe103/Leu100/Val6 recurrence.
- `ensemble_analysis.png` — three-panel figure.
- `ensemble_results.tex` — a drop-in LaTeX snippet with the numbers filled in.




