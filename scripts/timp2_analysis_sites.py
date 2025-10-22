#!/usr/bin/env python3
"""
TIMP2 Site-Aware ΔScore Analysis
--------------------------------
For merged files that contain both Site1 and Site2 data (columns: site, zinc_id, docking_score, smiles).

Computes site selectivity (ΔScore = Site2 − Site1) for ligands present at both binding sites,
adds a delta_score column, and exports:
  • full annotated CSV  -> *_analysis_full.csv
  • summary CSV         -> site_selectivity_summary.csv
"""

import pandas as pd
import numpy as np
from pathlib import Path
import argparse

def analyze_sites(infile: Path, outdir: Path):
    df = pd.read_csv(infile)
    print(f"Loaded {len(df)} rows from {infile}")

    required = {"id", "site", "docking_score"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    # Normalize site names (1/2)
    df["site"] = df["site"].astype(str).str.strip().str.replace("site", "", regex=False)

    # Base ligand identifier
    df["lig_base"] = df["id"].str.extract(r"(lig_\d{6,8})", expand=False)

    # Pivot to pair sites
    pivot = (
        df.pivot_table(index="lig_base", columns="site", values="docking_score", aggfunc="min")
        .dropna(subset=["1", "2"], how="any")
    )

    pivot["delta_score"] = pivot["2"] - pivot["1"]
    pivot = pivot.reset_index()

    # Merge back
    df = df.merge(pivot[["lig_base", "delta_score"]], on="lig_base", how="left")

    outdir.mkdir(parents=True, exist_ok=True)
    out_full = outdir / f"{infile.stem}_analysis_full.csv"
    out_summary = outdir / "site_selectivity_summary.csv"

    df.to_csv(out_full, index=False)
    pivot.sort_values("delta_score").to_csv(out_summary, index=False)

    print(f"Saved ΔScore analysis: {out_full}")
    print(f"Summary written to: {out_summary}")
    return df, pivot


def main():
    ap = argparse.ArgumentParser(description="Compute site selectivity ΔScore for TIMP2 merged files")
    ap.add_argument("infile", type=Path, help="Input merged CSV file with site data")
    ap.add_argument("-o", "--outdir", type=Path, required=True, help="Output directory")
    args = ap.parse_args()

    analyze_sites(args.infile, args.outdir)

if __name__ == "__main__":
    main()
