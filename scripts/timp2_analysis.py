#!/usr/bin/env python3
"""
TIMP2 Analysis Script — Site-Aware ΔScore Edition

This version computes site selectivity (ΔScore = Site2 − Site1) for ligands present at both binding sites,
adds the delta_score column to the output, and exports a site_selectivity_summary.csv for analysis.

Typical workflow:
  1. Run zinc_recovery.py to map ligand IDs to ZINC IDs
  2. Run this script on the ZINC-mapped CSV per chunk
  3. Follow with tier1_filter_admetai.py for ADMET/Tier assignment

Usage:
  python timp2_analysis.py input_with_zinc.csv -o output_directory
"""

import pandas as pd
import numpy as np
import argparse
from pathlib import Path
import re

# ------------------------------------------------------
# Helpers
# ------------------------------------------------------

def base_ligand_id(lig):
    if isinstance(lig, str):
        m = re.search(r'(lig_\d{8})', lig)
        return m.group(1) if m else lig
    return lig

# ------------------------------------------------------
# Core Analysis Logic
# ------------------------------------------------------

def analyze_file(infile: Path, outdir: Path):
    df = pd.read_csv(infile)
    print(f"Loaded {len(df)} entries from {infile}")

    # Basic filtering and descriptors (placeholder — assume precomputed)
    if 'docking_score' not in df.columns:
        raise ValueError('Input file must include docking_score column')

    # Compute base ligand ID and clean site labels
    if 'site' not in df.columns:
        # Try to infer site from ID pattern
        df['site'] = df['id'].apply(lambda x: '1' if 'site1' in str(x).lower() else ('2' if 'site2' in str(x).lower() else 'UNK'))

    df['lig_base'] = df['id'].apply(base_ligand_id)

    # ------------------------------------------------------
    # ΔScore (Site2 − Site1)
    # ------------------------------------------------------
    paired = (
        df.pivot_table(
            index='lig_base',
            columns='site',
            values='docking_score',
            aggfunc='min'
        )
        .dropna(subset=['1', '2'], how='any')
    )

    paired['delta_score'] = paired['2'] - paired['1']
    paired = paired.reset_index()

    df = df.merge(paired[['lig_base', 'delta_score']], on='lig_base', how='left')

    # ------------------------------------------------------
    # Exports
    # ------------------------------------------------------
    outdir.mkdir(parents=True, exist_ok=True)

    out_main = outdir / f"{infile.stem}_analysis_full.csv"
    out_summary = outdir / "site_selectivity_summary.csv"

    df.to_csv(out_main, index=False)
    paired.sort_values('delta_score').to_csv(out_summary, index=False)

    print(f"Saved analysis: {out_main}")
    print(f"ΔScore summary written: {out_summary}")

    return df, paired

# ------------------------------------------------------
# CLI
# ------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description='TIMP2 site-aware analysis with ΔScore computation')
    ap.add_argument('infile', type=Path, help='Input CSV with ZINC IDs and docking scores')
    ap.add_argument('-o', '--outdir', type=Path, required=True, help='Output directory for analysis results')
    args = ap.parse_args()

    analyze_file(args.infile, args.outdir)

if __name__ == '__main__':
    main()