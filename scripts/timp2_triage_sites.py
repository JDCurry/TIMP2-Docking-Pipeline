#!/usr/bin/env python3
"""
TIMP2 Site-Aware Triage (v2.1)

Purpose
-------
Preserve and analyze ligand docking results from multiple binding sites (e.g., Site 1, Site 2)
without collapsing entries by SMILES or clustering away duplicates.

Key Features
-------------
1. Reads scores + chunk SMILES (or direct hits CSV) just like timp2_triage.py.
2. Automatically detects or assigns 'site' labels (supports both 'site1'/'site2' and '_A'/'_B').
3. Computes RDKit descriptors, CNS MPO, and PAINS/Brenk filters.
4. Assigns each molecule to CNS/Peripheral lanes.
5. **Preserves duplicates** — no collapsing by SMILES or cluster diversity filtering.
6. Exports separate per-site and combined triage results for comparative analysis.

Usage Example
-------------
python timp2_triage_sites.py \
  --hits-csv strict_0109_hits.csv \
  --out triage_sites_out

Outputs
-------
out/
  - all_sites_full.csv   (complete dataset, all sites)
  - site1_hits.csv       (subset for Site 1)
  - site2_hits.csv       (subset for Site 2)
  - summary.md
"""

import argparse
from pathlib import Path
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import FilterCatalog

# ------------------------------------------------------
# Utilities
# ------------------------------------------------------

def extract_site_from_id(ligand_id: str) -> str:
    if not isinstance(ligand_id, str):
        return 'UNK'
    lid = ligand_id.lower()
    if 'site1' in lid or 'site_1' in lid:
        return '1'
    if 'site2' in lid or 'site_2' in lid:
        return '2'
    if lid.endswith('a'):
        return 'A'
    if lid.endswith('b'):
        return 'B'
    return 'UNK'

def build_name_to_smiles(chunks_dir: Path, pattern: str):
    mapping = {}
    for smi in sorted(chunks_dir.glob(pattern)):
        with open(smi, 'r', encoding='utf-8', errors='ignore') as f:
            for i, line in enumerate(f):
                line = line.strip()
                if not line:
                    continue
                parts = line.split()
                smi_txt = parts[0]
                name = parts[1] if len(parts) > 1 else f"{smi.stem}_{i:07d}"
                mapping[name] = smi_txt
    return mapping

# ------------------------------------------------------
# Descriptor + ADMET helpers
# ------------------------------------------------------

def compute_rdkit_props(smiles: str):
    m = Chem.MolFromSmiles(smiles)
    if not m:
        return None
    return {
        'mol': m,
        'mw': Descriptors.MolWt(m),
        'tpsa': rdMolDescriptors.CalcTPSA(m),
        'hbd': Lipinski.NumHDonors(m),
        'hba': Lipinski.NumHAcceptors(m),
        'rotb': Lipinski.NumRotatableBonds(m),
        'clogp': Crippen.MolLogP(m),
        'fsp3': rdMolDescriptors.CalcFractionCSP3(m),
        'scaffold': MurckoScaffold.MurckoScaffoldSmiles(mol=m)
    }

def pains_brenk_flags(m: Chem.Mol) -> str:
    params = FilterCatalog.FilterCatalogParams()
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)
    catalog = FilterCatalog.FilterCatalog(params)
    alerts = [entry.GetDescription() for entry in catalog.GetMatches(m)]
    return '; '.join(alerts)

def cns_mpo(mw, tpsa, hbd, clogp):
    s = 0
    s += 1 - min(max((mw - 360) / 140, 0), 1)
    s += 1 - min(max((tpsa - 60) / 60, 0), 1)
    s += 1 - min(max((hbd - 1) / 2, 0), 1)
    s += 1 - abs(clogp - 3) / 3
    return round(6 * s / 4, 2)

LANE_A_CNS = {'mw': (None, 450), 'tpsa': (None, 70), 'hbd': (None, 1), 'hba': (None, 7), 'rotb': (None, 8), 'clogp': (2.0, 3.5)}
LANE_B_PERIPH = {'mw': (None, 550), 'tpsa': (None, 100), 'hbd': (None, 2), 'hba': (None, 10), 'rotb': (None, 10), 'clogp': (1.0, 4.5)}

def passes_lane(props, lane):
    def in_range(val, lo, hi):
        if lo is not None and val < lo: return False
        if hi is not None and val > hi: return False
        return True
    return all(in_range(props[k], *lane[k]) for k in lane)

# ------------------------------------------------------
# Core Triage Logic
# ------------------------------------------------------

def triage_sites(df_hits: pd.DataFrame, outdir: Path):
    records = []
    for _, row in df_hits.iterrows():
        r = compute_rdkit_props(row['smiles'])
        if not r: continue
        m = r.pop('mol')
        alerts = pains_brenk_flags(m)
        mpo = cns_mpo(r['mw'], r['tpsa'], r['hbd'], r['clogp'])
        lane = 'CNS' if passes_lane(r, LANE_A_CNS) else ('PERIPH' if passes_lane(r, LANE_B_PERIPH) else 'REJECT')
        records.append({**row, **r, 'alerts': alerts, 'mpo_cns': mpo, 'lane': lane})

    df = pd.DataFrame(records)
    outdir.mkdir(parents=True, exist_ok=True)
    df.to_csv(outdir / 'all_sites_full.csv', index=False)

    for site in sorted(df['site'].unique()):
        df_site = df[df['site'] == site]
        df_site.to_csv(outdir / f'site{site}_hits.csv', index=False)

    with open(outdir / 'summary.md', 'w') as f:
        f.write('# TIMP2 Site-Aware Triage Summary\n\n')
        f.write(f'Total entries: {len(df)}\n\n')
        for site in sorted(df['site'].unique()):
            n = len(df[df['site']==site])
            cns = (df['site']==site) & (df['lane']=='CNS')
            per = (df['site']==site) & (df['lane']=='PERIPH')
            f.write(f'- Site {site}: {n} total ({cns.sum()} CNS, {per.sum()} Peripheral)\n')

    print(f"Exported site-aware triage results to: {outdir}")

# ------------------------------------------------------
# CLI
# ------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description='TIMP2 site-aware triage pipeline')
    ap.add_argument('--hits-csv', type=Path, help='CSV with id,smiles,docking_score (and optional site)')
    ap.add_argument('--scores-csv', type=Path, help='scores_all.csv (ligand,best_kcal_mol)')
    ap.add_argument('--chunks-dir', type=Path, help='Directory containing chunk .smi files')
    ap.add_argument('--pattern', default='strict_*.smi')
    ap.add_argument('--top', type=int, default=5000)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()

    if args.hits_csv and args.hits_csv.exists():
        df = pd.read_csv(args.hits_csv)
    elif args.scores_csv and args.chunks_dir:
        name2smi = build_name_to_smiles(args.chunks_dir, args.pattern)
        rows = []
        with open(args.scores_csv, 'r') as f:
            header = [h.strip() for h in f.readline().split(',')]
            # Flexible header mapping: supports 'ligand' or 'id', 'best_kcal_mol' or 'docking_score'
            colmap = {h.lower(): i for i, h in enumerate(header)}
            i_lig = colmap.get('ligand', colmap.get('id'))
            i_score = colmap.get('best_kcal_mol', colmap.get('docking_score'))
            if i_lig is None or i_score is None:
                raise ValueError('Could not find ligand/id or score column in CSV header')
            for line in f:
                parts = line.strip().split(',')
                if len(parts) <= max(i_lig, i_score): continue
                lig, score = parts[i_lig], float(parts[i_score])
                smi = name2smi.get(lig)
                if smi:
                    rows.append({'id': lig, 'smiles': smi, 'docking_score': score})
        df = pd.DataFrame(rows)
        if args.top > 0:
            df = df.nsmallest(args.top, 'docking_score')
    else:
        raise SystemExit('Provide either --hits-csv or --scores-csv + --chunks-dir')

    if 'site' not in df.columns:
        df['site'] = df['id'].apply(extract_site_from_id)

    triage_sites(df, args.out)

if __name__ == '__main__':
    main()
