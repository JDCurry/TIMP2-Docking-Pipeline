#!/usr/bin/env python3
"""
Add ZINC IDs to prioritized compound files
Maps lig_XXXXXXXX or other ligand names back to ZINC IDs using the original SMILES file
"""

import pandas as pd
import argparse
import sys
from pathlib import Path

def add_zinc_ids_to_prioritized(prioritized_csv, original_smi, output_csv):
    """
    Add ZINC IDs to a prioritized compounds CSV file
    
    Args:
        prioritized_csv: The prioritized compounds/singletons file
        original_smi: The original strict_XXXX.smi file with ZINC IDs
        output_csv: Output file with ZINC IDs added
    """
    
    print("=" * 60)
    print("ADD ZINC IDs TO PRIORITIZED COMPOUNDS")
    print("=" * 60)
    
    # Read the prioritized file
    print(f"\nReading {prioritized_csv}...")
    df = pd.read_csv(prioritized_csv)
    print(f"  Total entries: {len(df)}")
    
    # Read ZINC IDs from original SMILES file
    print(f"\nReading ZINC IDs from {original_smi}...")
    zinc_mapping = {}  # Maps lig_XXXXXXXX or compound name to ZINC ID
    smiles_mapping = {}  # Maps to SMILES
    
    with open(original_smi, 'r') as f:
        for idx, line in enumerate(f, 1):
            parts = line.strip().split()
            if len(parts) >= 2:
                smiles = parts[0]
                zinc_id = parts[1]
            else:
                zinc_id = f"mol_{idx:08d}"
                smiles = parts[0] if parts else ""
            
            # Create mapping for lig_00000001 format (corresponds to line number)
            lig_name = f"lig_{idx:08d}"
            zinc_mapping[lig_name] = zinc_id
            smiles_mapping[lig_name] = smiles
    
    print(f"  Loaded {len(zinc_mapping)} ZINC IDs")
    
    # Determine which column contains the ligand identifier
    id_column = None
    for col in ['id', 'ligand', 'compound_id', 'name']:
        if col in df.columns:
            id_column = col
            print(f"\nUsing column '{col}' as identifier")
            break
    
    if id_column is None:
        print("\nError: Could not find identifier column (tried: id, ligand, compound_id, name)")
        return None
    
    # Add ZINC IDs to dataframe
    print("\nMapping identifiers to ZINC IDs...")
    
    # Direct mapping first
    df['zinc_id'] = df[id_column].map(zinc_mapping)
    df['smiles'] = df[id_column].map(smiles_mapping)
    
    # Handle site-specific names (e.g., site1_lig_00000001)
    if df['zinc_id'].isna().any():
        print("  Attempting to extract lig_XXXXXXXX from complex names...")
        df['lig_base'] = df[id_column].str.extract(r'(lig_\d{8})')
        unmapped_before = df['zinc_id'].isna().sum()
        
        df.loc[df['zinc_id'].isna(), 'zinc_id'] = df.loc[df['zinc_id'].isna(), 'lig_base'].map(zinc_mapping)
        df.loc[df['smiles'].isna(), 'smiles'] = df.loc[df['smiles'].isna(), 'lig_base'].map(smiles_mapping)
        
        unmapped_after = df['zinc_id'].isna().sum()
        print(f"    Recovered {unmapped_before - unmapped_after} additional mappings")
        df = df.drop('lig_base', axis=1)
    
    # Handle ZINC IDs already in the identifier (e.g., ZINC001234567_1)
    if df['zinc_id'].isna().any():
        print("  Checking if identifiers already contain ZINC IDs...")
        zinc_pattern = df[id_column].str.extract(r'(ZINC\d+)')
        has_zinc = zinc_pattern[0].notna()
        if has_zinc.any():
            df.loc[df['zinc_id'].isna() & has_zinc, 'zinc_id'] = zinc_pattern.loc[has_zinc, 0]
            print(f"    Found {has_zinc.sum()} entries with ZINC IDs in identifier")
    
    # Report results
    mapped = df['zinc_id'].notna().sum()
    unmapped = df['zinc_id'].isna().sum()
    
    print("\n" + "=" * 60)
    print("RESULTS")
    print("=" * 60)
    print(f"Successfully mapped: {mapped} / {len(df)} entries ({100*mapped/len(df):.1f}%)")
    if unmapped > 0:
        print(f"Unable to map: {unmapped} entries")
        print("\nSample unmapped identifiers:")
        print(df[df['zinc_id'].isna()][id_column].head(5).tolist())
    
    # Reorder columns to put zinc_id near the front
    cols = df.columns.tolist()
    if 'zinc_id' in cols:
        cols.remove('zinc_id')
        # Insert after id column
        id_col_idx = cols.index(id_column)
        cols.insert(id_col_idx + 1, 'zinc_id')
    if 'smiles' in cols:
        cols.remove('smiles')
        cols.insert(cols.index('zinc_id') + 1, 'smiles')
    df = df[cols]
    
    # Save enhanced file
    df.to_csv(output_csv, index=False)
    print(f"\n✓ Output saved to: {output_csv}")
    print("=" * 60)
    
    return df

def main():
    parser = argparse.ArgumentParser(
        description="Add ZINC IDs to prioritized compound files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Add ZINC IDs to prioritized compounds
  python add_zinc_ids_to_prioritized.py prioritized_compounds.csv strict_0120.smi prioritized_compounds_with_zinc.csv
  
  # Add ZINC IDs to prioritized singletons
  python add_zinc_ids_to_prioritized.py prioritized_singletons.csv strict_0120.smi prioritized_singletons_with_zinc.csv
        """
    )
    
    parser.add_argument(
        "prioritized_csv",
        type=str,
        help="Path to prioritized compounds or singletons CSV file"
    )
    parser.add_argument(
        "original_smi",
        type=str,
        help="Path to original SMILES file with ZINC IDs (e.g., strict_0120.smi)"
    )
    parser.add_argument(
        "output_csv",
        type=str,
        help="Output CSV file with ZINC IDs added"
    )
    
    args = parser.parse_args()
    
    # Check if files exist
    if not Path(args.prioritized_csv).exists():
        print(f"Error: File not found: {args.prioritized_csv}")
        return 1
    
    if not Path(args.original_smi).exists():
        print(f"Error: File not found: {args.original_smi}")
        return 1
    
    # Add ZINC IDs
    result = add_zinc_ids_to_prioritized(args.prioritized_csv, args.original_smi, args.output_csv)
    
    if result is None:
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())
