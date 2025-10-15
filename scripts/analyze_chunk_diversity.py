#!/usr/bin/env python3
"""
Analyze chemical diversity across chunks by examining ZINC ID patterns and chemical properties.
This helps identify which chunks to prioritize for docking.
"""

import sys
from pathlib import Path
import random
from collections import defaultdict
import re

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors, AllChem
    from rdkit.Chem.Scaffolds import MurckoScaffold
    from rdkit import DataStructs
    RDKIT_AVAILABLE = True
except ImportError:
    print("Warning: RDKit not available. Install with: conda install -c conda-forge rdkit")
    RDKIT_AVAILABLE = False

def analyze_zinc_ids(smi_file, sample_size=1000):
    """
    Analyze ZINC ID patterns in a chunk.
    ZINC IDs often encode information about the source/synthesis route.
    """
    zinc_patterns = defaultdict(int)
    zinc_prefixes = defaultdict(int)
    
    with open(smi_file, 'r', encoding='utf-8', errors='ignore') as f:
        for i, line in enumerate(f):
            if i >= sample_size:
                break
            parts = line.strip().split()
            if len(parts) >= 2:
                zinc_id = parts[1]
                # Extract ZINC prefix patterns (first 6-8 digits often indicate source)
                if zinc_id.startswith('ZINC'):
                    # Get the numeric part
                    num_part = zinc_id[4:]
                    if len(num_part) >= 6:
                        prefix = num_part[:6]
                        zinc_prefixes[prefix] += 1
                        # Also track broader patterns (first 3 digits)
                        broad_prefix = num_part[:3]
                        zinc_patterns[broad_prefix] += 1
    
    return zinc_patterns, zinc_prefixes

def calculate_chemical_fingerprint(smi_file, sample_size=1000):
    """
    Calculate average chemical properties and scaffold diversity for a chunk.
    """
    if not RDKIT_AVAILABLE:
        return None
    
    properties = {
        'mw': [],
        'logp': [],
        'tpsa': [],
        'rotb': [],
        'hbd': [],
        'hba': [],
        'arom_rings': [],
        'fsp3': []
    }
    
    scaffolds = set()
    fingerprints = []
    
    with open(smi_file, 'r', encoding='utf-8', errors='ignore') as f:
        for i, line in enumerate(f):
            if i >= sample_size:
                break
            
            parts = line.strip().split()
            if len(parts) >= 1:
                smiles = parts[0]
                mol = Chem.MolFromSmiles(smiles)
                
                if mol:
                    # Calculate properties
                    properties['mw'].append(Descriptors.MolWt(mol))
                    properties['logp'].append(Descriptors.MolLogP(mol))
                    properties['tpsa'].append(rdMolDescriptors.CalcTPSA(mol))
                    properties['rotb'].append(rdMolDescriptors.CalcNumRotatableBonds(mol))
                    properties['hbd'].append(rdMolDescriptors.CalcNumHBD(mol))
                    properties['hba'].append(rdMolDescriptors.CalcNumHBA(mol))
                    properties['arom_rings'].append(rdMolDescriptors.CalcNumAromaticRings(mol))
                    properties['fsp3'].append(Descriptors.FractionCSP3(mol))
                    
                    # Get Murcko scaffold
                    try:
                        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
                        scaffolds.add(Chem.MolToSmiles(scaffold))
                    except:
                        pass
                    
                    # Generate Morgan fingerprint
                    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
                    fingerprints.append(fp)
    
    # Calculate averages
    avg_properties = {}
    for key, values in properties.items():
        if values:
            avg_properties[key] = sum(values) / len(values)
    
    return {
        'avg_properties': avg_properties,
        'scaffold_count': len(scaffolds),
        'fingerprints': fingerprints[:100]  # Keep subset for similarity
    }

def calculate_chunk_similarity(fp_list1, fp_list2):
    """
    Calculate Tanimoto similarity between two sets of fingerprints.
    """
    if not fp_list1 or not fp_list2:
        return 0.0
    
    similarities = []
    for fp1 in fp_list1[:50]:  # Sample to speed up
        for fp2 in fp_list2[:50]:
            sim = DataStructs.TanimotoSimilarity(fp1, fp2)
            similarities.append(sim)
    
    return sum(similarities) / len(similarities) if similarities else 0.0

def select_diverse_chunks(chunks_dir, n_select=6, sample_size=1000):
    """
    Analyze all chunks and select the most diverse subset.
    """
    chunks_dir = Path(chunks_dir)
    smi_files = sorted(chunks_dir.glob("strict_*.smi"))
    
    if not smi_files:
        print(f"No strict_*.smi files found in {chunks_dir}")
        return []
    
    print(f"Found {len(smi_files)} chunks to analyze\n")
    print("Analyzing chunks (this may take a few minutes)...")
    
    # Analyze each chunk
    chunk_data = {}
    for smi_file in smi_files:
        chunk_name = smi_file.stem
        print(f"  Analyzing {chunk_name}...", end=" ")
        
        # Get ZINC ID patterns
        zinc_patterns, zinc_prefixes = analyze_zinc_ids(smi_file, sample_size)
        
        # Get chemical fingerprint
        chem_data = calculate_chemical_fingerprint(smi_file, sample_size)
        
        chunk_data[chunk_name] = {
            'file': smi_file,
            'zinc_patterns': zinc_patterns,
            'zinc_prefixes': zinc_prefixes,
            'chem_data': chem_data
        }
        
        if chem_data and chem_data.get('avg_properties'):
            print(f"MW={chem_data['avg_properties']['mw']:.0f}, "
                  f"Scaffolds={chem_data['scaffold_count']}")
        else:
            print("Done")
    
    print("\n" + "="*80)
    print("CHUNK DIVERSITY ANALYSIS")
    print("="*80)
    
    # Strategy 1: Select based on ZINC ID distribution
    print("\n1. ZINC ID Pattern Distribution:")
    print("-"*40)
    
    # Group chunks by dominant ZINC prefix
    prefix_groups = defaultdict(list)
    for chunk_name, data in chunk_data.items():
        if data['zinc_patterns']:
            # Get most common prefix
            top_prefix = max(data['zinc_patterns'].items(), key=lambda x: x[1])[0]
            prefix_groups[top_prefix].append(chunk_name)
    
    for prefix, chunks in sorted(prefix_groups.items()):
        print(f"  ZINC**{prefix}***: {', '.join(chunks[:5])}{'...' if len(chunks) > 5 else ''}")
    
    # Strategy 2: Select based on chemical property spread
    if RDKIT_AVAILABLE:
        print("\n2. Chemical Property Ranges:")
        print("-"*40)
        
        # Find chunks with extreme property values
        property_extremes = defaultdict(lambda: {'min': (None, float('inf')), 
                                                 'max': (None, float('-inf'))})
        
        for chunk_name, data in chunk_data.items():
            if data['chem_data'] and data['chem_data'].get('avg_properties'):
                for prop, value in data['chem_data']['avg_properties'].items():
                    if value < property_extremes[prop]['min'][1]:
                        property_extremes[prop]['min'] = (chunk_name, value)
                    if value > property_extremes[prop]['max'][1]:
                        property_extremes[prop]['max'] = (chunk_name, value)
        
        for prop in ['mw', 'logp', 'tpsa', 'rotb']:
            if prop in property_extremes:
                min_chunk, min_val = property_extremes[prop]['min']
                max_chunk, max_val = property_extremes[prop]['max']
                print(f"  {prop.upper():5s}: {min_val:6.1f} ({min_chunk}) to {max_val:6.1f} ({max_chunk})")
    
    # Select diverse chunks
    print(f"\n3. Recommended {n_select} Diverse Chunks:")
    print("-"*40)
    
    selected = []
    
    # Method 1: Even spacing across chunk numbers (simple but effective)
    total_chunks = len(smi_files)
    step = max(1, total_chunks // (n_select - 1))
    for i in range(0, total_chunks, step):
        if len(selected) < n_select:
            selected.append(smi_files[i].stem)
    
    # Method 2: Select based on ZINC prefix diversity (more sophisticated)
    diverse_selected = []
    used_prefixes = set()
    
    # First pass: get one from each major prefix group
    for prefix, chunks in prefix_groups.items():
        if len(diverse_selected) < n_select and prefix not in used_prefixes:
            # Pick the middle chunk from this group
            chunk = chunks[len(chunks)//2]
            diverse_selected.append(chunk)
            used_prefixes.add(prefix)
    
    # Fill remaining slots with evenly spaced chunks
    if len(diverse_selected) < n_select:
        remaining = [c for c in selected if c not in diverse_selected]
        diverse_selected.extend(remaining[:n_select - len(diverse_selected)])
    
    print("\nRecommended chunks (evenly spaced):")
    for i, chunk in enumerate(selected[:n_select], 1):
        print(f"  {i}. {chunk}")
    
    print("\nAlternative selection (ZINC prefix diversity):")
    for i, chunk in enumerate(diverse_selected[:n_select], 1):
        data = chunk_data.get(chunk, {})
        if data.get('zinc_patterns'):
            top_pattern = max(data['zinc_patterns'].items(), key=lambda x: x[1])[0]
            print(f"  {i}. {chunk} (ZINC**{top_pattern}*** dominant)")
        else:
            print(f"  {i}. {chunk}")
    
    return selected[:n_select], diverse_selected[:n_select]

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Analyze chunk diversity and select subset for docking")
    parser.add_argument("chunks_dir", help="Directory containing strict_*.smi files")
    parser.add_argument("-n", "--num-chunks", type=int, default=6, 
                       help="Number of chunks to select (default: 6)")
    parser.add_argument("-s", "--sample-size", type=int, default=1000,
                       help="Number of molecules to sample per chunk (default: 1000)")
    args = parser.parse_args()
    
    selected, diverse = select_diverse_chunks(args.chunks_dir, args.num_chunks, args.sample_size)
    
    print("\n" + "="*80)
    print("NEXT STEPS:")
    print("="*80)
    print(f"1. Run docking on the {args.num_chunks} selected chunks")
    print(f"2. Expected yield: ~{args.num_chunks * 2360:,} high-quality candidates")
    print(f"3. Compute time: ~{args.num_chunks} weeks (vs 136 weeks for all)")
    print(f"4. Cost savings: {100 * (1 - args.num_chunks/136):.0f}%")

if __name__ == "__main__":
    main()