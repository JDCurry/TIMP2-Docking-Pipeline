#!/usr/bin/env python3
"""
Comprehensive analysis pipeline for TIMP2 virtual screening hits.
Performs filtering, clustering, ADMET prediction, and prioritization.
"""

import sys
import pandas as pd
import numpy as np
from pathlib import Path
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime

# RDKit imports
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Crippen, Lipinski, QED
    from rdkit.Chem.Scaffolds import MurckoScaffold
    from rdkit.Chem import AllChem, DataStructs
    from rdkit.ML.Cluster import Butina
    from rdkit.Chem import Draw
    from rdkit.Chem.Draw import IPythonConsole
    RDKIT_AVAILABLE = True
except ImportError:
    print("RDKit not installed. Install with: conda install -c conda-forge rdkit")
    RDKIT_AVAILABLE = False

# Optional: Advanced ADMET predictions
try:
    from chembl_structure_pipeline import checker, standardizer
    CHEMBL_AVAILABLE = True
except ImportError:
    print("ChEMBL structure pipeline not installed. Install with: pip install chembl_structure_pipeline")
    CHEMBL_AVAILABLE = False

class TIMP2CompoundAnalyzer:
    """Comprehensive analysis of TIMP2 virtual screening hits."""
    
    def __init__(self, input_csv, output_dir="timp2_analysis_results"):
        """
        Initialize analyzer with input CSV from triage.
        Expected columns: id, smiles, docking_score, mw, tpsa, hbd, hba, clogp, etc.
        """
        self.df = pd.read_csv(input_csv)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        print(f"Loaded {len(self.df)} compounds from {input_csv}")
        print(f"Columns available: {', '.join(self.df.columns)}")
        
        # Store results
        self.results = {}
        
    def calculate_additional_properties(self):
        """Calculate additional molecular properties for filtering."""
        print("\nCalculating additional molecular properties...")
        
        new_props = []
        for idx, row in self.df.iterrows():
            mol = Chem.MolFromSmiles(row['smiles'])
            if mol is None:
                new_props.append({})
                continue
                
            props = {
                'num_rings': Descriptors.RingCount(mol),
                'num_arom_rings': Descriptors.NumAromaticRings(mol),
                'num_heteroatoms': Descriptors.NumHeteroatoms(mol),
                'num_heavy_atoms': Descriptors.HeavyAtomCount(mol),
                'qed': QED.qed(mol),  # Quantitative Estimate of Drug-likeness
                'log_s': self.calculate_solubility(mol),  # Estimated solubility
                'psa': Descriptors.TPSA(mol),
                'fraction_sp3': Descriptors.FractionCSP3(mol)
            }
            
            # CNS MPO score (simplified version)
            props['cns_mpo'] = self.calculate_cns_mpo(mol, row)
            
            # Synthetic accessibility (simplified)
            props['sa_score'] = self.calculate_sa_score(mol)
            
            new_props.append(props)
            
            if (idx + 1) % 100 == 0:
                print(f"  Processed {idx + 1}/{len(self.df)} compounds...")
        
        # Add new properties to dataframe
        props_df = pd.DataFrame(new_props)
        self.df = pd.concat([self.df, props_df], axis=1)
        
        print(f"Added {len(props_df.columns)} new molecular properties")
        
    def calculate_solubility(self, mol):
        """Estimate aqueous solubility (LogS) using simple descriptors."""
        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        rb = Descriptors.NumRotatableBonds(mol)
        ap = len(mol.GetAromaticAtoms()) / mol.GetNumHeavyAtoms()
        
        # Simple ESOL-like equation
        log_s = 0.16 - 0.63 * logp - 0.0062 * mw + 0.066 * rb - 0.74 * ap
        return log_s
    
    def calculate_cns_mpo(self, mol, row):
        """
        Calculate CNS Multiparameter Optimization score.
        Based on Wager et al. (2016) ACS Chem. Neurosci.
        """
        scores = []
        
        # MW score
        mw = row.get('mw', Descriptors.MolWt(mol))
        if mw <= 360: scores.append(1.0)
        elif mw <= 500: scores.append(1.0 - (mw - 360) / 140)
        else: scores.append(0.0)
        
        # cLogP score  
        clogp = row.get('clogp', Crippen.MolLogP(mol))
        if clogp <= 3.0: scores.append(1.0)
        elif clogp <= 5.0: scores.append(1.0 - (clogp - 3.0) / 2.0)
        else: scores.append(0.0)
        
        # TPSA score
        tpsa = row.get('tpsa', Descriptors.TPSA(mol))
        if tpsa <= 75: scores.append(1.0)
        elif tpsa <= 120: scores.append(1.0 - (tpsa - 75) / 45)
        else: scores.append(0.0)
        
        # HBD score
        hbd = row.get('hbd', Lipinski.NumHDonors(mol))
        if hbd <= 1: scores.append(1.0)
        elif hbd <= 3: scores.append(1.0 - (hbd - 1) / 2.0)
        else: scores.append(0.0)
        
        # pKa - simplified, assume neutral
        scores.append(0.5)
        
        # HBA score
        hba = row.get('hba', Lipinski.NumHAcceptors(mol))
        if hba <= 5: scores.append(1.0)
        elif hba <= 8: scores.append(1.0 - (hba - 5) / 3.0)
        else: scores.append(0.0)
        
        return sum(scores)
    
    def calculate_sa_score(self, mol):
        """Simplified synthetic accessibility score (0=easy, 10=hard)."""
        # Very simplified - real SA score requires fragment database
        ring_count = Descriptors.RingCount(mol)
        complexity = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
        size_penalty = max(0, mol.GetNumHeavyAtoms() - 30) * 0.1
        
        return min(10, ring_count * 0.5 + complexity * 2 + size_penalty)
    
    def filter_compounds(self, strict=True):
        """Apply ADMET and drug-likeness filters."""
        print("\nApplying compound filters...")
        
        initial_count = len(self.df)
        
        if strict:
            # Strict CNS-focused filters
            filters = {
                'docking_score': ('<=', -7.5),
                'mw': ('between', [250, 450]),
                'clogp': ('between', [2.0, 3.5]),
                'tpsa': ('<=', 75),
                'hbd': ('<=', 1),
                'hba': ('<=', 7),
                'rotb': ('<=', 7),
                'qed': ('>=', 0.5),
                'log_s': ('>=', -5),
                'cns_mpo': ('>=', 4.0),
                'sa_score': ('<=', 5)
            }
        else:
            # Relaxed filters
            filters = {
                'docking_score': ('<=', -6.5),
                'mw': ('between', [200, 500]),
                'clogp': ('between', [1.5, 4.5]),
                'tpsa': ('<=', 90),
                'hbd': ('<=', 2),
                'hba': ('<=', 10),
                'qed': ('>=', 0.4),
                'log_s': ('>=', -6)
            }
        
        filtered_df = self.df.copy()
        filter_stats = {}
        
        for prop, (op, value) in filters.items():
            if prop not in filtered_df.columns:
                continue
                
            before = len(filtered_df)
            
            if op == '<=':
                filtered_df = filtered_df[filtered_df[prop] <= value]
            elif op == '>=':
                filtered_df = filtered_df[filtered_df[prop] >= value]
            elif op == 'between':
                filtered_df = filtered_df[
                    (filtered_df[prop] >= value[0]) & 
                    (filtered_df[prop] <= value[1])
                ]
            
            after = len(filtered_df)
            filter_stats[prop] = before - after
            
        self.filtered_df = filtered_df
        
        print(f"Filtering results:")
        print(f"  Initial: {initial_count} compounds")
        print(f"  After filters: {len(filtered_df)} compounds")
        print(f"  Removed: {initial_count - len(filtered_df)} compounds")
        print("\nCompounds removed by each filter:")
        for prop, count in filter_stats.items():
            if count > 0:
                print(f"  {prop}: {count}")
        
        return filtered_df
    
    def cluster_by_scaffold(self):
        """Group compounds by Murcko scaffold."""
        print("\nClustering by Murcko scaffold...")
        
        scaffolds = []
        for smiles in self.filtered_df['smiles']:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                scaffold = MurckoScaffold.GetScaffoldForMol(mol)
                scaffolds.append(Chem.MolToSmiles(scaffold))
            else:
                scaffolds.append('INVALID')
        
        self.filtered_df['scaffold'] = scaffolds
        
        # Count scaffold frequency
        scaffold_counts = self.filtered_df['scaffold'].value_counts()
        
        print(f"Unique scaffolds: {len(scaffold_counts)}")
        print(f"Top 10 scaffolds:")
        for scaffold, count in scaffold_counts.head(10).items():
            print(f"  {scaffold}: {count} compounds")
        
        # Select best compounds per scaffold
        best_per_scaffold = []
        for scaffold, group in self.filtered_df.groupby('scaffold'):
            # Take top 3 by docking score
            best = group.nsmallest(min(3, len(group)), 'docking_score')
            best_per_scaffold.append(best)
        
        self.scaffold_filtered = pd.concat(best_per_scaffold)
        print(f"\nSelected {len(self.scaffold_filtered)} compounds (best 3 per scaffold)")
        
        return self.scaffold_filtered
    
    def cluster_by_fingerprint(self, cutoff=0.6):
        """Cluster compounds by Tanimoto similarity."""
        print(f"\nClustering by fingerprint similarity (cutoff={cutoff})...")
        
        # Generate fingerprints
        fps = []
        valid_indices = []
        for idx, smiles in enumerate(self.scaffold_filtered['smiles']):
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
                fps.append(fp)
                valid_indices.append(idx)
        
        # Calculate distance matrix
        dists = []
        for i in range(len(fps)):
            for j in range(i):
                dists.append(1 - DataStructs.TanimotoSimilarity(fps[i], fps[j]))
        
        # Butina clustering
        clusters = Butina.ClusterData(dists, len(fps), cutoff, isDistData=True)
        
        print(f"Generated {len(clusters)} clusters")
        print(f"Cluster sizes: {[len(c) for c in clusters[:10]]}...")
        
        # Select representative from each cluster
        representatives = []
        for cluster in clusters:
            cluster_indices = [valid_indices[i] for i in cluster]
            cluster_df = self.scaffold_filtered.iloc[cluster_indices]
            # Select best by docking score
            best = cluster_df.nsmallest(1, 'docking_score')
            representatives.append(best)
        
        self.diverse_compounds = pd.concat(representatives)
        print(f"Selected {len(self.diverse_compounds)} diverse compounds")
        
        return self.diverse_compounds
    
    def prioritize_compounds(self, weights=None):
        """
        Create final prioritization score combining multiple factors.
        """
        print("\nPrioritizing compounds...")
        
        if weights is None:
            weights = {
                'docking_score': 0.35,
                'cns_mpo': 0.25,
                'qed': 0.15,
                'log_s': 0.10,
                'sa_score': 0.15
            }
        
        df = self.diverse_compounds.copy()
        
        # Normalize scores (0-1, higher is better)
        df['norm_docking'] = 1 - (df['docking_score'] - df['docking_score'].min()) / \
                                  (df['docking_score'].max() - df['docking_score'].min())
        df['norm_cns_mpo'] = df['cns_mpo'] / 6.0
        df['norm_qed'] = df['qed']
        df['norm_log_s'] = (df['log_s'] - df['log_s'].min()) / \
                           (df['log_s'].max() - df['log_s'].min())
        df['norm_sa'] = 1 - df['sa_score'] / 10.0
        
        # Calculate composite score
        df['priority_score'] = (
            weights['docking_score'] * df['norm_docking'] +
            weights['cns_mpo'] * df['norm_cns_mpo'] +
            weights['qed'] * df['norm_qed'] +
            weights['log_s'] * df['norm_log_s'] +
            weights['sa_score'] * df['norm_sa']
        )
        
        # Sort by priority score
        df = df.sort_values('priority_score', ascending=False)
        
        self.prioritized_compounds = df
        
        print(f"Top 10 compounds by priority score:")
        for idx, row in df.head(10).iterrows():
            print(f"  {row['id']}: Score={row['priority_score']:.3f}, "
                  f"Dock={row['docking_score']:.2f}, CNS_MPO={row['cns_mpo']:.2f}")
        
        return df
    
    def visualize_results(self):
        """Create comprehensive visualizations."""
        print("\nGenerating visualizations...")
        
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        
        # 1. Docking score distribution
        axes[0, 0].hist(self.df['docking_score'], bins=30, alpha=0.7, color='blue')
        axes[0, 0].axvline(self.filtered_df['docking_score'].median(), 
                           color='red', linestyle='--', label='Median (filtered)')
        axes[0, 0].set_xlabel('Docking Score (kcal/mol)')
        axes[0, 0].set_ylabel('Count')
        axes[0, 0].set_title('Docking Score Distribution')
        axes[0, 0].legend()
        
        # 2. MW vs cLogP with CNS box
        axes[0, 1].scatter(self.df['mw'], self.df['clogp'], 
                          alpha=0.5, s=20, label='All')
        axes[0, 1].scatter(self.prioritized_compounds['mw'], 
                          self.prioritized_compounds['clogp'],
                          color='red', s=50, label='Prioritized')
        axes[0, 1].add_patch(plt.Rectangle((250, 2.0), 200, 1.5, 
                                           fill=False, edgecolor='green', 
                                           linewidth=2, label='CNS space'))
        axes[0, 1].set_xlabel('Molecular Weight')
        axes[0, 1].set_ylabel('cLogP')
        axes[0, 1].set_title('Chemical Space')
        axes[0, 1].legend()
        
        # 3. TPSA vs HBD
        axes[0, 2].scatter(self.df['tpsa'], self.df['hbd'], 
                          alpha=0.5, s=20, label='All')
        axes[0, 2].scatter(self.prioritized_compounds['tpsa'], 
                          self.prioritized_compounds['hbd'],
                          color='red', s=50, label='Prioritized')
        axes[0, 2].axvline(70, color='green', linestyle='--', alpha=0.5)
        axes[0, 2].axhline(1, color='green', linestyle='--', alpha=0.5)
        axes[0, 2].set_xlabel('TPSA')
        axes[0, 2].set_ylabel('HBD')
        axes[0, 2].set_title('Polarity Profile')
        axes[0, 2].legend()
        
        # 4. CNS MPO distribution
        if 'cns_mpo' in self.filtered_df.columns:
            axes[1, 0].hist(self.filtered_df['cns_mpo'], bins=20, 
                           alpha=0.7, color='purple')
            axes[1, 0].axvline(4.0, color='red', linestyle='--', 
                              label='Desirable threshold')
            axes[1, 0].set_xlabel('CNS MPO Score')
            axes[1, 0].set_ylabel('Count')
            axes[1, 0].set_title('CNS MPO Distribution')
            axes[1, 0].legend()
        
        # 5. QED vs Solubility
        if 'qed' in self.filtered_df.columns and 'log_s' in self.filtered_df.columns:
            axes[1, 1].scatter(self.filtered_df['qed'], 
                              self.filtered_df['log_s'],
                              c=self.filtered_df['docking_score'],
                              cmap='viridis', alpha=0.6)
            axes[1, 1].set_xlabel('QED')
            axes[1, 1].set_ylabel('LogS')
            axes[1, 1].set_title('Drug-likeness vs Solubility')
            plt.colorbar(axes[1, 1].collections[0], ax=axes[1, 1], 
                        label='Docking Score')
        
        # 6. Priority score components
        if len(self.prioritized_compounds) > 0:
            top10 = self.prioritized_compounds.head(10)
            components = ['norm_docking', 'norm_cns_mpo', 'norm_qed', 
                         'norm_log_s', 'norm_sa']
            component_data = top10[components].values.T
            
            x = np.arange(len(top10))
            width = 0.15
            
            for i, comp in enumerate(components):
                axes[1, 2].bar(x + i*width, component_data[i], 
                              width, label=comp.replace('norm_', ''))
            
            axes[1, 2].set_xlabel('Top 10 Compounds')
            axes[1, 2].set_ylabel('Normalized Score')
            axes[1, 2].set_title('Priority Score Components')
            axes[1, 2].legend(fontsize=8)
            axes[1, 2].set_xticks(x + width * 2)
            axes[1, 2].set_xticklabels([f"#{i+1}" for i in range(10)])
        
        plt.tight_layout()
        plt.savefig(self.output_dir / 'compound_analysis.png', dpi=300)
        plt.show()
        
        print(f"Visualizations saved to {self.output_dir}/compound_analysis.png")
    
    def export_results(self):
        """Export all results to files."""
        print("\nExporting results...")
        
        # Export filtered compounds
        self.filtered_df.to_csv(
            self.output_dir / 'filtered_compounds.csv', index=False
        )
        
        # Export scaffold analysis
        scaffold_summary = self.filtered_df['scaffold'].value_counts().to_frame()
        scaffold_summary.to_csv(self.output_dir / 'scaffold_summary.csv')
        
        # Export prioritized compounds
        self.prioritized_compounds.to_csv(
            self.output_dir / 'prioritized_compounds.csv', index=False
        )
        
        # Export top 25 for manual inspection
        top25 = self.prioritized_compounds.head(25)[
            ['id', 'smiles', 'docking_score', 'mw', 'clogp', 'tpsa', 
             'hbd', 'cns_mpo', 'qed', 'priority_score']
        ]
        top25.to_csv(self.output_dir / 'top25_for_visualization.csv', index=False)
        
        # Generate SDF for top compounds (for Maestro/PyMOL)
        writer = Chem.SDWriter(str(self.output_dir / 'top25_compounds.sdf'))
        for idx, row in top25.iterrows():
            mol = Chem.MolFromSmiles(row['smiles'])
            if mol:
                mol.SetProp('_Name', str(row['id']))
                mol.SetProp('docking_score', str(row['docking_score']))
                mol.SetProp('priority_score', str(row['priority_score']))
                writer.write(mol)
        writer.close()
        
        # Generate summary report
        with open(self.output_dir / 'analysis_summary.txt', 'w') as f:
            f.write("TIMP2 Compound Analysis Summary\n")
            f.write("=" * 50 + "\n\n")
            f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            f.write("Pipeline Statistics:\n")
            f.write(f"  Initial compounds: {len(self.df)}\n")
            f.write(f"  After filtering: {len(self.filtered_df)}\n")
            f.write(f"  After scaffold selection: {len(self.scaffold_filtered)}\n")
            f.write(f"  After diversity selection: {len(self.diverse_compounds)}\n")
            f.write(f"  Final prioritized: {len(self.prioritized_compounds)}\n\n")
            
            f.write("Top 10 Compounds:\n")
            for idx, row in self.prioritized_compounds.head(10).iterrows():
                f.write(f"  {row['id']}: {row['smiles']}\n")
                f.write(f"    Docking: {row['docking_score']:.2f} kcal/mol\n")
                f.write(f"    Priority: {row['priority_score']:.3f}\n")
                f.write(f"    CNS MPO: {row.get('cns_mpo', 'N/A')}\n\n")
        
        print(f"Results exported to {self.output_dir}/")
        print(f"  - filtered_compounds.csv: All compounds passing filters")
        print(f"  - prioritized_compounds.csv: Final prioritized list")
        print(f"  - top25_for_visualization.csv: Top compounds for manual inspection")
        print(f"  - top25_compounds.sdf: SDF file for Maestro/PyMOL")
        print(f"  - analysis_summary.txt: Summary report")
    
    def run_full_pipeline(self, strict_filter=True):
        """Run the complete analysis pipeline."""
        print("\n" + "="*50)
        print("Running TIMP2 Compound Analysis Pipeline")
        print("="*50)
        
        # Calculate properties
        self.calculate_additional_properties()
        
        # Apply filters
        self.filter_compounds(strict=strict_filter)
        
        # Scaffold clustering
        self.cluster_by_scaffold()
        
        # Diversity selection
        self.cluster_by_fingerprint(cutoff=0.6)
        
        # Prioritization
        self.prioritize_compounds()
        
        # Visualizations
        self.visualize_results()
        
        # Export results
        self.export_results()
        
        print("\n" + "="*50)
        print("Pipeline complete!")
        print("="*50)
        
        return self.prioritized_compounds


def main():
    """Main execution function."""
    import argparse
    
    parser = argparse.ArgumentParser(description='Analyze TIMP2 virtual screening hits')
    parser.add_argument('input_csv', help='Path to CNS hits CSV from triage')
    parser.add_argument('-o', '--output', default='timp2_analysis_results',
                       help='Output directory (default: timp2_analysis_results)')
    parser.add_argument('--relaxed', action='store_true',
                       help='Use relaxed filters instead of strict CNS filters')
    args = parser.parse_args()
    
    # Run analysis
    analyzer = TIMP2CompoundAnalyzer(args.input_csv, args.output)
    results = analyzer.run_full_pipeline(strict_filter=not args.relaxed)
    
    print(f"\nTop compound for visualization:")
    if len(results) > 0:
        top = results.iloc[0]
        print(f"ID: {top['id']}")
        print(f"SMILES: {top['smiles']}")
        print(f"Docking Score: {top['docking_score']:.2f} kcal/mol")
        print(f"Priority Score: {top['priority_score']:.3f}")

if __name__ == "__main__":
    main()