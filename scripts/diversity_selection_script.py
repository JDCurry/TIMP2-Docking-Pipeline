#!/usr/bin/env python3
"""
Comprehensive Compound Diversity Selection Script for TIMP2 Virtual Screening
Selects diverse compounds from docking results balancing singleton and scaffold representation
"""

import pandas as pd
import numpy as np
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from rdkit.Chem.Scaffolds import MurckoScaffold
import argparse
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

class DiversitySelector:
    """
    Select diverse compounds from virtual screening results using a balanced approach
    between singleton scaffolds and populated scaffold families.
    """
    
    def __init__(self, singletons_file, non_singletons_file, output_dir="diversity_selection"):
        """
        Initialize with separate singleton and non-singleton files.
        """
        self.singletons_df = pd.read_csv(singletons_file)
        self.non_singletons_df = pd.read_csv(non_singletons_file)
        self.full_df = pd.concat([self.singletons_df, self.non_singletons_df], ignore_index=True)
        
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        print(f"Loaded {len(self.singletons_df)} singletons and {len(self.non_singletons_df)} non-singletons")
        print(f"Total compounds: {len(self.full_df)}")
        
        # Store selected compounds
        self.selected_compounds = None
        
    def calculate_tanimoto_similarity(self, smiles1, smiles2, radius=2, nBits=2048):
        """
        Calculate Tanimoto similarity between two molecules using Morgan fingerprints.
        """
        mol1 = Chem.MolFromSmiles(smiles1)
        mol2 = Chem.MolFromSmiles(smiles2)
        
        if mol1 is None or mol2 is None:
            return 0.0
        
        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, radius, nBits=nBits)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, radius, nBits=nBits)
        
        return DataStructs.TanimotoSimilarity(fp1, fp2)
    
    def get_scaffold_statistics(self):
        """
        Analyze scaffold distribution in the dataset.
        """
        # Get scaffold counts from non-singletons
        scaffold_counts = self.non_singletons_df.groupby('scaffold').size().sort_values(ascending=False)
        
        stats = {
            'total_compounds': len(self.full_df),
            'singleton_count': len(self.singletons_df),
            'singleton_percentage': (len(self.singletons_df) / len(self.full_df)) * 100,
            'unique_scaffolds': len(scaffold_counts) + len(self.singletons_df),
            'most_populated_scaffolds': scaffold_counts.head(10).to_dict(),
            'scaffold_size_distribution': {
                '1 (singleton)': len(self.singletons_df),
                '2-5': sum((scaffold_counts >= 2) & (scaffold_counts <= 5)),
                '6-10': sum((scaffold_counts >= 6) & (scaffold_counts <= 10)),
                '11-20': sum((scaffold_counts >= 11) & (scaffold_counts <= 20)),
                '>20': sum(scaffold_counts > 20)
            }
        }
        
        return stats
    
    def select_diverse_compounds(self, n_compounds=5, strategy='balanced', 
                                max_similarity=0.7, prioritize_score=True):
        """
        Select diverse compounds using specified strategy.
        
        Parameters:
        -----------
        n_compounds : int
            Number of compounds to select
        strategy : str
            'balanced' - mix of singletons and scaffold representatives
            'singleton_first' - prioritize unique scaffolds
            'scaffold_first' - prioritize validated scaffolds
            'pure_score' - select by score regardless of scaffold
        max_similarity : float
            Maximum Tanimoto similarity allowed between selected compounds
        prioritize_score : bool
            Whether to prioritize docking score in selection
        """
        
        selected = []
        selected_smiles = []
        
        if strategy == 'balanced':
            selected = self._select_balanced(n_compounds, max_similarity, prioritize_score)
            
        elif strategy == 'singleton_first':
            selected = self._select_singleton_first(n_compounds, max_similarity, prioritize_score)
            
        elif strategy == 'scaffold_first':
            selected = self._select_scaffold_first(n_compounds, max_similarity, prioritize_score)
            
        elif strategy == 'pure_score':
            selected = self._select_pure_score(n_compounds, max_similarity)
            
        else:
            raise ValueError(f"Unknown strategy: {strategy}")
        
        self.selected_compounds = pd.DataFrame(selected)
        return self.selected_compounds
    
    def _select_balanced(self, n_compounds, max_similarity, prioritize_score):
        """
        Balanced selection: mix of singletons and scaffold representatives.
        """
        selected = []
        selected_smiles = []
        
        # Step 1: Select top 2 singletons
        n_singletons = min(2, len(self.singletons_df), n_compounds)
        if prioritize_score:
            top_singletons = self.singletons_df.nsmallest(n_singletons * 2, 'docking_score')
        else:
            # Use composite score if available, otherwise use docking score
            if 'mpo_cns' in self.singletons_df.columns:
                # Create composite score
                temp_df = self.singletons_df.copy()
                temp_df['composite_score'] = (
                    -temp_df['docking_score'] * 0.5 + 
                    temp_df['mpo_cns'] * 0.3 +
                    (10 - temp_df.get('alerts', 0)) * 0.2
                )
                top_singletons = temp_df.nlargest(n_singletons * 2, 'composite_score')
            else:
                top_singletons = self.singletons_df.nsmallest(n_singletons * 2, 'docking_score')
        
        for _, compound in top_singletons.iterrows():
            if self._check_similarity(compound['smiles'], selected_smiles, max_similarity):
                selected.append(compound.to_dict())
                selected_smiles.append(compound['smiles'])
                if len(selected) >= n_singletons:
                    break
        
        # Step 2: Get scaffold families sorted by population
        scaffold_counts = self.non_singletons_df.groupby('scaffold').size().sort_values(ascending=False)
        
        # Step 3: Select best from top scaffold families
        remaining_slots = n_compounds - len(selected)
        scaffolds_checked = 0
        max_scaffolds_to_check = min(10, len(scaffold_counts))  # Don't check more than 10 scaffolds
        
        for scaffold, count in scaffold_counts.items():
            if len(selected) >= n_compounds:
                break
            if scaffolds_checked >= max_scaffolds_to_check:
                break
            scaffolds_checked += 1
            
            # Get best compound from this scaffold
            scaffold_compounds = self.non_singletons_df[self.non_singletons_df['scaffold'] == scaffold]
            
            if prioritize_score:
                best_from_scaffold = scaffold_compounds.nsmallest(3, 'docking_score')
            else:
                if 'mpo_cns' in scaffold_compounds.columns:
                    temp_df = scaffold_compounds.copy()
                    temp_df['composite_score'] = (
                        -temp_df['docking_score'] * 0.5 + 
                        temp_df['mpo_cns'] * 0.3 +
                        (10 - temp_df.get('alerts', 0)) * 0.2
                    )
                    best_from_scaffold = temp_df.nlargest(3, 'composite_score')
                else:
                    best_from_scaffold = scaffold_compounds.nsmallest(3, 'docking_score')
            
            for _, compound in best_from_scaffold.iterrows():
                if self._check_similarity(compound['smiles'], selected_smiles, max_similarity):
                    selected.append(compound.to_dict())
                    selected_smiles.append(compound['smiles'])
                    break
        
        # Step 4: If still need more, add more singletons
        if len(selected) < n_compounds:
            remaining_singletons = self.singletons_df[
                ~self.singletons_df['id'].isin([s['id'] for s in selected])
            ]
            
            if prioritize_score:
                additional = remaining_singletons.nsmallest(n_compounds - len(selected), 'docking_score')
            else:
                additional = remaining_singletons.sample(min(len(remaining_singletons), n_compounds - len(selected)))
            
            for _, compound in additional.iterrows():
                if self._check_similarity(compound['smiles'], selected_smiles, max_similarity):
                    selected.append(compound.to_dict())
                    selected_smiles.append(compound['smiles'])
                    if len(selected) >= n_compounds:
                        break
        
        return selected
    
    def _select_singleton_first(self, n_compounds, max_similarity, prioritize_score):
        """
        Prioritize singleton scaffolds for maximum novelty.
        """
        selected = []
        selected_smiles = []
        
        # Select as many singletons as possible
        if prioritize_score:
            candidates = self.singletons_df.nsmallest(min(len(self.singletons_df), n_compounds * 2), 'docking_score')
        else:
            candidates = self.singletons_df.sample(min(len(self.singletons_df), n_compounds * 2))
        
        for _, compound in candidates.iterrows():
            if self._check_similarity(compound['smiles'], selected_smiles, max_similarity):
                selected.append(compound.to_dict())
                selected_smiles.append(compound['smiles'])
                if len(selected) >= n_compounds:
                    break
        
        # Fill remaining with non-singletons if needed
        if len(selected) < n_compounds:
            remaining = n_compounds - len(selected)
            non_singleton_candidates = self.non_singletons_df.nsmallest(remaining * 2, 'docking_score')
            
            for _, compound in non_singleton_candidates.iterrows():
                if self._check_similarity(compound['smiles'], selected_smiles, max_similarity):
                    selected.append(compound.to_dict())
                    selected_smiles.append(compound['smiles'])
                    if len(selected) >= n_compounds:
                        break
        
        return selected
    
    def _select_scaffold_first(self, n_compounds, max_similarity, prioritize_score):
        """
        Prioritize compounds from validated scaffold families.
        """
        selected = []
        selected_smiles = []
        
        # Get scaffold families sorted by population
        scaffold_counts = self.non_singletons_df.groupby('scaffold').size().sort_values(ascending=False)
        
        # Select from most populated scaffolds
        for scaffold, count in scaffold_counts.items():
            if len(selected) >= n_compounds:
                break
            
            scaffold_compounds = self.non_singletons_df[self.non_singletons_df['scaffold'] == scaffold]
            best_compound = scaffold_compounds.nsmallest(1, 'docking_score').iloc[0]
            
            if self._check_similarity(best_compound['smiles'], selected_smiles, max_similarity):
                selected.append(best_compound.to_dict())
                selected_smiles.append(best_compound['smiles'])
        
        # Fill remaining with singletons if needed
        if len(selected) < n_compounds:
            remaining = n_compounds - len(selected)
            singleton_candidates = self.singletons_df.nsmallest(remaining * 2, 'docking_score')
            
            for _, compound in singleton_candidates.iterrows():
                if self._check_similarity(compound['smiles'], selected_smiles, max_similarity):
                    selected.append(compound.to_dict())
                    selected_smiles.append(compound['smiles'])
                    if len(selected) >= n_compounds:
                        break
        
        return selected
    
    def _select_pure_score(self, n_compounds, max_similarity):
        """
        Select purely based on docking score, ensuring diversity.
        """
        selected = []
        selected_smiles = []
        
        # Sort all compounds by docking score
        all_sorted = self.full_df.nsmallest(len(self.full_df), 'docking_score')
        
        for _, compound in all_sorted.iterrows():
            if self._check_similarity(compound['smiles'], selected_smiles, max_similarity):
                selected.append(compound.to_dict())
                selected_smiles.append(compound['smiles'])
                if len(selected) >= n_compounds:
                    break
        
        return selected
    
    def _check_similarity(self, smiles, selected_smiles, max_similarity):
        """
        Check if a compound is sufficiently different from already selected compounds.
        """
        if not selected_smiles:
            return True
        
        for selected_smile in selected_smiles:
            similarity = self.calculate_tanimoto_similarity(smiles, selected_smile)
            if similarity > max_similarity:
                return False
        
        return True
    
    def analyze_selection(self):
        """
        Analyze the diversity and properties of selected compounds.
        """
        if self.selected_compounds is None:
            print("No compounds selected yet. Run select_diverse_compounds() first.")
            return None
        
        analysis = {
            'n_selected': len(self.selected_compounds),
            'n_singletons': sum(self.selected_compounds['is_singleton'] == 'True'),
            'n_non_singletons': sum(self.selected_compounds['is_singleton'] == 'False'),
            'docking_scores': {
                'best': self.selected_compounds['docking_score'].min(),
                'worst': self.selected_compounds['docking_score'].max(),
                'mean': self.selected_compounds['docking_score'].mean(),
                'std': self.selected_compounds['docking_score'].std()
            }
        }
        
        # Add property statistics if available
        if 'mpo_cns' in self.selected_compounds.columns:
            analysis['cns_mpo'] = {
                'min': self.selected_compounds['mpo_cns'].min(),
                'max': self.selected_compounds['mpo_cns'].max(),
                'mean': self.selected_compounds['mpo_cns'].mean()
            }
        
        # Calculate pairwise similarities
        similarities = []
        for i in range(len(self.selected_compounds)):
            for j in range(i+1, len(self.selected_compounds)):
                sim = self.calculate_tanimoto_similarity(
                    self.selected_compounds.iloc[i]['smiles'],
                    self.selected_compounds.iloc[j]['smiles']
                )
                similarities.append(sim)
        
        if similarities:
            analysis['diversity'] = {
                'min_similarity': min(similarities),
                'max_similarity': max(similarities),
                'mean_similarity': np.mean(similarities)
            }
        
        return analysis
    
    def export_results(self, prefix="selected"):
        """
        Export selected compounds and analysis to files.
        """
        if self.selected_compounds is None:
            print("No compounds selected yet. Run select_diverse_compounds() first.")
            return
        
        # Export selected compounds
        output_file = self.output_dir / f"{prefix}_compounds.csv"
        self.selected_compounds.to_csv(output_file, index=False)
        print(f"Selected compounds saved to: {output_file}")
        
        # Export analysis
        analysis = self.analyze_selection()
        analysis_file = self.output_dir / f"{prefix}_analysis.txt"
        
        with open(analysis_file, 'w') as f:
            f.write(f"Diversity Selection Analysis\n")
            f.write(f"Generated: {datetime.now()}\n")
            f.write("="*60 + "\n\n")
            
            f.write(f"Selection Summary:\n")
            f.write(f"  Total selected: {analysis['n_selected']}\n")
            f.write(f"  Singletons: {analysis['n_singletons']}\n")
            f.write(f"  Non-singletons: {analysis['n_non_singletons']}\n\n")
            
            f.write(f"Docking Scores:\n")
            for key, value in analysis['docking_scores'].items():
                f.write(f"  {key}: {value:.3f}\n")
            f.write("\n")
            
            if 'cns_mpo' in analysis:
                f.write(f"CNS MPO Scores:\n")
                for key, value in analysis['cns_mpo'].items():
                    f.write(f"  {key}: {value:.3f}\n")
                f.write("\n")
            
            if 'diversity' in analysis:
                f.write(f"Chemical Diversity (Tanimoto):\n")
                for key, value in analysis['diversity'].items():
                    f.write(f"  {key}: {value:.3f}\n")
            
            f.write("\n" + "="*60 + "\n")
            f.write("Selected Compounds:\n\n")
            
            for idx, row in self.selected_compounds.iterrows():
                f.write(f"{idx+1}. {row['id']}\n")
                f.write(f"   SMILES: {row['smiles']}\n")
                f.write(f"   Docking Score: {row['docking_score']:.3f}\n")
                f.write(f"   Scaffold: {row.get('scaffold', 'N/A')}\n")
                f.write(f"   Singleton: {row.get('is_singleton', 'N/A')}\n")
                if 'mpo_cns' in row:
                    f.write(f"   CNS MPO: {row['mpo_cns']:.2f}\n")
                f.write("\n")
        
        print(f"Analysis saved to: {analysis_file}")
        
        # Export SDF if RDKit is available
        try:
            sdf_file = self.output_dir / f"{prefix}_compounds.sdf"
            writer = Chem.SDWriter(str(sdf_file))
            
            for idx, row in self.selected_compounds.iterrows():
                mol = Chem.MolFromSmiles(row['smiles'])
                if mol:
                    mol.SetProp('_Name', row['id'])
                    mol.SetProp('docking_score', str(row['docking_score']))
                    mol.SetProp('singleton', str(row.get('is_singleton', 'N/A')))
                    writer.write(mol)
            
            writer.close()
            print(f"SDF file saved to: {sdf_file}")
        except Exception as e:
            print(f"Could not export SDF: {e}")
    
    def compare_strategies(self, n_compounds=5, max_similarity=0.7):
        """
        Compare different selection strategies and return results.
        """
        strategies = ['balanced', 'singleton_first', 'scaffold_first', 'pure_score']
        results = {}
        
        for strategy in strategies:
            print(f"\nTesting strategy: {strategy}")
            selected = self.select_diverse_compounds(
                n_compounds=n_compounds,
                strategy=strategy,
                max_similarity=max_similarity
            )
            
            analysis = self.analyze_selection()
            results[strategy] = {
                'compounds': selected,
                'analysis': analysis
            }
            
            # Export each strategy's results
            self.export_results(prefix=f"{strategy}_{n_compounds}")
        
        # Create comparison report
        comparison_file = self.output_dir / f"strategy_comparison_{n_compounds}.txt"
        with open(comparison_file, 'w') as f:
            f.write("Strategy Comparison Report\n")
            f.write("="*60 + "\n\n")
            
            for strategy, data in results.items():
                f.write(f"\n{strategy.upper()} Strategy:\n")
                f.write("-"*40 + "\n")
                
                analysis = data['analysis']
                f.write(f"Singletons: {analysis['n_singletons']}/{analysis['n_selected']}\n")
                f.write(f"Best docking score: {analysis['docking_scores']['best']:.3f}\n")
                f.write(f"Mean docking score: {analysis['docking_scores']['mean']:.3f}\n")
                
                if 'diversity' in analysis:
                    f.write(f"Mean Tanimoto similarity: {analysis['diversity']['mean_similarity']:.3f}\n")
                
                f.write("\nCompound IDs:\n")
                for idx, row in data['compounds'].iterrows():
                    f.write(f"  - {row['id']} (score: {row['docking_score']:.2f})\n")
        
        print(f"\nComparison report saved to: {comparison_file}")
        return results


def main():
    parser = argparse.ArgumentParser(
        description='Select diverse compounds from virtual screening results'
    )
    parser.add_argument('singletons_csv', help='Path to singletons CSV file')
    parser.add_argument('non_singletons_csv', help='Path to non-singletons CSV file')
    parser.add_argument('-n', '--num-compounds', type=int, default=5,
                       help='Number of compounds to select (default: 5)')
    parser.add_argument('-s', '--strategy', default='balanced',
                       choices=['balanced', 'singleton_first', 'scaffold_first', 'pure_score'],
                       help='Selection strategy (default: balanced)')
    parser.add_argument('--max-similarity', type=float, default=0.7,
                       help='Maximum Tanimoto similarity between selected compounds (default: 0.7)')
    parser.add_argument('--compare', action='store_true',
                       help='Compare all strategies and export results')
    parser.add_argument('-o', '--output-dir', default='diversity_selection',
                       help='Output directory for results')
    
    args = parser.parse_args()
    
    # Initialize selector
    selector = DiversitySelector(
        args.singletons_csv,
        args.non_singletons_csv,
        output_dir=args.output_dir
    )
    
    # Print statistics
    print("\n" + "="*60)
    print("SCAFFOLD STATISTICS")
    print("="*60)
    stats = selector.get_scaffold_statistics()
    print(f"Total compounds: {stats['total_compounds']}")
    print(f"Singletons: {stats['singleton_count']} ({stats['singleton_percentage']:.1f}%)")
    print(f"Unique scaffolds: {stats['unique_scaffolds']}")
    print("\nScaffold size distribution:")
    for size_range, count in stats['scaffold_size_distribution'].items():
        print(f"  {size_range}: {count} scaffolds")
    
    if args.compare:
        # Compare all strategies
        print("\n" + "="*60)
        print("COMPARING SELECTION STRATEGIES")
        print("="*60)
        results = selector.compare_strategies(
            n_compounds=args.num_compounds,
            max_similarity=args.max_similarity
        )
    else:
        # Use specified strategy
        print("\n" + "="*60)
        print(f"SELECTING COMPOUNDS ({args.strategy} strategy)")
        print("="*60)
        
        selected = selector.select_diverse_compounds(
            n_compounds=args.num_compounds,
            strategy=args.strategy,
            max_similarity=args.max_similarity
        )
        
        print(f"\nSelected {len(selected)} compounds:")
        for idx, row in selected.iterrows():
            singleton_status = "singleton" if row.get('is_singleton') == 'True' else "non-singleton"
            print(f"  {idx+1}. {row['id']} (score: {row['docking_score']:.3f}, {singleton_status})")
        
        # Analyze and export
        analysis = selector.analyze_selection()
        print(f"\nDiversity metrics:")
        if 'diversity' in analysis:
            print(f"  Mean Tanimoto similarity: {analysis['diversity']['mean_similarity']:.3f}")
            print(f"  Max Tanimoto similarity: {analysis['diversity']['max_similarity']:.3f}")
        
        selector.export_results(prefix=f"{args.strategy}_{args.num_compounds}")
    
    print("\n" + "="*60)
    print("DONE! Check the output directory for results.")
    print("="*60)


if __name__ == "__main__":
    main()