#!/usr/bin/env python3
"""
Differentially Accessible Regions (DARs) and Differentially Expressed Genes (DEGs) Analysis Script

This script analyzes DARs (ATAC-seq) or DEGs (RNA-seq) data by:
1. Counting peaks/genes for each genotype in up/ and down/ folders
2. Creating bar plots showing peak/gene counts by genotype
3. Generating Upset plots for up/ and down/ regions

Date: 2025 (HS_ABA project)
"""

import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Literal
import logging
from collections import defaultdict
import matplotlib_venn as venn
import numpy as np
from upsetplot import UpSet, from_memberships
from itertools import combinations

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Type definitions
AnalysisType = Literal["DAR", "DEG"]


class DifferentialAnalyzer:
    """
    A class to analyze Differentially Accessible Regions (DARs) or Differentially Expressed Genes (DEGs) data.
    
    This class provides methods to:
    - Count peaks/genes for each genotype
    - Create visualizations including bar plots and Venn diagrams
    - Generate summary statistics
    """
    
    def __init__(self, data_folder: str, analysis_type: AnalysisType = "DAR"):
        """
        Initialize the DifferentialAnalyzer.
        
        Args:
            data_folder: Path to the folder containing up/ and down/ subfolders
            analysis_type: Type of analysis - "DAR" for ATAC-seq or "DEG" for RNA-seq
        """
        self.data_folder = Path(data_folder)
        self.analysis_type = analysis_type
        self.up_folder = self.data_folder / "up"
        self.down_folder = self.data_folder / "down"
        
        # Set terminology based on analysis type
        if analysis_type == "DAR":
            self.item_name = "peaks"
            self.item_name_singular = "peak"
            self.analysis_name = "DARs"
            self.full_name = "Differentially Accessible Regions"
        else:  # DEG
            self.item_name = "genes"
            self.item_name_singular = "gene"
            self.analysis_name = "DEGs"
            self.full_name = "Differentially Expressed Genes"
        
        # Validate folder structure
        if not self.data_folder.exists():
            raise FileNotFoundError(f"Data folder not found: {data_folder}")
        if not self.up_folder.exists():
            raise FileNotFoundError(f"Up folder not found: {self.up_folder}")
        if not self.down_folder.exists():
            raise FileNotFoundError(f"Down folder not found: {self.down_folder}")
        
        logger.info(f"Initialized {self.analysis_name} Analyzer with data folder: {data_folder}")
    
    def _parse_filename(self, filename: str) -> Tuple[str, str, str, str]:
        """
        Parse differential filename to extract genotype, treatment, and control information.
        
        Args:
            filename: Filename in format genotype_treatment_genotype_control_[u/d].txt
            
        Returns:
            Tuple of (genotype, treatment, control, direction)
        """
        # Remove .txt extension
        name = filename.replace('.txt', '')
        
        # Split by underscore
        parts = name.split('_')
        
        if len(parts) >= 4:
            genotype = parts[0]
            treatment = parts[1]
            control = parts[3]
            direction = parts[-1]  # 'u' for up, 'd' for down
            return genotype, treatment, control, direction
        else:
            logger.warning(f"Could not parse filename: {filename}")
            return "unknown", "unknown", "unknown", "unknown"
    
    def _read_items_from_file(self, file_path: Path) -> List[str]:
        """
        Read peak/gene regions from a differential file.
        
        Args:
            file_path: Path to the differential file
            
        Returns:
            List of peak/gene regions
        """
        try:
            with open(file_path, 'r') as f:
                lines = [line.strip() for line in f.readlines() if line.strip()]
            
            # Remove header line if it exists (usually 'x' or similar)
            if lines and lines[0] in ['x', 'X', 'peak_id', 'region', 'gene_id']:
                lines = lines[1:]
            
            # Return non-empty lines that contain peak/gene information
            items = [line for line in lines if line and not line.startswith('#')]
            
            return items
            
        except Exception as e:
            logger.error(f"Error reading file {file_path}: {e}")
            return []
    
    def _count_items_in_file(self, file_path: Path) -> int:
        """
        Count the number of peaks/genes in a differential file.
        
        Args:
            file_path: Path to the differential file
            
        Returns:
            Number of peaks/genes (non-empty lines minus header)
        """
        items = self._read_items_from_file(file_path)
        return len(items)
    
    def count_items_by_genotype(self) -> Dict[str, Dict[str, int]]:
        """
        Count peaks/genes for each genotype in both up/ and down/ folders.
        
        Returns:
            Dictionary with genotype as key and dict of {'up': count, 'down': count} as value
        """
        results = defaultdict(lambda: {'up': 0, 'down': 0})
        
        # Process up folder
        logger.info(f"Processing up folder for {self.item_name}...")
        for file_path in self.up_folder.glob("*.txt"):
            if file_path.name.startswith('.'):  # Skip hidden files
                continue
                
            genotype, treatment, control, direction = self._parse_filename(file_path.name)
            item_count = self._count_items_in_file(file_path)
            
            if genotype != "unknown":
                results[genotype]['up'] += item_count
                logger.info(f"Genotype {genotype}: {item_count} {self.item_name} in up folder")
        
        # Process down folder
        logger.info(f"Processing down folder for {self.item_name}...")
        for file_path in self.down_folder.glob("*.txt"):
            if file_path.name.startswith('.'):  # Skip hidden files
                continue
                
            genotype, treatment, control, direction = self._parse_filename(file_path.name)
            item_count = self._count_items_in_file(file_path)
            
            if genotype != "unknown":
                results[genotype]['down'] += item_count
                logger.info(f"Genotype {genotype}: {item_count} {self.item_name} in down folder")
        
        return dict(results)
    
    def create_bar_plot(self, item_counts: Dict[str, Dict[str, int]], 
                       output_path: str = None) -> None:
        """
        Create a bar plot showing peak/gene counts by genotype.
        
        Args:
            item_counts: Dictionary with genotype as key and dict of {'up': count, 'down': count} as value
            output_path: Path to save the plot (auto-generated if None)
        """
        if output_path is None:
            output_path = f"{self.analysis_name.lower()}_peak_counts.pdf"
        
        # Prepare data for plotting
        genotypes = list(item_counts.keys())
        up_counts = [item_counts[g]['up'] for g in genotypes]
        down_counts = [item_counts[g]['down'] for g in genotypes]
        
        # Create the plot
        fig, ax = plt.subplots(figsize=(12, 8))
        
        x = np.arange(len(genotypes))
        width = 0.35
        
        # Create bars
        up_bars = ax.bar(x - width/2, up_counts, width, label='Up-regulated', 
                        color='#2E8B57', alpha=0.8)
        down_bars = ax.bar(x + width/2, down_counts, width, label='Down-regulated', 
                          color='#DC143C', alpha=0.8)
        
        # Customize the plot
        ax.set_xlabel('Genotype', fontsize=14, fontweight='bold')
        ax.set_ylabel(f'Number of {self.item_name.capitalize()}', fontsize=14, fontweight='bold')
        ax.set_title(f'{self.analysis_name} {self.item_name.capitalize()} Counts by Genotype', fontsize=16, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(genotypes, fontsize=12)
        ax.legend(fontsize=12)
        ax.grid(True, alpha=0.3)
        
        # Add value labels on bars
        for bar in up_bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                   f'{int(height)}', ha='center', va='bottom', fontsize=10)
        
        for bar in down_bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                   f'{int(height)}', ha='center', va='bottom', fontsize=10)
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"Bar plot saved to: {output_path}")
    
    def _read_samples_from_files(self, folder_path: Path) -> Dict[str, set]:
        """
        Read all unique samples from files in a folder.
        
        Args:
            folder_path: Path to the folder containing sample files
            
        Returns:
            Dictionary with genotype as key and set of samples as value
        """
        samples_dict = {}
        
        for file_path in folder_path.glob("*.txt"):
            if file_path.name.startswith('.'):
                continue
                
            genotype, treatment, control, direction = self._parse_filename(file_path.name)
            if genotype != "unknown":
                samples = self._read_items_from_file(file_path)
                if genotype not in samples_dict:
                    samples_dict[genotype] = set()
                samples_dict[genotype].update(samples)
        
        return samples_dict
    
    def _create_boolean_matrix(self, samples_dict: Dict[str, set]) -> Tuple[pd.DataFrame, List[str]]:
        """
        Create a boolean matrix from samples dictionary.
        
        Args:
            samples_dict: Dictionary with genotype as key and set of samples as value
            
        Returns:
            Tuple of (boolean DataFrame, list of all unique samples)
        """
        # Get all unique samples across all genotypes
        all_samples = set()
        for samples in samples_dict.values():
            all_samples.update(samples)
        
        all_samples_list = sorted(list(all_samples))
        
        # Create boolean matrix
        data = []
        for sample in all_samples_list:
            row = {'sample': sample}
            for genotype in samples_dict.keys():
                # Ensure boolean values are properly set
                row[genotype] = sample in samples_dict[genotype]
            data.append(row)
        
        df = pd.DataFrame(data)
        df.set_index('sample', inplace=True)
        
        # Convert to boolean type to ensure consistency
        for col in df.columns:
            df[col] = df[col].astype(bool)
        
        return df, all_samples_list
    
    def _validate_upset_consistency(self, df: pd.DataFrame, samples_dict: Dict[str, set]) -> Dict[str, int]:
        """
        Validate that Upset plot counts match exported results.
        
        Args:
            df: Boolean DataFrame with samples as rows and genotypes as columns
            samples_dict: Dictionary with genotype as key and set of samples as value
            
        Returns:
            Dictionary with combination as key and count as value for validation
        """
        from itertools import combinations
        
        validation_counts = {}
        genotypes = list(samples_dict.keys())
        
        # Calculate counts using the SAME logic as Upset plot
        # For each sample, determine which combination it belongs to
        sample_combinations = {}
        
        for sample in df.index:
            # Get all genotypes this sample belongs to
            sample_genotypes = [name for name in genotypes if df.loc[sample, name]]
            
            if sample_genotypes:
                # Sort to ensure consistent naming
                combo_name = ' & '.join(sorted(sample_genotypes))
                sample_combinations[combo_name] = sample_combinations.get(combo_name, 0) + 1
        
        # Fill in all possible combinations (including those with 0 count)
        for r in range(1, len(genotypes) + 1):
            for combo in combinations(genotypes, r):
                combo_name = ' & '.join(sorted(combo))
                validation_counts[combo_name] = sample_combinations.get(combo_name, 0)
        
        return validation_counts
    
    def _export_upset_results(self, df: pd.DataFrame, samples_dict: Dict[str, set], 
                             output_path: Path, direction: str) -> None:
        """
        Export Upset plot results to CSV and TXT files.
        
        Args:
            df: Boolean DataFrame with samples as rows and genotypes as columns
            samples_dict: Dictionary with genotype as key and set of samples as value
            output_path: Base path for output files
            direction: 'up' or 'down' for file naming
        """
        # Export boolean matrix to CSV
        csv_path = output_path.parent / f"upset_{direction}_boolean_matrix.csv"
        df.to_csv(csv_path)
        logger.info(f"Boolean matrix saved to: {csv_path}")
        
        # Calculate and export intersection statistics using the same logic as Upset plot
        genotypes = list(samples_dict.keys())
        intersection_stats = []
        
        # Calculate all possible combinations using the SAME logic as Upset plot
        from itertools import combinations
        
        # First, calculate the actual counts using Upset plot logic
        sample_combinations = {}
        for sample in df.index:
            # Get all genotypes this sample belongs to
            sample_genotypes = [name for name in genotypes if df.loc[sample, name]]
            
            if sample_genotypes:
                # Sort to ensure consistent naming
                combo_name = ' & '.join(sorted(sample_genotypes))
                sample_combinations[combo_name] = sample_combinations.get(combo_name, 0) + 1
        
        # Now create statistics for all possible combinations
        for r in range(1, len(genotypes) + 1):
            for combo in combinations(genotypes, r):
                combo_name = ' & '.join(sorted(combo))
                count = sample_combinations.get(combo_name, 0)
                
                intersection_stats.append({
                    'combination': combo_name,
                    'count': count,
                    'percentage': (count / len(df)) * 100 if len(df) > 0 else 0
                })
        
        # Export intersection statistics to CSV
        stats_df = pd.DataFrame(intersection_stats)
        stats_csv_path = output_path.parent / f"upset_{direction}_intersection_stats.csv"
        stats_df.to_csv(stats_csv_path, index=False)
        logger.info(f"Intersection statistics saved to: {stats_csv_path}")
        
        # Export individual combination results to TXT files
        # For single genotypes, use the Upset plot logic to ensure consistency
        for combo in combinations(genotypes, 1):  # Single genotypes
            genotype = combo[0]
            # Get samples that belong to EXACTLY this genotype (Upset plot logic)
            samples = []
            for sample in df.index:
                sample_genotypes = [name for name in genotypes if df.loc[sample, name]]
                if sample_genotypes == [genotype]:  # Only this genotype
                    samples.append(sample)
            
            txt_path = output_path.parent / f"upset_{direction}_{genotype}_only.txt"
            with open(txt_path, 'w') as f:
                f.write(f"# {genotype} only {self.item_name}\n")
                f.write(f"# Count: {len(samples)}\n")
                for sample in sorted(samples):
                    f.write(f"{sample}\n")
            logger.info(f"{genotype} only {self.item_name} saved to: {txt_path}")
        
        # Export intersection combinations to TXT files
        # Use the SAME logic as Upset plot to ensure consistency
        for r in range(2, len(genotypes) + 1):
            for combo in combinations(genotypes, r):
                # Get samples that belong to EXACTLY this combination (Upset plot logic)
                intersection_samples = []
                combo_set = set(combo)
                
                for sample in df.index:
                    sample_genotypes = [name for name in genotypes if df.loc[sample, name]]
                    sample_genotypes_set = set(sample_genotypes)
                    
                    # Sample must belong to EXACTLY this combination
                    if sample_genotypes_set == combo_set:
                        intersection_samples.append(sample)
                
                if intersection_samples:
                    combo_name = '_'.join(sorted(combo))
                    txt_path = output_path.parent / f"upset_{direction}_{combo_name}_intersection.txt"
                    with open(txt_path, 'w') as f:
                        f.write(f"# {combo_name} combination (exact match)\n")
                        f.write(f"# Count: {len(intersection_samples)}\n")
                        for sample in sorted(intersection_samples):
                            f.write(f"{sample}\n")
                    logger.info(f"Combination {combo_name} saved to: {txt_path}")
                else:
                    # Log when no exact match is found
                    combo_name = '_'.join(sorted(combo))
                    logger.info(f"No exact match found for {combo_name}")
        
        # Validate consistency between exported results and Upset plot logic
        validation_counts = self._validate_upset_consistency(df, samples_dict)
        logger.info(f"Validation counts for {direction} direction:")
        for combo, count in validation_counts.items():
            logger.info(f"  {combo}: {count} {self.item_name}")
        
        # Export validation summary
        validation_df = pd.DataFrame(list(validation_counts.items()), columns=['Combination', 'Count'])
        validation_csv_path = output_path.parent / f"upset_{direction}_validation_summary.csv"
        validation_df.to_csv(validation_csv_path, index=False)
        logger.info(f"Validation summary saved to: {validation_csv_path}")
    
    def _create_upset_plot(self, samples_dict: Dict[str, set], title: str, 
                          output_path: Path, direction: str) -> None:
        """
        Create an Upset plot from samples data.
        
        Args:
            samples_dict: Dictionary with genotype as key and set of samples as value
            title: Title for the plot
            output_path: Path to save the plot
            direction: 'up' or 'down' for file naming
        """
        # Filter out empty sets
        non_empty_sets = {k: v for k, v in samples_dict.items() if len(v) > 0}
        
        if len(non_empty_sets) == 0:
            # If no samples, create a simple message plot
            fig, ax = plt.subplots(figsize=(10, 8))
            ax.text(0.5, 0.5, f'No {self.item_name} found', ha='center', va='center', 
                   transform=ax.transAxes, fontsize=16)
            ax.set_title(title, fontsize=16, fontweight='bold')
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            plt.tight_layout()
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
            return
        
        if len(non_empty_sets) == 1:
            # If only one genotype has samples, create a simple bar plot
            fig, ax = plt.subplots(figsize=(10, 8))
            genotype = list(non_empty_sets.keys())[0]
            count = len(non_empty_sets[genotype])
            ax.bar([genotype], [count], color='skyblue', alpha=0.7)
            ax.set_xlabel('Genotype')
            ax.set_ylabel(f'Number of {self.item_name.capitalize()}')
            ax.set_title(title)
            ax.text(0, count + count*0.1, f'{count}', ha='center', va='bottom')
            plt.tight_layout()
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
            return
        
        # Create boolean matrix
        df, all_samples = self._create_boolean_matrix(non_empty_sets)
        
        # Export results to CSV and TXT files
        self._export_upset_results(df, non_empty_sets, output_path, direction)
        
        # Create Upset plot for 2 or more genotypes
        try:
            # Use the same all_samples from boolean matrix creation for consistency
            # Build memberships list using the boolean matrix
            memberships = []
            for sample in all_samples:
                # Get membership based on boolean matrix to ensure consistency
                membership = [name for name in non_empty_sets.keys() if df.loc[sample, name]]
                memberships.append(membership)
            
            # Create UpSet data using from_memberships
            from upsetplot import from_memberships
            upset_data = from_memberships(memberships)
            
            # Create the plot
            fig = plt.figure(figsize=(16, 10))
            
            # Create UpSet object and plot with explicit subset_size and sorting
            upset = UpSet(upset_data, show_counts=True, sort_categories_by='cardinality', subset_size='count')
            upset.plot(fig=fig)
            
            # Add title
            plt.suptitle(title, fontsize=16, fontweight='bold')
            plt.tight_layout()
            plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
            plt.close()
            
            # Log success and validate consistency
            logger.info(f"Upset plot created successfully with {len(all_samples)} {self.item_name} and {len(non_empty_sets)} genotypes")
            
            # Additional validation: compare Upset plot counts with our validation counts
            validation_counts = self._validate_upset_consistency(df, non_empty_sets)
            logger.info("Upset plot validation - combination counts:")
            for combo, count in validation_counts.items():
                logger.info(f"  {combo}: {count} {self.item_name}")
            
        except Exception as e:
            logger.error(f"Error creating Upset plot: {e}")
            import traceback
            traceback.print_exc()
            
            # Fallback to simple bar plot
            fig, ax = plt.subplots(figsize=(12, 8))
            genotypes = list(non_empty_sets.keys())
            counts = [len(non_empty_sets[g]) for g in genotypes]
            
            bars = ax.bar(genotypes, counts, color='skyblue', alpha=0.7)
            ax.set_xlabel('Genotype', fontsize=12)
            ax.set_ylabel(f'Number of {self.item_name.capitalize()}', fontsize=12)
            ax.set_title(title, fontsize=14, fontweight='bold')
            
            # Add value labels on bars
            for bar, count in zip(bars, counts):
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2., height + max(counts)*0.01,
                       f'{count}', ha='center', va='bottom', fontsize=10)
            
            plt.tight_layout()
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()
    
    def create_upset_plots(self, item_counts: Dict[str, Dict[str, int]], 
                          output_dir: str = ".") -> None:
        """
        Create Upset plots for up/ and down/ regions.
        
        Args:
            item_counts: Dictionary with genotype as key and dict of {'up': count, 'down': count} as value
            output_dir: Directory to save the Upset plots
        """
        output_path = Path(output_dir)
        
        # Process up-regulated samples
        logger.info(f"Processing up-regulated {self.item_name} for Upset plot...")
        up_samples = self._read_samples_from_files(self.up_folder)
        if up_samples:
            self._create_upset_plot(up_samples, f"Up-regulated {self.analysis_name} by Genotype", 
                                  output_path / f"{self.analysis_name.lower()}_up_upset.pdf", "up")
            logger.info(f"Up-regulated Upset plot saved")
        else:
            logger.warning(f"No up-regulated {self.item_name} found")
        
        # Process down-regulated samples
        logger.info(f"Processing down-regulated {self.item_name} for Upset plot...")
        down_samples = self._read_samples_from_files(self.down_folder)
        if down_samples:
            self._create_upset_plot(down_samples, f"Down-regulated {self.analysis_name} by Genotype", 
                                  output_path / f"{self.analysis_name.lower()}_down_upset.pdf", "down")
            logger.info(f"Down-regulated Upset plot saved")
        else:
            logger.warning(f"No down-regulated {self.item_name} found")
    
    def create_summary_table(self, item_counts: Dict[str, Dict[str, int]], 
                           output_path: str = None) -> None:
        """
        Create a summary table with peak/gene counts.
        
        Args:
            item_counts: Dictionary with genotype as key and dict of {'up': count, 'down': count} as value
            output_path: Path to save the summary table (auto-generated if None)
        """
        if output_path is None:
            output_path = f"{self.analysis_name.lower()}_summary.csv"
        
        # Convert to DataFrame
        df = pd.DataFrame.from_dict(item_counts, orient='index')
        df.index.name = 'Genotype'
        df.columns = ['Up-regulated', 'Down-regulated']
        df['Total'] = df['Up-regulated'] + df['Down-regulated']
        
        # Save to CSV
        df.to_csv(output_path)
        logger.info(f"Summary table saved to: {output_path}")
        
        # Print summary
        print("\n" + "="*50)
        print(f"{self.analysis_name} ANALYSIS SUMMARY")
        print("="*50)
        print(df.to_string())
        print("="*50)
    
    def run_analysis(self, output_dir: str = "output") -> None:
        """
        Run the complete differential analysis.
        
        Args:
            output_dir: Directory to save all output files
        """
        # Create output directory
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)
        
        logger.info(f"Starting {self.analysis_name} analysis...")
        
        # Count items
        item_counts = self.count_items_by_genotype()
        
        # Create visualizations
        self.create_bar_plot(item_counts, output_path / f"{self.analysis_name.lower()}_peak_counts.pdf")
        self.create_upset_plots(item_counts, output_path)
        self.create_summary_table(item_counts, output_path / f"{self.analysis_name.lower()}_summary.csv")
        
        # Final validation: ensure consistency across all outputs
        self._final_validation_check(output_path)
        
        logger.info(f"{self.analysis_name} analysis completed successfully!")
    
    def _final_validation_check(self, output_path: Path) -> None:
        """
        Perform final validation to ensure consistency across all outputs.
        
        Args:
            output_path: Path to the output directory
        """
        logger.info("Performing final validation check...")
        
        # Check if validation files exist and compare with summary
        validation_files = list(output_path.glob("upset_*_validation_summary.csv"))
        
        if validation_files:
            logger.info("Validation summary files found:")
            for file_path in validation_files:
                try:
                    df = pd.read_csv(file_path)
                    direction = file_path.stem.split('_')[1]  # Extract 'up' or 'down'
                    logger.info(f"  {direction} direction validation:")
                    for _, row in df.iterrows():
                        logger.info(f"    {row['Combination']}: {row['Count']} {self.item_name}")
                except Exception as e:
                    logger.error(f"Error reading validation file {file_path}: {e}")
        else:
            logger.warning("No validation summary files found")
        
        logger.info("Final validation check completed")


def main():
    """
    Main function to run the differential analysis.
    """
    import argparse
    
    parser = argparse.ArgumentParser(description="Analyze DARs or DEGs data")
    parser.add_argument("data_folder", help="Path to the folder containing up/ and down/ subfolders")
    parser.add_argument("--output", "-o", default="output", help="Output directory (default: output)")
    parser.add_argument("--type", "-t", choices=["DAR", "DEG"], default="DAR", 
                       help="Analysis type: DAR for ATAC-seq or DEG for RNA-seq (default: DAR)")
    
    args = parser.parse_args()
    
    try:
        # Initialize analyzer
        analyzer = DifferentialAnalyzer(args.data_folder, args.type)
        
        # Run analysis
        analyzer.run_analysis(args.output)
        
    except Exception as e:
        logger.error(f"Error during analysis: {e}")
        raise


if __name__ == "__main__":
    main() 