#!/usr/bin/env python3
"""
Peaks Annotation Analysis Script (Fixed Version)

This script processes peak files and performs genomic annotation using R packages
GenomicFeatures and ChIPseeker. It provides comprehensive analysis including:
- Peak annotation with genomic features
- Distance to TSS analysis
- Proximal/distal peak classification
- Gene extraction and visualization

Author: Generated for HS_ABA analysis
Date: 2025
Fixed: Updated for modern rpy2 API compatibility
"""

import os
import sys
import glob
import logging
import argparse
import warnings
from pathlib import Path
from typing import List, Dict, Tuple, Optional
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime

# Suppress deprecation warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

# R integration with modern API
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri, Formula
from rpy2.robjects.packages import importr
from rpy2.robjects.conversion import localconverter

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class PeaksAnnotationAnalyzer:
    """
    A comprehensive analyzer for peaks annotation using R GenomicFeatures and ChIPseeker.
    
    This class handles the complete workflow of peaks annotation including:
    - Reading peak files in various formats
    - Creating genomic annotations from GFF files
    - Annotating peaks with genomic features
    - Classifying peaks as proximal or distal
    - Extracting associated genes
    - Generating visualizations
    """
    
    def __init__(self, output_dir: str = "output/peaks_annotation"):
        """
        Initialize the PeaksAnnotationAnalyzer.
        
        Args:
            output_dir: Directory to save output files
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize R packages
        self._setup_r_environment()
        
        # Store analysis results
        self.annotation_results: Dict = {}
        self.gene_lists: Dict = {}
        
    def _setup_r_environment(self) -> None:
        """Setup R environment and import required packages."""
        try:
            logger.info("Setting up R environment...")
            
            # Test basic R functionality first
            robjects.r('print("R is working")')
            
            # Import required R packages
            logger.info("Importing R packages...")
            self.genomic_features = importr('GenomicFeatures')
            self.chipseeker = importr('ChIPseeker')
            self.genomic_ranges = importr('GenomicRanges')
            self.base_r = importr('base')
            self.utils_r = importr('utils')
            
            logger.info("Successfully imported R packages")
            
        except Exception as e:
            logger.error(f"Failed to import R packages: {e}")
            logger.error("Please ensure R and required packages are installed:")
            logger.error("In R console, run:")
            logger.error("install.packages(c('BiocManager'))")
            logger.error("BiocManager::install(c('GenomicFeatures', 'ChIPseeker', 'GenomicRanges'))")
            sys.exit(1)
    
    def read_peaks_files(self, peaks_dir: str) -> Dict[str, pd.DataFrame]:
        """
        Read all peak files from the specified directory.
        
        Args:
            peaks_dir: Directory containing peak files
            
        Returns:
            Dictionary mapping file names to peak dataframes
        """
        peaks_data = {}
        peaks_path = Path(peaks_dir)
        
        if not peaks_path.exists():
            raise FileNotFoundError(f"Peaks directory not found: {peaks_dir}")
        
        # Look for common peak file extensions
        peak_files = []
        for pattern in ['*.narrowPeak', '*.bed', '*.txt']:
            peak_files.extend(glob.glob(str(peaks_path / pattern)))
        
        if not peak_files:
            raise ValueError(f"No peak files found in {peaks_dir}")
        
        logger.info(f"Found {len(peak_files)} peak files")
        
        for peak_file in peak_files:
            file_name = Path(peak_file).stem
            logger.info(f"Reading peak file: {file_name}")
            
            try:
                # Read peak file - assuming tab-separated format
                df = pd.read_csv(peak_file, sep='\t', header=0)
                
                # Handle different peak file formats
                if df.shape[1] >= 3:
                    # Standard format: chr, start, end, [additional columns]
                    peaks_df = pd.DataFrame({
                        'Chr': df.iloc[:, 0].astype(str),
                        'Start': df.iloc[:, 1].astype(int),
                        'End': df.iloc[:, 2].astype(int),
                        'Peak': [f"{file_name}_peak_{i+1}" for i in range(len(df))]
                    })
                elif df.shape[1] == 1 and 'Peaks' in df.columns:
                    # Cluster peaks format: single column with format "chr.start.end"
                    logger.info(f"Detected cluster peaks format in {file_name}")
                    peaks_list = []
                    
                    for peak_id in df['Peaks']:
                        try:
                            # Parse format like "chr2.8816585.8816956"
                            parts = str(peak_id).split('.')
                            if len(parts) >= 3:
                                chr_name = parts[0]
                                start = int(parts[1])
                                end = int(parts[2])
                                peaks_list.append({
                                    'Chr': chr_name,
                                    'Start': start,
                                    'End': end,
                                    'Peak': peak_id
                                })
                        except (ValueError, IndexError) as e:
                            logger.warning(f"Skipping invalid peak ID: {peak_id} ({e})")
                            continue
                    
                    if not peaks_list:
                        raise ValueError(f"No valid peaks found in cluster format: {peak_file}")
                    
                    peaks_df = pd.DataFrame(peaks_list)
                else:
                    raise ValueError(f"Invalid peak file format: {peak_file}")
                
                # Validate data
                if peaks_df.empty:
                    logger.warning(f"Empty peak file: {file_name}")
                    continue
                
                # Ensure proper chromosome naming (keep chr prefix to match GFF format)
                # peaks_df['Chr'] = peaks_df['Chr'].str.replace('chr', '', regex=False)
                
                peaks_data[file_name] = peaks_df
                logger.info(f"Successfully loaded {len(peaks_df)} peaks from {file_name}")
                
            except Exception as e:
                logger.error(f"Error reading peak file {peak_file}: {e}")
                continue
        
        if not peaks_data:
            raise ValueError("No valid peak files could be loaded")
        
        return peaks_data
    
    def create_txdb_from_gff(self, gff_file: str) -> robjects.RObject:
        """
        Create TxDb object from GFF file using GenomicFeatures.
        
        Args:
            gff_file: Path to GFF annotation file
            
        Returns:
            R TxDb object
        """
        if not Path(gff_file).exists():
            raise FileNotFoundError(f"GFF file not found: {gff_file}")
        
        logger.info(f"Creating TxDb from GFF file: {gff_file}")
        
        try:
            # Create TxDb from GFF file
            txdb = self.genomic_features.makeTxDbFromGFF(
                file=gff_file,
                format="gff3"
            )
            
            logger.info("Successfully created TxDb object")
            return txdb
            
        except Exception as e:
            logger.error(f"Error creating TxDb from GFF: {e}")
            raise
    
    def convert_peaks_to_granges(self, peaks_df: pd.DataFrame) -> robjects.RObject:
        """
        Convert peaks dataframe to R GRanges object using modern rpy2 API.
        
        Args:
            peaks_df: Pandas dataframe with Chr, Start, End, Peak columns
            
        Returns:
            R GRanges object
        """
        # Validate required columns
        required_cols = ['Chr', 'Start', 'End', 'Peak']
        if not all(col in peaks_df.columns for col in required_cols):
            raise ValueError(f"Peaks dataframe must contain columns: {required_cols}")
        
        try:
            # Convert to R dataframe using modern API
            with localconverter(robjects.default_converter + pandas2ri.converter):
                r_peaks_df = robjects.conversion.py2rpy(peaks_df)
            
            # Create GRanges object
            granges = self.genomic_ranges.makeGRangesFromDataFrame(
                r_peaks_df,
                keep_extra_columns=True
            )
            
            return granges
            
        except Exception as e:
            logger.error(f"Error converting peaks to GRanges: {e}")
            raise
    
    def annotate_peaks(self, granges: robjects.RObject, txdb: robjects.RObject) -> robjects.RObject:
        """
        Annotate peaks using ChIPseeker.
        
        Args:
            granges: R GRanges object containing peaks
            txdb: R TxDb object for annotation
            
        Returns:
            R annotated peak object
        """
        logger.info("Annotating peaks with ChIPseeker")
        
        try:
            # Define TSS region and annotation priority
            tss_region = robjects.IntVector([-1500, 1500])
            annotation_priority = robjects.StrVector([
                "Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"
            ])
            
            # Annotate peaks
            peak_anno = self.chipseeker.annotatePeak(
                peak=granges,
                TxDb=txdb,
                tssRegion=tss_region,
                genomicAnnotationPriority=annotation_priority
            )
            
            logger.info("Successfully annotated peaks")
            return peak_anno
            
        except Exception as e:
            logger.error(f"Error annotating peaks: {e}")
            raise
    
    def process_annotation_results(self, peak_anno: robjects.RObject, sample_name: str) -> pd.DataFrame:
        """
        Process peak annotation results and add distance transformation.
        
        Args:
            peak_anno: R annotated peak object
            sample_name: Name of the sample
            
        Returns:
            Processed annotation dataframe
        """
        logger.info(f"Processing annotation results for {sample_name}")
        
        try:
            # Extract annotation data
            anno_data = robjects.r('function(x) as.data.frame(x@anno)')(peak_anno)
            
            # Convert to pandas dataframe using modern API
            with localconverter(robjects.default_converter + pandas2ri.converter):
                anno_df = robjects.conversion.rpy2py(anno_data)
            
            # Add log10 distance to TSS
            if 'distanceToTSS' in anno_df.columns:
                anno_df['distanceToTSS_log10'] = np.log10(
                    (np.abs(anno_df['distanceToTSS']) / 1000) + 1
                )
            else:
                logger.warning("distanceToTSS column not found in annotation results")
                anno_df['distanceToTSS_log10'] = np.nan
            
            # Classify peaks as proximal or distal
            if 'annotation' in anno_df.columns:
                anno_df['peak_type'] = anno_df['annotation'].apply(
                    lambda x: 'distal' if 'Distal Intergenic' in str(x) else 'proximal'
                )
            else:
                logger.warning("annotation column not found in annotation results")
                anno_df['peak_type'] = 'unknown'
            
            # Add sample information
            anno_df['sample'] = sample_name
            
            logger.info(f"Processed {len(anno_df)} annotated peaks for {sample_name}")
            return anno_df
            
        except Exception as e:
            logger.error(f"Error processing annotation results: {e}")
            raise
    
    def extract_genes(self, anno_df: pd.DataFrame, sample_name: str) -> Dict[str, List[str]]:
        """
        Extract genes associated with all, distal, and proximal peaks.
        
        Args:
            anno_df: Processed annotation dataframe
            sample_name: Name of the sample
            
        Returns:
            Dictionary containing gene lists for different categories
        """
        logger.info(f"Extracting genes for {sample_name}")
        
        gene_lists = {}
        
        if 'geneId' not in anno_df.columns:
            logger.warning("geneId column not found in annotation results")
            return gene_lists
        
        # Remove NA values and get unique genes
        all_genes = anno_df['geneId'].dropna().unique().tolist()
        
        # Extract genes for different peak types
        distal_genes = anno_df[anno_df['peak_type'] == 'distal']['geneId'].dropna().unique().tolist()
        proximal_genes = anno_df[anno_df['peak_type'] == 'proximal']['geneId'].dropna().unique().tolist()
        
        gene_lists = {
            'all_peaks': all_genes,
            'distal_peaks': distal_genes,
            'proximal_peaks': proximal_genes
        }
        
        # Log statistics
        logger.info(f"Extracted genes - All: {len(all_genes)}, "
                   f"Distal: {len(distal_genes)}, Proximal: {len(proximal_genes)}")
        
        return gene_lists
    
    def save_gene_lists(self, gene_lists: Dict[str, List[str]], sample_name: str) -> None:
        """
        Save gene lists to separate files.
        
        Args:
            gene_lists: Dictionary containing gene lists
            sample_name: Name of the sample
        """
        for category, genes in gene_lists.items():
            if genes:
                output_file = self.output_dir / f"{sample_name}_{category}_genes.txt"
                
                with open(output_file, 'w') as f:
                    f.write(f"# Genes associated with {category} in {sample_name}\n")
                    f.write(f"# Total genes: {len(genes)}\n")
                    f.write("geneId\n")
                    for gene in sorted(genes):
                        f.write(f"{gene}\n")
                
                logger.info(f"Saved {len(genes)} genes to {output_file}")
    
    def _simplify_annotation(self, annotation: str) -> str:
        """
        Simplify annotation categories to 8 main types.
        
        Args:
            annotation: Original annotation string
            
        Returns:
            Simplified annotation category
        """
        annotation = str(annotation).lower()
        
        if 'promoter' in annotation:
            return 'Promoter'
        elif '5' in annotation and 'utr' in annotation:
            return '5UTR'
        elif '3' in annotation and 'utr' in annotation:
            return '3UTR'
        elif 'exon' in annotation:
            return 'Exon'
        elif 'intron' in annotation:
            return 'Intron'
        elif 'downstream' in annotation:
            return 'Downstream'
        elif 'distal intergenic' in annotation:
            return 'Distal intergenic'
        elif 'intergenic' in annotation:
            return 'Intergenic'
        else:
            return 'Other'

    def create_visualizations(self, all_results: Dict[str, pd.DataFrame]) -> None:
        """
        Create PDF visualizations for distance to TSS and annotation percentages.
        
        Args:
            all_results: Dictionary mapping sample names to annotation dataframes
        """
        logger.info("Creating visualizations")
        
        # Combine all results
        combined_df = pd.concat(all_results.values(), ignore_index=True)
        
        # Add simplified annotation column
        if 'annotation' in combined_df.columns:
            combined_df['simplified_annotation'] = combined_df['annotation'].apply(self._simplify_annotation)
        
        # Create figure with subplots
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        fig.suptitle('Peaks Annotation Analysis', fontsize=16, fontweight='bold')
        
        # 1. Distance to TSS (log10) distribution
        ax1 = axes[0, 0]
        if 'distanceToTSS_log10' in combined_df.columns:
            combined_df.boxplot(column='distanceToTSS_log10', by='sample', ax=ax1)
            ax1.set_title('Distance to TSS (log10) Distribution')
            ax1.set_xlabel('Sample')
            ax1.set_ylabel('Distance to TSS (log10)')
            ax1.tick_params(axis='x', rotation=45)
        
        # 2. Distance to TSS histogram
        ax2 = axes[0, 1]
        if 'distanceToTSS_log10' in combined_df.columns:
            for sample in combined_df['sample'].unique():
                sample_data = combined_df[combined_df['sample'] == sample]['distanceToTSS_log10'].dropna()
                ax2.hist(sample_data, alpha=0.6, label=sample, bins=30)
            ax2.set_title('Distance to TSS (log10) Histogram')
            ax2.set_xlabel('Distance to TSS (log10)')
            ax2.set_ylabel('Frequency')
            ax2.legend()
        
        # 3. Simplified annotation percentages by sample
        ax3 = axes[1, 0]
        if 'simplified_annotation' in combined_df.columns:
            # Calculate annotation percentages using simplified categories
            annotation_counts = combined_df.groupby(['sample', 'simplified_annotation']).size().unstack(fill_value=0)
            annotation_percentages = annotation_counts.div(annotation_counts.sum(axis=1), axis=0) * 100
            
            # Define consistent order for annotation categories
            category_order = ['Promoter', '5UTR', '3UTR', 'Exon', 'Intron', 'Downstream', 'Intergenic', 'Distal intergenic']
            
            # Reorder columns to match desired order (only include categories that exist)
            existing_categories = [cat for cat in category_order if cat in annotation_percentages.columns]
            if existing_categories:
                annotation_percentages = annotation_percentages[existing_categories]
            
            # Create stacked bar plot with consistent colors
            colors = plt.cm.Set3(np.linspace(0, 1, len(annotation_percentages.columns)))
            annotation_percentages.plot(kind='bar', stacked=True, ax=ax3, color=colors)
            ax3.set_title('Annotation Percentages by Sample (Simplified)')
            ax3.set_xlabel('Sample')
            ax3.set_ylabel('Percentage (%)')
            ax3.tick_params(axis='x', rotation=45)
            ax3.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        
        # 4. Peak type distribution
        ax4 = axes[1, 1]
        if 'peak_type' in combined_df.columns:
            peak_type_counts = combined_df.groupby(['sample', 'peak_type']).size().unstack(fill_value=0)
            peak_type_percentages = peak_type_counts.div(peak_type_counts.sum(axis=1), axis=0) * 100
            
            peak_type_percentages.plot(kind='bar', ax=ax4, color=['skyblue', 'lightcoral'])
            ax4.set_title('Peak Type Distribution (Proximal vs Distal)')
            ax4.set_xlabel('Sample')
            ax4.set_ylabel('Percentage (%)')
            ax4.tick_params(axis='x', rotation=45)
            ax4.legend()
        
        # Adjust layout and save
        plt.tight_layout()
        
        # Save as PDF
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        pdf_file = self.output_dir / f"peaks_annotation_analysis_{timestamp}.pdf"
        plt.savefig(pdf_file, format='pdf', dpi=300, bbox_inches='tight')
        
        # Also save as PNG for easy viewing
        png_file = self.output_dir / f"peaks_annotation_analysis_{timestamp}.png"
        plt.savefig(png_file, format='png', dpi=300, bbox_inches='tight')
        
        plt.close()
        
        logger.info(f"Saved visualizations to {pdf_file} and {png_file}")
    
    def run_analysis(self, peaks_dir: str, gff_file: str) -> Dict[str, pd.DataFrame]:
        """
        Run complete peaks annotation analysis.
        
        Args:
            peaks_dir: Directory containing peak files
            gff_file: Path to GFF annotation file
            
        Returns:
            Dictionary mapping sample names to annotation results
        """
        logger.info("Starting peaks annotation analysis")
        
        # Read peak files
        peaks_data = self.read_peaks_files(peaks_dir)
        
        # Create TxDb from GFF
        txdb = self.create_txdb_from_gff(gff_file)
        
        all_results = {}
        all_gene_lists = {}
        
        # Process each peak file
        for sample_name, peaks_df in peaks_data.items():
            logger.info(f"Processing sample: {sample_name}")
            
            try:
                # Convert to GRanges
                granges = self.convert_peaks_to_granges(peaks_df)
                
                # Annotate peaks
                peak_anno = self.annotate_peaks(granges, txdb)
                
                # Process results
                anno_df = self.process_annotation_results(peak_anno, sample_name)
                all_results[sample_name] = anno_df
                
                # Extract genes
                gene_lists = self.extract_genes(anno_df, sample_name)
                all_gene_lists[sample_name] = gene_lists
                
                # Save gene lists
                self.save_gene_lists(gene_lists, sample_name)
                
                # Save annotation results
                anno_file = self.output_dir / f"{sample_name}_annotation_results.csv"
                anno_df.to_csv(anno_file, index=False)
                logger.info(f"Saved annotation results to {anno_file}")
                
            except Exception as e:
                logger.error(f"Error processing sample {sample_name}: {e}")
                continue
        
        # Create visualizations
        if all_results:
            self.create_visualizations(all_results)
        
        # Save summary
        self._save_analysis_summary(all_results, all_gene_lists)
        
        logger.info("Peaks annotation analysis completed successfully")
        return all_results
    
    def _save_analysis_summary(self, all_results: Dict[str, pd.DataFrame], 
                              all_gene_lists: Dict[str, Dict[str, List[str]]]) -> None:
        """
        Save analysis summary to file.
        
        Args:
            all_results: Dictionary of annotation results
            all_gene_lists: Dictionary of gene lists
        """
        summary_file = self.output_dir / "analysis_summary.txt"
        
        with open(summary_file, 'w') as f:
            f.write("Peaks Annotation Analysis Summary\n")
            f.write("=" * 40 + "\n\n")
            f.write(f"Analysis completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Total samples processed: {len(all_results)}\n\n")
            
            for sample_name, anno_df in all_results.items():
                f.write(f"Sample: {sample_name}\n")
                f.write(f"  Total peaks: {len(anno_df)}\n")
                
                if 'peak_type' in anno_df.columns:
                    proximal_count = len(anno_df[anno_df['peak_type'] == 'proximal'])
                    distal_count = len(anno_df[anno_df['peak_type'] == 'distal'])
                    f.write(f"  Proximal peaks: {proximal_count}\n")
                    f.write(f"  Distal peaks: {distal_count}\n")
                
                # Add simplified annotation statistics
                if 'annotation' in anno_df.columns:
                    # Apply simplification to the dataframe
                    simplified_annotations = anno_df['annotation'].apply(self._simplify_annotation)
                    annotation_counts = simplified_annotations.value_counts()
                    
                    f.write("  Simplified annotation distribution:\n")
                    for annotation, count in annotation_counts.items():
                        percentage = (count / len(anno_df)) * 100
                        f.write(f"    {annotation}: {count} ({percentage:.1f}%)\n")
                
                if sample_name in all_gene_lists:
                    gene_lists = all_gene_lists[sample_name]
                    f.write(f"  Associated genes (all): {len(gene_lists.get('all_peaks', []))}\n")
                    f.write(f"  Associated genes (proximal): {len(gene_lists.get('proximal_peaks', []))}\n")
                    f.write(f"  Associated genes (distal): {len(gene_lists.get('distal_peaks', []))}\n")
                
                f.write("\n")
        
        logger.info(f"Saved analysis summary to {summary_file}")


def main():
    """Main function to run peaks annotation analysis."""
    parser = argparse.ArgumentParser(
        description="Peaks Annotation Analysis using GenomicFeatures and ChIPseeker (Fixed Version)"
    )
    parser.add_argument(
        "--peaks-dir", 
        required=True,
        help="Directory containing peak files"
    )
    parser.add_argument(
        "--gff-file",
        required=True,
        help="Path to GFF annotation file"
    )
    parser.add_argument(
        "--output-dir",
        default="output/peaks_annotation",
        help="Output directory for results"
    )
    
    args = parser.parse_args()
    
    try:
        # Initialize analyzer
        analyzer = PeaksAnnotationAnalyzer(output_dir=args.output_dir)
        
        # Run analysis
        results = analyzer.run_analysis(args.peaks_dir, args.gff_file)
        
        print(f"\nAnalysis completed successfully!")
        print(f"Results saved to: {args.output_dir}")
        print(f"Processed {len(results)} samples")
        
    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main() 