#!/usr/bin/env python3
"""
RNA-ATAC Overlapping Genes Correlation Analysis

This script analyzes the correlation between RNA-seq and ATAC-seq log2FoldChange expression
for overlapping genes between different clusters using various peak selection strategies.

Date: 2025 (HS_ABA project)
"""

import os
import sys
import logging
import argparse
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Union
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import pearsonr, spearmanr
import warnings
import re

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('rna_atac_overlapping_correlation_analysis.log'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

# Suppress warnings
warnings.filterwarnings('ignore')


class RNAATACOverlappingCorrelationAnalyzer:
    """
    Analyzer for RNA-seq and ATAC-seq overlapping genes correlation analysis.
    
    This class provides functionality to:
    - Load overlapping genes from specified directory with prefix filtering
    - Extract RNA-seq log2FoldChange values from DESeq2 files
    - Map genes to peaks using annotation data
    - Calculate correlations between RNA-seq and ATAC-seq data
    - Generate correlation heatmaps with significance indicators
    - Export data for R ComplexHeatmap customization
    """
    
    def __init__(self, 
                 overlapping_genes_dir: str,
                 rna_deseq2_dir: str,
                 atac_deseq2_dir: str,
                 peaks_annotation_dir: str,
                 dar_analysis_dir: str,
                 output_dir: str,
                 correlation_method: str = 'pearson',
                 gene_prefix: str = None):
        """
        Initialize the analyzer with data directories.
        
        Args:
            overlapping_genes_dir: Directory containing overlapping genes files
            rna_deseq2_dir: Directory containing RNA-seq DESeq2 files
            atac_deseq2_dir: Directory containing ATAC-seq DESeq2 files
            peaks_annotation_dir: Directory containing peaks annotation files
            dar_analysis_dir: Directory containing DAR analysis results
            output_dir: Directory for output files
            correlation_method: Correlation method ('pearson' or 'spearman')
            gene_prefix: Prefix to filter genes (e.g., 'Mp')
        """
        self.overlapping_genes_dir = Path(overlapping_genes_dir)
        self.rna_deseq2_dir = Path(rna_deseq2_dir)
        self.atac_deseq2_dir = Path(atac_deseq2_dir)
        self.peaks_annotation_dir = Path(peaks_annotation_dir)
        self.dar_analysis_dir = Path(dar_analysis_dir)
        self.output_dir = Path(output_dir)
        self.correlation_method = correlation_method.lower()
        self.gene_prefix = gene_prefix
        
        # Validate correlation method
        if self.correlation_method not in ['pearson', 'spearman']:
            raise ValueError("correlation_method must be 'pearson' or 'spearman'")
        
        # Create output directories
        self.output_dir.mkdir(parents=True, exist_ok=True)
        (self.output_dir / "correlation_results").mkdir(exist_ok=True)
        (self.output_dir / "extracted_data").mkdir(exist_ok=True)
        (self.output_dir / "R_data").mkdir(exist_ok=True)
        
        # Data storage
        self.rna_deseq2_data = {}
        self.atac_deseq2_data = {}
        self.peaks_annotation = {}
        self.overlapping_genes = {}
        self.dar_analysis_results = {}
        
        # Analysis results
        self.correlation_matrix = None
        self.pvalue_matrix = None
        self.cluster_pairs = []
        
        # Define expected clusters - will be dynamically determined from actual data
        self.expected_rna_clusters = []
        self.expected_atac_clusters = []
        
        # Define genotype groups - map RNA-seq sample names to ATAC-seq sample names
        self.genotype_groups = {
            'Tak1_heat_Tak1_NHS': 'Tak1_HS_Tak1_CK',
            'hsfa_heat_hsfa_NHS': 'hsfa_HS_hsfa_CK',
            'hsfb_heat_hsfb_NHS': 'hsfb_HS_hsfb_CK',
            'dko_heat_dko_NHS': 'dko_HS_dko_CK'
        }
        
        logger.info(f"RNA-ATAC Overlapping Correlation Analyzer initialized with {self.correlation_method} correlation")
        if self.gene_prefix:
            logger.info(f"Gene prefix filter: {self.gene_prefix}")
    
    def load_rna_deseq2_data(self) -> Dict[str, pd.DataFrame]:
        """Load RNA-seq DESeq2 data from CSV files."""
        logger.info("Loading RNA-seq DESeq2 data...")
        
        rna_data = {}
        for file_path in self.rna_deseq2_dir.glob("*.csv"):
            sample_name = file_path.stem
            try:
                df = pd.read_csv(file_path, index_col=0)
                if 'log2FoldChange' in df.columns:
                    rna_data[sample_name] = df[['log2FoldChange']]
                    logger.info(f"Loaded RNA-seq data for {sample_name}: {len(df)} genes")
                else:
                    logger.warning(f"No log2FoldChange column found in {file_path}")
            except Exception as e:
                logger.error(f"Error loading RNA-seq data from {file_path}: {e}")
        
        return rna_data
    
    def load_atac_deseq2_data(self) -> Dict[str, pd.DataFrame]:
        """Load ATAC-seq DESeq2 data from CSV files."""
        logger.info("Loading ATAC-seq DESeq2 data...")
        
        atac_data = {}
        for file_path in self.atac_deseq2_dir.glob("*.csv"):
            sample_name = file_path.stem
            try:
                df = pd.read_csv(file_path, index_col=0)
                if 'log2FoldChange' in df.columns:
                    atac_data[sample_name] = df[['log2FoldChange']]
                    logger.info(f"Loaded ATAC-seq data for {sample_name}: {len(df)} peaks")
                else:
                    logger.warning(f"No log2FoldChange column found in {file_path}")
            except Exception as e:
                logger.error(f"Error loading ATAC-seq data from {file_path}: {e}")
        
        return atac_data
    
    def load_peaks_annotation(self) -> Dict[str, pd.DataFrame]:
        """Load peaks annotation data from CSV files."""
        logger.info("Loading peaks annotation data...")
        
        peaks_annotation = {}
        for file_path in self.peaks_annotation_dir.glob("*_peaks_annotation_results.csv"):
            cluster_name = file_path.stem.replace('_peaks_annotation_results', '')
            try:
                df = pd.read_csv(file_path)
                required_cols = ['geneId', 'Peak', 'distanceToTSS', 'peak_type']
                if all(col in df.columns for col in required_cols):
                    peaks_annotation[cluster_name] = df[required_cols]
                    logger.info(f"Loaded annotation for {cluster_name}: {len(df)} gene-peak mappings")
                else:
                    logger.warning(f"Missing required columns in {file_path}")
            except Exception as e:
                logger.error(f"Error loading annotation from {file_path}: {e}")
        
        return peaks_annotation
    
    def load_dar_analysis_results(self) -> Dict[str, List[str]]:
        """Load DAR analysis results from text files."""
        logger.info("Loading DAR analysis results...")
        
        dar_results = {}
        for file_path in self.dar_analysis_dir.glob("*.txt"):
            cluster_name = file_path.stem
            try:
                with open(file_path, 'r') as f:
                    lines = f.readlines()
                    peaks = []
                    for line in lines:
                        line = line.strip()
                        if line and line != 'Peak' and line != 'Peaks' and not line.startswith('#'):
                            peaks.append(line)
                    
                dar_results[cluster_name] = peaks
                logger.info(f"Loaded DAR results for {cluster_name}: {len(peaks)} peaks")
            except Exception as e:
                logger.error(f"Error loading DAR results from {file_path}: {e}")
        
        return dar_results
    
    def extract_cluster_names(self, cluster_pair: str) -> Tuple[str, str]:
        """Extract RNA and ATAC cluster names from cluster pair string."""
        # Expected format: RNAcluster_vs_ATACcluster_genes_overlapping_genes
        if "_vs_" in cluster_pair and "_overlapping_genes" in cluster_pair:
            parts = cluster_pair.replace("_overlapping_genes", "").split("_vs_")
            if len(parts) == 2:
                # Remove "_genes" suffix if present
                rna_cluster = parts[0].replace("_genes", "")
                atac_cluster = parts[1].replace("_genes", "")
                return rna_cluster, atac_cluster
        
        logger.error(f"Invalid cluster pair format: {cluster_pair}")
        return None, None
    
    def discover_clusters(self):
        """Dynamically discover RNA and ATAC clusters from overlapping genes files."""
        logger.info("Discovering clusters from overlapping genes files...")
        
        rna_clusters = set()
        atac_clusters = set()
        
        for file_path in self.overlapping_genes_dir.glob("*.txt"):
            if "summary" in file_path.stem.lower():
                continue
            
            filename = file_path.stem
            if "_vs_" in filename and "_overlapping_genes" in filename:
                cluster_pair = filename
                rna_cluster, atac_cluster = self.extract_cluster_names(cluster_pair)
                
                if rna_cluster and atac_cluster:
                    rna_clusters.add(rna_cluster)
                    atac_clusters.add(atac_cluster)
        
        self.expected_rna_clusters = sorted(list(rna_clusters))
        self.expected_atac_clusters = sorted(list(atac_clusters))
        
        logger.info(f"Discovered {len(self.expected_rna_clusters)} RNA-seq clusters: {self.expected_rna_clusters}")
        logger.info(f"Discovered {len(self.expected_atac_clusters)} ATAC-seq clusters: {self.expected_atac_clusters}")
        
        return self.expected_rna_clusters, self.expected_atac_clusters
    
    def load_overlapping_genes(self) -> Dict[str, List[str]]:
        """Load overlapping genes from text files with prefix filtering."""
        logger.info("Loading overlapping genes...")
        
        overlapping_genes = {}
        for file_path in self.overlapping_genes_dir.glob("*.txt"):
            if "summary" in file_path.stem.lower():
                continue
            
            # Extract cluster pair from filename
            filename = file_path.stem
            if "_vs_" in filename and "_overlapping_genes" in filename:
                cluster_pair = filename
            else:
                logger.warning(f"Skipping file with unexpected naming format: {filename}")
                continue
            
            try:
                with open(file_path, 'r') as f:
                    lines = f.readlines()
                    
                genes = []
                for line in lines:
                    line = line.strip()
                    if (line and 
                        line != 'Gene' and 
                        line != 'GeneID' and 
                        not line.startswith('#') and
                        not line.startswith('gene') and
                        len(line) > 1 and
                        not line.isupper()):  # Skip header lines like "GeneID"
                        
                        # Apply prefix filtering if specified
                        if self.gene_prefix is None or line.startswith(self.gene_prefix):
                            genes.append(line)
                
                if genes:
                    overlapping_genes[cluster_pair] = genes
                    logger.info(f"Loaded {len(genes)} overlapping genes for {cluster_pair} from {file_path.name}")
                    if self.gene_prefix:
                        logger.info(f"Filtered by prefix: {self.gene_prefix}")
                else:
                    logger.warning(f"No genes found for {cluster_pair} after prefix filtering")
                
            except Exception as e:
                logger.error(f"Error loading overlapping genes from {file_path}: {e}")
        
        return overlapping_genes
    
    def select_peak_by_strategy(self, peaks_info: List[Dict], strategy: str) -> Optional[Dict]:
        """Select peak based on specified strategy."""
        if not peaks_info:
            return None
        
        if strategy == 'closest_to_tss':
            return min(peaks_info, key=lambda x: x['distance'])
        elif strategy == 'proximal_median':
            proximal_peaks = [p for p in peaks_info if p['peak_type'] == 'proximal']
            if proximal_peaks:
                proximal_peaks.sort(key=lambda x: x['distance'])
                median_idx = len(proximal_peaks) // 2
                return proximal_peaks[median_idx]
            else:
                return min(peaks_info, key=lambda x: x['distance'])
        elif strategy == 'distal_median':
            distal_peaks = [p for p in peaks_info if p['peak_type'] == 'distal']
            if distal_peaks:
                distal_peaks.sort(key=lambda x: x['distance'])
                median_idx = len(distal_peaks) // 2
                return distal_peaks[median_idx]
            else:
                return min(peaks_info, key=lambda x: x['distance'])
        elif strategy == 'all_median':
            all_peaks = sorted(peaks_info, key=lambda x: x['distance'])
            median_idx = len(all_peaks) // 2
            return all_peaks[median_idx]
        else:
            raise ValueError(f"Unknown strategy: {strategy}")
    
    def calculate_correlations(self, peak_strategy: str = 'all_median') -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Calculate correlations between RNA-seq and ATAC-seq data for overlapping genes.
        
        Args:
            peak_strategy: Strategy for selecting peaks ('all_median', 'closest_to_tss', 
                          'proximal_median', 'distal_median')
        
        Returns:
            Tuple of (correlation_matrix, pvalue_matrix)
        """
        logger.info("Starting correlation calculations...")
        
        # Debug: Print available keys in peaks_annotation
        logger.info(f"Available keys in peaks_annotation: {list(self.peaks_annotation.keys())}")
        
        # Initialize matrices
        rna_genotypes = list(self.genotype_groups.keys())
        atac_genotypes = list(self.genotype_groups.values())
        correlation_matrix = pd.DataFrame(index=rna_genotypes, columns=atac_genotypes, dtype=float)
        pvalue_matrix = pd.DataFrame(index=rna_genotypes, columns=atac_genotypes, dtype=float)
        
        # Calculate correlations for each cluster pair
        for cluster_pair, overlapping_genes in self.overlapping_genes.items():
            logger.info(f"\nProcessing cluster pair: {cluster_pair}")
            logger.info(f"Number of overlapping genes: {len(overlapping_genes)}")
            
            # Debug: Print the actual cluster_pair value
            logger.info(f"Cluster pair value: '{cluster_pair}'")
            
            # Extract RNA and ATAC cluster names
            if '_vs_' in cluster_pair:
                rna_cluster, atac_cluster = cluster_pair.split('_vs_')
            else:
                logger.warning(f"Could not parse cluster pair: {cluster_pair}")
                continue
            
            # Debug: Print the extracted cluster names
            logger.info(f"RNA cluster: '{rna_cluster}'")
            logger.info(f"ATAC cluster: '{atac_cluster}'")
            
            # Calculate correlations for each genotype
            for rna_genotype in rna_genotypes:
                atac_genotype = self.genotype_groups[rna_genotype]
                
                # Debug: Print genotype mapping
                logger.info(f"RNA genotype: {rna_genotype} -> ATAC genotype: {atac_genotype}")
                
                # Check if we have data for both genotypes
                if rna_genotype not in self.rna_deseq2_data:
                    logger.warning(f"Missing data for {rna_genotype}")
                    continue
                    
                if atac_genotype not in self.atac_deseq2_data:
                    logger.warning(f"Missing data for {atac_genotype}")
                    continue
                
                # Initialize lists for correlation calculation
                rna_values = []
                atac_values = []
                valid_genes = []
                genes_not_found = []
                genes_no_peaks = []
                
                # Process each overlapping gene
                for gene in overlapping_genes:
                    # Try to find gene in RNA-seq data with flexible matching
                    rna_val = None
                    if gene in self.rna_deseq2_data[rna_genotype].index:
                        rna_val = self.rna_deseq2_data[rna_genotype].loc[gene, 'log2FoldChange']
                    else:
                        # Try flexible matching
                        for idx in self.rna_deseq2_data[rna_genotype].index:
                            if (gene.lower() == idx.lower() or 
                                gene.replace('_', '') == idx.replace('_', '')):
                                rna_val = self.rna_deseq2_data[rna_genotype].loc[idx, 'log2FoldChange']
                                logger.debug(f"Found gene match: {gene} -> {idx}")
                                break
                    
                    if rna_val is None:
                        genes_not_found.append(gene)
                        continue
                    
                    # Find ATAC-seq data for this gene - ignore cluster names and search across all annotation clusters
                    gene_peaks_found = False
                    peaks_info = []
                    
                    # Search for the gene across all available peaks annotation clusters
                    for annotation_cluster, annotation_df in self.peaks_annotation.items():
                        gene_peaks = annotation_df[annotation_df['geneId'] == gene]
                        
                        if not gene_peaks.empty:
                            gene_peaks_found = True
                            logger.debug(f"Found gene {gene} in annotation cluster {annotation_cluster}")
                            
                            for _, row in gene_peaks.iterrows():
                                # Check if peak is in DAR analysis results
                                if annotation_cluster in self.dar_analysis_results:
                                    if row['Peak'] in self.dar_analysis_results[annotation_cluster]:
                                        peaks_info.append({
                                            'peak_id': row['Peak'],
                                            'distance': abs(row['distanceToTSS']),
                                            'peak_type': row['peak_type']
                                        })
                                
                                # If no peaks found in DAR results, use all peaks
                                if not peaks_info and row['Peak'] in self.atac_deseq2_data[atac_genotype].index:
                                    peaks_info.append({
                                        'peak_id': row['Peak'],
                                        'distance': abs(row['distanceToTSS']),
                                        'peak_type': row['peak_type']
                                    })
                            
                            # If we found peaks for this gene, no need to check other annotation clusters
                            break
                    
                    if gene_peaks_found and peaks_info:
                        selected_peak = self.select_peak_by_strategy(peaks_info, peak_strategy)
                        if selected_peak and selected_peak['peak_id'] in self.atac_deseq2_data[atac_genotype].index:
                            atac_val = self.atac_deseq2_data[atac_genotype].loc[selected_peak['peak_id'], 'log2FoldChange']
                            rna_values.append(rna_val)
                            atac_values.append(atac_val)
                            valid_genes.append(gene)
                            logger.debug(f"Successfully processed gene {gene} with peak {selected_peak['peak_id']}")
                        else:
                            genes_no_peaks.append(gene)
                            logger.debug(f"Gene {gene} has no valid peaks after selection")
                    else:
                        genes_no_peaks.append(gene)
                        logger.debug(f"Gene {gene} not found in any annotation cluster")
                
                # Log detailed information for debugging
                logger.info(f"{cluster_pair} - {rna_genotype}:")
                logger.info(f"  Total overlapping genes: {len(overlapping_genes)}")
                logger.info(f"  Genes not found in RNA-seq: {len(genes_not_found)}")
                logger.info(f"  Genes with no peaks: {len(genes_no_peaks)}")
                logger.info(f"  Valid genes for correlation: {len(valid_genes)}")
                if genes_not_found:
                    logger.info(f"  Examples of genes not found: {genes_not_found[:5]}")
                if genes_no_peaks:
                    logger.info(f"  Examples of genes with no peaks: {genes_no_peaks[:5]}")
                if valid_genes:
                    logger.info(f"  Examples of valid genes: {valid_genes[:5]}")
                
                # Calculate correlation if we have enough data
                if len(valid_genes) >= 3:
                    correlation, pvalue = spearmanr(rna_values, atac_values)
                    correlation_matrix.loc[rna_genotype, atac_genotype] = correlation
                    pvalue_matrix.loc[rna_genotype, atac_genotype] = pvalue
                    logger.info(f"  Correlation: {correlation:.4f}, p-value: {pvalue:.4f}")
                else:
                    logger.warning(f"Insufficient data for {cluster_pair} - {rna_genotype}: {len(valid_genes)} valid genes (need >=3)")
                    correlation_matrix.loc[rna_genotype, atac_genotype] = np.nan
                    pvalue_matrix.loc[rna_genotype, atac_genotype] = np.nan
        
        return correlation_matrix, pvalue_matrix
    
    def get_significance_symbol(self, p_value: float) -> str:
        """Get significance symbol based on p-value."""
        if np.isnan(p_value):
            return ""
        
        if p_value < 0.001:
            return "***"
        elif p_value < 0.01:
            return "**"
        elif p_value < 0.05:
            return "*"
        else:
            return ""
    
    def create_correlation_heatmap(self) -> plt.Figure:
        """Create correlation heatmap with significance indicators."""
        if self.correlation_matrix is None or self.pvalue_matrix is None:
            logger.error("No correlation data available for heatmap")
            return None
        
        # Create figure
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
        
        # Convert matrices to numeric
        corr_matrix_numeric = self.correlation_matrix.astype(float).fillna(0)
        pval_matrix_numeric = self.pvalue_matrix.astype(float).fillna(1)
        
        # 1. Correlation heatmap
        sns.heatmap(corr_matrix_numeric, 
                   annot=True, 
                   fmt='.3f',
                   cmap='RdBu_r',
                   center=0,
                   vmin=-1, 
                   vmax=1,
                   cbar_kws={'label': f'{self.correlation_method.capitalize()} Correlation'},
                   ax=ax1)
        
        # Add significance indicators
        for i in range(len(self.correlation_matrix.index)):
            for j in range(len(self.correlation_matrix.columns)):
                p_value = self.pvalue_matrix.iloc[i, j]
                if not np.isnan(p_value):
                    significance = self.get_significance_symbol(p_value)
                    if significance:
                        ax1.text(j + 0.5, i + 0.5, significance, 
                                ha='center', va='center', 
                                fontsize=12, fontweight='bold', color='black')
        
        ax1.set_title(f'RNA-ATAC Correlation Heatmap\n({self.correlation_method.capitalize()})', 
                      fontsize=14, fontweight='bold')
        ax1.set_xlabel('ATAC-seq Clusters', fontsize=12)
        ax1.set_ylabel('RNA-seq Clusters', fontsize=12)
        
        # 2. P-value heatmap
        sns.heatmap(pval_matrix_numeric, 
                   annot=True, 
                   fmt='.3e',
                   cmap='Reds_r',
                   vmin=0, 
                   vmax=0.05,
                   cbar_kws={'label': 'P-value'},
                   ax=ax2)
        
        ax2.set_title('P-value Heatmap', fontsize=14, fontweight='bold')
        ax2.set_xlabel('ATAC-seq Clusters', fontsize=12)
        ax2.set_ylabel('RNA-seq Clusters', fontsize=12)
        
        plt.tight_layout()
        return fig
    
    def save_correlation_results(self, peak_strategy: str):
        """Save correlation analysis results to files."""
        timestamp = pd.Timestamp.now().strftime("%Y%m%d_%H%M%S")
        results_dir = self.output_dir / "correlation_results"
        
        # Save correlation matrix
        corr_file = results_dir / f"correlation_matrix_{peak_strategy}_{timestamp}.csv"
        self.correlation_matrix.to_csv(corr_file)
        logger.info(f"Saved correlation matrix to: {corr_file}")
        
        # Save p-value matrix
        pval_file = results_dir / f"pvalue_matrix_{peak_strategy}_{timestamp}.csv"
        self.pvalue_matrix.to_csv(pval_file)
        logger.info(f"Saved p-value matrix to: {pval_file}")
        
        # Save heatmap data for R ComplexHeatmap
        heatmap_data = []
        for i in range(len(self.correlation_matrix.index)):
            for j in range(len(self.correlation_matrix.columns)):
                rna_cluster = self.correlation_matrix.index[i]
                atac_cluster = self.correlation_matrix.columns[j]
                correlation = self.correlation_matrix.iloc[i, j]
                p_value = self.pvalue_matrix.iloc[i, j]
                significance = self.get_significance_symbol(p_value)
                
                heatmap_data.append({
                    'RNA_cluster': rna_cluster,
                    'ATAC_cluster': atac_cluster,
                    'correlation': correlation,
                    'p_value': p_value,
                    'significance': significance
                })
        
        # Save as CSV for R
        r_data_file = results_dir / f"heatmap_data_for_R_{peak_strategy}_{timestamp}.csv"
        heatmap_df = pd.DataFrame(heatmap_data)
        heatmap_df.to_csv(r_data_file, index=False)
        logger.info(f"Saved R-compatible data to: {r_data_file}")
        
        # Save analysis summary
        summary_file = results_dir / f"analysis_summary_{peak_strategy}_{timestamp}.txt"
        with open(summary_file, 'w') as f:
            f.write("RNA-seq vs ATAC-seq Overlapping Genes Correlation Analysis Summary\n")
            f.write("=" * 70 + "\n\n")
            f.write(f"Analysis Date: {timestamp}\n")
            f.write(f"Correlation Method: {self.correlation_method.capitalize()}\n")
            f.write(f"Peak Selection Strategy: {peak_strategy}\n")
            if self.gene_prefix:
                f.write(f"Gene Prefix Filter: {self.gene_prefix}\n")
            f.write(f"Number of RNA-seq clusters: {len(self.correlation_matrix.index)}\n")
            f.write(f"Number of ATAC-seq clusters: {len(self.correlation_matrix.columns)}\n")
            f.write(f"Matrix size: {len(self.correlation_matrix.index)}x{len(self.correlation_matrix.columns)}\n\n")
            
            f.write("RNA-seq Clusters:\n")
            for cluster in self.expected_rna_clusters:
                f.write(f"  - {cluster}\n")
            f.write("\nATAC-seq Clusters:\n")
            for cluster in self.expected_atac_clusters:
                f.write(f"  - {cluster}\n")
            f.write("\n")
            
            f.write("Correlation Matrix:\n")
            f.write(self.correlation_matrix.to_string())
            f.write("\n\nP-value Matrix:\n")
            f.write(self.pvalue_matrix.to_string())
            
            f.write("\n\n")
            
            # Summary statistics
            try:
                corr_values = pd.to_numeric(self.correlation_matrix.values.flatten(), errors='coerce')
                valid_correlations = corr_values[~np.isnan(corr_values)]
                
                if len(valid_correlations) > 0:
                    f.write(f"Mean correlation: {np.mean(valid_correlations):.3f}\n")
                    f.write(f"Median correlation: {np.median(valid_correlations):.3f}\n")
                    f.write(f"Min correlation: {np.min(valid_correlations):.3f}\n")
                    f.write(f"Max correlation: {np.max(valid_correlations):.3f}\n")
                
                # Count significant correlations
                pval_values = pd.to_numeric(self.pvalue_matrix.values.flatten(), errors='coerce')
                valid_pvals = pval_values[~np.isnan(pval_values)]
                significant_pairs_001 = np.sum(valid_pvals < 0.001)
                significant_pairs_01 = np.sum(valid_pvals < 0.01)
                significant_pairs_05 = np.sum(valid_pvals < 0.05)
                
                f.write(f"\nP-value Significance:\n")
                f.write(f"Significant correlations (p < 0.001): {significant_pairs_001}\n")
                f.write(f"Significant correlations (p < 0.01): {significant_pairs_01}\n")
                f.write(f"Significant correlations (p < 0.05): {significant_pairs_05}\n")
                
            except Exception as e:
                f.write(f"Error calculating statistics: {e}\n")
        
        logger.info(f"Saved analysis summary to: {summary_file}")
    
    def generate_r_script(self, peak_strategy: str):
        """Generate R script for ComplexHeatmap customization."""
        timestamp = pd.Timestamp.now().strftime("%Y%m%d_%H%M%S")
        r_scripts_dir = self.output_dir / "R_data"
        
        r_script = f"""# RNA-ATAC Correlation Heatmap R Script
# Generated on: {timestamp}
# Peak Strategy: {peak_strategy}
# Correlation Method: {self.correlation_method.capitalize()}

# Load required libraries
library(ComplexHeatmap)
library(circlize)
library(dplyr)

# Set working directory (modify as needed)
# setwd("path/to/your/output/directory")

# Load correlation data
heatmap_data <- read.csv("correlation_results/heatmap_data_for_R_{peak_strategy}_{timestamp}.csv")

# Create correlation matrix
correlation_matrix <- heatmap_data %>%
  select(RNA_cluster, ATAC_cluster, correlation) %>%
  pivot_wider(names_from = ATAC_cluster, values_from = correlation) %>%
  column_to_rownames("RNA_cluster") %>%
  as.matrix()

# Create p-value matrix
pvalue_matrix <- heatmap_data %>%
  select(RNA_cluster, ATAC_cluster, p_value) %>%
  pivot_wider(names_from = ATAC_cluster, values_from = p_value) %>%
  column_to_rownames("RNA_cluster") %>%
  as.matrix()

# Create significance matrix
significance_matrix <- heatmap_data %>%
  select(RNA_cluster, ATAC_cluster, significance) %>%
  pivot_wider(names_from = ATAC_cluster, values_from = significance) %>%
  column_to_rownames("RNA_cluster") %>%
  as.matrix()

# Define color function for correlation values
col_fun <- colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B"))

# Create correlation heatmap
ht_corr <- Heatmap(
  correlation_matrix,
  name = "Correlation",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  cell_fun = function(j, i, x, y, width, height, fill) {{
    if (!is.na(correlation_matrix[i, j])) {{
      # Add correlation value
      grid.text(
        sprintf("%.3f", correlation_matrix[i, j]),
        x, y,
        gp = gpar(fontsize = 8, fontface = "bold")
      )
      
      # Add significance symbols
      if (significance_matrix[i, j] != "") {{
        grid.text(
          significance_matrix[i, j],
          x, y + height * 0.3,
          gp = gpar(fontsize = 12, fontface = "bold", col = "black")
        )
      }}
    }}
  }}
)

# Create p-value heatmap
ht_pval <- Heatmap(
  pvalue_matrix,
  name = "P-value",
  col = colorRamp2(c(0, 0.05), c("red", "white")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 10),
  cell_fun = function(j, i, x, y, width, height, fill) {{
    if (!is.na(pvalue_matrix[i, j])) {{
      grid.text(
        sprintf("%.3e", pvalue_matrix[i, j]),
        x, y,
        gp = gpar(fontsize = 6)
      )
    }}
  }}
)

# Combine heatmaps
combined_ht <- ht_corr + ht_pval

# Save heatmap
pdf("rna_atac_correlation_heatmap_{peak_strategy}_{timestamp}.pdf", width = 16, height = 10)
draw(combined_ht, heatmap_legend_side = "bottom")
dev.off()

# Save individual heatmaps
pdf("correlation_heatmap_{peak_strategy}_{timestamp}.pdf", width = 12, height = 8)
draw(ht_corr, heatmap_legend_side = "bottom")
dev.off()

pdf("pvalue_heatmap_{peak_strategy}_{timestamp}.pdf", width = 12, height = 8)
draw(ht_pval, heatmap_legend_side = "bottom")
dev.off()

cat("Heatmaps generated successfully!\\n")
cat("Files saved:\\n")
cat("- rna_atac_correlation_heatmap_{peak_strategy}_{timestamp}.pdf (combined)\\n")
cat("- correlation_heatmap_{peak_strategy}_{timestamp}.pdf (correlation only)\\n")
cat("- pvalue_heatmap_{peak_strategy}_{timestamp}.pdf (p-value only)\\n")
"""
        
        r_file = r_scripts_dir / f"complexheatmap_script_{peak_strategy}_{timestamp}.R"
        with open(r_file, 'w') as f:
            f.write(r_script)
        
        logger.info(f"Generated R script: {r_file}")
        
        # Also save the script content to a text file for easy viewing
        txt_file = r_scripts_dir / f"complexheatmap_script_{peak_strategy}_{timestamp}.txt"
        with open(txt_file, 'w') as f:
            f.write(r_script)
        
        logger.info(f"Saved R script content to: {txt_file}")
    
    def run_analysis(self, peak_strategy: str = 'all_median') -> None:
        """
        Run the complete RNA-ATAC overlapping correlation analysis.
        
        Args:
            peak_strategy: Strategy for selecting peaks ('all_median', 'closest_to_tss', 
                          'proximal_median', 'distal_median')
        """
        logger.info("Starting RNA-ATAC overlapping correlation analysis...")
        
        # Load all data
        logger.info("Loading data...")
        self.rna_deseq2_data = self.load_rna_deseq2_data()
        self.atac_deseq2_data = self.load_atac_deseq2_data()
        self.peaks_annotation = self.load_peaks_annotation()
        self.dar_analysis_results = self.load_dar_analysis_results()
        self.overlapping_genes = self.load_overlapping_genes()
        
        if not self.overlapping_genes:
            logger.error("No overlapping genes found. Analysis cannot proceed.")
            return
        
        # Calculate correlations
        logger.info("Calculating correlations...")
        correlation_matrix, pvalue_matrix = self.calculate_correlations(peak_strategy)
        
        # Store results as instance variables
        self.correlation_matrix = correlation_matrix
        self.pvalue_matrix = pvalue_matrix
        
        # Generate visualizations
        self.generate_r_script(peak_strategy)
        
        # Save results
        self.save_correlation_results(peak_strategy)
        
        logger.info("Analysis completed successfully!")


def main():
    """Main function to run the RNA-ATAC overlapping correlation analysis."""
    parser = argparse.ArgumentParser(
        description="Analyze RNA-seq and ATAC-seq overlapping genes correlations"
    )
    
    parser.add_argument(
        "--overlapping-genes-dir",
        type=str,
        default="output/upset_overlap_analysis/overlapping_genes",
        help="Directory containing overlapping genes files"
    )
    
    parser.add_argument(
        "--rna-deseq2-dir",
        type=str,
        default="data/HS_RNA_deseq2",
        help="Directory containing RNA-seq DESeq2 files"
    )
    
    parser.add_argument(
        "--atac-deseq2-dir",
        type=str,
        default="data/HS_ATAC_deseq2",
        help="Directory containing ATAC-seq DESeq2 files"
    )
    
    parser.add_argument(
        "--peaks-annotation-dir",
        type=str,
        default="output/HS_DAR_annotation",
        help="Directory containing peaks annotation files"
    )
    
    parser.add_argument(
        "--dar-analysis-dir",
        type=str,
        default="output/HS_DAR_analysis/all",
        help="Directory containing DAR analysis results"
    )
    
    parser.add_argument(
        "--output-dir",
        type=str,
        default="output/RNA_ATAC_overlapping_correlation_analysis",
        help="Directory for output files"
    )
    
    parser.add_argument(
        "--peak-strategy",
        type=str,
        choices=['closest_to_tss', 'proximal_median', 'distal_median', 'all_median'],
        default='closest_to_tss',
        help="Strategy for selecting peaks to represent genes"
    )
    
    parser.add_argument(
        "--correlation-method",
        type=str,
        choices=['pearson', 'spearman'],
        default='pearson',
        help="Correlation method to use ('pearson' or 'spearman')"
    )
    
    parser.add_argument(
        "--gene-prefix",
        type=str,
        default=None,
        help="Prefix to filter genes (e.g., 'Mp')"
    )
    
    args = parser.parse_args()
    
    # Create analyzer and run analysis
    analyzer = RNAATACOverlappingCorrelationAnalyzer(
        overlapping_genes_dir=args.overlapping_genes_dir,
        rna_deseq2_dir=args.rna_deseq2_dir,
        atac_deseq2_dir=args.atac_deseq2_dir,
        peaks_annotation_dir=args.peaks_annotation_dir,
        dar_analysis_dir=args.dar_analysis_dir,
        output_dir=args.output_dir,
        correlation_method=args.correlation_method,
        gene_prefix=args.gene_prefix
    )
    
    analyzer.run_analysis(peak_strategy=args.peak_strategy)


if __name__ == "__main__":
    main()
