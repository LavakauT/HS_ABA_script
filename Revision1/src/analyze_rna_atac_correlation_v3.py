#!/usr/bin/env python3
"""
RNA-seq and ATAC-seq Cluster Genes Overlapping Gene Expression Correlation Analysis V3

This script analyzes the correlation between RNA-seq and ATAC-seq log2FoldChange expression
for overlapping genes between different clusters using corrected data extraction methods.

Author: AI Assistant
Date: 2025-01-27
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
from statsmodels.stats.multitest import multipletests
import warnings

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('rna_atac_correlation_analysis_v3.log'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

# Suppress warnings
warnings.filterwarnings('ignore')


class RNAATACCorrelationAnalyzerV3:
    """
    Enhanced analyzer for RNA-seq and ATAC-seq cluster genes overlapping gene expression correlations.
    
    This class provides functionality to:
    - Load RNA-seq and ATAC-seq cluster results (genes and peaks)
    - Extract log2FoldChange values from DESeq2 files
    - Map genes to peaks using annotation data
    - Calculate correlations between RNA-seq and ATAC-seq data
    - Generate correlation heatmaps with significance indicators
    - Generate R scripts for further customization
    - Perform FDR correction for p-values
    - Create individual scatter plots for each cluster pair
    """
    
    def __init__(self, 
                 rna_deseq2_dir: str,
                 atac_deseq2_dir: str,
                 rna_cluster_results_dir: str,
                 atac_cluster_results_dir: str,
                 peaks_annotation_dir: str,
                 overlapping_genes_dir: str,
                 output_dir: str,
                 correlation_method: str = 'pearson',
                 fdr_method: str = 'fdr_bh'):
        """
        Initialize the analyzer with data directories.
        
        Args:
            rna_deseq2_dir: Directory containing RNA-seq DESeq2 files
            atac_deseq2_dir: Directory containing ATAC-seq DESeq2 files
            rna_cluster_results_dir: Directory containing RNA-seq cluster results
            atac_cluster_results_dir: Directory containing ATAC-seq cluster results
            peaks_annotation_dir: Directory containing peaks annotation files
            overlapping_genes_dir: Directory containing overlapping genes files
            output_dir: Directory for output files
            correlation_method: Correlation method ('pearson' or 'spearman')
            fdr_method: FDR correction method ('fdr_bh', 'fdr_by', 'bonferroni', 'holm')
        """
        self.rna_deseq2_dir = Path(rna_deseq2_dir)
        self.atac_deseq2_dir = Path(atac_deseq2_dir)
        self.rna_cluster_results_dir = Path(rna_cluster_results_dir)
        self.atac_cluster_results_dir = Path(atac_cluster_results_dir)
        self.peaks_annotation_dir = Path(peaks_annotation_dir)
        self.overlapping_genes_dir = Path(overlapping_genes_dir)
        self.output_dir = Path(output_dir)
        self.correlation_method = correlation_method.lower()
        self.fdr_method = fdr_method
        
        # Validate correlation method
        if self.correlation_method not in ['pearson', 'spearman']:
            raise ValueError("correlation_method must be 'pearson' or 'spearman'")
        
        # Validate FDR method
        valid_fdr_methods = ['fdr_bh', 'fdr_by', 'bonferroni', 'holm']
        if self.fdr_method not in valid_fdr_methods:
            raise ValueError(f"fdr_method must be one of {valid_fdr_methods}")
        
        # Create output directories
        self.output_dir.mkdir(parents=True, exist_ok=True)
        (self.output_dir / "extracted_data").mkdir(exist_ok=True)
        (self.output_dir / "peaks_to_genes").mkdir(exist_ok=True)
        (self.output_dir / "correlation_results").mkdir(exist_ok=True)
        (self.output_dir / "individual_correlations").mkdir(exist_ok=True)
        (self.output_dir / "R_scripts").mkdir(exist_ok=True)
        (self.output_dir / "scatter_plots").mkdir(exist_ok=True)
        
        # Data storage
        self.rna_deseq2_data = {}
        self.atac_deseq2_data = {}
        self.rna_cluster_genes = {}
        self.atac_cluster_peaks = {}
        self.peaks_annotation = {}
        self.overlapping_genes = {}
        
        # Analysis results
        self.correlation_matrix = None
        self.pvalue_matrix = None
        self.fdr_matrix = None  # New: FDR-corrected p-values
        self.cluster_pairs = []
        self.individual_correlations = {}  # Store individual correlation data
        
        # Define expected clusters (4x4 matrix)
        self.expected_rna_clusters = ['cluster_1', 'cluster_2', 'cluster_3', 'cluster_4']
        self.expected_atac_clusters = ['cluster_1', 'cluster_2', 'cluster_3', 'cluster_4']
        
        logger.info(f"RNA-ATAC Correlation Analyzer V3 initialized with {self.correlation_method} correlation")
        logger.info(f"FDR correction method: {self.fdr_method}")
        logger.info(f"Expected RNA clusters: {self.expected_rna_clusters}")
        logger.info(f"Expected ATAC clusters: {self.expected_atac_clusters}")
    
    def load_rna_deseq2_data(self) -> Dict[str, pd.DataFrame]:
        """
        Load RNA-seq DESeq2 data from CSV files.
        
        Returns:
            Dictionary mapping sample names to DESeq2 DataFrames
        """
        logger.info("Loading RNA-seq DESeq2 data...")
        
        rna_data = {}
        for file_path in self.rna_deseq2_dir.glob("*.csv"):
            sample_name = file_path.stem
            try:
                df = pd.read_csv(file_path, index_col=0)
                # Extract log2FoldChange column if it exists
                if 'log2FoldChange' in df.columns:
                    rna_data[sample_name] = df[['log2FoldChange']]
                    logger.info(f"Loaded RNA-seq data for {sample_name}: {len(df)} genes")
                else:
                    logger.warning(f"No log2FoldChange column found in {file_path}")
            except Exception as e:
                logger.error(f"Error loading RNA-seq data from {file_path}: {e}")
        
        return rna_data
    
    def load_atac_deseq2_data(self) -> Dict[str, pd.DataFrame]:
        """
        Load ATAC-seq DESeq2 data from CSV files.
        
        Returns:
            Dictionary mapping sample names to DESeq2 DataFrames
        """
        logger.info("Loading ATAC-seq DESeq2 data...")
        
        atac_data = {}
        for file_path in self.atac_deseq2_dir.glob("*.csv"):
            sample_name = file_path.stem
            try:
                df = pd.read_csv(file_path, index_col=0)
                # Extract log2FoldChange column if it exists
                if 'log2FoldChange' in df.columns:
                    atac_data[sample_name] = df[['log2FoldChange']]
                    logger.info(f"Loaded ATAC-seq data for {sample_name}: {len(df)} peaks")
                else:
                    logger.warning(f"No log2FoldChange column found in {file_path}")
            except Exception as e:
                logger.error(f"Error loading ATAC-seq data from {file_path}: {e}")
        
        return atac_data
    
    def load_cluster_genes(self) -> Dict[str, List[str]]:
        """
        Load RNA-seq cluster genes from text files.
        
        Returns:
            Dictionary mapping cluster names to lists of gene IDs
        """
        logger.info("Loading RNA-seq cluster genes...")
        
        cluster_genes = {}
        for file_path in self.rna_cluster_results_dir.glob("cluster_*_genes.txt"):
            cluster_name = file_path.stem.replace('_genes', '')
            try:
                with open(file_path, 'r') as f:
                    lines = f.readlines()
                    # Skip header lines and empty lines
                    genes = []
                    for line in lines:
                        line = line.strip()
                        if line and line != 'Gene' and line != 'GeneID' and not line.startswith('#'):
                            genes.append(line)
                    
                cluster_genes[cluster_name] = genes
                logger.info(f"Loaded {len(genes)} genes for {cluster_name}")
            except Exception as e:
                logger.error(f"Error loading genes from {file_path}: {e}")
        
        return cluster_genes
    
    def load_cluster_peaks(self) -> Dict[str, List[str]]:
        """
        Load ATAC-seq cluster peaks from text files.
        
        Returns:
            Dictionary mapping cluster names to lists of peak IDs
        """
        logger.info("Loading ATAC-seq cluster peaks...")
        
        cluster_peaks = {}
        for file_path in self.atac_cluster_results_dir.glob("cluster_*_peaks.txt"):
            cluster_name = file_path.stem.replace('_peaks', '')
            try:
                with open(file_path, 'r') as f:
                    lines = f.readlines()
                    # Skip header lines and empty lines
                    peaks = []
                    for line in lines:
                        line = line.strip()
                        if line and line != 'Peaks' and line != 'Peak' and not line.startswith('#'):
                            peaks.append(line)
                    
                cluster_peaks[cluster_name] = peaks
                logger.info(f"Loaded {len(peaks)} peaks for {cluster_name}")
            except Exception as e:
                logger.error(f"Error loading peaks from {file_path}: {e}")
        
        return cluster_peaks
    
    def load_peaks_annotation(self) -> Dict[str, pd.DataFrame]:
        """
        Load peaks annotation data from CSV files.
        
        Returns:
            Dictionary mapping cluster names to annotation DataFrames
        """
        logger.info("Loading peaks annotation data...")
        
        peaks_annotation = {}
        for file_path in self.peaks_annotation_dir.glob("*_peaks_annotation_results.csv"):
            cluster_name = file_path.stem.replace('_peaks_annotation_results', '')
            try:
                df = pd.read_csv(file_path)
                # Check for required columns
                required_cols = ['geneId', 'Peak', 'distanceToTSS', 'peak_type']
                if all(col in df.columns for col in required_cols):
                    peaks_annotation[cluster_name] = df[required_cols]
                    logger.info(f"Loaded annotation for {cluster_name}: {len(df)} gene-peak mappings")
                else:
                    logger.warning(f"Missing required columns in {file_path}")
            except Exception as e:
                logger.error(f"Error loading annotation from {file_path}: {e}")
        
        return peaks_annotation
    
    def load_overlapping_genes(self) -> Dict[str, List[str]]:
        """
        Load overlapping genes from text files with improved parsing.
        
        Returns:
            Dictionary mapping cluster pair names to lists of gene IDs
        """
        logger.info("Loading overlapping genes...")
        
        overlapping_genes = {}
        for file_path in self.overlapping_genes_dir.glob("*.txt"):
            if "summary" in file_path.stem.lower():
                continue
            
            cluster_pair = file_path.stem
            try:
                with open(file_path, 'r') as f:
                    lines = f.readlines()
                    
                # Parse genes, skipping headers and empty lines
                genes = []
                for line in lines:
                    line = line.strip()
                    # Skip empty lines, headers, and comment lines
                    if (line and 
                        line != 'Gene' and 
                        line != 'GeneID' and 
                        not line.startswith('#') and
                        not line.startswith('gene') and
                        len(line) > 1):  # Ensure it's not just a single character
                        genes.append(line)
                
                overlapping_genes[cluster_pair] = genes
                logger.info(f"Loaded {len(genes)} overlapping genes for {cluster_pair} from {file_path.name}")
                
                # Debug: show first few genes
                if genes:
                    logger.info(f"First 5 genes: {genes[:5]}")
                
            except Exception as e:
                logger.error(f"Error loading overlapping genes from {file_path}: {e}")
        
        return overlapping_genes
    
    def get_sample_mapping(self) -> Tuple[Dict[str, str], Dict[str, str]]:
        """
        Create mapping between RNA-seq and ATAC-seq sample names.
        
        Returns:
            Tuple of (sample_mapping, reverse_mapping)
        """
        # Map RNA-seq samples to ATAC-seq samples based on naming convention
        sample_mapping = {
            'dko_heat_dko_NHS': 'dko_HS_dko_CK',
            'hsfa_heat_hsfa_NHS': 'hsfa_HS_hsfa_CK',
            'hsfb_heat_hsfb_NHS': 'hsfb_HS_hsfb_CK',
            'Tak1_heat_Tak1_NHS': 'Tak1_HS_Tak1_CK'
        }
        
        # Create reverse mapping
        reverse_mapping = {v: k for k, v in sample_mapping.items()}
        
        return sample_mapping, reverse_mapping
    
    def get_ordered_sample_names(self) -> Tuple[List[str], List[str]]:
        """
        Get ordered sample names to ensure consistent ordering between RNA and ATAC data.
        
        Returns:
            Tuple of (ordered_rna_samples, ordered_atac_samples)
        """
        # Define the desired order: Tak1 -> hsfa -> hsfb -> dko
        ordered_rna_samples = [
            'Tak1_heat_Tak1_NHS',
            'hsfa_heat_hsfa_NHS', 
            'hsfb_heat_hsfb_NHS',
            'dko_heat_dko_NHS'
        ]
        
        # Get corresponding ATAC samples in the same order
        sample_mapping, _ = self.get_sample_mapping()
        ordered_atac_samples = [sample_mapping[rna_sample] for rna_sample in ordered_rna_samples]
        
        return ordered_rna_samples, ordered_atac_samples
    
    def extract_expression_data(self):
        """
        Extract expression data for each cluster and save to individual folders.
        """
        logger.info("Extracting expression data for each cluster...")
        
        # Create directories for extracted data
        extracted_dir = self.output_dir / "extracted_data"
        (extracted_dir / "RNA_clusters").mkdir(exist_ok=True)
        (extracted_dir / "ATAC_clusters").mkdir(exist_ok=True)
        
        # Extract RNA-seq expression data for each cluster
        for cluster_name, genes in self.rna_cluster_genes.items():
            cluster_data = {}
            
            for sample_name, deseq2_df in self.rna_deseq2_data.items():
                # Get log2FoldChange values for genes in this cluster
                cluster_genes_in_data = [gene for gene in genes if gene in deseq2_df.index]
                if cluster_genes_in_data:
                    cluster_data[sample_name] = deseq2_df.loc[cluster_genes_in_data, 'log2FoldChange']
            
            if cluster_data:
                # Combine all samples for this cluster
                combined_df = pd.DataFrame(cluster_data)
                combined_df.index.name = 'Gene'
                
                # Save to file
                output_file = extracted_dir / "RNA_clusters" / f"{cluster_name}_expression.csv"
                combined_df.to_csv(output_file)
                logger.info(f"Saved RNA expression data for {cluster_name}: {len(combined_df)} genes")
        
        # Extract ATAC-seq accessibility data for each cluster
        for cluster_name, peaks in self.atac_cluster_peaks.items():
            cluster_data = {}
            
            for sample_name, deseq2_df in self.atac_deseq2_data.items():
                # Get log2FoldChange values for peaks in this cluster
                cluster_peaks_in_data = [peak for peak in peaks if peak in deseq2_df.index]
                if cluster_peaks_in_data:
                    cluster_data[sample_name] = deseq2_df.loc[cluster_peaks_in_data, 'log2FoldChange']
            
            if cluster_data:
                # Combine all samples for this cluster
                combined_df = pd.DataFrame(cluster_data)
                combined_df.index.name = 'Peak'
                
                # Save to file
                output_file = extracted_dir / "ATAC_clusters" / f"{cluster_name}_accessibility.csv"
                combined_df.to_csv(output_file)
                logger.info(f"Saved ATAC accessibility data for {cluster_name}: {len(combined_df)} peaks")
    
    def process_peaks_to_genes(self, peak_strategy: str = 'closest_to_tss'):
        """
        Process peaks to genes mapping using specified strategy.
        
        Args:
            peak_strategy: Strategy for selecting peaks ('closest_to_tss', 'proximal_median', 
                         'distal_median', 'all_median')
        """
        logger.info(f"Processing peaks to genes using strategy: {peak_strategy}")
        
        peaks_to_genes_dir = self.output_dir / "peaks_to_genes"
        (peaks_to_genes_dir / peak_strategy).mkdir(exist_ok=True)
        
        # Process each cluster
        for cluster_name in self.atac_cluster_peaks.keys():
            if cluster_name not in self.peaks_annotation:
                logger.warning(f"No annotation data for {cluster_name}")
                continue
            
            annotation_df = self.peaks_annotation[cluster_name]
            peaks_data = {}
            
            # Group peaks by gene
            for _, row in annotation_df.iterrows():
                gene_id = row['geneId']
                peak_id = row['Peak']
                distance = abs(row['distanceToTSS'])
                peak_type = row['peak_type']
                
                if gene_id not in peaks_data:
                    peaks_data[gene_id] = []
                
                peaks_data[gene_id].append({
                    'peak_id': peak_id,
                    'distance': distance,
                    'peak_type': peak_type
                })
            
            # Process each gene
            processed_genes = {}
            for gene_id, peaks_info in peaks_data.items():
                if not peaks_info:
                    continue
                
                # Select peak based on strategy
                selected_peak = self.select_peak_by_strategy(peaks_info, peak_strategy)
                if selected_peak:
                    processed_genes[gene_id] = selected_peak
            
            # Save processed genes
            if processed_genes:
                output_file = peaks_to_genes_dir / peak_strategy / f"{cluster_name}_processed_genes.csv"
                with open(output_file, 'w') as f:
                    f.write("GeneID,SelectedPeak,DistanceToTSS,PeakType\n")
                    for gene_id, peak_info in processed_genes.items():
                        f.write(f"{gene_id},{peak_info['peak_id']},{peak_info['distance']},{peak_info['peak_type']}\n")
                
                logger.info(f"Saved processed genes for {cluster_name}: {len(processed_genes)} genes")
    
    def select_peak_by_strategy(self, peaks_info: List[Dict], strategy: str) -> Optional[Dict]:
        """
        Select peak based on specified strategy.
        
        Args:
            peaks_info: List of peak information dictionaries
            strategy: Peak selection strategy
            
        Returns:
            Selected peak information or None
        """
        if not peaks_info:
            return None
        
        if strategy == 'closest_to_tss':
            # Select peak closest to TSS
            return min(peaks_info, key=lambda x: x['distance'])
        
        elif strategy == 'proximal_median':
            # Select median of proximal peaks
            proximal_peaks = [p for p in peaks_info if p['peak_type'] == 'proximal']
            if proximal_peaks:
                proximal_peaks.sort(key=lambda x: x['distance'])
                median_idx = len(proximal_peaks) // 2
                return proximal_peaks[median_idx]
            else:
                logger.warning("No proximal peaks found, using closest peak")
                return min(peaks_info, key=lambda x: x['distance'])
        
        elif strategy == 'distal_median':
            # Select median of distal peaks
            distal_peaks = [p for p in peaks_info if p['peak_type'] == 'distal']
            if distal_peaks:
                distal_peaks.sort(key=lambda x: x['distance'])
                median_idx = len(distal_peaks) // 2
                return distal_peaks[median_idx]
            else:
                logger.warning("No distal peaks found, using closest peak")
                return min(peaks_info, key=lambda x: x['distance'])
        
        elif strategy == 'all_median':
            # Select median of all peaks
            all_peaks = sorted(peaks_info, key=lambda x: x['distance'])
            median_idx = len(all_peaks) // 2
            return all_peaks[median_idx]
        
        else:
            raise ValueError(f"Unknown strategy: {strategy}")
    
    def calculate_correlations(self, peak_strategy: str = 'closest_to_tss'):
        """
        Calculate correlations between RNA-seq and ATAC-seq data.
        
        Args:
            peak_strategy: Peak selection strategy used
        """
        logger.info(f"Calculating correlations using peak strategy: {peak_strategy}")
        
        # Get sample mapping
        sample_mapping, reverse_mapping = self.get_sample_mapping()
        
        # Initialize correlation matrices with expected 4x4 structure
        correlation_matrix = pd.DataFrame(index=self.expected_rna_clusters, 
                                        columns=self.expected_atac_clusters)
        pvalue_matrix = pd.DataFrame(index=self.expected_rna_clusters, 
                                   columns=self.expected_atac_clusters)
        
        # Initialize individual correlations storage
        self.individual_correlations = {}
        
        # Store all p-values for FDR correction
        all_pvalues = []
        pvalue_positions = []  # Track positions for FDR assignment
        
        # Process each expected cluster pair
        for rna_cluster in self.expected_rna_clusters:
            for atac_cluster in self.expected_atac_clusters:
                cluster_pair = f"{rna_cluster}_vs_{atac_cluster}_overlapping_genes"
                
                # Check if we have overlapping genes for this pair
                if cluster_pair not in self.overlapping_genes:
                    logger.warning(f"No overlapping genes file found for {cluster_pair}")
                    correlation_matrix.loc[rna_cluster, atac_cluster] = np.nan
                    pvalue_matrix.loc[rna_cluster, atac_cluster] = np.nan
                    continue
                
                overlapping_genes = self.overlapping_genes[cluster_pair]
                
                if not overlapping_genes:
                    logger.warning(f"No overlapping genes for {cluster_pair}")
                    correlation_matrix.loc[rna_cluster, atac_cluster] = np.nan
                    pvalue_matrix.loc[rna_cluster, atac_cluster] = np.nan
                    continue
                
                logger.info(f"Processing {cluster_pair}: {len(overlapping_genes)} overlapping genes")
                
                # Get RNA-seq expression data for overlapping genes
                rna_values = []
                atac_values = []
                gene_data = []  # Store individual gene data for output
                
                # Get ordered sample names to ensure consistent ordering
                ordered_rna_samples, ordered_atac_samples = self.get_ordered_sample_names()
                
                for gene in overlapping_genes:
                    # Get RNA-seq expression from all samples in correct order
                    gene_rna_values = []
                    for sample_name in ordered_rna_samples:
                        if gene in self.rna_deseq2_data[sample_name].index:
                            gene_rna_values.append(self.rna_deseq2_data[sample_name].loc[gene, 'log2FoldChange'])
                        else:
                            gene_rna_values.append(np.nan)
                    
                    # Get ATAC-seq accessibility for this gene in correct order
                    gene_atac_values = []
                    if atac_cluster in self.peaks_annotation:
                        annotation_df = self.peaks_annotation[atac_cluster]
                        gene_peaks = annotation_df[annotation_df['geneId'] == gene]
                        
                        if not gene_peaks.empty:
                            # Select peak based on strategy
                            peaks_info = []
                            for _, row in gene_peaks.iterrows():
                                peaks_info.append({
                                    'peak_id': row['Peak'],
                                    'distance': abs(row['distanceToTSS']),
                                    'peak_type': row['peak_type']
                                })
                            
                            selected_peak = self.select_peak_by_strategy(peaks_info, peak_strategy)
                            if selected_peak:
                                # Get ATAC-seq accessibility from all samples in correct order
                                for sample_name in ordered_atac_samples:
                                    if selected_peak['peak_id'] in self.atac_deseq2_data[sample_name].index:
                                        gene_atac_values.append(self.atac_deseq2_data[sample_name].loc[selected_peak['peak_id'], 'log2FoldChange'])
                                    else:
                                        gene_atac_values.append(np.nan)
                            else:
                                gene_atac_values = [np.nan] * len(ordered_atac_samples)
                        else:
                            gene_atac_values = [np.nan] * len(ordered_atac_samples)
                    else:
                        gene_atac_values = [np.nan] * len(ordered_atac_samples)
                    
                    # Store gene data for individual correlation output
                    gene_data.append({
                        'Gene': gene,
                        'RNA_values': gene_rna_values,
                        'ATAC_values': gene_atac_values,
                        'RNA_samples': ordered_rna_samples,
                        'ATAC_samples': ordered_atac_samples
                    })
                    
                    # Only include genes with valid data
                    if not (np.isnan(gene_rna_values).all() or np.isnan(gene_atac_values).all()):
                        rna_values.append(gene_rna_values)
                        atac_values.append(gene_atac_values)
                
                # Calculate correlation
                if len(rna_values) >= 3:
                    # Flatten the data for correlation calculation
                    flat_rna = [val for sublist in rna_values for val in sublist if not np.isnan(val)]
                    flat_atac = [val for sublist in atac_values for val in sublist if not np.isnan(val)]
                    
                    if len(flat_rna) >= 3 and len(flat_atac) >= 3:
                        try:
                            if self.correlation_method == 'pearson':
                                correlation, p_value = pearsonr(flat_rna, flat_atac)
                            else:
                                correlation, p_value = spearmanr(flat_rna, flat_atac)
                            
                            correlation_matrix.loc[rna_cluster, atac_cluster] = correlation
                            pvalue_matrix.loc[rna_cluster, atac_cluster] = p_value
                            
                            # Store p-value for FDR correction
                            all_pvalues.append(p_value)
                            pvalue_positions.append((rna_cluster, atac_cluster))
                            
                            logger.info(f"{cluster_pair}: correlation={correlation:.3f}, p-value={p_value:.3e}")
                            
                            # Store individual correlation data
                            self.individual_correlations[cluster_pair] = {
                                'correlation': correlation,
                                'p_value': p_value,
                                'gene_data': gene_data,
                                'flat_rna': flat_rna,
                                'flat_atac': flat_atac,
                                'n_genes': len(overlapping_genes),
                                'n_valid_genes': len(rna_values)
                            }
                            
                        except Exception as e:
                            logger.error(f"Error calculating correlation for {cluster_pair}: {e}")
                            correlation_matrix.loc[rna_cluster, atac_cluster] = np.nan
                            pvalue_matrix.loc[rna_cluster, atac_cluster] = np.nan
                    else:
                        correlation_matrix.loc[rna_cluster, atac_cluster] = np.nan
                        pvalue_matrix.loc[rna_cluster, atac_cluster] = np.nan
                else:
                    correlation_matrix.loc[rna_cluster, atac_cluster] = np.nan
                    pvalue_matrix.loc[rna_cluster, atac_cluster] = np.nan
        
        # Perform FDR correction
        if all_pvalues:
            logger.info(f"Performing FDR correction using {self.fdr_method} method...")
            try:
                # Perform multiple testing correction
                rejected, pvals_corrected, alphacSidak, alphacBonf = multipletests(
                    all_pvalues, 
                    alpha=0.05, 
                    method=self.fdr_method
                )
                
                # Create FDR matrix
                fdr_matrix = pd.DataFrame(index=self.expected_rna_clusters, 
                                         columns=self.expected_atac_clusters)
                fdr_matrix = fdr_matrix.fillna(np.nan)
                
                # Assign corrected p-values back to matrix
                for i, (rna_cluster, atac_cluster) in enumerate(pvalue_positions):
                    fdr_matrix.loc[rna_cluster, atac_cluster] = pvals_corrected[i]
                
                self.fdr_matrix = fdr_matrix
                logger.info(f"FDR correction completed. {sum(rejected)} significant correlations after correction.")
                
            except Exception as e:
                logger.error(f"Error performing FDR correction: {e}")
                self.fdr_matrix = pvalue_matrix.copy()  # Use original p-values if FDR fails
        
        self.correlation_matrix = correlation_matrix
        self.pvalue_matrix = pvalue_matrix
        
        return correlation_matrix, pvalue_matrix

    def create_individual_scatter_plots(self, peak_strategy: str):
        """
        Create individual scatter plots for each cluster pair overlapping results.
        
        Args:
            peak_strategy: Peak selection strategy used
        """
        logger.info("Creating individual scatter plots for each cluster pair...")
        
        scatter_dir = self.output_dir / "scatter_plots"
        (scatter_dir / peak_strategy).mkdir(exist_ok=True)
        
        # Get ordered sample names for consistent labeling
        ordered_rna_samples, ordered_atac_samples = self.get_ordered_sample_names()
        
        for cluster_pair, data in self.individual_correlations.items():
            if not data['gene_data']:
                continue
            
            # Create figure with subplots for each sample
            n_samples = len(ordered_rna_samples)
            fig, axes = plt.subplots(2, 2, figsize=(15, 12))
            axes = axes.flatten()
            
            # Plot each sample
            for i, (rna_sample, atac_sample) in enumerate(zip(ordered_rna_samples, ordered_atac_samples)):
                if i >= len(axes):
                    break
                
                ax = axes[i]
                
                # Extract data for this sample
                x_values = []
                y_values = []
                gene_labels = []
                
                for gene_info in data['gene_data']:
                    gene = gene_info['Gene']
                    rna_idx = ordered_rna_samples.index(rna_sample)
                    atac_idx = ordered_atac_samples.index(atac_sample)
                    
                    if (rna_idx < len(gene_info['RNA_values']) and 
                        atac_idx < len(gene_info['ATAC_values'])):
                        
                        rna_val = gene_info['RNA_values'][rna_idx]
                        atac_val = gene_info['ATAC_values'][atac_idx]
                        
                        if not (np.isnan(rna_val) or np.isnan(atac_val)):
                            x_values.append(rna_val)
                            y_values.append(atac_val)
                            gene_labels.append(gene)
                
                if x_values and y_values:
                    # Create scatter plot
                    scatter = ax.scatter(x_values, y_values, alpha=0.7, s=50, edgecolors='black', linewidth=0.5)
                    
                    # Add trend line
                    if len(x_values) > 1:
                        z = np.polyfit(x_values, y_values, 1)
                        p = np.poly1d(z)
                        ax.plot(x_values, p(x_values), "r--", alpha=0.8, linewidth=2)
                    
                    # Calculate correlation for this sample
                    if len(x_values) >= 3:
                        try:
                            if self.correlation_method == 'pearson':
                                corr, p_val = pearsonr(x_values, y_values)
                            else:
                                corr, p_val = spearmanr(x_values, y_values)
                            
                            # Add correlation info to plot
                            ax.text(0.05, 0.95, 
                                   f'r = {corr:.3f}\np = {p_val:.3e}', 
                                   transform=ax.transAxes, 
                                   bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
                                   fontsize=10, verticalalignment='top')
                        except:
                            pass
                    
                    # Customize plot
                    ax.set_xlabel(f'RNA-seq log2FoldChange ({rna_sample})', fontsize=10)
                    ax.set_ylabel(f'ATAC-seq log2FoldChange ({atac_sample})', fontsize=10)
                    ax.set_title(f'Sample: {rna_sample} vs {atac_sample}', fontsize=12, fontweight='bold')
                    ax.grid(True, alpha=0.3)
                    ax.axhline(y=0, color='black', linestyle='-', alpha=0.3)
                    ax.axvline(x=0, color='black', linestyle='-', alpha=0.3)
                    
                    # Add gene labels for top genes (if not too many)
                    if len(gene_labels) <= 20:  # Only label if reasonable number of genes
                        for j, gene in enumerate(gene_labels):
                            ax.annotate(gene, (x_values[j], y_values[j]), 
                                       xytext=(5, 5), textcoords='offset points',
                                       fontsize=8, alpha=0.7)
                else:
                    ax.text(0.5, 0.5, 'No valid data', 
                           transform=ax.transAxes, ha='center', va='center',
                           fontsize=12, style='italic')
                    ax.set_title(f'Sample: {rna_sample} vs {atac_sample}', fontsize=12, fontweight='bold')
            
            # Remove empty subplots if any
            for i in range(len(data['gene_data']), len(axes)):
                fig.delaxes(axes[i])
            
            # Add overall title
            fig.suptitle(f'Individual Scatter Plots: {cluster_pair}\nPeak Strategy: {peak_strategy}', 
                         fontsize=16, fontweight='bold')
            
            plt.tight_layout()
            
            # Save plot
            output_file = scatter_dir / peak_strategy / f"{cluster_pair}_scatter_plots.pdf"
            fig.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close(fig)
            
            logger.info(f"Saved scatter plots for {cluster_pair} to: {output_file}")
        
        # Create combined scatter plot with all samples
        self.create_combined_scatter_plot(peak_strategy, scatter_dir)

    def create_combined_scatter_plot(self, peak_strategy: str, scatter_dir: Path):
        """
        Create a combined scatter plot showing all samples together.
        
        Args:
            peak_strategy: Peak selection strategy used
            scatter_dir: Directory for scatter plots
        """
        logger.info("Creating combined scatter plot...")
        
        # Get ordered sample names
        ordered_rna_samples, ordered_atac_samples = self.get_ordered_sample_names()
        
        # Create figure
        fig, ax = plt.subplots(figsize=(12, 10))
        
        # Colors for different samples
        colors = plt.cm.Set1(np.linspace(0, 1, len(ordered_rna_samples)))
        
        # Plot each sample with different color
        for i, (rna_sample, atac_sample) in enumerate(zip(ordered_rna_samples, ordered_atac_samples)):
            all_x = []
            all_y = []
            
            # Collect data from all cluster pairs for this sample
            for cluster_pair, data in self.individual_correlations.items():
                if not data['gene_data']:
                    continue
                
                rna_idx = ordered_rna_samples.index(rna_sample)
                atac_idx = ordered_atac_samples.index(atac_sample)
                
                for gene_info in data['gene_data']:
                    if (rna_idx < len(gene_info['RNA_values']) and 
                        atac_idx < len(gene_info['ATAC_values'])):
                        
                        rna_val = gene_info['RNA_values'][rna_idx]
                        atac_val = gene_info['ATAC_values'][atac_idx]
                        
                        if not (np.isnan(rna_val) or np.isnan(atac_val)):
                            all_x.append(rna_val)
                            all_y.append(atac_val)
            
            if all_x and all_y:
                # Create scatter plot for this sample
                ax.scatter(all_x, all_y, alpha=0.6, s=40, 
                          color=colors[i], label=f'{rna_sample} vs {atac_sample}',
                          edgecolors='black', linewidth=0.3)
                
                # Add trend line
                if len(all_x) > 1:
                    z = np.polyfit(all_x, all_y, 1)
                    p = np.poly1d(z)
                    ax.plot(all_x, p(all_x), color=colors[i], alpha=0.8, 
                           linewidth=2, linestyle='--')
        
        # Customize plot
        ax.set_xlabel('RNA-seq log2FoldChange', fontsize=14, fontweight='bold')
        ax.set_ylabel('ATAC-seq log2FoldChange', fontsize=14, fontweight='bold')
        ax.set_title(f'Combined RNA-seq vs ATAC-seq Correlation\nPeak Strategy: {peak_strategy}', 
                     fontsize=16, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.axhline(y=0, color='black', linestyle='-', alpha=0.5)
        ax.axvline(x=0, color='black', linestyle='-', alpha=0.5)
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        
        # Calculate overall correlation
        all_rna_values = []
        all_atac_values = []
        for cluster_pair, data in self.individual_correlations.items():
            if data['flat_rna'] and data['flat_atac']:
                all_rna_values.extend(data['flat_rna'])
                all_atac_values.extend(data['flat_atac'])
        
        if all_rna_values and all_atac_values:
            try:
                if self.correlation_method == 'pearson':
                    overall_corr, overall_p = pearsonr(all_rna_values, all_atac_values)
                else:
                    overall_corr, overall_p = spearmanr(all_rna_values, all_atac_values)
                
                # Add overall correlation info
                ax.text(0.05, 0.95, 
                       f'Overall {self.correlation_method.capitalize()}:\nr = {overall_corr:.3f}\np = {overall_p:.3e}', 
                       transform=ax.transAxes, 
                       bbox=dict(boxstyle="round,pad=0.5", facecolor="white", alpha=0.9),
                       fontsize=12, verticalalignment='top', fontweight='bold')
            except:
                pass
        
        plt.tight_layout()
        
        # Save combined plot
        output_file = scatter_dir / peak_strategy / f"combined_scatter_plot_{peak_strategy}.pdf"
        fig.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close(fig)
        
        logger.info(f"Saved combined scatter plot to: {output_file}")

    def get_significance_symbol(self, p_value: float, fdr_value: float = None) -> str:
        """
        Get significance symbol based on p-value and FDR value.
        
        Args:
            p_value: P-value
            fdr_value: FDR-corrected p-value (optional)
            
        Returns:
            Significance symbol
        """
        if np.isnan(p_value):
            return ""
        
        # Use FDR value if available, otherwise use p-value
        test_value = fdr_value if fdr_value is not None and not np.isnan(fdr_value) else p_value
        
        if test_value < 0.001:
            return "***"
        elif test_value < 0.01:
            return "**"
        elif test_value < 0.05:
            return "*"
        else:
            return ""

    def create_correlation_heatmap(self) -> plt.Figure:
        """
        Create correlation heatmap with significance indicators using both p-values and FDR.
        
        Returns:
            Matplotlib figure object
        """
        if self.correlation_matrix is None or self.pvalue_matrix is None:
            logger.error("No correlation data available for heatmap")
            return None
        
        # Create figure
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))
        
        # Convert correlation matrix to numeric, replacing NaN with 0 for visualization
        corr_matrix_numeric = self.correlation_matrix.astype(float)
        corr_matrix_numeric = corr_matrix_numeric.fillna(0)
        
        # Create correlation heatmap with p-value significance
        sns.heatmap(corr_matrix_numeric, 
                   annot=True, 
                   fmt='.3f',
                   cmap='RdBu_r',
                   center=0,
                   vmin=-1, 
                   vmax=1,
                   cbar_kws={'label': f'{self.correlation_method.capitalize()} Correlation'},
                   ax=ax1)
        
        # Add p-value significance indicators
        for i in range(len(self.correlation_matrix.index)):
            for j in range(len(self.correlation_matrix.columns)):
                p_value = self.pvalue_matrix.iloc[i, j]
                significance_symbol = self.get_significance_symbol(p_value)
                if significance_symbol:
                    ax1.text(j + 0.5, i + 0.5, significance_symbol, 
                            ha='center', va='center', 
                            fontsize=12, fontweight='bold', color='black')
        
        ax1.set_title(f'RNA-seq vs ATAC-seq Cluster Correlation Heatmap\n({self.correlation_method.capitalize()}) - P-values', 
                      fontsize=14, fontweight='bold')
        ax1.set_xlabel('ATAC-seq Clusters', fontsize=12)
        ax1.set_ylabel('RNA-seq Clusters', fontsize=12)
        
        # Create correlation heatmap with FDR significance (if available)
        if self.fdr_matrix is not None:
            sns.heatmap(corr_matrix_numeric, 
                       annot=True, 
                       fmt='.3f',
                       cmap='RdBu_r',
                       center=0,
                       vmin=-1, 
                       vmax=1,
                       cbar_kws={'label': f'{self.correlation_method.capitalize()} Correlation'},
                       ax=ax2)
            
            # Add FDR significance indicators
            for i in range(len(self.correlation_matrix.index)):
                for j in range(len(self.correlation_matrix.columns)):
                    fdr_value = self.fdr_matrix.iloc[i, j]
                    p_value = self.pvalue_matrix.iloc[i, j]
                    significance_symbol = self.get_significance_symbol(p_value, fdr_value)
                    if significance_symbol:
                        ax2.text(j + 0.5, i + 0.5, significance_symbol, 
                                ha='center', va='center', 
                                fontsize=12, fontweight='bold', color='black')
            
            ax2.set_title(f'RNA-seq vs ATAC-seq Cluster Correlation Heatmap\n({self.correlation_method.capitalize()}) - FDR Corrected', 
                          fontsize=14, fontweight='bold')
            ax2.set_xlabel('ATAC-seq Clusters', fontsize=12)
            ax2.set_ylabel('RNA-seq Clusters', fontsize=12)
        else:
            ax2.text(0.5, 0.5, 'FDR correction not available', 
                     transform=ax2.transAxes, ha='center', va='center',
                     fontsize=14, style='italic')
            ax2.set_title('FDR Correction Not Available', fontsize=14, fontweight='bold')
        
        plt.tight_layout()
        return fig
    
    def save_individual_correlations(self, peak_strategy: str):
        """
        Save individual correlation data for each cluster pair.
        
        Args:
            peak_strategy: Peak selection strategy used
        """
        logger.info("Saving individual correlation data...")
        
        # Create directory for individual correlations
        individual_dir = self.output_dir / "individual_correlations" / peak_strategy
        individual_dir.mkdir(parents=True, exist_ok=True)
        
        timestamp = pd.Timestamp.now().strftime("%Y%m%d_%H%M%S")
        
        # Save individual correlation data for each cluster pair
        for cluster_pair, data in self.individual_correlations.items():
            if not data['gene_data']:
                continue
                
            # Create filename for this cluster pair
            cluster_pair_name = cluster_pair.replace(' vs ', '_vs_')
            data_file = individual_dir / f"{cluster_pair_name}_correlation_data.csv"
            
            # Prepare data for saving
            rows = []
            for gene_info in data['gene_data']:
                gene = gene_info['Gene']
                rna_values = gene_info['RNA_values']
                atac_values = gene_info['ATAC_values']
                
                # Get ordered sample names
                ordered_rna_samples, ordered_atac_samples = self.get_ordered_sample_names()
                
                # Create row for each gene
                # Extract cluster names from cluster_pair format: "cluster1_vs_cluster2_overlapping_genes"
                cluster_parts = cluster_pair.replace('_overlapping_genes', '').split('_vs_')
                if len(cluster_parts) >= 2:
                    rna_cluster = cluster_parts[0]
                    atac_cluster = cluster_parts[1]
                else:
                    rna_cluster = "unknown"
                    atac_cluster = "unknown"
                
                row = {
                    'Gene': gene,
                    'RNA_cluster': rna_cluster,
                    'ATAC_cluster': atac_cluster
                }
                
                # Add RNA values for each sample
                for i, sample in enumerate(ordered_rna_samples):
                    if i < len(rna_values):
                        row[f'RNA_{sample}'] = rna_values[i]
                    else:
                        row[f'RNA_{sample}'] = np.nan
                
                # Add ATAC values for each sample
                for i, sample in enumerate(ordered_atac_samples):
                    if i < len(atac_values):
                        row[f'ATAC_{sample}'] = atac_values[i]
                    else:
                        row[f'ATAC_{sample}'] = np.nan
                
                rows.append(row)
            
            # Save to CSV
            if rows:
                df = pd.DataFrame(rows)
                df.to_csv(data_file, index=False)
                logger.info(f"Saved individual correlation data for {cluster_pair} to: {data_file}")
        
        # Save summary of individual correlations
        summary_file = individual_dir / f"individual_correlations_summary_{timestamp}.csv"
        summary_rows = []
        
        for cluster_pair, data in self.individual_correlations.items():
            # Extract cluster names from cluster_pair format: "cluster1_vs_cluster2_overlapping_genes"
            cluster_parts = cluster_pair.replace('_overlapping_genes', '').split('_vs_')
            if len(cluster_parts) >= 2:
                rna_cluster = cluster_parts[0]
                atac_cluster = cluster_parts[1]
            else:
                rna_cluster = "unknown"
                atac_cluster = "unknown"
            
            summary_rows.append({
                'RNA_cluster': rna_cluster,
                'ATAC_cluster': atac_cluster,
                'correlation': data['correlation'],
                'p_value': data['p_value'],
                'n_genes': data['n_genes'],
                'n_valid_genes': data['n_valid_genes']
            })
        
        if summary_rows:
            summary_df = pd.DataFrame(summary_rows)
            summary_df.to_csv(summary_file, index=False)
            logger.info(f"Saved individual correlations summary to: {summary_file}")
    
    def save_results(self, peak_strategy: str):
        """
        Save analysis results to files.
        
        Args:
            peak_strategy: Peak selection strategy used
        """
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
        
        # Save FDR matrix if available
        if self.fdr_matrix is not None:
            fdr_file = results_dir / f"fdr_matrix_{peak_strategy}_{timestamp}.csv"
            self.fdr_matrix.to_csv(fdr_file)
            logger.info(f"Saved FDR matrix to: {fdr_file}")
        
        # Save heatmap data for R ComplexHeatmap
        heatmap_data = []
        for i in range(len(self.correlation_matrix.index)):
            for j in range(len(self.correlation_matrix.columns)):
                rna_cluster = self.correlation_matrix.index[i]
                atac_cluster = self.correlation_matrix.columns[j]
                correlation = self.correlation_matrix.iloc[i, j]
                p_value = self.pvalue_matrix.iloc[i, j]
                fdr_value = self.fdr_matrix.iloc[i, j] if self.fdr_matrix is not None else np.nan
                significance_p = self.get_significance_symbol(p_value)
                significance_fdr = self.get_significance_symbol(p_value, fdr_value)
                
                heatmap_data.append({
                    'RNA_cluster': rna_cluster,
                    'ATAC_cluster': atac_cluster,
                    'correlation': correlation,
                    'p_value': p_value,
                    'fdr_value': fdr_value,
                    'significance_p': significance_p,
                    'significance_fdr': significance_fdr
                })
        
        # Save as CSV for R
        r_data_file = results_dir / f"heatmap_data_for_R_{peak_strategy}_{timestamp}.csv"
        heatmap_df = pd.DataFrame(heatmap_data)
        heatmap_df.to_csv(r_data_file, index=False)
        logger.info(f"Saved R-compatible data to: {r_data_file}")
        
        # Save analysis summary
        summary_file = results_dir / f"analysis_summary_{peak_strategy}_{timestamp}.txt"
        with open(summary_file, 'w') as f:
            f.write("RNA-seq vs ATAC-seq Cluster Correlation Analysis Summary\n")
            f.write("=" * 60 + "\n\n")
            f.write(f"Analysis Date: {timestamp}\n")
            f.write(f"Correlation Method: {self.correlation_method.capitalize()}\n")
            f.write(f"FDR Correction Method: {self.fdr_method}\n")
            f.write(f"Peak Selection Strategy: {peak_strategy}\n")
            f.write(f"Number of RNA-seq clusters: {len(self.correlation_matrix.index)}\n")
            f.write(f"Number of ATAC-seq clusters: {len(self.correlation_matrix.columns)}\n")
            f.write(f"Expected matrix size: 4x4\n\n")
            
            f.write("Correlation Matrix:\n")
            f.write(self.correlation_matrix.to_string())
            f.write("\n\nP-value Matrix:\n")
            f.write(self.pvalue_matrix.to_string())
            
            if self.fdr_matrix is not None:
                f.write("\n\nFDR-corrected P-value Matrix:\n")
                f.write(self.fdr_matrix.to_string())
            
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
                
                # Count significant correlations (p-values)
                pval_values = pd.to_numeric(self.pvalue_matrix.values.flatten(), errors='coerce')
                valid_pvals = pval_values[~np.isnan(pval_values)]
                significant_pairs_001_p = np.sum(valid_pvals < 0.001)
                significant_pairs_01_p = np.sum(valid_pvals < 0.01)
                significant_pairs_05_p = np.sum(valid_pvals < 0.05)
                
                f.write(f"\nP-value Significance:\n")
                f.write(f"Significant correlations (p < 0.001): {significant_pairs_001_p}\n")
                f.write(f"Significant correlations (p < 0.01): {significant_pairs_01_p}\n")
                f.write(f"Significant correlations (p < 0.05): {significant_pairs_05_p}\n")
                
                # Count significant correlations (FDR)
                if self.fdr_matrix is not None:
                    fdr_values = pd.to_numeric(self.fdr_matrix.values.flatten(), errors='coerce')
                    valid_fdrs = fdr_values[~np.isnan(fdr_values)]
                    significant_pairs_001_fdr = np.sum(valid_fdrs < 0.001)
                    significant_pairs_01_fdr = np.sum(valid_fdrs < 0.01)
                    significant_pairs_05_fdr = np.sum(valid_fdrs < 0.05)
                    
                    f.write(f"\nFDR-corrected Significance:\n")
                    f.write(f"Significant correlations (FDR < 0.001): {significant_pairs_001_fdr}\n")
                    f.write(f"Significant correlations (FDR < 0.01): {significant_pairs_01_fdr}\n")
                    f.write(f"Significant correlations (FDR < 0.05): {significant_pairs_05_fdr}\n")
                
            except Exception as e:
                f.write(f"Error calculating statistics: {e}\n")
        
        logger.info(f"Saved analysis summary to: {summary_file}")
    
    def generate_r_script_for_heatmap(self, peak_strategy: str) -> str:
        """
        Generate R script for creating customizable correlation heatmap.
        
        Args:
            peak_strategy: Peak selection strategy used
            
        Returns:
            R script content as string
        """
        timestamp = pd.Timestamp.now().strftime("%Y%m%d_%H%M%S")
        
        r_script = f"""# RNA-seq vs ATAC-seq Cluster Correlation Heatmap R Script
# Generated on: {timestamp}
# Peak Strategy: {peak_strategy}
# Correlation Method: {self.correlation_method.capitalize()}

# Load required libraries
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(ggplot2)

# Set working directory (modify as needed)
# setwd("path/to/your/output/directory")

# Read correlation data
correlation_data <- read.csv("correlation_results/heatmap_data_for_R_{peak_strategy}_{timestamp}.csv")

# Create correlation matrix
corr_matrix <- correlation_data %>%
  select(RNA_cluster, ATAC_cluster, correlation) %>%
  pivot_wider(names_from = ATAC_cluster, values_from = correlation) %>%
  column_to_rownames("RNA_cluster") %>%
  as.matrix()

# Create p-value matrix
pval_matrix <- correlation_data %>%
  select(RNA_cluster, ATAC_cluster, p_value) %>%
  pivot_wider(names_from = ATAC_cluster, values_from = p_value) %>%
  column_to_rownames("RNA_cluster") %>%
  as.matrix()

# Create significance matrix
sig_matrix <- correlation_data %>%
  select(RNA_cluster, ATAC_cluster, significance) %>%
  pivot_wider(names_from = ATAC_cluster, values_from = significance) %>%
  column_to_rownames("RNA_cluster") %>%
  as.matrix()

# Define color palette
col_fun <- colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B"))

# Create correlation heatmap
ht <- Heatmap(
  corr_matrix,
  name = "{self.correlation_method.capitalize()} Correlation",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "left",
  column_names_side = "bottom",
  row_names_gp = gpar(fontsize = 12, fontface = "bold"),
  column_names_gp = gpar(fontsize = 12, fontface = "bold"),
  row_title = "RNA-seq Clusters",
  column_title = "ATAC-seq Clusters",
  row_title_gp = gpar(fontsize = 14, fontface = "bold"),
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  heatmap_legend_param = list(
    title_position = "topcenter",
    legend_direction = "horizontal",
    legend_width = unit(8, "cm")
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {{
    # Add correlation values
    if (!is.na(corr_matrix[i, j])) {{
      grid.text(
        sprintf("%.3f", corr_matrix[i, j]),
        x, y,
        gp = gpar(fontsize = 10, fontface = "bold")
      )
    }}
    
    # Add significance symbols
    if (sig_matrix[i, j] != "") {{
      grid.text(
        sig_matrix[i, j],
        x, y + height * 0.3,
        gp = gpar(fontsize = 14, fontface = "bold", col = "black")
      )
    }}
  }}
)

# Draw heatmap
pdf("correlation_heatmap_{peak_strategy}_{timestamp}_R.pdf", width = 10, height = 8)
draw(ht, heatmap_legend_side = "bottom")
dev.off()

# Alternative: Create heatmap with ggplot2 for more customization
# Convert data for ggplot2
plot_data <- correlation_data %>%
  mutate(
    correlation_rounded = round(correlation, 3),
    significance_label = ifelse(significance != "", significance, ""),
    fill_color = case_when(
      correlation >= 0.5 ~ "High Positive",
      correlation >= 0.1 ~ "Low Positive", 
      correlation >= -0.1 ~ "No Correlation",
      correlation >= -0.5 ~ "Low Negative",
      TRUE ~ "High Negative"
    )
  )

# Create ggplot2 heatmap
ggplot(plot_data, aes(x = ATAC_cluster, y = RNA_cluster, fill = correlation)) +
  geom_tile(color = "white", size = 0.5) +
  geom_text(aes(label = correlation_rounded), size = 4, fontface = "bold") +
  geom_text(aes(label = significance_label), 
            position = position_nudge(y = 0.3), 
            size = 6, fontface = "bold", color = "black") +
  scale_fill_gradient2(
    low = "#2166AC", 
    mid = "white", 
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1),
    name = "{self.correlation_method.capitalize()} Correlation"
  ) +
  labs(
    title = "RNA-seq vs ATAC-seq Cluster Correlation Heatmap",
    subtitle = paste("Peak Strategy:", "{peak_strategy}", "| Method:", "{self.correlation_method.capitalize()}"),
    x = "ATAC-seq Clusters",
    y = "RNA-seq Clusters"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  ) +
  coord_fixed()

# Save ggplot2 version
ggsave("correlation_heatmap_{peak_strategy}_{timestamp}_ggplot.pdf", 
       width = 10, height = 8, dpi = 300)

# Print summary statistics
cat("\\n=== Correlation Analysis Summary ===\\n")
cat("Peak Strategy:", "{peak_strategy}", "\\n")
cat("Correlation Method:", "{self.correlation_method.capitalize()}", "\\n")
cat("Matrix Size:", nrow(corr_matrix), "x", ncol(corr_matrix), "\\n\\n")

# Summary statistics
valid_correlations <- corr_matrix[!is.na(corr_matrix)]
if (length(valid_correlations) > 0) {{
  cat("Correlation Statistics:\\n")
  cat("Mean:", round(mean(valid_correlations), 3), "\\n")
  cat("Median:", round(median(valid_correlations), 3), "\\n")
  cat("Min:", round(min(valid_correlations), 3), "\\n")
  cat("Max:", round(max(valid_correlations), 3), "\\n\\n")
}}

# Significance summary
valid_pvals <- pval_matrix[!is.na(pval_matrix)]
if (length(valid_pvals) > 0) {{
  cat("Significance Summary:\\n")
  cat("p < 0.001:", sum(valid_pvals < 0.001), "\\n")
  cat("p < 0.01:", sum(valid_pvals < 0.01), "\\n")
  cat("p < 0.05:", sum(valid_pvals < 0.05), "\\n\\n")
}}

cat("Analysis completed! Check the generated PDF files.\\n")
"""
        
        return r_script

    def generate_r_script_for_individual_correlations(self, peak_strategy: str) -> str:
        """
        Generate R script for analyzing individual correlations.
        
        Args:
            peak_strategy: Peak selection strategy used
            
        Returns:
            R script content as string
        """
        timestamp = pd.Timestamp.now().strftime("%Y%m%d_%H%M%S")
        
        r_script = f"""# Individual RNA-ATAC Correlation Analysis R Script
# Generated on: {timestamp}
# Peak Strategy: {peak_strategy}
# Correlation Method: {self.correlation_method.capitalize()}

# Load required libraries
library(ggplot2)
library(dplyr)
library(gridExtra)
library(corrplot)

# Set working directory (modify as needed)
# setwd("path/to/your/output/directory")

# Function to analyze individual cluster pair correlations
analyze_cluster_pair <- function(cluster_pair_name) {{
  # Read individual correlation data
  data_file <- paste0("individual_correlations/{peak_strategy}/", cluster_pair_name, "_correlation_matrix.csv")
  
  if (!file.exists(data_file)) {{
    cat("Data file not found:", data_file, "\\n")
    return(NULL)
  }}
  
  # Read data
  correlation_data <- read.csv(data_file)
  
  # Extract RNA and ATAC values
  rna_values <- unlist(correlation_data$RNA_log2FC)
  atac_values <- unlist(correlation_data$ATAC_log2FC)
  
  # Remove NA values
  valid_data <- data.frame(
    RNA = rna_values[!is.na(rna_values) & !is.na(atac_values)],
    ATAC = atac_values[!is.na(rna_values) & !is.na(atac_values)]
  )
  
  if (nrow(valid_data) < 3) {{
    cat("Insufficient data for", cluster_pair_name, "\\n")
    return(NULL)
  }}
  
  # Calculate correlation
  cor_result <- cor.test(valid_data$RNA, valid_data$ATAC, method = "{self.correlation_method}")
  
  # Create scatter plot
  p1 <- ggplot(valid_data, aes(x = RNA, y = ATAC)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = "red") +
    labs(
      title = paste("RNA vs ATAC Correlation:", cluster_pair_name),
      subtitle = paste("r =", round(cor_result$estimate, 3), 
                      "p =", format.pval(cor_result$p.value, digits = 3)),
      x = "RNA-seq log2FoldChange",
      y = "ATAC-seq log2FoldChange"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.title = element_text(size = 12, face = "bold")
    )
  
  # Create correlation plot
  p2 <- ggplot(valid_data, aes(x = RNA, y = ATAC)) +
    stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
    scale_fill_viridis_c() +
    labs(
      title = "Density Plot",
      x = "RNA-seq log2FoldChange",
      y = "ATAC-seq log2FoldChange"
    ) +
    theme_minimal()
  
  # Combine plots
  combined_plot <- grid.arrange(p1, p2, ncol = 2)
  
  # Save plot
  ggsave(paste0("individual_correlations/{peak_strategy}/", cluster_pair_name, "_analysis.pdf"), 
         combined_plot, width = 12, height = 6, dpi = 300)
  
  # Return correlation results
  return(list(
    correlation = cor_result$estimate,
    p_value = cor_result$p.value,
    n_genes = nrow(valid_data),
    plot = combined_plot
  ))
}}

# Analyze all cluster pairs
cluster_pairs <- c(
  "cluster_1_vs_cluster_1_overlapping_genes",
  "cluster_1_vs_cluster_2_overlapping_genes",
  "cluster_1_vs_cluster_3_overlapping_genes",
  "cluster_1_vs_cluster_4_overlapping_genes",
  "cluster_2_vs_cluster_1_overlapping_genes",
  "cluster_2_vs_cluster_2_overlapping_genes",
  "cluster_2_vs_cluster_3_overlapping_genes",
  "cluster_2_vs_cluster_4_overlapping_genes",
  "cluster_3_vs_cluster_1_overlapping_genes",
  "cluster_3_vs_cluster_2_overlapping_genes",
  "cluster_3_vs_cluster_3_overlapping_genes",
  "cluster_3_vs_cluster_4_overlapping_genes",
  "cluster_4_vs_cluster_1_overlapping_genes",
  "cluster_4_vs_cluster_2_overlapping_genes",
  "cluster_4_vs_cluster_3_overlapping_genes",
  "cluster_4_vs_cluster_4_overlapping_genes"
)

# Process each cluster pair
results <- list()
for (pair in cluster_pairs) {{
  cat("Processing:", pair, "\\n")
  result <- analyze_cluster_pair(pair)
  if (!is.null(result)) {{
    results[[pair]] <- result
  }}
}}

# Create summary table
if (length(results) > 0) {{
  summary_data <- do.call(rbind, lapply(names(results), function(name) {{
    data.frame(
      Cluster_Pair = name,
      Correlation = results[[name]]$correlation,
      P_Value = results[[name]]$p_value,
      N_Genes = results[[name]]$n_genes
    )
  }}))
  
  # Save summary
  write.csv(summary_data, 
            paste0("individual_correlations/{peak_strategy}/individual_correlations_summary.csv"), 
            row.names = FALSE)
  
  # Print summary
  cat("\\n=== Individual Correlations Summary ===\\n")
  print(summary_data)
}}

cat("\\nIndividual correlation analysis completed!\\n")
"""
        
        return r_script

    def save_r_scripts(self, peak_strategy: str):
        """
        Save R scripts for further customization.
        
        Args:
            peak_strategy: Peak selection strategy used
        """
        logger.info("Generating R scripts for customization...")
        
        r_scripts_dir = self.output_dir / "R_scripts"
        timestamp = pd.Timestamp.now().strftime("%Y%m%d_%H%M%S")
        
        # Generate and save heatmap R script
        heatmap_r_script = self.generate_r_script_for_heatmap(peak_strategy)
        heatmap_r_file = r_scripts_dir / f"correlation_heatmap_{peak_strategy}_{timestamp}.R"
        with open(heatmap_r_file, 'w') as f:
            f.write(heatmap_r_script)
        logger.info(f"Saved heatmap R script to: {heatmap_r_file}")
        
        # Generate and save individual correlations R script
        individual_r_script = self.generate_r_script_for_individual_correlations(peak_strategy)
        individual_r_file = r_scripts_dir / f"individual_correlations_{peak_strategy}_{timestamp}.R"
        with open(individual_r_file, 'w') as f:
            f.write(individual_r_script)
        logger.info(f"Saved individual correlations R script to: {individual_r_file}")
        
        # Generate main analysis R script
        main_r_script = self.generate_main_r_script(peak_strategy, timestamp)
        main_r_file = r_scripts_dir / f"main_analysis_{peak_strategy}_{timestamp}.R"
        with open(main_r_file, 'w') as f:
            f.write(main_r_script)
        logger.info(f"Saved main analysis R script to: {main_r_file}")

    def generate_main_r_script(self, peak_strategy: str, timestamp: str) -> str:
        """
        Generate main R script that orchestrates the entire analysis.
        
        Args:
            peak_strategy: Peak selection strategy used
            timestamp: Timestamp for file naming
            
        Returns:
            Main R script content as string
        """
        r_script = f"""# Main RNA-ATAC Correlation Analysis R Script
# Generated on: {timestamp}
# Peak Strategy: {peak_strategy}
# Correlation Method: {self.correlation_method.capitalize()}

# This script provides a comprehensive analysis workflow
# Run sections as needed or run the entire script

# ============================================================================
# SETUP AND LIBRARIES
# ============================================================================

# Load required libraries
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(corrplot)
library(viridis)

# Set working directory (modify as needed)
# setwd("path/to/your/output/directory")

# Create output directories if they don't exist
dir.create("R_output", showWarnings = FALSE)
dir.create("R_output/plots", showWarnings = FALSE)
dir.create("R_output/tables", showWarnings = FALSE)

# ============================================================================
# DATA LOADING
# ============================================================================

# Load correlation data
correlation_data <- read.csv("correlation_results/heatmap_data_for_R_{peak_strategy}_{timestamp}.csv")

# Display data structure
cat("Data loaded successfully!\\n")
cat("Number of cluster pairs:", nrow(correlation_data), "\\n")
cat("Columns:", paste(colnames(correlation_data), collapse = ", "), "\\n\\n")

# ============================================================================
# CORRELATION HEATMAP
# ============================================================================

create_correlation_heatmap <- function() {{
  # Create correlation matrix
  corr_matrix <- correlation_data %>%
    select(RNA_cluster, ATAC_cluster, correlation) %>%
    pivot_wider(names_from = ATAC_cluster, values_from = correlation) %>%
    column_to_rownames("RNA_cluster") %>%
    as.matrix()
  
  # Create p-value matrix
  pval_matrix <- correlation_data %>%
    select(RNA_cluster, ATAC_cluster, p_value) %>%
    pivot_wider(names_from = ATAC_cluster, values_from = p_value) %>%
    column_to_rownames("RNA_cluster") %>%
    as.matrix()
  
  # Create significance matrix
  sig_matrix <- correlation_data %>%
    select(RNA_cluster, ATAC_cluster, significance) %>%
    pivot_wider(names_from = ATAC_cluster, values_from = significance) %>%
    column_to_rownames("RNA_cluster") %>%
    as.matrix()
  
  # Define color palette
  col_fun <- colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B"))
  
  # Create heatmap
  ht <- Heatmap(
    corr_matrix,
    name = "{self.correlation_method.capitalize()} Correlation",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_side = "left",
    column_names_side = "bottom",
    row_names_gp = gpar(fontsize = 12, fontface = "bold"),
    column_names_gp = gpar(fontsize = 12, fontface = "bold"),
    row_title = "RNA-seq Clusters",
    column_title = "ATAC-seq Clusters",
    row_title_gp = gpar(fontsize = 14, fontface = "bold"),
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    heatmap_legend_param = list(
      title_position = "topcenter",
      legend_direction = "horizontal",
      legend_width = unit(8, "cm")
    ),
    cell_fun = function(j, i, x, y, width, height, fill) {{
      # Add correlation values
      if (!is.na(corr_matrix[i, j])) {{
        grid.text(
          sprintf("%.3f", corr_matrix[i, j]),
          x, y,
          gp = gpar(fontsize = 10, fontface = "bold")
        )
      }}
      
      # Add significance symbols
      if (sig_matrix[i, j] != "") {{
        grid.text(
          sig_matrix[i, j],
          x, y + height * 0.3,
          gp = gpar(fontsize = 14, fontface = "bold", col = "black")
        )
      }}
    }}
  )
  
  # Save heatmap
  pdf("R_output/plots/correlation_heatmap_{peak_strategy}_{timestamp}.pdf", width = 10, height = 8)
  draw(ht, heatmap_legend_side = "bottom")
  dev.off()
  
  cat("Correlation heatmap saved to R_output/plots/\\n")
  
  return(list(corr_matrix = corr_matrix, pval_matrix = pval_matrix, sig_matrix = sig_matrix))
}}

# ============================================================================
# STATISTICAL ANALYSIS
# ============================================================================

perform_statistical_analysis <- function(corr_matrix, pval_matrix) {{
  # Convert to numeric and remove NA values
  corr_values <- as.numeric(corr_matrix)
  pval_values <- as.numeric(pval_matrix)
  
  valid_correlations <- corr_values[!is.na(corr_values)]
  valid_pvals <- pval_values[!is.na(pval_values)]
  
  # Calculate statistics
  stats_summary <- data.frame(
    Metric = c("Mean", "Median", "Min", "Max", "SD"),
    Value = c(
      mean(valid_correlations),
      median(valid_correlations),
      min(valid_correlations),
      max(valid_correlations),
      sd(valid_correlations)
    )
  )
  
  # Significance summary
  sig_summary <- data.frame(
    Threshold = c("p < 0.001", "p < 0.01", "p < 0.05"),
    Count = c(
      sum(valid_pvals < 0.001),
      sum(valid_pvals < 0.01),
      sum(valid_pvals < 0.05)
    )
  )
  
  # Save statistics
  write.csv(stats_summary, "R_output/tables/correlation_statistics.csv", row.names = FALSE)
  write.csv(sig_summary, "R_output/tables/significance_summary.csv", row.names = FALSE)
  
  # Print summary
  cat("\\n=== Statistical Summary ===\\n")
  print(stats_summary)
  cat("\\n=== Significance Summary ===\\n")
  print(sig_summary)
  
  return(list(stats = stats_summary, significance = sig_summary))
}}

# ============================================================================
# INDIVIDUAL CLUSTER PAIR ANALYSIS
# ============================================================================

analyze_individual_pairs <- function() {{
  # Get list of cluster pairs
  cluster_pairs <- unique(correlation_data$RNA_cluster)
  
  # Create individual plots for each pair
  for (rna_cluster in cluster_pairs) {{
    cluster_data <- correlation_data[correlation_data$RNA_cluster == rna_cluster, ]
    
    # Create bar plot for this RNA cluster
    p <- ggplot(cluster_data, aes(x = ATAC_cluster, y = correlation, fill = correlation)) +
      geom_bar(stat = "identity", color = "black", size = 0.5) +
      geom_text(aes(label = sprintf("%.3f", correlation)), 
                vjust = -0.5, size = 4, fontface = "bold") +
      geom_text(aes(label = significance), 
                vjust = 1.5, size = 6, fontface = "bold", color = "black") +
      scale_fill_gradient2(
        low = "#2166AC", 
        mid = "white", 
        high = "#B2182B",
        midpoint = 0,
        limits = c(-1, 1)
      ) +
      labs(
        title = paste("Correlations for", rna_cluster),
        subtitle = "vs ATAC-seq Clusters",
        x = "ATAC-seq Clusters",
        y = "Correlation Coefficient"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12),
        axis.title = element_text(size = 12, face = "bold"),
        legend.position = "none"
      )
    
    # Save plot
    ggsave(paste0("R_output/plots/", rna_cluster, "_correlations.pdf"), 
           p, width = 8, height = 6, dpi = 300)
  }}
  
  cat("Individual cluster pair plots saved to R_output/plots/\\n")
}}

# ============================================================================
# MAIN EXECUTION
# ============================================================================

# Run the complete analysis
cat("Starting RNA-ATAC correlation analysis...\\n\\n")

# Create heatmap
cat("1. Creating correlation heatmap...\\n")
heatmap_results <- create_correlation_heatmap()

# Perform statistical analysis
cat("\\n2. Performing statistical analysis...\\n")
stats_results <- perform_statistical_analysis(
  heatmap_results$corr_matrix, 
  heatmap_results$pval_matrix
)

# Analyze individual pairs
cat("\\n3. Analyzing individual cluster pairs...\\n")
analyze_individual_pairs()

cat("\\n=== Analysis Complete ===\\n")
cat("All outputs saved to R_output/ directory\\n")
cat("Check R_output/plots/ for visualizations\\n")
cat("Check R_output/tables/ for statistical results\\n")
"""
        
        return r_script

    def run_analysis(self, peak_strategy: str = 'closest_to_tss'):
        """
        Run complete correlation analysis.
        
        Args:
            peak_strategy: Strategy for selecting peaks
        """
        logger.info("Starting RNA-ATAC correlation analysis V3...")
        
        # Load all data
        logger.info("Loading data...")
        self.rna_deseq2_data = self.load_rna_deseq2_data()
        self.atac_deseq2_data = self.load_atac_deseq2_data()
        self.rna_cluster_genes = self.load_cluster_genes()
        self.atac_cluster_peaks = self.load_cluster_peaks()
        self.peaks_annotation = self.load_peaks_annotation()
        self.overlapping_genes = self.load_overlapping_genes()
        
        # Extract expression data
        logger.info("Extracting expression data...")
        self.extract_expression_data()
        
        # Process peaks to genes
        logger.info("Processing peaks to genes...")
        self.process_peaks_to_genes(peak_strategy)
        
        # Calculate correlations
        logger.info("Calculating correlations...")
        correlation_matrix, pvalue_matrix = self.calculate_correlations(peak_strategy)
        
        # Save individual correlations
        logger.info("Saving individual correlation data...")
        self.save_individual_correlations(peak_strategy)
        
        # Create individual scatter plots
        logger.info("Creating individual scatter plots...")
        self.create_individual_scatter_plots(peak_strategy)
        
        # Create heatmap
        logger.info("Creating heatmap...")
        fig = self.create_correlation_heatmap()
        
        # Save heatmap
        if fig:
            timestamp = pd.Timestamp.now().strftime("%Y%m%d_%H%M%S")
            heatmap_file = self.output_dir / f"correlation_heatmap_{peak_strategy}_{timestamp}.pdf"
            fig.savefig(heatmap_file, dpi=300, bbox_inches='tight')
            plt.close(fig)
            logger.info(f"Saved heatmap to: {heatmap_file}")
        
        # Save results
        logger.info("Saving results...")
        self.save_results(peak_strategy)
        
        # Generate R scripts
        logger.info("Generating R scripts for customization...")
        self.save_r_scripts(peak_strategy)
        
        logger.info("Analysis completed successfully!")


def main():
    """Main function to run the RNA-ATAC correlation analysis V3."""
    parser = argparse.ArgumentParser(
        description="Analyze RNA-seq and ATAC-seq cluster genes overlapping gene expression correlations V3"
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
        "--rna-cluster-results-dir",
        type=str,
        default="output/HS_DEG_hierarchical_clustering_4/cluster_results",
        help="Directory containing RNA-seq cluster results"
    )
    
    parser.add_argument(
        "--atac-cluster-results-dir",
        type=str,
        default="output/HS_DAR_hierarchical_clustering_4/cluster_results",
        help="Directory containing ATAC-seq cluster results"
    )
    
    parser.add_argument(
        "--peaks-annotation-dir",
        type=str,
        default="output/HS_DAR_annotation",
        help="Directory containing peaks annotation files"
    )
    
    parser.add_argument(
        "--overlapping-genes-dir",
        type=str,
        default="output/All_cluster_overlap_analysis/overlapping_genes",
        help="Directory containing overlapping genes files"
    )
    
    parser.add_argument(
        "--output-dir",
        type=str,
        default="output/RNA_ATAC_correlation_analysis_v3",
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
        "--fdr-method",
        type=str,
        choices=['fdr_bh', 'fdr_by', 'bonferroni', 'holm'],
        default='fdr_bh',
        help="FDR correction method to use"
    )
    
    args = parser.parse_args()
    
    # Create analyzer and run analysis
    analyzer = RNAATACCorrelationAnalyzerV3(
        rna_deseq2_dir=args.rna_deseq2_dir,
        atac_deseq2_dir=args.atac_deseq2_dir,
        rna_cluster_results_dir=args.rna_cluster_results_dir,
        atac_cluster_results_dir=args.atac_cluster_results_dir,
        peaks_annotation_dir=args.peaks_annotation_dir,
        overlapping_genes_dir=args.overlapping_genes_dir,
        output_dir=args.output_dir,
        correlation_method=args.correlation_method,
        fdr_method=args.fdr_method
    )
    
    analyzer.run_analysis(peak_strategy=args.peak_strategy)


if __name__ == "__main__":
    main()
