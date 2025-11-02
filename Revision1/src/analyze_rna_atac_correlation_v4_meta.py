#!/usr/bin/env python3
"""
RNA-seq and ATAC-seq Cluster Genes Overlapping Gene Expression Correlation Analysis V4
with Meta-Analysis and Heterogeneity Testing

This script analyzes the correlation between RNA-seq and ATAC-seq log2FoldChange expression
for overlapping genes between different clusters using per-genotype analysis and meta-analysis.

Date: 2025 (HS_ABA project)
"""

import os
import sys
import logging
import argparse
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Union, NamedTuple
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import pearsonr, spearmanr
from statsmodels.stats.multitest import multipletests
import warnings
from dataclasses import dataclass

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('rna_atac_correlation_analysis_v4_meta.log'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

# Suppress warnings
warnings.filterwarnings('ignore')


@dataclass
class GenotypeCorrelation:
    """Data class for storing per-genotype correlation results."""
    genotype: str
    correlation: float
    p_value: float
    n_genes: int
    rna_values: List[float]
    atac_values: List[float]
    valid_genes: List[str]


@dataclass
class MetaAnalysisResult:
    """Data class for storing meta-analysis results."""
    overall_correlation: float
    overall_p_value: float
    heterogeneity_q: float
    heterogeneity_p_value: float
    i2_statistic: float
    tau2_statistic: float
    direction_change_significant: bool
    genotype_results: List[GenotypeCorrelation]


class RNAATACCorrelationAnalyzerV4Meta:
    """
    Enhanced analyzer for RNA-seq and ATAC-seq cluster genes overlapping gene expression correlations
    with meta-analysis and heterogeneity testing.
    
    This class provides functionality to:
    - Load RNA-seq and ATAC-seq cluster results (genes and peaks)
    - Extract log2FoldChange values from DESeq2 files
    - Map genes to peaks using annotation data
    - Calculate per-genotype correlations between RNA-seq and ATAC-seq data
    - Perform meta-analysis to combine genotype-specific correlations
    - Test for heterogeneity between genotypes
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
        (self.output_dir / "meta_analysis").mkdir(exist_ok=True)
        (self.output_dir / "R_scripts").mkdir(exist_ok=True)
        (self.output_dir / "scatter_plots").mkdir(exist_ok=True)
        (self.output_dir / "statistical_analysis").mkdir(exist_ok=True)
        
        # Data storage
        self.rna_deseq2_data = {}
        self.atac_deseq2_data = {}
        self.rna_cluster_genes = {}
        self.atac_cluster_peaks = {}
        self.peaks_annotation = {}
        self.overlapping_genes = {}
        
        # Analysis results
        self.per_genotype_correlations = {}  # Store per-genotype correlations
        self.meta_analysis_results = {}      # Store meta-analysis results
        self.correlation_matrix = None       # Overall correlation matrix
        self.pvalue_matrix = None           # Overall p-value matrix
        self.fdr_matrix = None              # FDR-corrected p-values
        self.heterogeneity_matrix = None    # Heterogeneity test results
        self.cluster_pairs = []
        
        # Define expected clusters (4x4 matrix)
        self.expected_rna_clusters = ['cluster_1', 'cluster_2', 'cluster_3', 'cluster_4']
        self.expected_atac_clusters = ['cluster_1', 'cluster_2', 'cluster_3', 'cluster_4']
        
        # Define genotype groups for meta-analysis
        self.genotype_groups = {
            'Tak1': ['Tak1_heat_Tak1_NHS'],
            'hsfa': ['hsfa_heat_hsfa_NHS'],
            'hsfb': ['hsfb_heat_hsfb_NHS'],  # Fixed: was 'hsfb_heat_hsfa_NHS'
            'dko': ['dko_heat_dko_NHS']
        }
        
        logger.info(f"RNA-ATAC Correlation Analyzer V4 Meta initialized with {self.correlation_method} correlation")
        logger.info(f"FDR correction method: {self.fdr_method}")
        logger.info(f"Genotype groups: {list(self.genotype_groups.keys())}")
    
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
    
    def load_cluster_genes(self) -> Dict[str, List[str]]:
        """Load RNA-seq cluster genes from text files."""
        logger.info("Loading RNA-seq cluster genes...")
        
        cluster_genes = {}
        for file_path in self.rna_cluster_results_dir.glob("cluster_*_genes.txt"):
            cluster_name = file_path.stem.replace('_genes', '')
            try:
                with open(file_path, 'r') as f:
                    lines = f.readlines()
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
        """Load ATAC-seq cluster peaks from text files."""
        logger.info("Loading ATAC-seq cluster peaks...")
        
        cluster_peaks = {}
        for file_path in self.atac_cluster_results_dir.glob("cluster_*_peaks.txt"):
            cluster_name = file_path.stem.replace('_peaks', '')
            try:
                with open(file_path, 'r') as f:
                    lines = f.readlines()
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
    
    def load_overlapping_genes(self) -> Dict[str, List[str]]:
        """Load overlapping genes from text files with improved parsing."""
        logger.info("Loading overlapping genes...")
        
        overlapping_genes = {}
        for file_path in self.overlapping_genes_dir.glob("*.txt"):
            if "summary" in file_path.stem.lower():
                continue
            
            cluster_pair = file_path.stem
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
                        len(line) > 1):
                        genes.append(line)
                
                overlapping_genes[cluster_pair] = genes
                logger.info(f"Loaded {len(genes)} overlapping genes for {cluster_pair} from {file_path.name}")
                
            except Exception as e:
                logger.error(f"Error loading overlapping genes from {file_path}: {e}")
        
        return overlapping_genes
    
    def get_sample_mapping(self) -> Tuple[Dict[str, str], Dict[str, str]]:
        """Create mapping between RNA-seq and ATAC-seq sample names."""
        sample_mapping = {
            'dko_heat_dko_NHS': 'dko_HS_dko_CK',
            'hsfa_heat_hsfa_NHS': 'hsfa_HS_hsfa_CK',
            'hsfb_heat_hsfb_NHS': 'hsfb_HS_hsfb_CK',
            'Tak1_heat_Tak1_NHS': 'Tak1_HS_Tak1_CK'
        }
        
        reverse_mapping = {v: k for k, v in sample_mapping.items()}
        return sample_mapping, reverse_mapping
    
    def get_ordered_sample_names(self) -> Tuple[List[str], List[str]]:
        """Get ordered sample names to ensure consistent ordering between RNA and ATAC data."""
        ordered_rna_samples = [
            'Tak1_heat_Tak1_NHS',
            'hsfa_heat_hsfa_NHS', 
            'hsfb_heat_hsfb_NHS',
            'dko_heat_dko_NHS'
        ]
        
        sample_mapping, _ = self.get_sample_mapping()
        ordered_atac_samples = [sample_mapping[rna_sample] for rna_sample in ordered_rna_samples]
        
        return ordered_rna_samples, ordered_atac_samples
    
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
    
    def calculate_per_genotype_correlations(self, cluster_pair: str, 
                                         overlapping_genes: List[str],
                                         peak_strategy: str) -> Dict[str, GenotypeCorrelation]:
        """
        Calculate correlations for each genotype separately.
        
        Args:
            cluster_pair: Name of the cluster pair
            overlapping_genes: List of overlapping genes
            peak_strategy: Peak selection strategy
            
        Returns:
            Dictionary mapping genotype names to correlation results
        """
        logger.info(f"Calculating per-genotype correlations for {cluster_pair}")
        
        # Extract cluster names
        cluster_parts = cluster_pair.replace('_overlapping_genes', '').split('_vs_')
        if len(cluster_parts) < 2:
            logger.error(f"Invalid cluster pair format: {cluster_pair}")
            return {}
        
        rna_cluster = cluster_parts[0]
        atac_cluster = cluster_parts[1]
        
        genotype_correlations = {}
        
        # Process each genotype
        for genotype_name, sample_names in self.genotype_groups.items():
            if len(sample_names) != 1:
                logger.warning(f"Expected single sample per genotype, got {len(sample_names)} for {genotype_name}")
                continue
            
            rna_sample = sample_names[0]
            atac_sample = self.get_sample_mapping()[0].get(rna_sample)
            
            logger.info(f"Processing {genotype_name}: RNA={rna_sample}, ATAC={atac_sample}")
            
            if not atac_sample:
                logger.warning(f"No ATAC sample mapping found for {rna_sample}")
                continue
            
            # Get RNA-seq data for this genotype
            if rna_sample not in self.rna_deseq2_data:
                logger.warning(f"No RNA-seq data for {rna_sample}")
                continue
            
            # Get ATAC-seq data for this genotype
            if atac_sample not in self.atac_deseq2_data:
                logger.warning(f"No ATAC-seq data for {atac_sample}")
                continue
            
            # Extract expression data for overlapping genes
            rna_values = []
            atac_values = []
            valid_genes = []
            
            for gene in overlapping_genes:
                # Get RNA-seq expression
                if gene in self.rna_deseq2_data[rna_sample].index:
                    rna_val = self.rna_deseq2_data[rna_sample].loc[gene, 'log2FoldChange']
                else:
                    continue
                
                # Get ATAC-seq accessibility for this gene
                if atac_cluster in self.peaks_annotation:
                    annotation_df = self.peaks_annotation[atac_cluster]
                    gene_peaks = annotation_df[annotation_df['geneId'] == gene]
                    
                    if not gene_peaks.empty:
                        peaks_info = []
                        for _, row in gene_peaks.iterrows():
                            peaks_info.append({
                                'peak_id': row['Peak'],
                                'distance': abs(row['distanceToTSS']),
                                'peak_type': row['peak_type']
                            })
                        
                        selected_peak = self.select_peak_by_strategy(peaks_info, peak_strategy)
                        if selected_peak and selected_peak['peak_id'] in self.atac_deseq2_data[atac_sample].index:
                            atac_val = self.atac_deseq2_data[atac_sample].loc[selected_peak['peak_id'], 'log2FoldChange']
                            
                            # Store valid data
                            rna_values.append(rna_val)
                            atac_values.append(atac_val)
                            valid_genes.append(gene)
            
            # Calculate correlation for this genotype
            if len(rna_values) >= 3:
                try:
                    if self.correlation_method == 'pearson':
                        correlation, p_value = pearsonr(rna_values, atac_values)
                    else:
                        correlation, p_value = spearmanr(rna_values, atac_values)
                    
                    genotype_correlations[genotype_name] = GenotypeCorrelation(
                        genotype=genotype_name,
                        correlation=correlation,
                        p_value=p_value,
                        n_genes=len(valid_genes),
                        rna_values=rna_values,
                        atac_values=atac_values,
                        valid_genes=valid_genes
                    )
                    
                    logger.info(f"{cluster_pair} - {genotype_name}: r={correlation:.3f}, p={p_value:.3e}, n={len(valid_genes)}")
                    
                except Exception as e:
                    logger.error(f"Error calculating correlation for {cluster_pair} - {genotype_name}: {e}")
            else:
                logger.warning(f"Insufficient data for {cluster_pair} - {genotype_name}: {len(rna_values)} valid genes (need >=3)")
        
        logger.info(f"Completed per-genotype correlations for {cluster_pair}: {list(genotype_correlations.keys())}")
        return genotype_correlations
    
    def perform_meta_analysis(self, genotype_correlations: Dict[str, GenotypeCorrelation]) -> MetaAnalysisResult:
        """
        Perform meta-analysis to combine genotype-specific correlations.
        
        Args:
            genotype_correlations: Dictionary of genotype-specific correlation results
            
        Returns:
            MetaAnalysisResult object
        """
        if len(genotype_correlations) < 2:
            logger.warning("Need at least 2 genotypes for meta-analysis")
            return None
        
        # Extract correlation coefficients and sample sizes
        correlations = []
        sample_sizes = []
        genotype_names = []
        
        for genotype_name, result in genotype_correlations.items():
            correlations.append(result.correlation)
            sample_sizes.append(result.n_genes)
            genotype_names.append(genotype_name)
        
        # Convert correlations to Fisher's z-transformed values
        z_correlations = np.arctanh(correlations)
        
        # Calculate weights based on sample sizes
        weights = np.array(sample_sizes) - 3  # -3 for correlation degrees of freedom
        
        # Weighted mean of z-transformed correlations
        weighted_z_mean = np.average(z_correlations, weights=weights)
        
        # Convert back to correlation
        overall_correlation = np.tanh(weighted_z_mean)
        
        # Calculate standard error
        se_z = np.sqrt(1.0 / np.sum(weights))
        
        # Calculate z-statistic and p-value
        z_stat = weighted_z_mean / se_z
        overall_p_value = 2 * (1 - stats.norm.cdf(abs(z_stat)))
        
        # Test for heterogeneity (Q-statistic)
        q_statistic = np.sum(weights * (z_correlations - weighted_z_mean) ** 2)
        df = len(correlations) - 1
        heterogeneity_p_value = 1 - stats.chi2.cdf(q_statistic, df)
        
        # Calculate I² statistic
        if q_statistic > df:
            i2_statistic = (q_statistic - df) / q_statistic * 100
        else:
            i2_statistic = 0.0
        
        # Calculate tau² statistic
        if q_statistic > df:
            tau2_statistic = (q_statistic - df) / (np.sum(weights) - np.sum(weights**2) / np.sum(weights))
        else:
            tau2_statistic = 0.0
        
        # Check for significant direction changes
        positive_correlations = [c for c in correlations if c > 0]
        negative_correlations = [c for c in correlations if c < 0]
        
        direction_change_significant = False
        if len(positive_correlations) > 0 and len(negative_correlations) > 0:
            # Test if the difference between positive and negative correlations is significant
            try:
                t_stat, p_val = stats.ttest_ind(positive_correlations, negative_correlations)
                direction_change_significant = p_val < 0.05
            except:
                direction_change_significant = False
        
        return MetaAnalysisResult(
            overall_correlation=overall_correlation,
            overall_p_value=overall_p_value,
            heterogeneity_q=q_statistic,
            heterogeneity_p_value=heterogeneity_p_value,
            i2_statistic=i2_statistic,
            tau2_statistic=tau2_statistic,
            direction_change_significant=direction_change_significant,
            genotype_results=list(genotype_correlations.values())
        )
    
    def calculate_correlations_with_meta_analysis(self, peak_strategy: str = 'closest_to_tss'):
        """
        Calculate correlations with meta-analysis for each cluster pair.
        
        Args:
            peak_strategy: Peak selection strategy used
        """
        logger.info(f"Calculating correlations with meta-analysis using peak strategy: {peak_strategy}")
        
        # Initialize matrices
        correlation_matrix = pd.DataFrame(index=self.expected_rna_clusters, 
                                        columns=self.expected_atac_clusters)
        pvalue_matrix = pd.DataFrame(index=self.expected_rna_clusters, 
                                   columns=self.expected_atac_clusters)
        heterogeneity_matrix = pd.DataFrame(index=self.expected_rna_clusters, 
                                          columns=self.expected_atac_clusters)
        
        # Store all p-values for FDR correction
        all_pvalues = []
        pvalue_positions = []
        
        # Process each expected cluster pair
        for rna_cluster in self.expected_rna_clusters:
            for atac_cluster in self.expected_atac_clusters:
                cluster_pair = f"{rna_cluster}_vs_{atac_cluster}_overlapping_genes"
                
                # Check if we have overlapping genes for this pair
                if cluster_pair not in self.overlapping_genes:
                    logger.warning(f"No overlapping genes file found for {cluster_pair}")
                    correlation_matrix.loc[rna_cluster, atac_cluster] = np.nan
                    pvalue_matrix.loc[rna_cluster, atac_cluster] = np.nan
                    heterogeneity_matrix.loc[rna_cluster, atac_cluster] = np.nan
                    continue
                
                overlapping_genes = self.overlapping_genes[cluster_pair]
                
                if not overlapping_genes:
                    logger.warning(f"No overlapping genes for {cluster_pair}")
                    correlation_matrix.loc[rna_cluster, atac_cluster] = np.nan
                    pvalue_matrix.loc[rna_cluster, atac_cluster] = np.nan
                    heterogeneity_matrix.loc[rna_cluster, atac_cluster] = np.nan
                    continue
                
                logger.info(f"Processing {cluster_pair}: {len(overlapping_genes)} overlapping genes")
                
                # Calculate per-genotype correlations
                genotype_correlations = self.calculate_per_genotype_correlations(
                    cluster_pair, overlapping_genes, peak_strategy
                )
                
                # Store per-genotype results
                self.per_genotype_correlations[cluster_pair] = genotype_correlations
                
                # Perform meta-analysis
                if len(genotype_correlations) >= 2:
                    logger.info(f"Performing meta-analysis for {cluster_pair} with {len(genotype_correlations)} genotypes: {list(genotype_correlations.keys())}")
                    meta_result = self.perform_meta_analysis(genotype_correlations)
                    
                    if meta_result:
                        # Store meta-analysis results
                        self.meta_analysis_results[cluster_pair] = meta_result
                        
                        # Update matrices
                        correlation_matrix.loc[rna_cluster, atac_cluster] = meta_result.overall_correlation
                        pvalue_matrix.loc[rna_cluster, atac_cluster] = meta_result.overall_p_value
                        heterogeneity_matrix.loc[rna_cluster, atac_cluster] = meta_result.heterogeneity_p_value
                        
                        # Store p-value for FDR correction
                        all_pvalues.append(meta_result.overall_p_value)
                        pvalue_positions.append((rna_cluster, atac_cluster))
                        
                        logger.info(f"{cluster_pair} meta-analysis: r={meta_result.overall_correlation:.3f}, "
                                  f"p={meta_result.overall_p_value:.3e}, "
                                  f"heterogeneity p={meta_result.heterogeneity_p_value:.3e}")
                    else:
                        logger.warning(f"Meta-analysis failed for {cluster_pair}")
                        correlation_matrix.loc[rna_cluster, atac_cluster] = np.nan
                        pvalue_matrix.loc[rna_cluster, atac_cluster] = np.nan
                        heterogeneity_matrix.loc[rna_cluster, atac_cluster] = np.nan
                else:
                    logger.warning(f"Insufficient genotypes for meta-analysis in {cluster_pair}: {len(genotype_correlations)} genotypes (need >=2)")
                    correlation_matrix.loc[rna_cluster, atac_cluster] = np.nan
                    pvalue_matrix.loc[rna_cluster, atac_cluster] = np.nan
                    heterogeneity_matrix.loc[rna_cluster, atac_cluster] = np.nan
        
        # Perform FDR correction
        if all_pvalues:
            logger.info(f"Performing FDR correction using {self.fdr_method} method...")
            try:
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
                self.fdr_matrix = pvalue_matrix.copy()
        
        self.correlation_matrix = correlation_matrix
        self.pvalue_matrix = pvalue_matrix
        self.heterogeneity_matrix = heterogeneity_matrix
        
        return correlation_matrix, pvalue_matrix, heterogeneity_matrix

    def create_per_genotype_scatter_plots(self, peak_strategy: str):
        """Create scatter plots showing per-genotype correlations for each cluster pair."""
        logger.info("Creating per-genotype scatter plots...")
        
        scatter_dir = self.output_dir / "scatter_plots"
        (scatter_dir / peak_strategy).mkdir(exist_ok=True)
        
        for cluster_pair, genotype_correlations in self.per_genotype_correlations.items():
            if not genotype_correlations:
                continue
            
            # Create figure with subplots for each genotype
            n_genotypes = len(genotype_correlations)
            fig, axes = plt.subplots(2, 2, figsize=(15, 12))
            axes = axes.flatten()
            
            # Plot each genotype
            for i, (genotype_name, result) in enumerate(genotype_correlations.items()):
                if i >= len(axes):
                    break
                
                ax = axes[i]
                
                if result.rna_values and result.atac_values:
                    # Create scatter plot
                    scatter = ax.scatter(result.rna_values, result.atac_values, 
                                       alpha=0.7, s=50, edgecolors='black', linewidth=0.5)
                    
                    # Add trend line
                    if len(result.rna_values) > 1:
                        z = np.polyfit(result.rna_values, result.atac_values, 1)
                        p = np.poly1d(z)
                        ax.plot(result.rna_values, p(result.rna_values), "r--", alpha=0.8, linewidth=2)
                    
                    # Add correlation info
                    ax.text(0.05, 0.95, 
                           f'r = {result.correlation:.3f}\np = {result.p_value:.3e}\nn = {result.n_genes}', 
                           transform=ax.transAxes, 
                           bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
                           fontsize=10, verticalalignment='top')
                    
                    # Customize plot
                    ax.set_xlabel('RNA-seq log2FoldChange', fontsize=10)
                    ax.set_ylabel('ATAC-seq log2FoldChange', fontsize=10)
                    ax.set_title(f'{genotype_name.upper()}', fontsize=12, fontweight='bold')
                    ax.grid(True, alpha=0.3)
                    ax.axhline(y=0, color='black', linestyle='-', alpha=0.3)
                    ax.axvline(x=0, color='black', linestyle='-', alpha=0.3)
                else:
                    ax.text(0.5, 0.5, 'No valid data', 
                           transform=ax.transAxes, ha='center', va='center',
                           fontsize=12, style='italic')
                    ax.set_title(f'{genotype_name.upper()}', fontsize=12, fontweight='bold')
            
            # Remove empty subplots
            for i in range(n_genotypes, len(axes)):
                fig.delaxes(axes[i])
            
            # Add overall title
            fig.suptitle(f'Per-Genotype Correlations: {cluster_pair}\nPeak Strategy: {peak_strategy}', 
                         fontsize=16, fontweight='bold')
            
            plt.tight_layout()
            
            # Save plot
            output_file = scatter_dir / peak_strategy / f"{cluster_pair}_per_genotype_scatter.pdf"
            fig.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close(fig)
            
            logger.info(f"Saved per-genotype scatter plots for {cluster_pair} to: {output_file}")
    
    def create_meta_analysis_forest_plot(self, peak_strategy: str):
        """Create forest plot showing meta-analysis results for each cluster pair."""
        logger.info("Creating meta-analysis forest plots...")
        
        forest_dir = self.output_dir / "scatter_plots"
        (forest_dir / peak_strategy).mkdir(exist_ok=True)
        
        for cluster_pair, meta_result in self.meta_analysis_results.items():
            if not meta_result.genotype_results:
                continue
            
            # Create forest plot
            fig, ax = plt.subplots(figsize=(10, 8))
            
            # Extract data for forest plot
            genotypes = [result.genotype for result in meta_result.genotype_results]
            correlations = [result.correlation for result in meta_result.genotype_results]
            p_values = [result.p_value for result in meta_result.genotype_results]
            n_genes = [result.n_genes for result in meta_result.genotype_results]
            
            # Add overall result
            genotypes.append('Overall (Meta)')
            correlations.append(meta_result.overall_correlation)
            p_values.append(meta_result.overall_p_value)
            n_genes.append(sum(n_genes))
            
            # Create forest plot
            y_pos = np.arange(len(genotypes))
            colors = ['blue'] * (len(genotypes) - 1) + ['red']
            
            # Plot correlation coefficients
            ax.scatter(correlations, y_pos, c=colors, s=100, alpha=0.8, edgecolors='black')
            
            # Add error bars (95% CI approximation)
            for i, (corr, n) in enumerate(zip(correlations, n_genes)):
                if n > 3:
                    # Approximate 95% CI for correlation
                    se = np.sqrt((1 - corr**2) / (n - 3))
                    ci_95 = 1.96 * se
                    ax.plot([corr - ci_95, corr + ci_95], [i, i], 'k-', alpha=0.7)
            
            # Add vertical line at 0
            ax.axvline(x=0, color='black', linestyle='--', alpha=0.5)
            
            # Customize plot
            ax.set_yticks(y_pos)
            ax.set_yticklabels(genotypes)
            ax.set_xlabel('Correlation Coefficient', fontsize=12, fontweight='bold')
            ax.set_title(f'Forest Plot: {cluster_pair}\nPeak Strategy: {peak_strategy}', 
                         fontsize=14, fontweight='bold')
            ax.grid(True, alpha=0.3)
            
            # Add heterogeneity information
            if meta_result.heterogeneity_p_value < 0.05:
                ax.text(0.02, 0.98, f'Heterogeneity: p = {meta_result.heterogeneity_p_value:.3e}\nI² = {meta_result.i2_statistic:.1f}%', 
                       transform=ax.transAxes, 
                       bbox=dict(boxstyle="round,pad=0.5", facecolor="yellow", alpha=0.8),
                       fontsize=10, verticalalignment='top', fontweight='bold')
            
            # Add direction change information
            if meta_result.direction_change_significant:
                ax.text(0.02, 0.90, 'Direction Change: SIGNIFICANT', 
                       transform=ax.transAxes, 
                       bbox=dict(boxstyle="round,pad=0.5", facecolor="red", alpha=0.8),
                       fontsize=10, verticalalignment='top', fontweight='bold', color='white')
            
            plt.tight_layout()
            
            # Save plot
            output_file = forest_dir / peak_strategy / f"{cluster_pair}_forest_plot.pdf"
            fig.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close(fig)
            
            logger.info(f"Saved forest plot for {cluster_pair} to: {output_file}")
    
    def create_correlation_heatmap(self) -> plt.Figure:
        """Create correlation heatmap with significance indicators."""
        if self.correlation_matrix is None or self.pvalue_matrix is None:
            logger.error("No correlation data available for heatmap")
            return None
        
        # Create figure with three subplots
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(24, 8))
        
        # Convert matrices to numeric
        corr_matrix_numeric = self.correlation_matrix.astype(float).fillna(0)
        pval_matrix_numeric = self.pvalue_matrix.astype(float).fillna(1)
        het_matrix_numeric = self.heterogeneity_matrix.astype(float).fillna(1) if self.heterogeneity_matrix is not None else None
        
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
        
        ax1.set_title(f'Meta-Analysis Correlation Heatmap\n({self.correlation_method.capitalize()})', 
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
        
        # 3. Heterogeneity heatmap
        if het_matrix_numeric is not None:
            sns.heatmap(het_matrix_numeric, 
                       annot=True, 
                       fmt='.3e',
                       cmap='Blues_r',
                       vmin=0, 
                       vmax=0.05,
                       cbar_kws={'label': 'Heterogeneity P-value'},
                       ax=ax3)
            
            ax3.set_title('Heterogeneity Test Heatmap', fontsize=14, fontweight='bold')
            ax3.set_xlabel('ATAC-seq Clusters', fontsize=12)
            ax3.set_ylabel('RNA-seq Clusters', fontsize=12)
        else:
            ax3.text(0.5, 0.5, 'Heterogeneity data not available', 
                     transform=ax3.transAxes, ha='center', va='center',
                     fontsize=14, style='italic')
            ax3.set_title('Heterogeneity Test Not Available', fontsize=14, fontweight='bold')
        
        plt.tight_layout()
        return fig
    
    def get_significance_symbol(self, p_value: float, fdr_value: float = None) -> str:
        """Get significance symbol based on p-value and FDR value."""
        if np.isnan(p_value):
            return ""
        
        test_value = fdr_value if fdr_value is not None and not np.isnan(fdr_value) else p_value
        
        if test_value < 0.001:
            return "***"
        elif test_value < 0.01:
            return "**"
        elif test_value < 0.05:
            return "*"
        else:
            return ""
    
    def save_detailed_statistical_analysis(self, peak_strategy: str):
        """Save detailed statistical analysis results."""
        logger.info("Saving detailed statistical analysis...")
        
        stats_dir = self.output_dir / "statistical_analysis"
        timestamp = pd.Timestamp.now().strftime("%Y%m%d_%H%M%S")
        
        # 1. Per-genotype correlation summary
        per_genotype_summary = []
        for cluster_pair, genotype_correlations in self.per_genotype_correlations.items():
            for genotype_name, result in genotype_correlations.items():
                per_genotype_summary.append({
                    'Cluster_Pair': cluster_pair,
                    'Genotype': genotype_name,
                    'Correlation': result.correlation,
                    'P_Value': result.p_value,
                    'N_Genes': result.n_genes,
                    'Significant': result.p_value < 0.05
                })
        
        if per_genotype_summary:
            per_genotype_df = pd.DataFrame(per_genotype_summary)
            per_genotype_file = stats_dir / f"per_genotype_correlations_{peak_strategy}_{timestamp}.csv"
            per_genotype_df.to_csv(per_genotype_file, index=False)
            logger.info(f"Saved per-genotype correlations to: {per_genotype_file}")
        
        # 2. Meta-analysis summary
        meta_summary = []
        for cluster_pair, meta_result in self.meta_analysis_results.items():
            meta_summary.append({
                'Cluster_Pair': cluster_pair,
                'Overall_Correlation': meta_result.overall_correlation,
                'Overall_P_Value': meta_result.overall_p_value,
                'Heterogeneity_Q': meta_result.heterogeneity_q,
                'Heterogeneity_P_Value': meta_result.heterogeneity_p_value,
                'I2_Statistic': meta_result.i2_statistic,
                'Tau2_Statistic': meta_result.tau2_statistic,
                'Direction_Change_Significant': meta_result.direction_change_significant,
                'N_Genotypes': len(meta_result.genotype_results)
            })
        
        if meta_summary:
            meta_df = pd.DataFrame(meta_summary)
            meta_file = stats_dir / f"meta_analysis_summary_{peak_strategy}_{timestamp}.csv"
            meta_df.to_csv(meta_file, index=False)
            logger.info(f"Saved meta-analysis summary to: {meta_file}")
        
        # 3. Detailed statistical report
        report_file = stats_dir / f"detailed_statistical_report_{peak_strategy}_{timestamp}.txt"
        with open(report_file, 'w') as f:
            f.write("RNA-seq vs ATAC-seq Correlation Analysis - Detailed Statistical Report\n")
            f.write("=" * 80 + "\n\n")
            f.write(f"Analysis Date: {timestamp}\n")
            f.write(f"Correlation Method: {self.correlation_method.capitalize()}\n")
            f.write(f"FDR Correction Method: {self.fdr_method}\n")
            f.write(f"Peak Selection Strategy: {peak_strategy}\n\n")
            
            # Overall statistics
            f.write("OVERALL STATISTICS\n")
            f.write("-" * 40 + "\n")
            f.write(f"Total cluster pairs analyzed: {len(self.meta_analysis_results)}\n")
            f.write(f"Total genotypes: {len(self.genotype_groups)}\n")
            f.write(f"Genotypes: {', '.join(self.genotype_groups.keys())}\n\n")
            
            # Per-genotype analysis
            f.write("PER-GENOTYPE ANALYSIS\n")
            f.write("-" * 40 + "\n")
            for cluster_pair, genotype_correlations in self.per_genotype_correlations.items():
                f.write(f"\n{cluster_pair}:\n")
                for genotype_name, result in genotype_correlations.items():
                    f.write(f"  {genotype_name}: r={result.correlation:.3f}, p={result.p_value:.3e}, n={result.n_genes}\n")
            
            # Meta-analysis results
            f.write("\n\nMETA-ANALYSIS RESULTS\n")
            f.write("-" * 40 + "\n")
            for cluster_pair, meta_result in self.meta_analysis_results.items():
                f.write(f"\n{cluster_pair}:\n")
                f.write(f"  Overall correlation: {meta_result.overall_correlation:.3f}\n")
                f.write(f"  Overall p-value: {meta_result.overall_p_value:.3e}\n")
                f.write(f"  Heterogeneity Q: {meta_result.heterogeneity_q:.3f}\n")
                f.write(f"  Heterogeneity p-value: {meta_result.heterogeneity_p_value:.3e}\n")
                f.write(f"  I² statistic: {meta_result.i2_statistic:.1f}%\n")
                f.write(f"  Tau² statistic: {meta_result.tau2_statistic:.3f}\n")
                f.write(f"  Direction change significant: {meta_result.direction_change_significant}\n")
            
            # Heterogeneity interpretation
            f.write("\n\nHETEROGENEITY INTERPRETATION\n")
            f.write("-" * 40 + "\n")
            f.write("I² values:\n")
            f.write("  0-25%: Low heterogeneity\n")
            f.write("  25-50%: Moderate heterogeneity\n")
            f.write("  50-75%: Substantial heterogeneity\n")
            f.write("  75-100%: Considerable heterogeneity\n\n")
            
            # Direction change analysis
            f.write("DIRECTION CHANGE ANALYSIS\n")
            f.write("-" * 40 + "\n")
            significant_direction_changes = [pair for pair, result in self.meta_analysis_results.items() 
                                          if result.direction_change_significant]
            if significant_direction_changes:
                f.write("Significant direction changes detected in:\n")
                for pair in significant_direction_changes:
                    f.write(f"  {pair}\n")
            else:
                f.write("No significant direction changes detected.\n")
        
        logger.info(f"Saved detailed statistical report to: {report_file}")
    
    def generate_r_scripts(self, peak_strategy: str):
        """Generate comprehensive R scripts for analysis and visualization."""
        logger.info("Generating R scripts for analysis and visualization...")
        
        r_scripts_dir = self.output_dir / "R_scripts"
        timestamp = pd.Timestamp.now().strftime("%Y%m%d_%H%M%S")
        
        # Generate main analysis R script
        main_r_script = self.generate_main_analysis_r_script(peak_strategy, timestamp)
        main_r_file = r_scripts_dir / f"main_analysis_{peak_strategy}_{timestamp}.R"
        with open(main_r_file, 'w') as f:
            f.write(main_r_script)
        logger.info(f"Saved main analysis R script to: {main_r_file}")
        
        # Generate meta-analysis R script
        meta_r_script = self.generate_meta_analysis_r_script(peak_strategy, timestamp)
        meta_r_file = r_scripts_dir / f"meta_analysis_{peak_strategy}_{timestamp}.R"
        with open(meta_r_file, 'w') as f:
            f.write(meta_r_script)
        logger.info(f"Saved meta-analysis R script to: {meta_r_file}")
        
        # Generate visualization R script
        viz_r_script = self.generate_visualization_r_script(peak_strategy, timestamp)
        viz_r_file = r_scripts_dir / f"visualization_{peak_strategy}_{timestamp}.R"
        with open(viz_r_file, 'w') as f:
            f.write(viz_r_script)
        logger.info(f"Saved visualization R script to: {viz_r_file}")
    
    def generate_main_analysis_r_script(self, peak_strategy: str, timestamp: str) -> str:
        """Generate main R script for comprehensive analysis."""
        return f"""# Main RNA-ATAC Correlation Analysis with Meta-Analysis R Script
# Generated on: {timestamp}
# Peak Strategy: {peak_strategy}
# Correlation Method: {self.correlation_method.capitalize()}

# Load required libraries
library(metafor)
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(gridExtra)
library(viridis)

# Set working directory
# setwd("path/to/your/output/directory")

# Create output directories
dir.create("R_output", showWarnings = FALSE)
dir.create("R_output/plots", showWarnings = FALSE)
dir.create("R_output/tables", showWarnings = FALSE)
dir.create("R_output/meta_analysis", showWarnings = FALSE)

# Load data
per_genotype_data <- read.csv("statistical_analysis/per_genotype_correlations_{peak_strategy}_{timestamp}.csv")
meta_analysis_data <- read.csv("statistical_analysis/meta_analysis_summary_{peak_strategy}_{timestamp}.csv")

# Display data structure
cat("Data loaded successfully!\\n")
cat("Per-genotype correlations:", nrow(per_genotype_data), "\\n")
cat("Meta-analysis results:", nrow(meta_analysis_data), "\\n\\n")

# Perform additional statistical analyses
# ... (rest of the R script will be implemented in the next part)
"""
    
    def generate_meta_analysis_r_script(self, peak_strategy: str, timestamp: str) -> str:
        """Generate R script specifically for meta-analysis."""
        return f"""# Meta-Analysis R Script for RNA-ATAC Correlation
# Generated on: {timestamp}
# Peak Strategy: {peak_strategy}

library(metafor)
library(ggplot2)
library(dplyr)

# Load data
per_genotype_data <- read.csv("statistical_analysis/per_genotype_correlations_{peak_strategy}_{timestamp}.csv")
meta_analysis_data <- read.csv("statistical_analysis/meta_analysis_summary_{peak_strategy}_{timestamp}.csv")

# Function to perform meta-analysis for a cluster pair
perform_cluster_meta_analysis <- function(cluster_pair_data) {{
  # Convert correlations to Fisher's z
  z_correlations <- atanh(cluster_pair_data$Correlation)
  
  # Calculate standard errors
  se_z <- sqrt(1 / (cluster_pair_data$N_Genes - 3))
  
  # Perform random-effects meta-analysis
  meta_result <- rma(yi = z_correlations, sei = se_z, method = "REML")
  
  # Convert back to correlation scale
  overall_correlation <- tanh(meta_result$b)
  overall_ci_lower <- tanh(meta_result$ci.lb)
  overall_ci_upper <- tanh(meta_result$ci.ub)
  
  return(list(
    correlation = overall_correlation,
    ci_lower = overall_ci_lower,
    ci_upper = overall_ci_upper,
    heterogeneity_q = meta_result$Q,
    heterogeneity_p = meta_result$Qp,
    i2 = meta_result$I2,
    tau2 = meta_result$tau2
  ))
}}

# Perform meta-analysis for each cluster pair
cluster_pairs <- unique(per_genotype_data$Cluster_Pair)
meta_results <- list()

for (pair in cluster_pairs) {{
  pair_data <- per_genotype_data[per_genotype_data$Cluster_Pair == pair, ]
  meta_results[[pair]] <- perform_cluster_meta_analysis(pair_data)
}}

# Create forest plots
for (pair in names(meta_results)) {{
  pair_data <- per_genotype_data[per_genotype_data$Cluster_Pair == pair, ]
  
  # Create forest plot
  forest_plot <- ggplot(pair_data, aes(x = Correlation, y = Genotype)) +
    geom_point(aes(size = N_Genes), alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    geom_vline(xintercept = meta_results[[pair]]$correlation, 
               linetype = "solid", color = "blue", size = 1.2) +
    geom_errorbarh(aes(xmin = Correlation - 1.96 * sqrt(1/(N_Genes-3)), 
                       xmax = Correlation + 1.96 * sqrt(1/(N_Genes-3))), 
                   height = 0.2) +
    labs(
      title = paste("Forest Plot:", pair),
      subtitle = paste("Overall r =", round(meta_results[[pair]]$correlation, 3),
                      "I² =", round(meta_results[[pair]]$i2, 1), "%"),
      x = "Correlation Coefficient",
      y = "Genotype",
      size = "Number of Genes"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.title = element_text(size = 12, face = "bold")
    )
  
  # Save plot
  ggsave(paste0("R_output/meta_analysis/", gsub(" ", "_", pair), "_forest_plot.pdf"), 
         forest_plot, width = 10, height = 6, dpi = 300)
}}

cat("Meta-analysis completed! Check R_output/meta_analysis/ for forest plots.\\n")
"""
    
    def generate_visualization_r_script(self, peak_strategy: str, timestamp: str) -> str:
        """Generate R script for advanced visualizations."""
        return f"""# Advanced Visualization R Script for RNA-ATAC Correlation
# Generated on: {timestamp}
# Peak Strategy: {peak_strategy}

library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(circlize)
library(viridis)

# Load data
per_genotype_data <- read.csv("statistical_analysis/per_genotype_correlations_{peak_strategy}_{timestamp}.csv")
meta_analysis_data <- read.csv("statistical_analysis/meta_analysis_summary_{peak_strategy}_{timestamp}.csv")

# Create correlation matrix for heatmap
correlation_matrix <- per_genotype_data %>%
  group_by(Cluster_Pair, Genotype) %>%
  summarise(Correlation = mean(Correlation), .groups = 'drop') %>%
  pivot_wider(names_from = Genotype, values_from = Correlation) %>%
  column_to_rownames("Cluster_Pair") %>%
  as.matrix()

# Create ComplexHeatmap
col_fun <- colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B"))

ht <- Heatmap(
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
      grid.text(
        sprintf("%.3f", correlation_matrix[i, j]),
        x, y,
        gp = gpar(fontsize = 8, fontface = "bold")
      )
    }}
  }}
)

# Save heatmap
pdf("R_output/plots/correlation_heatmap_{peak_strategy}_{timestamp}.pdf", width = 10, height = 8)
draw(ht)
dev.off()

# Create heterogeneity analysis plot
heterogeneity_plot <- ggplot(meta_analysis_data, aes(x = I2_Statistic, y = Heterogeneity_P_Value)) +
  geom_point(aes(size = N_Genotypes, color = Direction_Change_Significant), alpha = 0.7) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 50, linetype = "dashed", color = "orange") +
  scale_color_manual(values = c("FALSE" = "blue", "TRUE" = "red")) +
  labs(
    title = "Heterogeneity Analysis",
    subtitle = "I² vs Heterogeneity P-value",
    x = "I² Statistic (%)",
    y = "Heterogeneity P-value",
    size = "Number of Genotypes",
    color = "Direction Change Significant"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    axis.title = element_text(size = 12, face = "bold")
  )

# Save heterogeneity plot
ggsave(paste0("R_output/plots/heterogeneity_analysis_{peak_strategy}_{timestamp}.pdf"), 
       heterogeneity_plot, width = 10, height = 6, dpi = 300)

cat("Visualization completed! Check R_output/plots/ for generated plots.\\n")
"""
    
    def run_analysis(self, peak_strategy: str = 'closest_to_tss'):
        """Run complete correlation analysis with meta-analysis."""
        logger.info("Starting RNA-ATAC correlation analysis V4 with meta-analysis...")
        
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
        self.extract_expression_data(peak_strategy)
        
        # Process peaks to genes
        logger.info("Processing peaks to genes...")
        self.process_peaks_to_genes(peak_strategy)
        
        # Calculate correlations with meta-analysis
        logger.info("Calculating correlations with meta-analysis...")
        correlation_matrix, pvalue_matrix, heterogeneity_matrix = self.calculate_correlations_with_meta_analysis(peak_strategy)
        
        # Save individual correlations
        logger.info("Saving individual correlation data...")
        self.save_individual_correlations(peak_strategy)
        
        # Save meta-analysis results
        logger.info("Saving meta-analysis results...")
        self.save_meta_analysis_results(peak_strategy)
        
        # Save correlation results
        logger.info("Saving correlation results...")
        self.save_correlation_results(peak_strategy)
        
        # Create visualizations
        logger.info("Creating visualizations...")
        self.create_per_genotype_scatter_plots(peak_strategy)
        self.create_meta_analysis_forest_plot(peak_strategy)
        
        # Create heatmap
        logger.info("Creating correlation heatmap...")
        fig = self.create_correlation_heatmap()
        
        # Save heatmap
        if fig:
            timestamp = pd.Timestamp.now().strftime("%Y%m%d_%H%M%S")
            heatmap_file = self.output_dir / f"correlation_heatmap_{peak_strategy}_{timestamp}.pdf"
            fig.savefig(heatmap_file, dpi=300, bbox_inches='tight')
            plt.close(fig)
            logger.info(f"Saved heatmap to: {heatmap_file}")
        
        # Save detailed statistical analysis
        logger.info("Saving detailed statistical analysis...")
        self.save_detailed_statistical_analysis(peak_strategy)
        
        # Generate R scripts
        logger.info("Generating R scripts for customization...")
        self.generate_r_scripts(peak_strategy)
        
        logger.info("Analysis completed successfully!")

    def extract_expression_data(self, peak_strategy: str):
        """
        Extract expression data for each cluster and save to individual folders.
        
        Args:
            peak_strategy: Peak selection strategy used
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
    
    def process_peaks_to_genes(self, peak_strategy: str):
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
        for cluster_pair, data in self.per_genotype_correlations.items():
            if not data:
                continue
                
            # Create filename for this cluster pair
            cluster_pair_name = cluster_pair.replace(' vs ', '_vs_')
            data_file = individual_dir / f"{cluster_pair_name}_correlation_data.csv"
            
            # Prepare data for saving
            rows = []
            for genotype_name, result in data.items():
                # Extract cluster names from cluster_pair format: "cluster1_vs_cluster2_overlapping_genes"
                cluster_parts = cluster_pair.replace('_overlapping_genes', '').split('_vs_')
                if len(cluster_parts) >= 2:
                    rna_cluster = cluster_parts[0]
                    atac_cluster = cluster_parts[1]
                else:
                    rna_cluster = "unknown"
                    atac_cluster = "unknown"
                
                row = {
                    'Gene': genotype_name,
                    'RNA_cluster': rna_cluster,
                    'ATAC_cluster': atac_cluster,
                    'Correlation': result.correlation,
                    'P_Value': result.p_value,
                    'N_Genes': result.n_genes,
                    'Significant': result.p_value < 0.05
                }
                
                # Add RNA and ATAC values
                if result.rna_values and result.atac_values:
                    row['RNA_values'] = ','.join(map(str, result.rna_values))
                    row['ATAC_values'] = ','.join(map(str, result.atac_values))
                    row['Valid_Genes'] = ','.join(result.valid_genes)
                else:
                    row['RNA_values'] = ''
                    row['ATAC_values'] = ''
                    row['Valid_Genes'] = ''
                
                rows.append(row)
            
            # Save to CSV
            if rows:
                df = pd.DataFrame(rows)
                df.to_csv(data_file, index=False)
                logger.info(f"Saved individual correlation data for {cluster_pair} to: {data_file}")
        
        # Save summary of individual correlations
        summary_file = individual_dir / f"individual_correlations_summary_{timestamp}.csv"
        summary_rows = []
        
        for cluster_pair, data in self.per_genotype_correlations.items():
            # Extract cluster names from cluster_pair format: "cluster1_vs_cluster2_overlapping_genes"
            cluster_parts = cluster_pair.replace('_overlapping_genes', '').split('_vs_')
            if len(cluster_parts) >= 2:
                rna_cluster = cluster_parts[0]
                atac_cluster = cluster_parts[1]
            else:
                rna_cluster = "unknown"
                atac_cluster = "unknown"
            
            for genotype_name, result in data.items():
                summary_rows.append({
                    'RNA_cluster': rna_cluster,
                    'ATAC_cluster': atac_cluster,
                    'Genotype': genotype_name,
                    'correlation': result.correlation,
                    'p_value': result.p_value,
                    'n_genes': result.n_genes,
                    'significant': result.p_value < 0.05
                })
        
        if summary_rows:
            summary_df = pd.DataFrame(summary_rows)
            summary_df.to_csv(summary_file, index=False)
            logger.info(f"Saved individual correlations summary to: {summary_file}")
    
    def save_meta_analysis_results(self, peak_strategy: str):
        """
        Save meta-analysis results to files.
        
        Args:
            peak_strategy: Peak strategy used
        """
        logger.info("Saving meta-analysis results...")
        
        meta_dir = self.output_dir / "meta_analysis"
        (meta_dir / peak_strategy).mkdir(exist_ok=True)
        timestamp = pd.Timestamp.now().strftime("%Y%m%d_%H%M%S")
        
        # Save detailed meta-analysis results for each cluster pair
        for cluster_pair, meta_result in self.meta_analysis_results.items():
            if not meta_result:
                continue
            
            # Create detailed results file
            cluster_pair_name = cluster_pair.replace(' vs ', '_vs_')
            detailed_file = meta_dir / peak_strategy / f"{cluster_pair_name}_detailed_meta_analysis.csv"
            
            # Prepare detailed data
            detailed_rows = []
            
            # Add genotype-specific results
            for genotype_result in meta_result.genotype_results:
                detailed_rows.append({
                    'Cluster_Pair': cluster_pair,
                    'Analysis_Type': 'Genotype_Specific',
                    'Genotype': genotype_result.genotype,
                    'Correlation': genotype_result.correlation,
                    'P_Value': genotype_result.p_value,
                    'N_Genes': genotype_result.n_genes,
                    'Significant': genotype_result.p_value < 0.05
                })
            
            # Add overall meta-analysis result
            detailed_rows.append({
                'Cluster_Pair': cluster_pair,
                'Analysis_Type': 'Meta_Analysis_Overall',
                'Genotype': 'Overall',
                'Correlation': meta_result.overall_correlation,
                'P_Value': meta_result.overall_p_value,
                'N_Genes': sum(r.n_genes for r in meta_result.genotype_results),
                'Significant': meta_result.overall_p_value < 0.05
            })
            
            # Save detailed results
            if detailed_rows:
                detailed_df = pd.DataFrame(detailed_rows)
                detailed_df.to_csv(detailed_file, index=False)
                logger.info(f"Saved detailed meta-analysis for {cluster_pair} to: {detailed_file}")
        
        # Save comprehensive meta-analysis summary
        summary_file = meta_dir / peak_strategy / f"comprehensive_meta_analysis_summary_{timestamp}.csv"
        summary_rows = []
        
        for cluster_pair, meta_result in self.meta_analysis_results.items():
            if not meta_result:
                continue
            
            # Extract cluster names
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
                'Cluster_Pair': cluster_pair,
                'Overall_Correlation': meta_result.overall_correlation,
                'Overall_P_Value': meta_result.overall_p_value,
                'Overall_Significant': meta_result.overall_p_value < 0.05,
                'Heterogeneity_Q': meta_result.heterogeneity_q,
                'Heterogeneity_P_Value': meta_result.heterogeneity_p_value,
                'Heterogeneity_Significant': meta_result.heterogeneity_p_value < 0.05,
                'I2_Statistic': meta_result.i2_statistic,
                'Tau2_Statistic': meta_result.tau2_statistic,
                'Direction_Change_Significant': meta_result.direction_change_significant,
                'N_Genotypes': len(meta_result.genotype_results),
                'Total_Genes': sum(r.n_genes for r in meta_result.genotype_results)
            })
        
        if summary_rows:
            summary_df = pd.DataFrame(summary_rows)
            summary_df.to_csv(summary_file, index=False)
            logger.info(f"Saved comprehensive meta-analysis summary to: {summary_file}")
        
        # Save heterogeneity interpretation guide
        guide_file = meta_dir / peak_strategy / f"heterogeneity_interpretation_guide_{timestamp}.txt"
        with open(guide_file, 'w') as f:
            f.write("Meta-Analysis Heterogeneity Interpretation Guide\n")
            f.write("=" * 50 + "\n\n")
            f.write("I² Statistic Interpretation:\n")
            f.write("- 0-25%: Low heterogeneity (genotypes are similar)\n")
            f.write("- 25-50%: Moderate heterogeneity\n")
            f.write("- 50-75%: Substantial heterogeneity\n")
            f.write("- 75-100%: Considerable heterogeneity (genotypes differ significantly)\n\n")
            
            f.write("Q-Statistic Interpretation:\n")
            f.write("- Q-statistic follows chi-square distribution\n")
            f.write("- Degrees of freedom = number of genotypes - 1\n")
            f.write("- P-value < 0.05 indicates significant heterogeneity\n\n")
            
            f.write("Tau² Statistic:\n")
            f.write("- Between-genotype variance component\n")
            f.write("- Higher values indicate more heterogeneity\n")
            f.write("- 0 indicates no heterogeneity\n\n")
            
            f.write("Direction Change Analysis:\n")
            f.write("- Tests if positive/negative correlations differ significantly\n")
            f.write("- Significant changes suggest genotype-specific regulatory patterns\n")
            f.write("- Important for understanding biological mechanisms\n")
        
        logger.info(f"Saved heterogeneity interpretation guide to: {guide_file}")
    
    def save_correlation_results(self, peak_strategy: str):
        """
        Save correlation analysis results to files.
        
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
        
        # Save heterogeneity matrix
        if self.heterogeneity_matrix is not None:
            het_file = results_dir / f"heterogeneity_matrix_{peak_strategy}_{timestamp}.csv"
            self.heterogeneity_matrix.to_csv(het_file)
            logger.info(f"Saved heterogeneity matrix to: {het_file}")
        
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
                heterogeneity_value = self.heterogeneity_matrix.iloc[i, j] if self.heterogeneity_matrix is not None else np.nan
                significance_p = self.get_significance_symbol(p_value)
                significance_fdr = self.get_significance_symbol(p_value, fdr_value)
                
                heatmap_data.append({
                    'RNA_cluster': rna_cluster,
                    'ATAC_cluster': atac_cluster,
                    'correlation': correlation,
                    'p_value': p_value,
                    'fdr_value': fdr_value,
                    'heterogeneity_p_value': heterogeneity_value,
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
            
            if self.heterogeneity_matrix is not None:
                f.write("\n\nHeterogeneity P-value Matrix:\n")
                f.write(self.heterogeneity_matrix.to_string())
            
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
                
                # Heterogeneity summary
                if self.heterogeneity_matrix is not None:
                    het_values = pd.to_numeric(self.heterogeneity_matrix.values.flatten(), errors='coerce')
                    valid_hets = het_values[~np.isnan(het_values)]
                    significant_heterogeneity_05 = np.sum(valid_hets < 0.05)
                    
                    f.write(f"\nHeterogeneity Analysis:\n")
                    f.write(f"Significant heterogeneity (p < 0.05): {significant_heterogeneity_05}\n")
                
            except Exception as e:
                f.write(f"Error calculating statistics: {e}\n")
        
        logger.info(f"Saved analysis summary to: {summary_file}")


def main():
    """Main function to run the RNA-ATAC correlation analysis V4 with meta-analysis."""
    parser = argparse.ArgumentParser(
        description="Analyze RNA-seq and ATAC-seq cluster genes overlapping gene expression correlations V4 with meta-analysis"
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
        default="output/RNA_ATAC_correlation_analysis_v4_meta",
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
    analyzer = RNAATACCorrelationAnalyzerV4Meta(
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
