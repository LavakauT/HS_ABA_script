#!/usr/bin/env python3
"""
Chromatin Enrichment Pipeline - Python Implementation
====================================================

Purpose:
    Reproducible, function-based Python pipeline for:
    (1) Preprocessing histone/ATAC peak count matrix for ChromHMM input (CPM)
    (2) Observed/Expected analysis of histone-mark overlaps (with chi-square via 2x2 table)
    (3) ChromHMM state × ATAC cluster enrichment (hypergeometric + heatmap)

Notes:
    - Each function is parameterized for flexible reuse.
    - Every block is step-annotated and documented.
    - Statistical tests are performed from explicit 2x2 contingency tables
      to make logic auditable and easy to cross-check.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Union
import warnings
from scipy import stats
from scipy.stats import chi2_contingency, fisher_exact, hypergeom
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

# Set plotting style
plt.style.use('default')
sns.set_palette("husl")


class ChromatinEnrichmentPipeline:
    """Main class for chromatin enrichment analysis pipeline."""
    
    def __init__(self, output_dir: str = "output"):
        """
        Initialize the pipeline.
        
        Args:
            output_dir: Directory for output files
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Default histone mark levels for plotting order
        self.histone_levels = [
            'H3K4me3', 'H3K36me3', 'H3K9ac', 'H3K14ac', 'H2AZ',
            'H2Aub', 'H3K9me1', 'H3K4me1', 'H3K27me1', 'H3K27me3'
        ]
    
    def ensure_dir(self, path: Union[str, Path]) -> Path:
        """
        Ensure a directory exists (create if missing).
        
        Args:
            path: Directory path
            
        Returns:
            Path object for the directory
        """
        path = Path(path)
        path.mkdir(parents=True, exist_ok=True)
        return path
    
    def contingency_and_tests(self, overlapping: int, total_cluster: int, 
                            total_state: int, total_universe: int) -> Dict:
        """
        Build a 2×2 contingency table for enrichment tests.
        
        Args:
            overlapping: |cluster ∩ state|
            total_cluster: |cluster|
            total_state: |state|
            total_universe: |universe| (background)
            
        Returns:
            Dictionary with contingency table and test results
        """
        a = int(overlapping)
        b = int(total_cluster - overlapping)
        c = int(total_state - overlapping)
        d = int(total_universe - total_cluster - total_state + overlapping)
        
        # Guard rails
        if any(x < 0 for x in [a, b, c, d]):
            raise ValueError(f"Invalid counts for contingency table: got negative cell(s). Check inputs. a={a}, b={b}, c={c}, d={d}")
        
        # Create contingency table
        contingency_table = np.array([[a, b], [c, d]])
        
        # Chi-square test with Yates' correction
        try:
            chi2_stat, chi2_p, dof, expected = chi2_contingency(contingency_table, correction=True)
        except:
            chi2_p = np.nan
        
        # Fisher's exact test (one-sided, enrichment: odds ratio > 1)
        try:
            fisher_p = fisher_exact(contingency_table, alternative='greater')[1]
        except:
            fisher_p = np.nan
        
        # Hypergeometric test (one-sided)
        try:
            # P[X >= a] = 1 - P[X <= a-1]
            hypergeo_p = 1 - hypergeom.cdf(a - 1, total_universe, total_state, total_cluster)
        except:
            hypergeo_p = np.nan
        
        return {
            'contingency_table': contingency_table,
            'chi2_p': chi2_p,
            'fisher_p': fisher_p,
            'hypergeo_p': hypergeo_p
        }
    
    def prepare_counts_for_chromhmm(self, input_counts_file: str, 
                                   output_counts_file: str,
                                   remove_threshold: int = 10000) -> Dict[str, pd.DataFrame]:
        """
        Prepare counts matrix for ChromHMM workflow.
        
        Steps:
        1) Read raw counts
        2) Combine Tak-1 and Upp5 by median (rounded) for specific marks (H3, H3K27me3)
        3) Copy Upp5_H2Aub into H2Aub_Tak1 (as per original logic)
        4) Rename/Drop columns to match ChromHMM inputs
        5) Filter rows where any sample > remove_threshold
        6) Compute CPM by column
        7) Drop rows with all-zero CPM
        
        Args:
            input_counts_file: Path to the counts table (tab-delimited with header)
            output_counts_file: Path to write the modified counts
            remove_threshold: Max read count allowed per row across samples
            
        Returns:
            Dictionary with filtered counts and CPM matrices
        """
        logger.info(f"[Preprocess] Reading counts: {input_counts_file}")
        
        # Read data
        data = pd.read_csv(input_counts_file, sep='\t', index_col=0)
        
        # Step 2: combine by median (rounded)
        if 'H3_Tak1' in data.columns and 'Upp5_H3' in data.columns:
            data['H3_Tak1'] = data[['H3_Tak1', 'Upp5_H3']].median(axis=1).round()
        
        if 'H3K27me3_Tak1' in data.columns and 'Upp5_H3K27me3' in data.columns:
            data['H3K27me3_Tak1'] = data[['H3K27me3_Tak1', 'Upp5_H3K27me3']].median(axis=1).round()
        
        # Step 3: copy Upp5_H2Aub to H2Aub_Tak1
        if 'Upp5_H2Aub' in data.columns:
            data['H2Aub_Tak1'] = data['Upp5_H2Aub']
        
        # Step 4: rename and drop
        data = data.rename(columns={'ATAC.seq_Tak1': 'ATAC_Tak1'})
        
        # Drop Upp5 columns and ATAC_Tak1
        columns_to_drop = ['Upp5_H2Aub', 'Upp5_H3K27me3', 'Upp5_H3', 'ATAC_Tak1']
        data = data.drop(columns=[col for col in columns_to_drop if col in data.columns])
        
        # Write modified counts (for record-keeping)
        output_path = Path(output_counts_file)
        self.ensure_dir(output_path.parent)
        data.to_csv(output_counts_file, sep='\t')
        
        logger.info(f"[Preprocess] Filtering rows with any sample > {remove_threshold}")
        counts_filtered = data[data.apply(lambda row: all(row <= remove_threshold), axis=1)]
        
        logger.info("[Preprocess] Computing CPM")
        # Compute CPM: (counts / total_counts) * 1e6
        cpm = counts_filtered.div(counts_filtered.sum(axis=0), axis=1) * 1e6
        
        # Drop rows with all-zero CPM
        cpm = cpm[cpm.sum(axis=1) > 0]
        
        logger.info(f"[Preprocess] CPM dimension: {cpm.shape[0]} x {cpm.shape[1]}")
        
        return {
            'counts_filtered': counts_filtered,
            'cpm': cpm
        }
    
    def compute_observed_expected(self, pos_dir: str, neg_dir: str, cluster_results_dir: str,
                                histone_levels: Optional[List[str]] = None,
                                out_table: str = "expect_observe.txt",
                                out_plot: str = "chromatin_state_ob_ex_correlation.pdf",
                                out_r_script: str = "chromatin_state_ob_ex_correlation.R") -> pd.DataFrame:
        """
        Observed/Expected enrichment for histone marks across clusters.
        
        Args:
            pos_dir: Directory containing subfolders per cluster with histone mark bed files
            neg_dir: Directory containing per-mark negative overlap bed files
            cluster_results_dir: Directory containing cluster bed files for total peak counts
            histone_levels: Character vector for plotting order
            out_table: Path to write the summary table (TSV)
            out_plot: Path to write the correlation PDF
            out_r_script: Path to write the R script for plot customization
            
        Returns:
            DataFrame with enrichment analysis results
        """
        if histone_levels is None:
            histone_levels = self.histone_levels
        
        logger.info(f"[Ob/Ex] Reading cluster results from: {cluster_results_dir}")
        
        # Read cluster bed files to get total peak counts per cluster
        cluster_path = Path(cluster_results_dir)
        cluster_counts = {}
        
        for bed_file in cluster_path.glob("*.bed"):
            cluster_name = bed_file.stem
            try:
                # Count lines in bed file (excluding header if any)
                with open(bed_file, 'r') as f:
                    line_count = sum(1 for line in f if line.strip() and not line.startswith('#'))
                cluster_counts[cluster_name] = line_count
                logger.info(f"[Ob/Ex] Cluster {cluster_name}: {line_count} peaks")
            except Exception as e:
                logger.warning(f"Error reading cluster file {bed_file}: {e}")
                continue
        
        logger.info(f"[Ob/Ex] Reading POS overlaps from: {pos_dir}")
        
        # Read positive overlaps from histonemark_intersect directory
        pos_path = Path(pos_dir)
        data_list = []
        
        for cluster_dir in pos_path.iterdir():
            if cluster_dir.is_dir():
                cluster_name = cluster_dir.name
                logger.info(f"[Ob/Ex] Processing cluster: {cluster_name}")
                
                # Get total peaks for this cluster
                total_cluster_peaks = cluster_counts.get(cluster_name, 0)
                if total_cluster_peaks == 0:
                    logger.warning(f"No peak count found for cluster {cluster_name}")
                    continue
                
                # Process each histone mark file in the cluster directory
                for histone_file in cluster_dir.glob("*_Tak1_*.bed"):
                    histone_name = histone_file.stem.replace(f'_Tak1_{cluster_name}', '')
                    
                    if histone_name in histone_levels:
                        try:
                            # Count lines in histone mark bed file
                            with open(histone_file, 'r') as f:
                                histone_count = sum(1 for line in f if line.strip() and not line.startswith('#'))
                            
                            data_list.append({
                                'cluster': cluster_name,
                                'histone': histone_name,
                                'total_cluster_peaks': total_cluster_peaks,
                                'histone_marked_peaks': histone_count,
                                'histone_unmarked_peaks': total_cluster_peaks - histone_count
                            })
                            
                            logger.info(f"[Ob/Ex] {cluster_name} - {histone_name}: {histone_count}/{total_cluster_peaks} peaks marked")
                            
                        except Exception as e:
                            logger.warning(f"Error reading histone file {histone_file}: {e}")
                            continue
        
        if not data_list:
            raise ValueError("No valid data found in positive directory")
        
        data = pd.DataFrame(data_list)
        
        logger.info(f"[Ob/Ex] Reading NEG overlaps from: {neg_dir}")
        
        # Read negative overlaps from histonemark_intersect_n directory
        neg_path = Path(neg_dir)
        neg_data_list = []
        
        # Define universe size (total background peaks)
        total_universe = 10092  # Default value, can be parameterized
        
        for histone_file in neg_path.glob("*_Tak1_TN.bed"):
            histone_name = histone_file.stem.replace('_Tak1_TN', '')
            
            if histone_name in histone_levels:
                try:
                    # Count lines in negative histone mark bed file
                    with open(histone_file, 'r') as f:
                        histone_count = sum(1 for line in f if line.strip() and not line.startswith('#'))
                    
                    neg_data_list.append({
                        'histone': histone_name,
                        'total_universe': total_universe,
                        'histone_marked_universe': histone_count,
                        'histone_unmarked_universe': total_universe - histone_count
                    })
                    
                    logger.info(f"[Ob/Ex] Background - {histone_name}: {histone_count}/{total_universe} peaks marked")
                    
                except Exception as e:
                    logger.warning(f"Error reading negative file {histone_file}: {e}")
                    continue
        
        neg_data = pd.DataFrame(neg_data_list)
        
        # Merge positive and negative data
        export = data.merge(neg_data, on='histone', how='inner')
        
        # Calculate percentages and observed/expected ratio
        export['percent_neg'] = (export['histone_marked_universe'] / export['total_universe']) * 100
        export['percent_pos'] = (export['histone_marked_peaks'] / export['total_cluster_peaks']) * 100
        export['ob_ex'] = export['percent_pos'] / export['percent_neg']
        
        # Chi-square via explicit 2×2 table for each row
        logger.info("[Ob/Ex] Running chi-square per (cluster, histone)")
        chi2_results = []
        
        for _, row in export.iterrows():
            # Construct 2×2 table:
            #   - POS group: total_cluster_peaks peaks, with histone_marked_peaks marked
            #   - NEG group: total_universe peaks, with histone_marked_universe marked
            contingency_table = np.array([
                [row['histone_marked_peaks'], row['histone_unmarked_peaks']],
                [row['histone_marked_universe'], row['histone_unmarked_universe']]
            ])
            
            try:
                chi2_stat, chi2_p, dof, expected = chi2_contingency(contingency_table, correction=True)
                chi2_results.append({
                    'chi2_statistic': chi2_stat,
                    'chi2_p_value': chi2_p,
                    'chi2_degrees_of_freedom': dof,
                    'expected_values': expected.tolist()
                })
            except Exception as e:
                logger.warning(f"Chi-square test failed for {row['cluster']}-{row['histone']}: {e}")
                chi2_results.append({
                    'chi2_statistic': np.nan,
                    'chi2_p_value': np.nan,
                    'chi2_degrees_of_freedom': np.nan,
                    'expected_values': [[np.nan, np.nan], [np.nan, np.nan]]
                })
        
        # Add chi-square results to export dataframe
        chi2_df = pd.DataFrame(chi2_results)
        export = pd.concat([export, chi2_df], axis=1)
        
        # Add significance levels
        export['sig'] = export['chi2_p_value'].apply(lambda x: 
            '****' if pd.notna(x) and x < 1e-4 else
            '***' if pd.notna(x) and x < 1e-3 else
            '**' if pd.notna(x) and x < 1e-2 else
            '*' if pd.notna(x) and x < 5e-2 else ''
        )
        
        # Convert to categorical for plotting
        export['cluster'] = pd.Categorical(export['cluster'], 
                                         categories=[f'C{i}' for i in range(1, 16)],
                                         ordered=True)
        export['histone'] = pd.Categorical(export['histone'], 
                                         categories=histone_levels,
                                         ordered=True)
        
        # Create plot
        self._create_ob_ex_plot(export, out_plot)
        
        # Generate R script for plot customization
        self._generate_r_script(export, out_r_script)
        
        # Save detailed table
        output_path = Path(out_table)
        self.ensure_dir(output_path.parent)
        export.to_csv(out_table, sep='\t', index=False)
        
        # Save summary statistics
        summary_path = Path(out_table).parent / "ob_ex_summary_statistics.txt"
        self._save_summary_statistics(export, summary_path)
        
        return export
    
    def _save_summary_statistics(self, data: pd.DataFrame, output_path: str):
        """Save detailed summary statistics for the observed/expected analysis."""
        with open(output_path, 'w') as f:
            f.write("Observed/Expected Analysis Summary Statistics\n")
            f.write("=" * 50 + "\n\n")
            
            # Overall statistics
            f.write(f"Total comparisons: {len(data)}\n")
            f.write(f"Number of clusters: {data['cluster'].nunique()}\n")
            f.write(f"Number of histone marks: {data['histone'].nunique()}\n\n")
            
            # Statistics by cluster
            f.write("Statistics by Cluster:\n")
            f.write("-" * 30 + "\n")
            for cluster in sorted(data['cluster'].unique()):
                cluster_data = data[data['cluster'] == cluster]
                f.write(f"\nCluster {cluster}:\n")
                f.write(f"  Total peaks: {cluster_data['total_cluster_peaks'].iloc[0]}\n")
                f.write(f"  Mean histone mark coverage: {cluster_data['percent_pos'].mean():.2f}%\n")
                f.write(f"  Mean O/E ratio: {cluster_data['ob_ex'].mean():.3f}\n")
                f.write(f"  Significant enrichments: {len(cluster_data[cluster_data['sig'] != ''])}\n")
            
            # Statistics by histone mark
            f.write("\n\nStatistics by Histone Mark:\n")
            f.write("-" * 30 + "\n")
            for histone in self.histone_levels:
                if histone in data['histone'].values:
                    histone_data = data[data['histone'] == histone]
                    f.write(f"\n{histone}:\n")
                    f.write(f"  Background coverage: {histone_data['percent_neg'].iloc[0]:.2f}%\n")
                    f.write(f"  Mean cluster coverage: {histone_data['percent_pos'].mean():.2f}%\n")
                    f.write(f"  Mean O/E ratio: {histone_data['ob_ex'].mean():.3f}\n")
                    f.write(f"  Significant enrichments: {len(histone_data[histone_data['sig'] != ''])}\n")
            
            # Top enrichments
            f.write("\n\nTop 10 Enrichments (by O/E ratio):\n")
            f.write("-" * 40 + "\n")
            top_enrichments = data.nlargest(10, 'ob_ex')[['cluster', 'histone', 'ob_ex', 'percent_pos', 'percent_neg', 'sig']]
            for _, row in top_enrichments.iterrows():
                f.write(f"{row['cluster']}-{row['histone']}: O/E={row['ob_ex']:.3f}, "
                       f"Cluster={row['percent_pos']:.2f}%, Background={row['percent_neg']:.2f}% {row['sig']}\n")
    
    def _generate_r_script(self, data: pd.DataFrame, output_path: str):
        """Generate R script for plot customization."""
        r_script = '''# R script for chromatin state observed/expected correlation plot customization
# Generated by Python chromatin enrichment pipeline

# Load required libraries
library(ggplot2)
library(dplyr)
library(viridis)

# Read the data
data <- read.table("expect_observe.txt", header=TRUE, sep="\\t", stringsAsFactors=FALSE)

# Convert to factors for proper ordering
data$cluster <- factor(data$cluster, levels=paste0("C", 1:15))
data$histone <- factor(data$histone, levels=c("H3K4me3", "H3K36me3", "H3K9ac", "H3K14ac", "H2AZ",
                                             "H2Aub", "H3K9me1", "H3K4me1", "H3K27me1", "H3K27me3"))

# Filter out H3 for plotting (if present)
plot_data <- data[data$histone != "H3", ]

# Create enhanced scatter plot
p1 <- ggplot(plot_data, aes(x=percent_pos, y=ob_ex, color=histone, size=histone_marked_peaks)) +
  geom_point(alpha=0.7) +
  geom_smooth(method="lm", se=TRUE, color="red", linetype="dashed") +
  facet_wrap(~cluster, scales="free_x", ncol=4) +
  scale_color_viridis(discrete=TRUE) +
  scale_size_continuous(range=c(2, 8)) +
  labs(title="Chromatin State Observed/Expected Correlation",
       subtitle="Size indicates number of marked peaks",
       x="Percentage of peaks with histone mark (%)",
       y="Observed/Expected ratio",
       color="Histone Mark",
       size="Marked Peaks") +
  theme_minimal() +
  theme(plot.title=element_text(size=16, face="bold"),
        plot.subtitle=element_text(size=12),
        axis.title=element_text(size=12, face="bold"),
        legend.position="bottom",
        strip.text=element_text(size=10, face="bold"))

# Save the plot
ggsave("chromatin_state_ob_ex_correlation_enhanced.pdf", p1, 
       width=16, height=12, dpi=300)

# Create heatmap of O/E ratios
p2 <- ggplot(plot_data, aes(x=cluster, y=histone, fill=ob_ex)) +
  geom_tile() +
  scale_fill_gradient2(low="blue", mid="white", high="red", 
                       midpoint=1, limits=c(0, max(plot_data$ob_ex, na.rm=TRUE))) +
  geom_text(aes(label=sprintf("%.2f", ob_ex)), size=3, color="black") +
  labs(title="Observed/Expected Ratio Heatmap",
       x="Cluster",
       y="Histone Mark",
       fill="O/E Ratio") +
  theme_minimal() +
  theme(plot.title=element_text(size=14, face="bold"),
        axis.text.x=element_text(angle=45, hjust=1),
        legend.position="right")

# Save the heatmap
ggsave("ob_ex_ratio_heatmap.pdf", p2, width=10, height=8, dpi=300)

# Print summary statistics
cat("\\nSummary Statistics:\\n")
cat("==================\\n")
cat("Total comparisons:", nrow(data), "\\n")
cat("Number of clusters:", length(unique(data$cluster)), "\\n")
cat("Number of histone marks:", length(unique(data$histone)), "\\n")

# Summary by cluster
cat("\\nSummary by Cluster:\\n")
cluster_summary <- data %>%
  group_by(cluster) %>%
  summarise(
    total_peaks = first(total_cluster_peaks),
    mean_coverage = mean(percent_pos, na.rm=TRUE),
    mean_oe_ratio = mean(ob_ex, na.rm=TRUE),
    significant_count = sum(sig != "", na.rm=TRUE)
  )
print(cluster_summary)

# Summary by histone mark
cat("\\nSummary by Histone Mark:\\n")
histone_summary <- data %>%
  group_by(histone) %>%
  summarise(
    background_coverage = first(percent_neg),
    mean_cluster_coverage = mean(percent_pos, na.rm=TRUE),
    mean_oe_ratio = mean(ob_ex, na.rm=TRUE),
    significant_count = sum(sig != "", na.rm=TRUE)
  )
print(histone_summary)
'''
        
        output_path = Path(output_path)
        self.ensure_dir(output_path.parent)
        with open(output_path, 'w') as f:
            f.write(r_script)
        
        logger.info(f"R script generated: {output_path}")
    
    def _create_ob_ex_plot(self, data: pd.DataFrame, out_plot: str):
        """Create the observed/expected correlation plot."""
        # Filter out H3 for plotting
        plot_data = data[data['histone'] != 'H3'].copy()
        
        # Create facet plot
        fig, axes = plt.subplots(2, 8, figsize=(20, 8))
        fig.suptitle('Chromatin State Observed/Expected Correlation', fontsize=16, fontweight='bold')
        
        # Flatten axes for easier iteration
        axes_flat = axes.flatten()
        
        for i, cluster in enumerate(sorted(plot_data['cluster'].unique())):
            if i < len(axes_flat):
                ax = axes_flat[i]
                cluster_data = plot_data[plot_data['cluster'] == cluster]
                
                # Scatter plot with size based on number of marked peaks
                scatter = ax.scatter(cluster_data['percent_pos'], cluster_data['ob_ex'], 
                                   c=cluster_data['histone'].cat.codes, 
                                   s=cluster_data['histone_marked_peaks']/10, alpha=0.7)
                
                # Add trend line
                if len(cluster_data) > 1:
                    z = np.polyfit(cluster_data['percent_pos'], cluster_data['ob_ex'], 1)
                    p = np.poly1d(z)
                    ax.plot(cluster_data['percent_pos'], p(cluster_data['percent_pos']), 
                           "r--", alpha=0.8)
                
                ax.set_title(f'Cluster {cluster}')
                ax.set_xlabel('Percentage (%)')
                ax.set_ylabel('Observed/Expected')
                ax.grid(True, alpha=0.3)
        
        # Hide unused subplots
        for i in range(len(plot_data['cluster'].unique()), len(axes_flat)):
            axes_flat[i].set_visible(False)
        
        plt.tight_layout()
        
        # Save plot
        output_path = Path(out_plot)
        self.ensure_dir(output_path.parent)
        plt.savefig(out_plot, dpi=300, bbox_inches='tight')
        plt.close()
    
    def chromhmm_state_enrichment(self, peaks_file: str, segments_bed: str,
                                 overlap_counts_dir: str,
                                 out_pdf: str = "state_ovl.pdf",
                                 out_r_script: str = "state_ovl.R") -> pd.DataFrame:
        """
        ChromHMM state × ATAC cluster enrichment with hypergeometric test and heatmap.
        
        Args:
            peaks_file: Path to ATAC DE peaks w/ cluster labels (expects column "cluster")
            segments_bed: ChromHMM segments BED (Chr, Start, End, state)
            overlap_counts_dir: Directory with files like C1_state_count.txt
            out_pdf: Path to write the ComplexHeatmap PDF
            out_r_script: Path to write the R script for plot customization
            
        Returns:
            DataFrame with enrichment analysis results
        """
        logger.info("[State×Cluster] Reading inputs")
        
        # Read peaks file
        peaks = pd.read_csv(peaks_file, sep='\t')
        
        # Read segments BED file
        segments = pd.read_csv(segments_bed, header=None, sep='\t', 
                              names=['Chr', 'Start', 'End', 'state'])
        
        # Calculate totals
        total_s = segments.groupby('state').size().reset_index(name='total_counts_state')
        total = peaks.groupby('cluster').size().reset_index(name='total_counts')
        
        # Read overlap count files
        overlap_path = Path(overlap_counts_dir)
        so_files = list(overlap_path.glob("*.txt"))
        
        sm_list = []
        for file_path in so_files:
            cluster = file_path.stem.replace('__state_count', '')
            
            try:
                # Read space-separated file without header
                df = pd.read_csv(file_path, header=None, sep='\s+', names=['counts', 'state'])
                df = df[df['state'] != 'total'].copy()
                
                # Clean state names - remove .bed extension and cluster suffix
                df['state'] = df['state'].str.replace(r'\.bed$', '', regex=True)
                df['state'] = df['state'].str.replace(r'_C\d+$', '', regex=True)
                df['cluster'] = cluster
                df['overlapping'] = df['counts']
                
                sm_list.append(df[['cluster', 'state', 'overlapping']])
                logger.info(f"[State×Cluster] Read {len(df)} states for cluster {cluster}")
            except Exception as e:
                logger.warning(f"Error reading overlap file {file_path}: {e}")
                continue
        
        if not sm_list:
            raise ValueError("No valid overlap data found")
        
        sm = pd.concat(sm_list, ignore_index=True)
        
        # Calculate universe size
        Total = total_s['total_counts_state'].sum()
        
        # Merge with totals and calculate enrichment
        sm2 = sm.merge(total, on='cluster').merge(total_s, on='state')
        
        sm2['cluster_ratio'] = (sm2['overlapping'] / sm2['total_counts']).round(2)
        
        # Hypergeometric test (one-sided, enrichment)
        sm2['hg'] = sm2.apply(lambda row: 
            1 - hypergeom.cdf(row['overlapping'] - 1, Total, 
                             row['total_counts_state'], row['total_counts']), axis=1)
        
        sm2['log10(hg)'] = -np.log10(sm2['hg'])
        sm2['mark'] = sm2['hg'].apply(lambda x: '*' if x < 0.05 else '')
        
        # Convert to categorical for plotting
        sm2['state'] = pd.Categorical(sm2['state'], 
                                     categories=[f'E{i}' for i in range(1, 13)],
                                     ordered=True)
        sm2['cluster'] = pd.Categorical(sm2['cluster'], 
                                       categories=[f'C{i}' for i in range(1, 15)],
                                       ordered=True)
        
        # Create heatmap
        self._create_state_enrichment_heatmap(sm2, out_pdf)
        
        # Generate R script for plot customization
        self._generate_state_enrichment_r_script(sm2, out_r_script)
        
        # Save detailed results
        detailed_output = Path(out_pdf).parent / "state_enrichment_detailed_results.txt"
        self._save_state_enrichment_statistics(sm2, detailed_output)
        
        return sm2
    
    def _generate_state_enrichment_r_script(self, data: pd.DataFrame, output_path: str):
        """Generate R script for state enrichment plot customization."""
        r_script = '''# R script for ChromHMM state × ATAC cluster enrichment plot customization
# Generated by Python chromatin enrichment pipeline
# This script creates plots similar to state_ovl.pdf

# Load required libraries
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(ggplot2)
library(viridis)
library(gridExtra)

# Set working directory to script location (optional)
# setwd(dirname(parent.frame(2)$ofile))

# Read the data - multiple format options
data_file <- "state_enrichment_detailed_results.txt"
if (file.exists(data_file)) {
  data <- read.table(data_file, header=TRUE, sep="\\t", stringsAsFactors=FALSE)
} else {
  # Try CSV format
  csv_file <- "state_enrichment_detailed_results.csv"
  if (file.exists(csv_file)) {
    data <- read.csv(csv_file, sep="\\t", stringsAsFactors=FALSE)
  } else {
    stop("No data file found. Please ensure state_enrichment_detailed_results.txt or .csv exists.")
  }
}

# Also try to read pre-made matrices if available
ratio_matrix_file <- "ratio_matrix.csv"
count_matrix_file <- "count_matrix.csv"
pvalue_matrix_file <- "pvalue_matrix.csv"
sig_matrix_file <- "significance_matrix.csv"

# Function to safely read matrix files
read_matrix_safe <- function(filename) {
  if (file.exists(filename)) {
    return(as.matrix(read.csv(filename, row.names=1)))
  } else {
    return(NULL)
  }
}

# Read matrices if available, otherwise create from data
ratio_matrix <- read_matrix_safe(ratio_matrix_file)
count_matrix <- read_matrix_safe(count_matrix_file)
pvalue_matrix <- read_matrix_safe(pvalue_matrix_file)
sig_matrix <- read_matrix_safe(sig_matrix_file)

# If matrices not available, create from data
if (is.null(ratio_matrix)) {
  ratio_matrix <- data %>%
    select(cluster, state, cluster_ratio) %>%
    pivot_wider(names_from=cluster, values_from=cluster_ratio) %>%
    column_to_rownames("state") %>%
    as.matrix()
}

if (is.null(count_matrix)) {
  count_matrix <- data %>%
    select(cluster, state, overlapping) %>%
    pivot_wider(names_from=cluster, values_from=overlapping) %>%
    column_to_rownames("state") %>%
    as.matrix()
}

if (is.null(pvalue_matrix)) {
  pvalue_matrix <- data %>%
    select(cluster, state, hg) %>%
    pivot_wider(names_from=cluster, values_from=hg) %>%
    column_to_rownames("state") %>%
    as.matrix()
}

if (is.null(sig_matrix)) {
  sig_matrix <- ifelse(pvalue_matrix < 0.05, "*", "")
}

# Ensure proper ordering of states and clusters
state_order <- paste0("E", 1:12)
cluster_order <- paste0("C", 1:14)

# Filter and order matrices
ratio_matrix <- ratio_matrix[state_order, cluster_order, drop=FALSE]
count_matrix <- count_matrix[state_order, cluster_order, drop=FALSE]
pvalue_matrix <- pvalue_matrix[state_order, cluster_order, drop=FALSE]
sig_matrix <- sig_matrix[state_order, cluster_order, drop=FALSE]

# Create enhanced heatmap similar to state_ovl.pdf
ht <- Heatmap(ratio_matrix,
              name="Cluster Ratio",
              col=colorRamp2(c(0, 0.25, 0.5), c("blue", "white", "red")),
              cluster_rows=FALSE,
              cluster_columns=FALSE,
              row_names_gp=gpar(fontsize=10, fontface="bold"),
              column_names_gp=gpar(fontsize=10, fontface="bold"),
              row_title="Chromatin State",
              column_title="ATAC Cluster",
              row_title_gp=gpar(fontsize=12, fontface="bold"),
              column_title_gp=gpar(fontsize=12, fontface="bold"),
              cell_fun=function(j, i, x, y, width, height, fill) {
                # Add ratio value
                grid.text(sprintf("%.3f", ratio_matrix[i, j]), x, y, gp=gpar(fontsize=8))
                # Add significance marker
                if(sig_matrix[i, j] != "") {
                  grid.text(sig_matrix[i, j], x, y+0.3, gp=gpar(fontsize=12, fontface="bold", col="black"))
                }
              },
              heatmap_legend_param=list(
                title_position="topcenter", 
                legend_height=unit(4, "cm"),
                title_gp=gpar(fontsize=12, fontface="bold")
              ))

# Save the main heatmap
pdf("state_enrichment_enhanced.pdf", width=14, height=10)
draw(ht)
dev.off()

# Create alternative heatmap with counts as annotation
ht_with_counts <- Heatmap(ratio_matrix,
                          name="Cluster Ratio",
                          col=colorRamp2(c(0, 0.25, 0.5), c("blue", "white", "red")),
                          cluster_rows=FALSE,
                          cluster_columns=FALSE,
                          row_names_gp=gpar(fontsize=10, fontface="bold"),
                          column_names_gp=gpar(fontsize=10, fontface="bold"),
                          row_title="Chromatin State",
                          column_title="ATAC Cluster",
                          row_title_gp=gpar(fontsize=12, fontface="bold"),
                          column_title_gp=gpar(fontsize=12, fontface="bold"),
                          cell_fun=function(j, i, x, y, width, height, fill) {
                            # Add count value
                            grid.text(count_matrix[i, j], x, y-0.3, gp=gpar(fontsize=8, col="darkblue"))
                            # Add ratio value
                            grid.text(sprintf("%.3f", ratio_matrix[i, j]), x, y, gp=gpar(fontsize=8))
                            # Add significance marker
                            if(sig_matrix[i, j] != "") {
                              grid.text(sig_matrix[i, j], x, y+0.3, gp=gpar(fontsize=12, fontface="bold", col="black"))
                            }
                          },
                          heatmap_legend_param=list(
                            title_position="topcenter", 
                            legend_height=unit(4, "cm"),
                            title_gp=gpar(fontsize=12, fontface="bold")
                          ))

# Save the count-annotated heatmap
pdf("state_enrichment_with_counts.pdf", width=14, height=10)
draw(ht_with_counts)
dev.off()

# Create alternative visualizations with ggplot2
# Scatter plot of enrichment
p1 <- ggplot(data, aes(x=cluster, y=state, size=overlapping, color=cluster_ratio)) +
  geom_point() +
  scale_size_continuous(range=c(2, 12)) +
  scale_color_viridis(option="plasma", direction=-1) +
  labs(title="ChromHMM State × ATAC Cluster Enrichment",
       subtitle="Size indicates overlap count, color indicates cluster ratio",
       x="ATAC Cluster",
       y="Chromatin State",
       size="Overlap Count",
       color="Cluster Ratio") +
  theme_minimal() +
  theme(plot.title=element_text(size=14, face="bold"),
        axis.text.x=element_text(angle=45, hjust=1),
        legend.position="right")

ggsave("state_enrichment_scatter.pdf", p1, width=12, height=8, dpi=300)

# Bar plot of cluster ratios by state
p2 <- ggplot(data, aes(x=cluster, y=cluster_ratio, fill=state)) +
  geom_bar(stat="identity", position="dodge") +
  scale_fill_viridis(discrete=TRUE, option="plasma") +
  labs(title="Cluster Ratios by Chromatin State",
       x="ATAC Cluster",
       y="Cluster Ratio",
       fill="Chromatin State") +
  theme_minimal() +
  theme(plot.title=element_text(size=14, face="bold"),
        axis.text.x=element_text(angle=45, hjust=1))

ggsave("state_enrichment_barplot.pdf", p2, width=14, height=8, dpi=300)

# Heatmap using ggplot2 (alternative to ComplexHeatmap)
p3 <- ggplot(data, aes(x=cluster, y=state, fill=cluster_ratio)) +
  geom_tile() +
  scale_fill_viridis(option="plasma") +
  geom_text(aes(label=sprintf("%.3f", cluster_ratio)), size=3) +
  geom_text(aes(label=ifelse(hg < 0.05, "*", "")), 
            position=position_nudge(y=0.3), size=6, fontface="bold") +
  labs(title="ChromHMM State × ATAC Cluster Enrichment (ggplot2)",
       x="ATAC Cluster",
       y="Chromatin State",
       fill="Cluster Ratio") +
  theme_minimal() +
  theme(plot.title=element_text(size=14, face="bold"),
        axis.text.x=element_text(angle=45, hjust=1))

ggsave("state_enrichment_ggplot2_heatmap.pdf", p3, width=14, height=8, dpi=300)

# Print summary statistics
cat("\\nState Enrichment Summary Statistics:\\n")
cat("=====================================\\n")
cat("Total comparisons:", nrow(data), "\\n")
cat("Number of clusters:", length(unique(data$cluster)), "\\n")
cat("Number of states:", length(unique(data$state)), "\\n")

# Summary by cluster
cat("\\nSummary by Cluster:\\n")
cluster_summary <- data %>%
  group_by(cluster) %>%
  summarise(
    total_peaks = first(total_counts),
    mean_ratio = mean(cluster_ratio, na.rm=TRUE),
    significant_states = sum(hg < 0.05, na.rm=TRUE),
    max_ratio = max(cluster_ratio, na.rm=TRUE),
    min_ratio = min(cluster_ratio, na.rm=TRUE)
  )
print(cluster_summary)

# Summary by state
cat("\\nSummary by State:\\n")
state_summary <- data %>%
  group_by(state) %>%
  summarise(
    total_segments = first(total_counts_state),
    mean_ratio = mean(cluster_ratio, na.rm=TRUE),
    significant_clusters = sum(hg < 0.05, na.rm=TRUE),
    max_ratio = max(cluster_ratio, na.rm=TRUE),
    min_ratio = min(cluster_ratio, na.rm=TRUE)
  )
print(state_summary)

# Top enrichments
cat("\\nTop 10 Enrichments (by cluster ratio):\\n")
top_enrichments <- data %>%
  arrange(desc(cluster_ratio)) %>%
  head(10) %>%
  select(cluster, state, cluster_ratio, overlapping, hg)
print(top_enrichments)

# Save summary statistics
write.csv(cluster_summary, "cluster_summary.csv", row.names=FALSE)
write.csv(state_summary, "state_summary.csv", row.names=FALSE)
write.csv(top_enrichments, "top_enrichments.csv", row.names=FALSE)

cat("\\nAll plots and data files have been saved.\\n")
cat("Main heatmap: state_enrichment_enhanced.pdf\\n")
cat("Count-annotated heatmap: state_enrichment_with_counts.pdf\\n")
cat("Scatter plot: state_enrichment_scatter.pdf\\n")
cat("Bar plot: state_enrichment_barplot.pdf\\n")
cat("ggplot2 heatmap: state_enrichment_ggplot2_heatmap.pdf\\n")
'''
        
        output_path = Path(output_path)
        self.ensure_dir(output_path.parent)
        with open(output_path, 'w') as f:
            f.write(r_script)
        
        logger.info(f"Enhanced state enrichment R script generated: {output_path}")
    
    def _save_state_enrichment_statistics(self, data: pd.DataFrame, output_path: str):
        """Save detailed statistics for the state enrichment analysis."""
        with open(output_path, 'w') as f:
            f.write("ChromHMM State × ATAC Cluster Enrichment Analysis\n")
            f.write("=" * 55 + "\n\n")
            
            # Overall statistics
            f.write(f"Total comparisons: {len(data)}\n")
            f.write(f"Number of clusters: {data['cluster'].nunique()}\n")
            f.write(f"Number of states: {data['state'].nunique()}\n\n")
            
            # Statistics by cluster
            f.write("Statistics by Cluster:\n")
            f.write("-" * 30 + "\n")
            for cluster in sorted(data['cluster'].unique()):
                cluster_data = data[data['cluster'] == cluster]
                f.write(f"\nCluster {cluster}:\n")
                f.write(f"  Total peaks: {cluster_data['total_counts'].iloc[0]}\n")
                f.write(f"  Mean cluster ratio: {cluster_data['cluster_ratio'].mean():.3f}\n")
                f.write(f"  Significant states: {len(cluster_data[cluster_data['hg'] < 0.05])}\n")
            
            # Statistics by state
            f.write("\n\nStatistics by State:\n")
            f.write("-" * 30 + "\n")
            for state in sorted(data['state'].unique()):
                state_data = data[data['state'] == state]
                f.write(f"\nState {state}:\n")
                f.write(f"  Total segments: {state_data['total_counts_state'].iloc[0]}\n")
                f.write(f"  Mean cluster ratio: {state_data['cluster_ratio'].mean():.3f}\n")
                f.write(f"  Significant clusters: {len(state_data[state_data['hg'] < 0.05])}\n")
            
            # Top enrichments
            f.write("\n\nTop 10 Enrichments (by cluster ratio):\n")
            f.write("-" * 40 + "\n")
            top_enrichments = data.nlargest(10, 'cluster_ratio')[['cluster', 'state', 'cluster_ratio', 'overlapping', 'hg']]
            for _, row in top_enrichments.iterrows():
                f.write(f"{row['cluster']}-{row['state']}: Ratio={row['cluster_ratio']:.3f}, "
                       f"Overlap={row['overlapping']}, p-value={row['hg']:.2e}\n")
            
            # Raw data table
            f.write("\n\nRaw Data Table (Tab-separated):\n")
            f.write("-" * 40 + "\n")
            f.write("cluster\tstate\toverlapping\ttotal_counts\ttotal_counts_state\tcluster_ratio\thg\tlog10(hg)\tmark\n")
            for _, row in data.iterrows():
                f.write(f"{row['cluster']}\t{row['state']}\t{row['overlapping']}\t"
                       f"{row['total_counts']}\t{row['total_counts_state']}\t"
                       f"{row['cluster_ratio']:.4f}\t{row['hg']:.2e}\t"
                       f"{row['log10(hg)']:.3f}\t{row['mark']}\n")
        
        # Also save as CSV for easier R import
        csv_output_path = str(output_path).replace('.txt', '.csv')
        data.to_csv(csv_output_path, sep='\t', index=False)
        logger.info(f"Raw data saved as CSV: {csv_output_path}")
        
        # Save matrices for R plotting
        self._save_matrices_for_r(data, output_path)
    
    def _save_matrices_for_r(self, data: pd.DataFrame, base_output_path: str):
        """Save data matrices in formats optimized for R plotting."""
        base_path = Path(base_output_path).parent
        
        # Create ratio matrix
        ratio_matrix = data.pivot(index='state', columns='cluster', values='cluster_ratio')
        ratio_matrix.to_csv(base_path / "ratio_matrix.csv")
        
        # Create count matrix
        count_matrix = data.pivot(index='state', columns='cluster', values='overlapping')
        count_matrix.to_csv(base_path / "count_matrix.csv")
        
        # Create p-value matrix
        pvalue_matrix = data.pivot(index='state', columns='cluster', values='hg')
        pvalue_matrix.to_csv(base_path / "pvalue_matrix.csv")
        
        # Create significance matrix
        sig_matrix = data.pivot(index='state', columns='cluster', values='mark')
        sig_matrix.to_csv(base_path / "significance_matrix.csv")
        
        # Create summary statistics
        summary_stats = {
            'total_comparisons': len(data),
            'num_clusters': data['cluster'].nunique(),
            'num_states': data['state'].nunique(),
            'total_peaks': data['total_counts'].sum(),
            'total_segments': data['total_counts_state'].sum()
        }
        
        with open(base_path / "summary_stats.txt", 'w') as f:
            for key, value in summary_stats.items():
                f.write(f"{key}: {value}\n")
        
        logger.info(f"Matrix files saved for R plotting in: {base_path}")
    
    def _create_state_enrichment_heatmap(self, data: pd.DataFrame, out_pdf: str):
        """Create the state enrichment heatmap."""
        # Prepare matrices for heatmap
        df_p = data.pivot(index='state', columns='cluster', values='cluster_ratio')
        df_c = data.pivot(index='state', columns='cluster', values='overlapping')
        df_s = data.pivot(index='state', columns='cluster', values='mark')
        
        # Create heatmap
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Create heatmap using seaborn
        sns.heatmap(df_p, annot=df_c, fmt='g', cmap='RdBu_r', center=0.25,
                   cbar_kws={'label': 'Cluster Ratio'}, ax=ax)
        
        # Add significance markers
        for i, state in enumerate(df_s.index):
            for j, cluster in enumerate(df_s.columns):
                if pd.notna(df_s.loc[state, cluster]) and df_s.loc[state, cluster] == '*':
                    ax.text(j + 0.5, i + 0.5, '*', ha='center', va='center', 
                           fontsize=12, fontweight='bold', color='black')
        
        ax.set_title('ChromHMM State × ATAC Cluster Enrichment', fontsize=14, fontweight='bold')
        ax.set_xlabel('ATAC Cluster', fontsize=12, fontweight='bold')
        ax.set_ylabel('Chromatin State', fontsize=12, fontweight='bold')
        
        plt.tight_layout()
        
        # Save plot
        output_path = Path(out_pdf)
        self.ensure_dir(output_path.parent)
        plt.savefig(out_pdf, dpi=300, bbox_inches='tight')
        plt.close()
    
    def run_pipeline(self, counts_in: str = "data/chromatin/histone_merge_peak_counts.txt",
                    counts_out: str = "output/chromatin/modified_histone_merge_peak_counts.txt",
                    pos_dir: str = "data/chromatin/histonemark_intersect",
                    neg_dir: str = "data/chromatin/histonemark_intersect_n",
                    cluster_results_dir: str = "data/chromatin/cluster_results",
                    obex_table: str = "output/chromatin/expect_observe.txt",
                    obex_plot: str = "output/chromatin/chromatin_state_ob_ex_correlation.pdf",
                    obex_r_script: str = "output/chromatin/chromatin_state_ob_ex_correlation.R",
                    peaks_file: str = "data/chromatin/DE_peaks_hclust.txt",
                    segments_bed: str = "data/chromatin/Tak1_12_segments.bed",
                    overlap_counts_dir: str = "data/chromatin/state_cluster_ovl",
                    state_plot: str = "output/chromatin/state_ovl.pdf",
                    state_r_script: str = "output/chromatin/state_ovl.R"):
        """
        Run the full pipeline.
        
        Args:
            counts_in: Input counts file path
            counts_out: Output counts file path
            pos_dir: Positive overlaps directory
            neg_dir: Negative overlaps directory
            cluster_results_dir: Cluster results directory with bed files
            obex_table: Observed/expected table output path
            obex_plot: Observed/expected plot output path
            obex_r_script: Observed/expected R script output path
            peaks_file: Peaks file path
            segments_bed: Segments BED file path
            overlap_counts_dir: Overlap counts directory
            state_plot: State enrichment plot output path
            state_r_script: State enrichment R script output path
        """
        try:
            # (1) Preprocess & CPM
            logger.info("Starting preprocessing and CPM calculation...")
            self.prepare_counts_for_chromhmm(counts_in, counts_out, remove_threshold=10000)
            
            # (2) Observed/Expected
            logger.info("Starting observed/expected analysis...")
            self.compute_observed_expected(
                pos_dir=pos_dir,
                neg_dir=neg_dir,
                cluster_results_dir=cluster_results_dir,
                out_table=obex_table,
                out_plot=obex_plot,
                out_r_script=obex_r_script
            )
            
            # (3) State × Cluster enrichment
            logger.info("Starting state enrichment analysis...")
            self.chromhmm_state_enrichment(
                peaks_file=peaks_file,
                segments_bed=segments_bed,
                overlap_counts_dir=overlap_counts_dir,
                out_pdf=state_plot,
                out_r_script=state_r_script
            )
            
            logger.info("Pipeline completed successfully!")
            
        except Exception as e:
            logger.error(f"Pipeline failed: {e}")
            raise


def main():
    """Main function to run the pipeline."""
    # Initialize pipeline
    pipeline = ChromatinEnrichmentPipeline()
    
    # Run with corrected paths
    pipeline.run_pipeline()


if __name__ == "__main__":
    main()
