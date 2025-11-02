#!/usr/bin/env python3
"""
RNA-seq and ATAC-seq Cluster Gene Overlap Enrichment Analysis

This script analyzes the overlap between RNA-seq cluster genes and ATAC-seq cluster genes,
calculates enrichment statistics using hypergeometric test, and generates heatmaps for visualization.

Author: Generated for HS_ABA project
Date: 2025
"""

import os
import sys
import argparse
import logging
from pathlib import Path
from typing import Dict, List, Set, Tuple, Optional, Union
from collections import defaultdict
import warnings

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import hypergeom
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as patches

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore', category=FutureWarning)


class ClusterOverlapAnalyzer:
    """
    Analyzer for RNA-seq and ATAC-seq cluster gene overlap enrichment analysis.
    
    This class handles loading cluster gene data, calculating overlaps,
    performing hypergeometric enrichment tests, and generating visualizations.
    """
    
    def __init__(self, output_dir: str = "output/cluster_overlap_analysis"):
        """
        Initialize the analyzer.
        
        Args:
            output_dir: Directory to save output files
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Data storage
        self.rna_clusters: Dict[str, Set[str]] = {}
        self.atac_clusters: Dict[str, Set[str]] = {}
        self.universal_genes: Set[str] = set()
        
        # Results storage
        self.overlap_matrix: Optional[pd.DataFrame] = None
        self.proportion_matrix: Optional[pd.DataFrame] = None
        self.pvalue_matrix: Optional[pd.DataFrame] = None
        self.enrichment_matrix: Optional[pd.DataFrame] = None
    
    def load_rna_clusters(self, rna_cluster_dir: str) -> None:
        """
        Load RNA-seq cluster gene files.
        
        Args:
            rna_cluster_dir: Directory containing RNA-seq cluster gene files
        """
        logger.info(f"Loading RNA-seq cluster genes from: {rna_cluster_dir}")
        
        cluster_dir = Path(rna_cluster_dir)
        if not cluster_dir.exists():
            raise FileNotFoundError(f"RNA cluster directory not found: {rna_cluster_dir}")
        
        # Find all cluster gene files
        gene_files = list(cluster_dir.glob("cluster_*_genes.txt"))
        if not gene_files:
            raise FileNotFoundError(f"No cluster gene files found in: {rna_cluster_dir}")
        
        for file_path in sorted(gene_files):
            cluster_name = file_path.stem.replace("_genes", "")
            logger.info(f"Loading RNA cluster: {cluster_name}")
            
            genes = self._load_gene_file(file_path, skip_header=True)
            self.rna_clusters[cluster_name] = genes
            self.universal_genes.update(genes)
        
        logger.info(f"Loaded {len(self.rna_clusters)} RNA-seq clusters")
        for cluster, genes in self.rna_clusters.items():
            logger.info(f"  {cluster}: {len(genes)} genes")
    
    def load_atac_clusters(self, atac_cluster_dir: str) -> None:
        """
        Load ATAC-seq cluster gene files.
        
        Args:
            atac_cluster_dir: Directory containing ATAC-seq cluster gene files
        """
        logger.info(f"Loading ATAC-seq cluster genes from: {atac_cluster_dir}")
        
        cluster_dir = Path(atac_cluster_dir)
        if not cluster_dir.exists():
            raise FileNotFoundError(f"ATAC cluster directory not found: {atac_cluster_dir}")
        
        # Try different file naming patterns
        file_patterns = [
            "cluster_*_peaks_all_peaks_genes.txt",      # all_peaks format
            "cluster_*_peaks_proximal_peaks_genes.txt", # proximal_peaks format
            "cluster_*_peaks_distal_peaks_genes.txt",   # distal_peaks format
            "cluster_*_peaks_genes.txt",                # generic format
            "cluster_*_genes.txt"                       # simple format
        ]
        
        gene_files = []
        detected_pattern = None
        
        for pattern in file_patterns:
            gene_files = list(cluster_dir.glob(pattern))
            if gene_files:
                detected_pattern = pattern
                logger.info(f"Detected file pattern: {pattern}")
                break
        
        if not gene_files:
            available_files = list(cluster_dir.glob("*.txt"))
            available_names = [f.name for f in available_files[:5]]  # Show first 5 files
            raise FileNotFoundError(
                f"No cluster gene files found in: {atac_cluster_dir}\n"
                f"Available .txt files: {available_names}\n"
                f"Expected patterns: {file_patterns}"
            )
        
        # Determine suffix to remove based on detected pattern
        if "_peaks_all_peaks_genes" in detected_pattern:
            suffix_to_remove = "_peaks_all_peaks_genes"
        elif "_peaks_proximal_peaks_genes" in detected_pattern:
            suffix_to_remove = "_peaks_proximal_peaks_genes"
        elif "_peaks_distal_peaks_genes" in detected_pattern:
            suffix_to_remove = "_peaks_distal_peaks_genes"
        elif "_peaks_genes" in detected_pattern:
            suffix_to_remove = "_peaks_genes"
        else:
            suffix_to_remove = "_genes"
        
        for file_path in sorted(gene_files):
            cluster_name = file_path.stem.replace(suffix_to_remove, "")
            logger.info(f"Loading ATAC cluster: {cluster_name}")
            
            genes = self._load_gene_file(file_path, skip_header=True, skip_comments=True)
            self.atac_clusters[cluster_name] = genes
            self.universal_genes.update(genes)
        
        logger.info(f"Loaded {len(self.atac_clusters)} ATAC-seq clusters")
        for cluster, genes in self.atac_clusters.items():
            logger.info(f"  {cluster}: {len(genes)} genes")
    
    def _load_gene_file(self, file_path: Path, skip_header: bool = True, 
                       skip_comments: bool = False) -> Set[str]:
        """
        Load genes from a file.
        
        Args:
            file_path: Path to the gene file
            skip_header: Whether to skip the first line
            skip_comments: Whether to skip comment lines starting with #
            
        Returns:
            Set of gene IDs
        """
        genes = set()
        
        with open(file_path, 'r') as f:
            lines = f.readlines()
        
        start_idx = 1 if skip_header else 0
        
        for line in lines[start_idx:]:
            line = line.strip()
            if not line:
                continue
            if skip_comments and line.startswith('#'):
                continue
            
            # Extract gene ID (assuming first column or entire line)
            gene_id = line.split('\t')[0] if '\t' in line else line
            if gene_id and gene_id != "Gene" and gene_id != "geneId":
                genes.add(gene_id)
        
        return genes
    
    def calculate_overlap_enrichment(self) -> None:
        """
        Calculate overlap counts, proportions, and enrichment p-values.
        """
        logger.info("Calculating overlap enrichment analysis...")
        
        if not self.rna_clusters or not self.atac_clusters:
            raise ValueError("Both RNA and ATAC clusters must be loaded first")
        
        # Initialize matrices
        rna_cluster_names = sorted(self.rna_clusters.keys())
        atac_cluster_names = sorted(self.atac_clusters.keys())
        
        overlap_data = []
        proportion_data = []
        pvalue_data = []
        enrichment_data = []
        
        universe_size = len(self.universal_genes)
        logger.info(f"Universal gene set size: {universe_size}")
        
        for rna_cluster in rna_cluster_names:
            overlap_row = []
            proportion_row = []
            pvalue_row = []
            enrichment_row = []
            
            rna_genes = self.rna_clusters[rna_cluster]
            rna_size = len(rna_genes)
            
            for atac_cluster in atac_cluster_names:
                atac_genes = self.atac_clusters[atac_cluster]
                atac_size = len(atac_genes)
                
                # Calculate overlap
                overlap_genes = rna_genes.intersection(atac_genes)
                overlap_count = len(overlap_genes)
                
                # Calculate proportion (overlap / RNA cluster size)
                proportion = overlap_count / rna_size if rna_size > 0 else 0
                
                # Hypergeometric test for enrichment
                # phyper(q, m, n, k, lower.tail=FALSE) in R
                # q: number of successes (overlap_count - 1 for upper tail)
                # m: number of success states in population (atac_size)
                # n: number of failure states in population (universe_size - atac_size)
                # k: number of draws (rna_size)
                
                if overlap_count > 0 and rna_size > 0 and atac_size > 0:
                    pvalue = hypergeom.sf(overlap_count - 1, universe_size, atac_size, rna_size)
                    # Calculate enrichment fold change
                    expected = (rna_size * atac_size) / universe_size
                    enrichment = overlap_count / expected if expected > 0 else 0
                else:
                    pvalue = 1.0
                    enrichment = 0.0
                
                overlap_row.append(overlap_count)
                proportion_row.append(proportion)
                pvalue_row.append(pvalue)
                enrichment_row.append(enrichment)
            
            overlap_data.append(overlap_row)
            proportion_data.append(proportion_row)
            pvalue_data.append(pvalue_row)
            enrichment_data.append(enrichment_row)
        
        # Create DataFrames
        self.overlap_matrix = pd.DataFrame(
            overlap_data, 
            index=rna_cluster_names, 
            columns=atac_cluster_names
        )
        
        self.proportion_matrix = pd.DataFrame(
            proportion_data, 
            index=rna_cluster_names, 
            columns=atac_cluster_names
        )
        
        self.pvalue_matrix = pd.DataFrame(
            pvalue_data, 
            index=rna_cluster_names, 
            columns=atac_cluster_names
        )
        
        self.enrichment_matrix = pd.DataFrame(
            enrichment_data, 
            index=rna_cluster_names, 
            columns=atac_cluster_names
        )
        
        logger.info("Overlap enrichment analysis completed")
    
    def generate_heatmaps(self, significance_threshold: float = 0.05) -> None:
        """
        Generate heatmaps for overlap analysis.
        
        Args:
            significance_threshold: P-value threshold for significance marking
        """
        logger.info("Generating heatmaps...")
        
        if self.overlap_matrix is None:
            raise ValueError("Must run calculate_overlap_enrichment() first")
        
        # Create output files
        overlap_heatmap_file = self.output_dir / "rna_atac_overlap_heatmap.pdf"
        enrichment_heatmap_file = self.output_dir / "rna_atac_enrichment_heatmap.pdf"
        
        # Generate first heatmap: Overlap counts with proportions as colors
        self._generate_overlap_heatmap(overlap_heatmap_file, significance_threshold)
        
        # Generate second heatmap: Enrichment analysis
        self._generate_enrichment_heatmap(enrichment_heatmap_file, significance_threshold)
        
        # Save matrices for R processing
        self._save_matrices_for_r()
        
        logger.info(f"Heatmaps saved to: {self.output_dir}")
    
    def _generate_overlap_heatmap(self, output_file: Path, significance_threshold: float) -> None:
        """Generate overlap count heatmap with proportion-based coloring."""
        
        with PdfPages(output_file) as pdf:
            fig, ax = plt.subplots(figsize=(12, 8))
            
            # Use proportion matrix for coloring
            sns.heatmap(
                self.proportion_matrix,
                annot=self.overlap_matrix,  # Annotate with overlap counts
                fmt='d',  # Integer format for annotations
                cmap='YlOrRd',
                cbar_kws={'label': 'Proportion of RNA-seq cluster genes'},
                ax=ax,
                annot_kws={'fontsize': 10}
            )
            
            # Add significance markers
            self._add_significance_markers(ax, significance_threshold)
            
            ax.set_title('RNA-seq vs ATAC-seq Cluster Gene Overlap\n'
                        '(Colors: Proportions, Numbers: Overlap Counts)', 
                        fontsize=14, fontweight='bold')
            ax.set_xlabel('ATAC-seq Clusters', fontsize=12)
            ax.set_ylabel('RNA-seq Clusters', fontsize=12)
            
            plt.tight_layout()
            pdf.savefig(fig, dpi=300, bbox_inches='tight')
            plt.close()
    
    def _generate_enrichment_heatmap(self, output_file: Path, significance_threshold: float) -> None:
        """Generate enrichment fold-change heatmap."""
        
        with PdfPages(output_file) as pdf:
            fig, ax = plt.subplots(figsize=(12, 8))
            
            # Use enrichment matrix for coloring
            sns.heatmap(
                self.enrichment_matrix,
                annot=True,
                fmt='.2f',
                cmap='RdBu_r',
                center=1,  # Center colormap at 1 (no enrichment)
                cbar_kws={'label': 'Enrichment Fold Change'},
                ax=ax,
                annot_kws={'fontsize': 10}
            )
            
            # Add significance markers
            self._add_significance_markers(ax, significance_threshold)
            
            ax.set_title('RNA-seq vs ATAC-seq Cluster Gene Enrichment Analysis\n'
                        '(Fold Change Enrichment with Significance)', 
                        fontsize=14, fontweight='bold')
            ax.set_xlabel('ATAC-seq Clusters', fontsize=12)
            ax.set_ylabel('RNA-seq Clusters', fontsize=12)
            
            plt.tight_layout()
            pdf.savefig(fig, dpi=300, bbox_inches='tight')
            plt.close()
    
    def _add_significance_markers(self, ax, significance_threshold: float) -> None:
        """Add asterisk markers for significant p-values."""
        
        for i in range(len(self.pvalue_matrix.index)):
            for j in range(len(self.pvalue_matrix.columns)):
                pval = self.pvalue_matrix.iloc[i, j]
                if pval < significance_threshold:
                    # Add asterisk marker
                    ax.text(j + 0.8, i + 0.2, '*', 
                           fontsize=16, fontweight='bold', 
                           color='black', ha='center', va='center')
    
    def _save_matrices_for_r(self) -> None:
        """Save matrices in formats suitable for R ComplexHeatmap processing."""
        
        # Save overlap matrix
        overlap_file = self.output_dir / "overlap_matrix_for_R.csv"
        self.overlap_matrix.to_csv(overlap_file)
        
        # Save proportion matrix
        proportion_file = self.output_dir / "proportion_matrix_for_R.csv"
        self.proportion_matrix.to_csv(proportion_file)
        
        # Save p-value matrix
        pvalue_file = self.output_dir / "pvalue_matrix_for_R.csv"
        self.pvalue_matrix.to_csv(pvalue_file)
        
        # Save enrichment matrix
        enrichment_file = self.output_dir / "enrichment_matrix_for_R.csv"
        self.enrichment_matrix.to_csv(enrichment_file)
        
        # Create R script template
        r_script_file = self.output_dir / "complexheatmap_template.R"
        self._create_r_script_template(r_script_file)
        
        logger.info("Matrices saved for R ComplexHeatmap processing")
    
    def _create_r_script_template(self, output_file: Path) -> None:
        """Create R script template for ComplexHeatmap visualization."""
        
        r_script = '''
# R script template for ComplexHeatmap visualization
# Load required libraries
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# Set working directory to the output folder
setwd("{output_dir}")

# Load data matrices
overlap_matrix <- read.csv("overlap_matrix_for_R.csv", row.names = 1)
proportion_matrix <- read.csv("proportion_matrix_for_R.csv", row.names = 1)
pvalue_matrix <- read.csv("pvalue_matrix_for_R.csv", row.names = 1)
enrichment_matrix <- read.csv("enrichment_matrix_for_R.csv", row.names = 1)

# Convert to matrices
overlap_mat <- as.matrix(overlap_matrix)
proportion_mat <- as.matrix(proportion_matrix)
pvalue_mat <- as.matrix(pvalue_matrix)
enrichment_mat <- as.matrix(enrichment_matrix)

# Create significance annotation
significance_mat <- ifelse(pvalue_mat < 0.05, "*", "")

# Define color functions
col_proportion <- colorRamp2(c(0, max(proportion_mat)), c("white", "red"))
col_enrichment <- colorRamp2(c(min(enrichment_mat), 1, max(enrichment_mat)), 
                            c("blue", "white", "red"))

# Create heatmap 1: Overlap counts with proportion coloring
pdf("complexheatmap_overlap.pdf", width = 10, height = 8)
ht1 <- Heatmap(
    proportion_mat,
    name = "Proportion",
    col = col_proportion,
    cell_fun = function(j, i, x, y, width, height, fill) {{
        grid.text(overlap_mat[i, j], x, y, gp = gpar(fontsize = 10))
        if(significance_mat[i, j] == "*") {{
            grid.text("*", x, y + unit(2, "mm"), gp = gpar(fontsize = 14, fontface = "bold"))
        }}
    }},
    column_title = "ATAC-seq Clusters",
    row_title = "RNA-seq Clusters",
    heatmap_legend_param = list(title = "Proportion of\\nRNA-seq cluster genes")
)
draw(ht1)
dev.off()

# Create heatmap 2: Enrichment analysis
pdf("complexheatmap_enrichment.pdf", width = 10, height = 8)
ht2 <- Heatmap(
    enrichment_mat,
    name = "Enrichment",
    col = col_enrichment,
    cell_fun = function(j, i, x, y, width, height, fill) {{
        grid.text(sprintf("%.2f", enrichment_mat[i, j]), x, y, gp = gpar(fontsize = 10))
        if(significance_mat[i, j] == "*") {{
            grid.text("*", x, y + unit(2, "mm"), gp = gpar(fontsize = 14, fontface = "bold"))
        }}
    }},
    column_title = "ATAC-seq Clusters",
    row_title = "RNA-seq Clusters",
    heatmap_legend_param = list(title = "Enrichment\\nFold Change")
)
draw(ht2)
dev.off()

cat("ComplexHeatmap visualizations saved to PDF files\\n")
'''.format(output_dir=str(self.output_dir.absolute()))
        
        with open(output_file, 'w') as f:
            f.write(r_script)
    
    def save_summary_report(self) -> None:
        """Save a summary report of the analysis."""
        
        report_file = self.output_dir / "overlap_analysis_summary.txt"
        
        with open(report_file, 'w') as f:
            f.write("RNA-seq and ATAC-seq Cluster Gene Overlap Analysis Summary\n")
            f.write("=" * 65 + "\n\n")
            
            # Basic statistics
            f.write("Dataset Information:\n")
            f.write(f"RNA-seq clusters: {len(self.rna_clusters)}\n")
            f.write(f"ATAC-seq clusters: {len(self.atac_clusters)}\n")
            f.write(f"Universal gene set size: {len(self.universal_genes)}\n\n")
            
            # Cluster sizes
            f.write("RNA-seq Cluster Sizes:\n")
            for cluster, genes in sorted(self.rna_clusters.items()):
                f.write(f"  {cluster}: {len(genes)} genes\n")
            
            f.write("\nATAC-seq Cluster Sizes:\n")
            for cluster, genes in sorted(self.atac_clusters.items()):
                f.write(f"  {cluster}: {len(genes)} genes\n")
            
            # Significant overlaps
            f.write("\nSignificant Overlaps (p < 0.05):\n")
            significant_pairs = []
            for i, rna_cluster in enumerate(self.overlap_matrix.index):
                for j, atac_cluster in enumerate(self.overlap_matrix.columns):
                    pval = self.pvalue_matrix.iloc[i, j]
                    if pval < 0.05:
                        overlap = self.overlap_matrix.iloc[i, j]
                        enrichment = self.enrichment_matrix.iloc[i, j]
                        significant_pairs.append((rna_cluster, atac_cluster, overlap, enrichment, pval))
            
            if significant_pairs:
                significant_pairs.sort(key=lambda x: x[4])  # Sort by p-value
                for rna, atac, overlap, enrichment, pval in significant_pairs:
                    f.write(f"  {rna} vs {atac}: {overlap} genes, "
                           f"enrichment={enrichment:.2f}, p={pval:.2e}\n")
            else:
                f.write("  No significant overlaps found (p < 0.05)\n")
        
        logger.info(f"Summary report saved to: {report_file}")
    
    def save_overlapping_genes(self) -> None:
        """
        Save overlapping genes between each RNA-seq and ATAC-seq cluster pair to individual text files.
        """
        logger.info("Saving overlapping genes to individual files...")
        
        if not self.rna_clusters or not self.atac_clusters:
            raise ValueError("Both RNA and ATAC clusters must be loaded first")
        
        # Create overlapping genes directory
        overlapping_genes_dir = self.output_dir / "overlapping_genes"
        overlapping_genes_dir.mkdir(exist_ok=True)
        
        # Get sorted cluster names for consistent ordering
        rna_cluster_names = sorted(self.rna_clusters.keys())
        atac_cluster_names = sorted(self.atac_clusters.keys())
        
        total_files_created = 0
        
        for rna_cluster in rna_cluster_names:
            rna_genes = self.rna_clusters[rna_cluster]
            
            for atac_cluster in atac_cluster_names:
                atac_genes = self.atac_clusters[atac_cluster]
                
                # Calculate overlap
                overlap_genes = rna_genes.intersection(atac_genes)
                overlap_count = len(overlap_genes)
                
                # Only create file if there are overlapping genes
                if overlap_count > 0:
                    # Create filename
                    filename = f"{rna_cluster}_vs_{atac_cluster}_overlapping_genes.txt"
                    filepath = overlapping_genes_dir / filename
                    
                    # Get statistics if matrices are available
                    proportion = None
                    pvalue = None
                    enrichment = None
                    
                    if (self.proportion_matrix is not None and 
                        self.pvalue_matrix is not None and 
                        self.enrichment_matrix is not None):
                        proportion = self.proportion_matrix.loc[rna_cluster, atac_cluster]
                        pvalue = self.pvalue_matrix.loc[rna_cluster, atac_cluster]
                        enrichment = self.enrichment_matrix.loc[rna_cluster, atac_cluster]
                    
                    # Write overlapping genes to file
                    with open(filepath, 'w') as f:
                        f.write(f"# Overlapping genes between {rna_cluster} (RNA-seq) and {atac_cluster} (ATAC-seq)\n")
                        f.write(f"# RNA-seq cluster size: {len(rna_genes)} genes\n")
                        f.write(f"# ATAC-seq cluster size: {len(atac_genes)} genes\n")
                        f.write(f"# Overlapping genes: {overlap_count} genes\n")
                        
                        if proportion is not None:
                            f.write(f"# Proportion of RNA-seq cluster: {proportion:.4f}\n")
                        if enrichment is not None:
                            f.write(f"# Enrichment fold change: {enrichment:.4f}\n")
                        if pvalue is not None:
                            f.write(f"# P-value (hypergeometric test): {pvalue:.2e}\n")
                        
                        f.write(f"# Generated on: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                        f.write("#\n")
                        f.write("# Gene list:\n")
                        f.write("GeneID\n")
                        
                        # Write genes in sorted order
                        for gene in sorted(overlap_genes):
                            f.write(f"{gene}\n")
                    
                    total_files_created += 1
                    logger.info(f"Created: {filename} ({overlap_count} genes)")
        
        # Create summary file
        summary_file = overlapping_genes_dir / "overlapping_genes_summary.txt"
        with open(summary_file, 'w') as f:
            f.write("Overlapping Genes Files Summary\n")
            f.write("=" * 40 + "\n\n")
            f.write(f"Total files created: {total_files_created}\n")
            f.write(f"RNA-seq clusters: {len(rna_cluster_names)}\n")
            f.write(f"ATAC-seq clusters: {len(atac_cluster_names)}\n")
            f.write(f"Generated on: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            f.write("Files created:\n")
            for rna_cluster in rna_cluster_names:
                for atac_cluster in atac_cluster_names:
                    overlap_genes = self.rna_clusters[rna_cluster].intersection(
                        self.atac_clusters[atac_cluster]
                    )
                    if len(overlap_genes) > 0:
                        filename = f"{rna_cluster}_vs_{atac_cluster}_overlapping_genes.txt"
                        f.write(f"- {filename}: {len(overlap_genes)} genes\n")
        
        logger.info(f"Saved {total_files_created} overlapping gene files to: {overlapping_genes_dir}")
        logger.info(f"Summary saved to: {summary_file}")


def main():
    """Main function for command-line interface."""
    
    parser = argparse.ArgumentParser(
        description="Analyze RNA-seq and ATAC-seq cluster gene overlaps with enrichment analysis"
    )
    
    parser.add_argument(
        "--rna-clusters", 
        required=True,
        help="Directory containing RNA-seq cluster gene files"
    )
    
    parser.add_argument(
        "--atac-clusters", 
        required=True,
        help="Directory containing ATAC-seq cluster gene files"
    )
    
    parser.add_argument(
        "--output-dir", 
        default="output/cluster_overlap_analysis",
        help="Output directory for results (default: output/cluster_overlap_analysis)"
    )
    
    parser.add_argument(
        "--significance-threshold", 
        type=float, 
        default=0.05,
        help="P-value threshold for significance (default: 0.05)"
    )
    
    args = parser.parse_args()
    
    try:
        # Initialize analyzer
        analyzer = ClusterOverlapAnalyzer(args.output_dir)
        
        # Load data
        analyzer.load_rna_clusters(args.rna_clusters)
        analyzer.load_atac_clusters(args.atac_clusters)
        
        # Perform analysis
        analyzer.calculate_overlap_enrichment()
        
        # Generate visualizations
        analyzer.generate_heatmaps(args.significance_threshold)
        
        # Save overlapping genes to individual files
        analyzer.save_overlapping_genes()
        
        # Save summary
        analyzer.save_summary_report()
        
        logger.info("Analysis completed successfully!")
        
    except Exception as e:
        logger.error(f"Analysis failed: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main() 