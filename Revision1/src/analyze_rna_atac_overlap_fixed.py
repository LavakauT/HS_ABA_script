#!/usr/bin/env python3
"""
RNA-seq and ATAC-seq Cluster Gene Overlap Enrichment Analysis (Fixed Version)

This script analyzes the overlap between RNA-seq cluster genes and ATAC-seq cluster genes,
calculates enrichment statistics using hypergeometric test, and generates heatmaps for visualization.

This is a fixed version that handles the actual file naming convention used in the project:
- RNA clusters: UU.txt, UN.txt, UD.txt, NU.txt, NN.txt, ND.txt, DU.txt, DN.txt, DD.txt
- ATAC clusters: Same naming convention in peak2gene directory

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
from statsmodels.stats.multitest import multipletests

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
    
    def __init__(self, output_dir: str = "output/cluster_overlap_analysis", total_genome_genes: Optional[int] = None):
        """
        Initialize the analyzer.
        
        Args:
            output_dir: Directory to save output files
            total_genome_genes: Total number of genes in the genome (if None, will use ATAC ∪ RNA genes as universe)
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Data storage
        self.rna_clusters: Dict[str, Set[str]] = {}
        self.atac_clusters: Dict[str, Set[str]] = {}
        self.universal_genes: Set[str] = set()
        self.total_genome_genes: Optional[int] = total_genome_genes
        
        # Results storage
        self.overlap_matrix: Optional[pd.DataFrame] = None
        self.proportion_matrix: Optional[pd.DataFrame] = None
        self.pvalue_matrix: Optional[pd.DataFrame] = None
        self.enrichment_matrix: Optional[pd.DataFrame] = None
        self.fdr_matrix: Optional[pd.DataFrame] = None  # Added FDR matrix
        
        # Detailed calculation storage for debugging
        self.calculation_details: Dict[Tuple[str, str], Dict] = {}
    
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
        
        # Find all cluster gene files with the actual naming convention
        # Expected files: UU.txt, UN.txt, UD.txt, NU.txt, NN.txt, ND.txt, DU.txt, DN.txt, DD.txt
        expected_clusters = ['UU', 'UN', 'UD', 'NU', 'NN', 'ND', 'DU', 'DN', 'DD']
        
        for cluster_name in expected_clusters:
            file_path = cluster_dir / f"{cluster_name}.txt"
            if file_path.exists():
                logger.info(f"Loading RNA cluster: {cluster_name}")
                genes = self._load_gene_file(file_path, skip_header=True)
                self.rna_clusters[cluster_name] = genes
                self.universal_genes.update(genes)
            else:
                logger.warning(f"RNA cluster file not found: {file_path}")
        
        if not self.rna_clusters:
            raise FileNotFoundError(f"No cluster gene files found in: {rna_cluster_dir}")
        
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
        
        # Find all cluster gene files with the actual naming convention
        # Expected files: UU.txt, UN.txt, UD.txt, NU.txt, ND.txt, DU.txt, DN.txt, DD.txt
        expected_clusters = ['UU', 'UN', 'UD', 'NU', 'NN', 'ND', 'DU', 'DN', 'DD']
        
        for cluster_name in expected_clusters:
            file_path = cluster_dir / f"{cluster_name}.txt"
            if file_path.exists():
                logger.info(f"Loading ATAC cluster: {cluster_name}")
                genes = self._load_gene_file(file_path, skip_header=True)
                self.atac_clusters[cluster_name] = genes
                self.universal_genes.update(genes)
            else:
                logger.warning(f"ATAC cluster file not found: {file_path}")
        
        if not self.atac_clusters:
            raise FileNotFoundError(f"No cluster gene files found in: {atac_cluster_dir}")
        
        logger.info(f"Loaded {len(self.atac_clusters)} ATAC-seq clusters")
        for cluster, genes in self.atac_clusters.items():
            logger.info(f"  {cluster}: {len(genes)} genes")
    
    def _load_gene_file(self, file_path: Path, skip_header: bool = True) -> Set[str]:
        """
        Load genes from a text file.
        
        Args:
            file_path: Path to the gene file
            skip_header: Whether to skip the first line (header)
            
        Returns:
            Set of gene identifiers
        """
        genes = set()
        
        try:
            with open(file_path, 'r') as f:
                lines = f.readlines()
                
            if skip_header and lines:
                lines = lines[1:]  # Skip header line
                
            for line in lines:
                line = line.strip()
                if line and not line.startswith('#'):
                    genes.add(line)
                    
        except Exception as e:
            logger.error(f"Error reading file {file_path}: {e}")
            raise
            
        return genes
    
    def _estimate_total_genome_size(self) -> int:
        """
        Estimate the total genome size based on the loaded data.
        
        Returns:
            Estimated total number of genes in the genome
        """
        if self.total_genome_genes is not None:
            return self.total_genome_genes
        
        # FIXED: Use ATAC ∪ RNA genes as universe (union approach)
        # This represents all genes that are present in either ATAC or RNA datasets
        atac_all_genes = set()
        rna_all_genes = set()
        
        # Collect all ATAC genes
        for genes in self.atac_clusters.values():
            atac_all_genes.update(genes)
        
        # Collect all RNA genes
        for genes in self.rna_clusters.values():
            rna_all_genes.update(genes)
        
        # Universe = ATAC ∪ RNA genes (union)
        universe_genes = atac_all_genes | rna_all_genes
        universe_size = len(universe_genes)
        
        logger.info(f"ATAC total genes: {len(atac_all_genes)}")
        logger.info(f"RNA total genes: {len(rna_all_genes)}")
        logger.info(f"Universe (ATAC ∪ RNA): {universe_size} genes")
        
        if universe_size == 0:
            raise ValueError("No genes found in ATAC or RNA datasets. Check data consistency.")
        
        return universe_size
    
    def calculate_overlaps(self) -> None:
        """
        Calculate overlaps between RNA-seq and ATAC-seq clusters.
        """
        logger.info("Calculating cluster overlaps...")
        
        # Get the correct total genome size
        N_total = self._estimate_total_genome_size()
        logger.info(f"Using total genome size: {N_total}")
        
        # Initialize matrices
        rna_clusters = sorted(self.rna_clusters.keys())
        atac_clusters = sorted(self.atac_clusters.keys())
        
        self.overlap_matrix = pd.DataFrame(0, index=rna_clusters, columns=atac_clusters)
        self.proportion_matrix = pd.DataFrame(0.0, index=rna_clusters, columns=atac_clusters)
        self.pvalue_matrix = pd.DataFrame(1.0, index=rna_clusters, columns=atac_clusters)
        self.enrichment_matrix = pd.DataFrame(0.0, index=rna_clusters, columns=atac_clusters)
        
        # Calculate overlaps and enrichment
        for rna_cluster in rna_clusters:
            rna_genes = self.rna_clusters[rna_cluster]
            
            for atac_cluster in atac_clusters:
                atac_genes = self.atac_clusters[atac_cluster]
                
                # Calculate overlap
                overlap = len(rna_genes & atac_genes)
                self.overlap_matrix.loc[rna_cluster, atac_cluster] = overlap
                
                # Calculate proportion
                if len(rna_genes) > 0:
                    proportion = overlap / len(rna_genes)
                    self.proportion_matrix.loc[rna_cluster, atac_cluster] = proportion
                
                # Calculate hypergeometric test p-value and enrichment
                if len(rna_genes) > 0 and len(atac_genes) > 0:
                    # FIXED: Use correct total genome size for hypergeometric test
                    # N = total genes in genome (not m + n)
                    N = N_total
                    m = len(atac_genes)            # ATAC cluster size (white balls)
                    n = N - m                      # Non-ATAC genes (black balls)
                    k = len(rna_genes)             # RNA cluster size (balls drawn)
                    q = overlap                    # Observed overlap (white balls drawn)
                    
                    # Calculate expected overlap under null hypothesis
                    expected = (m * k) / N
                    
                    if expected > 0:
                        # FIXED: Calculate p-value using R phyper equivalent
                        # R: phyper(q, m, n, k, lower.tail=FALSE, log.p=FALSE)
                        # This gives P(X > q) = 1 - P(X <= q)
                        if q > 0:
                            # For P(X > q), use 1 - P(X <= q)
                            p_value = 1 - hypergeom.cdf(q, N, m, k)
                        else:
                            p_value = 1.0  # No overlap observed
                        
                        self.pvalue_matrix.loc[rna_cluster, atac_cluster] = p_value
                        
                        # FIXED: Calculate enrichment using log2 fold change
                        # This should properly reflect the biological significance
                        if q > 0 and expected > 0:
                            # Calculate fold change
                            fold_change = q / expected
                            
                            # Calculate log2 enrichment
                            if fold_change > 0:
                                enrichment = np.log2(fold_change)
                            else:
                                enrichment = -np.inf  # No overlap, very low enrichment
                            
                            # Additional validation: if p-value is significant, enrichment should reflect it
                            if p_value < 0.05 and q > expected:
                                # For significant positive enrichment, ensure enrichment is positive
                                if enrichment <= 0:
                                    logger.warning(f"Significant p-value ({p_value:.6f}) but low enrichment ({enrichment:.4f}) for {rna_cluster} vs {atac_cluster}")
                                    logger.warning(f"q={q}, expected={expected:.4f}, fold_change={fold_change:.4f}")
                        else:
                            enrichment = -np.inf  # No overlap or invalid expected value
                        
                        self.enrichment_matrix.loc[rna_cluster, atac_cluster] = enrichment
                        
                        # Store detailed calculation info for debugging
                        self._store_calculation_details(
                            rna_cluster, atac_cluster, 
                            N, m, n, k, q, expected, p_value, enrichment
                        )
        
        logger.info("Overlap calculation completed")
        
        # Apply FDR correction
        self._apply_fdr_correction()
    
    def _apply_fdr_correction(self) -> None:
        """
        Apply FDR correction using Benjamini-Hochberg method.
        """
        logger.info("Applying FDR correction...")
        
        # Flatten p-values for correction
        pvalues_flat = self.pvalue_matrix.values.flatten()
        pvalues_flat = pvalues_flat[~np.isnan(pvalues_flat)]  # Remove NaN values
        
        if len(pvalues_flat) > 0:
            # Apply FDR correction
            rejected, pvals_corrected, alphacSidak, alphacBonf = multipletests(
                pvalues_flat, 
                alpha=0.05, 
                method='fdr_bh'
            )
            
            # Reshape corrected p-values back to matrix
            self.fdr_matrix = self.pvalue_matrix.copy()
            mask = ~np.isnan(self.pvalue_matrix.values)
            self.fdr_matrix.values[mask] = pvals_corrected
            
            logger.info(f"FDR correction applied to {len(pvalues_flat)} p-values")
            logger.info(f"Significant after FDR correction (alpha=0.05): {np.sum(rejected)}")
        else:
            logger.warning("No valid p-values found for FDR correction")
            self.fdr_matrix = self.pvalue_matrix.copy()
    
    def _store_calculation_details(self, rna_cluster: str, atac_cluster: str, 
                                 N: int, m: int, n: int, k: int, q: int, 
                                 expected: float, p_value: float, enrichment: float) -> None:
        """
        Store detailed calculation parameters for debugging and verification.
        
        Args:
            rna_cluster: RNA cluster name
            atac_cluster: ATAC cluster name
            N: Total number of genes
            m: Number of ATAC genes (white balls)
            n: Number of non-ATAC genes (black balls)
            k: Number of RNA genes (balls drawn)
            q: Observed overlap (white balls drawn)
            expected: Expected overlap under null hypothesis
            p_value: Calculated p-value
            enrichment: Calculated enrichment
        """
        key = (rna_cluster, atac_cluster)
        self.calculation_details[key] = {
            'N': N,                    # Total genes
            'm': m,                    # ATAC genes (white balls)
            'n': n,                    # Non-ATAC genes (black balls)
            'k': k,                    # RNA genes (balls drawn)
            'q': q,                    # Observed overlap
            'expected': expected,      # Expected overlap
            'p_value': p_value,       # P-value
            'enrichment': enrichment,  # Log2 enrichment
            'r_phyper_equivalent': f"phyper({q}, {m}, {n}, {k}, lower.tail=FALSE, log.p=FALSE)" if q > 0 else "N/A"
        }
    
    def generate_heatmaps(self) -> None:
        """
        Generate heatmap visualizations for the analysis results.
        """
        logger.info("Generating heatmap visualizations...")
        
        # Create a PDF with multiple heatmaps
        pdf_path = self.output_dir / "cluster_overlap_heatmaps.pdf"
        
        with PdfPages(pdf_path) as pdf:
            # 1. Overlap count heatmap
            self._create_heatmap(
                self.overlap_matrix, 
                "Gene Overlap Counts Between RNA-seq and ATAC-seq Clusters",
                "Overlap Count",
                pdf
            )
            
            # 2. Proportion heatmap
            self._create_heatmap(
                self.proportion_matrix, 
                "Proportion of RNA-seq Genes Overlapping with ATAC-seq Clusters",
                "Proportion",
                pdf
            )
            
            # 3. P-value heatmap (log-transformed)
            pvalue_log = -np.log10(self.pvalue_matrix)
            pvalue_log = pvalue_log.replace([np.inf, -np.inf], 0)
            self._create_heatmap(
                pvalue_log, 
                "Statistical Significance of Overlaps (-log10 P-value)",
                "-log10 P-value",
                pdf
            )
            
            # 4. FDR-corrected P-value heatmap (log-transformed)
            if self.fdr_matrix is not None:
                fdr_log = -np.log10(self.fdr_matrix)
                fdr_log = fdr_log.replace([np.inf, -np.inf], 0)
                self._create_heatmap(
                    fdr_log, 
                    "FDR-Corrected Statistical Significance (-log10 FDR)",
                    "-log10 FDR",
                    pdf
                )
            
            # 5. Enrichment heatmap (log2 fold change)
            # Handle -inf values for visualization
            enrichment_clean = self.enrichment_matrix.replace([-np.inf, np.inf], np.nan)
            enrichment_clean = enrichment_clean.fillna(-10)  # Replace -inf with -10 for visualization
            
            self._create_heatmap(
                enrichment_clean, 
                "Enrichment of Overlaps (log2 Fold Change)",
                "log2 Fold Change",
                pdf
            )
        
        logger.info(f"Heatmaps saved to: {pdf_path}")
    
    def _create_heatmap(self, data: pd.DataFrame, title: str, cbar_label: str, pdf: PdfPages) -> None:
        """
        Create a single heatmap and add it to the PDF.
        
        Args:
            data: Data matrix for the heatmap
            title: Title for the heatmap
            cbar_label: Label for the colorbar
            pdf: PDF object to add the plot to
        """
        plt.figure(figsize=(10, 8))
        
        # Create heatmap
        sns.heatmap(
            data, 
            annot=True, 
            fmt='.3f', 
            cmap='RdBu_r', 
            center=0 if 'Enrichment' in title else None,
            cbar_kws={'label': cbar_label}
        )
        
        plt.title(title, fontsize=14, fontweight='bold')
        plt.xlabel('ATAC-seq Clusters', fontsize=12)
        plt.ylabel('RNA-seq Clusters', fontsize=12)
        plt.tight_layout()
        
        # Add to PDF
        pdf.savefig()
        plt.close()
    
    def save_results(self) -> None:
        """
        Save analysis results to CSV files.
        """
        logger.info("Saving analysis results...")
        
        # Save matrices
        self.overlap_matrix.to_csv(self.output_dir / "overlap_counts.csv")
        self.proportion_matrix.to_csv(self.output_dir / "overlap_proportions.csv")
        self.pvalue_matrix.to_csv(self.output_dir / "overlap_pvalues.csv")
        
        # Save FDR-corrected p-values
        if self.fdr_matrix is not None:
            self.fdr_matrix.to_csv(self.output_dir / "overlap_fdr_corrected.csv")
        
        # Handle -inf values in enrichment matrix for CSV output
        enrichment_clean = self.enrichment_matrix.replace([-np.inf, np.inf], "NA")
        enrichment_clean.to_csv(self.output_dir / "overlap_enrichment.csv")
        
        # Create summary statistics
        summary_data = []
        for rna_cluster in self.rna_clusters.keys():
            for atac_cluster in self.atac_clusters.keys():
                enrichment_val = self.enrichment_matrix.loc[rna_cluster, atac_cluster]
                # Handle -inf values in summary
                if np.isinf(enrichment_val):
                    enrichment_val = "NA" if enrichment_val < 0 else "Inf"
                
                fdr_val = self.fdr_matrix.loc[rna_cluster, atac_cluster] if self.fdr_matrix is not None else "NA"
                
                summary_data.append({
                    'RNA_Cluster': rna_cluster,
                    'ATAC_Cluster': atac_cluster,
                    'Overlap_Count': self.overlap_matrix.loc[rna_cluster, atac_cluster],
                    'Proportion': self.proportion_matrix.loc[rna_cluster, atac_cluster],
                    'P_value': self.pvalue_matrix.loc[rna_cluster, atac_cluster],
                    'FDR_Corrected': fdr_val,
                    'neg_log10_P_value': -np.log10(self.pvalue_matrix.loc[rna_cluster, atac_cluster]) if self.pvalue_matrix.loc[rna_cluster, atac_cluster] > 0 else np.inf,
                    'neg_log10_FDR': -np.log10(fdr_val) if isinstance(fdr_val, (int, float)) and fdr_val > 0 else np.inf
                })
        
        summary_df = pd.DataFrame(summary_data)
        summary_df.to_csv(self.output_dir / "overlap_summary.csv", index=False)
        
        logger.info(f"Results saved to: {self.output_dir}")
    
    def save_detailed_calculations(self) -> None:
        """
        Save detailed calculation parameters for verification.
        """
        logger.info("Saving detailed calculation parameters...")
        
        # Create detailed calculation report
        detailed_data = []
        for (rna_cluster, atac_cluster), details in self.calculation_details.items():
            fdr_val = self.fdr_matrix.loc[rna_cluster, atac_cluster] if self.fdr_matrix is not None else "NA"
            
            row = {
                'RNA_Cluster': rna_cluster,
                'ATAC_Cluster': atac_cluster,
                'N_Total_Genes': details['N'],
                'm_ATAC_Genes': details['m'],
                'n_NonATAC_Genes': details['n'],
                'k_RNA_Genes': details['k'],
                'q_Observed_Overlap': details['q'],
                'P_value': details['p_value'],
                'FDR_Corrected': fdr_val,
                'neg_log10_P_value': -np.log10(details['p_value']) if details['p_value'] > 0 else np.inf,
                'neg_log10_FDR': -np.log10(fdr_val) if isinstance(fdr_val, (int, float)) and fdr_val > 0 else np.inf,
                'R_phyper_Equivalent': details['r_phyper_equivalent']
            }
            detailed_data.append(row)
        
        # Save as CSV
        detailed_df = pd.DataFrame(detailed_data)
        detailed_df.to_csv(self.output_dir / "detailed_calculation_parameters.csv", index=False)
        
        # Save as text report
        report_path = self.output_dir / "detailed_calculation_report.txt"
        with open(report_path, 'w') as f:
            f.write("RNA-ATAC Cluster Overlap Analysis - Detailed Calculation Report\n")
            f.write("=" * 80 + "\n\n")
            
            f.write("Hypergeometric Test Parameters (R phyper format):\n")
            f.write("phyper(q, m, n, k, lower.tail = FALSE, log.p = FALSE)\n\n")
            
            f.write("Where:\n")
            f.write("- q = observed overlap (white balls drawn)\n")
            f.write("- m = ATAC cluster size (white balls in urn)\n")
            f.write("- n = non-ATAC genes (black balls in urn)\n")
            f.write("- k = RNA cluster size (balls drawn)\n")
            f.write("- N = total genes (m + n)\n\n")
            
            f.write("FDR Correction: Benjamini-Hochberg method (alpha = 0.05)\n\n")
            
            f.write("Detailed Results:\n")
            f.write("-" * 80 + "\n")
            
            for (rna_cluster, atac_cluster), details in self.calculation_details.items():
                fdr_val = self.fdr_matrix.loc[rna_cluster, atac_cluster] if self.fdr_matrix is not None else "NA"
                
                f.write(f"\nRNA Cluster: {rna_cluster} vs ATAC Cluster: {atac_cluster}\n")
                f.write(f"  N (Total genes): {details['N']}\n")
                f.write(f"  m (ATAC genes): {details['m']}\n")
                f.write(f"  n (Non-ATAC genes): {details['n']}\n")
                f.write(f"  k (RNA genes): {details['k']}\n")
                f.write(f"  q (Observed overlap): {details['q']}\n")
                f.write(f"  P-value: {details['p_value']:.6f}\n")
                f.write(f"  FDR-corrected: {fdr_val}\n")
                f.write(f"  -log10(P-value): {-np.log10(details['p_value']):.4f}\n")
                f.write(f"  -log10(FDR): {-np.log10(fdr_val):.4f}\n" if isinstance(fdr_val, (int, float)) and fdr_val > 0 else f"  -log10(FDR): N/A\n")
                f.write(f"  R phyper equivalent: {details['r_phyper_equivalent']}\n")
        
        logger.info(f"Detailed calculations saved to: {self.output_dir}")
    
    def run_analysis(self, rna_cluster_dir: str, atac_cluster_dir: str) -> None:
        """
        Run the complete analysis pipeline.
        
        Args:
            rna_cluster_dir: Directory containing RNA-seq cluster gene files
            atac_cluster_dir: Directory containing ATAC-seq cluster gene files
        """
        try:
            # Load data
            self.load_rna_clusters(rna_cluster_dir)
            self.load_atac_clusters(atac_cluster_dir)
            
            # Calculate overlaps
            self.calculate_overlaps()
            
            # Generate visualizations
            self.generate_heatmaps()
            
            # Save results
            self.save_results()
            
            # Save detailed calculations
            self.save_detailed_calculations()
            
            logger.info("Analysis completed successfully!")
            
        except Exception as e:
            logger.error(f"Analysis failed: {e}")
            raise


def main():
    """
    Main function to run the cluster overlap analysis.
    """
    parser = argparse.ArgumentParser(
        description="Analyze overlap between RNA-seq and ATAC-seq cluster genes"
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
        "--total-genome-genes",
        type=int,
        help="Total number of genes in the genome (recommended for accurate statistics)"
    )
    
    args = parser.parse_args()
    
    # Create analyzer and run analysis
    analyzer = ClusterOverlapAnalyzer(
        args.output_dir, 
        total_genome_genes=args.total_genome_genes
    )
    analyzer.run_analysis(args.rna_clusters, args.atac_clusters)


if __name__ == "__main__":
    main()
