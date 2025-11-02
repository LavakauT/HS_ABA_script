#!/usr/bin/env python3
"""
Hierarchical Clustering Analysis for VST Transformed Counts Data

This script performs hierarchical clustering analysis on VST transformed counts data
filtered by DEG/DAR files. It includes z-score normalization, bootstrap analysis
for cluster stability, and outputs comprehensive results with heatmaps.

Date: 2025 (HS_ABA project)
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Set
import logging
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import pdist
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score
import warnings
warnings.filterwarnings('ignore')

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class HierarchicalClusteringAnalyzer:
    """
    A comprehensive hierarchical clustering analyzer for VST transformed counts data.
    
    This class handles the complete pipeline from data loading to clustering
    analysis and visualization, including bootstrap analysis for cluster stability.
    """
    
    def __init__(self, vst_data_path: str, deg_dar_folder_path: str, output_dir: str):
        """
        Initialize the hierarchical clustering analyzer.
        
        Args:
            vst_data_path: Path to the VST transformed counts data file
            deg_dar_folder_path: Path to the folder containing DEG/DAR files
            output_dir: Directory to save output files
        """
        self.vst_data_path = Path(vst_data_path)
        self.deg_dar_folder_path = Path(deg_dar_folder_path)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize data containers
        self.vst_data = None
        self.deg_dar_genes = set()
        self.filtered_data = None
        self.z_score_data = None
        self.cluster_labels = None
        self.linkage_matrix = None
        
        # Detect if this is DAR analysis based on folder name
        self.is_dar_analysis = self._detect_dar_analysis()
        self.term_for_features = "peaks" if self.is_dar_analysis else "genes"
        
        logger.info(f"Initialized analyzer with VST data: {self.vst_data_path}")
        logger.info(f"DEG/DAR folder: {self.deg_dar_folder_path}")
        logger.info(f"Output directory: {self.output_dir}")
        logger.info(f"Analysis type: {'DAR' if self.is_dar_analysis else 'DEG'}")
        logger.info(f"Feature terminology: {self.term_for_features}")
    
    def _detect_dar_analysis(self) -> bool:
        """
        Detect if this is a DAR analysis based on folder name and content.
        
        Returns:
            True if DAR analysis, False if DEG analysis
        """
        # Check folder name for DAR indicators
        folder_name = self.deg_dar_folder_path.name.lower()
        if 'dar' in folder_name or 'peak' in folder_name or 'atac' in folder_name:
            return True
        
        # Check if folder contains DAR-related files
        try:
            for file_path in self.deg_dar_folder_path.rglob("*.txt"):
                file_name = file_path.name.lower()
                if 'dar' in file_name or 'peak' in file_name:
                    return True
        except Exception:
            pass
        
        return False
    
    def load_vst_data(self) -> pd.DataFrame:
        """
        Load VST transformed counts data from file.
        
        Returns:
            DataFrame containing the VST transformed counts data
        """
        logger.info("Loading VST transformed counts data...")
        
        try:
            # Load the data, assuming tab-separated format
            self.vst_data = pd.read_csv(self.vst_data_path, sep='\t', index_col=0)
            logger.info(f"Loaded VST data with shape: {self.vst_data.shape}")
            logger.info(f"Sample columns: {list(self.vst_data.columns[:5])}")
            logger.info(f"Sample genes: {list(self.vst_data.index[:5])}")
            
            return self.vst_data
            
        except Exception as e:
            logger.error(f"Error loading VST data: {e}")
            raise
    
    def load_deg_dar_files(self) -> Set[str]:
        """
        Load all DEG/DAR files from the specified folder and extract gene names.
        
        Returns:
            Set of unique gene names from all DEG/DAR files
        """
        logger.info("Loading DEG/DAR files...")
        
        deg_dar_genes = set()
        
        # Walk through all subdirectories
        for root, dirs, files in os.walk(self.deg_dar_folder_path):
            for file in files:
                if file.endswith('.txt'):
                    file_path = Path(root) / file
                    logger.info(f"Processing file: {file_path}")
                    
                    try:
                        # Read the file and extract gene names
                        with open(file_path, 'r') as f:
                            lines = f.readlines()
                        
                        # Skip header if present (first line contains 'x')
                        if lines and lines[0].strip() == 'x':
                            lines = lines[1:]
                        
                        # Extract gene names
                        file_genes = {line.strip() for line in lines if line.strip()}
                        deg_dar_genes.update(file_genes)
                        
                        logger.info(f"Added {len(file_genes)} genes from {file}")
                        
                    except Exception as e:
                        logger.warning(f"Error processing file {file_path}: {e}")
                        continue
        
        self.deg_dar_genes = deg_dar_genes
        logger.info(f"Total unique genes from DEG/DAR files: {len(self.deg_dar_genes)}")
        
        return self.deg_dar_genes
    
    def filter_data_by_deg_dar(self) -> pd.DataFrame:
        """
        Filter VST data to include only genes present in DEG/DAR files.
        
        Returns:
            Filtered DataFrame containing only DEG/DAR genes
        """
        logger.info("Filtering VST data by DEG/DAR genes...")
        
        if self.vst_data is None:
            raise ValueError("VST data not loaded. Call load_vst_data() first.")
        
        if not self.deg_dar_genes:
            raise ValueError("DEG/DAR genes not loaded. Call load_deg_dar_files() first.")
        
        # Find intersection of genes
        available_genes = set(self.vst_data.index)
        common_genes = self.deg_dar_genes.intersection(available_genes)
        
        logger.info(f"{self.term_for_features.capitalize()} in VST data: {len(available_genes)}")
        logger.info(f"{self.term_for_features.capitalize()} in DEG/DAR files: {len(self.deg_dar_genes)}")
        logger.info(f"Common {self.term_for_features}: {len(common_genes)}")
        
        if len(common_genes) == 0:
            raise ValueError(f"No common {self.term_for_features} found between VST data and DEG/DAR files.")
        
        # Filter the data
        self.filtered_data = self.vst_data.loc[list(common_genes)]
        logger.info(f"Filtered data shape: {self.filtered_data.shape}")
        
        return self.filtered_data
    
    def perform_z_score_normalization(self) -> pd.DataFrame:
        """
        Perform row-wise z-score normalization on the filtered data.
        
        Returns:
            Z-score normalized DataFrame
        """
        logger.info("Performing row-wise z-score normalization...")
        
        if self.filtered_data is None:
            raise ValueError("Filtered data not available. Call filter_data_by_deg_dar() first.")
        
        # Transpose for row-wise normalization (genes as features)
        data_for_scaling = self.filtered_data.T
        
        # Perform z-score normalization
        scaler = StandardScaler()
        z_score_scaled = scaler.fit_transform(data_for_scaling)
        
        # Convert back to DataFrame with original structure
        self.z_score_data = pd.DataFrame(
            z_score_scaled.T,
            index=self.filtered_data.index,
            columns=self.filtered_data.columns
        )
        
        logger.info(f"Z-score data shape: {self.z_score_data.shape}")
        logger.info(f"Z-score statistics - Mean: {self.z_score_data.values.mean():.4f}, "
                   f"Std: {self.z_score_data.values.std():.4f}")
        
        return self.z_score_data
    
    def perform_hierarchical_clustering(self, method: str = 'ward', 
                                     metric: str = 'euclidean',
                                     n_clusters: Optional[int] = None,
                                     distance_threshold: Optional[float] = None) -> np.ndarray:
        """
        Perform hierarchical clustering on the z-score normalized data.
        
        Args:
            method: Linkage method ('ward', 'complete', 'average', 'single')
            metric: Distance metric for linkage
            n_clusters: Number of clusters to form (if None, use distance_threshold)
            distance_threshold: Distance threshold for clustering (if None, use n_clusters)
        
        Returns:
            Array of cluster labels
        """
        logger.info(f"Performing hierarchical clustering with method: {method}, metric: {metric}")
        
        if self.z_score_data is None:
            raise ValueError("Z-score data not available. Call perform_z_score_normalization() first.")
        
        # Prepare data for clustering
        data_for_clustering = self.z_score_data.values
        
        # Calculate linkage matrix
        if method == 'ward' and metric != 'euclidean':
            logger.warning("Ward method requires euclidean distance. Using euclidean metric.")
            metric = 'euclidean'
        
        self.linkage_matrix = linkage(data_for_clustering, method=method, metric=metric)
        
        # Determine number of clusters
        if n_clusters is not None:
            self.cluster_labels = fcluster(self.linkage_matrix, n_clusters, criterion='maxclust')
        elif distance_threshold is not None:
            self.cluster_labels = fcluster(self.linkage_matrix, distance_threshold, criterion='distance')
        else:
            # Auto-determine number of clusters using silhouette analysis
            n_clusters = self._auto_determine_clusters(data_for_clustering)
            self.cluster_labels = fcluster(self.linkage_matrix, n_clusters, criterion='maxclust')
        
        logger.info(f"Clustering completed. Number of clusters: {len(set(self.cluster_labels))}")
        logger.info(f"Cluster sizes: {np.bincount(self.cluster_labels)}")
        
        return self.cluster_labels
    
    def _auto_determine_clusters(self, data: np.ndarray, max_clusters: int = 10) -> int:
        """
        Automatically determine optimal number of clusters using silhouette analysis.
        
        Args:
            data: Data matrix for clustering
            max_clusters: Maximum number of clusters to test
        
        Returns:
            Optimal number of clusters
        """
        logger.info("Auto-determining optimal number of clusters...")
        
        best_score = -1
        best_n_clusters = 2
        
        for n_clusters in range(2, min(max_clusters + 1, len(data) // 2)):
            try:
                labels = fcluster(self.linkage_matrix, n_clusters, criterion='maxclust')
                if len(set(labels)) > 1:  # Ensure we have at least 2 clusters
                    score = silhouette_score(data, labels)
                    if score > best_score:
                        best_score = score
                        best_n_clusters = n_clusters
            except Exception as e:
                logger.warning(f"Error calculating silhouette score for {n_clusters} clusters: {e}")
                continue
        
        logger.info(f"Optimal number of clusters: {best_n_clusters} (silhouette score: {best_score:.4f})")
        return best_n_clusters
    
    def bootstrap_cluster_stability(self, n_bootstrap: int = 100, 
                                  sample_fraction: float = 0.8) -> Dict:
        """
        Perform bootstrap analysis to assess cluster stability.
        
        Args:
            n_bootstrap: Number of bootstrap iterations
            sample_fraction: Fraction of data to sample in each bootstrap
        
        Returns:
            Dictionary containing bootstrap stability metrics
        """
        logger.info(f"Performing bootstrap analysis with {n_bootstrap} iterations...")
        
        if self.z_score_data is None:
            raise ValueError("Z-score data not available.")
        
        data = self.z_score_data.values
        n_samples = len(data)
        sample_size = int(n_samples * sample_fraction)
        
        cluster_agreements = []
        
        for i in range(n_bootstrap):
            # Sample data
            indices = np.random.choice(n_samples, sample_size, replace=False)
            sample_data = data[indices]
            
            # Perform clustering on sample
            sample_linkage = linkage(sample_data, method='ward', metric='euclidean')
            sample_labels = fcluster(sample_linkage, len(set(self.cluster_labels)), criterion='maxclust')
            
            # Calculate agreement with original clustering
            original_sample_labels = self.cluster_labels[indices]
            
            # Simple agreement metric (can be improved with more sophisticated measures)
            agreement = np.mean(sample_labels == original_sample_labels)
            cluster_agreements.append(agreement)
            
            if (i + 1) % 20 == 0:
                logger.info(f"Bootstrap iteration {i + 1}/{n_bootstrap}")
        
        stability_metrics = {
            'mean_agreement': np.mean(cluster_agreements),
            'std_agreement': np.std(cluster_agreements),
            'min_agreement': np.min(cluster_agreements),
            'max_agreement': np.max(cluster_agreements),
            'agreement_scores': cluster_agreements
        }
        
        logger.info(f"Bootstrap stability - Mean agreement: {stability_metrics['mean_agreement']:.4f}")
        logger.info(f"Standard deviation: {stability_metrics['std_agreement']:.4f}")
        
        return stability_metrics
    
    def save_clustering_results(self) -> None:
        """
        Save comprehensive clustering results to files.
        """
        logger.info("Saving clustering results...")
        
        if self.cluster_labels is None:
            raise ValueError("Clustering not performed. Call perform_hierarchical_clustering() first.")
        
        # Create results DataFrame
        results_df = pd.DataFrame({
            'Feature': self.z_score_data.index,
            'Cluster': self.cluster_labels
        })
        
        # Save individual cluster results
        cluster_dir = self.output_dir / "cluster_results"
        cluster_dir.mkdir(exist_ok=True)
        
        unique_clusters = sorted(set(self.cluster_labels))
        cluster_summary = []
        
        for cluster_id in unique_clusters:
            cluster_features = results_df[results_df['Cluster'] == cluster_id]['Feature'].tolist()
            
            # Save individual cluster file
            cluster_file = cluster_dir / f"cluster_{cluster_id}_{self.term_for_features}.txt"
            with open(cluster_file, 'w') as f:
                f.write(f"{self.term_for_features.capitalize()}\n")
                for feature in cluster_features:
                    f.write(f"{feature}\n")
            
            cluster_summary.append({
                'Cluster_ID': cluster_id,
                'Size': len(cluster_features),
                'Features': cluster_features
            })
        
        # Save comprehensive results
        results_file = self.output_dir / "clustering_results_comprehensive.txt"
        with open(results_file, 'w') as f:
            f.write(f"{self.term_for_features.capitalize()}\tCluster\n")
            for _, row in results_df.iterrows():
                f.write(f"{row['Feature']}\t{row['Cluster']}\n")
        
        # Save cluster summary
        summary_file = self.output_dir / "cluster_summary.txt"
        with open(summary_file, 'w') as f:
            f.write(f"Cluster_ID\tSize\t{self.term_for_features.capitalize()}\n")
            for cluster_info in cluster_summary:
                features_str = ",".join(cluster_info['Features'])
                f.write(f"{cluster_info['Cluster_ID']}\t{cluster_info['Size']}\t{features_str}\n")
        
        # Save z-score data for further analysis
        zscore_file = self.output_dir / "zscore_normalized_data.txt"
        self.z_score_data.to_csv(zscore_file, sep='\t')
        
        logger.info(f"Results saved to {self.output_dir}")
        logger.info(f"Individual cluster files saved to {cluster_dir}")
    
    def create_heatmap(self, figsize: Tuple[int, int] = (12, 10),
                      cmap: str = 'RdBu_r',
                      save_pdf: bool = True,
                      save_r_compatible: bool = True) -> None:
        """
        Create and save heatmap visualization of the clustering results.
        
        Args:
            figsize: Figure size for the heatmap
            cmap: Colormap for the heatmap
            save_pdf: Whether to save as PDF
            save_r_compatible: Whether to save data compatible with R ComplexHeatmap
        """
        logger.info("Creating heatmap visualization...")
        
        if self.z_score_data is None or self.cluster_labels is None:
            raise ValueError("Z-score data or cluster labels not available.")
        
        # Sort data by cluster labels
        cluster_order = np.argsort(self.cluster_labels)
        sorted_data = self.z_score_data.iloc[cluster_order]
        sorted_labels = self.cluster_labels[cluster_order]
        
        # Create figure
        fig, ax = plt.subplots(figsize=figsize)
        
        # Create heatmap
        sns.heatmap(sorted_data, 
                   cmap=cmap, 
                   center=0,
                   xticklabels=True,
                   yticklabels=False,  # Don't show gene names for clarity
                   cbar_kws={'label': 'Z-score'},
                   ax=ax)
        
        # Add cluster boundaries
        unique_clusters = np.unique(sorted_labels)
        cluster_boundaries = []
        current_pos = 0
        
        for cluster_id in unique_clusters:
            cluster_size = np.sum(sorted_labels == cluster_id)
            cluster_boundaries.append(current_pos + cluster_size)
            current_pos += cluster_size
        
        # Draw cluster boundaries
        for boundary in cluster_boundaries[:-1]:
            ax.axhline(y=boundary, color='black', linewidth=2)
        
        ax.set_title('Hierarchical Clustering Heatmap (Z-score normalized)')
        ax.set_xlabel('Samples')
        ax.set_ylabel(self.term_for_features.capitalize())
        
        # Save PDF
        if save_pdf:
            pdf_file = self.output_dir / "hierarchical_clustering_heatmap.pdf"
            plt.savefig(pdf_file, dpi=300, bbox_inches='tight')
            logger.info(f"Heatmap saved as PDF: {pdf_file}")
        
        plt.close()
        
        # Save R-compatible data for ComplexHeatmap
        if save_r_compatible:
            self._save_r_compatible_data(sorted_data, sorted_labels)
    
    def _save_r_compatible_data(self, sorted_data: pd.DataFrame, 
                               sorted_labels: np.ndarray) -> None:
        """
        Save data in format compatible with R ComplexHeatmap.
        
        Args:
            sorted_data: Z-score data sorted by cluster
            sorted_labels: Cluster labels sorted by cluster
        """
        # Save z-score data for R
        r_data_file = self.output_dir / "zscore_data_for_R.txt"
        sorted_data.to_csv(r_data_file, sep='\t')
        
        # Save cluster annotations for R
        cluster_df = pd.DataFrame({
            'Feature': sorted_data.index,
            'Cluster': sorted_labels
        })
        cluster_df.to_csv(self.output_dir / "cluster_annotations_for_R.txt", sep='\t', index=False)
        
        # Create R script template
        r_script_file = self.output_dir / "complexheatmap_script.R"
        with open(r_script_file, 'w') as f:
            f.write("""# R script for ComplexHeatmap visualization with replicate median calculation
# Load required libraries
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)

# Read data
zscore_data <- read.table("zscore_data_for_R.txt", header=TRUE, row.names=1, sep="\\t")
cluster_anno <- read.table("cluster_annotations_for_R.txt", header=TRUE, sep="\\t")

# Function to calculate median for replicates
calculate_replicate_median <- function(data_matrix) {
  # Convert to data frame for easier manipulation
  df <- as.data.frame(data_matrix)
  df$Feature <- rownames(df)
  
  # Reshape data to long format
  df_long <- df %>%
    gather(key = "Sample", value = "Z_score", -Feature)
  
  # Extract condition and replicate information from sample names
  # Using improved logic for various naming patterns
  df_long <- df_long %>%
    mutate(
      Condition = sapply(Sample, function(sample_name) {
        # Method 1: Remove trailing numbers (replicate numbers)
        condition <- sub("-[0-9]+$", "", sample_name)
        
        # Method 2: If Method 1 doesn't work well, try splitting by hyphen
        if(condition == sample_name) {
          parts <- strsplit(sample_name, "-")[[1]]
          if(length(parts) >= 2) {
            # Check if last part is a number
            if(grepl("^[0-9]+$", parts[length(parts)])) {
              paste(parts[1:(length(parts)-1)], collapse = "-")
            } else {
              sample_name
            }
          } else {
            sample_name
          }
        } else {
          condition
        }
      }),
      Replicate = sapply(Sample, function(sample_name) {
        # Extract replicate number
        parts <- strsplit(sample_name, "-")[[1]]
        if(length(parts) >= 2 && grepl("^[0-9]+$", parts[length(parts)])) {
          parts[length(parts)]
        } else {
          "1"
        }
      })
    )
  
  # Calculate median for each feature-condition combination
  df_median <- df_long %>%
    group_by(Feature, Condition) %>%
    summarise(
      Median_Z_score = median(Z_score, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Reshape back to wide format
  df_median_wide <- df_median %>%
    spread(key = Condition, value = Median_Z_score)
  
  # Convert back to matrix
  median_matrix <- as.matrix(df_median_wide[, -1])
  rownames(median_matrix) <- df_median_wide$Feature
  
  return(median_matrix)
}

# Calculate median for replicates
cat("Calculating median values for replicates...\\n")
zscore_median <- calculate_replicate_median(zscore_data)

# Save median data
write.table(zscore_median, "zscore_median_data.txt", sep="\\t", quote=FALSE, col.names=NA)

# Create cluster annotation
cluster_colors <- rainbow(length(unique(cluster_anno$Cluster)))
names(cluster_colors) <- unique(cluster_anno$Cluster)
cluster_ha <- HeatmapAnnotation(
    cluster = cluster_anno$Cluster,
    col = list(cluster = cluster_colors),
    which = "row"
)

# Create original heatmap (with all replicates)
cat("Creating original heatmap with all replicates...\\n")
ht_original <- Heatmap(
    as.matrix(zscore_data),
    name = "Z-score (All Replicates)",
    col = colorRamp2(c(-1.5, 0, 1.5), c("#5ab4ac", "white", "#d8b365")),
    cluster_rows = FALSE,  # Use our hierarchical clustering
    cluster_columns = TRUE,
    show_row_names = FALSE,
    show_column_names = TRUE,
    left_annotation = cluster_ha,
    heatmap_legend_param = list(
        title = "Z-score",
        title_position = "leftcenter-rot"
    )
)

# Create median heatmap
cat("Creating median heatmap...\\n")
ht_median <- Heatmap(
    as.matrix(zscore_median),
    name = "Z-score (Median)",
    col = colorRamp2(c(-1.5, 0, 1.5), c("#5ab4ac", "white", "#d8b365")),
    cluster_rows = FALSE,  # Use our hierarchical clustering
    cluster_columns = TRUE,
    show_row_names = FALSE,
    show_column_names = TRUE,
    left_annotation = cluster_ha,
    heatmap_legend_param = list(
        title = "Z-score",
        title_position = "leftcenter-rot"
    )
)

# Draw original heatmap
pdf("complexheatmap_original_replicates.pdf", width=4, height=4)
draw(ht_original)
dev.off()

# Draw median heatmap
pdf("complexheatmap_median_replicates.pdf", width=4, height=4)
draw(ht_median)
dev.off()

# Create combined heatmap (side by side)
cat("Creating combined heatmap...\\n")
pdf("complexheatmap_combined.pdf", width=8, height=4)
draw(ht_original + ht_median)
dev.off()

# Print summary statistics
cat("\\n=== SUMMARY STATISTICS ===\\n")
cat("Original data dimensions:", dim(zscore_data), "\\n")
cat("Median data dimensions:", dim(zscore_median), "\\n")
cat("Number of clusters:", length(unique(cluster_anno$Cluster)), "\\n")
cat("Cluster sizes:", table(cluster_anno$Cluster), "\\n")

# Calculate and print correlation between original and median data
if(ncol(zscore_data) == ncol(zscore_median)) {
  correlations <- sapply(1:ncol(zscore_data), function(i) {
    cor(zscore_data[, i], zscore_median[, i], use = "complete.obs")
  })
  cat("Correlation between original and median data:", round(mean(correlations, na.rm=TRUE), 4), "\\n")
}

cat("\\nComplexHeatmap visualization completed!")
cat("\\nFiles generated:")
cat("\\n- complexheatmap_original_replicates.pdf: Original data with all replicates")
cat("\\n- complexheatmap_median_replicates.pdf: Median values for each condition")
cat("\\n- complexheatmap_combined.pdf: Side-by-side comparison")
cat("\\n- zscore_median_data.txt: Median data matrix")
""")
        
        logger.info(f"R-compatible data saved to {self.output_dir}")
        logger.info(f"R script template created: {r_script_file}")
        
        # Create additional R script with flexible replicate handling
        r_script_advanced_file = self.output_dir / "complexheatmap_advanced_script.R"
        with open(r_script_advanced_file, 'w') as f:
            f.write("""# Advanced R script for ComplexHeatmap visualization with flexible replicate handling
# Load required libraries
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)

# Read data
zscore_data <- read.table("zscore_data_for_R.txt", header=TRUE, row.names=1, sep="\\t")
cluster_anno <- read.table("cluster_annotations_for_R.txt", header=TRUE, sep="\\t")

# Function to calculate median for replicates with flexible naming patterns
calculate_replicate_median_flexible <- function(data_matrix, pattern_type = "auto") {
  # Convert to data frame for easier manipulation
  df <- as.data.frame(data_matrix)
  df$Feature <- rownames(df)
  
  # Reshape data to long format
  df_long <- df %>%
    gather(key = "Sample", value = "Z_score", -Feature)
  
  # Function to extract condition and replicate based on pattern type
  extract_condition_replicate <- function(sample_names, pattern_type) {
    if(pattern_type == "auto") {
      # Try to detect pattern automatically
      # Look for common patterns like: condition_replicate, condition-replicate, etc.
      samples <- unique(sample_names)
      if(any(grepl("_", samples))) {
        # Pattern: condition_replicate (e.g., "dko_ABA1", "hsfa_Mock2")
        condition <- sapply(strsplit(samples, "_"), function(x) {
          if(length(x) >= 2) {
            paste(x[1:(length(x)-1)], collapse = "_")
          } else {
            samples[1]
          }
        })
        replicate <- sapply(strsplit(samples, "_"), function(x) {
          if(length(x) >= 2) {
            x[length(x)]
          } else {
            "1"
          }
        })
      } else if(any(grepl("-", samples))) {
        # Pattern: condition-replicate (e.g., "dko-ABA1", "hsfa-Mock2")
        condition <- sapply(strsplit(samples, "-"), function(x) {
          if(length(x) >= 2) {
            paste(x[1:(length(x)-1)], collapse = "-")
          } else {
            samples[1]
          }
        })
        replicate <- sapply(strsplit(samples, "-"), function(x) {
          if(length(x) >= 2) {
            x[length(x)]
          } else {
            "1"
          }
        })
      } else {
        # No clear pattern, treat each sample as unique condition
        condition <- samples
        replicate <- rep("1", length(samples))
      }
    } else if(pattern_type == "underscore") {
      # Pattern: condition_replicate
      condition <- sapply(strsplit(sample_names, "_"), function(x) {
        if(length(x) >= 2) {
          paste(x[1:(length(x)-1)], collapse = "_")
        } else {
          sample_names[1]
        }
      })
      replicate <- sapply(strsplit(sample_names, "_"), function(x) {
        if(length(x) >= 2) {
          x[length(x)]
        } else {
          "1"
        }
      })
    } else if(pattern_type == "hyphen") {
      # Pattern: condition-replicate
      condition <- sapply(strsplit(sample_names, "-"), function(x) {
        if(length(x) >= 2) {
          paste(x[1:(length(x)-1)], collapse = "-")
        } else {
          sample_names[1]
        }
      })
      replicate <- sapply(strsplit(sample_names, "-"), function(x) {
        if(length(x) >= 2) {
          x[length(x)]
        } else {
          "1"
        }
      })
    } else {
      # Custom pattern or no replicates
      condition <- sample_names
      replicate <- rep("1", length(sample_names))
    }
    
    return(list(condition = condition, replicate = replicate))
  }
  
  # Extract condition and replicate information
  sample_info <- extract_condition_replicate(df_long$Sample, pattern_type)
  
  # Add condition and replicate columns
  df_long$Condition <- sample_info$condition[match(df_long$Sample, unique(df_long$Sample))]
  df_long$Replicate <- sample_info$replicate[match(df_long$Sample, unique(df_long$Sample))]
  
  # Calculate median for each feature-condition combination
  df_median <- df_long %>%
    group_by(Feature, Condition) %>%
    summarise(
      Median_Z_score = median(Z_score, na.rm = TRUE),
      N_replicates = n(),
      .groups = 'drop'
    )
  
  # Reshape back to wide format
  df_median_wide <- df_median %>%
    select(Feature, Condition, Median_Z_score) %>%
    spread(key = Condition, value = Median_Z_score)
  
  # Convert back to matrix
  median_matrix <- as.matrix(df_median_wide[, -1])
  rownames(median_matrix) <- df_median_wide$Feature
  
  return(list(
    median_matrix = median_matrix,
    replicate_info = df_median
  ))
}

# Function to create condition annotation
create_condition_annotation <- function(data_matrix) {
  # Extract condition names from column names
  conditions <- colnames(data_matrix)
  
  # Create color palette for conditions
  condition_colors <- rainbow(length(unique(conditions)))
  names(condition_colors) <- unique(conditions)
  
  # Create annotation
  condition_ha <- HeatmapAnnotation(
    condition = conditions,
    col = list(condition = condition_colors),
    which = "column"
  )
  
  return(condition_ha)
}

# Calculate median for replicates with different pattern options
cat("Calculating median values for replicates...\\n")

# Try different pattern types
pattern_types <- c("auto", "underscore", "hyphen", "none")
median_results <- list()

for(pattern in pattern_types) {
  cat("Trying pattern type:", pattern, "\\n")
  tryCatch({
    result <- calculate_replicate_median_flexible(zscore_data, pattern)
    median_results[[pattern]] <- result
    cat("Success with pattern:", pattern, "\\n")
    cat("Median matrix dimensions:", dim(result$median_matrix), "\\n")
    cat("Number of conditions:", ncol(result$median_matrix), "\\n")
    break
  }, error = function(e) {
    cat("Failed with pattern:", pattern, "-", e$message, "\\n")
  })
}

# Use the first successful result
if(length(median_results) > 0) {
  selected_pattern <- names(median_results)[1]
  zscore_median <- median_results[[selected_pattern]]$median_matrix
  replicate_info <- median_results[[selected_pattern]]$replicate_info
  
  cat("Using pattern type:", selected_pattern, "\\n")
} else {
  cat("Warning: Could not calculate median. Using original data.\\n")
  zscore_median <- zscore_data
  replicate_info <- NULL
}

# Save median data
write.table(zscore_median, "zscore_median_data.txt", sep="\\t", quote=FALSE, col.names=NA)

# Save replicate information
if(!is.null(replicate_info)) {
  write.table(replicate_info, "replicate_info.txt", sep="\\t", quote=FALSE, row.names=FALSE)
}

# Create cluster annotation
cluster_colors <- rainbow(length(unique(cluster_anno$Cluster)))
names(cluster_colors) <- unique(cluster_anno$Cluster)
cluster_ha <- HeatmapAnnotation(
    cluster = cluster_anno$Cluster,
    col = list(cluster = cluster_colors),
    which = "row"
)

# Create condition annotation for median data
condition_ha <- create_condition_annotation(zscore_median)

# Create original heatmap (with all replicates)
cat("Creating original heatmap with all replicates...\\n")
ht_original <- Heatmap(
    as.matrix(zscore_data),
    name = "Z-score (All Replicates)",
    col = colorRamp2(c(-1.5, 0, 1.5), c("#5ab4ac", "white", "#d8b365")),
    cluster_rows = FALSE,  # Use our hierarchical clustering
    cluster_columns = TRUE,
    show_row_names = FALSE,
    show_column_names = TRUE,
    left_annotation = cluster_ha,
    heatmap_legend_param = list(
        title = "Z-score",
        title_position = "leftcenter-rot"
    )
)

# Create median heatmap
cat("Creating median heatmap...\\n")
ht_median <- Heatmap(
    as.matrix(zscore_median),
    name = "Z-score (Median)",
    col = colorRamp2(c(-1.5, 0, 1.5), c("#5ab4ac", "white", "#d8b365")),
    cluster_rows = FALSE,  # Use our hierarchical clustering
    cluster_columns = TRUE,
    show_row_names = FALSE,
    show_column_names = TRUE,
    left_annotation = cluster_ha,
    top_annotation = condition_ha,
    heatmap_legend_param = list(
        title = "Z-score",
        title_position = "leftcenter-rot"
    )
)

# Draw original heatmap
pdf("complexheatmap_original_replicates.pdf", width=4, height=4)
draw(ht_original)
dev.off()

# Draw median heatmap
pdf("complexheatmap_median_replicates.pdf", width=4, height=4)
draw(ht_median)
dev.off()

# Create combined heatmap (side by side)
cat("Creating combined heatmap...\\n")
pdf("complexheatmap_combined.pdf", width=8, height=4)
draw(ht_original + ht_median)
dev.off()

# Create cluster-specific median plots
if(!is.null(replicate_info)) {
  cat("Creating cluster-specific median plots...\\n")
  
  # Merge cluster information with replicate info
  cluster_replicate_info <- replicate_info %>%
    left_join(cluster_anno, by = "Feature")
  
  # Create plots for each cluster
  unique_clusters <- unique(cluster_replicate_info$Cluster)
  
  for(cluster_id in unique_clusters) {
    cluster_data <- cluster_replicate_info %>%
      filter(Cluster == cluster_id)
    
    if(nrow(cluster_data) > 0) {
      # Create boxplot for this cluster
      pdf(paste0("cluster_", cluster_id, "_median_boxplot.pdf"), width=10, height=6)
      
      # Reshape data for plotting
      plot_data <- cluster_data %>%
        select(Feature, Condition, Median_Z_score) %>%
        spread(key = Condition, value = Median_Z_score)
      
      # Create boxplot
      boxplot(plot_data[, -1], 
              main = paste("Cluster", cluster_id, "Median Z-scores"),
              ylab = "Median Z-score",
              col = rainbow(ncol(plot_data) - 1))
      
      dev.off()
    }
  }
}

# Print summary statistics
cat("\\n=== SUMMARY STATISTICS ===\\n")
cat("Original data dimensions:", dim(zscore_data), "\\n")
cat("Median data dimensions:", dim(zscore_median), "\\n")
cat("Number of clusters:", length(unique(cluster_anno$Cluster)), "\\n")
cat("Cluster sizes:", table(cluster_anno$Cluster), "\\n")

if(!is.null(replicate_info)) {
  cat("Replicate information:\\n")
  print(table(replicate_info$Condition, replicate_info$N_replicates))
}

# Calculate and print correlation between original and median data
if(ncol(zscore_data) == ncol(zscore_median)) {
  correlations <- sapply(1:ncol(zscore_data), function(i) {
    cor(zscore_data[, i], zscore_median[, i], use = "complete.obs")
  })
  cat("Correlation between original and median data:", round(mean(correlations, na.rm=TRUE), 4), "\\n")
}

cat("\\nComplexHeatmap visualization completed!")
cat("\\nFiles generated:")
cat("\\n- complexheatmap_original_replicates.pdf: Original data with all replicates")
cat("\\n- complexheatmap_median_replicates.pdf: Median values for each condition")
cat("\\n- complexheatmap_combined.pdf: Side-by-side comparison")
cat("\\n- zscore_median_data.txt: Median data matrix")
if(!is.null(replicate_info)) {
  cat("\\n- replicate_info.txt: Detailed replicate information")
  cat("\\n- cluster_*_median_boxplot.pdf: Cluster-specific median plots")
}
""")
        
        logger.info(f"Advanced R script created: {r_script_advanced_file}")
        
        # Create simplified R script for median calculation
        r_script_simple_file = self.output_dir / "median_heatmap_script.R"
        with open(r_script_simple_file, 'w') as f:
            f.write("""# Simple R script for median calculation and heatmap visualization
# Load required libraries
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)

# Read data
zscore_data <- read.table("zscore_data_for_R.txt", header=TRUE, row.names=1, sep="\\t")
cluster_anno <- read.table("cluster_annotations_for_R.txt", header=TRUE, sep="\\t")

# Function to safely read data with error checking
safe_read_data <- function(file_path, expected_cols = NULL) {
  if(!file.exists(file_path)) {
    stop(paste("File not found:", file_path))
  }
  
  data <- tryCatch({
    read.table(file_path, header=TRUE, row.names=1, sep="\\t", check.names=FALSE)
  }, error = function(e) {
    stop(paste("Error reading", file_path, ":", e$message))
  })
  
  cat("Successfully read", file_path, "with dimensions:", dim(data), "\\n")
  cat("Sample names:", paste(colnames(data), collapse=", "), "\\n")
  return(data)
}

# Function to extract condition from sample name by removing last digit
extract_condition <- function(sample_names) {
  cat("Extracting conditions from sample names...\\n")
  cat("Sample names:", paste(sample_names, collapse=", "), "\\n")
  
  # Remove the last digit from each sample name
  conditions <- sub("[0-9]$", "", sample_names)
  
  cat("Extracted conditions:", paste(conditions, collapse=", "), "\\n")
  
  # Get unique conditions
  unique_conditions <- unique(conditions)
  cat("Unique conditions found:", paste(unique_conditions, collapse=", "), "\\n")
  
  return(conditions)
}

# Function to calculate median for replicates with error checking
calculate_median_by_condition <- function(data_matrix) {
  cat("Calculating median values...\\n")
  
  # Convert to data frame
  df <- as.data.frame(data_matrix)
  df$Feature <- rownames(df)
  
  cat("Original data dimensions:", dim(df), "\\n")
  cat("Sample columns:", paste(colnames(df)[1:min(5, ncol(df))], collapse=", "), "\\n")
  
  # Reshape to long format
  df_long <- df %>%
    gather(key = "Sample", value = "Z_score", -Feature)
  
  cat("Long format data rows:", nrow(df_long), "\\n")
  
  # Extract condition from sample name using improved logic
  df_long$Condition <- extract_condition(df_long$Sample)
  
  # Check for unique conditions
  unique_conditions <- unique(df_long$Condition)
  cat("Unique conditions found:", paste(unique_conditions, collapse=", "), "\\n")
  
  # Verify we have multiple conditions
  if(length(unique_conditions) <= 1) {
    cat("WARNING: Only one condition found. This may indicate a problem with condition extraction.\\n")
    cat("Using original sample names as conditions...\\n")
    df_long$Condition <- df_long$Sample
    unique_conditions <- unique(df_long$Condition)
  }
  
  # Calculate median for each feature-condition combination
  df_median <- df_long %>%
    group_by(Feature, Condition) %>%
    summarise(
      Median_Z_score = median(Z_score, na.rm = TRUE),
      N_replicates = n(),
      .groups = 'drop'
    )
  
  cat("Median calculation completed. Rows:", nrow(df_median), "\\n")
  cat("Conditions in median data:", paste(unique(df_median$Condition), collapse=", "), "\\n")
  
  # Reshape back to wide format
  df_median_wide <- df_median %>%
    select(Feature, Condition, Median_Z_score) %>%
    spread(key = Condition, value = Median_Z_score)
  
  cat("Wide format dimensions:", dim(df_median_wide), "\\n")
  cat("Wide format columns:", paste(colnames(df_median_wide), collapse=", "), "\\n")
  
  # Convert to matrix
  median_matrix <- as.matrix(df_median_wide[, -1])
  rownames(median_matrix) <- df_median_wide$Feature
  
  # IMPORTANT: Maintain original feature order from input data
  # Get the original feature order from the input data
  original_features <- rownames(data_matrix)
  
  # Reorder median matrix to match original feature order
  median_matrix <- median_matrix[original_features, ]
  
  cat("Maintaining original feature order for cluster compatibility\\n")
  cat("Final median matrix dimensions:", dim(median_matrix), "\\n")
  cat("Final median matrix columns:", paste(colnames(median_matrix), collapse=", "), "\\n")
  
  return(list(
    median_matrix = median_matrix,
    replicate_info = df_median
  ))
}

# Read data with error checking
cat("Reading data files...\\n")
zscore_data <- safe_read_data("zscore_data_for_R.txt")
cluster_anno <- read.table("cluster_annotations_for_R.txt", header=TRUE, sep="\\t", stringsAsFactors=FALSE)
cat("Successfully read cluster_annotations_for_R.txt with dimensions:", dim(cluster_anno), "\\n")

# Calculate median with error handling
cat("\\n=== MEDIAN CALCULATION ===\\n")
result <- tryCatch({
  calculate_median_by_condition(zscore_data)
}, error = function(e) {
  cat("Error in median calculation:", e$message, "\\n")
  cat("Using original data as fallback...\\n")
  list(
    median_matrix = zscore_data,
    replicate_info = NULL
  )
})

zscore_median <- result$median_matrix
replicate_info <- result$replicate_info

# Save median data
write.table(zscore_median, "zscore_median_data.txt", sep="\\t", quote=FALSE, col.names=NA)
if(!is.null(replicate_info)) {
  write.table(replicate_info, "replicate_info.txt", sep="\\t", quote=FALSE, row.names=FALSE)
}

# Ensure cluster annotations match the median data
cat("\\n=== CLUSTER ANNOTATION PROCESSING ===\\n")
cat("Original cluster annotations:", nrow(cluster_anno), "features\\n")
cat("Median data features:", nrow(zscore_median), "features\\n")

# Filter cluster annotations to only include features present in median data
cluster_anno_filtered <- cluster_anno[cluster_anno$Feature %in% rownames(zscore_median), ]
cat("Features in both cluster annotations and median data:", nrow(cluster_anno_filtered), "\\n")

if(nrow(cluster_anno_filtered) == 0) {
  cat("ERROR: No features match between cluster annotations and median data!\\n")
  cat("This will cause errors in the heatmap.\\n")
  cat("Cluster annotation features (first 10):", paste(head(cluster_anno$Feature, 10), collapse=", "), "\\n")
  cat("Median data features (first 10):", paste(head(rownames(zscore_median), 10), collapse=", "), "\\n")
  stop("No matching features found between cluster annotations and median data")
}

# Sort cluster annotations to match median data order
cluster_anno_filtered <- cluster_anno_filtered[match(rownames(zscore_median), cluster_anno_filtered$Feature), ]

# Check for any NA values
na_check <- is.na(cluster_anno_filtered$Cluster)
if(any(na_check)) {
  cat("WARNING: Found", sum(na_check), "features with NA cluster assignments\\n")
  # Remove features with NA clusters
  cluster_anno_filtered <- cluster_anno_filtered[!na_check, ]
  zscore_median <- zscore_median[!na_check, ]
}

# Create cluster annotation
unique_clusters <- unique(cluster_anno_filtered$Cluster)
cat("Unique clusters:", length(unique_clusters), "\\n")
cat("Cluster sizes:", table(cluster_anno_filtered$Cluster), "\\n")

if(length(unique_clusters) == 0) {
  cat("ERROR: No valid clusters found!\\n")
  stop("No valid clusters found in cluster annotations")
}

cluster_colors <- rainbow(length(unique_clusters))
names(cluster_colors) <- unique_clusters

cluster_ha <- HeatmapAnnotation(
    cluster = cluster_anno_filtered$Cluster,
    col = list(cluster = cluster_colors),
    which = "row"
)

# Create condition annotation
unique_conditions <- unique(colnames(zscore_median))
condition_colors <- rainbow(length(unique_conditions))
names(condition_colors) <- unique_conditions

condition_ha <- HeatmapAnnotation(
    condition = colnames(zscore_median),
    col = list(condition = condition_colors),
    which = "column"
)

# Create median heatmap with error checking
cat("\\n=== CREATING HEATMAP ===\\n")
ht_median <- tryCatch({
  Heatmap(
      as.matrix(zscore_median),
    name = "Z-score (Median)",
    col = colorRamp2(c(-1.5, 0, 1.5), c("#5ab4ac", "white", "#d8b365")),
    cluster_rows = FALSE,  # Use hierarchical clustering from Python
    cluster_columns = TRUE,
    show_row_names = FALSE,
    show_column_names = TRUE,
    left_annotation = cluster_ha,
    top_annotation = condition_ha,
    heatmap_legend_param = list(
        title = "Z-score",
        title_position = "leftcenter-rot"
    )
)

# Draw median heatmap
pdf("median_heatmap.pdf", width=4, height=4)
draw(ht_median)
dev.off()

# Create cluster-specific plots
cat("Creating cluster-specific plots...\\n")
cluster_replicate_info <- replicate_info %>%
  left_join(cluster_anno_filtered, by = "Feature")

unique_clusters <- unique(cluster_replicate_info$Cluster)

for(cluster_id in unique_clusters) {
  cluster_data <- cluster_replicate_info %>%
    filter(Cluster == cluster_id)
  
  if(nrow(cluster_data) > 0) {
    # Create boxplot
    pdf(paste0("cluster_", cluster_id, "_boxplot.pdf"), width=10, height=6)
    
    plot_data <- cluster_data %>%
      select(Feature, Condition, Median_Z_score) %>%
      spread(key = Condition, value = Median_Z_score)
    
    # Ensure we have data to plot
    if(ncol(plot_data) > 1) {
      boxplot(plot_data[, -1], 
              main = paste("Cluster", cluster_id, "Median Z-scores"),
              ylab = "Median Z-score",
              col = rainbow(ncol(plot_data) - 1))
    } else {
      plot(1, 1, type = "n", main = paste("Cluster", cluster_id, "- No data to plot"))
      text(1, 1, "No data available for this cluster")
    }
    
    dev.off()
    
    # Create heatmap for this cluster
    cluster_features <- cluster_data$Feature
    cluster_matrix <- zscore_median[cluster_features, , drop = FALSE]
    
    if(nrow(cluster_matrix) > 0) {
      pdf(paste0("cluster_", cluster_id, "_heatmap.pdf"), width=8, height=6)
      
      ht_cluster <- Heatmap(
        as.matrix(cluster_matrix),
        name = paste("Cluster", cluster_id),
        col = colorRamp2(c(-1.5, 0, 1.5), c("#5ab4ac", "white", "#d8b365")),
        cluster_rows = FALSE,
        cluster_columns = TRUE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        top_annotation = condition_ha,
        heatmap_legend_param = list(
          title = "Z-score",
          title_position = "leftcenter-rot"
        )
      )
      
      draw(ht_cluster)
      dev.off()
    }
  }
}

# Print summary
cat("\\n=== SUMMARY ===\\n")
cat("Original samples:", ncol(zscore_data), "\\n")
cat("Conditions (after median):", ncol(zscore_median), "\\n")
cat("Features:", nrow(zscore_median), "\\n")
cat("Clusters:", length(unique_clusters), "\\n")
cat("Cluster sizes:", table(cluster_anno_filtered$Cluster), "\\n")

cat("\\nFiles generated:")
cat("\\n- median_heatmap.pdf: Main median heatmap")
cat("\\n- zscore_median_data.txt: Median data matrix")
cat("\\n- replicate_info.txt: Replicate information")
cat("\\n- cluster_*_boxplot.pdf: Cluster-specific boxplots")
cat("\\n- cluster_*_heatmap.pdf: Cluster-specific heatmaps")
""")
        
        logger.info(f"Simple median R script created: {r_script_simple_file}")
        
        # Create improved R script with better error handling
        r_script_improved_file = self.output_dir / "median_heatmap_improved_script.R"
        with open(r_script_improved_file, 'w') as f:
            f.write("""# Improved R script for median calculation and heatmap visualization
# Load required libraries
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)

# Function to safely read data with error checking
safe_read_data <- function(file_path, expected_cols = NULL) {
  if(!file.exists(file_path)) {
    stop(paste("File not found:", file_path))
  }
  
  data <- tryCatch({
    read.table(file_path, header=TRUE, row.names=1, sep="\\t", check.names=FALSE)
  }, error = function(e) {
    stop(paste("Error reading", file_path, ":", e$message))
  })
  
  cat("Successfully read", file_path, "with dimensions:", dim(data), "\\n")
  return(data)
}

# Read data with error checking
cat("Reading data files...\\n")
zscore_data <- safe_read_data("zscore_data_for_R.txt")
cluster_anno <- read.table("cluster_annotations_for_R.txt", header=TRUE, sep="\\t", stringsAsFactors=FALSE)
cat("Successfully read cluster_annotations_for_R.txt with dimensions:", dim(cluster_anno), "\\n")

# Function to extract condition from sample name
extract_condition <- function(sample_names) {
  cat("Extracting conditions from sample names...\\n")
  cat("Sample names:", paste(sample_names, collapse=", "), "\\n")
  
  # Method 1: Remove trailing numbers (replicate numbers)
  conditions <- sub("-[0-9]+$", "", sample_names)
  
  # Method 2: If Method 1 doesn't work well, try splitting by hyphen
  # Remove the last digit from each sample name
  conditions <- sub("[0-9]$", "", sample_names)
  
  cat("Extracted conditions:", paste(conditions, collapse=", "), "\\n")
  
  # Get unique conditions
  unique_conditions <- unique(conditions)
  cat("Unique conditions found:", paste(unique_conditions, collapse=", "), "\\n")
  
  return(conditions)
}

# Function to calculate median for replicates with error checking
calculate_median_by_condition <- function(data_matrix) {
  cat("Calculating median values...\\n")
  
  # Convert to data frame
  df <- as.data.frame(data_matrix)
  df$Feature <- rownames(df)
  
  cat("Original data dimensions:", dim(df), "\\n")
  cat("Sample columns:", paste(colnames(df)[1:min(5, ncol(df))], collapse=", "), "\\n")
  
  # Reshape to long format
  df_long <- df %>%
    gather(key = "Sample", value = "Z_score", -Feature)
  
  cat("Long format data rows:", nrow(df_long), "\\n")
  
  # Extract condition from sample name using improved logic
  df_long$Condition <- extract_condition(df_long$Sample)
  
  # Check for unique conditions
  unique_conditions <- unique(df_long$Condition)
  cat("Unique conditions found:", paste(unique_conditions, collapse=", "), "\\n")
  
  # Verify we have multiple conditions
  if(length(unique_conditions) <= 1) {
    cat("WARNING: Only one condition found. This may indicate a problem with condition extraction.\\n")
    cat("Using original sample names as conditions...\\n")
    df_long$Condition <- df_long$Sample
    unique_conditions <- unique(df_long$Condition)
  }
  
  # Calculate median for each feature-condition combination
  df_median <- df_long %>%
    group_by(Feature, Condition) %>%
    summarise(
      Median_Z_score = median(Z_score, na.rm = TRUE),
      N_replicates = n(),
      .groups = 'drop'
    )
  
  cat("Median calculation completed. Rows:", nrow(df_median), "\\n")
  
  # Reshape back to wide format
  df_median_wide <- df_median %>%
    select(Feature, Condition, Median_Z_score) %>%
    spread(key = Condition, value = Median_Z_score)
  
  cat("Wide format dimensions:", dim(df_median_wide), "\\n")
  
  # Convert to matrix
  median_matrix <- as.matrix(df_median_wide[, -1])
  rownames(median_matrix) <- df_median_wide$Feature
  
  # IMPORTANT: Maintain original feature order from input data
  # Get the original feature order from the input data
  original_features <- rownames(data_matrix)
  
  # Reorder median matrix to match original feature order
  median_matrix <- median_matrix[original_features, ]
  
  cat("Maintaining original feature order for cluster compatibility\\n")
  cat("Final median matrix dimensions:", dim(median_matrix), "\\n")
  
  return(list(
    median_matrix = median_matrix,
    replicate_info = df_median
  ))
}

# Calculate median with error handling
cat("\\n=== MEDIAN CALCULATION ===\\n")
result <- tryCatch({
  calculate_median_by_condition(zscore_data)
}, error = function(e) {
  cat("Error in median calculation:", e$message, "\\n")
  cat("Using original data as fallback...\\n")
  list(
    median_matrix = zscore_data,
    replicate_info = NULL
  )
})

zscore_median <- result$median_matrix
replicate_info <- result$replicate_info

# Save median data
write.table(zscore_median, "zscore_median_data.txt", sep="\\t", quote=FALSE, col.names=NA)
if(!is.null(replicate_info)) {
  write.table(replicate_info, "replicate_info.txt", sep="\\t", quote=FALSE, row.names=FALSE)
}

# Ensure cluster annotations match the median data
cat("\\n=== CLUSTER ANNOTATION PROCESSING ===\\n")
cat("Original cluster annotations:", nrow(cluster_anno), "features\\n")
cat("Median data features:", nrow(zscore_median), "features\\n")

# Filter cluster annotations to only include features present in median data
cluster_anno_filtered <- cluster_anno[cluster_anno$Feature %in% rownames(zscore_median), ]
cat("Features in both cluster annotations and median data:", nrow(cluster_anno_filtered), "\\n")

if(length(cluster_anno_filtered) > 0 && nrow(cluster_anno_filtered) == 0) {
  cat("WARNING: No features match between cluster annotations and median data!\\n")
  cat("This will cause errors in the heatmap.\\n")
  cat("Cluster annotation features:", paste(head(cluster_anno$Feature, 10), collapse=", "), "\\n")
  cat("Median data features:", paste(head(rownames(zscore_median), 10), collapse=", "), "\\n")
}

# Sort cluster annotations to match median data order
cluster_anno_filtered <- cluster_anno_filtered[match(rownames(zscore_median), cluster_anno_filtered$Feature), ]

# Check for any NA values
na_check <- is.na(cluster_anno_filtered$Cluster)
if(any(na_check)) {
  cat("WARNING: Found", sum(na_check), "features with NA cluster assignments\\n")
  # Remove features with NA clusters
  cluster_anno_filtered <- cluster_anno_filtered[!na_check, ]
  zscore_median <- zscore_median[!na_check, ]
}

# Create cluster annotation
unique_clusters <- unique(cluster_anno_filtered$Cluster)
cat("Unique clusters:", length(unique_clusters), "\\n")
cat("Cluster sizes:", table(cluster_anno_filtered$Cluster), "\\n")

cluster_colors <- rainbow(length(unique_clusters))
names(cluster_colors) <- unique_clusters

cluster_ha <- HeatmapAnnotation(
    cluster = cluster_anno_filtered$Cluster,
    col = list(cluster = cluster_colors),
    which = "row"
)

# Create condition annotation
unique_conditions <- unique(colnames(zscore_median))
condition_colors <- rainbow(length(unique_conditions))
names(condition_colors) <- unique_conditions

condition_ha <- HeatmapAnnotation(
    condition = colnames(zscore_median),
    col = list(condition = condition_colors),
    which = "column"
)

# Create median heatmap with error checking
cat("\\n=== CREATING HEATMAP ===\\n")
ht_median <- tryCatch({
  Heatmap(
      as.matrix(zscore_median),
      name = "Z-score (Median)",
      col = colorRamp2(c(-1.5, 0, 1.5), c("#5ab4ac", "white", "#d8b365")),
      cluster_rows = FALSE,  # Use hierarchical clustering from Python
      cluster_columns = TRUE,
      show_row_names = FALSE,
      show_column_names = TRUE,
      left_annotation = cluster_ha,
      top_annotation = condition_ha,
      heatmap_legend_param = list(
          title = "Z-score",
          title_position = "leftcenter-rot"
      )
  )
}, error = function(e) {
  cat("Error creating heatmap:", e$message, "\\n")
  cat("Trying simplified heatmap without annotations...\\n")
  Heatmap(
      as.matrix(zscore_median),
      name = "Z-score (Median)",
      col = colorRamp2(c(-1.5, 0, 1.5), c("#5ab4ac", "white", "#d8b365")),
      cluster_rows = FALSE,
      cluster_columns = TRUE,
      show_row_names = FALSE,
      show_column_names = TRUE
  )
})

# Draw median heatmap
cat("Drawing median heatmap...\\n")
pdf("median_heatmap.pdf", width=4, height=4)
draw(ht_median)
dev.off()
cat("Median heatmap saved as median_heatmap.pdf\\n")

# Create cluster-specific plots only if we have replicate info
if(!is.null(replicate_info) && nrow(cluster_anno_filtered) > 0) {
  cat("\\n=== CREATING CLUSTER-SPECIFIC PLOTS ===\\n")
  cluster_replicate_info <- replicate_info %>%
    left_join(cluster_anno_filtered, by = "Feature")
  
  unique_clusters <- unique(cluster_replicate_info$Cluster)
  cat("Creating plots for clusters:", paste(unique_clusters, collapse=", "), "\\n")
  
  for(cluster_id in unique_clusters) {
    cluster_data <- cluster_replicate_info %>%
      filter(Cluster == cluster_id)
    
    cat("Processing cluster", cluster_id, "with", nrow(cluster_data), "features\\n")
    
    if(nrow(cluster_data) > 0) {
      # Create boxplot
      pdf(paste0("cluster_", cluster_id, "_boxplot.pdf"), width=10, height=6)
      
      plot_data <- cluster_data %>%
        select(Feature, Condition, Median_Z_score) %>%
        spread(key = Condition, value = Median_Z_score)
      
      # Ensure we have data to plot
      if(ncol(plot_data) > 1) {
        boxplot(plot_data[, -1], 
                main = paste("Cluster", cluster_id, "Median Z-scores"),
                ylab = "Median Z-score",
                col = rainbow(ncol(plot_data) - 1))
      } else {
        plot(1, 1, type = "n", main = paste("Cluster", cluster_id, "- No data to plot"))
        text(1, 1, "No data available for this cluster")
      }
      
      dev.off()
      
      # Create heatmap for this cluster
      cluster_features <- cluster_data$Feature
      cluster_matrix <- zscore_median[cluster_features, , drop = FALSE]
      
      if(nrow(cluster_matrix) > 0) {
        pdf(paste0("cluster_", cluster_id, "_heatmap.pdf"), width=8, height=6)
        
        ht_cluster <- Heatmap(
          as.matrix(cluster_matrix),
          name = paste("Cluster", cluster_id),
          col = colorRamp2(c(-1.5, 0, 1.5), c("#5ab4ac", "white", "#d8b365")),
          cluster_rows = FALSE,
          cluster_columns = TRUE,
          show_row_names = TRUE,
          show_column_names = TRUE,
          top_annotation = condition_ha,
          heatmap_legend_param = list(
            title = "Z-score",
            title_position = "leftcenter-rot"
          )
        )
        
        draw(ht_cluster)
        dev.off()
      }
    }
  }
} else {
  cat("Skipping cluster-specific plots (no replicate info or cluster data)\\n")
}

# Print comprehensive summary
cat("\\n=== FINAL SUMMARY ===\\n")
cat("Original data dimensions:", dim(zscore_data), "\\n")
cat("Median data dimensions:", dim(zscore_median), "\\n")
cat("Cluster annotations processed:", nrow(cluster_anno_filtered), "\\n")
cat("Unique clusters:", length(unique(cluster_anno_filtered$Cluster)), "\\n")
if(!is.null(replicate_info)) {
  cat("Replicate information available: Yes\\n")
} else {
  cat("Replicate information available: No\\n")
}

cat("\\nFiles generated:")
cat("\\n- median_heatmap.pdf: Main median heatmap")
cat("\\n- zscore_median_data.txt: Median data matrix")
if(!is.null(replicate_info)) {
  cat("\\n- replicate_info.txt: Replicate information")
  cat("\\n- cluster_*_boxplot.pdf: Cluster-specific boxplots")
  cat("\\n- cluster_*_heatmap.pdf: Cluster-specific heatmaps")
}
cat("\\n\\nAnalysis completed successfully!")
""")
        
        logger.info(f"Improved median R script created: {r_script_improved_file}")
    
    def run_complete_analysis(self, clustering_method: str = 'ward',
                             clustering_metric: str = 'euclidean',
                             n_clusters: Optional[int] = None,
                             distance_threshold: Optional[float] = None,
                             perform_bootstrap: bool = True,
                             n_bootstrap: int = 100) -> Dict:
        """
        Run the complete hierarchical clustering analysis pipeline.
        
        Args:
            clustering_method: Method for hierarchical clustering
            clustering_metric: Distance metric for clustering
            n_clusters: Number of clusters (if None, auto-determine)
            distance_threshold: Distance threshold for clustering
            perform_bootstrap: Whether to perform bootstrap analysis
            n_bootstrap: Number of bootstrap iterations
        
        Returns:
            Dictionary containing analysis results
        """
        logger.info("Starting complete hierarchical clustering analysis...")
        
        # Step 1: Load data
        self.load_vst_data()
        self.load_deg_dar_files()
        
        # Step 2: Filter and normalize
        self.filter_data_by_deg_dar()
        self.perform_z_score_normalization()
        
        # Step 3: Perform clustering
        self.perform_hierarchical_clustering(
            method=clustering_method,
            metric=clustering_metric,
            n_clusters=n_clusters,
            distance_threshold=distance_threshold
        )
        
        # Step 4: Bootstrap analysis (optional)
        bootstrap_results = None
        if perform_bootstrap:
            bootstrap_results = self.bootstrap_cluster_stability(n_bootstrap=n_bootstrap)
        
        # Step 5: Save results
        self.save_clustering_results()
        
        # Step 6: Create visualizations
        self.create_heatmap()
        
        # Compile results
        results = {
            'n_genes_original': len(self.vst_data),
            'n_genes_filtered': len(self.filtered_data),
            'n_clusters': len(set(self.cluster_labels)),
            'cluster_sizes': np.bincount(self.cluster_labels).tolist(),
            'bootstrap_results': bootstrap_results
        }
        
        logger.info("Analysis completed successfully!")
        return results


def main():
    """Main function to run the hierarchical clustering analysis."""
    parser = argparse.ArgumentParser(
        description="Hierarchical Clustering Analysis for VST Transformed Counts Data"
    )
    
    parser.add_argument(
        "--vst_data", 
        required=True,
        help="Path to VST transformed counts data file"
    )
    
    parser.add_argument(
        "--deg_dar_folder", 
        required=True,
        help="Path to folder containing DEG/DAR files"
    )
    
    parser.add_argument(
        "--output_dir", 
        required=True,
        help="Output directory for results"
    )
    
    parser.add_argument(
        "--clustering_method", 
        default="ward",
        choices=["ward", "complete", "average", "single"],
        help="Hierarchical clustering method"
    )
    
    parser.add_argument(
        "--clustering_metric", 
        default="euclidean",
        choices=["euclidean", "manhattan", "cosine"],
        help="Distance metric for clustering"
    )
    
    parser.add_argument(
        "--n_clusters", 
        type=int,
        help="Number of clusters (if not specified, auto-determine)"
    )
    
    parser.add_argument(
        "--distance_threshold", 
        type=float,
        help="Distance threshold for clustering"
    )
    
    parser.add_argument(
        "--no_bootstrap", 
        action="store_true",
        help="Skip bootstrap analysis"
    )
    
    parser.add_argument(
        "--n_bootstrap", 
        type=int,
        default=100,
        help="Number of bootstrap iterations"
    )
    
    args = parser.parse_args()
    
    # Create analyzer and run analysis
    analyzer = HierarchicalClusteringAnalyzer(
        vst_data_path=args.vst_data,
        deg_dar_folder_path=args.deg_dar_folder,
        output_dir=args.output_dir
    )
    
    try:
        results = analyzer.run_complete_analysis(
            clustering_method=args.clustering_method,
            clustering_metric=args.clustering_metric,
            n_clusters=args.n_clusters,
            distance_threshold=args.distance_threshold,
            perform_bootstrap=not args.no_bootstrap,
            n_bootstrap=args.n_bootstrap
        )
        
        # Print summary
        print("\n" + "="*50)
        print("ANALYSIS SUMMARY")
        print("="*50)
        print(f"Original genes: {results['n_genes_original']}")
        print(f"Filtered genes: {results['n_genes_filtered']}")
        print(f"Number of clusters: {results['n_clusters']}")
        print(f"Cluster sizes: {results['cluster_sizes']}")
        
        if results['bootstrap_results']:
            print(f"Bootstrap mean agreement: {results['bootstrap_results']['mean_agreement']:.4f}")
            print(f"Bootstrap std agreement: {results['bootstrap_results']['std_agreement']:.4f}")
        
        print(f"\nResults saved to: {args.output_dir}")
        print("="*50)
        
    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main() 