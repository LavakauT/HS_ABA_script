#!/bin/bash

# Compare Different Clustering Numbers for HS RNA Analysis
# This script runs hierarchical clustering with different numbers of clusters

set -e  # Exit on any error

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
PYTHON_SCRIPT="$PROJECT_DIR/src/hierarchical_clustering_analysis.py"

# Data paths
VST_DATA="data/HS_counts/RNA/vst_norm_counts_hs_RNA.txt"
DEG_DAR_FOLDER="data/HS_DEG"
BASE_OUTPUT_DIR="output/HS_clustering_comparison"

# Clustering numbers to test
CLUSTER_NUMBERS=(2 3 4 5 6 7 8)

echo "="*60
echo "COMPARING DIFFERENT CLUSTERING NUMBERS FOR HS RNA ANALYSIS"
echo "="*60
echo "VST Data: $VST_DATA"
echo "DEG/DAR Folder: $DEG_DAR_FOLDER"
echo "Base Output Directory: $BASE_OUTPUT_DIR"
echo "Testing cluster numbers: ${CLUSTER_NUMBERS[@]}"
echo "="*60

# Create base output directory
mkdir -p "$BASE_OUTPUT_DIR"

# Function to run clustering with specified number of clusters
run_clustering() {
    local n_clusters=$1
    local output_dir="$BASE_OUTPUT_DIR/clusters_${n_clusters}"
    
    echo ""
    echo "Testing with $n_clusters clusters..."
    echo "Output directory: $output_dir"
    
    # Run the analysis
    python3 "$PYTHON_SCRIPT" \
        --vst_data "$VST_DATA" \
        --deg_dar_folder "$DEG_DAR_FOLDER" \
        --output_dir "$output_dir" \
        --clustering_method ward \
        --clustering_metric euclidean \
        --n_clusters "$n_clusters" \
        --no_bootstrap  # Skip bootstrap for speed
    
    # Check if analysis was successful
    if [[ $? -eq 0 ]]; then
        echo "✅ Analysis completed successfully for $n_clusters clusters"
        
        # Get cluster sizes
        if [[ -f "$output_dir/cluster_summary.txt" ]]; then
            echo "Cluster sizes:"
            tail -n +2 "$output_dir/cluster_summary.txt" | cut -f2 | paste -sd ' ' -
        fi
    else
        echo "❌ Analysis failed for $n_clusters clusters"
    fi
}

# Run clustering for each number
for n_clusters in "${CLUSTER_NUMBERS[@]}"; do
    run_clustering "$n_clusters"
done

# Create summary
echo ""
echo "="*60
echo "SUMMARY OF ALL CLUSTERING RESULTS"
echo "="*60

for n_clusters in "${CLUSTER_NUMBERS[@]}"; do
    output_dir="$BASE_OUTPUT_DIR/clusters_${n_clusters}"
    summary_file="$output_dir/cluster_summary.txt"
    
    if [[ -f "$summary_file" ]]; then
        echo ""
        echo "Clusters: $n_clusters"
        echo "Output: $output_dir"
        echo "Files generated:"
        ls -la "$output_dir" | grep -E "\.(txt|pdf)$" | wc -l
        echo "Cluster sizes:"
        tail -n +2 "$summary_file" | cut -f2 | paste -sd ' ' -
    fi
done

echo ""
echo "="*60
echo "COMPARISON COMPLETED!"
echo "="*60
echo "All results saved to: $BASE_OUTPUT_DIR"
echo ""
echo "To view specific results:"
for n_clusters in "${CLUSTER_NUMBERS[@]}"; do
    echo "  $n_clusters clusters: $BASE_OUTPUT_DIR/clusters_${n_clusters}/"
done
echo ""
echo "To create R visualizations:"
for n_clusters in "${CLUSTER_NUMBERS[@]}"; do
    echo "  Rscript $BASE_OUTPUT_DIR/clusters_${n_clusters}/complexheatmap_script.R"
done 