#!/bin/bash

# Hierarchical Clustering Analysis Runner Script
# This script runs the hierarchical clustering analysis on VST transformed counts data

set -e  # Exit on any error

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
PYTHON_SCRIPT="$PROJECT_DIR/src/hierarchical_clustering_analysis.py"

# Default parameters
VST_DATA_PATH=""
DEG_DAR_FOLDER=""
OUTPUT_DIR=""
CLUSTERING_METHOD="ward"
CLUSTERING_METRIC="euclidean"
N_CLUSTERS=""
DISTANCE_THRESHOLD=""
NO_BOOTSTRAP=false
N_BOOTSTRAP=100

# Function to print usage
print_usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Required options:"
    echo "  --vst-data PATH          Path to VST transformed counts data file"
    echo "  --deg-dar-folder PATH    Path to folder containing DEG/DAR files"
    echo "  --output-dir PATH        Output directory for results"
    echo ""
    echo "Optional options:"
    echo "  --clustering-method METHOD   Clustering method (ward, complete, average, single) [default: ward]"
    echo "  --clustering-metric METRIC   Distance metric (euclidean, manhattan, cosine) [default: euclidean]"
    echo "  --n-clusters N              Number of clusters (auto-determine if not specified)"
    echo "  --distance-threshold THRESH  Distance threshold for clustering"
    echo "  --no-bootstrap              Skip bootstrap analysis"
    echo "  --n-bootstrap N             Number of bootstrap iterations [default: 100]"
    echo "  --help                      Show this help message"
    echo ""
    echo "Examples:"
    echo "  $0 --vst-data data/ABA_counts/RNA/vst_norm_counts_ABA_RNA.txt \\"
    echo "       --deg-dar-folder data/ABA_DEG \\"
    echo "       --output-dir output/ABA_hierarchical_clustering"
    echo ""
    echo "  $0 --vst-data data/HS_counts/RNA/vst_norm_counts_hs_RNA.txt \\"
    echo "       --deg-dar-folder data/HS_DEG \\"
    echo "       --output-dir output/HS_hierarchical_clustering \\"
    echo "       --clustering-method complete --n-clusters 5"
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --vst-data)
            VST_DATA_PATH="$2"
            shift 2
            ;;
        --deg-dar-folder)
            DEG_DAR_FOLDER="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --clustering-method)
            CLUSTERING_METHOD="$2"
            shift 2
            ;;
        --clustering-metric)
            CLUSTERING_METRIC="$2"
            shift 2
            ;;
        --n-clusters)
            N_CLUSTERS="$2"
            shift 2
            ;;
        --distance-threshold)
            DISTANCE_THRESHOLD="$2"
            shift 2
            ;;
        --no-bootstrap)
            NO_BOOTSTRAP=true
            shift
            ;;
        --n-bootstrap)
            N_BOOTSTRAP="$2"
            shift 2
            ;;
        --help)
            print_usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            print_usage
            exit 1
            ;;
    esac
done

# Check required parameters
if [[ -z "$VST_DATA_PATH" ]]; then
    echo "Error: --vst-data is required"
    print_usage
    exit 1
fi

if [[ -z "$DEG_DAR_FOLDER" ]]; then
    echo "Error: --deg-dar-folder is required"
    print_usage
    exit 1
fi

if [[ -z "$OUTPUT_DIR" ]]; then
    echo "Error: --output-dir is required"
    print_usage
    exit 1
fi

# Check if files exist
if [[ ! -f "$VST_DATA_PATH" ]]; then
    echo "Error: VST data file not found: $VST_DATA_PATH"
    exit 1
fi

if [[ ! -d "$DEG_DAR_FOLDER" ]]; then
    echo "Error: DEG/DAR folder not found: $DEG_DAR_FOLDER"
    exit 1
fi

# Check if Python script exists
if [[ ! -f "$PYTHON_SCRIPT" ]]; then
    echo "Error: Python script not found: $PYTHON_SCRIPT"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Build command
CMD="python \"$PYTHON_SCRIPT\""
CMD="$CMD --vst_data \"$VST_DATA_PATH\""
CMD="$CMD --deg_dar_folder \"$DEG_DAR_FOLDER\""
CMD="$CMD --output_dir \"$OUTPUT_DIR\""
CMD="$CMD --clustering_method \"$CLUSTERING_METHOD\""
CMD="$CMD --clustering_metric \"$CLUSTERING_METRIC\""

if [[ -n "$N_CLUSTERS" ]]; then
    CMD="$CMD --n_clusters $N_CLUSTERS"
fi

if [[ -n "$DISTANCE_THRESHOLD" ]]; then
    CMD="$CMD --distance_threshold $DISTANCE_THRESHOLD"
fi

if [[ "$NO_BOOTSTRAP" == true ]]; then
    CMD="$CMD --no_bootstrap"
fi

CMD="$CMD --n_bootstrap $N_BOOTSTRAP"

# Print configuration
echo "="*60
echo "HIERARCHICAL CLUSTERING ANALYSIS"
echo "="*60
echo "VST Data: $VST_DATA_PATH"
echo "DEG/DAR Folder: $DEG_DAR_FOLDER"
echo "Output Directory: $OUTPUT_DIR"
echo "Clustering Method: $CLUSTERING_METHOD"
echo "Clustering Metric: $CLUSTERING_METRIC"
if [[ -n "$N_CLUSTERS" ]]; then
    echo "Number of Clusters: $N_CLUSTERS"
fi
if [[ -n "$DISTANCE_THRESHOLD" ]]; then
    echo "Distance Threshold: $DISTANCE_THRESHOLD"
fi
echo "Bootstrap Analysis: $([ "$NO_BOOTSTRAP" == true ] && echo "Disabled" || echo "Enabled ($N_BOOTSTRAP iterations)")"
echo "="*60

# Run the analysis
echo "Starting hierarchical clustering analysis..."
echo "Command: $CMD"
echo ""

eval $CMD

# Check if analysis was successful
if [[ $? -eq 0 ]]; then
    echo ""
    echo "="*60
    echo "ANALYSIS COMPLETED SUCCESSFULLY!"
    echo "="*60
    echo "Results saved to: $OUTPUT_DIR"
    echo ""
    echo "Generated files:"
    ls -la "$OUTPUT_DIR"
    echo ""
    echo "For R visualization, you have several options:"
    echo "1. Final corrected median calculation and heatmap (recommended):"
    echo "   Rscript $OUTPUT_DIR/median_heatmap_final_fix.R"
    echo "2. Fixed median calculation and heatmap (for troubleshooting):"
    echo "   Rscript $OUTPUT_DIR/median_heatmap_fixed_script.R"
    echo "3. Improved median calculation and heatmap:"
    echo "   Rscript $OUTPUT_DIR/median_heatmap_improved_script.R"
    echo "4. Basic median calculation and heatmap:"
    echo "   Rscript $OUTPUT_DIR/median_heatmap_script.R"
    echo "5. Advanced visualization with flexible replicate handling:"
    echo "   Rscript $OUTPUT_DIR/complexheatmap_advanced_script.R"
    echo "6. Original visualization (all replicates):"
    echo "   Rscript $OUTPUT_DIR/complexheatmap_script.R"
    echo ""
    echo "The median scripts will calculate replicate medians and generate:"
    echo "- Median heatmaps"
    echo "- Cluster-specific plots"
    echo "- Replicate information files"
    echo ""
    echo "Note: The final corrected script (median_heatmap_final_fix.R) uses the"
    echo "corrected condition extraction logic that properly removes replicate numbers."
    echo "If you encounter issues with other scripts, try the final corrected version."
else
    echo ""
    echo "="*60
    echo "ANALYSIS FAILED!"
    echo "="*60
    exit 1
fi 