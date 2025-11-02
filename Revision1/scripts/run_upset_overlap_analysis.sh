#!/bin/bash
# Script to run RNA-seq and ATAC-seq upset file overlap analysis
# This script analyzes the overlap between HS_DEG_analysis and HS_DAR_analysis upset files

set -e  # Exit on any error

echo "Starting RNA-seq and ATAC-seq upset file overlap analysis..."

# Set paths
SCRIPT_DIR="src"
RNA_DATA_DIR="output/HS_DEG_analysis"
ATAC_DATA_DIR="output/HS_DAR_analysis"
OUTPUT_DIR="output/upset_overlap_analysis"

# Check if required directories exist
if [ ! -d "$RNA_DATA_DIR" ]; then
    echo "Error: RNA data directory not found: $RNA_DATA_DIR"
    exit 1
fi

if [ ! -d "$ATAC_DATA_DIR" ]; then
    echo "Error: ATAC data directory not found: $ATAC_DATA_DIR"
    exit 1
fi

# Check if up/ and down/ subdirectories exist
if [ ! -d "$RNA_DATA_DIR/up" ] && [ ! -d "$RNA_DATA_DIR/down" ]; then
    echo "Error: Neither 'up' nor 'down' directory found in RNA data directory"
    exit 1
fi

if [ ! -d "$ATAC_DATA_DIR/up" ] && [ ! -d "$ATAC_DATA_DIR/down" ]; then
    echo "Error: Neither 'up' nor 'down' directory found in ATAC data directory"
    exit 1
fi

echo "RNA data directory: $RNA_DATA_DIR"
echo "ATAC data directory: $ATAC_DATA_DIR"
echo "Output directory: $OUTPUT_DIR"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Run the analysis
echo "Running upset overlap analysis..."
python "$SCRIPT_DIR/analyze_rna_atac_overlap_upset.py" \
    --rna-clusters "$RNA_DATA_DIR" \
    --atac-clusters "$ATAC_DATA_DIR" \
    --output-dir "$OUTPUT_DIR" \
    --significance-threshold 0.05

if [ $? -eq 0 ]; then
    echo "✅ Analysis completed successfully!"
    echo "Results saved to: $OUTPUT_DIR"
    
    # List output files
    echo ""
    echo "Generated files:"
    find "$OUTPUT_DIR" -type f | sort
    
    echo ""
    echo "Summary report: $OUTPUT_DIR/upset_overlap_analysis_summary.txt"
    echo "Overlapping genes: $OUTPUT_DIR/overlapping_genes/"
    echo "Heatmaps: $OUTPUT_DIR/*.pdf"
    echo "R matrices: $OUTPUT_DIR/*.csv"
    
else
    echo "❌ Analysis failed!"
    exit 1
fi
