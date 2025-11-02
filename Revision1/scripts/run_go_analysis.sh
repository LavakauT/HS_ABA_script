#!/bin/bash

# GO Term Analysis Runner Script
# This script runs the GO term analysis on the specified data folders

echo "Starting GO Term Analysis..."
echo "=============================="

# Set the project root directory
PROJECT_ROOT="/Users/user/Documents/HS_ABA/Revision1"
cd "$PROJECT_ROOT"

# Check if Python script exists
if [ ! -f "src/go_term_analysis.py" ]; then
    echo "Error: GO analysis script not found at src/go_term_analysis.py"
    exit 1
fi

# Check if required data files exist
if [ ! -f "data/Mp_Ath_Ortho/MpTak_v6.1r1.protein__v__Athaliana_447_Araport11.protein.csv" ]; then
    echo "Error: Orthology mapping file not found"
    exit 1
fi

if [ ! -f "data/Mp_Ath_Ortho/Athaliana_GO.txt" ]; then
    echo "Error: GO annotation file not found"
    exit 1
fi

# Check if required directories exist
REQUIRED_DIRS=(
    "output/HS_DEG_hierarchical_clustering_4/cluster_results"
    "output/HS_DAR_hierarchical_clustering_4/cluster_results"
    "output/HS_DAR_analysis/up"
    "output/HS_DAR_analysis/down"
    "output/HS_DEG_analysis/up"
    "output/HS_DEG_analysis/down"
)

echo "Checking required directories..."
for dir in "${REQUIRED_DIRS[@]}"; do
    if [ ! -d "$dir" ]; then
        echo "Warning: Directory $dir not found, skipping..."
    else
        echo "✓ Found: $dir"
    fi
done

echo ""
echo "Running GO Term Analysis..."
echo "=============================="

# Run the Python script
python src/go_term_analysis.py

# Check if the script ran successfully
if [ $? -eq 0 ]; then
    echo ""
    echo "✓ GO Term Analysis completed successfully!"
    echo ""
    echo "Results are available in:"
    echo "  - Tables: output/go_analysis/tables/"
    echo "  - Plots: output/go_analysis/plots/"
    echo "  - Summary: output/go_analysis/go_enrichment_summary.xlsx"
    echo ""
else
    echo ""
    echo "✗ GO Term Analysis failed!"
    echo "Check the error messages above for details."
    exit 1
fi

echo "Analysis complete!"
