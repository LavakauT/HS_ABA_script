#!/bin/bash

# Differential Analysis Script (DARs and DEGs)
# This script provides an easy way to run DARs (ATAC-seq) or DEGs (RNA-seq) analysis on different datasets

set -e  # Exit on any error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Check if Python is available
if ! command -v python &> /dev/null; then
    print_error "Python is not installed or not in PATH"
    exit 1
fi

# Check if required packages are installed
print_status "Checking required packages..."
python -c "import pandas, matplotlib, seaborn, matplotlib_venn, numpy, upsetplot" 2>/dev/null || {
    print_warning "Some required packages are missing. Installing..."
    pip install pandas matplotlib seaborn matplotlib-venn numpy upsetplot
}

# Function to run analysis
run_analysis() {
    local data_folder=$1
    local output_folder=$2
    local analysis_type=$3
    
    if [ ! -d "$data_folder" ]; then
        print_error "Data folder not found: $data_folder"
        return 1
    fi
    
    print_status "Analyzing $data_folder for $analysis_type..."
    python src/analyze_differential.py "$data_folder" --output "$output_folder" --type "$analysis_type"
    
    if [ $? -eq 0 ]; then
        print_status "$analysis_type analysis completed successfully for $data_folder"
        print_status "Results saved to: $output_folder"
    else
        print_error "$analysis_type analysis failed for $data_folder"
        return 1
    fi
}

# Main script
main() {
    print_status "Starting Differential Analysis"
    print_status "============================="
    
    # Create output directory
    mkdir -p output
    
    # Analyze ABA_DAR data (ATAC-seq)
    if [ -d "data/ABA_DAR" ]; then
        run_analysis "data/ABA_DAR" "output/ABA_DAR_analysis" "DAR"
    else
        print_warning "ABA_DAR data folder not found"
    fi
    
    # Analyze HS_DAR data (ATAC-seq)
    if [ -d "data/HS_DAR" ]; then
        run_analysis "data/HS_DAR" "output/HS_DAR_analysis" "DAR"
    else
        print_warning "HS_DAR data folder not found"
    fi
    
    # Analyze ABA_DEG data (RNA-seq)
    if [ -d "data/ABA_DEG" ]; then
        run_analysis "data/ABA_DEG" "output/ABA_DEG_analysis" "DEG"
    else
        print_warning "ABA_DEG data folder not found"
    fi
    
    # Analyze HS_DEG data (RNA-seq)
    if [ -d "data/HS_DEG" ]; then
        run_analysis "data/HS_DEG" "output/HS_DEG_analysis" "DEG"
    else
        print_warning "HS_DEG data folder not found"
    fi
    
    print_status "All analyses completed!"
    print_status "Check the output/ directory for results"
}

# Run main function
main "$@" 