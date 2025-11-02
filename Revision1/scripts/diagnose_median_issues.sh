#!/bin/bash

# Diagnostic script for median analysis issues
# This script helps identify common problems with the median calculation

set -e

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    local color=$1
    local message=$2
    echo -e "${color}${message}${NC}"
}

# Function to check if file exists
check_file() {
    local file_path=$1
    local description=$2
    
    if [[ -f "$file_path" ]]; then
        print_status $GREEN "✓ $description: $file_path"
        return 0
    else
        print_status $RED "✗ $description: $file_path (NOT FOUND)"
        return 1
    fi
}

# Function to check R package installation
check_r_package() {
    local package=$1
    local result=$(Rscript -e "if(require($package, quietly=TRUE)) cat('OK') else cat('MISSING')" 2>/dev/null)
    
    if [[ "$result" == "OK" ]]; then
        print_status $GREEN "✓ R package: $package"
        return 0
    else
        print_status $RED "✗ R package: $package (NOT INSTALLED)"
        return 1
    fi
}

# Function to analyze sample naming patterns
analyze_sample_names() {
    local data_file=$1
    local description=$2
    
    if [[ ! -f "$data_file" ]]; then
        print_status $RED "Cannot analyze sample names: $data_file not found"
        return 1
    fi
    
    print_status $BLUE "Analyzing sample names in $description..."
    
    # Extract sample names and analyze patterns
    local sample_names=$(head -1 "$data_file" | cut -f2- | tr '\t' '\n')
    local total_samples=$(echo "$sample_names" | wc -l)
    
    echo "Total samples: $total_samples"
    echo "Sample names:"
    echo "$sample_names" | head -10
    
    # Check for naming patterns
    local underscore_count=$(echo "$sample_names" | grep -c "_" || true)
    local hyphen_count=$(echo "$sample_names" | grep -c "-" || true)
    local no_separator_count=$(echo "$sample_names" | grep -v "_" | grep -v "-" | wc -l)
    
    echo "Naming pattern analysis:"
    echo "  - Underscore separator: $underscore_count samples"
    echo "  - Hyphen separator: $hyphen_count samples"
    echo "  - No separator: $no_separator_count samples"
    
    # Suggest appropriate script
    if [[ $underscore_count -gt 0 ]]; then
        print_status $GREEN "Recommended: Use underscore pattern detection"
    elif [[ $hyphen_count -gt 0 ]]; then
        print_status $GREEN "Recommended: Use hyphen pattern detection"
    else
        print_status $YELLOW "Warning: No clear naming pattern detected"
    fi
}

# Function to check data dimensions
check_data_dimensions() {
    local data_file=$1
    local description=$2
    
    if [[ ! -f "$data_file" ]]; then
        print_status $RED "Cannot check dimensions: $data_file not found"
        return 1
    fi
    
    local rows=$(wc -l < "$data_file")
    local cols=$(head -1 "$data_file" | tr '\t' '\n' | wc -l)
    
    echo "$description dimensions: ${rows} rows x ${cols} columns"
    
    if [[ $rows -lt 10 ]]; then
        print_status $YELLOW "Warning: Very few genes ($rows)"
    fi
    
    if [[ $cols -lt 3 ]]; then
        print_status $YELLOW "Warning: Very few samples ($cols)"
    fi
}

# Function to test R script execution
test_r_script() {
    local script_file=$1
    local description=$2
    
    if [[ ! -f "$script_file" ]]; then
        print_status $RED "✗ $description: $script_file (NOT FOUND)"
        return 1
    fi
    
    print_status $BLUE "Testing $description..."
    
    # Create a temporary test directory
    local test_dir=$(mktemp -d)
    cd "$test_dir"
    
    # Copy required files if they exist
    local required_files=("zscore_data_for_R.txt" "cluster_annotations_for_R.txt")
    local files_copied=0
    
    for file in "${required_files[@]}"; do
        if [[ -f "$PROJECT_DIR/output/*/$file" ]]; then
            cp "$PROJECT_DIR/output/*/$file" . 2>/dev/null || true
            files_copied=$((files_copied + 1))
        fi
    done
    
    if [[ $files_copied -eq 0 ]]; then
        print_status $YELLOW "Warning: No test data found, skipping script test"
        cd - > /dev/null
        rm -rf "$test_dir"
        return 1
    fi
    
    # Test script syntax
    local syntax_check=$(Rscript -e "source('$script_file', echo=FALSE)" 2>&1)
    
    if [[ $? -eq 0 ]]; then
        print_status $GREEN "✓ $description: Syntax OK"
        cd - > /dev/null
        rm -rf "$test_dir"
        return 0
    else
        print_status $RED "✗ $description: Syntax error"
        echo "Error details:"
        echo "$syntax_check" | head -5
        cd - > /dev/null
        rm -rf "$test_dir"
        return 1
    fi
}

# Main diagnostic function
main() {
    print_status $BLUE "=========================================="
    print_status $BLUE "MEDIAN ANALYSIS DIAGNOSTIC TOOL"
    print_status $BLUE "=========================================="
    echo ""
    
    # Check for output directories
    print_status $BLUE "Checking output directories..."
    local output_dirs=($(find "$PROJECT_DIR/output" -name "*hierarchical*" -type d 2>/dev/null || true))
    
    if [[ ${#output_dirs[@]} -eq 0 ]]; then
        print_status $RED "No hierarchical clustering output directories found"
        echo "Please run hierarchical clustering analysis first:"
        echo "  ./scripts/run_hierarchical_clustering.sh"
        exit 1
    fi
    
    # Check each output directory
    for output_dir in "${output_dirs[@]}"; do
        echo ""
        print_status $BLUE "Analyzing: $output_dir"
        
        # Check required files
        check_file "$output_dir/zscore_data_for_R.txt" "Z-score data"
        check_file "$output_dir/cluster_annotations_for_R.txt" "Cluster annotations"
        check_file "$output_dir/median_heatmap_fixed_script.R" "Fixed median script"
check_file "$output_dir/median_heatmap_script.R" "Basic median script"
check_file "$output_dir/median_heatmap_improved_script.R" "Improved median script"
        check_file "$output_dir/complexheatmap_advanced_script.R" "Advanced script"
        
        # Analyze sample naming patterns
        if [[ -f "$output_dir/zscore_data_for_R.txt" ]]; then
            echo ""
            analyze_sample_names "$output_dir/zscore_data_for_R.txt" "Z-score data"
            check_data_dimensions "$output_dir/zscore_data_for_R.txt" "Z-score data"
        fi
        
        # Check cluster annotations
        if [[ -f "$output_dir/cluster_annotations_for_R.txt" ]]; then
            echo ""
            check_data_dimensions "$output_dir/cluster_annotations_for_R.txt" "Cluster annotations"
        fi
        
        # Test R scripts
        echo ""
        print_status $BLUE "Testing R scripts..."
        test_r_script "$output_dir/median_heatmap_fixed_script.R" "Fixed median script"
test_r_script "$output_dir/median_heatmap_improved_script.R" "Improved median script"
test_r_script "$output_dir/median_heatmap_script.R" "Basic median script"
    done
    
    # Check R package installation
    echo ""
    print_status $BLUE "Checking R package installation..."
    local required_packages=("ComplexHeatmap" "circlize" "dplyr" "tidyr")
    local missing_packages=0
    
    for package in "${required_packages[@]}"; do
        if ! check_r_package "$package"; then
            missing_packages=$((missing_packages + 1))
        fi
    done
    
    if [[ $missing_packages -gt 0 ]]; then
        echo ""
        print_status $YELLOW "Missing R packages detected. Install with:"
        echo "Rscript -e \"install.packages(c('ComplexHeatmap', 'circlize', 'dplyr', 'tidyr'))\""
    fi
    
    # Summary and recommendations
    echo ""
    print_status $BLUE "=========================================="
    print_status $BLUE "RECOMMENDATIONS"
    print_status $BLUE "=========================================="
    echo ""
    echo "1. If you encounter issues with the basic median script, try the improved version:"
    echo "   Rscript median_heatmap_fixed_script.R"
    echo "   Rscript median_heatmap_improved_script.R"
    echo ""
    echo "2. For complex naming patterns, use the advanced script:"
    echo "   Rscript complexheatmap_advanced_script.R"
    echo ""
    echo "3. If you see 'No genes match' warnings, check your sample naming convention"
    echo ""
    echo "4. For debugging, run the improved script which provides detailed output"
    echo ""
    print_status $GREEN "Diagnostic completed!"
}

# Run main function
main "$@" 