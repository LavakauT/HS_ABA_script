#!/bin/bash

# Peaks Annotation Analysis Runner Script
# This script provides easy execution of peaks annotation analysis

# Set script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Default values
PEAKS_DIR=""
GFF_FILE=""
OUTPUT_DIR="output/peaks_annotation"
PYTHON_ENV=""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Function to show usage
show_usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  -p, --peaks-dir DIR     Directory containing peak files (required)"
    echo "  -g, --gff-file FILE     Path to GFF annotation file (required)"
    echo "  -o, --output-dir DIR    Output directory (default: output/peaks_annotation)"
    echo "  -e, --env ENV           Python environment to activate (optional)"
    echo "  -h, --help              Show this help message"
    echo ""
    echo "Examples:"
    echo "  # Basic usage"
    echo "  $0 -p data/ABA_narrowpeak -g data/genome/MpTak_v6.1r2.gff"
    echo ""
    echo "  # With custom output directory"
    echo "  $0 -p data/HS_narrowpeak -g data/genome/MpTak_v6.1r2.gff -o output/HS_peaks_annotation"
    echo ""
    echo "  # With conda environment"
    echo "  $0 -p data/ABA_narrowpeak -g data/genome/MpTak_v6.1r2.gff -e my_conda_env"
}

# Function to check if file/directory exists
check_path() {
    local path="$1"
    local type="$2"
    
    if [ "$type" = "file" ]; then
        if [ ! -f "$path" ]; then
            print_error "File not found: $path"
            return 1
        fi
    elif [ "$type" = "dir" ]; then
        if [ ! -d "$path" ]; then
            print_error "Directory not found: $path"
            return 1
        fi
    fi
    return 0
}

# Function to check R packages
check_r_packages() {
    print_info "Checking R packages..."
    
    # Check if R is installed
    if ! command -v R &> /dev/null; then
        print_error "R is not installed or not in PATH"
        return 1
    fi
    
    # Check required R packages
    local packages=("GenomicFeatures" "ChIPseeker" "GenomicRanges")
    local missing_packages=()
    
    for package in "${packages[@]}"; do
        if ! R --slave -e "if (!require('$package', quietly=TRUE)) quit(status=1)" &> /dev/null; then
            missing_packages+=("$package")
        fi
    done
    
    if [ ${#missing_packages[@]} -gt 0 ]; then
        print_error "Missing R packages: ${missing_packages[*]}"
        print_info "Install missing packages in R:"
        echo "  if (!require('BiocManager', quietly = TRUE))"
        echo "      install.packages('BiocManager')"
        echo "  BiocManager::install(c('${missing_packages[*]// /', '}'))"
        return 1
    fi
    
    print_success "All required R packages are installed"
    return 0
}

# Function to activate Python environment
activate_env() {
    local env_name="$1"
    
    if [ -n "$env_name" ]; then
        print_info "Activating Python environment: $env_name"
        
        # Try conda first
        if command -v conda &> /dev/null; then
            source "$(conda info --base)/etc/profile.d/conda.sh"
            conda activate "$env_name" 2>/dev/null
            if [ $? -eq 0 ]; then
                print_success "Activated conda environment: $env_name"
                return 0
            fi
        fi
        
        # Try virtualenv
        if [ -d "$env_name" ]; then
            source "$env_name/bin/activate"
            if [ $? -eq 0 ]; then
                print_success "Activated virtual environment: $env_name"
                return 0
            fi
        fi
        
        print_warning "Could not activate environment: $env_name"
        print_info "Continuing with current Python environment"
    fi
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -p|--peaks-dir)
            PEAKS_DIR="$2"
            shift 2
            ;;
        -g|--gff-file)
            GFF_FILE="$2"
            shift 2
            ;;
        -o|--output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -e|--env)
            PYTHON_ENV="$2"
            shift 2
            ;;
        -h|--help)
            show_usage
            exit 0
            ;;
        *)
            print_error "Unknown option: $1"
            show_usage
            exit 1
            ;;
    esac
done

# Check required arguments
if [ -z "$PEAKS_DIR" ] || [ -z "$GFF_FILE" ]; then
    print_error "Missing required arguments"
    show_usage
    exit 1
fi

# Convert to absolute paths if relative
if [[ ! "$PEAKS_DIR" = /* ]]; then
    PEAKS_DIR="$PROJECT_ROOT/$PEAKS_DIR"
fi

if [[ ! "$GFF_FILE" = /* ]]; then
    GFF_FILE="$PROJECT_ROOT/$GFF_FILE"
fi

if [[ ! "$OUTPUT_DIR" = /* ]]; then
    OUTPUT_DIR="$PROJECT_ROOT/$OUTPUT_DIR"
fi

# Validate inputs
print_info "Validating inputs..."
check_path "$PEAKS_DIR" "dir" || exit 1
check_path "$GFF_FILE" "file" || exit 1

# Check for peak files in directory
peak_count=$(find "$PEAKS_DIR" -name "*.narrowPeak" -o -name "*.bed" -o -name "*.txt" | wc -l)
if [ "$peak_count" -eq 0 ]; then
    print_error "No peak files found in $PEAKS_DIR"
    print_info "Looking for files with extensions: .narrowPeak, .bed, .txt"
    exit 1
fi

print_success "Found $peak_count peak files in $PEAKS_DIR"

# Activate Python environment if specified
activate_env "$PYTHON_ENV"

# Check Python dependencies
print_info "Checking Python dependencies..."
python -c "import pandas, numpy, matplotlib, rpy2" 2>/dev/null
if [ $? -ne 0 ]; then
    print_error "Missing Python dependencies"
    print_info "Install with: pip install -r requirement/requirements_peaks_annotation.txt"
    exit 1
fi

# Check R packages
check_r_packages || exit 1

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Run the analysis
print_info "Starting peaks annotation analysis..."
print_info "Peaks directory: $PEAKS_DIR"
print_info "GFF file: $GFF_FILE"
print_info "Output directory: $OUTPUT_DIR"

cd "$PROJECT_ROOT"

python src/peaks_annotation.py \
    --peaks-dir "$PEAKS_DIR" \
    --gff-file "$GFF_FILE" \
    --output-dir "$OUTPUT_DIR"

if [ $? -eq 0 ]; then
    print_success "Peaks annotation analysis completed successfully!"
    print_info "Results saved to: $OUTPUT_DIR"
    
    # List output files
    if [ -d "$OUTPUT_DIR" ]; then
        print_info "Generated files:"
        find "$OUTPUT_DIR" -type f -name "*.txt" -o -name "*.csv" -o -name "*.pdf" -o -name "*.png" | sort | while read file; do
            echo "  - $(basename "$file")"
        done
    fi
else
    print_error "Peaks annotation analysis failed!"
    exit 1
fi 