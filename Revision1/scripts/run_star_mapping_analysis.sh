#!/bin/bash

#!/bin/bash

# STAR Mapping Analysis Script
# This script runs the STAR mapping analysis to process Log.final.out files
# and generate CSV reports and PDF visualizations
#
# Usage:
#   bash scripts/run_star_mapping_analysis.sh
#   bash scripts/run_star_mapping_analysis.sh --input /path/to/input --output /path/to/output
#   bash scripts/run_star_mapping_analysis.sh -i /path/to/input -o /path/to/output

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
if ! command -v python3 &> /dev/null; then
    print_error "Python3 is not installed or not in PATH"
    exit 1
fi

# Check if required Python packages are installed
print_status "Checking Python dependencies..."
python -c "import pandas, matplotlib, seaborn" 2>/dev/null || {
    print_error "Required Python packages (pandas, matplotlib, seaborn) are not installed"
    print_status "Please install them using: pip install pandas matplotlib seaborn"
    exit 1
}

# Parse command line arguments
INPUT_DIR="data/RNA_STAR"
OUTPUT_DIR="output/STAR_mapping_analysis"

while [[ $# -gt 0 ]]; do
    case $1 in
        --input|-i)
            INPUT_DIR="$2"
            shift 2
            ;;
        --output|-o)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --help|-h)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --input, -i DIR     Input directory containing Log.final.out files"
            echo "  --output, -o DIR    Output directory for results"
            echo "  --help, -h          Show this help message"
            echo ""
            echo "Examples:"
            echo "  $0"
            echo "  $0 --input data/RNA_STAR --output output/my_analysis"
            echo "  $0 -i data/RNA_STAR -o output/my_analysis"
            exit 0
            ;;
        *)
            print_error "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Check if input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    print_error "Input directory $INPUT_DIR does not exist"
    exit 1
fi

# Check if there are any Log.final.out files
LOG_FILES=$(find "$INPUT_DIR" -name "*_Log.final.out" 2>/dev/null | wc -l)
if [ "$LOG_FILES" -eq 0 ]; then
    print_error "No Log.final.out files found in $INPUT_DIR"
    exit 1
fi

print_status "Found $LOG_FILES Log.final.out files"
print_status "Input directory: $INPUT_DIR"
print_status "Output directory: $OUTPUT_DIR"

# Create output directory
mkdir -p "$OUTPUT_DIR"

print_status "Starting STAR mapping analysis..."

# Run the Python script with arguments
python src/analyze_star_mapping.py --input "$INPUT_DIR" --output "$OUTPUT_DIR"

if [ $? -eq 0 ]; then
    print_status "STAR mapping analysis completed successfully!"
    print_status "Check the output directory: $OUTPUT_DIR"
else
    print_error "STAR mapping analysis failed"
    exit 1
fi

print_status "Analysis complete!" 