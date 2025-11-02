#!/bin/bash

# Peak Analysis Script
# This script runs the peak analysis tool on narrowPeak files

set -e  # Exit on any error

# Default values
INPUT_DIR=""
OUTPUT_DIR=""
BED_FILE=""
VERBOSE=""

# Function to print usage
print_usage() {
    echo "Usage: $0 -i <input_dir> -o <output_dir> [-b <bed_file>] [-v]"
    echo ""
    echo "Options:"
    echo "  -i <input_dir>    Directory containing narrowPeak files (required)"
    echo "  -o <output_dir>   Output directory for results (required)"
    echo "  -b <bed_file>     Optional merged bed file for reference line"
    echo "  -v                Enable verbose logging"
    echo ""
    echo "Examples:"
    echo "  $0 -i data/ABA_narrowpeak/ -o output/peaks/"
    echo "  $0 -i data/ABA_narrowpeak/ -o output/peaks/ -b data/ABA_narrowpeak/ABAmerged.bed"
    echo ""
}

# Parse command line arguments
while getopts "i:o:b:vh" opt; do
    case $opt in
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        b) BED_FILE="$OPTARG" ;;
        v) VERBOSE="-v" ;;
        h) print_usage; exit 0 ;;
        *) print_usage; exit 1 ;;
    esac
done

# Check required arguments
if [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Error: Input directory (-i) and output directory (-o) are required"
    print_usage
    exit 1
fi

# Check if input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory '$INPUT_DIR' does not exist"
    exit 1
fi

# Check if narrowPeak files exist in input directory
if ! ls "$INPUT_DIR"/*.narrowPeak >/dev/null 2>&1; then
    echo "Error: No narrowPeak files found in '$INPUT_DIR'"
    exit 1
fi

# Check if bed file exists (if provided)
if [ -n "$BED_FILE" ] && [ ! -f "$BED_FILE" ]; then
    echo "Error: Bed file '$BED_FILE' does not exist"
    exit 1
fi

echo "Starting peak analysis..."
echo "Input directory: $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
if [ -n "$BED_FILE" ]; then
    echo "Bed file: $BED_FILE"
fi
echo ""

# Run the analysis
if [ -n "$BED_FILE" ]; then
    python src/analyze_peaks.py -i "$INPUT_DIR" -o "$OUTPUT_DIR" -b "$BED_FILE" $VERBOSE
else
    python src/analyze_peaks.py -i "$INPUT_DIR" -o "$OUTPUT_DIR" $VERBOSE
fi

echo ""
echo "Analysis completed successfully!"
echo "Results saved to: $OUTPUT_DIR" 