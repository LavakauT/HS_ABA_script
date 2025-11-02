#!/bin/bash

source activate featurecounts
ml  SAMtools/1.21-GCC-13.3.0
ml  BEDTools/2.31.0-GCC-12.3.0



# FRiP (Fraction of Reads in Peaks) Calculator
# This script calculates FRiP for all BAM files in a given directory
# against a merged BED file containing all peaks

# Don't exit on error, handle errors manually
# set -e  # Exit on any error

# Function to display usage
usage() {
    echo "Usage: $0 <bam_folder_path> <bed_file_path> [output_file]"
    echo ""
    echo "Arguments:"
    echo "  bam_folder_path  : Path to folder containing BAM files"
    echo "  bed_file_path    : Path to the merged BED file"
    echo "  output_file      : Output file path (optional, default: FRiP_results.txt)"
    echo ""
    echo "Example:"
    echo "  $0 /path/to/bam/files /path/to/merged_peaks.bed"
    echo "  $0 /path/to/bam/files /path/to/merged_peaks.bed my_FRiP_results.txt"
    exit 1
}

# Check if required tools are available
check_dependencies() {
    local missing_tools=()
    
    if ! command -v samtools &> /dev/null; then
        missing_tools+=("samtools")
    fi
    
    if ! command -v bedtools &> /dev/null; then
        missing_tools+=("bedtools")
    fi
    
    if [ ${#missing_tools[@]} -ne 0 ]; then
        echo "Error: Missing required tools: ${missing_tools[*]}"
        echo "Please install the missing tools and try again."
        exit 1
    fi
}

# Function to calculate FRiP for a single BAM file
calculate_FRiP() {
    local bam_file="$1"
    local bed_file="$2"
    local sample_name=$(basename "$bam_file" .bam)
    
    echo "  Processing: $sample_name" >&2
    
    # Check if BAM file is valid
    if ! samtools quickcheck "$bam_file" 2>/dev/null; then
        echo "  Warning: BAM file validation failed for $sample_name, but continuing..." >&2
        # Don't return 1, continue processing
    fi
    
    # Get total mapped reads (excluding unmapped reads)
    local total_mapped_reads
    total_mapped_reads=$(samtools view -c -F 4 "$bam_file" 2>/dev/null || echo "0")
    
    if [ -z "$total_mapped_reads" ] || [ "$total_mapped_reads" = "0" ]; then
        echo "  Warning: No mapped reads found for $sample_name" >&2
        total_mapped_reads=0
    fi
    
    # Get reads overlapping peaks
    local overlapping_reads
    overlapping_reads=$(bedtools intersect -a "$bam_file" -b "$bed_file" -u 2>/dev/null | samtools view -c 2>/dev/null || echo "0")
    
    if [ -z "$overlapping_reads" ] || [ "$overlapping_reads" = "0" ]; then
        echo "  Warning: No overlapping reads found for $sample_name" >&2
        overlapping_reads=0
    fi
    
    # Calculate FRiP
    local frip=0
    if [ "$total_mapped_reads" -gt 0 ]; then
        frip=$(echo "scale=6; $overlapping_reads / $total_mapped_reads" | bc -l 2>/dev/null || echo "0")
        if [ -z "$frip" ]; then
            frip=0
        fi
    fi
    
    echo "  Results: Total reads=$total_mapped_reads, Overlapping reads=$overlapping_reads, FRiP=$frip" >&2
    echo "$sample_name,$total_mapped_reads,$overlapping_reads,$frip"
}

# Main script
main() {
    # Check arguments
    if [ $# -lt 2 ] || [ $# -gt 3 ]; then
        usage
    fi
    
    local bam_folder="$1"
    local bed_file="$2"
    local output_file="${3:-FRiP_results.txt}"
    
    # Check if directories and files exist
    if [ ! -d "$bam_folder" ]; then
        echo "Error: BAM folder '$bam_folder' does not exist"
        exit 1
    fi
    
    if [ ! -f "$bed_file" ]; then
        echo "Error: BED file '$bed_file' does not exist"
        exit 1
    fi
    
    # Check dependencies
    check_dependencies
    
    echo "=== FRiP Calculator ==="
    echo "BAM folder: $bam_folder"
    echo "BED file: $bed_file"
    echo "Output file: $output_file"
    echo ""
    
    # Create output file with header
    echo "Sample,Total_Mapped_Reads,Reads_in_Peaks,FRiP" > "$output_file"
    
    # Find all BAM files in the directory
    echo "Searching for BAM files in: $bam_folder"
    
    # Use a more robust method to find BAM files
    local bam_files=()
    while IFS= read -r -d '' file; do
        bam_files+=("$file")
    done < <(find "$bam_folder" -name "*.bam" -type f -print0 2>/dev/null)
    
    if [ ${#bam_files[@]} -eq 0 ]; then
        echo "Error: No BAM files found in '$bam_folder'"
        echo "Please check:"
        echo "  1. The directory path is correct"
        echo "  2. BAM files have .bam extension"
        echo "  3. You have read permissions for the directory"
        echo ""
        echo "Files found in directory:"
        ls -la "$bam_folder" 2>/dev/null || echo "Cannot list directory contents"
        exit 1
    fi
    
    echo "Found ${#bam_files[@]} BAM file(s):"
    for bam_file in "${bam_files[@]}"; do
        echo "  - $(basename "$bam_file")"
    done
    echo ""
    
    # Process each BAM file
    local processed=0
    for bam_file in "${bam_files[@]}"; do
        if [ -f "$bam_file" ] && [ -r "$bam_file" ]; then
            echo "Processing file $((processed + 1)) of ${#bam_files[@]}: $(basename "$bam_file")"
            if calculate_FRiP "$bam_file" "$bed_file" >> "$output_file"; then
                ((processed++))
            else
                echo "Warning: Failed to process BAM file: $bam_file" >&2
            fi
        else
            echo "Warning: Cannot read BAM file: $bam_file"
        fi
    done
    
    echo ""
    echo "=== Summary ==="
    echo "Total BAM files found: ${#bam_files[@]}"
    echo "Successfully processed: $processed BAM file(s)"
    echo "Results saved to: $output_file"
    echo ""
    
    # Show a preview of the results
    if [ $processed -gt 0 ]; then
        echo "Preview of results:"
        echo "=================="
        head -5 "$output_file"
        echo ""
        echo "FRiP calculation completed successfully!"
    else
        echo "No BAM files were successfully processed."
        echo "Please check the error messages above."
        exit 1
    fi
}

# Run main function with all arguments
main "$@" 