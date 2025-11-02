#!/bin/bash

# Script to run samtools flagstat on all BAM files in a directory
# and output results as a formatted table with column names and sample information

# Check if directory path is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <directory_path>"
    echo "Example: $0 /path/to/bam/files"
    echo "Default: Using /Users/user/Documents/HS_ABA/f1/bam"
    BAM_DIR="/Users/user/Documents/HS_ABA/f1/bam"
else
    BAM_DIR="$1"
fi

OUTPUT_FILE="flagstat_results.txt"

# Check if directory exists
if [ ! -d "$BAM_DIR" ]; then
    echo "Error: Directory '$BAM_DIR' does not exist"
    exit 1
fi

# Check if samtools is available
if ! command -v samtools &> /dev/null; then
    echo "Warning: samtools is not installed or not in PATH"
    echo "You may need to load a module or install samtools"
    echo "Try: module load SAMtools/1.21-GCC-13.3.0"
    exit 1
fi

echo "Running samtools flagstat on BAM files in: $BAM_DIR"
echo "Output will be saved to: $OUTPUT_FILE"
echo ""

# Create header for the output file
cat > "$OUTPUT_FILE" << EOF
# Samtools Flagstat Results
# Generated on: $(date)
# Directory: $BAM_DIR
#
# Columns:
# Sample: Sample name (derived from BAM filename)
# Total_Reads: Total number of reads
# QC_Passed: Number of reads that passed QC
# QC_Failed: Number of reads that failed QC
# Secondary: Number of secondary alignments
# Supplementary: Number of supplementary alignments
# Duplicates: Number of duplicate reads
# Mapped: Number of mapped reads
# Paired: Number of paired reads
# Read1: Number of first reads in pair
# Read2: Number of second reads in pair
# Properly_Paired: Number of properly paired reads
# Both_Mapped: Number of reads where both mates mapped
# Singletons: Number of singleton reads
# Mate_Other_Chr: Number of reads with mate mapped to different chromosome
# Mate_Other_Chr_Qual: Number of reads with mate mapped to different chromosome with mapping quality >= 5
#
Sample	Total_Reads	QC_Passed	QC_Failed	Secondary	Supplementary	Duplicates	Mapped	Paired	Read1	Read2	Properly_Paired	Both_Mapped	Singletons	Mate_Other_Chr	Mate_Other_Chr_Qual
EOF

# Counter for processed files
processed=0
errors=0

# Process each BAM file
for bam_file in "$BAM_DIR"/*.bam; do
    # Check if file exists (in case no .bam files found)
    if [ ! -f "$bam_file" ]; then
        echo "No BAM files found in $BAM_DIR"
        exit 1
    fi
    
    # Extract sample name from filename (remove .bam extension and path)
    sample_name=$(basename "$bam_file" .bam)
    
    echo "Processing: $sample_name"
    
    # Run samtools flagstat and capture output
    if samtools flagstat "$bam_file" > /tmp/flagstat_temp.txt 2>/dev/null; then
        # Parse the flagstat output and extract values
        # Total reads and QC status
        total_line=$(grep "in total" /tmp/flagstat_temp.txt)
        total_reads=$(echo "$total_line" | awk '{print $1}')
        qc_passed=$(echo "$total_line" | awk '{print $1}')
        qc_failed=$(echo "$total_line" | awk '{print $3}')
        
        # Secondary alignments
        secondary=$(grep "secondary" /tmp/flagstat_temp.txt | awk '{print $1}')
        
        # Supplementary alignments
        supplementary=$(grep "supplementary" /tmp/flagstat_temp.txt | awk '{print $1}')
        
        # Duplicates
        duplicates=$(grep "duplicates" /tmp/flagstat_temp.txt | head -1 | awk '{print $1}')
        
        # Mapped reads
        mapped=$(grep "mapped" /tmp/flagstat_temp.txt | head -1 | awk '{print $1}')
        
        # Paired reads
        paired=$(grep "paired in sequencing" /tmp/flagstat_temp.txt | awk '{print $1}')
        
        # Read1 and Read2
        read1=$(grep "read1" /tmp/flagstat_temp.txt | awk '{print $1}')
        read2=$(grep "read2" /tmp/flagstat_temp.txt | awk '{print $1}')
        
        # Properly paired
        properly_paired=$(grep "properly paired" /tmp/flagstat_temp.txt | awk '{print $1}')
        
        # Both mapped
        both_mapped=$(grep "with itself and mate mapped" /tmp/flagstat_temp.txt | awk '{print $1}')
        
        # Singletons
        singletons=$(grep "singletons" /tmp/flagstat_temp.txt | awk '{print $1}')
        
        # Mate mapped to different chromosome
        mate_other_chr=$(grep "with mate mapped to a different chr" /tmp/flagstat_temp.txt | head -1 | awk '{print $1}')
        
        # Mate mapped to different chromosome with quality >= 5
        mate_other_chr_qual=$(grep "with mate mapped to a different chr (mapQ>=5)" /tmp/flagstat_temp.txt | awk '{print $1}')
        
        # Write results to output file
        echo -e "$sample_name\t$total_reads\t$qc_passed\t$qc_failed\t$secondary\t$supplementary\t$duplicates\t$mapped\t$paired\t$read1\t$read2\t$properly_paired\t$both_mapped\t$singletons\t$mate_other_chr\t$mate_other_chr_qual" >> "$OUTPUT_FILE"
        
        ((processed++))
        echo "  ✓ Completed"
    else
        echo "  ✗ Error processing $sample_name"
        ((errors++))
    fi
done

# Clean up temporary file
rm -f /tmp/flagstat_temp.txt

echo ""
echo "Processing completed!"
echo "Files processed: $processed"
echo "Errors encountered: $errors"
echo "Results saved to: $OUTPUT_FILE"

# Display summary
if [ $processed -gt 0 ]; then
    echo ""
    echo "Summary of results:"
    echo "==================="
    tail -n +18 "$OUTPUT_FILE" | head -5
    if [ $processed -gt 5 ]; then
        echo "..."
        echo "(showing first 5 samples, see $OUTPUT_FILE for complete results)"
    fi
fi 