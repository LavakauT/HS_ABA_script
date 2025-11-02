#!/bin/bash

# Debug script to help identify issues with BAM file processing
# This script will examine a directory and show what BAM files are found

set -e

echo "=== BAM File Debugging Tool ==="

if [ $# -eq 0 ]; then
    echo "Usage: $0 <bam_folder_path>"
    echo ""
    echo "This script will examine the specified directory and show:"
    echo "  1. All files in the directory"
    echo "  2. BAM files found"
    echo "  3. File permissions and sizes"
    echo "  4. Test if files are valid BAM files"
    exit 1
fi

BAM_FOLDER="$1"

echo "Examining directory: $BAM_FOLDER"
echo ""

# Check if directory exists
if [ ! -d "$BAM_FOLDER" ]; then
    echo "❌ Error: Directory '$BAM_FOLDER' does not exist"
    exit 1
fi

# Show all files in directory
echo "📁 All files in directory:"
ls -la "$BAM_FOLDER"
echo ""

# Find BAM files using different methods
echo "🔍 Searching for BAM files..."

# Method 1: Simple find
echo "Method 1 - Simple find:"
find "$BAM_FOLDER" -name "*.bam" -type f 2>/dev/null | while read -r file; do
    echo "  Found: $file"
done

# Method 2: Find with null termination
echo ""
echo "Method 2 - Find with null termination:"
find "$BAM_FOLDER" -name "*.bam" -type f -print0 2>/dev/null | while IFS= read -r -d '' file; do
    echo "  Found: $file"
done

# Method 3: Using ls
echo ""
echo "Method 3 - Using ls:"
ls "$BAM_FOLDER"/*.bam 2>/dev/null || echo "  No .bam files found with ls"

# Count BAM files
BAM_COUNT=$(find "$BAM_FOLDER" -name "*.bam" -type f 2>/dev/null | wc -l)
echo ""
echo "📊 Summary:"
echo "  Total BAM files found: $BAM_COUNT"

# Test each BAM file
if [ "$BAM_COUNT" -gt 0 ]; then
    echo ""
    echo "🔬 Testing BAM files:"
    find "$BAM_FOLDER" -name "*.bam" -type f 2>/dev/null | while read -r bam_file; do
        echo "  Testing: $(basename "$bam_file")"
        
        # Check file permissions
        if [ -r "$bam_file" ]; then
            echo "    ✅ Readable"
        else
            echo "    ❌ Not readable"
        fi
        
        # Check file size
        SIZE=$(stat -f%z "$bam_file" 2>/dev/null || stat -c%s "$bam_file" 2>/dev/null || echo "unknown")
        echo "    📏 Size: $SIZE bytes"
        
        # Test if it's a valid BAM file (if samtools is available)
        if command -v samtools &> /dev/null; then
            if samtools quickcheck "$bam_file" 2>/dev/null; then
                echo "    ✅ Valid BAM file"
                
                # Get basic stats
                TOTAL_READS=$(samtools view -c "$bam_file" 2>/dev/null || echo "ERROR")
                MAPPED_READS=$(samtools view -c -F 4 "$bam_file" 2>/dev/null || echo "ERROR")
                echo "    📊 Total reads: $TOTAL_READS"
                echo "    📊 Mapped reads: $MAPPED_READS"
            else
                echo "    ❌ Invalid or corrupted BAM file"
            fi
        else
            echo "    ⚠️  samtools not available for validation"
        fi
        echo ""
    done
else
    echo ""
    echo "❌ No BAM files found in directory"
    echo ""
    echo "Possible issues:"
    echo "  1. Files don't have .bam extension"
    echo "  2. Files are not regular files (might be symlinks)"
    echo "  3. Directory permissions issue"
    echo "  4. Files are in subdirectories"
    echo ""
    echo "Files with similar extensions:"
    find "$BAM_FOLDER" -name "*bam*" -o -name "*BAM*" 2>/dev/null | head -10
fi

echo ""
echo "=== Debug Complete ===" 