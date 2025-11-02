#!/bin/bash

# Test script to verify FRiP calculator processes multiple files
# This script creates test BAM files and runs the FRiP calculation

set -e

echo "=== Testing FRiP Calculator with Multiple Files ==="

# Create temporary directory
TEMP_DIR=$(mktemp -d)
echo "Created temporary directory: $TEMP_DIR"

# Create test BED file
cat > "$TEMP_DIR/test_peaks.bed" << 'EOF'
chr1	1000	2000	peak1	100	.
chr1	3000	4000	peak2	100	.
chr1	5000	6000	peak3	100	.
chr2	1000	2000	peak4	100	.
chr2	3000	4000	peak5	100	.
EOF

echo "Created test BED file with 5 peaks"

# Create test BAM files (simulated)
# Note: These are not real BAM files, just for testing the script logic
for i in {1..5}; do
    cat > "$TEMP_DIR/sample${i}.bam" << EOF
# This is a fake BAM file for testing
# In real usage, these would be actual BAM files
EOF
    echo "Created test BAM file: sample${i}.bam"
done

# Test the script
echo ""
echo "Testing FRiP calculator..."
echo "BAM directory: $TEMP_DIR"
echo "BED file: $TEMP_DIR/test_peaks.bed"

# Run the FRiP calculator
if ./scripts/calculate_FRiP.sh "$TEMP_DIR" "$TEMP_DIR/test_peaks.bed" "$TEMP_DIR/test_results.txt"; then
    echo ""
    echo "✅ Script executed successfully!"
    echo ""
    echo "Results file contents:"
    cat "$TEMP_DIR/test_results.txt"
else
    echo ""
    echo "❌ Script failed!"
fi

# Clean up
echo ""
echo "Cleaning up temporary files..."
rm -rf "$TEMP_DIR"
echo "Test completed!" 