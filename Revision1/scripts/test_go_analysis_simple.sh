#!/bin/bash

# Simple GO Analysis Test Script
# This script tests the GO analysis with a small subset of data

echo "Testing GO Analysis with Real Data..."
echo "====================================="

# Set the project root directory
PROJECT_ROOT="/Users/user/Documents/HS_ABA/Revision1"
cd "$PROJECT_ROOT"

# Create a test directory with a small subset of data
TEST_DIR="test_go_analysis"
mkdir -p "$TEST_DIR"

echo "Creating test data directory: $TEST_DIR"

# Copy a small subset of data for testing
if [ -d "output/HS_DEG_hierarchical_clustering_4/cluster_results" ]; then
    echo "Copying test data from HS_DEG clusters..."
    cp output/HS_DEG_hierarchical_clustering_4/cluster_results/cluster_1_genes.txt "$TEST_DIR/"
    cp output/HS_DEG_hierarchical_clustering_4/cluster_results/cluster_2_genes.txt "$TEST_DIR/"
    echo "✓ Copied 2 cluster files"
else
    echo "Warning: HS_DEG cluster directory not found"
fi

# Check if we have test data
if [ -f "$TEST_DIR/cluster_1_genes.txt" ]; then
    echo ""
    echo "Running GO analysis on test data..."
    echo "=================================="
    
    # Create a simple test script
    cat > test_go_simple.py << 'EOF'
#!/usr/bin/env python3
"""
Simple test script for GO analysis
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from go_term_analysis import GOTermAnalyzer

def main():
    # Configuration
    ortholog_file = "data/Mp_Ath_Ortho/MpTak_v6.1r1.protein__v__Athaliana_447_Araport11.protein.csv"
    go_file = "data/Mp_Ath_Ortho/Athaliana_GO.txt"
    output_dir = "test_go_analysis/output"
    
    # Test with just one folder
    input_folders = ["test_go_analysis"]
    group_names = ["Test_Cluster"]
    
    try:
        # Initialize analyzer
        print("Initializing GO Term Analyzer...")
        analyzer = GOTermAnalyzer(ortholog_file, go_file, output_dir)
        
        # Perform analysis
        print("Starting GO term analysis...")
        results = analyzer.analyze_gene_lists(input_folders, group_names)
        
        # Show results
        if results:
            print(f"Analysis complete! Found results for {len(results)} groups")
            for group, df in results.items():
                if not df.empty:
                    print(f"  {group}: {len(df)} GO terms")
                else:
                    print(f"  {group}: No significant terms")
        else:
            print("No results generated")
            
    except Exception as e:
        print(f"Error during analysis: {e}")
        raise

if __name__ == "__main__":
    main()
EOF

    # Run the test
    python test_go_simple.py
    
    if [ $? -eq 0 ]; then
        echo ""
        echo "✓ Test completed successfully!"
        echo ""
        echo "Test results available in: test_go_analysis/output/"
        
        # Clean up test files
        echo "Cleaning up test files..."
        rm -rf "$TEST_DIR"
        rm test_go_simple.py
        
        echo "✓ Test cleanup completed"
    else
        echo ""
        echo "✗ Test failed!"
        exit 1
    fi
    
else
    echo "Error: No test data available"
    exit 1
fi

echo ""
echo "Test completed!"
