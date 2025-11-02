#!/usr/bin/env python3
"""
Extract Gene IDs from Annotation Results

This script extracts gene IDs from peak annotation CSV files and creates
gene list files for GO term analysis.

Date: 2025 (HS_ABA project)
"""

import pandas as pd
import argparse
from pathlib import Path
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def extract_genes_from_annotation(csv_file: str, output_file: str, gene_id_column: str = 'geneId'):
    """
    Extract gene IDs from annotation CSV file
    
    Args:
        csv_file: Path to annotation CSV file
        output_file: Path to output gene list file
        gene_id_column: Name of column containing gene IDs
    """
    try:
        # Read CSV file
        df = pd.read_csv(csv_file)
        
        # Check if gene ID column exists
        if gene_id_column not in df.columns:
            logger.error(f"Column '{gene_id_column}' not found in {csv_file}")
            logger.info(f"Available columns: {list(df.columns)}")
            return False
        
        # Extract unique gene IDs
        gene_ids = df[gene_id_column].dropna().unique()
        
        # Filter out empty strings and 'nan' values
        gene_ids = [gid for gid in gene_ids if str(gid).strip() and str(gid).lower() != 'nan']
        
        logger.info(f"Extracted {len(gene_ids)} unique gene IDs from {csv_file}")
        
        # Write to output file
        with open(output_file, 'w') as f:
            f.write("GeneID\n")  # Header
            for gene_id in sorted(gene_ids):
                f.write(f"{gene_id}\n")
        
        logger.info(f"Gene IDs saved to {output_file}")
        return True
        
    except Exception as e:
        logger.error(f"Error processing {csv_file}: {e}")
        return False

def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description="Extract gene IDs from peak annotation CSV files"
    )
    
    parser.add_argument(
        "--input-file", 
        required=True,
        help="Input annotation CSV file"
    )
    
    parser.add_argument(
        "--output-file", 
        required=True,
        help="Output gene list file"
    )
    
    parser.add_argument(
        "--gene-id-column",
        default="geneId",
        help="Name of column containing gene IDs (default: geneId)"
    )
    
    args = parser.parse_args()
    
    # Extract genes
    success = extract_genes_from_annotation(
        args.input_file, 
        args.output_file, 
        args.gene_id_column
    )
    
    if success:
        logger.info("Gene extraction completed successfully!")
    else:
        logger.error("Gene extraction failed!")
        exit(1)

if __name__ == "__main__":
    main()
