#!/usr/bin/env python3
"""
Extract Genes for Peak Files from HS_DAR_Annotation

This script reads peak files (_intersection.txt or _only.txt) and extracts
corresponding genes from HS_DAR_annotation files, avoiding duplicates.

Date: 2025 (HS_ABA project)
"""

import pandas as pd
import argparse
from pathlib import Path
import logging
import re
from typing import Dict, Set, List, Tuple

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def load_peaks_annotation(annotation_dir: Path) -> Dict[str, Set[str]]:
    """
    Load all peaks annotation files and create a mapping from peak coordinates to gene IDs
    
    Args:
        annotation_dir: Directory containing peaks annotation files
        
    Returns:
        Dictionary mapping peak coordinates to set of gene IDs
    """
    peak_to_genes: Dict[str, Set[str]] = {}
    
    # Find all peaks annotation files
    annotation_files = list(annotation_dir.glob("*_peaks_annotation_results.csv"))
    
    for file_path in annotation_files:
        logger.info(f"Loading annotation file: {file_path.name}")
        try:
            df = pd.read_csv(file_path)
            
            # Check if required columns exist
            if 'Peak' not in df.columns or 'geneId' not in df.columns:
                logger.warning(f"Skipping {file_path.name} - missing required columns")
                continue
            
            # Process each peak
            for _, row in df.iterrows():
                peak_coord = str(row['Peak']).strip()
                gene_id = str(row['geneId']).strip()
                
                # Skip invalid entries
                if pd.isna(peak_coord) or pd.isna(gene_id) or peak_coord == 'nan' or gene_id == 'nan':
                    continue
                
                # Initialize set if peak doesn't exist
                if peak_coord not in peak_to_genes:
                    peak_to_genes[peak_coord] = set()
                
                # Add gene ID to the set
                peak_to_genes[peak_coord].add(gene_id)
                
        except Exception as e:
            logger.error(f"Error loading {file_path}: {e}")
            continue
    
    logger.info(f"Loaded {len(peak_to_genes)} unique peaks with gene annotations")
    return peak_to_genes

def extract_genes_for_peak_file(peak_file: Path, peak_to_genes: Dict[str, Set[str]]) -> Set[str]:
    """
    Extract genes for a specific peak file
    
    Args:
        peak_file: Path to peak file
        peak_to_genes: Mapping from peak coordinates to gene IDs
        
    Returns:
        Set of unique gene IDs
    """
    genes: Set[str] = set()
    
    try:
        with open(peak_file, 'r') as f:
            for line in f:
                line = line.strip()
                
                # Skip comments and empty lines
                if line.startswith('#') or not line:
                    continue
                
                # Extract peak coordinates (format: chr1.10215491.10216093)
                if re.match(r'chr\d+\.\d+\.\d+', line):
                    if line in peak_to_genes:
                        genes.update(peak_to_genes[line])
                    else:
                        logger.debug(f"Peak {line} not found in annotation data")
                        
    except Exception as e:
        logger.error(f"Error reading {peak_file}: {e}")
    
    return genes

def process_peak_files(up_dir: Path, down_dir: Path, annotation_dir: Path, output_dir: Path):
    """
    Process all peak files and extract corresponding genes
    
    Args:
        up_dir: Directory containing up-regulated peak files
        down_dir: Directory containing down-regulated peak files
        annotation_dir: Directory containing peaks annotation files
        output_dir: Output directory for gene files
    """
    # Create output directories
    up_output_dir = output_dir / "up_genes"
    down_output_dir = output_dir / "down_genes"
    
    up_output_dir.mkdir(parents=True, exist_ok=True)
    down_output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load peaks annotation data
    logger.info("Loading peaks annotation data...")
    peak_to_genes = load_peaks_annotation(annotation_dir)
    
    # Process up-regulated peak files
    logger.info("Processing up-regulated peak files...")
    up_files = list(up_dir.glob("*_intersection.txt")) + list(up_dir.glob("*_only.txt"))
    
    for peak_file in up_files:
        logger.info(f"Processing up file: {peak_file.name}")
        genes = extract_genes_for_peak_file(peak_file, peak_to_genes)
        
        if genes:
            # Create output filename (replace .txt with _genes.txt)
            output_filename = peak_file.stem + "_genes.txt"
            output_path = up_output_dir / output_filename
            
            # Write genes to file
            with open(output_path, 'w') as f:
                f.write(f"# Genes extracted from {peak_file.name}\n")
                f.write(f"# Total genes: {len(genes)}\n")
                f.write("GeneID\n")
                for gene_id in sorted(genes):
                    f.write(f"{gene_id}\n")
            
            logger.info(f"Saved {len(genes)} genes to {output_path}")
        else:
            logger.warning(f"No genes found for {peak_file.name}")
    
    # Process down-regulated peak files
    logger.info("Processing down-regulated peak files...")
    down_files = list(down_dir.glob("*_intersection.txt")) + list(down_dir.glob("*_only.txt"))
    
    for peak_file in down_files:
        logger.info(f"Processing down file: {peak_file.name}")
        genes = extract_genes_for_peak_file(peak_file, peak_to_genes)
        
        if genes:
            # Create output filename (replace .txt with _genes.txt)
            output_filename = peak_file.stem + "_genes.txt"
            output_path = down_output_dir / output_filename
            
            # Write genes to file
            with open(output_path, 'w') as f:
                f.write(f"# Genes extracted from {peak_file.name}\n")
                f.write(f"# Total genes: {len(genes)}\n")
                f.write("GeneID\n")
                for gene_id in sorted(genes):
                    f.write(f"{gene_id}\n")
            
            logger.info(f"Saved {len(genes)} genes to {output_path}")
        else:
            logger.warning(f"No genes found for {peak_file.name}")

def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description="Extract genes for peak files from HS_DAR_annotation"
    )
    
    parser.add_argument(
        "--up-dir",
        default="output/HS_DAR_analysis/up",
        help="Directory containing up-regulated peak files (default: output/HS_DAR_analysis/up)"
    )
    
    parser.add_argument(
        "--down-dir",
        default="output/HS_DAR_analysis/down",
        help="Directory containing down-regulated peak files (default: output/HS_DAR_analysis/down)"
    )
    
    parser.add_argument(
        "--annotation-dir",
        default="output/HS_DAR_annotation",
        help="Directory containing peaks annotation files (default: output/HS_DAR_annotation)"
    )
    
    parser.add_argument(
        "--output-dir",
        default="output/HS_DAR_genes_extracted",
        help="Output directory for extracted gene files (default: output/HS_DAR_genes_extracted)"
    )
    
    args = parser.parse_args()
    
    # Convert to Path objects
    up_dir = Path(args.up_dir)
    down_dir = Path(args.down_dir)
    annotation_dir = Path(args.annotation_dir)
    output_dir = Path(args.output_dir)
    
    # Check if directories exist
    if not up_dir.exists():
        logger.error(f"Up directory does not exist: {up_dir}")
        exit(1)
    
    if not down_dir.exists():
        logger.error(f"Down directory does not exist: {down_dir}")
        exit(1)
    
    if not annotation_dir.exists():
        logger.error(f"Annotation directory does not exist: {annotation_dir}")
        exit(1)
    
    # Process peak files
    process_peak_files(up_dir, down_dir, annotation_dir, output_dir)
    
    logger.info("Gene extraction completed successfully!")

if __name__ == "__main__":
    main()
