#!/usr/bin/env python3
"""
STAR Mapping Results Analysis Script

This script analyzes STAR mapping results from Log.final.out files and generates:
1. A comprehensive CSV report with all mapping statistics
2. PDF visualizations of key mapping metrics

"""

import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import logging
from datetime import datetime

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class STARMappingAnalyzer:
    """
    Analyzer for STAR mapping results from Log.final.out files.
    
    This class provides functionality to:
    - Parse STAR Log.final.out files
    - Extract mapping statistics
    - Generate CSV reports
    - Create PDF visualizations
    """
    
    def __init__(self, input_dir: str, output_dir: str):
        """
        Initialize the STAR mapping analyzer.
        
        Args:
            input_dir: Directory containing Log.final.out files
            output_dir: Directory to save output files
        """
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Define fields to exclude from CSV output
        self.excluded_fields = {
            'Started job on',
            'Started mapping on', 
            'Finished on',
            'Mapping speed, Million of reads per hour'
        }
        
        # Define fields for visualization
        self.visualization_fields = [
            'Number of input reads',
            'Average input read length',
            'Uniquely mapped reads number',
            'Uniquely mapped reads %',
            'Average mapped length',
            'Number of reads unmapped: too many mismatches',
            '% of reads unmapped: too many mismatches',
            'Number of reads unmapped: too short',
            '% of reads unmapped: too short',
            'Number of reads unmapped: other',
            '% of reads unmapped: other'
        ]
    
    def find_log_files(self) -> List[Path]:
        """
        Find all Log.final.out files in the input directory.
        
        Returns:
            List of Path objects pointing to Log.final.out files
        """
        log_files = list(self.input_dir.rglob("*_Log.final.out"))
        logger.info(f"Found {len(log_files)} Log.final.out files")
        return log_files
    
    def extract_sample_name(self, file_path: Path) -> str:
        """
        Extract sample name from file path by removing '_Log.final.out' suffix.
        
        Args:
            file_path: Path to the Log.final.out file
            
        Returns:
            Sample name without the suffix
        """
        filename = file_path.name
        if filename.endswith('_Log.final.out'):
            return filename[:-len('_Log.final.out')]
        return file_path.stem
    
    def parse_log_file(self, file_path: Path) -> Dict[str, str]:
        """
        Parse a single Log.final.out file and extract mapping statistics.
        
        Args:
            file_path: Path to the Log.final.out file
            
        Returns:
            Dictionary containing mapping statistics
        """
        stats = {}
        
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if '|' in line:
                        # Split by | and clean up whitespace
                        parts = [part.strip() for part in line.split('|')]
                        if len(parts) >= 2:
                            field_name = parts[0].strip()
                            field_value = parts[1].strip()
                            
                            # Skip excluded fields
                            if field_name not in self.excluded_fields:
                                stats[field_name] = field_value
        
        except Exception as e:
            logger.error(f"Error parsing file {file_path}: {e}")
            return {}
        
        return stats
    
    def convert_value_to_numeric(self, value: str) -> float:
        """
        Convert string value to numeric, handling percentages and special formats.
        
        Args:
            value: String value from log file
            
        Returns:
            Numeric value
        """
        if not value or value == '':
            return 0.0
        
        # Remove percentage signs and convert to float
        if '%' in value:
            return float(value.replace('%', ''))
        
        # Handle comma-separated numbers
        if ',' in value:
            return float(value.replace(',', ''))
        
        # Try to convert to float
        try:
            return float(value)
        except ValueError:
            return 0.0
    
    def process_all_files(self) -> pd.DataFrame:
        """
        Process all Log.final.out files and create a DataFrame.
        
        Returns:
            DataFrame containing all mapping statistics
        """
        log_files = self.find_log_files()
        all_data = []
        
        for file_path in log_files:
            sample_name = self.extract_sample_name(file_path)
            stats = self.parse_log_file(file_path)
            
            if stats:
                stats['Sample'] = sample_name
                all_data.append(stats)
                logger.info(f"Processed {sample_name}")
            else:
                logger.warning(f"Failed to parse {file_path}")
        
        if not all_data:
            logger.error("No data was extracted from log files")
            return pd.DataFrame()
        
        # Create DataFrame
        df = pd.DataFrame(all_data)
        
        # Convert numeric columns
        for col in df.columns:
            if col != 'Sample':
                df[col] = df[col].apply(self.convert_value_to_numeric)
        
        return df
    
    def save_csv_report(self, df: pd.DataFrame) -> str:
        """
        Save the mapping statistics to a CSV file.
        
        Args:
            df: DataFrame containing mapping statistics
            
        Returns:
            Path to the saved CSV file
        """
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        csv_path = self.output_dir / f"star_mapping_summary_{timestamp}.csv"
        
        df.to_csv(csv_path, index=False)
        logger.info(f"CSV report saved to {csv_path}")
        
        return str(csv_path)
    
    def create_visualizations(self, df: pd.DataFrame) -> List[str]:
        """
        Create individual PDF visualizations for each key mapping metric.
        
        Args:
            df: DataFrame containing mapping statistics
            
        Returns:
            List of paths to the saved PDF files
        """
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        pdf_paths = []
        
        # Set up the plotting style
        plt.style.use('default')
        sns.set_palette("husl")
        
        for field in self.visualization_fields:
            if field in df.columns:
                # Create individual figure for each metric
                fig, ax = plt.subplots(figsize=(12, 8))
                
                # Create bar plot
                bars = ax.bar(range(len(df)), df[field], color='skyblue', alpha=0.7)
                ax.set_title(f'{field}', fontsize=14, fontweight='bold', pad=20)
                ax.set_xlabel('Samples', fontsize=12)
                ax.set_ylabel('Value', fontsize=12)
                
                # Rotate x-axis labels for better readability
                ax.set_xticks(range(len(df)))
                ax.set_xticklabels(df['Sample'], rotation=45, ha='right', fontsize=10)
                
                # Remove value labels on bars to avoid overlapping with many samples
                
                # Adjust layout
                plt.tight_layout()
                
                # Save individual PDF
                safe_field_name = field.replace(' ', '_').replace(',', '').replace('%', 'percent').replace(':', '')
                pdf_path = self.output_dir / f"star_mapping_{safe_field_name}_{timestamp}.pdf"
                plt.savefig(pdf_path, dpi=300, bbox_inches='tight')
                plt.close()
                
                pdf_paths.append(str(pdf_path))
                logger.info(f"Visualization for {field} saved to {pdf_path}")
        
        return pdf_paths
    
    def run_analysis(self) -> Tuple[str, List[str]]:
        """
        Run the complete STAR mapping analysis.
        
        Returns:
            Tuple of (csv_path, pdf_paths)
        """
        logger.info("Starting STAR mapping analysis...")
        
        # Process all files
        df = self.process_all_files()
        
        if df.empty:
            logger.error("No data to analyze")
            return "", []
        
        # Save CSV report
        csv_path = self.save_csv_report(df)
        
        # Create visualizations
        pdf_paths = self.create_visualizations(df)
        
        logger.info("STAR mapping analysis completed successfully")
        return csv_path, pdf_paths


def main():
    """
    Main function to run the STAR mapping analysis.
    """
    import argparse
    
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description="STAR Mapping Analysis Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python src/analyze_star_mapping.py
  python src/analyze_star_mapping.py --input data/RNA_STAR --output output/my_analysis
  python src/analyze_star_mapping.py -i data/RNA_STAR -o output/my_analysis
        """
    )
    
    parser.add_argument(
        "--input", "-i",
        type=str,
        default="data/RNA_STAR",
        help="Input directory containing Log.final.out files (default: data/RNA_STAR)"
    )
    
    parser.add_argument(
        "--output", "-o",
        type=str,
        default="output/STAR_mapping_analysis",
        help="Output directory for results (default: output/STAR_mapping_analysis)"
    )
    
    # Parse arguments
    args = parser.parse_args()
    
    # Validate input directory
    if not os.path.exists(args.input):
        print(f"Error: Input directory '{args.input}' does not exist.")
        return
    
    # Create analyzer and run analysis
    analyzer = STARMappingAnalyzer(args.input, args.output)
    csv_path, pdf_paths = analyzer.run_analysis()
    
    if csv_path and pdf_paths:
        print(f"Analysis completed successfully!")
        print(f"CSV report: {csv_path}")
        print(f"PDF visualizations ({len(pdf_paths)} files):")
        for pdf_path in pdf_paths:
            print(f"  - {pdf_path}")
    else:
        print("Analysis failed. Please check the logs for details.")


if __name__ == "__main__":
    main() 