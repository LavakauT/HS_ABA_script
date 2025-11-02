#!/usr/bin/env python3
"""
Peak Analysis and Visualization Script

This script analyzes narrowPeak files to count peaks and visualize their distribution.
It generates plots for peak counts and peak length distributions, with optional
reference line from a merged bed file.

"""

import os
import sys
import argparse
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle
import warnings

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')


class PeakAnalyzer:
    """
    A class to analyze narrowPeak files and generate visualizations.
    
    This class provides functionality to:
    - Count peaks in individual narrowPeak files
    - Analyze peak length distributions
    - Generate visualizations with optional reference lines
    """
    
    def __init__(self, input_dir: str, output_dir: str, bed_file: Optional[str] = None):
        """
        Initialize the PeakAnalyzer.
        
        Args:
            input_dir: Directory containing narrowPeak files
            output_dir: Directory to save output files
            bed_file: Optional merged bed file for reference line
        """
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.bed_file = Path(bed_file) if bed_file else None
        
        # Create output directory if it doesn't exist
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Set up plotting style
        plt.style.use('default')
        sns.set_palette("husl")
        
        logger.info(f"Initialized PeakAnalyzer with input_dir: {self.input_dir}")
        logger.info(f"Output directory: {self.output_dir}")
        if self.bed_file:
            logger.info(f"Reference bed file: {self.bed_file}")
    
    def find_narrowpeak_files(self) -> List[Path]:
        """
        Find all narrowPeak files in the input directory.
        
        Returns:
            List of paths to narrowPeak files
        """
        narrowpeak_files = list(self.input_dir.glob("*.narrowPeak"))
        
        if not narrowpeak_files:
            logger.error(f"No narrowPeak files found in {self.input_dir}")
            raise FileNotFoundError(f"No narrowPeak files found in {self.input_dir}")
        
        logger.info(f"Found {len(narrowpeak_files)} narrowPeak files")
        return narrowpeak_files
    
    def extract_sample_name(self, file_path: Path) -> str:
        """
        Extract sample name from narrowPeak file path.
        
        Args:
            file_path: Path to narrowPeak file
            
        Returns:
            Sample name (filename without _peaks.narrowPeak suffix)
        """
        filename = file_path.stem  # Remove .narrowPeak extension
        if filename.endswith('_peaks'):
            return filename[:-6]  # Remove '_peaks' suffix
        return filename
    
    def count_peaks_in_file(self, file_path: Path) -> int:
        """
        Count the number of peaks in a narrowPeak file.
        
        Args:
            file_path: Path to narrowPeak file
            
        Returns:
            Number of peaks in the file
        """
        try:
            with open(file_path, 'r') as f:
                # Count non-empty lines (each line represents one peak)
                peak_count = sum(1 for line in f if line.strip())
            return peak_count
        except Exception as e:
            logger.error(f"Error reading file {file_path}: {e}")
            return 0
    
    def get_peak_lengths(self, file_path: Path) -> List[int]:
        """
        Extract peak lengths from a narrowPeak file.
        
        Args:
            file_path: Path to narrowPeak file
            
        Returns:
            List of peak lengths
        """
        lengths = []
        try:
            with open(file_path, 'r') as f:
                for line in f:
                    if line.strip():
                        parts = line.strip().split('\t')
                        if len(parts) >= 3:
                            start = int(parts[1])
                            end = int(parts[2])
                            length = end - start
                            lengths.append(length)
        except Exception as e:
            logger.error(f"Error reading peak lengths from {file_path}: {e}")
        
        return lengths
    
    def count_peaks_in_bed_file(self, bed_file: Path) -> int:
        """
        Count the number of peaks in a bed file.
        
        Args:
            bed_file: Path to bed file
            
        Returns:
            Number of peaks in the bed file
        """
        try:
            with open(bed_file, 'r') as f:
                peak_count = sum(1 for line in f if line.strip())
            return peak_count
        except Exception as e:
            logger.error(f"Error reading bed file {bed_file}: {e}")
            return 0
    
    def analyze_peaks(self) -> Tuple[Dict[str, int], Dict[str, List[int]]]:
        """
        Analyze all narrowPeak files in the input directory.
        
        Returns:
            Tuple of (peak_counts, peak_lengths) dictionaries
        """
        narrowpeak_files = self.find_narrowpeak_files()
        
        peak_counts = {}
        peak_lengths = {}
        
        for file_path in narrowpeak_files:
            sample_name = self.extract_sample_name(file_path)
            peak_count = self.count_peaks_in_file(file_path)
            lengths = self.get_peak_lengths(file_path)
            
            peak_counts[sample_name] = peak_count
            peak_lengths[sample_name] = lengths
            
            if lengths:
                logger.info(f"Sample {sample_name}: {peak_count} peaks, "
                           f"length range: {min(lengths)}-{max(lengths)} bp")
            else:
                logger.info(f"Sample {sample_name}: {peak_count} peaks, no length data")
        
        return peak_counts, peak_lengths
    
    def create_peak_count_plot(self, peak_counts: Dict[str, int], 
                              bed_peak_count: Optional[int] = None) -> None:
        """
        Create a bar plot of peak counts.
        
        Args:
            peak_counts: Dictionary mapping sample names to peak counts
            bed_peak_count: Optional peak count from bed file for reference line
        """
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Sort samples by peak count for better visualization
        sorted_items = sorted(peak_counts.items(), key=lambda x: x[1], reverse=True)
        samples = [item[0] for item in sorted_items]
        counts = [item[1] for item in sorted_items]
        
        # Create bar plot
        bars = ax.bar(range(len(samples)), counts, 
                     color=sns.color_palette("husl", len(samples)))
        
        # Add reference line if bed file is provided
        if bed_peak_count is not None:
            ax.axhline(y=bed_peak_count, color='red', linestyle='--', 
                      linewidth=2, label=f'Merged peaks ({bed_peak_count})')
            ax.legend()
        
        # Customize plot
        ax.set_xlabel('Samples', fontsize=12, fontweight='bold')
        ax.set_ylabel('Number of Peaks', fontsize=12, fontweight='bold')
        ax.set_title('Peak Counts by Sample', fontsize=14, fontweight='bold')
        ax.set_xticks(range(len(samples)))
        ax.set_xticklabels(samples, rotation=45, ha='right')
        
        # Add value labels on bars
        for i, (bar, count) in enumerate(zip(bars, counts)):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + max(counts)*0.01,
                   f'{count:,}', ha='center', va='bottom', fontweight='bold')
        
        # Adjust layout
        plt.tight_layout()
        
        # Save plot
        output_file = self.output_dir / "peak_counts.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"Peak count plot saved to {output_file}")
        
        # Also save as PDF
        pdf_file = self.output_dir / "peak_counts.pdf"
        plt.savefig(pdf_file, bbox_inches='tight')
        logger.info(f"Peak count plot saved to {pdf_file}")
        
        plt.close()
    
    def create_length_distribution_plot(self, peak_lengths: Dict[str, List[int]]) -> None:
        """
        Create a violin plot of peak length distributions.
        
        Args:
            peak_lengths: Dictionary mapping sample names to lists of peak lengths
        """
        # Prepare data for plotting
        plot_data = []
        sample_names = []
        
        for sample_name, lengths in peak_lengths.items():
            if lengths:  # Only include samples with peaks
                plot_data.extend(lengths)
                sample_names.extend([sample_name] * len(lengths))
        
        if not plot_data:
            logger.warning("No peak length data available for plotting")
            return
        
        # Create DataFrame for easier plotting
        df = pd.DataFrame({
            'Sample': sample_names,
            'Peak Length (bp)': plot_data
        })
        
        # Create violin plot
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
        
        # Violin plot
        sns.violinplot(data=df, x='Sample', y='Peak Length (bp)', ax=ax1)
        ax1.set_title('Peak Length Distribution by Sample', fontsize=14, fontweight='bold')
        ax1.set_xlabel('Samples', fontsize=12, fontweight='bold')
        ax1.set_ylabel('Peak Length (bp)', fontsize=12, fontweight='bold')
        ax1.tick_params(axis='x', rotation=45)
        
        # Box plot
        sns.boxplot(data=df, x='Sample', y='Peak Length (bp)', ax=ax2)
        ax2.set_title('Peak Length Distribution (Box Plot)', fontsize=14, fontweight='bold')
        ax2.set_xlabel('Samples', fontsize=12, fontweight='bold')
        ax2.set_ylabel('Peak Length (bp)', fontsize=12, fontweight='bold')
        ax2.tick_params(axis='x', rotation=45)
        
        plt.tight_layout()
        
        # Save plots
        output_file = self.output_dir / "peak_length_distribution.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"Peak length distribution plot saved to {output_file}")
        
        # Also save as PDF
        pdf_file = self.output_dir / "peak_length_distribution.pdf"
        plt.savefig(pdf_file, bbox_inches='tight')
        logger.info(f"Peak length distribution plot saved to {pdf_file}")
        
        plt.close()
    
    def create_summary_statistics(self, peak_counts: Dict[str, int], 
                                peak_lengths: Dict[str, List[int]]) -> None:
        """
        Create summary statistics and save to CSV.
        
        Args:
            peak_counts: Dictionary mapping sample names to peak counts
            peak_lengths: Dictionary mapping sample names to lists of peak lengths
        """
        summary_data = []
        
        for sample_name in peak_counts.keys():
            count = peak_counts[sample_name]
            lengths = peak_lengths[sample_name]
            
            if lengths:
                summary_data.append({
                    'Sample': sample_name,
                    'Peak_Count': count,
                    'Mean_Length': np.mean(lengths),
                    'Median_Length': np.median(lengths),
                    'Std_Length': np.std(lengths),
                    'Min_Length': np.min(lengths),
                    'Max_Length': np.max(lengths),
                    'Q25_Length': np.percentile(lengths, 25),
                    'Q75_Length': np.percentile(lengths, 75)
                })
            else:
                summary_data.append({
                    'Sample': sample_name,
                    'Peak_Count': count,
                    'Mean_Length': 0,
                    'Median_Length': 0,
                    'Std_Length': 0,
                    'Min_Length': 0,
                    'Max_Length': 0,
                    'Q25_Length': 0,
                    'Q75_Length': 0
                })
        
        # Create DataFrame and save
        df = pd.DataFrame(summary_data)
        output_file = self.output_dir / "peak_summary_statistics.csv"
        df.to_csv(output_file, index=False)
        logger.info(f"Summary statistics saved to {output_file}")
        
        # Also save as Excel file
        excel_file = self.output_dir / "peak_summary_statistics.xlsx"
        df.to_excel(excel_file, index=False)
        logger.info(f"Summary statistics saved to {excel_file}")
        
        return df
    
    def run_analysis(self) -> None:
        """
        Run the complete peak analysis pipeline.
        """
        logger.info("Starting peak analysis...")
        
        # Analyze peaks
        peak_counts, peak_lengths = self.analyze_peaks()
        
        # Get bed file peak count if provided
        bed_peak_count = None
        if self.bed_file and self.bed_file.exists():
            bed_peak_count = self.count_peaks_in_bed_file(self.bed_file)
            logger.info(f"Merged bed file contains {bed_peak_count} peaks")
        
        # Create visualizations
        self.create_peak_count_plot(peak_counts, bed_peak_count)
        self.create_length_distribution_plot(peak_lengths)
        
        # Create summary statistics
        summary_df = self.create_summary_statistics(peak_counts, peak_lengths)
        
        # Print summary
        logger.info("\n" + "="*50)
        logger.info("PEAK ANALYSIS SUMMARY")
        logger.info("="*50)
        for sample, count in sorted(peak_counts.items()):
            logger.info(f"{sample}: {count:,} peaks")
        
        if bed_peak_count:
            logger.info(f"Merged peaks: {bed_peak_count:,}")
        
        logger.info("="*50)
        logger.info("Analysis completed successfully!")


def main():
    """
    Main function to run the peak analysis.
    """
    parser = argparse.ArgumentParser(
        description="Analyze narrowPeak files and generate visualizations",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python analyze_peaks.py -i data/narrowpeak/ -o output/peaks/
  python analyze_peaks.py -i data/narrowpeak/ -o output/peaks/ -b merged_peaks.bed
        """
    )
    
    parser.add_argument(
        '-i', '--input_dir',
        required=True,
        help='Directory containing narrowPeak files'
    )
    
    parser.add_argument(
        '-o', '--output_dir',
        required=True,
        help='Output directory for results'
    )
    
    parser.add_argument(
        '-b', '--bed_file',
        help='Optional merged bed file for reference line'
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose logging'
    )
    
    args = parser.parse_args()
    
    # Set logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    try:
        # Create analyzer and run analysis
        analyzer = PeakAnalyzer(args.input_dir, args.output_dir, args.bed_file)
        analyzer.run_analysis()
        
    except Exception as e:
        logger.error(f"Error during analysis: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main() 