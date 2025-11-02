#!/usr/bin/env python3
"""
FRiP (Fraction of Reads in Peaks) Visualization Tool

This script creates comprehensive visualizations for FRiP analysis results,
including bar plots, box plots, correlation analysis, and statistical summaries.

Author: AI Assistant
Date: 2024
"""

import argparse
import logging
import os
import sys
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.figure import Figure
from scipy import stats

# Suppress warnings
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Set style for better plots
plt.style.use('seaborn-v0_8')
sns.set_palette("husl")


class FRiPVisualizer:
    """
    A comprehensive FRiP visualization tool that creates various plots
    and statistical analyses from FRiP calculation results.
    """
    
    def __init__(self, input_file: str, output_dir: str):
        """
        Initialize the FRiP visualizer.
        
        Args:
            input_file: Path to the FRiP results CSV file
            output_dir: Directory to save visualization outputs
        """
        self.input_file = input_file
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Load data
        self.data = self._load_data()
        
        # Set up color palette
        self.colors = sns.color_palette("husl", 8)
        
    def _load_data(self) -> pd.DataFrame:
        """
        Load and validate FRiP data from CSV file.
        
        Returns:
            DataFrame containing FRiP results
            
        Raises:
            FileNotFoundError: If input file doesn't exist
            ValueError: If data format is invalid
        """
        if not os.path.exists(self.input_file):
            raise FileNotFoundError(f"Input file not found: {self.input_file}")
        
        try:
            data = pd.read_csv(self.input_file)
            logger.info(f"Loaded data with {len(data)} samples")
            
            # Validate required columns
            required_columns = ['Sample', 'Total_Mapped_Reads', 'Reads_in_Peaks', 'FRiP']
            missing_columns = [col for col in required_columns if col not in data.columns]
            
            if missing_columns:
                raise ValueError(f"Missing required columns: {missing_columns}")
            
            # Convert FRiP to numeric, handling any non-numeric values
            data['FRiP'] = pd.to_numeric(data['FRiP'], errors='coerce')
            
            # Remove any rows with NaN values
            data = data.dropna()
            
            if len(data) == 0:
                raise ValueError("No valid data found after cleaning")
            
            logger.info(f"Valid data for {len(data)} samples")
            return data
            
        except Exception as e:
            raise ValueError(f"Error loading data: {str(e)}")
    
    def _extract_sample_groups(self) -> Dict[str, List[str]]:
        """
        Extract sample groups based on sample names.
        
        Returns:
            Dictionary mapping group names to sample lists
        """
        groups = {}
        
        for sample in self.data['Sample']:
            # Extract group from sample name (e.g., "dko-ABA1" -> "dko")
            if '-' in sample:
                group = sample.split('-')[0]
            else:
                group = 'Unknown'
            
            if group not in groups:
                groups[group] = []
            groups[group].append(sample)
        
        return groups
    
    def create_frip_bar_plot(self) -> Figure:
        """
        Create a bar plot showing FRiP values for each sample.
        
        Returns:
            Matplotlib figure object
        """
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Sort data by FRiP values
        sorted_data = self.data.sort_values('FRiP', ascending=False)
        
        # Create bar plot
        bars = ax.bar(range(len(sorted_data)), sorted_data['FRiP'], 
                     color=self.colors[0], alpha=0.7, edgecolor='black')
        
        # Customize plot
        ax.set_xlabel('Samples', fontsize=12, fontweight='bold')
        ax.set_ylabel('FRiP (Fraction of Reads in Peaks)', fontsize=12, fontweight='bold')
        ax.set_title('FRiP Values by Sample', fontsize=14, fontweight='bold', pad=20)
        
        # Set x-axis labels
        ax.set_xticks(range(len(sorted_data)))
        ax.set_xticklabels(sorted_data['Sample'], rotation=45, ha='right')
        
        # Add value labels on bars
        for i, (bar, value) in enumerate(zip(bars, sorted_data['FRiP'])):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.001,
                   f'{value:.4f}', ha='center', va='bottom', fontsize=8)
        
        # Add grid
        ax.grid(True, alpha=0.3, axis='y')
        
        # Adjust layout
        plt.tight_layout()
        
        return fig
    
    def create_frip_box_plot(self) -> Figure:
        """
        Create a box plot showing FRiP distribution by sample groups.
        
        Returns:
            Matplotlib figure object
        """
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Add group information to data
        groups = self._extract_sample_groups()
        group_data = []
        group_labels = []
        
        for group, samples in groups.items():
            group_frip = self.data[self.data['Sample'].isin(samples)]['FRiP']
            if len(group_frip) > 0:
                group_data.append(group_frip.values)
                group_labels.append(group)
        
        # Create box plot
        bp = ax.boxplot(group_data, labels=group_labels, patch_artist=True)
        
        # Color boxes
        for patch, color in zip(bp['boxes'], self.colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        
        # Customize plot
        ax.set_xlabel('Sample Groups', fontsize=12, fontweight='bold')
        ax.set_ylabel('FRiP (Fraction of Reads in Peaks)', fontsize=12, fontweight='bold')
        ax.set_title('FRiP Distribution by Sample Groups', fontsize=14, fontweight='bold', pad=20)
        
        # Add grid
        ax.grid(True, alpha=0.3, axis='y')
        
        # Rotate x-axis labels if needed
        plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
        
        plt.tight_layout()
        
        return fig
    
    def create_reads_correlation_plot(self) -> Figure:
        """
        Create a scatter plot showing correlation between total reads and reads in peaks.
        
        Returns:
            Matplotlib figure object
        """
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Create scatter plot
        scatter = ax.scatter(self.data['Total_Mapped_Reads'], 
                           self.data['Reads_in_Peaks'], 
                           c=self.data['FRiP'], 
                           cmap='viridis', 
                           s=100, 
                           alpha=0.7,
                           edgecolors='black')
        
        # Add colorbar
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('FRiP Value', fontsize=12, fontweight='bold')
        
        # Calculate correlation
        correlation, p_value = stats.pearsonr(self.data['Total_Mapped_Reads'], 
                                            self.data['Reads_in_Peaks'])
        
        # Add correlation line
        z = np.polyfit(self.data['Total_Mapped_Reads'], self.data['Reads_in_Peaks'], 1)
        p = np.poly1d(z)
        ax.plot(self.data['Total_Mapped_Reads'], p(self.data['Total_Mapped_Reads']), 
               "r--", alpha=0.8, linewidth=2)
        
        # Customize plot
        ax.set_xlabel('Total Mapped Reads', fontsize=12, fontweight='bold')
        ax.set_ylabel('Reads in Peaks', fontsize=12, fontweight='bold')
        ax.set_title(f'Correlation: Total Reads vs Reads in Peaks\n'
                    f'Pearson r = {correlation:.3f}, p = {p_value:.3e}', 
                    fontsize=14, fontweight='bold', pad=20)
        
        # Add grid
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        return fig
    
    def create_frip_vs_reads_plot(self) -> Figure:
        """
        Create a scatter plot showing FRiP vs total reads.
        
        Returns:
            Matplotlib figure object
        """
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Create scatter plot
        ax.scatter(self.data['Total_Mapped_Reads'], 
                  self.data['FRiP'], 
                  color=self.colors[1], 
                  s=100, 
                  alpha=0.7,
                  edgecolors='black')
        
        # Calculate correlation
        correlation, p_value = stats.pearsonr(self.data['Total_Mapped_Reads'], 
                                            self.data['FRiP'])
        
        # Add trend line
        z = np.polyfit(self.data['Total_Mapped_Reads'], self.data['FRiP'], 1)
        p = np.poly1d(z)
        ax.plot(self.data['Total_Mapped_Reads'], p(self.data['Total_Mapped_Reads']), 
               "r--", alpha=0.8, linewidth=2)
        
        # Customize plot
        ax.set_xlabel('Total Mapped Reads', fontsize=12, fontweight='bold')
        ax.set_ylabel('FRiP (Fraction of Reads in Peaks)', fontsize=12, fontweight='bold')
        ax.set_title(f'FRiP vs Total Reads\n'
                    f'Pearson r = {correlation:.3f}, p = {p_value:.3e}', 
                    fontsize=14, fontweight='bold', pad=20)
        
        # Add grid
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        return fig
    
    def create_summary_statistics(self) -> Figure:
        """
        Create a summary statistics table as a figure.
        
        Returns:
            Matplotlib figure object
        """
        fig, ax = plt.subplots(figsize=(12, 8))
        ax.axis('tight')
        ax.axis('off')
        
        # Calculate summary statistics
        stats_data = {
            'Metric': ['Number of Samples', 'Mean FRiP', 'Median FRiP', 'Std FRiP', 
                      'Min FRiP', 'Max FRiP', 'Mean Total Reads', 'Mean Reads in Peaks'],
            'Value': [
                len(self.data),
                f"{self.data['FRiP'].mean():.4f}",
                f"{self.data['FRiP'].median():.4f}",
                f"{self.data['FRiP'].std():.4f}",
                f"{self.data['FRiP'].min():.4f}",
                f"{self.data['FRiP'].max():.4f}",
                f"{self.data['Total_Mapped_Reads'].mean():,.0f}",
                f"{self.data['Reads_in_Peaks'].mean():,.0f}"
            ]
        }
        
        # Create table
        col_labels = list(stats_data.keys())
        table = ax.table(cellText=pd.DataFrame(stats_data).values,
                        colLabels=col_labels,
                        cellLoc='center',
                        loc='center',
                        colWidths=[0.6, 0.4])
        
        # Style the table
        table.auto_set_font_size(False)
        table.set_fontsize(12)
        table.scale(1.2, 1.5)
        
        # Color header
        for i in range(len(col_labels)):
            table[(0, i)].set_facecolor('#4CAF50')
            table[(0, i)].set_text_props(weight='bold', color='white')
        
        # Color alternating rows
        for i in range(1, len(stats_data['Metric']) + 1):
            if i % 2 == 0:
                for j in range(len(col_labels)):
                    table[(i, j)].set_facecolor('#f0f0f0')
        
        ax.set_title('FRiP Analysis Summary Statistics', 
                    fontsize=16, fontweight='bold', pad=20)
        
        return fig
    
    def create_all_visualizations(self) -> None:
        """
        Create all visualizations and save them to files.
        """
        logger.info("Creating FRiP visualizations...")
        
        # Create individual plots
        plots = {
            'frip_bar_plot': self.create_frip_bar_plot(),
            'frip_box_plot': self.create_frip_box_plot(),
            'reads_correlation': self.create_reads_correlation_plot(),
            'frip_vs_reads': self.create_frip_vs_reads_plot(),
            'summary_statistics': self.create_summary_statistics()
        }
        
        # Save individual plots
        for name, fig in plots.items():
            output_path = self.output_dir / f"{name}.pdf"
            fig.savefig(output_path, dpi=300, bbox_inches='tight', format='pdf')
            logger.info(f"Saved {name} to {output_path}")
            plt.close(fig)
        
        # Create combined PDF report
        pdf_path = self.output_dir / "FRiP_analysis_report.pdf"
        with PdfPages(pdf_path) as pdf:
            for name, fig in plots.items():
                pdf.savefig(fig, bbox_inches='tight')
                plt.close(fig)
        
        logger.info(f"Saved combined report to {pdf_path}")
        
        # Save summary statistics to CSV
        summary_path = self.output_dir / "FRiP_summary_statistics.csv"
        summary_stats = {
            'Metric': ['Number of Samples', 'Mean FRiP', 'Median FRiP', 'Std FRiP', 
                      'Min FRiP', 'Max FRiP', 'Mean Total Reads', 'Mean Reads in Peaks'],
            'Value': [
                len(self.data),
                self.data['FRiP'].mean(),
                self.data['FRiP'].median(),
                self.data['FRiP'].std(),
                self.data['FRiP'].min(),
                self.data['FRiP'].max(),
                self.data['Total_Mapped_Reads'].mean(),
                self.data['Reads_in_Peaks'].mean()
            ]
        }
        
        pd.DataFrame(summary_stats).to_csv(summary_path, index=False)
        logger.info(f"Saved summary statistics to {summary_path}")
        
        # Save detailed data
        detailed_path = self.output_dir / "FRiP_detailed_results.csv"
        self.data.to_csv(detailed_path, index=False)
        logger.info(f"Saved detailed results to {detailed_path}")


def main():
    """Main function to run the FRiP visualizer."""
    parser = argparse.ArgumentParser(
        description="Visualize FRiP (Fraction of Reads in Peaks) analysis results",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python visualize_FRiP.py FRiP_results.txt output/
  python visualize_FRiP.py FRiP_results.txt output/ --output-dir my_results/
        """
    )
    
    parser.add_argument(
        'input_file',
        help='Path to the FRiP results CSV file'
    )
    
    parser.add_argument(
        'output_dir',
        help='Directory to save visualization outputs'
    )
    
    parser.add_argument(
        '--output-dir',
        dest='output_dir_alt',
        help='Alternative output directory (overrides positional argument)'
    )
    
    args = parser.parse_args()
    
    # Use alternative output directory if provided
    output_dir = args.output_dir_alt if args.output_dir_alt else args.output_dir
    
    try:
        # Create visualizer and generate plots
        visualizer = FRiPVisualizer(args.input_file, output_dir)
        visualizer.create_all_visualizations()
        
        print(f"\n✅ FRiP visualization completed successfully!")
        print(f"📁 Results saved to: {output_dir}")
        print(f"📊 Generated visualizations:")
        print(f"   - FRiP bar plot")
        print(f"   - FRiP box plot by groups")
        print(f"   - Reads correlation plot")
        print(f"   - FRiP vs reads plot")
        print(f"   - Summary statistics")
        print(f"   - Combined PDF report")
        
    except Exception as e:
        logger.error(f"Error during visualization: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main() 