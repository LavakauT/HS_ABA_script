#!/usr/bin/env python3
"""
Samtools Flagstat Results Visualizer

This script reads samtools flagstat results from a tab-separated file and creates
various visualizations including bar charts, heatmaps, and summary statistics.

Usage:
    python visualize_flagstat.py <flagstat_results_file>

Example:
    python visualize_flagstat.py flagstat_results.txt
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
import sys
from typing import Optional, Tuple, List
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Set style for better looking plots
plt.style.use('default')
sns.set_palette("husl")


def load_flagstat_data(file_path: str) -> pd.DataFrame:
    """
    Load flagstat results from a tab-separated file.
    
    Args:
        file_path: Path to the flagstat results file
        
    Returns:
        DataFrame containing the flagstat results
        
    Raises:
        FileNotFoundError: If the file doesn't exist
        ValueError: If the file format is invalid
    """
    try:
        # Read the file, skipping comment lines
        df = pd.read_csv(file_path, sep='\t', comment='#')
        
        if df.empty:
            raise ValueError("No data found in the file")
        
        # Clean column names by removing ".sorted" from all column names
        # Handle both ".sorted" and ".sorted." patterns
        df.columns = df.columns.str.replace('.sorted', '', regex=False)
        # Also remove any trailing dots that might be left
        df.columns = df.columns.str.rstrip('.')
        
        # Clean sample names in the first column by removing ".sorted"
        if 'Sample' in df.columns:
            df['Sample'] = df['Sample'].str.replace('.sorted', '', regex=False)
            
        logger.info(f"Loaded {len(df)} samples from {file_path}")
        return df
        
    except FileNotFoundError:
        logger.error(f"File not found: {file_path}")
        raise
    except Exception as e:
        logger.error(f"Error reading file: {e}")
        raise ValueError(f"Invalid file format: {e}")


def create_mapping_stats_plot(df: pd.DataFrame, output_dir: str) -> None:
    """
    Create a bar plot showing mapping statistics for each sample.
    
    Args:
        df: DataFrame containing flagstat results
        output_dir: Directory to save the plot
    """
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('Samtools Flagstat - Mapping Statistics', fontsize=16, fontweight='bold')
    
    # Calculate percentages
    df['Mapping_Rate'] = (df['Mapped'] / df['Total_Reads'] * 100).round(2)
    df['Properly_Paired_Rate'] = (df['Properly_Paired'] / df['Paired'] * 100).round(2)
    df['Duplicate_Rate'] = (df['Duplicates'] / df['Total_Reads'] * 100).round(2)
    df['Singleton_Rate'] = (df['Singletons'] / df['Paired'] * 100).round(2)
    
    # Plot 1: Mapping Rate
    axes[0, 0].bar(df['Sample'], df['Mapping_Rate'], color='skyblue', alpha=0.7)
    axes[0, 0].set_title('Mapping Rate (%)')
    axes[0, 0].set_ylabel('Mapping Rate (%)')
    axes[0, 0].tick_params(axis='x', rotation=45)
    axes[0, 0].grid(True, alpha=0.3)
    
    # Plot 2: Properly Paired Rate
    axes[0, 1].bar(df['Sample'], df['Properly_Paired_Rate'], color='lightgreen', alpha=0.7)
    axes[0, 1].set_title('Properly Paired Rate (%)')
    axes[0, 1].set_ylabel('Properly Paired Rate (%)')
    axes[0, 1].tick_params(axis='x', rotation=45)
    axes[0, 1].grid(True, alpha=0.3)
    
    # Plot 3: Duplicate Rate
    axes[1, 0].bar(df['Sample'], df['Duplicate_Rate'], color='salmon', alpha=0.7)
    axes[1, 0].set_title('Duplicate Rate (%)')
    axes[1, 0].set_ylabel('Duplicate Rate (%)')
    axes[1, 0].tick_params(axis='x', rotation=45)
    axes[1, 0].grid(True, alpha=0.3)
    
    # Plot 4: Singleton Rate
    axes[1, 1].bar(df['Sample'], df['Singleton_Rate'], color='gold', alpha=0.7)
    axes[1, 1].set_title('Singleton Rate (%)')
    axes[1, 1].set_ylabel('Singleton Rate (%)')
    axes[1, 1].tick_params(axis='x', rotation=45)
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    output_path = Path(output_dir) / 'mapping_statistics.pdf'
    plt.savefig(output_path, dpi=300, bbox_inches='tight', format='pdf')
    plt.close()
    logger.info(f"Mapping statistics plot saved to: {output_path}")


def create_read_distribution_plot(df: pd.DataFrame, output_dir: str) -> None:
    """
    Create a stacked bar plot showing read distribution.
    
    Args:
        df: DataFrame containing flagstat results
        output_dir: Directory to save the plot
    """
    # Prepare data for stacked bar plot
    plot_data = df[['Sample', 'Mapped', 'Secondary', 'Supplementary']].copy()
    plot_data['Unmapped'] = df['Total_Reads'] - df['Mapped']
    
    # Create stacked bar plot
    fig, ax = plt.subplots(figsize=(12, 8))
    
    bottom = np.zeros(len(df))
    colors = ['#2E8B57', '#FF6B6B', '#4ECDC4', '#45B7D1']
    labels = ['Mapped', 'Unmapped', 'Secondary', 'Supplementary']
    
    for i, (col, color) in enumerate(zip(['Mapped', 'Unmapped', 'Secondary', 'Supplementary'], colors)):
        ax.bar(df['Sample'], plot_data[col], bottom=bottom, label=labels[i], color=color, alpha=0.8)
        bottom += plot_data[col]
    
    ax.set_title('Read Distribution by Sample', fontsize=14, fontweight='bold')
    ax.set_xlabel('Sample')
    ax.set_ylabel('Number of Reads')
    ax.legend()
    ax.tick_params(axis='x', rotation=45)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    output_path = Path(output_dir) / 'read_distribution.pdf'
    plt.savefig(output_path, dpi=300, bbox_inches='tight', format='pdf')
    plt.close()
    logger.info(f"Read distribution plot saved to: {output_path}")


def create_heatmap(df: pd.DataFrame, output_dir: str) -> None:
    """
    Create a heatmap showing correlation between different metrics.
    
    Args:
        df: DataFrame containing flagstat results
        output_dir: Directory to save the plot
    """
    # Calculate correlation matrix for numerical columns
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    correlation_matrix = df[numeric_cols].corr()
    
    # Create heatmap
    plt.figure(figsize=(12, 10))
    mask = np.triu(np.ones_like(correlation_matrix, dtype=bool))
    
    sns.heatmap(correlation_matrix, 
                mask=mask,
                annot=True, 
                cmap='coolwarm', 
                center=0,
                square=True,
                fmt='.2f',
                cbar_kws={"shrink": .8})
    
    plt.title('Correlation Matrix of Flagstat Metrics', fontsize=14, fontweight='bold')
    plt.tight_layout()
    output_path = Path(output_dir) / 'correlation_heatmap.pdf'
    plt.savefig(output_path, dpi=300, bbox_inches='tight', format='pdf')
    plt.close()
    logger.info(f"Correlation heatmap saved to: {output_path}")


def create_summary_statistics(df: pd.DataFrame, output_dir: str) -> None:
    """
    Create a summary statistics table and save as CSV.
    
    Args:
        df: DataFrame containing flagstat results
        output_dir: Directory to save the summary
    """
    # Calculate summary statistics
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    summary_stats = df[numeric_cols].describe()
    
    # Add percentage calculations
    percentage_stats = pd.DataFrame({
        'Mapping_Rate_Mean': (df['Mapped'] / df['Total_Reads'] * 100).mean(),
        'Mapping_Rate_Std': (df['Mapped'] / df['Total_Reads'] * 100).std(),
        'Properly_Paired_Rate_Mean': (df['Properly_Paired'] / df['Paired'] * 100).mean(),
        'Properly_Paired_Rate_Std': (df['Properly_Paired'] / df['Paired'] * 100).std(),
        'Duplicate_Rate_Mean': (df['Duplicates'] / df['Total_Reads'] * 100).mean(),
        'Duplicate_Rate_Std': (df['Duplicates'] / df['Total_Reads'] * 100).std(),
    }, index=[0])
    
    # Save summary statistics
    summary_path = Path(output_dir) / 'summary_statistics.csv'
    summary_stats.to_csv(summary_path)
    percentage_stats.to_csv(Path(output_dir) / 'percentage_statistics.csv')
    
    logger.info(f"Summary statistics saved to: {summary_path}")
    
    # Print summary to console
    print("\n" + "="*60)
    print("SUMMARY STATISTICS")
    print("="*60)
    print(f"Number of samples: {len(df)}")
    print(f"Average mapping rate: {percentage_stats['Mapping_Rate_Mean'].iloc[0]:.2f}% ± {percentage_stats['Mapping_Rate_Std'].iloc[0]:.2f}%")
    print(f"Average properly paired rate: {percentage_stats['Properly_Paired_Rate_Mean'].iloc[0]:.2f}% ± {percentage_stats['Properly_Paired_Rate_Std'].iloc[0]:.2f}%")
    print(f"Average duplicate rate: {percentage_stats['Duplicate_Rate_Mean'].iloc[0]:.2f}% ± {percentage_stats['Duplicate_Rate_Std'].iloc[0]:.2f}%")
    print("="*60)


def create_quality_metrics_plot(df: pd.DataFrame, output_dir: str) -> None:
    """
    Create a comprehensive quality metrics plot.
    
    Args:
        df: DataFrame containing flagstat results
        output_dir: Directory to save the plot
    """
    # Calculate quality metrics
    df['QC_Pass_Rate'] = (df['QC_Passed'] / df['Total_Reads'] * 100).round(2)
    df['Mapping_Rate'] = (df['Mapped'] / df['Total_Reads'] * 100).round(2)
    df['Properly_Paired_Rate'] = (df['Properly_Paired'] / df['Paired'] * 100).round(2)
    df['Duplicate_Rate'] = (df['Duplicates'] / df['Total_Reads'] * 100).round(2)
    
    # Create subplots
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Quality Metrics Overview', fontsize=16, fontweight='bold')
    
    # Plot 1: QC Pass Rate
    axes[0, 0].bar(df['Sample'], df['QC_Pass_Rate'], color='lightblue', alpha=0.8)
    axes[0, 0].set_title('QC Pass Rate')
    axes[0, 0].set_ylabel('QC Pass Rate (%)')
    axes[0, 0].tick_params(axis='x', rotation=45)
    axes[0, 0].axhline(y=95, color='red', linestyle='--', alpha=0.7, label='95% threshold')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    # Plot 2: Mapping Rate
    axes[0, 1].bar(df['Sample'], df['Mapping_Rate'], color='lightgreen', alpha=0.8)
    axes[0, 1].set_title('Mapping Rate')
    axes[0, 1].set_ylabel('Mapping Rate (%)')
    axes[0, 1].tick_params(axis='x', rotation=45)
    axes[0, 1].axhline(y=70, color='red', linestyle='--', alpha=0.7, label='70% threshold')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)
    
    # Plot 3: Properly Paired Rate
    axes[1, 0].bar(df['Sample'], df['Properly_Paired_Rate'], color='lightcoral', alpha=0.8)
    axes[1, 0].set_title('Properly Paired Rate')
    axes[1, 0].set_ylabel('Properly Paired Rate (%)')
    axes[1, 0].tick_params(axis='x', rotation=45)
    axes[1, 0].axhline(y=70, color='red', linestyle='--', alpha=0.7, label='70% threshold')
    axes[1, 0].legend()
    axes[1, 0].grid(True, alpha=0.3)
    
    # Plot 4: Duplicate Rate
    axes[1, 1].bar(df['Sample'], df['Duplicate_Rate'], color='gold', alpha=0.8)
    axes[1, 1].set_title('Duplicate Rate')
    axes[1, 1].set_ylabel('Duplicate Rate (%)')
    axes[1, 1].tick_params(axis='x', rotation=45)
    axes[1, 1].axhline(y=20, color='red', linestyle='--', alpha=0.7, label='20% threshold')
    axes[1, 1].legend()
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    output_path = Path(output_dir) / 'quality_metrics.pdf'
    plt.savefig(output_path, dpi=300, bbox_inches='tight', format='pdf')
    plt.close()
    logger.info(f"Quality metrics plot saved to: {output_path}")


def main():
    """Main function to run the flagstat visualizer."""
    parser = argparse.ArgumentParser(
        description='Visualize samtools flagstat results',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python visualize_flagstat.py flagstat_results.txt
    python visualize_flagstat.py results.txt --output-dir plots/
        """
    )
    
    parser.add_argument('input_file', help='Path to the flagstat results file')
    parser.add_argument('--output-dir', default='flagstat_plots', 
                       help='Output directory for plots (default: flagstat_plots)')
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)
    
    try:
        # Load data
        logger.info(f"Loading data from: {args.input_file}")
        df = load_flagstat_data(args.input_file)
        
        # Create visualizations
        logger.info("Creating visualizations...")
        
        create_mapping_stats_plot(df, output_dir)
        create_read_distribution_plot(df, output_dir)
        create_heatmap(df, output_dir)
        create_quality_metrics_plot(df, output_dir)
        create_summary_statistics(df, output_dir)
        
        logger.info(f"All visualizations completed! Results saved to: {output_dir}")
        
    except Exception as e:
        logger.error(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main() 