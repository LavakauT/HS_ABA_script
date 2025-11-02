#!/usr/bin/env python3
"""
Create a comprehensive summary report comparing all peak selection strategies
"""

import os
import pandas as pd
from pathlib import Path
import glob

def create_summary_report():
    """Create a comprehensive summary report for all strategies."""
    
    base_dir = "output/test_all_strategies"
    strategies = ["closest_to_tss", "proximal_median", "distal_median", "all_median"]
    
    print("Creating comprehensive summary report...")
    print("="*60)
    
    # Collect results from all strategies
    results = {}
    
    for strategy in strategies:
        strategy_dir = Path(base_dir) / strategy
        
        if not strategy_dir.exists():
            print(f"Warning: Directory for {strategy} not found")
            continue
            
        # Find the latest analysis summary
        summary_files = list(strategy_dir.glob("analysis_summary_*.txt"))
        if not summary_files:
            print(f"Warning: No summary files found for {strategy}")
            continue
            
        # Get the most recent file
        latest_file = max(summary_files, key=lambda x: x.stat().st_mtime)
        
        print(f"\nStrategy: {strategy}")
        print(f"Latest file: {latest_file.name}")
        
        # Read the summary
        with open(latest_file, 'r') as f:
            content = f.read()
            
        # Extract key metrics
        lines = content.split('\n')
        correlation_matrix_start = None
        pvalue_matrix_start = None
        
        for i, line in enumerate(lines):
            if "Correlation Matrix:" in line:
                correlation_matrix_start = i + 1
            elif "P-value Matrix:" in line:
                pvalue_matrix_start = i + 1
                break
        
        if correlation_matrix_start and pvalue_matrix_start:
            # Extract correlation matrix
            corr_lines = lines[correlation_matrix_start:pvalue_matrix_start-1]
            corr_data = []
            for line in corr_lines:
                if line.strip() and not line.startswith('='):
                    parts = line.split()
                    if len(parts) >= 5:
                        corr_data.append([float(x) for x in parts[1:]])
            
            if corr_data:
                corr_df = pd.DataFrame(corr_data, 
                                     columns=['cluster_1', 'cluster_2', 'cluster_3', 'cluster_4'],
                                     index=['cluster_1', 'cluster_2', 'cluster_3', 'cluster_4'])
                
                # Calculate statistics
                all_correlations = corr_df.values.flatten()
                mean_corr = all_correlations.mean()
                median_corr = pd.Series(all_correlations).median()
                min_corr = all_correlations.min()
                max_corr = all_correlations.max()
                
                results[strategy] = {
                    'correlation_matrix': corr_df,
                    'mean_correlation': mean_corr,
                    'median_correlation': median_corr,
                    'min_correlation': min_corr,
                    'max_correlation': max_corr,
                    'file': latest_file.name
                }
                
                print(f"  Mean correlation: {mean_corr:.3f}")
                print(f"  Median correlation: {median_corr:.3f}")
                print(f"  Min correlation: {min_corr:.3f}")
                print(f"  Max correlation: {max_corr:.3f}")
    
    # Create comparison table
    print("\n" + "="*60)
    print("COMPARISON SUMMARY")
    print("="*60)
    
    comparison_data = []
    for strategy in strategies:
        if strategy in results:
            data = results[strategy]
            comparison_data.append({
                'Strategy': strategy,
                'Mean Correlation': f"{data['mean_correlation']:.3f}",
                'Median Correlation': f"{data['median_correlation']:.3f}",
                'Min Correlation': f"{data['min_correlation']:.3f}",
                'Max Correlation': f"{data['max_correlation']:.3f}",
                'File': data['file']
            })
    
    if comparison_data:
        comparison_df = pd.DataFrame(comparison_data)
        print(comparison_df.to_string(index=False))
        
        # Save comparison table
        output_file = "output/test_all_strategies/strategy_comparison_summary.csv"
        comparison_df.to_csv(output_file, index=False)
        print(f"\nComparison table saved to: {output_file}")
        
        # Find best strategy
        best_strategy = max(results.keys(), key=lambda x: results[x]['mean_correlation'])
        print(f"\nBest strategy by mean correlation: {best_strategy}")
        print(f"Mean correlation: {results[best_strategy]['mean_correlation']:.3f}")
        
    else:
        print("No results found for any strategy")
    
    print("\n" + "="*60)
    print("Analysis completed!")

if __name__ == "__main__":
    create_summary_report() 