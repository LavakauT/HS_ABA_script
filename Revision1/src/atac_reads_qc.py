#!/usr/bin/env python3
"""
ATAC-seq QC Report Processing Script

This script is used to process HTML and JSON files from ATAC-seq fastq completion,
extract important QC metrics and generate complete report tables and visualization charts.

Date: 2025 (HS_ABA project)
"""

import os
import json
import re
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from bs4 import BeautifulSoup
import logging
from datetime import datetime

# 設定中文字體
plt.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'SimHei', '標楷體']
plt.rcParams['axes.unicode_minus'] = False

# 設定日誌
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


@dataclass
class QCData:
    """QC data structure"""
    sample_name: str
    source_file: str
    mean_length_before_filtering: Optional[float] = None
    mean_length_after_filtering: Optional[float] = None
    duplication_rate: Optional[float] = None
    insert_size_peak: Optional[float] = None
    before_filtering_total_reads: Optional[int] = None
    before_filtering_q30_bases: Optional[int] = None
    before_filtering_gc_content: Optional[float] = None
    after_filtering_total_reads: Optional[int] = None
    after_filtering_q30_bases: Optional[int] = None
    after_filtering_gc_content: Optional[float] = None
    reads_passed_filters: Optional[int] = None
    reads_with_low_quality: Optional[int] = None
    reads_with_too_many_n: Optional[int] = None
    reads_too_short: Optional[int] = None


class ATACQCProcessor:
    """ATAC-seq QC report processor"""
    
    def __init__(self, input_dir: str):
        """
        Initialize processor
        
        Args:
            input_dir: Directory path containing QC report files
        """
        self.input_dir = Path(input_dir)
        self.qc_data_list: List[QCData] = []
        
    def find_qc_files(self) -> List[Path]:
        """
        Find all QC report files recursively
        
        Returns:
            List containing all QC file paths
        """
        qc_files = []
        
        # Find HTML files recursively
        html_files = list(self.input_dir.rglob("*.html"))
        qc_files.extend(html_files)
        
        # Find JSON files recursively
        json_files = list(self.input_dir.rglob("*.json"))
        qc_files.extend(json_files)
        
        logger.info(f"Found {len(qc_files)} QC files")
        return qc_files
    
    def extract_sample_name(self, file_path: Path) -> str:
        """
        Extract sample name from file name
        
        Args:
            file_path: File path
            
        Returns:
            Sample name
        """
        # Remove file extension
        name = file_path.stem
        
        # Remove common file extensions that might be part of the filename
        # Remove .fastp, .json, .html extensions from the stem
        extensions_to_remove = ['.fastp', '.json', '.html', '.qc', '.report', '.stats']
        for ext in extensions_to_remove:
            if name.endswith(ext):
                name = name[:-len(ext)]
        
        # Try to extract sample name from file name
        # Common patterns: sample_name_xxx.html or sample_name_xxx.json
        patterns = [
            r'(.+?)_(fastqc|qc|report|stats|fastp)',
            r'(.+?)_(before|after)',
            r'(.+?)_(filtered|raw)',
            r'(.+?)_(\d+)',
        ]
        
        for pattern in patterns:
            match = re.search(pattern, name, re.IGNORECASE)
            if match:
                return match.group(1)
        
        # If no pattern matches, return the cleaned name
        return name
    
    def parse_html_file(self, file_path: Path) -> Optional[QCData]:
        """
        Parse HTML file and extract QC data
        
        Args:
            file_path: HTML file path
            
        Returns:
            QCData object or None
        """
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                content = f.read()
            
            soup = BeautifulSoup(content, 'html.parser')
            sample_name = self.extract_sample_name(file_path)
            qc_data = QCData(sample_name=sample_name, source_file=str(file_path))
            
            # Extract various QC metrics
            self._extract_basic_metrics(soup, qc_data)
            self._extract_length_metrics(soup, qc_data)
            self._extract_quality_metrics(soup, qc_data)
            self._extract_content_metrics(soup, qc_data)
            self._extract_filtering_metrics(soup, qc_data)
            
            return qc_data
            
        except Exception as e:
            logger.error(f"Error parsing HTML file {file_path}: {e}")
            return None
    
    def parse_json_file(self, file_path: Path) -> Optional[QCData]:
        """
        Parse JSON file and extract QC data
        
        Args:
            file_path: JSON file path
            
        Returns:
            QCData object or None
        """
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                data = json.load(f)
            
            sample_name = self.extract_sample_name(file_path)
            qc_data = QCData(sample_name=sample_name, source_file=str(file_path))
            
            # Extract data from JSON
            self._extract_json_metrics(data, qc_data)
            
            return qc_data
            
        except Exception as e:
            logger.error(f"Error parsing JSON file {file_path}: {e}")
            return None
    
    def _extract_basic_metrics(self, soup: BeautifulSoup, qc_data: QCData) -> None:
        """Extract basic metrics"""
        # Find duplication rate
        duplication_patterns = [
            r'duplication.*?(\d+\.?\d*)%',
            r'duplicate.*?(\d+\.?\d*)%',
            r'(\d+\.?\d*)%.*?duplicate'
        ]
        
        text = soup.get_text()
        for pattern in duplication_patterns:
            match = re.search(pattern, text, re.IGNORECASE)
            if match:
                qc_data.duplication_rate = float(match.group(1))
                break
    
    def _extract_length_metrics(self, soup: BeautifulSoup, qc_data: QCData) -> None:
        """Extract length-related metrics"""
        # Find average length
        length_patterns = [
            r'mean.*?length.*?(\d+\.?\d*)',
            r'average.*?length.*?(\d+\.?\d*)',
            r'(\d+\.?\d*).*?bp.*?average'
        ]
        
        text = soup.get_text()
        for pattern in length_patterns:
            match = re.search(pattern, text, re.IGNORECASE)
            if match:
                length = float(match.group(1))
                if qc_data.mean_length_before_filtering is None:
                    qc_data.mean_length_before_filtering = length
                else:
                    qc_data.mean_length_after_filtering = length
                break
    
    def _extract_quality_metrics(self, soup: BeautifulSoup, qc_data: QCData) -> None:
        """Extract quality-related metrics"""
        # Find Q30 bases
        q30_patterns = [
            r'Q30.*?(\d+)',
            r'(\d+).*?Q30',
            r'quality.*?30.*?(\d+)'
        ]
        
        text = soup.get_text()
        for pattern in q30_patterns:
            match = re.search(pattern, text, re.IGNORECASE)
            if match:
                q30_bases = int(match.group(1))
                if qc_data.before_filtering_q30_bases is None:
                    qc_data.before_filtering_q30_bases = q30_bases
                else:
                    qc_data.after_filtering_q30_bases = q30_bases
                break
    
    def _extract_content_metrics(self, soup: BeautifulSoup, qc_data: QCData) -> None:
        """Extract content-related metrics"""
        # Find GC content
        gc_patterns = [
            r'GC.*?content.*?(\d+\.?\d*)%',
            r'(\d+\.?\d*)%.*?GC',
            r'GC.*?(\d+\.?\d*)'
        ]
        
        text = soup.get_text()
        for pattern in gc_patterns:
            match = re.search(pattern, text, re.IGNORECASE)
            if match:
                gc_content = float(match.group(1))
                if qc_data.before_filtering_gc_content is None:
                    qc_data.before_filtering_gc_content = gc_content
                else:
                    qc_data.after_filtering_gc_content = gc_content
                break
    
    def _extract_filtering_metrics(self, soup: BeautifulSoup, qc_data: QCData) -> None:
        """Extract filtering-related metrics"""
        # Find total reads
        reads_patterns = [
            r'total.*?reads.*?(\d+)',
            r'(\d+).*?total.*?reads',
            r'reads.*?(\d+)'
        ]
        
        text = soup.get_text()
        for pattern in reads_patterns:
            match = re.search(pattern, text, re.IGNORECASE)
            if match:
                reads = int(match.group(1))
                if qc_data.before_filtering_total_reads is None:
                    qc_data.before_filtering_total_reads = reads
                else:
                    qc_data.after_filtering_total_reads = reads
                break
    
    def _extract_json_metrics(self, data: Dict[str, Any], qc_data: QCData) -> None:
        """Extract metrics from JSON data (fastp format)"""
        try:
            # Handle fastp JSON format
            if 'summary' in data and isinstance(data['summary'], dict):
                summary = data['summary']
                
                # Extract before filtering metrics
                if 'before_filtering' in summary and isinstance(summary['before_filtering'], dict):
                    before_data = summary['before_filtering']
                    if 'total_reads' in before_data:
                        qc_data.before_filtering_total_reads = int(before_data['total_reads'])
                    if 'q30_bases' in before_data:
                        qc_data.before_filtering_q30_bases = int(before_data['q30_bases'])
                    if 'gc_content' in before_data:
                        qc_data.before_filtering_gc_content = float(before_data['gc_content'])
                    if 'read1_mean_length' in before_data:
                        qc_data.mean_length_before_filtering = float(before_data['read1_mean_length'])
                
                # Extract after filtering metrics
                if 'after_filtering' in summary and isinstance(summary['after_filtering'], dict):
                    after_data = summary['after_filtering']
                    if 'total_reads' in after_data:
                        qc_data.after_filtering_total_reads = int(after_data['total_reads'])
                    if 'q30_bases' in after_data:
                        qc_data.after_filtering_q30_bases = int(after_data['q30_bases'])
                    if 'gc_content' in after_data:
                        qc_data.after_filtering_gc_content = float(after_data['gc_content'])
                    if 'read1_mean_length' in after_data:
                        qc_data.mean_length_after_filtering = float(after_data['read1_mean_length'])
            
            # Extract duplication rate
            if 'duplication' in data and isinstance(data['duplication'], dict):
                if 'rate' in data['duplication']:
                    qc_data.duplication_rate = float(data['duplication']['rate']) * 100  # Convert to percentage
            
            # Extract insert size peak
            if 'insert_size' in data and isinstance(data['insert_size'], dict):
                if 'peak' in data['insert_size']:
                    qc_data.insert_size_peak = float(data['insert_size']['peak'])
            
            # Extract filtering statistics
            if 'filtering_result' in data and isinstance(data['filtering_result'], dict):
                filtering_data = data['filtering_result']
                if 'passed_filter_reads' in filtering_data:
                    qc_data.reads_passed_filters = int(filtering_data['passed_filter_reads'])
                if 'low_quality_reads' in filtering_data:
                    qc_data.reads_with_low_quality = int(filtering_data['low_quality_reads'])
                if 'too_many_N_reads' in filtering_data:
                    qc_data.reads_with_too_many_n = int(filtering_data['too_many_N_reads'])
                if 'too_short_reads' in filtering_data:
                    qc_data.reads_too_short = int(filtering_data['too_short_reads'])
            
            # Try alternative key names for backward compatibility
            metrics_mapping = {
                'duplication_rate': ['duplicate_rate', 'duplicates'],
                'mean_length': ['mean_length', 'avg_length', 'length'],
                'total_reads': ['total_reads', 'reads', 'total'],
                'q30_bases': ['q30_bases', 'q30', 'quality_30'],
                'gc_content': ['gc_content', 'gc', 'gc_percent'],
                'insert_size_peak': ['insert_peak', 'peak_size']
            }
            
            for attr, keys in metrics_mapping.items():
                for key in keys:
                    if key in data:
                        value = data[key]
                        if attr == 'duplication_rate' and qc_data.duplication_rate is None:
                            qc_data.duplication_rate = float(value)
                        elif attr == 'mean_length' and qc_data.mean_length_before_filtering is None:
                            qc_data.mean_length_before_filtering = float(value)
                        elif attr == 'total_reads' and qc_data.before_filtering_total_reads is None:
                            qc_data.before_filtering_total_reads = int(value)
                        elif attr == 'q30_bases' and qc_data.before_filtering_q30_bases is None:
                            qc_data.before_filtering_q30_bases = int(value)
                        elif attr == 'gc_content' and qc_data.before_filtering_gc_content is None:
                            qc_data.before_filtering_gc_content = float(value)
                        elif attr == 'insert_size_peak' and qc_data.insert_size_peak is None:
                            qc_data.insert_size_peak = float(value)
                        break
                        
        except Exception as e:
            logger.error(f"Error extracting JSON metrics: {e}")
            return None
    
    def process_all_files(self) -> None:
        """Process all QC files"""
        qc_files = self.find_qc_files()
        
        for file_path in qc_files:
            logger.info(f"Processing file: {file_path}")
            
            if file_path.suffix.lower() == '.html':
                qc_data = self.parse_html_file(file_path)
            elif file_path.suffix.lower() == '.json':
                qc_data = self.parse_json_file(file_path)
            else:
                logger.warning(f"Unsupported file format: {file_path}")
                continue
            
            if qc_data:
                self.qc_data_list.append(qc_data)
        
        logger.info(f"Successfully processed {len(self.qc_data_list)} files")
    
    def create_summary_table(self) -> pd.DataFrame:
        """
        Create QC data summary table
        
        Returns:
            DataFrame containing all QC data
        """
        if not self.qc_data_list:
            logger.warning("No QC data available to create table")
            return pd.DataFrame()
        
        # Convert QCData objects to dictionary list
        data_dicts = []
        for qc_data in self.qc_data_list:
            data_dict = {
                'Sample Name': qc_data.sample_name,
                'Source File': qc_data.source_file,
                'Mean Length Before Filtering': qc_data.mean_length_before_filtering,
                'Mean Length After Filtering': qc_data.mean_length_after_filtering,
                'Duplication Rate (%)': qc_data.duplication_rate,
                'Insert Size Peak': qc_data.insert_size_peak,
                'Total Reads Before Filtering': qc_data.before_filtering_total_reads,
                'Q30 Bases Before Filtering': qc_data.before_filtering_q30_bases,
                'GC Content Before Filtering (%)': qc_data.before_filtering_gc_content,
                'Total Reads After Filtering': qc_data.after_filtering_total_reads,
                'Q30 Bases After Filtering': qc_data.after_filtering_q30_bases,
                'GC Content After Filtering (%)': qc_data.after_filtering_gc_content,
                'Reads Passed Filters': qc_data.reads_passed_filters,
                'Reads with Low Quality': qc_data.reads_with_low_quality,
                'Reads with Too Many N': qc_data.reads_with_too_many_n,
                'Reads Too Short': qc_data.reads_too_short
            }
            data_dicts.append(data_dict)
        
        df = pd.DataFrame(data_dicts)
        return df
    
    def create_comparison_plot(self, output_dir: str = "output") -> None:
        """
        Create before and after filtering comparison plots
        
        Args:
            output_dir: Output directory
        """
        if not self.qc_data_list:
            logger.warning("No data available to create plots")
            return
        
        # Create output directory
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)
        
        # Prepare data
        samples = []
        before_reads = []
        after_reads = []
        before_q30 = []
        after_q30 = []
        before_gc = []
        after_gc = []
        
        for qc_data in self.qc_data_list:
            if (qc_data.before_filtering_total_reads is not None and 
                qc_data.after_filtering_total_reads is not None):
                samples.append(qc_data.sample_name)
                before_reads.append(qc_data.before_filtering_total_reads)
                after_reads.append(qc_data.after_filtering_total_reads)
                before_q30.append(qc_data.before_filtering_q30_bases or 0)
                after_q30.append(qc_data.after_filtering_q30_bases or 0)
                before_gc.append(qc_data.before_filtering_gc_content or 0)
                after_gc.append(qc_data.after_filtering_gc_content or 0)
        
        if not samples:
            logger.warning("Insufficient data to create comparison plots")
            return
        
        # Create plots
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        fig.suptitle('ATAC-seq QC Before vs After Filtering Comparison', fontsize=16, fontweight='bold')
        
        # 1. Total reads comparison
        ax1 = axes[0, 0]
        x = range(len(samples))
        width = 0.35
        
        ax1.bar([i - width/2 for i in x], before_reads, width, 
                label='Before Filtering', color='lightblue', alpha=0.8)
        ax1.bar([i + width/2 for i in x], after_reads, width, 
                label='After Filtering', color='lightcoral', alpha=0.8)
        
        ax1.set_xlabel('Sample')
        ax1.set_ylabel('Total Reads')
        ax1.set_title('Total Reads Comparison')
        ax1.set_xticks(x)
        ax1.set_xticklabels(samples, rotation=45)
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # 2. Q30 bases comparison
        ax2 = axes[0, 1]
        ax2.bar([i - width/2 for i in x], before_q30, width, 
                label='Before Filtering', color='lightblue', alpha=0.8)
        ax2.bar([i + width/2 for i in x], after_q30, width, 
                label='After Filtering', color='lightcoral', alpha=0.8)
        
        ax2.set_xlabel('Sample')
        ax2.set_ylabel('Q30 Bases')
        ax2.set_title('Q30 Bases Comparison')
        ax2.set_xticks(x)
        ax2.set_xticklabels(samples, rotation=45)
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        # 3. GC content comparison
        ax3 = axes[1, 0]
        ax3.bar([i - width/2 for i in x], before_gc, width, 
                label='Before Filtering', color='lightblue', alpha=0.8)
        ax3.bar([i + width/2 for i in x], after_gc, width, 
                label='After Filtering', color='lightcoral', alpha=0.8)
        
        ax3.set_xlabel('Sample')
        ax3.set_ylabel('GC Content (%)')
        ax3.set_title('GC Content Comparison')
        ax3.set_xticks(x)
        ax3.set_xticklabels(samples, rotation=45)
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        
        # 4. Retention rate
        ax4 = axes[1, 1]
        retention_rates = [(after/before)*100 for before, after in zip(before_reads, after_reads)]
        ax4.bar(samples, retention_rates, color='lightgreen', alpha=0.8)
        ax4.set_xlabel('Sample')
        ax4.set_ylabel('Retention Rate (%)')
        ax4.set_title('Post-Filtering Retention Rate')
        ax4.set_xticklabels(samples, rotation=45)
        ax4.grid(True, alpha=0.3)
        
        # Add value labels on bars
        for i, rate in enumerate(retention_rates):
            ax4.text(i, rate + 1, f'{rate:.1f}%', ha='center', va='bottom')
        
        plt.tight_layout()
        plt.savefig(output_path / 'atac_qc_comparison.png', dpi=300, bbox_inches='tight')
        plt.savefig(output_path / 'atac_qc_comparison.pdf', bbox_inches='tight')
        plt.show()
        
        logger.info(f"Comparison plots saved to {output_path}")
    
    def calculate_retention_statistics(self) -> pd.DataFrame:
        """
        Calculate filtering retention statistics
        
        Returns:
            DataFrame containing retention statistics
        """
        stats_data = []
        
        for qc_data in self.qc_data_list:
            if (qc_data.before_filtering_total_reads is not None and 
                qc_data.after_filtering_total_reads is not None):
                
                before_reads = qc_data.before_filtering_total_reads
                after_reads = qc_data.after_filtering_total_reads
                retention_rate = (after_reads / before_reads) * 100
                
                stats_dict = {
                    'Sample Name': qc_data.sample_name,
                    'Reads Before Filtering': before_reads,
                    'Reads After Filtering': after_reads,
                    'Retained Reads': after_reads,
                    'Filtered Reads': before_reads - after_reads,
                    'Retention Rate (%)': retention_rate,
                    'Filtering Rate (%)': 100 - retention_rate
                }
                stats_data.append(stats_dict)
        
        return pd.DataFrame(stats_data)
    
    def save_results(self, output_dir: str = "output") -> None:
        """
        Save all results
        
        Args:
            output_dir: Output directory
        """
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)
        
        # Save summary table
        summary_df = self.create_summary_table()
        if not summary_df.empty:
            summary_df.to_csv(output_path / 'atac_qc_summary.csv', index=False, encoding='utf-8-sig')
            summary_df.to_excel(output_path / 'atac_qc_summary.xlsx', index=False)
            logger.info(f"Summary table saved to {output_path}")
        
        # Save retention statistics
        retention_df = self.calculate_retention_statistics()
        if not retention_df.empty:
            retention_df.to_csv(output_path / 'atac_retention_stats.csv', index=False, encoding='utf-8-sig')
            retention_df.to_excel(output_path / 'atac_retention_stats.xlsx', index=False)
            logger.info(f"Retention statistics saved to {output_path}")
        
        # Create comparison plots
        self.create_comparison_plot(output_dir)
        
        # Create detailed report
        self._create_detailed_report(output_path)
    
    def _create_detailed_report(self, output_path: Path) -> None:
        """
        Create detailed report
        
        Args:
            output_path: Output path
        """
        report_content = []
        report_content.append("# ATAC-seq QC Report")
        report_content.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        report_content.append(f"Files processed: {len(self.qc_data_list)}")
        report_content.append("")
        
        # Add detailed information for each sample
        for i, qc_data in enumerate(self.qc_data_list, 1):
            report_content.append(f"## Sample {i}: {qc_data.sample_name}")
            report_content.append(f"Source file: {qc_data.source_file}")
            report_content.append("")
            
            if qc_data.mean_length_before_filtering:
                report_content.append(f"- Mean length before filtering: {qc_data.mean_length_before_filtering:.2f} bp")
            if qc_data.mean_length_after_filtering:
                report_content.append(f"- Mean length after filtering: {qc_data.mean_length_after_filtering:.2f} bp")
            if qc_data.duplication_rate:
                report_content.append(f"- Duplication rate: {qc_data.duplication_rate:.2f}%")
            if qc_data.insert_size_peak:
                report_content.append(f"- Insert size peak: {qc_data.insert_size_peak:.2f} bp")
            
            if qc_data.before_filtering_total_reads:
                report_content.append(f"- Total reads before filtering: {qc_data.before_filtering_total_reads:,}")
            if qc_data.after_filtering_total_reads:
                report_content.append(f"- Total reads after filtering: {qc_data.after_filtering_total_reads:,}")
                
                if qc_data.before_filtering_total_reads:
                    retention = (qc_data.after_filtering_total_reads / qc_data.before_filtering_total_reads) * 100
                    report_content.append(f"- Retention rate: {retention:.2f}%")
            
            report_content.append("")
        
        # Save report
        with open(output_path / 'atac_qc_report.md', 'w', encoding='utf-8') as f:
            f.write('\n'.join(report_content))
        
        logger.info(f"Detailed report saved to {output_path / 'atac_qc_report.md'}")


def main():
    """Main function"""
    import argparse
    
    parser = argparse.ArgumentParser(description='ATAC-seq QC Report Processing Tool')
    parser.add_argument('input_dir', nargs='?', default='data', help='Directory path containing QC report files (default: data)')
    parser.add_argument('-o', '--output', default='output', help='Output directory (default: output)')
    parser.add_argument('-v', '--verbose', action='store_true', help='Show detailed logs')
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Check input directory
    if not os.path.exists(args.input_dir):
        logger.error(f"Input directory does not exist: {args.input_dir}")
        logger.info("Please provide a valid directory path containing QC report files")
        logger.info("Example: python src/atac_reads_qc.py data/HS_fastp")
        return
    
    # Create processor and process files
    processor = ATACQCProcessor(args.input_dir)
    processor.process_all_files()
    
    if processor.qc_data_list:
        # Save results
        processor.save_results(args.output)
        
        # Display summary
        summary_df = processor.create_summary_table()
        print("\n=== QC Data Summary ===")
        print(summary_df.to_string(index=False))
        
        retention_df = processor.calculate_retention_statistics()
        if not retention_df.empty:
            print("\n=== Retention Statistics ===")
            print(retention_df.to_string(index=False))
    else:
        logger.warning("No valid QC data found")
        logger.info("Please check that your input directory contains valid QC report files")
        logger.info("Supported formats: HTML and JSON files from fastp or other QC tools")


if __name__ == "__main__":
    main() 