
#predict DEG by using ATAC-seq data.
import os
import pandas as pd
from pybedtools import BedTool
import pyBigWig
from pyfaidx import Fasta
import numpy as np
from Bio.Seq import Seq
import argparse
from icecream import ic
import ml_util
#from sklearn.feature_extraction.text import CountVectorizer
import joblib

# Desired fixed length for sequences
UPSTREAM=1000
DOWNSTREAM=500
FIXED_LENGTH = UPSTREAM + DOWNSTREAM

def extract_promoter_sequence(row, genome, fixed_length=FIXED_LENGTH):
    chrom = row['chrom']
    start = row['promoter_start']
    end = row['promoter_end']
    strand = row['strand']
    seq = genome[chrom][start:end].seq
    if strand == '-':
        seq = str(Seq(seq).reverse_complement())
    if len(seq) < fixed_length:
        seq += 'N' * (fixed_length - len(seq))
    return seq[:fixed_length].upper()


def main(args):
    '''
    gene_file = '../input_data/Mp_gene.bed'
    genome_file = '../input_data/Mp_genome.fa'
    deg_file = '../input_data/Mp_DEG_HS.txt'
    query = 'Mp_DEG_HS'

    gene_file = "MSUv7_gene_primary.bed"
    #atac_test_file = "ATAC_Os_HS_peaks.narrowPeak"
    #atac_control_file = "ATAC_Os_Control_peaks.narrowPeak"
    genome_file = "MSUv7.fa"
    deg_file = "DEG_Os_HS.txt"
    model_name = 'At_exp'
    '''
    gene_file = args.gene_file
    genome_file = args.genome_file
    tp_file = args.tp
    tn_file = args.tn
    query = args.query_name
    targets_list = args.model_list

    # Load tp and tn file
    tp = pd.read_csv(tp_file, header=None)
    tn = pd.read_csv(tn_file, header=None)

    # Load genome fasta file
    genome = Fasta(genome_file)
    sequence_lengths = {}
    for seq_id in genome.keys():
        sequence_lengths[seq_id] = len(genome[seq_id])

    # Load gene bed and identify promoter regions
    genes = BedTool(gene_file)
    promoter_df = genes.to_dataframe()
    promoter_df['chrom'] = promoter_df['chrom'].astype(str)
    # Calculate promoter start and end positions with boundaries
    promoter_df['promoter_start'] = promoter_df.apply(
        lambda row: max(1, row['start'] - UPSTREAM) if row['strand'] == '+' else max(1, row['end'] - DOWNSTREAM), axis=1)
    promoter_df['promoter_end'] = promoter_df.apply(
        lambda row: min(sequence_lengths[row['chrom']], row['start'] + DOWNSTREAM) if row['strand'] == '+' else min(sequence_lengths[row['chrom']], row['end'] + UPSTREAM), axis=1)
    promoter_df['sequence'] = promoter_df.apply(lambda row: extract_promoter_sequence(row, genome), axis=1)

    # Merge with tp/tn data
    promoter_df['deg'] = np.where(promoter_df['name'].isin(tp[0]), 1, np.where(promoter_df['name'].isin(tn[0]), 0, 2))
    promoter_df = promoter_df[promoter_df['deg'] != 2].reset_index(drop=True)
    ic(promoter_df['deg'].value_counts())

    #seq_df = promoter_df['sequence']
    #y = np.array(promoter_df['deg'].tolist())
    targets = targets_list.split(',')
    
    results = []
    for target in targets:
        target_df = promoter_df
        target_path = f'model/{target}'
        #generate Kmer count matrix from sequence based on target features and scaled
        target_scaler = f'{target_path}/Feature_scaler.pkl'
        scaler = joblib.load(target_scaler)
        selected_features_name = scaler.get_feature_names_out().tolist()
        X_df_selected = ml_util.count_listed_kmers(selected_features_name, target_df['sequence'], count_rc=True)
        X_scaled = scaler.transform(X_df_selected)
        X_scaled_df = pd.DataFrame(X_scaled, columns=selected_features_name)
        y = np.array(target_df['deg'].tolist())

        result = ml_util.ml_testing(X_scaled_df, y, query, target)
        results.append(result)
    out_summary = f'testing_result_{query}.tsv'
    results_df = pd.concat(results)
    results_df.to_csv(out_summary, sep='\t', index=False)
    #ml_util.ml_training_and_testing(X_df, y, model_name)


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Prepare input data for deep learning model.')
    parser.add_argument('-gene_file', required=True, help='Path to the gene BED file.')
    parser.add_argument('-genome_file', required=True, help='Path to the genome FASTA file.')
    parser.add_argument('-tp', required=True, help='TP genes list.')
    parser.add_argument('-tn', required=True, help='TN genes list.')
    parser.add_argument('-model_list', required=True, help='List models for testing. Separated by comma "," ')
    parser.add_argument('-query_name', required=True, help='Output prefix.')
    #parser.add_argument('--output_file', required=True, help='Path to save the output npz file.')

    args = parser.parse_args()
    main(args)

