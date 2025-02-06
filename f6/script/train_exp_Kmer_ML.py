

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
from sklearn.feature_extraction.text import CountVectorizer
import joblib
from tqdm import tqdm

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


class KmerTransformer:
    from sklearn.feature_extraction.text import CountVectorizer
    def __init__(self, k):
        self.k = k
        self.vectorizer = CountVectorizer(analyzer=self._generate_kmers)

    def _generate_kmers(self, sequence):
        kmers = [sequence[i:i+self.k] for i in range(len(sequence) - self.k + 1)]
        return [kmer for kmer in kmers if all(char in 'ATCG' for char in kmer)]
    
    def fit(self, sequences):
        self.vectorizer.fit(sequences)
        self.kmer_features = self.vectorizer.get_feature_names_out()
        return self
    
    def transform(self, sequences):
        return self.vectorizer.transform(sequences)
    
    def get_kmer_features(self):
        return self.kmer_features


def kmer_task(k, sequences, filter_size):
    # Helper function for reverse complement
    def reverse_complement(kmer):
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return ''.join(complement[base] for base in reversed(kmer))

    # Helper function to get the lexicographically smaller of kmer and its reverse complement
    def get_canonical_kmer(kmer):
        rc_kmer = reverse_complement(kmer)
        return min(kmer, rc_kmer)
    
    # Helper function to check if kmer consists of only one type of character
    def is_not_homopolymer(kmer):
        return len(set(kmer)) != 1
    
    # Generate Kmers from input sequence
    kmer_transformer = KmerTransformer(k)
    kmer_transformer.fit(sequences)
    kmer_features = kmer_transformer.get_kmer_features()
    kmer_matrix = kmer_transformer.transform(sequences)

    canonical_kmer_to_index = {}
    #index_to_canonical_kmer = {}
    for i, kmer in enumerate(kmer_features):
        canonical_kmer = get_canonical_kmer(kmer)
        #index_to_canonical_kmer[i] = canonical_kmer
        if canonical_kmer not in canonical_kmer_to_index:
            canonical_kmer_to_index[canonical_kmer] = []
        canonical_kmer_to_index[canonical_kmer].append(i)

    # Calculate canonical kmer count and filter out low count kmer and homopolymer 
    kmer_counts = kmer_matrix.getnnz(axis=0) #showing how many genes has count
    canonical_kmer_counts = {}
    filtered_canonical_kmers = []
    
    for canonical_kmer, indices in canonical_kmer_to_index.items():
        # Sum up the counts for all indices corresponding to the canonical k-mer
        total_count = sum(kmer_counts[i] for i in indices)
        canonical_kmer_counts[canonical_kmer] = total_count
        if total_count > filter_size and canonical_kmer and is_not_homopolymer(canonical_kmer):
            filtered_canonical_kmers.append(canonical_kmer)
    
    def sum_fr_kmers(batch_start):
        batch_canonical_matrix = []
        batch_end = min(batch_start + batch_size, len(filtered_canonical_kmers))
        batch_canonical_kmers = filtered_canonical_kmers[batch_start:batch_end]
        batch_indices = sorted(set(
            idx for canonical_kmer in batch_canonical_kmers
            for idx in canonical_kmer_to_index[canonical_kmer]
        ))

        # Convert only the necessary columns to dense format
        batch_matrix = kmer_matrix[:, batch_indices].toarray()
        
        # Sum K-mer values between K-mer and its reverse complement
        index_mapping = {index: pos for pos, index in enumerate(batch_indices)}
        for canonical_kmer in batch_canonical_kmers:
            indices = canonical_kmer_to_index[canonical_kmer]

        #for canonical_kmer, indices in filtered_canonical_kmers[batch_start:batch_end]:
            if len(indices) > 1:
                combined_column = np.sum(batch_matrix[:, [index_mapping[idx] for idx in indices]], axis=1)
            else:
                combined_column = batch_matrix[:, index_mapping[indices[0]]]
            batch_canonical_matrix.append(combined_column)
        
        batch_kmer_df = pd.DataFrame(np.array(batch_canonical_matrix).T, columns=batch_canonical_kmers)
        return batch_kmer_df

    batch_size = 1000
    kmer_df_list = []
    # Process K-mers in batches
    for batch_start in tqdm(range(0, len(filtered_canonical_kmers), batch_size), desc=f"Combining F/R kmers {k}"):
        kmer_df_list.append(sum_fr_kmers(batch_start))

    kmer_df = pd.concat(kmer_df_list, axis = 1)
    num_kmer = kmer_df.shape[1]
    print(f'k={k} has {num_kmer} kmers')

    return kmer_df


def gen_kmer_table(sequences, k_min, k_max, filter_size):
    from concurrent.futures import ProcessPoolExecutor, as_completed
    merged_kmer = pd.DataFrame()
    with ProcessPoolExecutor() as executor:
        futures = [executor.submit(kmer_task, k, sequences, filter_size) for k in range(k_min, k_max)]
        for future in as_completed(futures):
            kmer_df = future.result()
            if merged_kmer.empty:
                merged_kmer = kmer_df
            else:
                merged_kmer = pd.concat([merged_kmer, kmer_df], axis=1)

    return merged_kmer


def kmer_distribution(sequences, kmer, bs = 10):
    # Calculate the distribution of kmer and its reverse complement in bins along the sequences.
    def reverse_complement(kmer):
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return ''.join(complement[base] for base in reversed(kmer))

    def find_positions(sequence, kmer):
        """Find all positions of kmer and its reverse complement in the sequence."""
        rc_kmer = reverse_complement(kmer)
        positions = set()
        
        # Find positions of kmer
        pos = sequence.find(kmer)
        while pos != -1:
            positions.add(pos)
            pos = sequence.find(kmer, pos + 1)
        
        # Find positions of reverse complement
        pos = sequence.find(rc_kmer)
        while pos != -1:
            positions.add(pos)
            pos = sequence.find(rc_kmer, pos + 1)
        
        return positions

    kmer_counts = []
    
    for sequence in sequences:
        positions = find_positions(sequence, kmer)
        bin_counts = [0] * ((len(sequence) + bs - 1) // bs)  # Initialize bins
        
        for pos in positions:
            bin_index = pos // bs
            bin_counts[bin_index] += 1
        
        kmer_counts.append(bin_counts)
    
    return np.sum(kmer_counts, axis = 0) / len(sequences) / bs * len(kmer)


def main(args):
    ic(args)
    gene_file = args.gene_file
    genome_file = args.genome_file
    tp_file = args.tp
    tn_file = args.tn
    filter = args.filter
    model_name = args.model_name

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
    #promoter_df = promoter_df.merge(gene_pro_df, how='inner', left_on='name', right_on='name')
    #promoter_df = merged_df.drop(columns=['gene_id'])

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

    #generate Kmer count matrix from sequence
    k_min = 5
    k_max = 13
    filter_size = len(promoter_df) * filter
    X_df = gen_kmer_table(promoter_df['sequence'], k_min, k_max, filter_size)
    y = np.array(promoter_df['deg'].tolist())

    ml_util.ml_training_and_testing(X_df, y, model_name, n_feature=100)

    # Generate Kmer distribution
    kmers_dis = []
    kmers_names = []
    target_path = f'model/{model_name}'
    target_scaler = f'{target_path}/Feature_scaler.pkl'
    scaler = joblib.load(target_scaler)
    selected_features_name = scaler.get_feature_names_out().tolist()
    for kmer in selected_features_name:
        kmers_dis.append(kmer_distribution(promoter_df[promoter_df['deg'] == 1].sequence, kmer, bs = 20))
        kmers_names.append(f'{kmer}_1')
        kmers_dis.append(kmer_distribution(promoter_df[promoter_df['deg'] == 0].sequence, kmer, bs = 20))
        kmers_names.append(f'{kmer}_0')
    kmer_dist_df = pd.DataFrame(np.array(kmers_dis), index = kmers_names)
    out_kmer_dist = f'{target_path}/kmer_distribution.tsv'
    kmer_dist_df.to_csv(out_kmer_dist, header=False, sep='\t')
    out_positive = f'{target_path}/Positive_set.tsv'
    promoter_df[promoter_df['deg'] == 1].iloc[:,:7].to_csv(out_positive, header=False, sep='\t')

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Prepare input data for deep learning model.')
    parser.add_argument('-gene_file', required=True, help='Path to the gene BED file.')
    parser.add_argument('-genome_file', required=True, help='Path to the genome FASTA file.')
    parser.add_argument('-tp', required=True, help='TP genes list.')
    parser.add_argument('-tn', required=True, help='TN genes list.')
    parser.add_argument('-model_name', required=True, help='Output model prefix.')
    parser.add_argument('-filter', required=False, help='Filtering out the Kmer in less than potion of genes. (Defalue 0.05)', default=0.05, type=float)
    args = parser.parse_args()
    main(args)

