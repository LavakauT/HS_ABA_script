#!/bin/bash
#SBATCH -J Consensus
#SBATCH -o Consensus_mp.out
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 12  
#SBATCH --mem-per-cpu 1G
#SBATCH --time 00-12:00

date
echo "job submitted"
echo "requiring tasks under this K-mer user"

echo "Environment loading....."

cd /RAID1/working/R425/lavakau/kmer/mp
echo "python 3.10 envs activate: R4.3.2"
source /RAID1/working/R425/lavakau/miniconda3/etc/profile.d/conda.sh
conda activate R4.3.2

Rscript consensus_motif.R
Rscript consensus_motif_tfbm.R