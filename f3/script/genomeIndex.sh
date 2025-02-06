#!/bin/bash
#SBATCH -J STAR_MB
#SBATCH -o STAR_model_build.out
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 15  
#SBATCH --mem-per-cpu 1G
#SBATCH --time 00-12:00

echo "###### Job submitted ######"
date
echo "Requiring 1 task with 10 cpus in total 10G Memory under this DAP-seq user"
echo

echo "###### Environment loading..... ######"
echo "module load Miniconda3/4.9.2 STAR/2.7.9a-GCC-11.2.0"
module load Miniconda3/4.9.2
module load STAR/2.7.9a-GCC-11.2.0
echo
echo "source activate base"
source activate base

path=/RAID1/working/R425/lavakau/rna-seq/marchantia/20240530_la
cd ${path}
echo "###### Current working directory: $path ######"
echo


file_path=${path}/gtf_genome

STAR \
    --runMode genomeGenerate \
    --genomeDir genomeIndex \
    --genomeFastaFiles ${file_path}/MpTak_v6.1r1.genome.fasta \
    --sjdbGTFfile ${file_path}/MpTak_v6.1r1.gtf \
    --runThreadN 10 --genomeSAindexNbases 11