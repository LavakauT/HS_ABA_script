#!/bin/bash
#SBATCH -J FeatureCounts
#SBATCH -o FeatureCounts_mp.out
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 10  
#SBATCH --mem-per-cpu 1G
#SBATCH --time 01-00:00


# Function to print a separator line
print_separator() {
    printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' -
}

echo "###### Job submitted ######"
date
echo "Requiring $SLURM_NTASKS task with $SLURM_CPUS_PER_TASK cpus in total $SLURM_MEM_PER_CPU Memory under this RNA-seq user"
print_separator

echo "###### Environment loading..... ######"
echo "module load Miniconda3/4.9.2"
module load Miniconda3/4.9.2
module load STAR/2.7.9a-GCC-11.2.0
print_separator
echo "source activate featurecounts"
source activate featurecounts

path=/RAID1/working/R425/lavakau/rna-seq/marchantia/20240530_la
cd ${path}

echo "###### Current working directory: $path ######"
print_separator

echo "###### Export paired-end reads to all their overlapping meta-features ######"
name=counts
indir=/RAID1/working/R425/lavakau/rna-seq/marchantia/20240530_la/fastp/bam
genome=/RAID1/working/R425/lavakau/rna-seq/marchantia/20240530_la/gtf_genome

# -t is the feature
# -g is the corresponding information
# please adjust them based on your gtf or gff files

featureCounts -p -O -T 10 \
    -t gene \
    -g ID \
    -a ${genome}/MpTak_v6.1r1.gff \
    -o ${name}.txt ${indir}/*.bam
echo "###### Finished FeatureCounts ######"
###########################################