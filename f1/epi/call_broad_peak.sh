#!/bin/bash
#SBATCH -J CallPeak
#SBATCH -o CallPeak_macs2_%a.out
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 20  
#SBATCH --mem-per-cpu 1G
#SBATCH --array=1-13
#SBATCH --time 01-00:00


# Function to print a separator line
print_separator() {
    printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' -
}

date
echo "###### Job submitted ######"
echo "Requiring 1 task with 15 cpus in total 15G Memory under this DAP-seq user"
print_separator

echo "###### Environment loading..... ######"
echo "module load Miniconda3/4.9.2, SAMtools/1.14-GCC-11.2.0, BEDTools/2.30.0-GCC-10.2.0, FastQC/0.11.9-Java-11, Bowtie2/2.4.4-GCC-11.2.0 Homer/4.11-GCC-11.2.0"
module load Miniconda3/4.9.2
module load SAMtools/1.14-GCC-11.2.0
module load BEDTools/2.30.0-GCC-10.2.0
module load FastQC/0.11.9-Java-11
module load Bowtie2/2.4.4-GCC-11.2.0
module load Homer/4.11-GCC-11.2.0
echo "source activate base"
source activate base

path=/RAID1/working/R425/lavakau/chromatin
cd ${path};
echo "###### Current working directory: /RAID1/working/R425/lavakau/chromatin ######"
print_separator



###### set up array for parallel running ######
echo "###### set up array for parallel running ######"
# Specify the path to the config file
config=${path}/range_callpeak.txt
# Extract the sample name for the current $SLURM_ARRAY_TASK_ID
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
# extsize=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)
print_separator
###########################################


###### macs2 peak calling ######
echo "###### macs2 peak calling ######"

indir=${path}/bowtie2/merge_bam
outdir=${path}/bowtie2/macs2


mkdir ${outdir}

# -t: ChIP-seq treatment file. If multiple files are given
# as '-t A B C', then they will all be read and pooled
# together. REQUIRED.

# -c: Control file. If multiple files are given as '-c A B
# C', they will be pooled to estimate ChIP-seq
# background noise.

# -g: Effective genome size. for mp is 1.1e8 (larger than Arabidopsis 1.6e7)
# -B: Store all details
# -m: Select the regions within MFOLD range of low-high
# confidence enrichment ratio against background to
# build model.
# -q: Minimum FDR (q-value) cutoff for peak detection.

# background file: NA
# specific condition: broad range peak calling

mkdir ${outdir}/${sample}

macs2 callpeak -t ${indir}/${sample}.merged_mapped_sorted.bam \
    -f BAM --outdir ${outdir}/${sample} -g 1.6e7 -n ${sample} -B -q 0.05 \
    --nomodel --broad --broad-cutoff 0.05 \
# --extsize ${extsize}, macs2 can automatically detect the extsize. 
print_separator
###########################################




# ###### HOMER peak calling ######
# echo "###### HOMER peak calling ######"

# outdir2=${path}/bowtie2/homer

# mkdir ${outdir2}

# # Create a "tag directory" for ChIP-seq experiment using the makeTagDirectory command.
# makeTagDirectory ${outdir2}/${sample} ${indir}/${sample}.merged_mapped_sorted.bam


# # histone (histone modification ChIP-Seq, region based, uses -region -size 500 -L 0, regions.txt)
# findPeaks ${outdir2}/${sample}/ \
#     -style histone -o ${outdir2}/${sample}/${sample}.peaks -region -size 1000 -minDist 2500
# print_separator
# ###########################################

echo "FInished calling peaks by Macs2 and HOMER"