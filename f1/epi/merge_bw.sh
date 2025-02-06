#!/bin/bash
#SBATCH -J MergeBW
#SBATCH -o MergeBW_%a.out
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
echo "module load Miniconda3/4.9.2, SAMtools/1.14-GCC-11.2.0, BEDTools/2.30.0-GCC-10.2.0, FastQC/0.11.9-Java-11, Bowtie2/2.4.4-GCC-11.2.0"
module load Miniconda3/4.9.2
module load SAMtools/1.14-GCC-11.2.0
module load BEDTools/2.30.0-GCC-10.2.0
module load FastQC/0.11.9-Java-11
module load Bowtie2/2.4.4-GCC-11.2.0
echo "source activate base"
source activate base

path=/RAID1/working/R425/lavakau/chromatin
cd ${path};
echo "###### Current working directory: /RAID1/working/R425/lavakau/chromatin ######"
print_separator



###### set up array for parallel running ######
echo "###### set up array for parallel running ######"
# Specify the path to the config file
config=${path}/range_merge.txt
# Extract the sample name for the current $SLURM_ARRAY_TASK_ID
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
print_separator
###########################################


indir=${path}/bowtie2/bam_sort
outdir=${path}/bowtie2/merge_bam

mkdir ${outdir}

###### Merged bam files from biological replicates ######
echo "###### samtools merge and index bam ######"
# merge two replicates
samtools merge ${outdir}/${sample}.merged_mapped_sorted.bam ${indir}/${sample}*.mapped_sorted.bam
# samtools index bam 
samtools index ${outdir}/${sample}.merged_mapped_sorted.bam
############################################


###### Visualization all siganls by all mapped read ######
echo "###### Visualization all signals by all mapped read ######"
echo "Normalized method: RPKM, Reads Per Kilobase per Million mapped reads"
echo "conda ENV: deeptool"
source activate deeptool
echo
outdir2=${path}/bowtie2/bigwig
mkdir ${outdir2}
bamCoverage -b ${outdir}/${sample}.merged_mapped_sorted.bam \
    --normalizeUsing RPKM \
    -p 20 --extendReads 300 \
    -o ${outdir2}/${sample}_merged.bw
print_separator
###########################################


echo "Finished merged bam and export bigwig files!"