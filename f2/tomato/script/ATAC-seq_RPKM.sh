#!/bin/bash
#SBATCH -J ATAC_seq_plot_RPKM
#SBATCH -o ATAC_seq_plot_RPKM_%a.out
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 15  
#SBATCH --mem-per-cpu 1G
#SBATCH --array=1-16%8
#SBATCH --time 00-12:00

# Function to print a separator line
print_separator() {
    printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' -
}

date
echo "###### Job submitted ######"
echo "Requiring total $SLURM_NTASKS task(s) with $SLURM_CPUS_PER_TASK cpus in total $SLURM_MEM_PER_CPU Memory under this RNA-seq user"
print_separator
echo "###### Environment loading..... ######"
echo "module load Miniconda3 FastQC SAMtools BEDTools Bowtie2"
module load Miniconda3/4.9.2
module load FastQC/0.11.9-Java-11
module load SAMtools/1.14-GCC-11.2.0
module load BEDTools/2.30.0-GCC-10.2.0
module load Bowtie2/2.4.4-GCC-11.2.0
echo "source activate deeptool"
source activate deeptool

path=/RAID1/working/R425/lavakau/atac-seq/marchantia/time_course
cd ${path}
echo "###### Current working directory: ${path} ######"
print_separator



###### set up array for parallel running ######
echo "###### Specify the path to the config file ######"
config=${path}/range3.txt

# Extract the sample name for the current $SLURM_ARRAY_TASK_ID
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
print_separator
###########################################



###### bamCoverage ######
echo "###### bamCoverage ######"
# convert bam into bigwig
# file transformation should take 15-20 minutes
indir=${path}/merge_flt_bam

outdir=${path}/merge_flt_bam_RPKM_bw
outdir2=${path}/merge_flt_bam_BPM_bw

mkdir ${outdir}

file=${indir}/${sample}.flt.bam
# bamCoverage -b ${file} \
#     --binSize 10 \
#     --normalizeUsing RPKM \
#     -o ${outdir}/${sample}.merge_RPKM.bw

bamCoverage -b ${file} \
    --binSize 10 \
    --normalizeUsing BPM \
    -o ${outdir2}/${sample}.merge.bw
print_separator
###########################################