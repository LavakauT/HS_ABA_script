#!/bin/bash
#SBATCH -J ATAC-seq-plot_RPKM
#SBATCH -o ATAC-seq_plot_RPKM.out
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 20  
#SBATCH --mem-per-cpu 1G
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

path=/RAID1/working/R425/lavakau/atac-seq/marchantia/20230626
cd ${path}
echo "###### Current working directory: ${path} ######"
print_separator


###### bamCoverage ######
echo "###### bamCoverage ######"
# convert bam into bigwig
# file transformation should take 15-20 minutes
indir=${path}

outdir=${path}/merge_bam_RPKM_bw

mkdir ${outdir}

for file in ${indir}/*-merge.bam
do
    b=$(basename $file -merge.bam)
    bamCoverage -b ${file} \
        --binSize 10 \
        --normalizeUsing RPKM \
        -o ${outdir}/${b}.merge_RPKM.bw 
done
print_separator
###########################################