#!/bin/bash
#SBATCH -J RNA_bw_export
#SBATCH -o RNA_bw_export_%a.out
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 20  
#SBATCH --mem-per-cpu 1G
#SBATCH --array=1-8
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

path=/RAID1/working/R425/lavakau/rna-seq/marchantia/20230801_la
cd ${path}
echo "###### Current working directory: ${path} ######"
print_separator


###### set up array for parallel running ######
echo "###### Specify the path to the config file ######"
config=${path}/range_rna_bam.txt

# Extract the sample name for the current $SLURM_ARRAY_TASK_ID
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
sample1=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)
sample2=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $config)
sample3=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $5}' $config)
print_separator
###########################################



###### directory setting ######
echo "###### directory setting ######"

indir=${path}/bam
outdir=${path}/merge_bam

mkdir ${outdir}

print_separator
###########################################




###### samtools merge and index bam ######
echo "###### samtools index and filter bam ######"
end=Aligned.sortedByCoord.out
# merge two replicates
samtools merge ${outdir}/${sample}.bam ${indir}/${sample1}${end}.bam ${indir}/${sample2}${end}.bam ${indir}/${sample3}${end}.bam
# samtools index bam 
samtools index ${outdir}/${sample}.bam

print_separator
###########################################



###### bamCoverage ######
echo "###### bamCoverage ######"
# convert bam into bigwig
# file transformation should take 15-20 minutes
indir=${path}/merge_bam

outdir=${path}/merge_bam_bw

mkdir ${outdir}

for file in ${indir}/*.bam
do
    b=$(basename $file .bam)
    bamCoverage -b ${file} \
        --normalizeUsing RPKM \
        -o ${outdir}/${b}.bw 
done
print_separator
###########################################