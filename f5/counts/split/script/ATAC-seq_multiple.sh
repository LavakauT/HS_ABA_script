#!/bin/bash
#SBATCH -J ATAC-mapping
#SBATCH -o ATAC-seq_mapping.out
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 15  
#SBATCH --mem-per-cpu 1G
#SBATCH --array=10
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
echo "source activate fastp"
source activate fastp

path=/RAID1/working/R425/lavakau/atac-seq/marchantia/20230928
cd ${path}
echo "###### Current working directory: ${path} ######"
print_separator


###### set up array for parallel running ######
echo "###### Specify the path to the config file ######"
config=${path}/range.txt

# Extract the sample name for the current $SLURM_ARRAY_TASK_ID
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
print_separator
###########################################



###### bowtie2 mapping ######
echo "###### bowtie2 mapping ######"
echo "geonme index: mar (If you do not have. please build it from bowtie2 manual)"
echo "Warning: It takes some time............."

indir=${path}/fastp/clean_fq

outdir=${path}/sorted_bam
outdir2=${path}/merge_sorted_markdup_bam
outdir3=${path}/merge_flt_bam
outdir4=${path}/merge_log
# outdir5=${path}/rm_organelle
outdir6=${path}/merge_sorted_bam

mkdir ${outdir}
mkdir ${outdir2}
mkdir ${outdir3}
mkdir ${outdir4}
# mkdir outdir5
mkdir ${outdir6}

# replicate 1
bowtie2 -p 10 -x mar/mar \
-1 ${indir}/${sample}1_R1.clean.fq.gz \
-2 ${indir}/${sample}1_R2.clean.fq.gz 2> ${outdir4}/${sample}1.log\
| samtools sort -@ 20 -O bam -o ${outdir}/${sample}1.sorted.bam -
samtools index ${outdir}/${sample}1.sorted.bam

# replicate 2
bowtie2 -p 10 -x mar/mar \
-1 ${indir}/${sample}2_R1.clean.fq.gz \
-2 ${indir}/${sample}2_R2.clean.fq.gz 2> ${outdir4}/${sample}2.log\
| samtools sort -@ 20 -O bam -o ${outdir}/${sample}2.sorted.bam -
samtools index ${outdir}/${sample}2.sorted.bam
print_separator
###########################################



###### samtools index and filter bam ######
echo "###### samtools index and filter bam ######"
# merge two replicates
samtools merge ${outdir6}/${sample}.sorted.bam ${outdir4}/${sample}1.sorted.bam ${outdir4}/${sample}2.sorted.bam
# samtools index bam 
samtools index ${outdir6}/${sample}.sorted.bam

# sambamba mark the duplicated reads
sambamba markdup ${outdir6}/${sample}.sorted.bam ${outdir2}/${sample}.sorted.markdup.bam 2> ${outdir4}/${sample}.log

# remove low-quality mapped and duplicated reads
samtools view -@ 20 -bF 1804 -q 20 ${outdir2}/${sample}.sorted.markdup.bam -o ${outdir3}/${sample}.flt.bam
samtools index -@ 20 ${outdir3}/${sample}.flt.bam

sambamba markdup ${outdir3}/${sample}.sorted.bam ${outdir2}/${sample}.sorted.markdup.bam 2> ${outdir3}/${sample}.log;
samtools view -@ 20 -bF 1804 -q 20 ${outdir2}/${sample}.sorted.markdup.bam -o ${outdir3}/${sample}.flt.bam;
samtools index -@ 20 ${outdir3}/${sample}.flt.bam

# bedtools remove ChrC and ChrM reads
# The bed.file (-b) here is the region we want to remove
# If, there is no region you want to remove, then just skip following code:
# bedtools intersect -abam ${outdir3}/${sample}.flt.bam -b region you want to remove -v > ${outdir5}/${sample}.rm_organelle.bam
# samtools index ${outdir5}/${sample}.rm_organelle.bam
print_separator
###########################################



###### macs2 peak calling ######
echo "###### macs2 peak calling ######"
source activate ml
macs2 callpeak -t ${outdir3}/${sample}.flt.bam \
    -f BAMPE \
    -n ${sample} \
    -g 1.1e8 \
    --cutoff-analysis --nomodel --shift -100 --extsize 200 \
    --outdir ${sample}.peakcall
print_separator 
###########################################


###### count macs2 peak calling ######
# wc -l *.narrowPeak
###########################################