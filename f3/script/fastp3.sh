#!/bin/bash
#SBATCH -J FastpQC
#SBATCH -o Fastp_mp_3.out
#SBATCH --ntasks 10
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task 5  
#SBATCH --mem-per-cpu 2G
#SBATCH --time 00-12:00

date
echo "###### Job submitted ######"
echo "Requiring total $SLURM_NTASKS task(s) with $SLURM_CPUS_PER_TASK cpus in total $SLURM_MEM_PER_CPU Memory under this RNA-seq user"
echo
echo "###### Environment loading..... ######"
echo "module load Miniconda3/4.9.2 FastQC/0.11.9-Java-11"
module load Miniconda3/4.9.2
module load FastQC/0.11.9-Java-11
echo "source activate fastp"
source activate fastp

cd /RAID1/working/R425/lavakau/rna-seq/marchantia/20240530_la
echo "###### Current working directory: /RAID1/working/R425/lavakau/rna-seq/marchantia/20240530_la/rawData ######"
echo


###### set up output folder and sequencing numbers ######
dir=/RAID1/working/R425/lavakau/rna-seq/marchantia/20240530_la/rawData
out=/RAID1/working/R425/lavakau/rna-seq/marchantia/20240530_la/fastp
outdir=/RAID1/working/R425/lavakau/rna-seq/marchantia/20240530_la/fastp/clean_fq
outdir1=/RAID1/working/R425/lavakau/rna-seq/marchantia/20240530_la/fastp/clean_fq_json
outdir2=/RAID1/working/R425/lavakau/rna-seq/marchantia/20240530_la/fastp/clean_fq_html
outdir3=/RAID1/working/R425/lavakau/rna-seq/marchantia/20240530_la/fastp/clean_fq_qc
seq=223C37LT4

mkdir ${out}
mkdir ${outdir}
mkdir ${outdir1}
mkdir ${outdir2}
mkdir ${outdir3}
###########################################


###### set up array for parallel running ######
num_tasks=$SLURM_NTASKS
# Specify the path to the config file
config=/RAID1/working/R425/lavakau/rna-seq/marchantia/20240530_la/rawData/range.txt
###########################################

for ((i=21; i>num_tasks; i++))
do
    # Extract the sample name for the current $SLURM_PROCID
    index=${i}
    sample=$(awk -v ArrayTaskID=$index '$1==ArrayTaskID {print $2}' $config)
    echo "Task $index is processing $sample"

    ###### Trimming with fastp ######
    echo "###### auto detection of adpter under Novaseq 6000 platform and trimming ######"
    date

    fastp \
    -i ${dir}/${sample}_${seq}_R1.fastq.gz \
    -I ${dir}/${sample}_${seq}_R2.fastq.gz \
    -o ${outdir}/${sample}_R1.clean.fq.gz \
    -O ${outdir}/${sample}_R2.clean.fq.gz \
    -j ${outdir1}/${sample}.fastp.json \
    -h ${outdir2}/${sample}.fastp.html
    ###########################################

    ###### QC Trimming file ######
    echo "###### QC Trimming file ######"
    date

    fastqc ${outdir}/${sample}_R1.clean.fq.gz -o ${outdir3}
    fastqc ${outdir}/${sample}_R2.clean.fq.gz -o ${outdir3}
    ###########################################


    echo "Task $SLURM_PROCID has finished processing $sample"
done

# wait for finishing all the samples
wait

###### Multi-QC Trimming file as one reports ######
echo "###### Multi-QC Trimming file as one reports ######"
echo
outdir3=/RAID1/working/R425/lavakau/rna-seq/marchantia/20240530_la/fastp/clean_fq_qc
multiqc .${ourdir3}
###########################################