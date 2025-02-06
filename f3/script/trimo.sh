#!/bin/bash
#SBATCH -J TMQC
#SBATCH -o TMQC_mp_%A-%a.out
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 5  
#SBATCH --mem-per-cpu 2G
#SBATCH --array=1-43
#SBATCH --time 00-12:00

date
echo "###### Job submitted ######"
echo "Requiring 1 task with 10 cpus in total 10G Memory under this DAP-seq user"

echo "###### Environment loading..... ######"
echo "module load Miniconda3/4.9.2 FastQC/0.11.9-Java-11 Trimmomatic/0.39-Java-11"
module load Miniconda3/4.9.2
module load FastQC/0.11.9-Java-11
module load Trimmomatic/0.39-Java-11
echo
echo "source activate fastp"
echo "To execute Trimmomatic run: java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar"
source activate fastp

cd /RAID1/working/R425/lavakau/rna-seq/marchantia/20240530_la/rowData
echo "###### Current working directory: /RAID1/working/R425/lavakau/rna-seq/marchantia/20240530_la/rowData ######"
echo

cd /home;
cd /RAID1/working/R425/lavakau/rna-seq/marchantia/20240530_la/rowData;
/RAID1/working/R425/lavakau/rna-seq/marchantia/20240530_la/rowData
###### set up array for parallel running ######
# Specify the path to the config file
config=/RAID1/working/R425/lavakau/rna-seq/marchantia/20240530_la/rowData/range.txt
# Extract the sample name for the current $SLURM_ARRAY_TASK_ID
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
###########################################


###### Trimming with fastp ######
# echo "###### auto detection of adpter under Novaseq 6000 platform and trimming ######"
# date

# dir=/RAID1/working/R425/lavakau/rna-seq/marchantia/20240530_la/rowData/hsfb_test/rawData
# outdir=/RAID1/working/R425/lavakau/rna-seq/marchantia/20240530_la/rowData/hsfb_test/clean_fq
# outdir1=/RAID1/working/R425/lavakau/rna-seq/marchantia/20240530_la/rowData/hsfb_test/clean_fq_json
# outdir2=/RAID1/working/R425/lavakau/rna-seq/marchantia/20240530_la/rowData/hsfb_test/clean_fq_html

# fastp \
# -i ${dir}/${sample}_2237NNLT4_R1.fastq.gz \
# -I ${dir}/${sample}_2237NNLT4_R2.fastq.gz \
# -o ${outdir}/${sample}_R1.clean.fq.gz \
# -O ${outdir}/${sample}_R2.clean.fq.gz \
# -j ${outdir1}/${sample}.fastp.json \
# -h ${outdir2}/${sample}.fastp.html
###########################################

###### QC Trimming file ######
# echo "###### QC Trimming file ######"
# date
# outdir3=/RAID1/working/R425/lavakau/rna-seq/marchantia/20240530_la/rowData/hsfb_test/clean_fq_qc
# mkdir ${outdir3}
# fastqc ${outdir}/${sample}_R1.clean.fq.gz -o ${outdir3}
# fastqc ${outdir}/${sample}_R2.clean.fq.gz -o ${outdir3}
###########################################

###### Multi-QC Trimming file as one reports ######
# echo "###### Multi-QC Trimming file as one reports ######"
# echo
# outdir3=/RAID1/working/R425/lavakau/rna-seq/marchantia/20240530_la/rowData/hsfb_test/clean_fq_qc
# multiqc .${ourdir3}
###########################################


## REMEMBER TO CHANGE THE LOG NAMES!!!
###### Trimmomatic method ######
echo 
echo "###### auto detection of adpter under Novaseq 6000 platform and trimming by Trimmomatic ######"
dir=/RAID1/working/R425/lavakau/rna-seq/marchantia/20240530_la/rowData
outdir=/RAID1/working/R425/lavakau/rna-seq/marchantia/20240530_la/trimo_clean_fq
seq=223C37LT4
mkdir ${outdir}

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar  \
    PE -phred33 ${dir}/${sample}_${seq}_R1.fastq.gz ${dir}/${sample}_${seq}_R2.fastq.gz ${outdir}/${sample}_R1.clean.fq.gz ${outdir}/${sample}_R1.clean_unpaired.fq.gz ${outdir}/${sample}_R2.clean.fq.gz ${outdir}/${sample}_R2.clean_unpaired.fq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 
###########################################

###### QC Trimming file ######
echo "###### QC Trimming file ######"
outdir3=/RAID1/working/R425/lavakau/rna-seq/marchantia/20240530_la/trimo_clean_fq_qc
mkdir ${outdir3}
fastqc ${outdir}/${sample}_R1.clean.fq.gz -o ${outdir3}
fastqc ${outdir}/${sample}_R2.clean.fq.gz -o ${outdir3}
###########################################

###### Multi-QC Trimming file as one reports ######
echo "###### Multi-QC Trimming file as one reports ######"
outdir3=/RAID1/working/R425/lavakau/rna-seq/marchantia/20240530_la/trimo_clean_fq_qc
multiqc .${ourdir3}
###########################################