#!/bin/bash
#SBATCH -J ChIP-TMQC
#SBATCH -o ChIP-TMQC2_%a.out
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 5  
#SBATCH --mem-per-cpu 2G
#SBATCH --array=1-28
#SBATCH --time 00-12:00


date
echo "###### Job submitted ######"
echo "Requiring 1 task with 10 cpus in total 10G Memory under this ChIP-seq user"
echo

echo "###### Environment loading..... ######"
echo "module load Miniconda3/4.9.2 FastQC/0.11.9-Java-11 Trimmomatic/0.39-Java-11"
module load Miniconda3/4.9.2
module load FastQC/0.11.9-Java-11
module load Trimmomatic/0.39-Java-11
echo
echo "source activate fastp"
echo "To execute Trimmomatic run: java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar"
source activate fastp

path=/RAID1/working/R425/lavakau/chromatin
cd ${path};
echo "###### Current working directory: /RAID1/working/R425/lavakau/chromatin ######"
echo




###### set up array for parallel running ######
# Specify the path to the config file
config=${path}/range_meta.txt
# Extract the sample name for the current $SLURM_ARRAY_TASK_ID
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
input_dir=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)
###########################################


##### Trimming with fastp ######
echo "###### auto detection of adpter under Novaseq 6000 platform and trimming ######"
date

dir=${path}/rawData
outdir=${path}/clean_fq
outdir1=${path}/clean_fq_json
outdir2=${path}/clean_fq_html

mkdir ${outdir}
mkdir ${outdir1}
mkdir ${outdir2}

fastp \
-i ${dir}/${sample}_R1.fastq.gz \
-I ${dir}/${sample}_R2.fastq.gz \
-o ${outdir}/${sample}_R1.clean.fq.gz \
-O ${outdir}/${sample}_R2.clean.fq.gz \
-j ${outdir1}/${sample}.fastp.json \
-h ${outdir2}/${sample}.fastp.html
##########################################

##### QC Trimming file ######
echo "###### QC Trimming file ######"
date
outdir3=${path}/clean_fq_qc
mkdir ${outdir3}
fastqc ${outdir}/${sample}_R1.clean.fq.gz -o ${outdir3}
fastqc ${outdir}/${sample}_R2.clean.fq.gz -o ${outdir3}
##########################################


##### Multi-QC Trimming file as one reports ######
echo "###### Multi-QC Trimming file as one reports ######"
echo
outdir3=${path}/clean_fq_qc

if [[ $SLURM_ARRAY_TASK_ID -eq 28 ]]
then
    multiqc .${ourdir3}
else
    echo "ArrayTaskID is not 28. Skipping execution."
fi
##########################################

echo "Finished!"