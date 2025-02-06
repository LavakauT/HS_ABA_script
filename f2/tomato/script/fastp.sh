#!/bin/bash
#SBATCH -J FastpQC
#SBATCH -o Fastp_tomato_%a.out
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 20  
#SBATCH --mem-per-cpu 1G
#SBATCH --array=1-6
#SBATCH --time 00-12:00


# Function to print a separator line
print_separator() {
    printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' -
}

date
echo "###### Job submitted ######"
echo "Requiring total $SLURM_NTASKS task(s) with $SLURM_CPUS_PER_TASK cpus in total $SLURM_MEM_PER_CPU Memory under this ATAC-seq user"
echo
echo "###### Environment loading..... ######"
echo "module load Miniconda3/4.9.2 FastQC/0.11.9-Java-11"
module load Miniconda3/4.9.2
module load FastQC/0.11.9-Java-11
echo "source activate fastp"
source activate fastp

path=/RAID1/working/R425/huanchi/atac/tomato
cd ${path}
echo "###### Current working directory: ${path}/rowData ######"


###### set up output folder and sequencing numbers ######
dir=${path}/rawData
out=${path}/fastp
outdir=${path}/fastp/clean_fq
outdir1=${path}/fastp/clean_fq_json
outdir2=${path}/fastp/clean_fq_html
outdir3=${path}/fastp/clean_fq_qc
# seq=H3NCKDSX7

mkdir ${out}
mkdir ${outdir}
mkdir ${outdir1}
mkdir ${outdir2}
mkdir ${outdir3}
###########################################


###### set up array for parallel running ######
# Specify the path to the config file
config=${path}/tomato_hsfa1a_atac.txt
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)
###########################################




###### Trimming with fastp ######
echo "###### auto detection of adpter under Novaseq 6000 platform and trimming ######"
echo "Task ${SLURM_ARRAY_TASK_ID} is processing ${sample}"

fastp \
-i ${dir}/${sample}_R1.fastq.gz \
-I ${dir}/${sample}_R2.fastq.gz \
-o ${outdir}/${sample}_R1.clean.fq.gz \
-O ${outdir}/${sample}_R2.clean.fq.gz \
-j ${outdir1}/${sample}.fastp.json \
-h ${outdir2}/${sample}.fastp.html
###########################################

###### QC Trimming file ######
echo "###### QC Trimming file ######"

fastqc ${outdir}/${sample}_R1.clean.fq.gz -o ${outdir3}
fastqc ${outdir}/${sample}_R2.clean.fq.gz -o ${outdir3}
###########################################




###### Multi-QC Trimming file as one reports ######
if [[ ${SLURM_ARRAY_TASK_ID} -eq 6 ]]; then
    echo "###### Multi-QC Trimming file as one reports ######"
    multiqc .${ourdir3}
else
    echo "Array ${SLURM_ARRAY_TASK_ID} and Job ${sample} is finished."
fi
print_separator
###########################################

echo "###### Finished Fastp Trimming ######"