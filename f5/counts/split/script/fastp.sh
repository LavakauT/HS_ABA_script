#!/bin/bash
#SBATCH -J FastpQC
#SBATCH -o Fastp_mp.out
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

path=/RAID1/working/R425/lavakau/atac-seq/marchantia/20230928
cd ${path}
echo "###### Current working directory: ${path}/rowData ######"
echo


###### set up output folder and sequencing numbers ######
dir=${path}/rawData
out=${path}/fastp
outdir=${path}/fastp/clean_fq
outdir1=${path}/fastp/clean_fq_json
outdir2=${path}/fastp/clean_fq_html
outdir3=${path}/fastp/clean_fq_qc
seq=22757NLT3_L4

mkdir ${out}
mkdir ${outdir}
mkdir ${outdir1}
mkdir ${outdir2}
mkdir ${outdir3}
###########################################


###### set up array for parallel running ######
num_tasks=$SLURM_NTASKS
# Specify the path to the config file
config=${path}/rowData/range.txt
###########################################

for ((i=1; i<num_tasks; i++))
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


    echo "Task $index has finished processing $sample"
done

# wait for finishing all the samples
wait

echo "###### Finished Fastp Trimming ######"