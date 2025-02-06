#!/bin/bash
#SBATCH -J Compress
#SBATCH -o Compress_%a.out
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 5  
#SBATCH --mem-per-cpu 2G
#SBATCH --array=1-28
#SBATCH --time 00-12:00


date
echo "###### Job submitted ######"
echo "Requiring 1 task with 5 cpus in total 10G Memory under this ChIP-seq user"
echo


path=/RAID1/working/R425/lavakau/chromatin
cd ${path};


###### set up array for parallel running ######
# Specify the path to the config file
config=${path}/range_meta.txt
# Extract the sample name for the current $SLURM_ARRAY_TASK_ID
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
echo "sample name: "${sample}
input_dir=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)
echo "sample directory: "${input_dir}
###########################################

# make output directory
output_dir=${path}/rawData
mkdir ${output_dir}

# copy the files to output directory
# there are paired-end files (R1, R2) 
cp ${input_dir}/${sample}_R1.fastq ${output_dir}
cp ${input_dir}/${sample}_R2.fastq ${output_dir}


# compress the copied files
# there are paired-end files (R1, R2) 
gzip ${output_dir}/${sample}_R1.fastq
gzip ${output_dir}/${sample}_R2.fastq

echo "Finished!"