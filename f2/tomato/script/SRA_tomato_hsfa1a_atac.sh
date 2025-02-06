#!/bin/sh
#SBATCH -J SRA
#SBATCH -o SRA_%a.out
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 5  
#SBATCH --mem-per-cpu 2G
#SBATCH --array=1-6
#SBATCH --time 01-00:00

###### set up array for parallel running ######
# Specify the path to the config file
config=/RAID1/working/R425/lavakau/chromatin/range.txt

# Extract the sample name for the current $SLURM_ARRAY_TASK_ID
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
output=/RAID1/working/R425/huanchi/atac/tomato/rawData
###########################################


# parameter:
# -O: option to specify the output directory
# -e <threads>: Number of threads to use.
# --split-files: Split the output into separate files for paired-end reads.
# --gzip: Compress the output files using gzip.




echo ">>> Job submitted <<<"
echo ">>> Platform: sratoolkit.3.0.6-centos_linux64"

cd /RAID1/working/R425/huanchi/sratoolkit.3.0.6-centos_linux64/bin

./fasterq-dump -O ${output} -e 8 --split-files ${sample}

# rememeber to compress the files after downloading all the files
echo ">>> Job finished <<<"