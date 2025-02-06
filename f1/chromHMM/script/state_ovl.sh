#!/bin/bash
#SBATCH -J Merge_state
#SBATCH -o Merge_state_%a.out
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 5  
#SBATCH --mem-per-cpu 2G
#SBATCH --array=1-15
#SBATCH --time 03-00:00

# Function to print a separator line
print_separator() {
    printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' -
}


date
echo "job submitted"
echo "requiring tasks under this KMER user"

echo "Environment loading....."
echo "module load Miniconda3/4.9.2 SAMtools/1.14-GCC-11.2.0 BEDTools/2.30.0-GCC-10.2.0"
module load Miniconda3/4.9.2
module load SAMtools/1.14-GCC-11.2.0
module load BEDTools/2.30.0-GCC-10.2.0
echo "Conda environment: deeptool"
source activate deeptool 


path=/RAID1/working/R425/lavakau/ChromHMM/marpol_chrom_state12
cd ${path}
echo "Current working directory: /RAID1/working/R425/lavakau/ChromHMM/marpol_chrom_state12"

###### set up array for parallel running ######
# Specify the path to the config file
config=/RAID1/working/R425/lavakau/pCRE/mp_acr/range.txt

# Extract the sample name for the current $SLURM_ARRAY_TASK_ID
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
print_separator
###########################################


###### overlapping peaks ######
a_bed_path=${path}/TP_bed
a_bed=${a_bed_path}/${sample}.bed
b_bed_path=${path}/bed

output=${a_bed_path}/${sample}
mkdir ${output}

for i in {1..12}
do
    state="E"${i}
    b_bed=${b_bed_path}/${state}.bed

    bedtools intersect -wa -a ${a_bed} -b ${b_bed} > ${output}/${sample}_${state}.bed
done


cd ${output}
wc -l *.bed > ${sample}_state_count.txt
###########################################