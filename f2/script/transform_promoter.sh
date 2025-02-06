#!/bin/bash
#SBATCH -J Transform
#SBATCH -o Transform_%a.out
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 12  
#SBATCH --mem-per-cpu 1G
#SBATCH --array=1-15
#SBATCH --time 03-00:00

date
echo "job submitted"
echo "requiring tasks under this KMER user"
echo

echo "Environment loading....."
echo "module load Miniconda3/4.9.2"
module load Miniconda3/4.9.2
source activate R4.3.2
echo

cd /RAID1/working/R425/lavakau/pCRE
echo "Current working directory: /RAID1/working/R425/lavakau/pCRE"
echo


###### set up array for parallel running ######
# Specify the path to the config file
config=/RAID1/working/R425/lavakau/pCRE/mp_acr/range.txt

# Extract the sample name for the current $SLURM_ARRAY_TASK_ID
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
###########################################



###### TRANSFORM ######
echo "###### TRANSFORM ######"
echo "Running: "${sample}
dir=/RAID1/working/R425/lavakau/pCRE/mp_acr
indir=${dir}/peak_promoter
outdir=${dir}/peak_promoter_fasta


ends='.pcre_df_p0.01.txt'
Rscript ${dir}/pos_neg_conversion.R ${outdir} ${sample} ${ends}


dir2=/RAID1/working/R425/lavakau/ML_master
outdir2=${dir2}/mp_acr_promoter

mkdir ${outdir2}
mkdir ${outdir2}/${sample}

Rscript ${dir}/pcc_filtering.R ${outdir} ${sample} ${outdir2}/${sample}
###########################################