#!/bin/bash
#SBATCH -J KMER
#SBATCH -o KMER_FET_promoter_%a.out
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
echo

echo "Environment loading....."
echo "module load Miniconda3/4.9.2"
module load Miniconda3/4.9.2
source activate ml
echo

cd /RAID1/working/R425/lavakau/pCRE
echo "Current working directory: /RAID1/working/R425/lavakau/pCRE"
print_separator


###### set up array for parallel running ######
# Specify the path to the config file
config=/RAID1/working/R425/lavakau/pCRE/mp_acr/range.txt

# Extract the sample name for the current $SLURM_ARRAY_TASK_ID
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
###########################################




###### get fasta ######
dir=/RAID1/working/R425/lavakau/pCRE/mp_acr
# for inp in ${dir}/peak_promoter2/*.txt
# do
#     python FastaManager_modified.py \
#         -f getseq2 \
#         -fasta ${dir}/peak.coord.fa \
#         -name ${inp}
# done

# for i in {1..15}
# do
#     cluster=C${i}
#     for inp in ${dir}/peak_promoter2/${cluster}/*.txt
#     do
#         python FastaManager_modified.py \
#             -f getseq2 \
#             -fasta ${dir}/peak.coord.fa \
#             -name ${inp}
#     done
# done
# print_separator
###########################################


###### finding kmer with FET ######
echo "###### finding kmer with FET ######"
echo "Running: "${sample}
dir=/RAID1/working/R425/lavakau/pCRE/mp_acr
indir=${dir}/peak_promoter2
outdir=${dir}/peak_promoter2_fasta

# mkdir ${outdir}
# mkdir ${outdir}/${sample}

# for input_neg in ${indir}/${sample}/*.fa
# do
#     b=$(basename ${input_neg} .txt.fa)
#     python pCRE_Finding_FET.py \
#         -pos ${indir}/${sample}.txt.fa  \
#         -neg ${input_neg} \
#         -k ${dir}/6mer.txt \
#         -FDR Y \
#         -save ${outdir}/${sample}/${b}.pcre
# done


conda deactivate
source activate R4.3.2
ends='.pcre_df_p0.01.txt'
Rscript ${dir}/pos_neg_conversion.R ${outdir} ${sample} ${ends}


dir2=/RAID1/working/R425/lavakau/ML_master
outdir2=${dir2}/mp_acr_promoter2

mkdir ${outdir2}
mkdir ${outdir2}/${sample}


Rscript ${dir}/pcc_filtering.R ${outdir} ${sample} ${outdir2}/${sample}
print_separator
###########################################
echo "Finished!"