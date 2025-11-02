#!/bin/bash
#SBATCH -J KMER
#SBATCH -o KMER_FET_proximal_%a.out
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 5
#SBATCH --mem-per-cpu 2G
#SBATCH --array=1-4
#SBATCH --time 03-00:00

set -Eeuo pipefail

# Function to print a separator line
print_separator() {
    local cols
    cols=$( (tput cols 2>/dev/null) || echo 80 )
    printf '%*s\n' "${cols}" '' | tr ' ' -
}

date
echo "job submitted"
echo "requiring tasks under this KMER user"
echo

echo "Environment loading....."
echo "module load Miniconda3/4.9.2"
# Ensure PS1 is defined to avoid errors under 'set -u'
export PS1=${PS1:-noninteractive}
module load Miniconda3/4.9.2
# Activate Python env for fasta & FET steps
set +u
source activate ml
set -u
echo

cd /RAID1/working/R425/lavakau/pCRE
echo "Current working directory: /RAID1/working/R425/lavakau/pCRE"
print_separator

###### set up array for parallel running ######
# Specify the path to the config file
config=/RAID1/working/R425/lavakau/pCRE/mp_acr_new/array.txt

# Extract the sample name for the current $SLURM_ARRAY_TASK_ID
sample=$(awk -v ArrayTaskID="${SLURM_ARRAY_TASK_ID}" '($1==ArrayTaskID){print $2; exit}' "${config}")
if [[ -z "${sample}" ]]; then
    echo "ERROR: Failed to resolve sample for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} from ${config}" >&2
    exit 2
fi
###########################################

###### get fasta ######
dir=/RAID1/working/R425/lavakau/pCRE/mp_acr_new

# Convert all top-level .txt to .fa
for inp in "${dir}"/peak_proximal/*.txt; do
    [[ -e "$inp" ]] || continue
    python FastaManager_modified2.py \
        -f getseq2 \
        -fasta "${dir}/peak.coord.fa" \
        -name "${inp}"
done

# Convert per-cluster .txt to .fa (clusters 1..4)
for i in {1..4}; do
    cluster="cluster_${i}"
    for inp in "${dir}/peak_proximal/${cluster}"/*.txt; do
        [[ -e "$inp" ]] || continue
        python FastaManager_modified2.py \
            -f getseq2 \
            -fasta "${dir}/peak.coord.fa" \
            -name "${inp}"
    done
done
print_separator
###########################################

###### finding kmer with FET ######
echo "###### finding kmer with FET ######"
echo "Running: ${sample}"
indir="${dir}/peak_proximal"
outdir="${dir}/peak_proximal_fasta"

mkdir -p "${outdir}/${sample}"

# Positive set (combined) should be ${indir}/${sample}.txt.fa
pos_fa="${indir}/${sample}.txt.fa"
if [[ ! -s "${pos_fa}" ]]; then
    echo "ERROR: Positive FASTA not found: ${pos_fa}" >&2
    exit 3
fi

for input_neg in "${indir}/${sample}"/*.fa; do
    [[ -e "$input_neg" ]] || continue
    b=$(basename "${input_neg}" .txt.fa)
    python pCRE_Finding_FET.py \
        -pos "${pos_fa}" \
        -neg "${input_neg}" \
        -k "${dir}/6mer.txt" \
        -FDR Y \
        -save "${outdir}/${sample}/${b}.pcre"
done

# Switch to R env
set +u
conda deactivate || true
source activate R4.3.2
set -u

ends='.pcre_df_p0.01.txt'
Rscript "${dir}/pos_neg_conversion.R" "${outdir}" "${sample}" "${ends}"

outdir2="${dir}/mp_acr_proximal"
mkdir -p "${outdir2}/${sample}"

Rscript "${dir}/pcc_filtering.R" "${outdir}" "${sample}" "${outdir2}/${sample}"
print_separator
###########################################
echo "Finished!"
