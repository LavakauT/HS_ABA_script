#!/bin/bash
#SBATCH -J merge_peak_featurecount
#SBATCH -o merge_peak_featurecount.out
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 20  
#SBATCH --mem-per-cpu 1G
#SBATCH --time 00-12:00

# Function to print a separator line
print_separator() {
    printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' -
}

date
echo "###### Job submitted ######"
echo "Requiring total $SLURM_NTASKS task(s) with $SLURM_CPUS_PER_TASK cpus in total $SLURM_MEM_PER_CPU Memory under this ATAC-seq user"
print_separator
echo "###### Environment loading..... ######"
echo "module load Miniconda3 FastQC SAMtools BEDTools Bowtie2"
module load Miniconda3/4.9.2
module load FastQC/0.11.9-Java-11
module load SAMtools/1.14-GCC-11.2.0
module load BEDTools/2.30.0-GCC-10.2.0
module load Bowtie2/2.4.4-GCC-11.2.0
echo "source activate featurecounts"
source activate featurecounts

path=/RAID1/working/R425/huanchi/atac/tomato
cd ${path}
echo "###### Current working directory: ${path} ######"
print_separator


indir1=${path}/narrowpeak_copy
indir2=${path}/merge_narrowpeak_copy
outdir1=${path}/flt_bam
outdir2=${path}/merge_flt_bam

###### merge peaks ######
# cat bed1.bed bed2.bed bedN.bed (...) | sort -k1,1 -k2,2n | bedtools merge -i - > merged.bed
# cp merged.bed to the folder of all bam files
# FeatureCounts with the bed file
# You can get the counts files to peak

cat ${indir1}/*.narrowPeak | sort -k1,1 -k2,2n | bedtools merge -i - > ${indir1}/merged.bed
cat ${indir2}/*.narrowPeak | sort -k1,1 -k2,2n | bedtools merge -i - > ${indir2}/merged.bed

cp ${indir1}/merged.bed ${outdir1}
cp ${indir2}/merged.bed ${outdir2}

# converting bed file to saf for featureCounts
# awk 'OFS="\t" {print $1"."$2"."$3, $1, $2, $3, "."}' merged.bed > merged.saf
awk 'OFS="\t" {print $1"."$2"."$3, $1, $2, $3, "."}' ${outdir1}/merged.bed > ${outdir1}/merged.saf
awk 'OFS="\t" {print $1"."$2"."$3, $1, $2, $3, "."}' ${outdir2}/merged.bed > ${outdir2}/merged.saf


# featureCounts -T 10 -a merged.saf -F SAF -o atac_merge_peak_counts.txt -p -B -C *.bam
# featureCounts -T 20 -a ${outdir1}/merged.saf -F SAF -o ${outdir1}/atac_peak_counts.txt -p -B -C ${outdir1}/*.flt.bam
# featureCounts -T 20 -a ${outdir2}/merged.saf -F SAF -o ${outdir2}/atac_merge_peak_counts.txt -p -B -C ${outdir2}/*.flt.bam
featureCounts -T 20 -a ${outdir2}/merged.saf -F SAF -o ${outdir1}/atac_merge_peak_counts.txt -p -B -C ${outdir1}/*.flt.bam
#######################