#!/bin/bash
#SBATCH -J Fig_exprot
#SBATCH -o Fig_export_merge.out
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
echo "source activate deeptools"
source activate deeptools


path=/RAID1/working/R425/huanchi/atac/tomato
cd ${path}
echo "###### Current working directory: ${path} ######"
print_separator

###### bamCoverage ######
echo "###### bamCoverage ######"
# convert bam into bigwig
# file transformation should take 15-20 minutes
indir=${path}/merge_flt_bam

outdir=${path}/merge_flt_bam_RPKM_bw

mkdir ${outdir}

for file in ${indir}/*.flt.bam
do
    b=$(basename $file .flt.bam)
    bamCoverage -b ${file} \
        --binSize 10 \
        --normalizeUsing RPKM \
        -o ${outdir}/${b}.merge.bw 
done
print_separator
###########################################

indir=${path}/merge_flt_bam_RPKM_bw
indir2=${path}/bed
outdir=${path}/plot
mkdir ${outdir}



computeMatrix reference-point \
    -S ${indir}/ATAC_M82_0h.merge.bw ${indir}/ATAC_M82_1h.merge.bw ${indir}/ATAC_M82_6h.merge.bw \
    -R ${indir2}/C1.bed ${indir2}/C2.bed ${indir2}/C3.bed ${indir2}/C4.bed ${indir2}/C5.bed  \
    --referencePoint center \
    -b 3000 -a 3000 \
    --skipZeros -o ${outdir}/matrix_reference_point.scale.gz \
############


plotHeatmap \
 -m ${outdir}/matrix_reference_point.scale.gz\
 -out ${outdir}/heatmap_bed.pdf \
 --colorMap Blues \
 --heatmapHeight 5 \
 --heatmapWidth 2 \
 --refPointLabel "center" \
 --whatToShow 'heatmap and colorbar' \
 --regionsLabel C1 C2 C3 C4 C5 \
 --samplesLabel 0H 1H 6H \
############