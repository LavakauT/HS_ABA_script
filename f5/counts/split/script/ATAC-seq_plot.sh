#!/bin/bash
#SBATCH -J ATAC-seq-plot
#SBATCH -o ATAC-seq_plot.out
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
echo "Requiring total $SLURM_NTASKS task(s) with $SLURM_CPUS_PER_TASK cpus in total $SLURM_MEM_PER_CPU Memory under this RNA-seq user"
print_separator
echo "###### Environment loading..... ######"
echo "module load Miniconda3 FastQC SAMtools BEDTools Bowtie2"
module load Miniconda3/4.9.2
module load FastQC/0.11.9-Java-11
module load SAMtools/1.14-GCC-11.2.0
module load BEDTools/2.30.0-GCC-10.2.0
module load Bowtie2/2.4.4-GCC-11.2.0
echo "source activate deeptool"
source activate deeptool

path=/RAID1/working/R425/lavakau/atac-seq/marchantia/20230928
cd ${path}
echo "###### Current working directory: ${path} ######"
print_separator


###### bamCoverage ######
echo "###### bamCoverage ######"
# convert bam into bigwig
# file transformation should take 15-20 minutes
indir=${path}/flt_bam

outdir=${path}/flt_bam_bw

mkdir ${outdir}

for file in ${indir}/*.flt.bam
do
    b=$(basename $file .flt.bam)
    bamCoverage -b ${file} \
        --binSize 10 \
        --normalizeUsing BPM \
        -o ${outdir}/${b}.flt.bw 
done
print_separator
###########################################


###### computeMatrix and Plot ######
echo "###### computeMatrix ######"
# matrix computation takes quite long time
# please update to deeptools3.5.2 for preventing dumpy error
indir=${path}/flt_bam

outdir2=${path}/matrix_scale_region
outdir3=${path}/matrix_reference_point

mkdir outdir2
mkdir outdir3

computeMatrix scale-regions \
    -S ${indir}/*.flt.bw \
    -R MpTak_v6.1_whole_gene_for_deeptools.bed \
    --regionBodyLength 2000 \
    --beforeRegionStartLength 3000 --afterRegionStartLength 3000 \
    --skipZeros -o ${outdir2}/matrix_scale_region.scale.gz

computeMatrix reference-point \
    -S ${indir}/*.flt.bw \
    -R MpTak_v6.1_whole_gene_for_deeptools.bed \
    --referencePoint TSS \
    -b 3000 -a 3000 \
    --skipZeros -o ${outdir3}/matrix_reference_point.reference.gz


plotProfile -m ${outdir2}/matrix_scale_region.scale.gz \
    -out ${outdir2}/scale_region.pdf \
    --perGroup;

plotProfile -m ${outdir2}/matrix_scale_region.scale.gz \
    -out ${outdir2}/scale_region_persample.pdf \
    --numPlotsPerRow 4;

plotProfile -m ${outdir3}/matrix_reference_point.reference.gz \
    -out ${outdir3}/reference_point_region.pdf \
    --perGroup;

plotProfile -m ${outdir3}/matrix_reference_point.reference.gz \
    -out ${outdir3}/reference_point_region_persample.pdf \
    --numPlotsPerRow 4;


plotHeatmap \
 -m ${outdir2}/matrix_scale_region.scale.gz\
 -out ${outdir2}/scale_ATAC_hp.pdf \
 --heatmapHeight 15  \
 --refPointLabel TSS \
 --regionsLabel promoters \
 --plotTitle 'ATAC-seq signal'

plotHeatmap \
 -m ${outdir3}/matrix_reference_point.reference.gz\
 -out ${outdir3}/refe_ATAC_hp.pdf \
 --heatmapHeight 15  \
 --refPointLabel TSS \
 --regionsLabel promoters \
 --plotTitle 'ATAC-seq signal'

print_separator
###########################################


###### multiBigwigSummary ######
echo "###### multiBigwigSummary ######"

outdir4=bw_summary

mkdir ${outdir4}

multiBigwigSummary bins -b ${indir}/*.flt.bw -o ${outdir4}/multibw_results.npz


plotCorrelation -in ${outdir4}/multibw_results.npz \
--corMethod spearman --skipZeros \
--whatToPlot scatterplot \
--plotTitle "Spearman Correlation" \
--removeOutliers \
--plotFile ${outdir4}/correlation_spearman_bwscore_scatterplot.pdf

plotCorrelation -in ${outdir4}/multibw_results.npz \
--corMethod spearman --skipZeros \
--whatToPlot heatmap \
--plotTitle "Spearman Correlation" \
--removeOutliers \
--plotNumbers \
--plotFile ${outdir4}/correlation_spearman_bwscore_heatmapplot.pdf
print_separator
###########################################

## .deb file needed in computeMatrix
# awk '$3=="gene"' Slycopersicum_691_ITAG4.0.gene.gff3 | awk 'BEGIN {OFS="\t"} {print $1,$4,$5,$9,$6,$7}' | sed -E 's/ID=(Mg.{7});.*\t\./\1\t\./' > Slycopersicum_691_ITAG4_gene_for_deeptools.bed

# collectinsertsizes
# java -jar picard.jar CollectInsertSizeMetrics -I rna-seq/WT-EtOH_2Aligned.sortedByCoord.out.bam -O rna-seq/WT-EtOH_2.txt -H rna-seq/WT-EtOH_2.pdf -M 0.5
# 
# java -jar picard.jar CollectInsertSizeMetrics -I rna-seq/Fy-3-mM-FAAligned.sortedByCoord.out.bam -O rna-seq/Fy-3-mM-FA.txt -H rna-seq/Fy-3-mM-FA.pdf -M 0.5