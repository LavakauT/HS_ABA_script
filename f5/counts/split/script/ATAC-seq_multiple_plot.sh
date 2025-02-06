#!/bin/bash
#SBATCH -J Fig_exprot
#SBATCH -o Fig_export.out
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16  
#SBATCH --mem-per-cpu 1G
#SBATCH --time 00-12:00

date
echo "job submitted"
echo "requiring 1 tasks under this ATAC-seq user"
cd /RAID1/working/R425/lavakau
echo "python 3.10.11 envs activate: tobias"
source /RAID1/working/R425/lavakau/miniconda3/etc/profile.d/conda.sh
conda activate atac

echo " modules loading....."
ml SAMtools/1.14-GCC-11.2.0
ml BEDTools/2.30.0-GCC-10.2.0
ml Miniconda3/4.9.2
ml FastQC/0.11.9-Java-11
ml Bowtie2/2.4.4-GCC-11.2.0

echo  "bowtie2 mapping in sensitive model"
echo "please update to deeptools3.5.2 for removal of dumpy error"
echo  "index:mar (mar: Marchantia polymorpha)"

cd /RAID1/working/R425/lavakau/atac-seq/marchantia/20230928

# computeMatrix scale-regions \
# -S *.merge.bw \
# -R MpTak_v6.1_whole_gene_for_deeptools.bed \
# --regionBodyLength 2000 \
# --beforeRegionStartLength 3000 --afterRegionStartLength 3000 \
# --skipZeros -o matrix_scale_region.scale.all.gz

# multiple beds:scale-regions
# computeMatrix scale-regions \
# -S *.merge.bw \
# -R bed_mp/C1.bed bed_mp/C3.bed bed_mp/C4.bed bed_mp/C5.bed \
# --regionBodyLength 2000 \
# --beforeRegionStartLength 3000 --afterRegionStartLength 3000 \
# --skipZeros -o matrix_scale_region.scale.1.gz

# multiple beds:scale-regions
# computeMatrix scale-regions \
# -S *.merge.bw \
# -R bed_mp/C7.bed bed_mp/C8.bed bed_mp/C14.bed bed_mp/C15.bed \
# --regionBodyLength 2000 \
# --beforeRegionStartLength 3000 --afterRegionStartLength 3000 \
# --skipZeros -o matrix_scale_region.scale.2.gz

# computeMatrix reference-point \
# -S *.merge.bw \
# -R MpTak_v6.1_whole_gene_for_deeptools.bed \
# --referencePoint TSS \
# -b 3000 -a 3000 \
# --skipZeros -o matrix_reference_point.reference.all.gz

# multiple beds:reference-point
# computeMatrix reference-point \
# -S *.merge.bw \
# -R bed_mp/C1.bed bed_mp/C3.bed bed_mp/C4.bed bed_mp/C5.bed \
# --referencePoint TSS \
# -b 3000 -a 3000 \
# --skipZeros -o matrix_reference_point.reference.1.gz

# multiple beds:reference-point
# computeMatrix reference-point \
# -S *.merge.bw \
# -R bed_mp/C7.bed bed_mp/C8.bed bed_mp/C14.bed bed_mp/C15.bed \
# --referencePoint TSS \
# -b 3000 -a 3000 \
# --skipZeros -o matrix_reference_point.reference.2.gz

# plotHeatmap \
#  -m matrix_scale_region.scale.all.gz\
#  -out scale_ATAC_hp.pdf \
#  --colorList '#ffffff,orange,#000000' \
#  --heatmapHeight 15  \
#  --refPointLabel TSS \
#  --regionsLabel "All DEGs" \
#  --plotTitle 'ATAC-seq signal' 

# plotHeatmap \
#  -m matrix_reference_point.reference.all.gz\
#  -out refe_ATAC_hp.pdf \
#  --colorList '#ffffff,orange,#000000' \
#  --heatmapHeight 15  \
#  --refPointLabel TSS \
#  --regionsLabel "All DEGs" \
#  --plotTitle 'ATAC-seq signal' 

# plotHeatmap \
#  -m matrix_scale_region.scale.2.gz\
#  -out scale_ATAC_aba_C1-C5.pdf \
#  --colorList '#ffffff,orange,#000000' \
#  --heatmapHeight 15  \
#  --refPointLabel TSS \
#  --regionsLabel C1 C3 C4 C5 \
#  --plotTitle 'ATAC-seq signal' 

# plotHeatmap \
#  -m matrix_reference_point.reference.2.gz\
#  -out refe_ATAC_aba_C1-C5.pdf \
#  --colorList '#ffffff,orange,#000000' \
#  --heatmapHeight 15  \
#  --refPointLabel TSS \
#  --regionsLabel C1 C3 C4 C5 \
#  --plotTitle 'ATAC-seq signal'

# plotHeatmap \
#  -m matrix_scale_region.scale.2.gz\
#  -out scale_ATAC_aba_C7-C15.pdf \
#  --colorList '#ffffff,orange,#000000' \
#  --heatmapHeight 15  \
#  --refPointLabel TSS \
#  --regionsLabel C7 C8 C14 C15 \
#  --plotTitle 'ATAC-seq signal' 

# plotHeatmap \
#  -m matrix_reference_point.reference.2.gz\
#  -out refe_ATAC_aba_C7-C15.pdf \
#  --colorList '#ffffff,orange,#000000' \
#  --heatmapHeight 15  \
#  --refPointLabel TSS \
#  --regionsLabel C7 C8 C14 C15 \
#  --plotTitle 'ATAC-seq signal'