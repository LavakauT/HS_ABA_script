#!/bin/bash
#SBATCH -J ATACsignal
#SBATCH -o ATACsignal__%A-%a.out
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 10  
#SBATCH --mem-per-cpu 1G
#SBATCH --array=1-10
#SBATCH --time 00-12:00

date
echo "job submitted"
echo "requiring 1 tasks under this ATAC-seq user"
cd /RAID1/working/R425/lavakau
echo "python 3.10.11 envs activate: "
source /RAID1/working/R425/lavakau/miniconda3/etc/profile.d/conda.sh
conda activate atac

echo " modules loading....."
ml SAMtools/1.14-GCC-11.2.0
ml BEDTools/2.30.0-GCC-10.2.0
ml Miniconda3/4.9.2
ml FastQC/0.11.9-Java-11
ml Bowtie2/2.4.4-GCC-11.2.0

echo  "bowtie2 mapping in sensitive model"
echo "please update to deeptools3.5.2 for dumpy error"
echo  "index:mar"

cd /RAID1/working/R425/lavakau/atac-seq/marchantia/20230626
# Specify the path to the config file
config=/RAID1/working/R425/lavakau/atac-seq/marchantia/20230626/array.txt

# Extract the sample name for the current $SLURM_ARRAY_TASK_ID
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)


computeMatrix scale-regions \
-S Tak1-CK.merge.bw hsfb-CK.merge.bw \
-R bed/${sample}.bed \
--regionBodyLength 2000 \
--scale 1 \
--samplesLabel Tak1 hsfb \
--beforeRegionStartLength 3000 --afterRegionStartLength 3000 \
--skipZeros -o signal/${sample}.msrs.gz

computeMatrix reference-point \
-S Tak1-CK.merge.bw hsfb-CK.merge.bw \
-R bed/${sample}.bed \
--referencePoint TSS \
--scale 1 \
--samplesLabel Tak1 hsfb \
-b 3000 -a 3000 \
--skipZeros -o signal/${sample}.mrpr.gz

plotProfile -m signal/${sample}.msrs.gz -out signal/msrs/${sample}.msrs.pdf --colors black red --perGroup;
plotProfile -m signal/${sample}.msrs.gz -out signal/per_msrs/${sample}.persample.msrs.pdf --colors black red --numPlotsPerRow 2

plotProfile -m signal/${sample}.mrpr.gz -out signal/mrpr/${sample}.mrpr.pdf --colors black red --perGroup;
plotProfile -m signal/${sample}.mrpr.gz -out signal/per_mrpr/${sample}.persample.mrpr.pdf --colors black red --numPlotsPerRow 2

plotHeatmap \
 -m signal/${sample}.msrs.gz\
 -out signal/${sample}.scale_ATAC_hp.pdf \
 --heatmapHeight 15  \
 --refPointLabel TSS \
 --regionsLabel signal \
 --plotTitle ${sample}_ATAC-seq_signal \

plotHeatmap \
 -m signal/${sample}.mrpr.gz\
 -out signal/${sample}.refe_ATAC_hp.pdf \
 --heatmapHeight 15  \
 --refPointLabel TSS \
 --regionsLabel signal \
 --plotTitle ${sample}_ATAC-seq_signal \
