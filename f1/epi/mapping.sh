#!/bin/bash
#SBATCH -J Mapping
#SBATCH -o Mapping_%a.out
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 12  
#SBATCH --mem-per-cpu 1G
#SBATCH --array=1-28
#SBATCH --time 01-00:00


date
echo "###### Job submitted ######"
echo "Requiring 1 task with 15 cpus in total 15G Memory under this DAP-seq user"
echo

echo "###### Environment loading..... ######"
echo "module load Miniconda3/4.9.2, SAMtools/1.14-GCC-11.2.0, BEDTools/2.30.0-GCC-10.2.0, FastQC/0.11.9-Java-11, Bowtie2/2.4.4-GCC-11.2.0"
module load Miniconda3/4.9.2
module load SAMtools/1.14-GCC-11.2.0
module load BEDTools/2.30.0-GCC-10.2.0
module load FastQC/0.11.9-Java-11
module load Bowtie2/2.4.4-GCC-11.2.0
echo
echo "source activate base"
source activate base

path=/RAID1/working/R425/lavakau/chromatin
cd ${path};
echo "###### Current working directory: /RAID1/working/R425/lavakau/chromatin ######"
echo



###### set up array for parallel running ######
echo
echo "###### set up array for parallel running ######"
# Specify the path to the config file
config=${path}/range_meta.txt
# Extract the sample name for the current $SLURM_ARRAY_TASK_ID
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
###########################################

###### bowtie2 mapping ######
# GUIDELINE FROM https://www.cell.com/current-biology/fulltext/S0960-9822(19)31610-0?dgcid=raven_jbs_aip_email#secsectitle0075
# CUT&RUN and ChIp-seq reads were mapped to the Tak-1 v6.1 genome using Bowtie2 v2.1.0
# and further processed using Samtools v1.3 and Bedtools v2.17.0.

# MAPQ less than ten and duplicates were removed with Samtools.

# Inserts less than 150 bp were removed from further analyses, as these fragments are 
# sub-nucleosomal in size and likely represent noise when profiling histones and histone modifications.

# Deduplicated reads from 2-4 biological replicates were merged.

# We called peaks for chromatin marks using HOMER and 
# considered a gene associated with a mark if at least 50% of the gene length overlapped with peaks.
# settings: -style histone -size 250 -minDist 500.

# Bigwig files were made using deepTools.

echo
echo "###### bowtie2 mapping ######"
echo "please refer the manual of bowtie2 for building specie genome index file"
echo "Ex: bowtie2-build genome.fasta genome_index"
echo "here, the marchantia polymorpha index is in folder mar"
echo

cd ${path}
indir=${path}/clean_fq

outdir=${path}/bowtie2
outdir1=${path}/bowtie2/log
outdir2=${path}/bowtie2/sam_file
outdir3=${path}/bowtie2/bam_sort

mkdir ${outdir}
mkdir ${outdir1}
mkdir ${outdir2}
mkdir ${outdir3}

# -p number of alignment threads to launch
# -x Index filename prefix
# -1 -2 Files with #1/#2 mates, paired with files in <m2>/<m1> Could be gzip'ed
# -L length of seed substrings; must be >3, <32 (22, default)
# --mp Mismatch penalty. The sequence error rate is approximately: {.75 * exp[-log(4) * B/A]} (6, default)
# --dovetail Concordant reads (Properly aligned reads)

bowtie2 -p 5 --mp 6 -L 22 --dovetail -x ${path}/mar/mar \
-1 ${indir}/${sample}_R1.clean.fq.gz \
-2 ${indir}/${sample}_R2.clean.fq.gz -S ${outdir2}/${sample}.sorted.sam 2> ${outdir1}/${sample}.log

# filter out unmapped reads and sort mapped reads
# parsed to discard poorly mapped reads (more than three mismatches, mapping quality below 30 and suboptimal alignments)
# remove PCR duplicates before being converted to BAM file

samtools view -@ 20 -h -F 4 -q 30 -u -S ${outdir2}/${sample}.sorted.sam \
 | samtools sort -@ 20 -o ${outdir3}/${sample}.sorted.temp
samtools rmdup -s ${outdir3}/${sample}.sorted.temp ${outdir3}/${sample}.mapped_sorted.bam

# Filter BAM file to include only reads with insert size >= 150 bp
# we also keep the no filtered bam files
samtools view -h ${outdir3}/${sample}.mapped_sorted.bam | \
awk 'BEGIN {OFS="\t"} /^@/ {print} $9 >= 150 || $9 <= -150 {print}' | \
samtools view -Sb - > ${outdir3}/${sample}.mapped_sorted_150bp.bam

samtools index ${outdir3}/${sample}.mapped_sorted.bam
samtools index ${outdir3}/${sample}.mapped_sorted_150bp.bam
###########################################


echo
echo "Bowtie2 mapping done!"