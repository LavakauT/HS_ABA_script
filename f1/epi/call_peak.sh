#!/bin/bash
#SBATCH -J DAP-CP
#SBATCH -o DAP-CP2_mp_%A-%a.out
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 15  
#SBATCH --mem-per-cpu 1G
#SBATCH --array=1-6
#SBATCH --time 01-00:00


date
echo "###### Job submitted ######"
echo "Requiring 1 task with 15 cpus in total 15G Memory under this DAP-seq user"
echo

echo "###### Environment loading..... ######"
echo "module load Miniconda3/4.9.2, SAMtools/1.14-GCC-11.2.0, BEDTools/2.30.0-GCC-10.2.0, FastQC/0.11.9-Java-11, Bowtie2/2.4.4-GCC-11.2.0"
echo "module load Java/11.0.2 GEMmotif/3.4-Java-11"
module load Miniconda3/4.9.2
module load SAMtools/1.14-GCC-11.2.0
module load BEDTools/2.30.0-GCC-10.2.0
module load FastQC/0.11.9-Java-11
module load Bowtie2/2.4.4-GCC-11.2.0
module load Java/11.0.2
module load GEMmotif/3.4-Java-11 
echo "To execute GEM run: java -jar $EBROOTGEMMOTIF/gem.jar"
echo
echo "source activate macs2"
source activate macs2

cd /RAID1/working/R425/huanchi/dap
echo
echo "###### Current working directory: /RAID1/working/R425/huanchi/dap ######"
echo

cd /home;
cd /RAID1/working/R425/huanchi/dap


###### set up array for parallel running ######
echo
echo "###### set up array for parallel running ######"
# Specify the path to the config file
config=/RAID1/working/R425/huanchi/dap/hsfb_test/range.txt
# Extract the sample name for the current $SLURM_ARRAY_TASK_ID
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
###########################################


# ###### macs2 peak calling ######
# echo
# echo "###### macs2 peak calling ######"
# echo

# indir=/RAID1/working/R425/huanchi/dap/hsfb_test/bowtie2/sam_sort

# outdir=/RAID1/working/R425/huanchi/dap/hsfb_test/bowtie2/macs2

# mkdir ${outdir}
# # -t: ChIP-seq treatment file. If multiple files are given
# # as '-t A B C', then they will all be read and pooled
# # together. REQUIRED.

# # -c: Control file. If multiple files are given as '-c A B
# # C', they will be pooled to estimate ChIP-seq
# # background noise.

# # -g: Effective genome size. for mp is 1.1e8 (larger than Arabidopsis 1.6e7)
# # -B: Store all details
# # -m: Select the regions within MFOLD range of high-
# # confidence enrichment ratio against background to
# # build model.
# # -q: Minimum FDR (q-value) cutoff for peak detection.

# # background file are Y100, Y50, Y30, E30, No100 

# # 20240513:
# # WARNING @ Mon, 13 May 2024 12:29:09: Too few paired peaks (34) so I can not build the model! Broader your MFOLD range parameter may erase this error. If it still can't build the model, we suggest to use --nomodel and --extsize 147 or other fixed number instead. 
# # WARNING @ Mon, 13 May 2024 12:29:09: Process for pairing-model is terminated! 
# # original:${indir}/Y100.mapped_sorted.bam ${indir}/Y50.mapped_sorted.bam ${indir}/Y30.mapped_sorted.bam ${indir}/E30.mapped_sorted.bam ${indir}/No100.mapped_sorted.bam
# # modified: ${indir}/Y100.mapped_sorted.bam ${indir}/E30.mapped_sorted.bam ${indir}/No100.mapped_sorted.bam

# macs2 callpeak -t ${indir}/${sample}.mapped_sorted.bam \
#     -c ${indir}/Y100.mapped_sorted.bam ${indir}/E30.mapped_sorted.bam ${indir}/No100.mapped_sorted.bam \
#     -f BAM --outdir ${outdir} -g 1.6e7 -n ${sample} -B -q 0.05 -m 2 50
# ###########################################


# ###### GEM peak calling and motif model creation ######
# # mkdir /RAID1/working/R425/huanchi/dap/genome
# # put Read_Distribution_default.txt and Marpol_genome_size.txt into genome folder
# # mkdir /RAID1/working/R425/huanchi/dap/genome/mar 
# # put genome fasta file into mar folder
# echo
# echo "###### GEM peak calling and motif model creation ######"
# echo

# indir=/RAID1/working/R425/huanchi/dap/hsfb_test/bowtie2/sam_sort

# outdir1=/RAID1/working/R425/huanchi/dap/hsfb_test/bowtie2/gem

# mkdir ${outdir1}

# # GEM tool for peak calling and motif model creation - use BAM file
# # force q value to be 0.001, q = 3, 0.01, q =2
# # Only to use GPS data here

# ## The p-value of the motif occurence. 
# # The p-value is the probability of a random sequence of the same 
# # length as the motif matching that position of the sequence with a 
# # score at least as good.

# ## The q-vlavlue of the motif occurence. 
# # The q-value is the estimated false discovery rate if the occurrence 
# # is accepted as significant.

# ## The width (bp) to smooth the read distribution. 
# # If it is set to -1, there will be no smoothing (default=30). 
# # For ChIP-exo, use a smaller number (say 3 or 5) to achieve higher spatial accuracy.

# ## Fold (IP/Control) cutoff to filter predicted events (default=3)
# ## Depending on the sequencing depth, the sample read
# ##################################### Proc. Natl Acad. Sci. USA (2003) 100:9440–9445
# # --t: maximum numbers of threads
# # --d: read spatial distribution file (supported by GEM)
# # --g: genome chrom.sizes file with chr name/length pairs (chrom \t total_length)
# # please generat genome chrom.sizes filee Marchantia polymorpha
# # --genome: the path to the genome sequence directory, for motif finding
# # --expt: aligned read file for experimental treatment
# # --ctrl: aligned reads file for control
# # --k: length of the k-mer for motif finding, use --k or (--kmin & --kmax)
# # -k_seqs: number of binding events to use for motif discovery (default=5000)


# java -jar -Xmx15G $EBROOTGEMMOTIF/gem.jar --t 8 \
#         --f BAM \
#         --d /RAID1/working/R425/huanchi/dap/genome/Read_Distribution_default.txt \
#         --g /RAID1/working/R425/huanchi/dap/genome/Marpol_genome_size.txt \
#         --genome /RAID1/working/R425/huanchi/dap/genome \
#         --expt ${indir}/${sample}.mapped_sorted.bam \
#         --ctrl ${indir}/Y100.mapped_sorted.bam ${indir}/Y50.mapped_sorted.bam ${indir}/Y30.mapped_sorted.ba ${indir}/No100.mapped_sorted.bam ${indir}/E30.mapped_sorted.bam \
#         --outBED \
#         --out ${outdir1}/${sample} \
#         --k_min 6 --kmax 20 --k_seqs 600 --k_neg_dinu_shuffle
# ###########################################


# echo
# date
# echo "Macs2 and GEM peak calling done!"



# REMEMBER TO CHANGE LOG NAME
###### macs2 peak calling ######
echo
echo "###### macs2 peak calling ######"
echo

indir=/RAID1/working/R425/huanchi/dap/hsfb_test/bowtie2/sam_sort_trimo

outdir=/RAID1/working/R425/huanchi/dap/hsfb_test/bowtie2/macs2_trimo

mkdir ${outdir}
# -t: ChIP-seq treatment file. If multiple files are given
# as '-t A B C', then they will all be read and pooled
# together. REQUIRED.

# -c: Control file. If multiple files are given as '-c A B
# C', they will be pooled to estimate ChIP-seq
# background noise.

# -g: Effective genome size. for mp is 1.1e8 (larger than Arabidopsis 1.6e7)
# -B: Store all details
# -m: Select the regions within MFOLD range of high-
# confidence enrichment ratio against background to
# build model.
# -q: Minimum FDR (q-value) cutoff for peak detection.

# background file are Y100, Y50, Y30, E30, No100 

# 20240513:
# WARNING @ Mon, 13 May 2024 12:29:09: Too few paired peaks (34) so I can not build the model! Broader your MFOLD range parameter may erase this error. If it still can't build the model, we suggest to use --nomodel and --extsize 147 or other fixed number instead. 
# WARNING @ Mon, 13 May 2024 12:29:09: Process for pairing-model is terminated! 
# original:${indir}/Y100.mapped_sorted.bam ${indir}/Y50.mapped_sorted.bam ${indir}/Y30.mapped_sorted.bam ${indir}/E30.mapped_sorted.bam ${indir}/No100.mapped_sorted.bam
# modified: ${indir}/Y100.mapped_sorted.bam ${indir}/E30.mapped_sorted.bam ${indir}/No100.mapped_sorted.bam

macs2 callpeak -t ${indir}/${sample}.mapped_sorted.bam \
    -c ${indir}/Y100.mapped_sorted.bam ${indir}/E30.mapped_sorted.bam ${indir}/No100.mapped_sorted.bam \
    -f BAM --outdir ${outdir} -g 1.6e7 -n ${sample} -B -q 0.05 -m 2 50
###########################################


###### GEM peak calling and motif model creation ######
# mkdir /RAID1/working/R425/huanchi/dap/genome
# put Read_Distribution_default.txt and Marpol_genome_size.txt into genome folder
# mkdir /RAID1/working/R425/huanchi/dap/genome/mar 
# put genome fasta file into mar folder
echo
echo "###### GEM peak calling and motif model creation ######"
echo

indir=/RAID1/working/R425/huanchi/dap/hsfb_test/bowtie2/sam_sort

outdir1=/RAID1/working/R425/huanchi/dap/hsfb_test/bowtie2/gem

mkdir ${outdir1}

# GEM tool for peak calling and motif model creation - use BAM file
# force q value to be 0.001, q = 3, 0.01, q =2
# Only to use GPS data here

## The p-value of the motif occurence. 
# The p-value is the probability of a random sequence of the same 
# length as the motif matching that position of the sequence with a 
# score at least as good.

## The q-vlavlue of the motif occurence. 
# The q-value is the estimated false discovery rate if the occurrence 
# is accepted as significant.

## The width (bp) to smooth the read distribution. 
# If it is set to -1, there will be no smoothing (default=30). 
# For ChIP-exo, use a smaller number (say 3 or 5) to achieve higher spatial accuracy.

## Fold (IP/Control) cutoff to filter predicted events (default=3)
## Depending on the sequencing depth, the sample read
##################################### Proc. Natl Acad. Sci. USA (2003) 100:9440–9445
# --t: maximum numbers of threads
# --d: read spatial distribution file (supported by GEM)
# --g: genome chrom.sizes file with chr name/length pairs (chrom \t total_length)
# please generat genome chrom.sizes filee Marchantia polymorpha
# --genome: the path to the genome sequence directory, for motif finding
# --expt: aligned read file for experimental treatment
# --ctrl: aligned reads file for control
# --k: length of the k-mer for motif finding, use --k or (--kmin & --kmax)
# -k_seqs: number of binding events to use for motif discovery (default=5000)


java -jar -Xmx15G $EBROOTGEMMOTIF/gem.jar --t 8 \
        --f BAM \
        --d /RAID1/working/R425/huanchi/dap/genome/Read_Distribution_default.txt \
        --g /RAID1/working/R425/huanchi/dap/genome/Marpol_genome_size.txt \
        --genome /RAID1/working/R425/huanchi/dap/genome \
        --expt ${indir}/${sample}.mapped_sorted.bam \
        --ctrl ${indir}/Y100.mapped_sorted.bam ${indir}/Y50.mapped_sorted.bam ${indir}/Y30.mapped_sorted.ba ${indir}/No100.mapped_sorted.bam ${indir}/E30.mapped_sorted.bam \
        --outBED \
        --out ${outdir1}/${sample} \
        --k_min 6 --kmax 20 --k_seqs 600 --k_neg_dinu_shuffle
###########################################


echo
date
echo "Macs2 and GEM peak calling done!"