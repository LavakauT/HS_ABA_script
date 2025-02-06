#!/bin/bash
#SBATCH -J TOBIAS_pd
#SBATCH -o TOBIAS_pd.out
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 20  
#SBATCH --mem-per-cpu 1G
#SBATCH --time 03-00:00


# Function to print a separator line
print_separator() {
    printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' -
}

date
echo "###### Job submitted ######"
echo "Requiring total $SLURM_NTASKS task(s) with $SLURM_CPUS_PER_TASK cpus in total $SLURM_MEM_PER_CPU Memory under this ATAC-seq user"
print_separator
echo "###### Environment loading..... ######"
echo "module load Miniconda3"
module load Miniconda3/4.9.2
echo "source activate base"
source activate base

path=/RAID1/working/R425/huanchi/atac/tomato
cd ${path}
echo "###### Current working directory: ${path} ######"
print_separator



###### TOBIAS Corection ######
# echo "###### TOBIAS Corection ######"

indir=${path}/merge_flt_bam

outdir=${path}/ATACorrect_DEP

mkdir ${outdir}

for bam in ${indir}/*.flt.bam
do
    b=$(basename $bam .flt.bam)
    TOBIAS ATACorrect \
        --core 20 \
        --bam ${bam}  \
        --genome Slycopersicum_796_ITAG5.0.fa \
        --peaks DE_peaks_hclust_available.bed \
        --outdir ${outdir}/${b}
done
###########################################






###### TOBIAS FootprintScores ######
echo "###### TOBIAS FootprintScores ######"

indir2=${path}/ATACorrect_DEP

outdir2=${indir2}/FootprintScores

mkdir ${outdir2}

for dir in ${indir2}/*/
do
    b=$(basename $dir)
    TOBIAS FootprintScores \
    --cores 20 \
    --signal ${indir2}/${b}/*_corrected.bw \
    --regions DE_peaks_hclust_available.bed \
    --output ${outdir2}/${b}.footprint.bw
done
###########################################






###### TOBIAS BINDetect ######
echo "###### TOBIAS BINDetect ######"

outdir3=${indir2}/BINDetect_pd
outdir4=${indir2}/BINDetect_pd/proximal
outdir5=${indir2}/BINDetect_pd/distal

mkdir ${outdir3}
mkdir ${outdir4}
mkdir ${outdir5}


# promoter kmer
for a_full in ${outdir2}/*.footprint.bw
do
    if [[ "$a_full" != *_0h.footprint.bw ]]; then
        a_g=$(basename $a_full .footprint.bw)
        echo "a_g is "${a_g}

        for b_full in ${outdir2}/*.footprint.bw
        do
            if [[ "$b_full" == *_0h.footprint.bw ]]; then
                b_g=$(basename $b_full .footprint.bw)
                echo "b_g is "${b_g}

                if [[ "$a_g" != "$b_g" ]]; then
                    echo "Processing: $a_full vs $b_full"
                    TOBIAS BINDetect --motifs /RAID1/working/R425/lavakau/kmer/sly_acr_pro/proximal_pfm.meme  \
                        --signals ${a_full} ${b_full} \
                        --genome Slycopersicum_796_ITAG5.0.fa \
                        --peaks DE_peaks_hclust_available_proximal.bed \
                        --peak_header DE_peaks_hclust_header.txt \
                        --outdir ${outdir4}/${a_g}_vs_${b_g} \
                        --cond_names HS CK \
                        --cores 20
                    echo "Finished processing : $a_full vs $b_full"
                else
                    continue
                fi
            else
                continue
            fi
        done
    else
        continue
    fi
done


# distal kmer
for a_full in ${outdir2}/*.footprint.bw
do
    if [[ "$a_full" != *_0h.footprint.bw ]]; then
        a_g=$(basename $a_full .footprint.bw)
        echo "a_g is "${a_g}

        for b_full in ${outdir2}/*.footprint.bw
        do
            if [[ "$b_full" == *_0h.footprint.bw ]]; then
                b_g=$(basename $b_full .footprint.bw)
                echo "b_g is "${b_g}

                if [[ "$a_g" != "$b_g" ]]; then
                    echo "Processing: $a_full vs $b_full"
                    TOBIAS BINDetect --motifs /RAID1/working/R425/lavakau/kmer/sly_acr_dis/distal_pfm.meme \
                        --signals ${a_full} ${b_full} \
                        --genome Slycopersicum_796_ITAG5.0.fa \
                        --peaks DE_peaks_hclust_available_distal.bed \
                        --peak_header DE_peaks_hclust_header.txt \
                        --outdir ${outdir5}/${a_g}_vs_${b_g} \
                        --cond_names HS CK \
                        --cores 20
                    echo "Finished processing : $a_full vs $b_full"
                else
                    continue
                fi
            else
                continue
            fi
        done
    else
        continue
    fi
done


echo "Job finished!"
###########################################