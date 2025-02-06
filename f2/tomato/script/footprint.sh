#!/bin/bash
#SBATCH -J TOBIAS
#SBATCH -o TOBIAS.out
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
echo "source activate tobias"
source activate tobias

path=/RAID1/working/R425/lavakau/atac-seq/marchantia/20230626
cd ${path}
echo "###### Current working directory: ${path} ######"
print_separator





###### TOBIAS Corection ######
# echo "###### TOBIAS Corection ######"

indir=${path}

outdir=${path}/ATACorrect_DEP

mkdir ${outdir}

# for bam in ${indir}/*-merge.bam
# do
#     b=$(basename $bam .bam)
#     TOBIAS ATACorrect \
#         --core 12 \
#         --bam ${bam}  \
#         --genome MpTak_v6.1r2.genome.fasta \
#         --peaks DE_peaks_hclust_available.bed \
#         --outdir ${outdir}/${b}
# done
###########################################






###### TOBIAS FootprintScores ######
echo "###### TOBIAS FootprintScores ######"

indir2=${path}/ATACorrect_DEP

outdir2=${indir2}/FootprintScores

mkdir ${outdir2}

# for dir in ${indir2}/*/
# do
#     b=$(basename $dir)
#     TOBIAS FootprintScores \
#     --cores 20 \
#     --signal ${indir2}/${b}/*_corrected.bw \
#     --regions DE_peaks_hclust_available.bed \
#     --output ${outdir2}/${b}.footprint.bw
# done
###########################################






###### TOBIAS BINDetect ######
echo "###### TOBIAS BINDetect ######"

# available chromosome 'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chrU', 'chrV',
# 'unplaced-scaffold_010', 'unplaced-scaffold_056', 'unplaced-scaffold_078', 'unplaced-scaffold_086', 'unplaced-scaffold_098',
# 'unplaced-scaffold_131', 'unplaced-scaffold_145', 'unplaced-scaffold_162', 'unplaced-scaffold_190', 'unplaced-scaffold_202',
# 'unplaced-scaffold_221', 'unplaced-scaffold_250', 'unplaced-scaffold_257', 'unplaced-scaffold_267', 'unplaced-scaffold_279',
# 'unplaced-scaffold_281', 'unplaced-scaffold_315', 'unplaced-scaffold_349', 'unplaced-scaffold_362', 'unplaced-scaffold_406',
# 'unplaced-scaffold_431', 'unplaced-scaffold_432', 'unplaced-scaffold_433', 'unplaced-scaffold_435', 'unplaced-scaffold_436',
# 'unplaced-scaffold_439', 'unplaced-scaffold_440', 'unplaced-scaffold_441', 'unplaced-scaffold_445'
# therefore, 5 peaks were removed as DE_peaks_hclust_available.bed

outdir3=${indir2}/BINDetect
outdir4=${indir2}/BINDetect/promoter
outdir5=${indir2}/BINDetect/distal

# mkdir ${outdir3}
# mkdir ${outdir4}
# mkdir ${outdir5}

# promoter kmer
for a_full in ${outdir2}/*.footprint.bw
do
    if [[ "$a_full" == *-HS-merge.footprint.bw ]]; then
        a_g=$(basename $a_full -HS-merge.footprint.bw)
        echo "a_g is "${a_g}

        for b_full in ${outdir2}/*.footprint.bw
        do
            if [[ "$b_full" == *-CK-merge.footprint.bw ]]; then
                b_g=$(basename $b_full -CK-merge.footprint.bw)
                echo "b_g is "${b_g}

                if [[ "$a_g" == "$b_g" ]]; then
                    echo "Processing: $a_full vs $b_full"
                    TOBIAS BINDetect --motifs /RAID1/working/R425/lavakau/kmer/mp_acr_promoter/promoter_pfm.meme  \
                        --signals ${a_full} ${b_full} \
                        --genome MpTak_v6.1r2.genome.fasta \
                        --peaks DE_peaks_hclust_available2.bed \
                        --peak_header DE_peaks_hclust_header.txt \
                        --outdir ${outdir4}/${a_g}_HS_vs_${b_g}_CK \
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
# for a_full in ${outdir2}/*.footprint.bw
# do
#     if [[ "$a_full" == *-HS-merge.footprint.bw ]]; then
#         a_g=$(basename $a_full -HS-merge.footprint.bw)
#         echo "a_g is "${a_g}

#         for b_full in ${outdir2}/*.footprint.bw
#         do
#             if [[ "$b_full" == *-CK-merge.footprint.bw ]]; then
#                 b_g=$(basename $b_full -CK-merge.footprint.bw)
#                 echo "b_g is "${b_g}

#                 if [[ "$a_g" == "$b_g" ]]; then
#                     echo "Processing: $a_full vs $b_full"
#                     TOBIAS BINDetect --motifs /RAID1/working/R425/lavakau/kmer/mp_acr_distal/distal_pfm.meme  \
#                         --signals ${a_full} ${b_full} \
#                         --genome MpTak_v6.1r2.genome.fasta \
#                         --peaks DE_peaks_hclust_available2.bed \
#                         --peak_header DE_peaks_hclust_header.txt \
#                         --outdir ${outdir5}/${a_g}_HS_vs_${b_g}_CK \
#                         --cond_names HS CK \
#                         --cores 20
#                     echo "Finished processing : $a_full vs $b_full"
#                 else
#                     continue
#                 fi
#             else
#                 continue
#             fi
#         done
#     else
#         continue
#     fi
# done

echo "Job finished!"
###########################################