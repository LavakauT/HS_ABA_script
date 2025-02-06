#!/bin/bash

# conda config --env --set subdir osx-64
# conda install bioconda::bedtools
# conda install bioconda::bioawk
# conda install conda-forge::coreutils

conda activate bedtool
p=/Users/user/Desktop/tomato_atac/merge_narrowpeak_copy
out1=/Users/user/Desktop/tomato_atac/merge_narrowpeak_copy/bed

cd ${p}
mkdir ${out1}
mkdir ${out1}/intersect
mkdir ${out1}/gain
mkdir ${out1}/loss



for h in 1h 6h
do
    # assign data
    data1=${p}/ATAC_M82_0h_peaks.narrowPeak
    data2=${p}/ATAC_M82_${h}_peaks.narrowPeak

    # make bed
    bioawk -c bed '{print $1, $2, $3}' ${data1} > ${out1}/ATAC_M82_0h.bed
    bioawk -c bed '{print $1, $2, $3}' ${data2} > ${out1}/ATAC_M82_${h}.bed

    # assign files
    file1=${out1}/ATAC_M82_0h.bed
    file2=${out1}/ATAC_M82_${h}.bed

    # overlapping peaks
    # -a is heat stress
    # -b is control
    bedtools intersect -a ${file2}\
        -b ${file1} -wa > ${out1}/intersect/intersect_ATAC_M82_${h}_vs_0h.bed
    
    # gain
    bedtools intersect -a ${file2}\
        -b ${file1} -v > ${out1}/gain/gain_ATAC_M82_${h}_vs_0h.bed
    
    # loss
    bedtools intersect -a ${file1}\
        -b ${file2} -v > ${out1}/loss/loss_ATAC_M82_${h}_vs_0h.bed
done



for h in 6h
do
    # assign data
    data1=${p}/ATAC_M82_1h_peaks.narrowPeak
    data2=${p}/ATAC_M82_${h}_peaks.narrowPeak

    # make bed
    bioawk -c bed '{print $1, $2, $3}' ${data1} > ${out1}/ATAC_M82_1h.bed
    bioawk -c bed '{print $1, $2, $3}' ${data2} > ${out1}/ATAC_M82_${h}.bed

    # assign files
    file1=${out1}/ATAC_M82_1h.bed
    file2=${out1}/ATAC_M82_${h}.bed

    # overlapping peaks
    # -a is heat stress
    # -b is control
    bedtools intersect -a ${file2}\
        -b ${file1} -wa > ${out1}/intersect/intersect_ATAC_M82_${h}_vs_1h.bed
    
    # gain
    bedtools intersect -a ${file2}\
        -b ${file1} -v > ${out1}/gain/gain_ATAC_M82_${h}_vs_1h.bed
    
    # loss
    bedtools intersect -a ${file1}\
        -b ${file2} -v > ${out1}/loss/loss_ATAC_M82_${h}_vs_1h.bed
done


wc -l ${p}/*.narrowPeak > ${p}/all_peaks.txt
wc -l ${out1}/gain/gain*.bed > ${p}/gain_peaks.txt
wc -l ${out1}/loss/loss*.bed > ${p}/loss_peaks.txt
wc -l ${out1}/intersect/intersect*.bed > ${p}/intersect_peaks.txt
#############################################

echo "Finished all the bedtools intersection analysis!"