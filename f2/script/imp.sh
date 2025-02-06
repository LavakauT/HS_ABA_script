#!/bin/bash
#SBATCH -J IMP
#SBATCH -o IMP.out
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 5  
#SBATCH --mem-per-cpu 2G
#SBATCH --time 03-00:00

date
echo "job submitted"
echo "requiring tasks under this KMER user"
echo

echo "Environment loading....."
echo "module load Miniconda3/4.9.2"
module load Miniconda3/4.9.2
source activate ml
echo

###### COPY IMP FILES ######
path=/RAID1/working/R425/lavakau/ML_master
for dir in mp_acr mp_acr_strict mp_acr_promoter mp_acr_distal
do
    for i in {1..15}
    do
        for alg in SVM RF LogReg
        do
            cluster=C${i}
            mkdir ${path}/${dir}/imp
            mkdir ${path}/${dir}/imp/${alg}
            mkdir ${path}/${dir}/imp/${alg}/${cluster}
            for inp in ${path}/${dir}/${cluster}/*_${alg}_imp
            do
                cp ${inp} ${path}/${dir}/imp/${alg}/${cluster}
            done
        done
    done
done
###########################################





###### summary the IMP ######
path=/RAID1/working/R425/lavakau/ML_master
for dir in mp_acr mp_acr_strict mp_acr_promoter mp_acr_distal
do
    for i in {1..15}
    do
        for alg in SVM RF LogReg
        do
            cluster=C${i}
            indir=${path}/${dir}/imp
            indir2=${path}/${dir}/imp/${alg}
            indir3=${path}/${dir}/imp/${alg}/${cluster}
            for inp in ${path}/${dir}/${cluster}/*_${alg}_imp
            do
                python ${path}/get_kmer-imp_overlap_source.py \
                    ${indir3} \
                    ${indir2}/${cluster}_${alg}_imp.txt
            done
        done
    done
done
###########################################