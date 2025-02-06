#!/bin/bash
#SBATCH -J FEATURE
#SBATCH -o FEATURE.out
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
            mkdir ${path}/${dir}/feature
            mkdir ${path}/${dir}/feature/${alg}
            mkdir ${path}/${dir}/feature/${alg}/${cluster}
            for inp in ${path}/${dir}/${cluster}/*_distinct_pcc_enriched_kmer.txt
            do
                cp ${inp} ${path}/${dir}/feature/${alg}/${cluster}
            done
        done
    done
done
###########################################





###### summary the KMERS ######
path=/RAID1/working/R425/lavakau/ML_master
for dir in mp_acr mp_acr_strict mp_acr_promoter mp_acr_distal
do
    for i in {1..15}
    do
        for alg in SVM RF LogReg
        do
            cluster=C${i}
            indir=${path}/${dir}/feature
            indir2=${path}/${dir}/feature/${alg}
            indir3=${path}/${dir}/feature/${alg}/${cluster}
            python ${path}/get_kmer_overlap_source.py ${indir3} ${indir2}/${cluster}_${alg}_kmer.txt
        done
    done
done
###########################################