#!/bin/bash
#SBATCH -J feature
#SBATCH -o feature.out
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 20 
#SBATCH --mem-per-cpu 1G
#SBATCH --time 03-00:00

date
echo "job submitted"
echo "requiring tasks under this KMER user"
echo

echo "Environment loading....."
echo "module load Miniconda3/4.9.2"
module load Miniconda3/4.9.2
source activate ml


path=/RAID1/working/R425/lavakau/ML_master
cd ${path}

# copy features
file="/RAID1/working/R425/lavakau/ML_master/sly_acr/dirs.txt"
out=${path}/feature_raw

mkdir ${out}

while IFS= read -r line
do
    echo "Processing directory: $line"
    indir=${path}/${line}

    for i in {1..5}
    do
        cluster="C"${i}
        mkdir ${out}/${line}
        mkdir ${out}/${line}/${cluster}

        for file in ${indir}/${cluster}/*_distinct_pcc_enriched_kmer.txt
        do
            cp ${file} ${out}/${line}/${cluster}
        done
    done
done < ${file}




# counts overlapping features
file="/RAID1/working/R425/lavakau/ML_master/sly_acr/dirs.txt"
out3=${path}/feature
indir=${path}/feature_raw

mkdir ${out3}

while IFS= read -r line
do
    echo "Processing directory: "${line}
    echo "Input directory: "${indir}/${line} 

    for i in {1..5}
    do
        cluster="C"${i}
        mkdir ${out3}/${line}

        python get_kmer_overlap_source.py ${indir}/${line}/${cluster} ${out3}/${line}/${cluster}_kmer.txt
    done
done < ${file}


echo "Finished"