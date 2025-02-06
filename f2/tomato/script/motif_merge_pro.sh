#!/bin/bash
#SBATCH -J MMer
#SBATCH -o MMer_pro.out
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 20  
#SBATCH --mem-per-cpu 1G
#SBATCH --time 00-12:00

echo "Job submitted"
echo "Requiring tasks under this KMER user"

echo ">>>ENVIRONMENT LAODING<<<"
echo " Module load Miniconda3/4.9.2 "
echo " Conda environment: R4.3.2 "
module load Miniconda3/4.9.2
source activate R4.3.2

cd /RAID1/working/R425/lavakau
echo " Current working directory: /RAID1/working/R425/lavakau "



### timer
start_time=$(date +%s)
echo ">>> START RUNNING MOTIF SIMILARITY <<<"

### motif files, upper folder path
get_folder_path=/RAID1/working/R425/lavakau/kmer
i_path=${get_folder_path}/sly_acr_pro
echo "RUNNING DIRECTORY: "${i_path}
echo "RUNNING RANGE: proximal"
Rscript ${i_path}/kmer_pool_pro.R ${i_path} "proximal"

end_time=$(date +%s)
timer1=$((end_time - start_time))
echo "RUNNING TIME: ${timer1}"
echo ">>> FINISH! <<<"
#########################################