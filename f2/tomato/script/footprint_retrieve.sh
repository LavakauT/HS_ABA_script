#!/bin/bash
#SBATCH -J BinRet
#SBATCH -o BinRet.out
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 20  
#SBATCH --mem-per-cpu 1G
#SBATCH --time 00-12:00

echo "Job submitted"
echo "Requiring tasks under this TOBIAS user"

echo ">>>ENVIRONMENT LAODING<<<"
echo " Module load Miniconda3/4.9.2 "
echo " Conda environment: R4.3.2 "
module load Miniconda3/4.9.2
source activate R4.3.2

path=/RAID1/working/R425/huanchi/atac/tomato
cd ${path}
echo " Current working directory: ${path} "



### timer
start_time=$(date +%s)
echo ">>> START RUNNING MOTIF SIMILARITY <<<"

### motif files, upper folder path
get_folder_path=${path}/ATACorrect_DEP
i_path=${get_folder_path}/BINDetect_pd
echo "RUNNING DIRECTORY: "${i_path}
echo "RUNNING RANGE: proximal"
Rscript ${i_path}/footprint_retrieve.R ${i_path}

end_time=$(date +%s)
timer1=$((end_time - start_time))
echo "RUNNING TIME: ${timer1}"
echo ">>> FINISH! <<<"
#########################################