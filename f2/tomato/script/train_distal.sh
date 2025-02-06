#!/bin/bash
#SBATCH -J ML
#SBATCH -o ML_acr_dis_%a.out
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 15  
#SBATCH --mem-per-cpu 1G
#SBATCH --array=1-5
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

###### set up array for parallel running ######
# Specify the path to the config file
config=/RAID1/working/R425/lavakau/pCRE/sly_acr/range.txt

# Extract the sample name for the current $SLURM_ARRAY_TASK_ID
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
###########################################


path=/RAID1/working/R425/lavakau/ML_master
cd ${path}
indir=${path}/sly_acr_dis

###### CLEAN DATA ######
# WARNING: please change the 1,0 to pos,neg after pCRE
# cd ML_master
# python ML_preprocess.py -df [cluster_aba/hsfb_dwn_10.txt] -na_method median -onehot t

for inp in ${indir}/${sample}/*df_p0.01.txt
do
    python ML_preprocess.py -df ${inp} -na_method median -onehot t
done
###########################################


###### DEFINE TEST SET ######
for inp in ${indir}/${sample}/*df_p0.01_mod.txt
do
    python test_set.py -df ${inp} -use pos,neg -type c -p 0.1 -save ${inp}.test
done
###########################################


###### SELECT FEATURES BY LASSO ######
for inp in ${indir}/${sample}/*df_p0.01_mod.txt
do
    python Feature_Selection.py -df ${inp} -cl_train pos,neg -type c -alg lasso -p 0.01 -save ${inp}.top_feat_lasso
done
###########################################




echo "Preprocess finished! Pass files to training procedure."




###### TRAIN WITH LASSO ######
# IF IT IS UNFEASIBLE TO DO LASSO, IT WILL TRAIN WITH NO LASSO FILE.
echo ">>> SVM training <<<"
for inp in ${indir}/${sample}/*df_p0.01_mod.txt
do 
    if [ -s "${inp}.top_feat_lasso" ]; then
        echo "Processing "${inp}" in LogReg with lasso features"
        python ML_classification_modified.py \
            -df ${inp} -test ${inp}.test -feat ${inp}.top_feat_lasso -cl_train pos,neg -alg SVM -tag lasso -apply all -n 10 -cv_num 10 -gs_n 10 -n_jobs 30 -plots T
    else
        echo "Processing "${inp}" in LogReg without lasso features"
        python ML_classification_modified.py \
            -df ${inp}  -test ${inp}.test -cl_train pos,neg -alg SVM -tag lasso -apply all -n 10 -cv_num 10 -gs_n 10 -n_jobs 30 -plots T
    fi
done


echo ">>> LogReg training <<<"
for inp in ${indir}/${sample}/*df_p0.01_mod.txt
do 
    if [ -s "${inp}.top_feat_lasso" ]; then
        echo "Processing "${inp}" in LogReg with lasso features"
        python ML_classification_modified.py \
            -df ${inp} -test ${inp}.test -feat ${inp}.top_feat_lasso -cl_train pos,neg -alg LogReg -tag lasso -apply all -n 10 -cv_num 10 -gs_n 10 -n_jobs 30 -plots T
    else
        echo "Processing "${inp}" in LogReg without lasso features"
        python ML_classification_modified.py \
            -df ${inp}  -test ${inp}.test -cl_train pos,neg -alg LogReg -tag lasso -apply all -n 10 -cv_num 10 -gs_n 10 -n_jobs 30 -plots T
    fi
done


echo ">>> RF training <<<"
for inp in ${indir}/${sample}/*df_p0.01_mod.txt
do 
    if [ -s "${inp}.top_feat_lasso" ]; then
        echo "Processing "${inp}" in LogReg with lasso features"
        python ML_classification_modified.py \
            -df ${inp} -test ${inp}.test -feat ${inp}.top_feat_lasso -cl_train pos,neg -alg RF -tag lasso -apply all -n 10 -cv_num 10 -gs_n 10 -n_jobs 30 -plots T
    else
        echo "Processing "${inp}" in LogReg without lasso features"
        python ML_classification_modified.py \
            -df ${inp}  -test ${inp}.test -cl_train pos,neg -alg RF -tag lasso -apply all -n 10 -cv_num 10 -gs_n 10 -n_jobs 30 -plots T
    fi
done
###########################################





###### TRAIN WITHOUT LASSO ######
# IF IT IS UNFEASIBLE TO DO LASSO, IT WILL TRAIN WITH NO LASSO FILE.
echo ">>> SVM training <<<"
for inp in ${indir}/${sample}/*df_p0.01_mod.txt
do 
    echo "Processing "${inp}" in LogReg without lasso features"
    python ML_classification_modified.py \
        -df ${inp}  -test ${inp}.test -cl_train pos,neg -alg SVM -tag no_lasso -apply all -n 10 -cv_num 10 -gs_n 10 -n_jobs 30 -plots T
done


echo ">>> LogReg training <<<"
for inp in ${indir}/${sample}/*df_p0.01_mod.txt
do 
    echo "Processing "${inp}" in LogReg without lasso features"
    python ML_classification_modified.py \
        -df ${inp}  -test ${inp}.test -cl_train pos,neg -alg LogReg -tag no_lasso -apply all -n 10 -cv_num 10 -gs_n 10 -n_jobs 30 -plots T
done


echo ">>> RF training <<<"
for inp in ${indir}/${sample}/*df_p0.01_mod.txt
do 
    echo "Processing "${inp}" in LogReg without lasso features"
    python ML_classification_modified.py \
        -df ${inp}  -test ${inp}.test -cl_train pos,neg -alg RF -tag no_lasso -apply all -n 10 -cv_num 10 -gs_n 10 -n_jobs 30 -plots T
done
###########################################

echo "Job finished!"