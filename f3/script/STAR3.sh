#!/bin/bash
#SBATCH -J STAR_MAP
#SBATCH -o STAR_mapping_mp_3.out
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 20  
#SBATCH --mem-per-cpu 1G
#SBATCH --time 01-00:00


echo "###### Job submitted ######"
date
echo "Requiring $SLURM_NTASKS task with $SLURM_CPUS_PER_TASK cpus in total $SLURM_MEM_PER_CPU Memory under this RNA-seq user"
echo

echo "###### Environment loading..... ######"
echo "module load Miniconda3/4.9.2 STAR/2.7.9a-GCC-11.2.0"
module load Miniconda3/4.9.2
module load STAR/2.7.9a-GCC-11.2.0
echo
echo "source activate base"
source activate base

path=/RAID1/working/R425/lavakau/rna-seq/marchantia/20240530_la
cd ${path}
echo "###### Current working directory: $path ######"
echo

# Function to print a separator line
print_separator() {
    printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' -
}


###### set up array for parallel running ######
# Specify the path to the config file
config=/RAID1/working/R425/lavakau/rna-seq/marchantia/20240530_la/rawData/range.txt
print_separator
###########################################


############# MAPPING LOOP
# loop mapping from file 21~30

indir=/RAID1/working/R425/lavakau/rna-seq/marchantia/20240530_la/fastp/clean_fq
outdir=/RAID1/working/R425/lavakau/rna-seq/marchantia/20240530_la/fastp/bam

mkdir ${outdir}

for i in {21..30}
do
    # Extract the sample name for the current $SLURM_PROCID
    index=${i}
    sample=$(awk -v ArrayTaskID=$index '$1==ArrayTaskID {print $2}' $config)
    echo "STAR mapping task $index is processing $sample"

    ###### MAPPING ######
    STAR --genomeDir genomeIndex \
        --runThreadN 10 \
        --readFilesIn ${indir}/${sample}_R1.clean.fq.gz ${indir}/${sample}_R2.clean.fq.gz \
        --outFileNamePrefix ${outdir}/${sample} \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped None \
        --outFilterMismatchNmax 3 \
        --outFilterMultimapNmax 1 \
        --outSAMattributes All --limitBAMsortRAM 1222620308
    echo "STAR mapping task $index finished processing $sample"
    print_separator
done

# wait for finishing all the samples
wait

echo "###### Finished STAR mapping ######"
###########################################