#!/bin/bash
#SBATCH -J ChromHMM
#SBATCH -o ChromHMM.out
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 20  
#SBATCH --mem-per-cpu 1G
#SBATCH --time 01-00:00

module load Miniconda3/4.9.2
module load Java/11.0.2
source activate gtf2gpd

cd /RAID1/working/R425/lavakau/ChromHMM
# gtfToGenePred -genePredExt -ignoreGroupsWithoutExons -geneNameAsName2 marpol/MpTak_v6.1r1.gtf marpol/MpTak_v6.1r1.gpd
# java -mx4000M -jar ChromHMM.jar ConvertGeneTable -l marpol/chromosomelengthfile.txt  marpol/MpTak_v6.1r1.gpd marpol marpol
# java -mx4000M -jar ChromHMM.jar BinarizeBam marpol/chromosomelengthfile.txt histone_ocr/ marpol/cellmarktable.txt marpol_chrom/ 


# java -mx4000M -jar ChromHMM.jar LearnModel marpol_chrom/ marpol_chrom_state9/ 9 marpol
java -mx4000M -jar ChromHMM.jar LearnModel marpol_chrom/ marpol_chrom_state12/ 12 marpol
# java -mx4000M -jar ChromHMM.jar LearnModel marpol_chrom/ marpol_chrom_state15/ 15 marpol



# chose 12 states
java -mx4000M -jar ChromHMM.jar OverlapEnrichment marpol_chrom_state12/Tak1_12_segments.bed COORDS/marpol marpol_12
echo "Finished ChromHMM!"