#!/bin/bash
#****************************************
# To test the repository, use phoenix:

# Login to the cluster:

# ssh ntavakoli6@login-phoenix.pace.gatech.edu 
# ssh ntavakoli6@login-phoenix-3.pace.gatech.edu

# pace-quota
# To submit an interactive job:

#  qsub -I -q inferno -A GT-saluru8-CODA20 -l nodes=1:ppn=24,mem=100gb,walltime=96:00:00
#*****************************************

# # Change this line according to your projcet directory
# cd /storage/coda1/p-saluru8/0/ntavakoli6/hged
# #****  NOTE:  This part needs to be done only once to download vcf and the required reference
# chmod +x download_vcf_and_ref.sh
# ./download vcf_and_ref.sh
# chmod +x download_sw_dependencies.sh
# ./download_sw_dependencies.sh
# ***************************************************************************************

# Change this line according to your projcet directory
cd /storage/coda1/p-saluru8/0/ntavakoli6/hged
#----------------------------------------------------------------

project_dir=$(pwd)  #project top-level directory

cd ${project_dir}/data
DATA=$(pwd)

cd ${project_dir}/software
software_dir=$(pwd)

cd ${project_dir}
bcftools=${software_dir}/bcftools-1.9/bcftools
vcftools=${software_dir}/vcftools-0.1.16/bin/vcftools
tabix=${software_dir}/htslib-1.12/tabix

cd ${DATA}
# ***************************************************************************************
# extract required data to make the graph

for id in $(seq 1 22; echo X; echo Y)
do
    MyVariants=${DATA}/chr${id}.vcf.gz
   
    #  saving SNP postions to files
    $bcftools view -H  $MyVariants | grep VT=SNP  | awk '{print $2"\t"$4"\t"$5}'  > t1_chr${id}.txt

    
    #  saving indel postions to files
    $bcftools view -H  $MyVariants | grep VT=INDEL  |awk '{print $2"\t"$4"\t"$5}'  > t2_chr${id}.txt


    sort -n t1_chr${id}.txt  t2_chr${id}.txt > chr${id}_snps_indel_POS_REF_ALT.txt

done

