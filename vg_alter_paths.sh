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
samtools=${software_dir}/samtools-1.12/samtools
cd ${DATA}
# ***************************************************************************************

id=22
sample=HG00096
ref= hs37d5  

# Input: chr22_snps_indel_POS_REF_ALT.txt
    # cat chr22_snps_indel_POS_REF_ALT.txt | head -10
    # 16050075	A	G
    # 16050115	G	A
    # 16050213	C	T
    # 16050319	C	T
    # 16050527	C	A
    # 16050568	C	A
    # 16050607	G	A
    # 16050627	G	T
    # 16050646	G	T
    # 16050655	G	A


    # cat chr22_snps_indel_POS_REF_ALT.txt | tail -10
    # 51241101	A	T
    # 51241102	T	C
    # 51241285	T	G
    # 51241298	T	G
    # 51241309	C	T
    # 51241342	C	A
    # 51241386	C	G
    # 51244163	A	G
    # 51244205	C	T
    # 51244237	C	T


# Get the linear backbone for the genome graph for chr22  
# ***************************************************************************************
cat chr${id}_snps_indel_POS_REF_ALT.txt | head -1
cat chr${id}_snps_indel_POS_REF_ALT.txt | tail -1
$samtools faidx hs37d5.fa 22:16050075-51244237 
