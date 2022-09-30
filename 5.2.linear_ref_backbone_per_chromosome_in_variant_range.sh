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
|
cd ${project_dir}
bcftools=${software_dir}/bcftools-1.9/bcftools
vcftools=${software_dir}/vcftools-0.1.16/bin/vcftools
tabix=${software_dir}/htslib-1.12/tabix
samtools=${software_dir}/samtools-1.12/samtools
cd ${DATA}
# ***************************************************************************************

# id=22

# Get the linear backbone for the genome graph for all choromosomes 
# ***************************************************************************************
for id in $(seq 1 22; echo X; echo Y)
do
    start=$(cat variant_positions_snps_indels_chr${id}.txt | head -1)
    end=$(cat variant_positions_snps_indels_chr${id}.txt | tail -1)
    $samtools faidx linear_bc_chr${id}.fa ${id}:${start}-${end} > linear_bc_chr${id}_in_variant_range.fa 
done

# Sample output
# cat linear_bc_chr22.fa | head -3
# >22:16050075-51244237
# AGTGGGCCTAAGTGCCTCCTCTCGGGACTGGTATGGGGACGGTCATGCAATCTGGACAAC
# ATTCACCTTTAAAAGTTTATTGATCTTTTGTGACATGCACGTGGGTTCCCAGTAGCAAGA
# 

# ***************************************************************************************

# Important note: to index the linear backbone, we can use the following command
# # indexing the linear backbonefile
#  $samtools faidx linear_bc_chr22.fa  22:16050075-51244237:1-20
#     AGTGGGCCTAAGTGCCTCCT


