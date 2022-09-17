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

id=22
sample=HG00096
ref=hs37d5  

# ******* Input files:
# *
# * linear_bc_chr22.fa
# * chr22_snps_indel_POS_REF_ALT.txt

# Get the first variant position 
adjust=$(cat chr22_snps_indel_POS_REF_ALT.txt | head -1 | cut -f1)
${adjust}
# 16050075

# Adjust postions in POS REF ALT to start from 0
cat chr22_POS_REF_ALT.txt | awk -v s=${adjust} '{print $1-s"\t"$2"\t"$3}' > chr22_adjustedPOS_REF_ALT.txt

# note that  ${#mystring} returns the string length
# https://stackoverflow.com/questions/71502887/how-to-retrieve-a-character-in-a-string-at-a-index-in-bash
mystring=$(cat linear_bc_chr22.fa)
${#mystring:1:5}

echo ${mystring}  | awk '{split($0,a,""); print a[1]}'
# print a[1] This is for Index Character a[2] ... a[n] you can do

# {split($0,a,"") This is for character to split "" meaning split every character

# index the fasta file of choromosome 22 
$samtools faidx linear_bc_chr22.fa
#22:16050075-51244237	35194163	22	60	61

cat linear_bc_chr22.fa | head -10
# >22:16050075-51244237
# AGTGGGCCTAAGTGCCTCCTCTCGGGACTGGTATGGGGACGGTCATGCAATCTGGACAAC
# ATTCACCTTTAAAAGTTTATTGATCTTTTGTGACATGCACGTGGGTTCCCAGTAGCAAGA
# AACTAAAGGGTCGCAGGCCGGTTTCTGCTAATTTCTTTAATTCCAAGACAGTCTCAAATA
# TTTTCTTATTAACTTCCTGGAGGGAGGCTTATCATTCTCTCTTTTGGATGATTCTAAGTA
# CCAGCTAAAATACAGCTATCATTCATTTTCCTTGATTTGGGAGCCTAATTTCTTTAATTT
# AGTATGCAAGAAAACCAATTTGGAAATATCAACTGTTTTGGAAACCTTAGACCTAGGTCA
# TCCTTAGTAAGATCTTCCCATTTATATAAATACTTGCAAGTAGTAGTGCCATAATTACCA
# AACATAAAGCCAACTGAGATGCCCAAAGGGGGCCACTCTCCTTGCTTTTCCTCCTTTTTA
# GAGGATTTATTTCCCATTTTTCTTAAAAAGGAAGAACAAACTGTGCCCTAGGGTTTACTG

# indexing the linear backbonefile
$samtools faidx linear_bc_chr22.fa  22:16050075-51244237:1-20
#AGTGGGCCTAAGTGCCTCCT

$samtools faidx hs37d5.fa 22:16050075-51244237