# Sanity check

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
bgzip=${software_dir}/htslib-1.12/bgzip
samtools=${software_dir}/samtools-1.12/samtools

cd ${DATA} 

id=22
sample=HG00096
ref=hs37d5
MyVariants=${DATA}/chr${id}.vcf.gz

cat chr22_sample.txt | head -10
# HG00096
# HG00097
# HG00099
# HG00100
# HG00101
# HG00102
# HG00103
# HG00105
# HG00106
# HG00107

# Get variant positions for chr22
cat chr22_indels.frq.count | grep 51214426
# 22	51214426	2	5008	CGT:4994	C:14

$bcftools view -H  $MyVariants  | grep 51214426


# get the substring of length alpha from the second haplotype of the sample HG00096
$samtools faidx hs37d5.fa 22:51214426-51214430 | $bcftools consensus -s HG00096 -H 2 chr${id}_sample_22.vcf.gz.gz >  chr${id}_sample_${sample}.fa

