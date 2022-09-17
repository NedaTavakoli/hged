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
# ***************************************************************************************


# List of samples  #printing a list of samples from a VCF:
# Each sample shows the information of two halplotypes 1 and 2
$bcftools query -l chr${id}.vcf > chr${id}_samples.txt
    #cat chr22_sample.txt | head -10
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

id=22
sample=HG00096
ref= hs37d5    
#***************************************************************************************

# 1: Generate vcf file for a given sample (haplotype), -oz is used to get gz format

$bcftools view -s ${sample} -Oz chr${id}.vcf.gz > chr${id}_${sample}.vcf.gz

#***************************************************************************************

# 2. Index vcf file (required for bcftools consensus)

$bcftools index chr${id}_${sample}.vcf.gz
#***************************************************************************************

# 3. Convert vcf file to fasta, for a halplotype in a region

    # Get the consensus for one region. The fasta header lines are then expected
    # in the form ">chr:from-to".

    # To get the range for variant positions: for example
    # cat variant_positions_snps_indels_chr21.txt | head -10
    # 9411239
    # 9411245
    # 9411264
    # 9411267
    # 9411302
    # 9411313
    # 9411332
    # 9411347
    # 9411356

    # cat variant_positions_snps_indels_chr22.txt | head -10
    # 16050075
    # 16050115
    # 16050213
    # 16050319
    # 16050527
    # 16050568
    # 16050607
    # 16050627
    # 16050646
    # 16050655

# Substring of lengh alpha for the first haplotype of the sample    
$samtools faidx ${ref}.fa 21:9412076-9412080 | $bcftools consensus -s ${sample} -H 1  chr${id}_${sample}.vcf.gz >  chr${id}_${sample}.fa

# Substring of lengh alpha for the second haplotype of the sample   
$samtools faidx hs37d5.fa 21:9412076-9412080 | $bcftools consensus -s ${sample} -H 2 chr${id}_${sample}.vcf.gz >  chr${id}_${sample}.fa 


