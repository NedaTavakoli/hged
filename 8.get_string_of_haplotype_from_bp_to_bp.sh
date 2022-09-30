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

#***************************************************************************************
# Generate vcf files for each samples and chromosomes and index them
#***************************************************************************************

# for id in $(seq 1 22; echo X; echo Y)
for id in {1,22}
do
    while read sample; do
    echo "${sample}"

    # 1: Generate vcf file for a given sample (haplotype), -oz is used to get gz format
    $bcftools view -s ${sample} -Oz chr${id}.vcf.gz > chr_${id}_${sample}.vcf.gz
    echo "done vcf for $sample"

    # 2. Index vcf file (required for bcftools consensus)
    $bcftools index chr_${id}_${sample}.vcf.gz
    done < chr22_sample.txt 
done

#***************************************************************************************
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

id=22
sample=HG00096
ref=hs37d5 
#***************************************************************************************
# Get the substring of length alpha for the first haplotype of the sample
#***************************************************************************************
# Substring of lengh alpha for the first haplotype of the sample  

# $samtools faidx ${ref}.fa 21:9412076-9412080 | $bcftools consensus -s ${sample} -H 1  chr_${id}_${sample}.vcf.gz >  chr${id}_${sample}.fa
$samtools faidx ${ref}.fa 22:16050075-16050175| $bcftools consensus -s ${sample} -H 1  chr_${id}_${sample}.vcf.gz >  chr${id}_${sample}.fa

# Substring of lengh alpha for the second haplotype of the sample   
# $samtools faidx hs37d5.fa 21:9412076-9412080 | $bcftools consensus -s ${sample} -H 2 chr_${id}_${sample}.vcf.gz >  chr${id}_${sample}.fa 
$samtools faidx ${ref}.fa 22:16050075-16050175| $bcftools consensus -s ${sample} -H 1  chr_${id}_${sample}.vcf.gz >  chr${id}_${sample}.fa


#-----------------------------------------------------

id=22
# Get the variant positions from the file 
while read column
do 
echo "$column"
done< file.txt

# Get the variabt positions from the text file and store them into an array 
#   Input: variant_positions_snps_indels_chr${id}.txt
#   output: array variant_positions
variant_positions=($(cut -d ',' -f2 variant_positions_snps_indels_chr${id}.txt))
# printf "%s\n" "${variant_positions[0]}"

#**************************************************



v=16577044
alpha=150
id=22
sample=HG00096
ref=hs37d5 


for alpha in {150, 1000, 5000, 50000}
do
    for id in {1,22}    # for each chromosome
    do
        # For each variant position 
        for v in "${variant_positions[@]}"
        do
            echo ${v}
            while read sample; 
            do
                echo "${sample}"

                # Check if we can ignore haplotypes, if GT=0, we can ignore haplotypes for that position
           
                # Get GT field for the sample     
                GT=$($bcftools view -H chr_${id}_${sample}.vcf.gz | grep ${v}| cut -f 10)
                echo ${GT}

                # only get the information for the first haplotype
                h1=${GT:0:1}

                # only get the information for the second haplotype
                h2=${GT:2:3}

                # Possible GT values: 0|0, 0|1, 0|2, 0|3, 1|0, 1|1, 1|2, 1|3, 2|0, 2|1, 2|2, 2|3, 3|0, 3|1, 3|2, 3|3

                # if GT field for the haplotype is not zero
              
                if [ ${h1} -ne 0 ]; then
                    # Substring of lengh alpha for the first haplotype of the sample 
                    # $samtools faidx ${ref}.fa 22:16577044-16577050 | $bcftools consensus -s ${sample} -H 1  chr_${id}_${sample}.vcf.gz
                    $samtools faidx ${ref}.fa ${id}:${v}-$((${v} + ${alpha}-1)) | $bcftools consensus -s ${sample} -H 1  chr_${id}_${sample}.vcf.gz
                fi

                if [ ${h2} -ne 0 ]; then
                    # Substring of lengh alpha for the second haplotype of the sample   
                    # $samtools faidx hs37d5.fa 21:9412076-9412080 | $bcftools consensus -s ${sample} -H 2 chr${id}_${sample}.vcf.gz >  chr${id}_${sample}.fa 
                    $samtools faidx ${ref}.fa ${id}:${v}-$((${v} + ${alpha}-1)) | $bcftools consensus -s ${sample} -H 2 chr_${id}_${sample}.vcf.gz
                fi
                echo "done consenses for $sample"
            done < chr22_sample.txt 
        done
    done
done








