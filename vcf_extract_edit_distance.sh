#!/bin/bash
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
MyVariants=${DATA}/chr${id}.vcf.gz
#***************************************************************************************************************************
#*   In the below code, we use BCFtools to count the number of SNPs and INDELs at each vcf file for each choromosome
#*   For each chromosome extract SNPs only and INDELs only
#*   In the log files # of SNPs and INDELs are kept, as well as the time that vcftools is used to extract these data
#*
#***************************************************************************************************************************
for id in $(seq 1 22; echo X; echo Y)
do
    MyVariants=${DATA}/chr${id}.vcf.gz

    # Count the total number of variants
    #$bcftools view -H  $MyVariants | wc -l 

    # # count number of SNPs
   if grep -q VT=SNP "$bcftools view -H  $MyVariants "; then print $2; fi

    # # count number of INDLEs
    $bcftools view -H  $MyVariants | grep VT=INDEL | wc -l
        
done
#***************************************************************************************************************************
# *
# *
#***************************************************************************************************************************
 # get list of variant positions (INDELS+ SNPs)
# Saving SNP and INDEL positions in two speparrte files
#***************************************************************************************************************************

for id in $(seq 1 22; echo X; echo Y)
do
    MyVariants=${DATA}/chr${id}.vcf.gz
   
    #  saving SNP postions to files
    $bcftools view -H  $MyVariants | grep VT=SNP  | awk '{print $2}'  > snps_positions_chr${id}.txt

    
    #  saving indel postions to files
    $bcftools view -H  $MyVariants | grep VT=INDEL  | awk '{print $2}'  > indels_positions_chr${id}.txt


    sort -n snps_positions_chr${id}.txt  indels_positions_chr${id}.txt > variant_positions_snps_indels_chr${id}.txt

done
#***************************************************************************************************************************
#* Compute statistics on the variant positions ( how far they are from each other)
#***************************************************************************************************************************
# Find average distance between values of a colomn 
for id in $(seq 1 22; echo X; echo Y)
do
    variant_positions_snps_indels_chr${id}.txt





done











