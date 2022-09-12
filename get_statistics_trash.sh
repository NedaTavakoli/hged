
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
for id in $(seq 1 22; echo X; echo Y)
do
    MyVariants=${DATA}/chr${id}.vcf.gz
    echo("Total Variants", "SNPs", "Indels") > chr${id}_numVariants.txt

    # Count the total number of variants
    #$bcftools view -H  $MyVariants | wc -l 

    # # count number of SNPs
    $bcftools view -H  $MyVariants | grep VT=SNP | wc -l

    # # count number of INDLEs
    $bcftools view -H  $MyVariants | grep VT=INDEL | wc -l