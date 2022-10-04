#!/bin/bash

# Change this line according to your projcet directory
cd /storage/coda1/p-saluru8/0/ntavakoli6/hged
id=22
alpha=3
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
# v=51214426
# ***************************************************************************************
ref=hs37d5
variant_positions=($(cut -d ',' -f2 variant_positions_snps_indels_chr${id}.txt))
samples=($($bcftools query -l chr${id}.vcf)) # array of samples, index from 0

# For each variant position 
for v in "${variant_positions[@]}"
do
    a=($v)
    arr=($($bcftools view -H -r 22:${v} chr22.vcf.gz | awk -F"\t" '{split($0, header, "\t");} \
        {for (i=10; i<=NF; i++) {if (gsub(/0\|1|1\|0|0\/1|1\/0/, "", $(i))==1) {printf header[i]",";printf i-10"\t"} if (i==NF) {printf "\n"}}}'))  
    for gt in "${arr[@]}"
    do
        h1=${gt:0:1}
        h2=${gt:2:2}  # it has comma at the end 
        h2=($(echo $h2| sed 's/\(.*\),/\1 /'))  # remove the comma at the end
        sample_index=${gt:4:5}
        sample=${samples[sample_index]}
        if [ ${h1} -ne 0 ]; then
            substring1=($($samtools faidx ${ref}.fa ${id}:${v}-$((${v} + ${alpha}-1)) | $bcftools consensus -s ${sample} -H 1  chr${id}.vcf.gz | tail -1))
        fi
        a+=($substring1)
        if [ ${h2} -ne 0 ]; then
            substring2=($($samtools faidx ${ref}.fa ${id}:${v}-$((${v} + ${alpha}-1)) | $bcftools consensus -s ${sample} -H 2  chr${id}.vcf.gz | tail -1))
        fi
        a+=($substring2)
    done
    echo ${a[@]}
    uniq_s=$(printf "%s\n" "${a[@]}" | sort -u)
    echo $uniq_s >> chr${id}_pos_substrings.txt
done    








