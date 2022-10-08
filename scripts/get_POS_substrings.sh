#!/bin/bash
if [ $# -eq 0 ];
then
    echo "$0: Missing arguments"
    exit 1
elif [ $# -gt 2 ];
then
    echo "$0: Too many arguments: $@"
    exit 1
else
    echo "Construct edge-labeled varion graph ( args: chr id, length of substrings"
    echo "==========================="
    echo "Number of arguments.: $#"
    echo "List of arguments...: $@"
    echo "Arg #1 chrId (Ex: 22)..................................: $1"
    echo "Arg #2 alpha (leng of substring) ......................: $2"
    echo "==========================="
fi

id=$1
alpha=$2

cd ..  
project_dir=$(pwd)
cd data
DATA=$(pwd)
cd ../software
software_dir=$(pwd)
bcftools=${software_dir}/bcftools-1.9/bcftools
samtools=${software_dir}/samtools-1.12/samtools
cd ${DATA}

graph=${DATA}/chr${id}_graph_alpha_${alpha}

MyVariants=chr${id}_snps_indels.vcf.gz
ref=hs37d5

# variant_positions=($(cut -d ',' -f2 variant_positions_snps_indels_chr${id}.txt))
# samples=($($bcftools query -l chr${id}.vcf.gz)) # array of samples, index from 0
# # For each variant position 
# for v in "${variant_positions[@]}"
# do
#     arr=($v)
#     arr+=($($bcftools view -H -r 22:${v} chr${id}.vcf.gz| awk -F"\t" '{split($0, header, "\t");} \
#         {for (i=10; i<=NF; i++) {if (gsub(/0\|1|1\|0|0\/1|1\/0/, "", $(i))==1) {printf header[i]",";printf i-10"\t"} if (i==NF) {printf "\n"}}}'))  
#     echo ${arr[@]}    
# done >> chr${id}_POS_substrings_len_${alpha}_unsorted2.txt


cut -d ',' -f1 variant_positions_snps_indels_chr${id}.txt > variant_positions.txt
$bcftools query -l chr${id}.vcf.gz > samples.txt

variant_positions_file=variant_positions.txt
samples_file=samples.txt
samtools_loc=${samtools}
bcf_loc=${bcftools}
pos_sub_unsorted_file=chr${id}_POS_substrings_len_${alpha}_unsorted2.txt

cd ../src
python extract_substrings.py ${ref} ${id} ${alpha} ${MyVariants} ${variant_positions_file}\
     ${samples_file} ${samtools_loc} ${bcf_loc}

# python extract_substrings.py ${ref} ${id} ${alpha} ${MyVariants} ${variant_positions_file}\
#      ${samples_file} ${pos_sub_unsorted_file} ${samtools_loc} ${bcf_loc}