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

#cd ..  # TODO: UNCOMENT IN ON REPO
project_dir=$(pwd)
cd data
DATA=$(pwd)
cd ../build
software_dir=$(pwd)
bcftools=${software_dir}/bcftools-1.9/bcftools
samtools=${software_dir}/samtools-1.12/samtools
cd ${DATA}

graph=${DATA}/chr${id}_graph_alpha_${alpha}

id=22
alpha=150

MyVariants=chr${id}_snps_indels.vcf.gz
ref=hs37d5
variant_positions=($(cut -d ',' -f2 variant_positions_snps_indels_chr${id}.txt))
samples=($($bcftools query -l chr${id}.vcf.gz)) # array of samples, index from 0

# For each variant position 
for v in "${variant_positions[@]}"
do
    a=($v)
    arr=($($bcftools view -H -r 22:${v} chr${id}.vcf.gz| awk -F"\t" '{split($0, header, "\t");} \
        {for (i=10; i<=NF; i++) {if (gsub(/0\|1|1\|0|0\/1|1\/0/, "", $(i))==1) {printf header[i]",";printf i-10"\t"} if (i==NF) {printf "\n"}}}'))  
    for gt in "${arr[@]}"
    do
        h1=${gt:0:1}
        h2=${gt:2:2}  # it has comma at the end 
        h2=($(echo $h2| sed 's/\(.*\),/\1 /'))  # remove the comma at the end
        sample_index=${gt:4:5}
        sample=${samples[sample_index]}
        if [ ${h1} -ne 0 ]; then
            substring1=($($samtools faidx ${ref}.fa ${id}:${v}-$((${v} + ${alpha})) | $bcftools consensus -s ${sample} -H 1 ${MyVariants}))
            s1=''
            for i in "${substring1[@]:1}"
            do
                s1+=$i
            done
            t=$((${alpha}))
            len_s1=${#s1}
            min_v1=$((t<len_s1? t : len_s1))
            a+=(${s1:0:min_v1})
        fi
        if [ ${h2} -ne 0 ]; then
            substring2=($($samtools faidx ${ref}.fa ${id}:${v}-$((${v} + ${alpha})) | $bcftools consensus -s ${sample} -H 2 ${MyVariants}))
            s2=''
            for i in "${substring2[@]:1}"
            do
                s2+=$i
            done
            t=$((${alpha}))
            len_s2=${#s2}
            min_v2=$((t<len_s2? t : len_s2))
            a+=(${s2:0:min_v2})
        fi
    done
    # echo ${a[@]}
    uniq_s=$(printf "%s\n" "${a[@]}" | sort -u)
    echo $uniq_s
done >> ${graph}/chr${id}_POS_substrings_len_${alpha}_2.txt






