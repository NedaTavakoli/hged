#!/bin/bash
id=22
alpha=50

MyVariants=chr${id}_snps_indels.vcf.gz
ref=hs37d5
variant_positions=($(cut -f2 chr${id}_snps_indel_POS_REF_ALT.txt))
samples=($(bcftools query -l chr${id}_snps_indels.vcf.gz)) # array of samples, index from 0

# For each variant position 
#for v in "${variant_positions[@]}"
for v in "${variant_positions:0:499}" # specify a range of variant position just for testing
do
    a=($v)
    arr=($(bcftools view -H -r 22:${v} chr${id}_snps_indels.vcf.gz |  awk -F"\t" '{split($0, header, "\t");} \
        {for (i=10; i<=NF; i++) {if ((gsub(/0\|0|0\/0|/, "", $(i)) !=1))  {printf header[i]",";printf i-10"\t"} if (i==NF) {printf "\n"}}}'))   

    for gt in "${arr[@]}"
    do
        h1=${gt:0:1} 
        h2=${gt:2:2}  # it has comma at the end 
        h2=($(echo $h2| sed 's/\(.*\),/\1 /'))  # remove the comma at the end
        sample_index=${gt:4:5}
        sample=${samples[sample_index]}
        if [ ${h1} -ne 0 ]; then
            s1=$(samtools faidx ${ref}.fa ${id}:${v}-$((${v} + ${alpha}-1)) | bcftools consensus -s ${sample} -H 1 ${MyVariants} |  sed '1d' | tr -d "[:space:]")
            len_s1=$(echo ${s1[@]} | wc -c)
            t=$((${alpha}))
            min_v1=$((t<len_s1? t : len_s1))
            a+=(${s1:0:min_v1})
            rm -f tmp_*
        fi
        if [ ${h2} -ne 0 ]; then
            s2=$(samtools faidx ${ref}.fa ${id}:${v}-$((${v} + ${alpha}-1)) | bcftools consensus -s ${sample} -H 2 ${MyVariants} |  sed '1d' | tr -d "[:space:]")
            len_s2=$(echo ${s2[@]} | wc -c)
            t=$((${alpha}))
            min_v2=$((t<len_s2? t : len_s2))
            a+=(${s2:0:min_v2})
        fi
    done
    # echo ${a[@]} 
    uniq_s=$(printf "%s\n" "${a[@]}" | sort -u)
    echo $uniq_s
done >> chr${id}_POS_substrings_len_${alpha}.txt








