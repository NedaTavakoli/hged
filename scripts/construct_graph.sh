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

project_dir=$(pwd)
cd data
DATA=$(pwd)
cd ../build
# cd ../software
software_dir=$(pwd)
bcftools=${software_dir}/bcftools-1.9/bcftools
samtools=${software_dir}/samtools-1.12/samtools
cd ${DATA}

mkdir -p chr${id}_graph_alpha_${alpha}
graph=${DATA}/chr${id}_graph_alpha_${alpha}
graph=$(echo $graph| tr -d '\r') # remove -r from the directory

cd ${DATA}

start=$(cat variant_positions_snps_indels_chr${id}.txt | head -1)
end=$(cat variant_positions_snps_indels_chr${id}.txt |tail -1)
REF_end=$(cat chr${id}_snps_indel_POS_REF_ALT.txt | tail -1 | awk '{print $2}')
len_REF_end=${#REF_end}
end=$(($end+$len_REF_end))
$samtools faidx hs37d5.fa ${id}:${start}-${end} > ${graph}/linear_bc_chr${id}_alpha_${alpha}.fa
REF_l=$($samtools faidx hs37d5.fa ${id}:${start}-${end} |  sed '1d' | tr -d "[:space:]") # remove the first linea and remove whitespaces
leng=$(echo ${REF_l[@]} | wc -c)
# echo ${REF_l:0:20}
num_vertice_linear_bc=$((${leng}+1))
num_alt_vertices=$(cat chr${id}_snps_indel_POS_REF_ALT.txt | awk ' length($3)>1'| awk '{sum = 0; t = split($3, arr, ","); for(i = 1; i <= t; i++) sum+=(length(arr[i])-1);print sum}' | awk '{ sum_t += $1 } END { print sum_t }')
Total_vertices=$(($num_vertice_linear_bc + $num_alt_vertices))
e=$(($start+ $Total_vertices))
seq $start $e > ${graph}/chr${id}_vertices_alpha_${alpha}.txt
#num_edges_for_alpha=150 : 35252189
#num_edges_for_alpha=50 : 35252188


variant_pos_file=variant_positions_snps_indels_chr${id}.txt
variant_POS_ALT_REF_file=chr${id}_snps_indel_POS_REF_ALT.txt

linear_bc_file=linear_bc_chr${id}_alpha_${alpha}.fa
linear_edges=chr${id}_linear_edges_alpha_${alpha}.txt
alt_edges=chr${id}_alt_edges2_alpha_${alpha}.txt



cd ../src
python get_edges_chr.py ${start} ${end} ${num_vertice_linear_bc}\
    ${DATA}/${variant_pos_file} ${DATA}/${variant_POS_ALT_REF_file} ${graph}/${linear_bc_file}\
     ${id} ${linear_edges} ${alt_edges} 

cat chr${id}_linear_edges_alpha_${alpha}.txt chr${id}_alt_edges2_alpha_${alpha}.txt | sort > ${graph}/chr${id}_vg_edges_alpha_${alpha}.txt

rm chr${id}_linear_edges_alpha_${alpha}.txt 
rm chr${id}_alt_edges2_alpha_${alpha}.txt