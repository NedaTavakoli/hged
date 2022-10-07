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
    echo "Arg #3 subrange .......................................: $3"
    echo "==========================="
fi

id=$1
alpha=$2
subrange=$3

#cd ..  # TODO: UNCOMENT IN ON REPO
project_dir=$(pwd)
cd data
DATA=$(pwd)
cd ../software
software_dir=$(pwd)
bcftools=${software_dir}/bcftools-1.9/bcftools
samtools=${software_dir}/samtools-1.12/samtools
cd ${DATA}
mkdir graph_alpha_${alpha}_subrange_${subrange}
subgraph=${DATA}/graph_alpha_${alpha}_subrange_${subrange}

# ****** Changes for the subrange starts from here ****************
start=$(cat variant_positions_snps_indels_chr${id}.txt | head -1)
end=$(cat variant_positions_snps_indels_chr${id}.txt | head -${subrange} |tail -1)
REF_end=$(cat chr${id}_snps_indel_POS_REF_ALT.txt | head -${subrange} |tail -1 | awk '{print $2}')
len_REF_end=${#REF_end}
t=$((${end}))
min_last=$((len_REF_end<t? len_REF_end: t))
end=$(($end+$min_last))
$samtools faidx hs37d5.fa ${id}:${start}-${end} > ${subgraph}/linear_bc_chr${id}_alpha_${alpha}_subrange_${subrange}.fa
REF_l=$(cat linear_bc_chr${id}_alpha_${alpha}_subrange_${subrange}.fa)
leng=${#REF_l} 
num_vertice_linear_bc=$((${leng}+1))
num_alt_vertices=$(cat chr${id}_snps_indel_POS_REF_ALT.txt | head -${subrange} | awk 'length($3)>1'| awk '{ sum += (length($3)-1); } END { print sum; }')
Total_vertices=$(($num_vertice_linear_bc + $num_alt_vertices))
e=$(($start+ $Total_vertices))
seq $start $e > ${subgraph}/chr${id}_vertices_subrange_${subrange}.txt

cat variant_positions_snps_indels_chr${id}.txt | head -${subrange} > ${subgraph}/variant_positions_snps_indels_chr${id}_head_${subrange}.txt
cat chr${id}_snps_indel_POS_REF_ALT.txt | head -${subrange} > ${subgraph}/chr${id}_snps_indel_POS_REF_ALT_head_${subrange}.txt

variant_pos_file=variant_positions_snps_indels_chr${id}_head_${subrange}.txt
variant_POS_ALT_REF_file=chr${id}_snps_indel_POS_REF_ALT_head_${subrange}.txt
linear_bc_file=linear_bc_chr${id}_alpha_${alpha}_subrange_${subrange}.fa
linear_edges=chr${id}_linear_edges_alpha_${alpha}_subrange_${subrange}.txt
alt_edges=chr${id}_alt_edges2_alpha_${alpha}_subrange_${subrange}.txt

cd ../src
python get_edges_chr.py ${start} ${end} ${num_vertice_linear_bc}\
    ${subgraph}/${variant_pos_file} ${subgraph}/${variant_POS_ALT_REF_file} ${subgraph}/${linear_bc_file}\
     ${id} ${linear_edges} ${alt_edges} 

cat chr${id}_linear_edges_alpha_${alpha}_subrange_${subrange}.txt chr${id}_alt_edges2_alpha_${alpha}_subrange_${subrange}.txt | sort > ${subgraph}/chr${id}_vg_edges_alpha_${alpha}_subrange_${subrange}.txt

rm chr${id}_alt_edges2_alpha_${alpha}_subrange_${subrange}.txt
rm chr${id}_linear_edges_alpha_${alpha}_subrange_${subrange}.txt

