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
    echo "We got some argument(s)"
    echo "==========================="
    echo "Number of arguments.: $#"
    echo "List of arguments...: $@"
    echo "Arg #1 chrId (Ex: 22)..............:  $1"
    echo "Arg #2 BCFtools_path  ..............: $2"
    echo "Arg #3 SAMtools_path  ..............: $3"
    echo "Arg #4 VCF_file  ...................: $4"
    echo "Arg #5 ref (without fa) ......................: $5"
    echo "==========================="
fi

# ***************************************************************************************
# inputs
samtools=$2
bcftools=$3
MyVariants=$3
ref=$5 
# ***************************************************************************************


# get list of variant positions_indel_snps
# ***************************************************************************************
$bcftools view -H  $MyVariants | grep VT=SNP  | awk '{print $2}'  > snps_positions_chr${id}.txt
$bcftools view -H  $MyVariants | grep VT=INDEL  | awk '{print $2}'  > indels_positions_chr${id}.txt
sort -n snps_positions_chr${id}.txt  indels_positions_chr${id}.txt > variant_positions_snps_indels_chr${id}.txt
# ***************************************************************************************


# get vcf" POSITON REF ALT
# ***************************************************************************************
$bcftools view -H  $MyVariants | grep VT=SNP  | awk '{print $2"\t"$4"\t"$5}'  > t1_chr${id}.txt   
$bcftools view -H  $MyVariants | grep VT=INDEL  |awk '{print $2"\t"$4"\t"$5}'  > t2_chr${id}.txt
sort -n t1_chr${id}.txt  t2_chr${id}.txt > chr${id}_snps_indel_POS_REF_ALT.txt
# ***************************************************************************************

# get number of variant position
# ***************************************************************************************
Num_variant=$(cat chr${id}_snps_indel_POS_REF_ALT.txt | wc -l)
print("Number of variant position is:",Num_variant)
# ***************************************************************************************

# get linar bc
# ***************************************************************************************
start=$(cat variant_positions_snps_indels_chr${id}.txt | head -1)
end=$(cat variant_positions_snps_indels_chr${id}.txt | tail -1)
$samtools faidx hs37d5.fa ${id}:${start}-${end} > linear_bc_chr${id}.fa
# ***************************************************************************************

# get linar bc within the range
# ***************************************************************************************
$samtools faidx linear_bc_chr${id}.fa ${id}:${start}-${end} > linear_bc_chr${id}_in_variant_range.fa 
# ***************************************************************************************


#    Get the list of vertices
# ***************************************************************************************
start=$(cat variant_positions_snps_indels_chr${id}.txt | head -1)
end=$(cat variant_positions_snps_indels_chr${id}.txt | tail -1)
mystring=$(cat linear_bc_chr${id}.fa)
leng=${#mystring} 
num_vertice_linear_bc=$((${leng}+1))
num_alt_vertices=$(cat chr${id}_snps_indel_POS_REF_ALT.txt | awk ' length($3)>1'|  awk '{ sum += (length($3)-1); } END { print sum; }')
Total_vertices=$(($num_vertice_linear_bc + $num_alt_vertices))
e=$(($start+ $Total_vertices))
seq $start $e > chr${id}_vertices.txt
# ***************************************************************************************

# get list of samples
# ***************************************************************************************
$bcftools query -l chr${id}.vcf > chr${id}_samples.txt
# ***************************************************************************************


# get number of haplotypes
# ***************************************************************************************
num_samples=$(cat chr${id}_samples.txt | wc -l)
num_haplotypes=$(($num_samples*2))
echo 'Number of haplotypes:'$num_haplotypes
# ***************************************************************************************

# get list of edges
# **************************************************************************************
# after running the following python file
python 7.1.vg.edges.py --chr ${id} --start ${start}\
    --end ${end} --bc ${num_vertice_linear_bc}
cat chr${id}_edges.txt chr${id}_edges_2.txt | sort > chr${id}_vg_edges.txt


