#!/bin/bash
#****************************************
# To test the repository, use phoenix:

# Login to the cluster:

# ssh ntavakoli6@login-phoenix.pace.gatech.edu 
# ssh ntavakoli6@login-phoenix-3.pace.gatech.edu

# pace-quota
# To submit an interactive job:

#  qsub -I -q inferno -A GT-saluru8-CODA20 -l nodes=1:ppn=1,mem=300gb,walltime=96:00:00
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
samtools=${software_dir}/samtools-1.12/samtools
bgzip=${software_dir}/htslib-1.12/bin/bgzip
cd ${DATA}
id=22
# ***************************************************************************************



# ******* Input files:
# *
# * linear_bc_chr22.fa
# * chr22_snps_indel_POS_REF_ALT.txt

# 


# get the first varinat position
start=$(cat variant_positions_snps_indels_chr${id}.txt | head -1)
echo $start
#16050075
# Chr1: 10177

# get the last variant position
end=$(cat variant_positions_snps_indels_chr${id}.txt | tail -1)
echo $end
#51244237
# for chr1 249240543

# get the length of linear backbone
mystring=$(cat linear_bc_chr${id}.fa)
leng=${#mystring} # note that  ${#mystring} returns the string length
echo ${leng}
#*****************************************************************
#***** generate vertices of edges in the graph
#*****************************************************************

# The coornidate numbers has the 1. linear bc and 2.the length of ALT 
#       1.number of graph vertices is equal to the number of chracters plus 1
num_vertice_linear_bc=$((${leng}+1))
echo ${num_vertice_linear_bc}
#   for the chromosome 22 it is 35780755
#   for the chromosome 1 it is 253384226
#       2.the length of ALT colomn -1, we need to ignore ALT=1, because they won't add any new vertices
num_alt_vertices=$(cat chr${id}_snps_indel_POS_REF_ALT.txt | awk ' length($3)>1'|  awk '{ sum += (length($3)-1); } END { print sum; }')
echo ${num_alt_vertices}
# for chr22 it is 71228
# for chr1 it is 347952

# Total number of vertices for the edge-labled variation graph
Total_vertices=$(($num_vertice_linear_bc + $num_alt_vertices))
echo ${Total_vertices}
# 35893981 for chr22
e=$(($start+ $Total_vertices))
seq $start $e > chr${id}_vertices.txt

#*****************************************************************
#***** generate list of edges in the graph
#*****************************************************************

# #                   1. edges accociated with the linear backbone

# # indexing the linear backbonefile, to get the chracter string
# #for i in $(seq $start $end) # edges for the linear backbone
# end1=$(($end-$start+1))  # we need to index the linear bc to the length of the end1
# echo $end1
# index=${start}
# #for i in $(seq 1 $end1) # edges for the linear backbone


# #for i in $(seq 1 4) # edges for the linear backbone
# for ((i=1;i<=${end1};i++)); 
# do
#     # get the character for that specific range
#     label=$($samtools faidx linear_bc_chr${id}.fa  $id:$start-$end:$i-$(($i)) | tail -1)
#     echo "($index,$(($index+1)),$label)" #>> chr${id}_edges1.txt
#     ((index++))  # increament the index by 1
# done >> chr${id}_edges1.txt


# #               *****************************************************************


# #          2. construct edges associated to paths


# # Get the variant positions from the text file and store them into an array 
# #   Input: variant_positions_snps_indels_chr${id}.txt
# #   output: array variant_positions
# variant_positions=($(cut -d ',' -f2 variant_positions_snps_indels_chr${id}.txt))

# # new edges are started from the last linear bc
# new_v=$(($num_vertice_linear_bc+1))
# echo $new_v
# # 35780756 for chr22

# for POS in "${variant_positions[@]}"
# do
#     REF=$(cat chr${id}_snps_indel_POS_REF_ALT.txt | awk -F " " -v m="$POS" '$1 == m{print $2}')
#     ALT=$(cat chr${id}_snps_indel_POS_REF_ALT.txt | awk -F " " -v m="$POS" '$1 == m{print $3}')
  
#     ref_len=${#REF}
#     # echo $ref_len
#     end_POS_bc=$(($POS+$ref_len))
#     # echo $end_POS_bc

#     # Parse ALT to characters 
#     IFS=', ' read -r -a array <<< "$ALT" 
#     for element in "${array[@]}"
#     do
#         # echo "$element" 
#         # echo "${#element}"

#         # "when ALT has only one character"
#         if [ ${#element} -eq 1 ];
#         then
#         # echo "when ALT has only one character"
#         echo "($POS,$end_POS_bc,$element)" >> chr${id}_edges2.txt
#         fi

#         # when the length of ALT is >1
#         if [[ ${#element} -ne 1 ]];
#         then
#        # echo "when the length of ALT is >1"

#         for ((i=0; i<$((${#element}-1)); i++));
#         do
        
          
#             arr[$i]="${element:$i:1}"
#             label2=${arr[$i]}

#             echo "(${new_v},$(($new_v+1)),$label2)" >> chr${id}_edges2.txt
#             ((new_v=new_v+1))  # increament the vertex by 1
           
#         done
       
#         # "The last character of non-one element"
#         label2=${arr[$((${#element}-1))]}
#         # echo "The last character of non-one element"
#         echo "($new_v,${end_POS_bc},$label2)" >> chr${id}_edges2.txt
#         # echo "$edge2"
#         fi
#     done       
# done 

# Combine two files together and keep them sorted


# cat chr${id}_edges.txt chr${id}_edges2.txt | sort > chr${id}_vg_edges.txt





