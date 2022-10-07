# Haplotype_aware Variation Selection in Genome Graphs under Edit Disatance

The overall workflow is:

```sh
git clone https://github.com/NedaTavakoli/hged
cd hged
project_dir=$(pwd)  #project top-level directory
# download data and softwares
chmod +x dependencies.sh
./dependencies.sh
# cosntruct edge_labeled variation graph 
chmod +x script/construct_graph.sh
./script/construct_graph.sh
# construct pos and list of substrings of length alpha from halpotypes
chmod +x script/get_POS_substrings.sh
./script/get_POS_substrings.sh
# OPTIONAL: Construct fasta file per each haplotype
chmod +x script/get_fa_for_each_haplotype.sh
# ./script/get_fa_for_each_haplotype.sh
```

List of vertices

List of edges (start, end, symbol, variant number)

List of variant positions

Number of variants: 
    chr22: 1102765  
    chr1: 6463830

Number of haplotypes: 5008
