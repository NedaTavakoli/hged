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


## This repository is to solve variant selection in genome graphs under edit disatnce
in other words for a given complete variation graph, it creates a reduced variation graph in which 
some variants are removed subject to some constraints. The constraints are for every substring of length 
alpha observed in haplotypes, the reduced varaition graph guarantees to preserve those substrings with
at most delta errors (i.e., edit distance of delta among alpha-long substrings of haplotypes in complete variation graph with those of reduced variation graph).

The project has the following folder structure:
```
hged
|___script # scripts to extract data from genome graphs
|   |___construct_graph.sh # construct complete edge-labeled variation graph
|   |___get_POS_substrings.sh # construct haplotypes substrings of length alpha for each variant position
|
|___src  
    |___get_edges_chr.py # to construct edges of variation graph, used in src/construct_graph.sh
    |___main.py # code to construct ILPs 
...
```




