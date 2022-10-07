<p align="center">
 <span style="font-size: 26px">Haplotype_aware Variation Selection in Genome Graphs under Edit Disatance</span>
</p>

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

This repository is used to solve variant selection in genome graphs under edit disatnce
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

The algorithm has the following steps:
```
1- Loading data: 
    Inputs: edges_file_name, location_substring_file_name
    outputs: E, locations, substrings, num_variants

2- Cosntruct graph  
    Inputs: E
    Outputs: G

3- Create global ILP
    Inputs: G, locations, substrings, number_variants, alpha, delta
    Outputs: model  

        |___Create global ILP 
        |   |___G_ind = reachable_subgraph(G, pos, alpha + delta) #bFind the reachable subgraph of distance d from a pos
        |   |___G_a, start_v, end_v = create_alignment_graph(G_ind, pos, S) # igonre source, pos is the top-left vertex
        |   |___G_a_pruned = prune_alignment_graph(G_a, start_v, end_v, delta) 
        |   |___|___G_no_dup = remove_multiedges(G) # remove multi-edges keeping the ones with the lowest weight, sort edges and remove duplicates with largest wei
        |   |___model = create_sub_ILP(model, G_a_pruned, start_v, end_v, delta, index_offset, global_var)  # create sub ILPS
 ```     
  


