## Haplotype_aware Variation Selection in Genome Graphs under Edit Distance

## Dependencies
- A C++ compiler with c++11 support, e.g., GNU g++ (version 5+)
- [samtools](https://vcftools.github.io/)
- [bcftools](https://vcftools.github.io/)
- [Gurobi](https://www.gurobi.com)


## Installation
The above dependencies can be handled by running script `dependencies.sh`.

### To use Gurobi:
Get your free Gurobi academic license code by registering here: https://www.gurobi.com/downloads/end-user-license-agreement-academic"
Add your licence key by using build/gurobi910/linux64/bin/grbgetkey tool"

The overall workflow is:

```sh
git clone https://github.com/NedaTavakoli/hged
cd hged
project_dir=$(pwd)  #project top-level directory
chr_id=22 #* change this numbers according to your needs
alpha=150 #* change this numbers according to your needs
# download data and softwares
chmod +x dependencies.sh
./dependencies.sh ${chr_id} 
```

After a successful compilation, expect executables named as `ilp_snp_indels` in a directory named `build`.

## Usage
All the executables implement a variety of algorithms to achieve variant graph size reduction, but they all have a similar interface.
```
SYNOPSIS

        ilp_snp_indels    -a <alpha> -d <delta>  -chr <id> -vcf <file1>  -fa <file2>  -pos <file3> [-prefix <file4>]


OPTIONS
        <alpha>     path length in variation graph (e.g., 500)
        <delta>     differences allowed (e.g., 10)
        <id>        chromosome id (e.g., 22), make it consistent with vcf file
        <file1>     compressed vcf file (something.vcf.gz)
        <file2>     reference genome fasta file (something.fa)
        <file3>     variant position file for SNPs and INDELs (something.txt)
        <file4>     filename to optionally save input and output variants
```

This repository is used to solve variant selection in genome graphs under edit disatnce
in other words for a given complete variation graph for each chromosome, it creates a reduced variation graph in which 
some variants are removed subject to some constraints. The constraints are for every substring of length 
alpha observed in haplotypes, the reduced varaition graph guarantees to preserve those substrings with
at most delta errors (i.e., edit distance of delta among alpha-long substrings of haplotypes in complete variation graph with those of reduced variation graph).

The project has the following folder structure:
```
hged
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
  


