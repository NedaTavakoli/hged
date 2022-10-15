
# 16050075 16654125
ssh ntavakoli6@login-phoenix-3.pace.gatech.edu
qsub -I -q inferno -A GT-saluru8-CODA20 -l nodes=1:ppn=24,mem=100gb,walltime=96:00:00

samtools=/storage/home/hcoda1/6/ntavakoli6/p-saluru8-0/hged/software/samtools-1.12/samtools
bcftools=/storage/home/hcoda1/6/ntavakoli6/p-saluru8-0/hged/software/bcftools-1.9/bcftools
#***********************************************************
# get dependencies
$bcftools view -v 'snps,indels' -Oz chr22.vcf.gz > chr22_snps_indels.vcf.gz
$bcftools index chr22_snps_indels.vcf.gz
$bcftools view -H -v 'snps,indels' chr22.vcf.gz | awk '{print $2}' > variant_positions_snps_indels_chr22.txt
$bcftools view -H -v 'snps,indels' chr22.vcf.gz | awk '{print $2"\t"$4"\t"$5}' > chr22_snps_indel_POS_REF_ALT.txt
$bcftools query -l chr22_snps_indels.vcf.gz > samples.txt

#***********************************************************
# get variant retained vcf file

python src/solution_analyzer.py graph_chr22_10k_100.txt \
    ILP_sol_vectors/ILP_sol_100_1.txt retained_variants_ed_ILP_100_1.txt    

# do only once
cat chr22_snps_indels.vcf | grep '^#' > header_chr22.vcf
cat chr22_snps_indels.vcf | grep  -vE '^#' > non_header_chr22.txt

# do for each config
cp header_chr22.vcf retained_variants_ed_ILP_100_1_chr22.vcf
awk 'NR == FNR {a[$0]; next } $2 in a {print $0} ' retained_variants_ed_ILP_100_1.txt  non_header_chr22.txt >> retained_variants_ed_ILP_100_1_chr22.vcf

# Generate vcf.gz file and its index file vcf.gz.tbi
cp retained_variants_ed_ILP_100_1_chr22.vcf retained_variants_ed_ILP_100_1_chr22_copy.vcf
$bgzip -c retained_variants_ed_ILP_100_1_chr22.vcf > retained_variants_ed_ILP_100_1_chr22.vcf.gz
$tabix -p vcf retained_variants_ed_ILP_100_1_chr22.vcf.gz
#***********************************************************

# READ MAPPING startrs from here

# construct list of random;y chosen haplotypes at range and location
# location: select a random number from a variant position list 

num_reads=50

cat data/variant_positions_snps_indels_chr22.txt | head -${num_reads} > data/variant_positions_snps_indels_chr22_head_${num_reads}.txt
variant_positions=($(cut -d ',' -f2 data/variant_positions_snps_indels_chr22_head_${num_reads}.txt))
samples=($($bcftools query -l data/chr22_snps_indels.vcf.gz)) # array of samples, index from 0

count=0
for v in "${variant_positions[@]}"
do
    # echo "Location"$v
    sample=${samples[$count]} # random sample
    count=$(($count+1))
    # echo "Haplotype: "$sample
    s1=$($samtools faidx data/hs37d5.fa 22:${v}-$((${v}+100-1)) | $bcftools consensus -s ${sample} -H 2 data/chr22_snps_indels.vcf.gz |  sed '1d' | tr -d "[:space:]")
    # echo $s1
    a=($v $sample $s1) 
    printf '%s\n' "${a[*]}" >> pos.hap.reads.chr22.txt  
done
cat pos.hap.reads.chr22.txt  | awk '{print $3}' > list.reads.chr22_${num_reads}.txt 


#----------
# Map on complete graph

# construct graph from snps and indels
./vg construct -r data/hs37d5.fa -v chr22_snps_indels.vcf.gz -a -f -m 32 > graph.chr22.vg
#index the graph
./vg index -x graph.chr22.xg -g graph.chr22.gcsa graph.chr22.vg

-bash-4.2$ cat pos_substrings_chr22_10k_100.txt  | head -1
16050075 GGTGGGCCTAAGTGCCTCCTCTCGGGACTGGTATGGGGACGGTCATGCAA

# to map a single read
./vg map -s GGTGGGCCTAAGTGCCTCCTCTCGGGACTGGTATGGGGACGGTCATGCAA -x graph.chr22.xg -g graph.chr22.gcsa > read.chr22.gam
./vg surject -x graph.chr22.xg -b read.chr22.gam > read.chr22.gam.bam
$samtools flagstat read.chr22.gam.bam


# to map lists of reads
./vg map -T list.reads.chr22.txt -x graph.chr22.xg -g graph.chr22.gcsa > list.read.chr22.gam
./vg surject -x graph.chr22.xg -b list.read.chr22.gam > list.read.chr22.gam.bam
$samtools flagstat list.read.chr22.gam.bam

#----------
# on reduced graph
./vg construct -r data/hs37d5.fa -v retained_variants_ed_ILP_100_1_chr22.vcf.gz -a -f -m 32 > reduced_100_1.graph.chr22.vg
# warning:[vg::Constructor] Unsupported IUPAC ambiguity codes found in 3; coercing to N.
# warning:[vg::Constructor] Lowercase characters found in hs37d5; coercing to uppercase.
# warning:[vg::Constructor] Unsupported IUPAC ambiguity codes found in hs37d5; coercing to N.

#index the graph
./vg index -x  reduced_100_1.graph.chr22.xg -g  reduced_100_1.graph.chr22.gcsa  reduced_100_1.graph.chr22.vg

# to map a single read
./vg map -s GGTGGGCCTAAGTGCCTCCTCTCGGGACTGGTATGGGGACGGTCATGCAA -x reduced_100_1.graph.chr22.xg -g reduced_100_1.graph.chr22.gcsa > reduced_100_1.read.chr22.gam
./vg surject -x reduced_100_1.graph.chr22.xg -b reduced_100_1.read.chr22.gam > reduced_100_1.read.chr22.gam.bam
$samtools flagstat reduced_100_1.read.chr22.gam.bam


# to map lists of reads
./vg map -T list.reads.chr22.txt -x reduced_100_1.graph.chr22.xg -g reduced_100_1.graph.chr22.gcsa > reduced_100_1.list.read.chr22.gam
./vg surject -x reduced_100_1.graph.chr22.xg -b reduced_100_1.list.read.chr22.gam > reduced_100_1.list.read.chr22.gam.bam
$samtools flagstat reduced_100_1.list.read.chr22.gam.bam


