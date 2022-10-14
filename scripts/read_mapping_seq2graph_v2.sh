
mason=/storage/home/hcoda1/6/ntavakoli6/p-saluru8-0/hged/software/mason-0.1.2-Linux-x86_64/bin/mason
jellyfish=/storage/home/hcoda1/6/ntavakoli6/p-saluru8-0/hged/software/mason-0.1.2-Linux-x86_64/bin/jellyfish-2.2.6/bin
k8=/storage/home/hcoda1/6/ntavakoli6/p-saluru8-0/hged/software/k8/k8-linux
FORGe=/storage/home/hcoda1/6/ntavakoli6/p-saluru8-0/hged/software/FORGe
tabix=/storage/home/hcoda1/6/ntavakoli6/p-saluru8-0/software/tabix/tabix
bgzip=/storage/home/hcoda1/6/ntavakoli6/p-saluru8-0/software/tabix/bgzip
VF=storage/home/hcoda1/6/ntavakoli6/p-saluru8-0/hged/software/VF/build
hisat2=/storage/home/hcoda1/6/ntavakoli6/p-saluru8-0/hged/software/hisat-genotype-top/hisat2
GraphAligner=/storage/home/hcoda1/6/ntavakoli6/p-saluru8-0/software/GraphAligner/bin
samtools=/storage/home/hcoda1/6/ntavakoli6/p-saluru8-0/hged/software/samtools-1.12/samtools
bcftools=/storage/home/hcoda1/6/ntavakoli6/p-saluru8-0/hged/software/bcftools-1.9/bcftools
htsfile=/storage/home/hcoda1/6/ntavakoli6/p-saluru8-0/hged/software/samtools-1.12/htslib-1.12/htsfile

REF=/storage/home/hcoda1/6/ntavakoli6/p-saluru8-0/hged/data/hs37d5.fa
CHRID=chr22
VARFILE_SNP_INDELs_zipped=/storage/home/hcoda1/6/ntavakoli6/p-saluru8-0/hged/data/chr22_snps_indels.vcf.gz



# get list of retained varinat positions
python src/solution_analyzer.py graph_chr22_10k_100.txt \
    ILP_sol_vectors/ILP_sol_100_1.txt retained_variants_ed_ILP_100_1.txt  

#get variant retained vcf file
# do only once
cat chr22_snps_indels.vcf | grep '^#' > header_chr22.vcf
cat chr22_snps_indels.vcf | grep  -vE '^#' > non_header_chr22.txt

# do for each config
cp header_chr22.vcf retained_variants_ed_ILP_100_1_chr22.vcf
awk 'NR == FNR {a[$0]; next } $2 in a {print $0} ' retained_variants_ed_ILP_100_1.txt  non_header_chr22.txt >> retained_variants_ed_ILP_100_1_chr22.vcf

reduced_vcf=retained_variants_ed_ILP_100_1_chr22.vcf
# to check the content
# $bcftools view retained_variants_ed_ILP_100_10_chr22.vcf

# get haplotype substrings
cat pos_substrings_chr22_10k_100.txt | head -100 | awk '{printf $2}' > reads.${CHRID}.fastq 



#*************************************************
# Map on comptete vcf files

# construct graph from snps and indels
./vg construct -r $REF -v chr22_snps_indels.vcf.gz -a -f -m 32 > graph.${CHRID}.vg
#index the graph
./vg index -x graph_${CHRID}.xg -g graph_${CHRID}.gcsa graph.${CHRID}.vg


# Map on complete variation graph and vcf 
#Run FORGe
mkdir -p  chr22_hisat_index

$FORGe/src/vcf_to_1ksnp.py  --reference $REF --vcf chr22_snps_indels.vcf --out ${CHRID}_snp.1ksnp 
$FORGe/src/rank.py --method popcov --reference $REF  --vars ${CHRID}_snp.1ksnp   --window-size 100  --output ordered_${CHRID}_snp.txt
$FORGe/src/build.py --reference  $REF --vars ${CHRID}_snp.1ksnp --window-size 100 --hisat hisat_input_${CHRID}_snp.snp --sorted ordered_${CHRID}_snp.txt --pct 10
$hisat/hisat2-build --snp hisat_input_${CHRID}_snp.snp $REF ${CHRID}_hisat_index/index_${CHRID}

#map using HISAT2: haplotype reads to graph indexes
$hisat/hisat2 -f -x $FORGe/${CHRID}_hisat_index/index_${CHRID} -U reads.${CHRID}.fastq -S hisat_alignment_${CHRID}_snp.sam #reads.${CHRID}.fastq are haplotype reads

#map by graph aligner

$GraphAligner/GraphAligner -g graph.${CHRID}.vg -f reads.${CHRID}.fastq -a GraphAligner.${CHRID}.aligned.gam -x vg
./vg surject -s GraphAligner.${CHRID}.aligned.gam -x graph.${CHRID}.xg > GraphAligner.${CHRID}.aligned.sam #linearize

#*************************************************
# Map on reduced vcf files

# Generate vcf.gz file and its index file vcf.gz.tbi
cp retained_variants_ed_ILP_100_1_chr22.vcf retained_variants_ed_ILP_100_1_chr22_copy.vcf
$bgzip -c retained_variants_ed_ILP_100_1_chr22.vcf > retained_variants_ed_ILP_100_1_chr22.vcf.gz
$tabix -p vcf retained_variants_ed_ILP_100_1_chr22.vcf.gz

# construct graph from snps and indels
./vg construct -r $REF -v retained_variants_ed_ILP_100_1_chr22_copy.vcf -a -f -m 32 > reduced_graph.${CHRID}_100_1.vg
#index the graph
./vg index -x graph_${CHRID}.xg -g graph_${CHRID}.gcsa graph.${CHRID}.vg

 #Run our hged framework for Chr #ID

mkdir -p  chr22_hisat_index_hged
$FORGe/src/vcf_to_1ksnp.py  --reference $REF --vcf retained_variants_ed_ILP_100_1_chr22.vcf --out ${CHRID}_snp_VF.1ksnp 
$FORGe/src/rank.py --method popcov-blowup --reference $REF --vars ${CHRID}_snp_VF.1ksnp --window-size 100 --output ordered_${CHRID}_snp_VF.txt
$FORGe/src/build.py --reference  $REF --vars ${CHRID}_snp_VF.1ksnp --window-size 100 --hisat hisat_input_${CHRID}_snp_VF.snp --sorted ordered_${CHRID}_snp_VF.txt --pct 100
$hisat2/hisat2-build --snp hisat_input_${CHRID}_snp_VF.snp $REF ${CHRID}_hisat_index_hged/chr22_index_hged
echo "FORGe with pct=100 for VF variants is done"

#map using HISAT2: simulated reads to graph indexes
$hisat2/hisat2 -f -x ${CHRID}_hisat_index_hged/chr22_index_hged -U reads_${CHRID}_${LEN}.fasta -S hisat_alignment_${CHRID}_snp_hged.sam 
echo "HISAT2 mapping for VF is done"


#map by graph aligner
$GraphAligner/GraphAligner -g reduced_graph.${CHRID}_100_1.vg -f reads.${CHRID}.fastq -a GraphAligner.${CHRID}.aligned.gam -x vg
./vg surject -s GraphAligner.${CHRID}.aligned.gam -x graph.${CHRID}.xg > GraphAligner.${CHRID}.aligned.sam #linearize

#*************************************************

# get statistics
$samtools flagstat GraphAligner.${CHRID}.aligned.sam #linearize