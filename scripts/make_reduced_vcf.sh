# get variant retained vcf file
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