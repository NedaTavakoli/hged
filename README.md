# hged
Haplotype_aware variation graphs_edit disatance

------------------------------------------------------------------
To test the repository, use phoenix:

Login to the cluster:

ssh ntavakoli6@login-phoenix.pace.gatech.edu 
OR
ssh ntavakoli6@login-phoenix-3.pace.gatech.edu

To submit an interactive job:

 qsub -I -q inferno -A GT-saluru8-CODA20 -l nodes=1:ppn=24,mem=100gb,walltime=96:00:00

 # Change this line according to your projcet directory
cd /storage/coda1/p-saluru8/0/ntavakoli6/hged

# Download vcf files and Human referene genome
chmod +x download_vcf_and_ref.sh
./download_vcf_and_ref.sh

# Download software dependencies
chmod +x download_sw_dependencies.sh
./download_sw_dependencies.sh

# Extract infotmation from vcf files: variant positions for SNPs and INDELs
chmod +x vcf_extract_edit_distance.sh
./vcf_extract_edit_distance.sh

# Note : git reset HEAD^ --hard

# get POS REF ALT per variant position (required to construct graph)
chmod +x get_vcfinfo_of_snps_indels_pos_ref_alt.sh
./get_vcfinfo_of_snps_indels_pos_ref_alt.sh
 
# Construt edge-lable variation graph
chmod vg_alter_paths.sh
./vg_alter_paths.sh

# get the substring of length alpha for haplotype h1/or h2 for a given sample
chmod +x get_string_of_haplotype_from_bp_to_bp.sh
./get_string_of_haplotype_from_bp_to_bp.sh
