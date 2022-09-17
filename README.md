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
chmod +x 1.download_vcf_and_ref.sh
./1.download_vcf_and_ref.sh

# Download software dependencies
chmod +x 2.download_sw_dependencies.sh
./2.download_sw_dependencies.sh

# Extract infotmation from vcf files: variant positions for SNPs and INDELs
chmod +x 3.variant_positions_INDEL_SNPs.sh
./3.variant_positions_INDEL_SNPs

# Note : git reset HEAD^ --hard

# get POS REF ALT per variant position (required to construct graph)
chmod +x 4.vcf_of_snps_indels_POS_REF_ALT.sh
./4.vcf_of_snps_indels_POS_REF_ALT.sh

 # Get the linear backbone per chromosome (required to construct the graph)
 chmod +x 5.linear_ref_backbone_per_chromosome.sh
 ./5.linear_ref_backbone_per_chromosome.sh

# Construt edge-lable variation graph
chmod vg_alter_paths.sh
./vg_alter_paths.sh

# get the substring of length alpha for haplotype h1/or h2 for a given sample
chmod +x 7.get_string_of_haplotype_from_bp_to_bp.sh
./7.get_string_of_haplotype_from_bp_to_bp.sh
