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

# Extract infotmation from vcf files and ref and construct edge-labled variation graph
chmod +x vcf_extract_edit_distance.sh
./vcf_extract_edit_distance.sh



