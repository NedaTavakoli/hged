#!/bin/bash
#!/bin/bash
#****************************************
# To test the repository, use phoenix:

# Login to the cluster:

# ssh ntavakoli6@login-phoenix.pace.gatech.edu 
# ssh ntavakoli6@login-phoenix-3.pace.gatech.edu

# pace-quota
# To submit an interactive job:

#  qsub -I -q inferno -A GT-saluru8-CODA20 -l nodes=1:ppn=24,mem=100gb,walltime=96:00:00
#*****************************************

# Change this line according to your projcet directory
cd /storage/coda1/p-saluru8/0/ntavakoli6/hged
#----------------------------------------------------------------

project_dir=$(pwd)  #project top-level directory

cd ${project_dir}/data
DATA=$(pwd)

cd ${project_dir}/software
software_dir=$(pwd)

cd ${project_dir}
bcftools=${software_dir}/bcftools-1.9/bcftools
vcftools=${software_dir}/vcftools-0.1.16/bin/vcftools
tabix=${software_dir}/htslib-1.12/tabix

cd ${DATA}

#***************************************************************************************************************************
#*** first we normalize vcf files using bcftools
#* 
for id in $(seq 1 22; echo X; echo Y)
do
    $bcftools norm -m-any chr${id}.vcf.gz -Ov > chr${id}_norm.vcf
done
#***************************************************************************************************************************

#***************************************************************************************************************************
#--------------------------- JUST FOR TEST
# Just to test if bcftools works or not
#printing variants without header:
id=22
MyVariants=${DATA}/chr${id}.vcf.gz
$bcftools view -H  $MyVariants | head -2
# To count the number of SNPs in each chr
$bcftools view -H  $MyVariants | grep VT=SNP | wc -l

# https://www.biostars.org/p/298361/
# https://gist.github.com/elowy01/93922762e131d7abd3c7e8e166a74a0b
# https://www.biostars.org/p/48204/

#***************************************************************************************************************************



#***************************************************************************************************************************
#*   In the below code, we use VCFtools to count the number of SNPs and INDELs at each vcf file for each choromosome
#*   For each chromosome extract SNPs only and INDELs only
#*   In the log files # of SNPs and INDELs are kept, as well as the time that vcftools is used to extract these data
#*
#***************************************************************************************************************************
# *** Scenario 1
# *
for id in $(seq 1 22; echo X; echo Y)
do
    MyVariants=${DATA}/chr${id}.vcf.gz

    #output  out.frq.count, keep only SNPs
    $vcftools --vcf chr${id}.vcf --counts --remove-indels --out chr${id}_snps_only
    # mv out.frq.count chr${id}_snp_only.out.frq.count

    #output  out.frq.count, keep only INDEL, for these options "indel" means any variant that alters the length of the REF allele.
    # INDEl refers to insertion, deletion, or insertion and deletion of nucleotides in genomic DNA
    $vcftools --vcf chr${id}.vcf --counts --keep-only-indels --out chr${id}_indels
    #mv out.frq.count chr${id}_indels.out.frq.count
done
#                       ----------------------------------------------------------------------
# Note that: I received other types of variants as 0, so I used the following commands to extract the required fieild using shell and bcftools
# *** Scenario 2
# *
for id in $(seq 1 22; echo X; echo Y)
do
    MyVariants=${DATA}/chr${id}.vcf.gz
    echo("Total Variants", "SNPs", "Indels") > chr${id}_numVariants.txt

    # Count the total number of variants
    #$bcftools view -H  $MyVariants | wc -l 

    # # count number of SNPs
    $bcftools view -H  $MyVariants | grep VT=SNP | wc -l

    # # count number of INDLEs
    $bcftools view -H  $MyVariants | grep VT=INDEL | wc -l

    
    
    # # Not accurave because we might have combinations of insertions and deletions
    # # compute number of DEL, if len(ALT)-len(REF) < 0
    #             $bcftools view -H  $MyVariants | grep VT=INDEL | awk -F"\t" '{if(length($5)<length($4)) print $0 }' | wc -l 

    #             # compute number of DEL, if len(ALT)-len(REF) < 0
    #             $bcftools view -H  $MyVariants | grep VT=INDEL | awk -F"\t" '{if(length($5)==length($4)) print $0 }' | wc -l 
                
    #             # Find the max size of ALT and REF
    #             $bcftools view -H  $MyVariants | grep VT=INDEL | awk -F"\t" '{print length($5);} {print length($4)}' > out.txt | sort -r | head -n1 | awk '{print $2}'
            

done
#                       ----------------------------------------------------------------------
# *** Scenario 3
# *
#***** Some notes to test
# * https://gist.github.com/elowy01/93922762e131d7abd3c7e8e166a74a0b
# ** selecting snps from file:
# * $bcftools view -v snps $MyVariants | wc -l
 #*  is equivalent to $bcftools view -H  $MyVariants| $bcftools view -v snps $MyVariants | wc -l
     $bcftools view -v indels $MyVariants | wc -l 

$bcftools view -H -i 'VT="SNP"' $MyVariants | wc -l 

# selcting only SNPs in vcffiles using bcftools
$bcftools view -v snps  $MyVariants  | wc -l 

# selcting only SNPs in vcffiles using bcftools
$bcftools view -v indels  $MyVariants  | wc -l 
 $bcftools view -H -v snps  $MyVariants  | wc -l 

 #select only biallelic (excluding multiallelic) snps
  $bcftools view -m2 -M2 -v snps  $MyVariants  | wc -l 

# selcting only INDels in vcffiles using bcftools
$bcftools view -H -v indels  $MyVariants  | wc -l 



#***************************************************************************************************************************




#************************************
 # get list of variant positions (INDELS+ SNPs)
for id in $(seq 1 22; echo X; echo Y)
do
    # list of snps posiions
    cat  chr${id}_snps_only.frq.count | cut -f2 >chr${id}_snps_positions.txt
    # remove the first line (POS)
    sed '1d'  chr${id}_snps_positions.txt  > tmpfile; mv tmpfile chr${id}_snps_positions.txt

    # list of INDEL posiions
    cat  chr${id}_indels.frq.count | cut -f2 >chr${id}_indels_positions.txt
   # remove the first line (POS)    
    sed '1d'  chr${id}_indels_positions.txt  > tmpfile; mv tmpfile chr${id}_indels_positions.txt

    # Combine two colomns of positions for SNPs and INDELs
    sort  chr${id}_snps_positions.txt chr${id}_indels_positions.txt | uniq  > chr${id}_snps_indels_positions.txt
done

#**********************************
#* Compute statistics on the variant positions ( how far they are from each other)







##############$$$$$

id=22
#---------------------------
##### Code starts here
MyVariants=${DATA}/chr${id}.vcf.gz

$bcftools view $MyVariants |\
    awk -F"\t" 'BEGIN {print "CHR\tPOS\tID\tREF\tALT"} \
      !/^#/ {print $1"\t"$2"\t"$3"\t"$4"\t"$5}'


#output  out.frq.count, keep only SNPs
$vcftools --vcf chr${id}.vcf --counts --remove-indels 
mv out.frq.count chr${id}_snp_only.out.frq.count


#output  out.frq.count, keep only INDEL
$vcftools --vcf chr${id}.vcf --counts --keep-only-indels
mv out.frq.count chr${id}_indels.out.frq.count


#*******************************************************************************************************
# The below code creates a file of this format
#The obtained file has these colomns: 
#  $1   $2  $3  $4   $5     $6          $7                $8                $9            $10
# CHR  POS  ID  REF  ALT  N_ALLELE  {REF:count}  {Allele1:count}  {Allele2:count}   {Allele2:count}
#********************************************************************************************************

# Concatenate two files and sort according to positions
paste <($bcftools view $MyVariants |\
    awk -F"\t" 'BEGIN {print "CHR\tPOS\tID\tREF\tALT"} \
      !/^#/ {print $1"\t"$2"\t"$3"\t"$4"\t"$5}'\
    \
    <(cat chr${id}._indels.out.frq.count |\
        awk -F"\t"  '{print $3"\t"$5"\t"$6"\t"$7"\t"$8}') \
    <(cat chr${id}._snp_only.out.frq.count |\
        awk -F"\t"  '{print $3"\t"$5"\t"$6"\t"$7"\t"$8}') \
    | sed 's/,\t/\t/g' | sed 's/,$//g'   > chr${id}_concat.txt # Replacing all the occurrence of the pattern in a line 




# Concatenate two files and sort according to positions
paste <($bcftools view $MyVariants |\
    awk -F"\t" 'BEGIN {print "CHR\tPOS\tID\tREF\tALT"} \
      !/^#/ {print $1"\t"$2"\t"$3"\t"$4"\t"$5}'\
    \
    <(cat chr${id}._indels.out.frq.count |\
        awk -F"\t"  '{print $3"\t"$5"\t"$6"\t"$7"\t"$8}') \
    <(cat chr${id}._snp_only.out.frq.count |\
        awk -F"\t"  '{print $3"\t"$5"\t"$6"\t"$7"\t"$8}') \
    | sed 's/,\t/\t/g' | sed 's/,$//g'   > chr${id}_concat.txt # Replacing all the occurrence of the pattern in a line 






