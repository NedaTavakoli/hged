#!/bin/bash

if [ $# -eq 0 ];
then
    echo "$0: Missing arguments"
    exit 1
elif [ $# -gt 1 ];
then
    echo "$0: Too many arguments: $@"
    exit 1
else
    echo "Construct edge-labeled varion graph ( args: chr id, length of substrings"
    echo "==========================="
    echo "Number of arguments.: $#"
    echo "List of arguments...: $@"
    echo "Arg #1 chrId (Ex: 22)..................................: $1"
    echo "==========================="
fi

id=$1

project_dir=$(pwd)
mkdir -p data && cd data
DATA=$(pwd)

# Get the reference
REFERENCE=hs37d5.fa
rm -f ${REFERENCE} ${REFERENCE}.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/${REFERENCE}.gz
gunzip ${REFERENCE}

# Get the phasings for chromosomes 
PREFIX=ALL.chr
SPREFIX=chr
SUFFIX=.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
SSUFFIX=.vcf.gz
for i in ${id}
do
  NAME=${PREFIX}${i}${SUFFIX}
  SHORT=${SPREFIX}${i}${SSUFFIX}
  rm -f ${SHORT} ${SHORT}.tbi
  wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/${NAME}
  wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/${NAME}.tbi
  mv ${NAME} ${SHORT}
  mv ${NAME}.tbi ${SHORT}.tbi
done

# Download softwares
cd ..
project_dir=$(pwd)  #project top-level directory
mkdir -p build && cd build
software_dir=$(pwd)

# get bcftools
echo "downloading bcftools"
cd ${software_dir}
wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2
tar -vxjf bcftools-1.9.tar.bz2
cd bcftools-1.9
./configure  --prefix=${software_dir}/bcftools
make
bcftools=$(pwd)/bcftools
echo "bcftools download and compilation finished"

#get htslib
#Note: tabix, htslib and bgzip2 will be installed in the bin directory and the main directory
echo "downloading htslib"
cd ${software_dir}
wget https://github.com/samtools/htslib/releases/download/1.12/htslib-1.12.tar.bz2
tar xvjf htslib-1.12.tar.bz2
cd htslib-1.12
autoreconf -i  # Build the configure script and install files it uses
./configure  --prefix=${software_dir}/htslib-1.12  # Optional but recommended, for choosing extra functionality, --disable-lzma
make
make install
htslib=$(pwd)
bgzip2=$(pwd)
tabix=$(pwd)
rm -f "htslib-1.12.tar.bz2"
echo "htslib download and compilation finished"

#get samtools
echo "downloading SAMtools"
cd ${software_dir}
wget https://github.com/samtools/samtools/releases/download/1.12/samtools-1.12.tar.bz2
tar -xvf samtools-1.12.tar.bz2
cd samtools-1.12
autoheader            # Build config.h.in (this may generate a warning about # AC_CONFIG_SUBDIRS - please ignore it).
autoconf -Wno-syntax  # Generate the configure script
./configure  --prefix=${software_dir}/samtools-1.12        # Needed for choosing optional functionality, --without-curses for compbio
make
make install
samtools=$(pwd) # samtools=$(pwd)/samtools
rm -f "samtools-1.12.tar.bz2"
echo "samtools download and compilation finished"

#get Gurobi
echo "downloading gurobi"
cd ${software_dir}
wget  https://packages.gurobi.com/9.1/gurobi9.1.0_linux64.tar.gz
tar xzf gurobi9.1.0_linux64.tar.gz
cd gurobi910/linux64/
make -j -C gurobi910/linux64/src/build #re-compile gurobi cpp files using user's c++ compiler
cp gurobi910/linux64/src/build/libgurobi_c++.a gurobi910/linux64/lib
# GUROBI_HOME=$(pwd)
# python setup.py install #re-compile gurobi python files using user's python
# cp gurobi910/linux64/src/build/python${ver}_utf32 gurobi910/linux64/lib
# export PYTHONPATH=$GUROBI_HOME/lib/python${ver}_utf32:$PYTHONPATH
# gurobi=${software_dir}/gurobi910/linux64/lib/python${ver}_utf32
rm -f "gurobi9.1.0_linux64.tar.gz"
echo "gurobi download and compilation finished"

cd ${DATA}
$bcftools view -v 'snps,indels' -Oz chr${id}.vcf.gz > chr${id}_snps_indels.vcf.gz
$bcftools index chr${id}_snps_indels.vcf.gz
$bcftools view -H -v 'snps,indels' chr${id}.vcf.gz | awk '{print $2}' > variant_positions_snps_indels_chr${id}.txt
$bcftools view -H -v 'snps,indels' chr${id}.vcf.gz | awk '{print $2"\t"$4"\t"$5}' > chr${id}_snps_indel_POS_REF_ALT.txt
Num_variant=$(cat chr${id}_snps_indel_POS_REF_ALT.txt | wc -l)
echo "Number of variant position is: "$Num_variant
$bcftools query -l chr${id}.vcf.gz > samples.txt
num_samples=$(cat samples.txt | wc -l)
num_haplotypes=$(($num_samples*2))
echo 'Number of haplotypes:'$num_haplotypes
echo "Looks like it went okay, now run  chmode +x ./scripts/construct_graph.sh && ./scripts/construct_graph.sh"
