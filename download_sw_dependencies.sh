#!/bin/bash
#****************************************
# To test the repository, use phoenix:

# Login to the cluster:

# ssh ntavakoli6@login-phoenix.pace.gatech.edu 
# ssh ntavakoli6@login-phoenix-3.pace.gatech.edu

# To submit an interactive job:

#  qsub -I -q inferno -A GT-saluru8-CODA20 -l nodes=1:ppn=24,mem=100gb,walltime=96:00:00
#*****************************************

# Change this line according to your projcet directory
cd /storage/coda1/p-saluru8/0/ntavakoli6/hged
#----------------------------------------------------------------

project_dir=$(pwd)  #project top-level directory
mkdir -p software && cd software
software_dir=$(pwd)

# get vcftools 
echo "downloading vcftools"
wget https://github.com/vcftools/vcftools/releases/download/v0.1.16/vcftools-0.1.16.tar.gz
tar xzf vcftools-0.1.16.tar.gz && cd vcftools-0.1.16
./autogen.sh
./configure --prefix=$(pwd)
make -j -s && make install
vcftools=$(pwd)
rm -f "vcftools-0.1.16.tar.gz"
echo "vcftools download and compilation finished"


# get bcftools
wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2
tar -vxjf bcftools-1.9.tar.bz2
cd bcftools-1.9
./configure  --prefix=${software_dir}/bcftools
make
bcftools=$(pwd)/bcftools
echo "bcftools download and compilation finished"
#--------------

#get htslib
#Note: tabix, htslib and bgzip2 will be installed in the bin directory and the main directory
echo "downloading htslib"
cd ${software_dir}
wget https://github.com/samtools/htslib/releases/download/1.12/htslib-1.12.tar.bz2
tar xvjf htslib-1.12.tar.bz2
cd htslib-1.12
autoreconf -i  # Build the configure script and install files it uses
./configure  --prefix=${software_dir}/htslib-1.12  # Optional but recommended, for choosing extra functionality
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
./configure   --prefix=$softwarewd/samtools-1.12        # Needed for choosing optional functionality
make
make install
samtools=$(pwd)
rm -f "samtools-1.12.tar.bz2"
echo "samtools download and compilation finished"

#get Gurobi
echo "downloading gurobi"
cd ${software_dir}
wget  https://packages.gurobi.com/9.1/gurobi9.1.0_linux64.tar.gz
tar xzf gurobi9.1.0_linux64.tar.gz
make -j -C gurobi910/linux64/src/build #re-compile gurobi cpp files using user's c++ compiler
cp gurobi910/linux64/src/build/libgurobi_c++.a gurobi910/linux64/lib
gurobi=$(pwd)
rm -f "gurobi9.1.0_linux64.tar.gz"
echo "gurobi download and compilation finished"

#check if appropriate files exist
if [ ! -f gurobi910/linux64/include/gurobi_c++.h ]; then echo "gurobi download failed"; fi
if [ ! -f gurobi910/linux64/src/build/libgurobi_c++.a ]; then echo "gurobi compilation failed"; fi
if [ ! -f vcftools-0.1.16/bin/vcftools ]; then echo "vcftools compilation failed"; fi


#--------------------------Download softwares done---------------------------------------------
#-------------------------- Dependecies done---------------------------------------------------