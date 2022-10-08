GUROBI_INSTALL=$(shell pwd)/build/gurobi910/linux64
BCFTOOLS_INSTALL=$(shell pwd)/build/bcftools-1.9/bcftools
SAMTOOLS_INSTALL=$(shell pwd)/build/samtools-1.12/samtools
VCFTOOLS_INSTAL=$(shell pwd)/build/vcftools-0.1.16/bin/vcftools
TABIX_INSTAL=$(shell pwd)/build/htslib-1.12/tabix
TARGET_DIR=$(shell pwd)/build
CPPFLAGS= -g -std=c++11 -DNDEBUG -O3 

all:
	mkdir -p build
	$(CXX) $(CPPFLAGS) -D VCFTOOLSPATH=$(VCFTOOLS_INSTALL) -m64 -o $(TARGET_DIR)/ilp_snp_indels -I $(GUROBI_INSTALL)/include/ -L  $(GUROBI_INSTALL)/lib/ src/ilp_snp_indels.cpp -lgurobi_c++ $(GUROBI_INSTALL)/lib/libgurobi91.so -lm
	@echo "check executables in build directory"


clean:
	rm -f $(TARGET_DIR)/ilp_snp_indels
