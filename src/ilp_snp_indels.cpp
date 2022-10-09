#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector> 
#include <algorithm>
#include <ctime>
#include <cstring>
#include <chrono>
#include <numeric>
#include <cassert>
#include <unordered_map>
#include "ext/prettyprint.hpp"
#include "common.hpp"
#include <boost/graph/graph_traits.hpp>
#include "gurobi_c++.h"
#include "ext/subprocess.hpp"

/********* Helper functions ******/
// using namespace subprocess;

/**
 * @brief   Extract substrings of length alpha at each variant position. 
 *          This function outputs substrings of length alpha at each variant position.
 */
void extract_pos_substring (const std::string &vcf_file, const std::string &fasta_file, const std::string &pos_file, const int &alpha, const int &chr, std::vector<int> &variant_pos, std::vector<int> &samples)
{
   // Extract variant positions from variant position file
   // seed random generator by time in seconds (this may create issue if two instances are launched at the same time)
    srand(time(0)); int random = rand() % 100000;  
    std::string tmp_file = ".hged." + std::to_string(random) + ".txt";
    std::string cmd = " cut -f2 " + pos_file + " >  " + tmp_file;
    std::cout << "INFO, hged::main, extracting pos from variant position file using command: " << cmd << std::endl;
    std::system(cmd.c_str());

    std::ifstream file (tmp_file);
    std::string line;
    while (std::getline(file, line))
    {
      int col1;

      std::istringstream iss(line);
      iss >> col1;
      variant_pos.push_back(col1);
    }

   // Extract samples from vcf file
    std::string tmp_file2 = "sample." + std::to_string(random) + ".txt";
    std::string cmd2 = std::string(TOSTRING(BCFTOOLSPATH)) + " query -l " + vcf_file + " >  " + tmp_file2;
    // std::cout << "INFO, hged::main, extracting samples from vcf file using command: " << cmd2 << std::endl;
    std::system(cmd2.c_str());

    std::ifstream file2 (tmp_file2);
    std::string line2;
    while (std::getline(file2, line2))
    {
      int col1;

      std::istringstream iss(line2);
      iss >> col1;
      samples.push_back(col1);
    }


      int i =1;
      //  arr=($($bcftools view -H -r 22:${v} chr${id}.vcf.gz| awk -F"\t" '{split($0, header, "\t");} \
        {for (i=10; i<=NF; i++) {if (gsub(/0\|1|1\|0|0\/1|1\/0/, "", $(i))==1) {printf header[i]",";printf i-10"\t"} if (i==NF) {printf "\n"}}}'))  
      srand(time(0)); int random3 = rand() % 100000;  
      std::string tmp_file3 = "pos_substrings_unsorted." + std::to_string(random3) + ".txt";

      std::string cmd3 = std::string(TOSTRING(BCFTOOLSPATH)) + " view  -H -r " + std::string(TOSTRING(chr)) + ":" + std::string(TOSTRING(variant_pos[i])) + " " + vcf_file + " >  " + tmp_file3;
      std::system(cmd3.c_str());

      // std::string cmd3 = std::string(TOSTRING(BCFTOOLSPATH)) + " view  -H -r " + std::string(TOSTRING(chr)) + ":" + std::string(TOSTRING(variant_pos[i])) + " " + vcf_file + \
      //   '|' + "awk -F'\t' '{split($0, header,'\t');} \
      //    {for (i=10; i<=NF; i++) {if (gsub(/0'\'|1|1'\'|0|0'\'/1|1'\'/0/, "", $(i))==1) {printf header[i]'\t';printf i-10'\t'} if (i==NF) {printf '\n'}}}'" + " >  " + tmp_file3;
      // std::system(cmd3.c_str());

    
 
}

int main(int argc, char **argv) {

  //parse command line arguments
  Parameters parameters;
  parseandSave_ILP(argc, argv, parameters);

  std::vector<int> variant_pos; 
  std::vector<int> samples; 
  extract_pos_substring (parameters.vcffile, parameters.fasta_ref_file, parameters.pos_file, parameters.alpha, parameters.chr, variant_pos, samples);
  assert (std::is_sorted(variant_pos.begin(), variant_pos.end())); //must be sorted in ascending order

  //variant positions (unique values)
  std::vector<int> pos_u;
  pos_u.insert (pos_u.end(), variant_pos.begin(), variant_pos.end());
  std::sort (pos_u.begin(), pos_u.end());
  pos_u.erase(std::unique(pos_u.begin(), pos_u.end()), pos_u.end() );

  std::cout<< "INFO, hged::main, count of variant containing positions = " << pos_u.size() << "\n";
  std::cout<< "INFO, hged::main, number of samples in the chromosome = " << samples.size() << "\n";

  return 0;
}

  
  
  
  
  
  
  
  