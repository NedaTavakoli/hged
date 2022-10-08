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
#include "gurobi_c++.h"

/********* Helper functions ******/

/**
 * @brief   Extract substrings of length alpha at each variant position. 
 *          This function outputs substrings of length alpha at each variant position.
 */

void extract_pos_substring (const std::string &vcf_file, const std::string &fasta_file, const std::string &pos_file, const int &alpha, std::vector<int> &ivariant_pos)
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
    
}

int main(int argc, char **argv) {

  //parse command line arguments
  Parameters parameters;
  parseandSave_ILP(argc, argv, parameters);

  std::vector<int> variant_pos; 
  extract_pos_substring (parameters.vcffile, parameters.fasta_ref_file, parameters.pos_file, parameters.alpha, variant_pos)
  assert (std::is_sorted(variant_pos.begin(), variant_pos.end())); //must be sorted in ascending order

  //variant positions (unique values)
  std::vector<int> pos_u;
  pos_u.insert (pos_u.end(), variant_pos.begin(), variant_pos.end());
  std::sort (pos_u.begin(), pos_u.end());
  pos_u.erase(std::unique(pos_u.begin(), pos_u.end()), pos_u.end() );

  std::cout<< "INFO, hged::main, count of variant containing positions = " << pos_u.size() << "\n";
  std::cout<< "INFO, hged::main, count of indels = " << indelpos.size() << "\n";
  std::cout<< "INFO, hged::main, count of SNP variants = " << std::accumulate(snpcount.begin(), snpcount.end(), 0) << "\n";

  
  
  
  
  
  
  
  