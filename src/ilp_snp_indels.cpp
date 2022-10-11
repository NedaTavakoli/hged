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

namespace sp = subprocess;
  
/**
 * @brief   Extract substrings of length alpha at each variant position. 
 *          This function outputs substrings of length alpha at each variant position.
 */
void extract_pos_substring (const std::string &vcf_file, const std::string &fasta_file, const std::string &pos_file, const int &alpha, const int &chr, std::vector<int> &variant_pos, std::vector<std::string> &samples)
{
   // Extract variant positions from variant position file
   // seed random generator by time in seconds (this may create issue if two instances are launched at the same time)
    srand(time(0)); int random = rand() % 100000;  
    std::string tmp_file = ".variant_positions." + std::to_string(random) + ".txt";
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
      std::string col1;
      std::istringstream iss(line2);
      iss >> col1;
      samples.push_back(col1);
    }

    // Extract substrings of lengh alpha for each variant position and save the results in pos_subtring_alpha_chr
    srand(time(0)); int random3 = rand() % 100000;
    std::ofstream fout("pos_substring_alpha_" + std::to_string(alpha) + "_chr_" + std::to_string(chr) + "." + std::to_string(random3) + ".txt");
      
    std::string tmp_file3 = "pos_substrings_unsorted." + std::to_string(random3) + ".txt";
    for (std::size_t i = 0; i <= variant_pos.size(); i++) // vector in C++ is 0-indexed
    { 
    //int i = 1;
   
      // read a line of vcf file at the variant position i, save to tmp_file3, then save to the vector named "vecOfStrs"
   
      //
      std::string cmd3 = std::string(TOSTRING(BCFTOOLSPATH)) + " view -H -r " + std::to_string(chr) + ":" + std::to_string(variant_pos[i]) + " " + vcf_file + " >  " + tmp_file3;
      std::system(cmd3.c_str());

      std::ifstream file3 (tmp_file3);
      std::string line3;

      std::vector<std::string> vecOfStrs;
      while (std::getline(file3, line3))
      {
        std::istringstream ss(line3);
        vecOfStrs.push_back(line3);
      }

      // split the line of each variant position by space and save the results to the vector named "vstrings"
      std::stringstream ss(vecOfStrs[0]);  // vecOfStrs[0] is the whole line of vcf file corresponds to a variant position 
      std::istream_iterator<std::string> begin(ss);
      std::istream_iterator<std::string> end;
      std::vector<std::string> vstrings(begin, end);

      // extract sequences correspond to each haplotype in a vector named list_haplotypes
      std::vector<std::string> list_haplotypes;  
      std::vector<int> sample_index;  

      // get the samples correspond to each GT number, note the colomn of GT starts from the colomn 9 (0-indexed)    
      for (std::size_t j = 9; j <= vstrings.size(); j++)  // j is the sample index
      {
        if (vstrings[j].compare("0|0") != 0)  // only extract haplotype sequences of GT != 0|0
        {
          sample_index.push_back(j-9);

          int h1 = 1, h2 = 2 ;

          // vstrings[i][0] represents the first element of GT: ex for 1|0, it shows 1
          // vstrings[i][2] represents the third element of GT: ex for 1|0, it shows 0 ( the second element is |)
          if ( vstrings[j][0] != 0)  
          {  int t = variant_pos[i] + alpha -1;

            auto sam = sp::Popen({std::string(TOSTRING(SAMTOOLSPATH)), "faidx", fasta_file , std::to_string(chr) + ":" + std::to_string(variant_pos[i]) + "-" + std::to_string(t)}, sp::output{sp::PIPE});
            auto consensus = sp::Popen({std::string(TOSTRING(BCFTOOLSPATH)), "consensus", "-s", samples[j-9], "-H", std::to_string(h1), vcf_file}, sp::input{sam.output()}, sp::output{sp::PIPE});
            auto res = consensus.communicate().first;
            // std::cout << "haplotype 1" << res.buf.data() << std::endl;
            list_haplotypes.push_back(res.buf.data());
          }

          if ( vstrings[j][2] != 0)
          {  int t = variant_pos[i] + alpha -1;

            auto sam = sp::Popen({std::string(TOSTRING(SAMTOOLSPATH)), "faidx", fasta_file , std::to_string(chr) + ":" + std::to_string(variant_pos[i]) + "-" + std::to_string(t)}, sp::output{sp::PIPE});
            auto consensus = sp::Popen({std::string(TOSTRING(BCFTOOLSPATH)), "consensus", "-s", samples[j-9], "-H", std::to_string(h2), vcf_file}, sp::input{sam.output()}, sp::output{sp::PIPE});
            auto res = consensus.communicate().first;
            //std::cout << "haplotype 2" << res.buf.data() << std::endl;
            list_haplotypes.push_back(res.buf.data());
          }

        }

      }

      
      // // From the PIP results, ignore the first line, only keep unique values, remove spaces and make sure the length of the substring is alpha
      // //list of haplotypes1 (unique values)
      std::vector<std::string> hap_u;
      hap_u.insert (hap_u.end(), list_haplotypes.begin(), list_haplotypes.end());
      std::sort (hap_u.begin(), hap_u.end());
      hap_u.erase(std::unique(hap_u.begin(), hap_u.end()), hap_u.end() );

      // // save the results of std::cout into file named coutbuf
      
      fout << variant_pos[i];
      for (std::size_t k = 0; k < hap_u.size(); k++)
      {
        std::size_t X = hap_u[k].find_first_of("ACGT");
        //remove(hap_u[k].substr(X).begin(), hap_u[k].substr(X).end(), ' ');  //remove spaces
        //std::cout << variant_pos[i] << " " << hap_u[k].substr(X) << std::endl; // from the first char to the end
        fout <<  " " + hap_u[k].substr(X);
        
      } 
      fout << std::endl;
    }   

  }

    




int main(int argc, char **argv) {

  //parse command line arguments
  Parameters parameters;
  parseandSave_ILP(argc, argv, parameters);

  std::vector<int> variant_pos; 
  std::vector<std::string> samples; 
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

  
     
        

     // std::size_t X = hap_u[i].find_first_of("ACGT");   // to skip the first numbers in the line
          // std::cout << "The first position of the String" << X << std::endl;
          // auto end2 = hap_u[i].end();
  
  
  
  