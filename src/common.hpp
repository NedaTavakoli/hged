#include <cassert>
#include "ext/clipp.h"

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

struct Parameters
{
  int alpha;
  int delta;
  int chr;
  std::string vcffile;
  std::string fasta_ref_file;
  std::string pos_file;
  std::string prefix;
};

inline bool exists (const std::string& filename) {
  std::ifstream f(filename.c_str());
  return f.good();
}

/**
 * @brief  parse and print command line arguments
 */
void parseandSave(int argc, char** argv, Parameters &param)
{
    param.pos = false; //default

  //define all arguments
  auto cli =
    (
     clipp::required("-a") & clipp::value("alpha", param.alpha).doc("path length in variation graph (e.g., 500)"),
     clipp::required("-d") & clipp::value("delta", param.delta).doc("differences allowed (e.g., 10)"),
     clipp::required("-chr") & clipp::value("chr", param.chr).doc("choromosome id (e.g., 22), make it consistent with vcf file"),
     clipp::required("-vcf") & clipp::value("file1", param.vcffile).doc("compressed vcf file (something.vcf.gz)"),
     clipp::required("-fa") & clipp::value("file2", param.fasta_ref_file).doc("reference genome fasta file (something.fa)"),
     clipp::required("-fa") & clipp::value("file2", param.pos).doc("reference genome fasta file (something.fa)"),
     clipp::option("-prefix") & clipp::value("file3", param.prefix).doc("filename to optionally save input and output variants")
    );

  if(!clipp::parse(argc, argv, cli))
  {
    //print help page
    clipp::operator<<(std::cout, clipp::make_man_page(cli, argv[0])) << std::endl;
    exit(1);
  }

  std::cout << "INFO, hged::parseandSave, alpha = " << param.alpha << std::endl;
  std::cout << "INFO, hged::parseandSave, delta = " << param.delta << std::endl;
  std::cout << "INFO, hged::parseandSave, chromosome id = " << param.chr << std::endl;
  std::cout << "INFO, hged::parseandSave, vcf file = " << param.vcffile << std::endl;
  std::cout << "INFO, hged::parseandSave, fasta file = " << param.fasta_ref_file << std::endl;
  if (param.prefix.length() > 0) std::cout << "INFO, VF::parseandSave, prefix = " << param.prefix << std::endl;

  if (! exists(param.vcffile))
  {
    std::cerr << "ERROR, hged::parseandSave, vcf file cannot be opened" << std::endl;
    exit(1);
  }
}

/**
 * @brief  parse and print command line arguments (modified for ILP)
 */
void parseandSave_ILP(int argc, char** argv, Parameters &param)
{
    param.pos = false; //default

  //define all arguments
  auto cli =
    (
     clipp::required("-a") & clipp::value("alpha", param.alpha).doc("path length in variation graph (e.g., 500)"),
     clipp::required("-d") & clipp::value("delta", param.delta).doc("differences allowed (e.g., 10)"),
     clipp::required("-chr") & clipp::value("id", param.chr).doc("chromosome id (e.g., 22 ), make it consistent with vcf file"),
     clipp::required("-vcf") & clipp::value("file1", param.vcffile).doc("compressed vcf file (something.vcf.gz)"),
     clipp::required("-fa") & clipp::value("file2", param.fasta_ref_file).doc("reference genome fasta file (something.fa)"),
     clipp::option("-prefix") & clipp::value("file3", param.prefix).doc("filename to optionally save input and output variants"),
    );

  if(!clipp::parse(argc, argv, cli))
  {
    //print help page
    clipp::operator<<(std::cout, clipp::make_man_page(cli, argv[0])) << std::endl;
    exit(1);
  }

  std::cout << "INFO, hged::parseandSave, alpha = " << param.alpha << std::endl;
  std::cout << "INFO, hged::parseandSave, delta = " << param.delta << std::endl;
  std::cout << "INFO, hged::parseandSave, chromosome id = " << param.chr << std::endl;
  std::cout << "INFO, hged::parseandSave, vcf file = " << param.vcffile << std::endl;
  std::cout << "INFO, hged::parseandSave, fasta file = " << param.fasta_ref_file << std::endl;
  if (param.prefix.length() > 0) std::cout << "INFO, hged::parseandSave, prefix = " << param.prefix << std::endl;

  if (! exists(param.vcffile))
  {
    std::cerr << "ERROR, hged::parseandSave, vcf file cannot be opened" << std::endl;
    exit(1);
  }
}

