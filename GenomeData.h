#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <bitset>
#include <sstream>

using namespace std;

class GenomeData {
    int individual_count;
    int snp_count;
    static const int MAX_INDIVIDUAL = 50000;
    int pheno_mask[MAX_INDIVIDUAL];
    map<int, int> id_to_pos;
    ifstream* bed_file;
    char* bin_genotype_data = NULL;
    bool remainder_exists;
    int total_bytes;
    vector<string> individual_data;
    vector<string> snp_data;
  public:
    GenomeData();
    ~GenomeData();
    bool init_individual_data(string filename);
    bool init_snp_data(string filename);
    bool init_binary_genotype_data(string filename);
    bool init_phenotype_file(string filename);
    void count_snp_alleles(int snp_num, int allele_counts[4]);
    static double chi_sq_val(int counts[4]);
    int get_individual_count();
    int get_snp_count();
    static void process_byte(bitset<8> b, int n, int counts[8],
                             int pat_no, int pheno_mask[]);
};
