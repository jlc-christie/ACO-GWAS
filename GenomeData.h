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
    char* bin_genotype_data = nullptr;
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
    static int genotype_val(int, int);
    void count_snp_alleles(int snp_num, int allele_counts[4]);
    void count_pairwise_genotypes(int snp1, int snp2, int counts[32]);
    static double chi_sq_val(int counts[4]);
    static double chi_sq_val_32(int counts[32]);
    int get_individual_count();
    int get_snp_count();
    void roulette_wheel_select(int n, int snps[], double pheremone_vals[]);
    void tournament_select(int n, int snps[], int tourn_size, double pheremone_vals[]);
    static void process_byte(bitset<8> b, int n, int counts[8],
                             int pat_no, int pheno_mask[]);
    static void process_byte_pair(bitset<8> b1, bitset<8> b2, int n,
                                  int counts[18], int pat_no, int pheno_mask[]);
};
