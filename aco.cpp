#include <iostream>
#include <boost/math/distributions/chi_squared.hpp>
#include "GenomeData.h"

using namespace std;

class ACO {
    int n_ants;
    int n_snps;
    int tournament_size;
    //int pheremone_vals[];
    float initial_pheremone;
    float evap_rate;

  public:
    ACO(int);
    ACO(int, int, int, float, float);
    int getNSnps();
    void init_pheremone();
};

ACO::ACO(int n_snps) {
  n_ants = 100;
  ACO::n_snps = n_snps;
  tournament_size = 50;
  initial_pheremone = 1;
  evap_rate = 0.01;
}

ACO::ACO(int n_ants, int n_snps, int tournament_size,
         float initial_pheremone, float evap_rate) {
  ACO::n_ants = n_ants;
  ACO::n_snps = n_snps;
  ACO::tournament_size = tournament_size;
  ACO::initial_pheremone = initial_pheremone;
  ACO::evap_rate = evap_rate;
}

int ACO::getNSnps() {
  return n_snps;
}

int main(void) {
  GenomeData gd = GenomeData();
  gd.init_individual_data("data/extend_directly_genotyped_binary.fam");
  gd.init_snp_data("data/extend_directly_genotyped_binary.bim");
  gd.init_binary_genotype_data("data/extend_directly_genotyped_binary.bed");
  gd.init_phenotype_file("data/phenotypes.txt");

  int n_snps = gd.get_snp_count();

  int counts[4];
  for (int i = 0; i < n_snps; i++) {
    cout << "\rsnp: " << i;
    gd.count_snp_alleles(i, counts);
    cout << endl;
    cout << counts[0] << " " << counts[1] << endl;
    cout << counts[2] << " " << counts[3] << endl;
    float crit = gd.chi_sq_val(counts);
    cout << endl << "chi-sq: " << crit << endl;
    int dof = 1; // 2x2 contingency table, hence (2-1)*(2-1) dof
    boost::math::chi_squared distribution(1);
    float p = 1 - boost::math::cdf(distribution, crit);
    cout << "p-val: " << p << endl;
    int z;
    cin >> z;
  }

  ACO aco = ACO(n_snps);
  cout << aco.getNSnps() << endl;
}
