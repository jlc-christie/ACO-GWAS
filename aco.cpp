#include <iostream>
#include <random>
#include <boost/math/distributions/chi_squared.hpp>
#include "GenomeData.h"

using namespace std;

class ACO {
    int n_ants;
    int n_snps;
    int n_individuals;
    int tournament_size;
    double* pheremone_vals = nullptr;
    float initial_pheremone;
    float evap_rate;
    GenomeData* data = nullptr;

  public:
    ACO();
    ACO(int);
    ACO(int, int, int, float, float);
    ~ACO();
    int get_n_snps();
    void set_n_snps(int);
    int get_n_individuals();
    void set_n_individuals(int);
    void init_from_gd();
    void attach_gd_obj(GenomeData &obj);
    void init_pheremone(double);
    void evap_pheremone();
    void run(int);
};

ACO::ACO() {
  n_ants = 100;
  tournament_size = 50;
  initial_pheremone = 1;
  evap_rate = 0.99;
}

ACO::ACO(int n_snps) {
  n_ants = 100;
  ACO::n_snps = n_snps;
  tournament_size = 50;
  initial_pheremone = 1;
  evap_rate = 0.99;
}

ACO::ACO(int n_ants, int n_snps, int tournament_size,
         float initial_pheremone, float evap_rate) {
  ACO::n_ants = n_ants;
  ACO::n_snps = n_snps;
  ACO::tournament_size = tournament_size;
  ACO::initial_pheremone = initial_pheremone;
  ACO::evap_rate = evap_rate;
}

ACO::~ACO() {
  if (ACO::pheremone_vals != nullptr) {
    delete[] pheremone_vals;
  }
}

int ACO::get_n_snps() {
  return n_snps;
}

void ACO::set_n_snps(int n) {
  ACO::n_snps = n;
}

int ACO::get_n_individuals() {
  return n_individuals;
}

void ACO::set_n_individuals(int n) {
  ACO::n_individuals = n;
}

void ACO::init_from_gd() {
  if (ACO::data != nullptr) {
    ACO::set_n_snps(ACO::data->get_snp_count());
    ACO::set_n_individuals(ACO::data->get_individual_count());
  } else {
    cerr << "No GenomeData object initialised" << endl;
  }
}

void ACO::attach_gd_obj(GenomeData &obj) {
  ACO::data = &obj;
}

void ACO::init_pheremone(double d) {
  ACO::pheremone_vals = new double[ACO::n_snps];
  fill_n(ACO::pheremone_vals, ACO::n_snps, d);
}

void ACO::evap_pheremone() {
  for (int i = 0; i < ACO::n_snps; i++) {
    ACO::pheremone_vals[i] *= ACO::evap_rate;
  }
}

/*
This function should be fragmented up so that only the core algorithm is
implemented here, the fitness function should be called, not inline within
this method.
*/
void ACO::run(int max_iter) {
  for (int iter = 0; iter < max_iter; iter++) {
    // Calculate sum of pheremone over all ant paths
    double pheremone_sum = 0;
    for (int i = 0; i < ACO::n_snps; i++) {
      pheremone_sum += ACO::pheremone_vals[i];
    }
    cout << "Pheremone Sum: " << pheremone_sum << endl;

    // Generate random number from uniform distribution to represent a SNP
    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_real_distribution<double> uni(0, pheremone_sum);

    for (int ant = 0; ant < n_ants; ant++) {
      double rand_doub = uni(rng);

      for (int i = 0; i < ACO::n_snps; i++) {
        if (rand_doub <= 0) {
          int counts[4];
          ACO::data->count_snp_alleles(i, counts);
          double crit = ACO::data->chi_sq_val(counts) / 50;
          // Deposit pheremone
          ACO::pheremone_vals[i] += crit;
          break;
        } else {
          rand_doub -= ACO::pheremone_vals[i];
        }
      }

    }

    ACO::evap_pheremone();
  }
  for (int i = 0; i < ACO::n_snps; i++) {
    if (ACO::pheremone_vals[i] > 0.1) {
      cout << "pheremone_vals[" << i << "] = " << pheremone_vals[i] << endl; 
    }
  }
}

int main(void) {
  GenomeData gd = GenomeData();
  gd.init_individual_data("data/extend_directly_genotyped_binary.fam");
  gd.init_snp_data("data/extend_directly_genotyped_binary.bim");
  gd.init_binary_genotype_data("data/extend_directly_genotyped_binary.bed");
  gd.init_phenotype_file("data/phenotypes.txt");

  // int counts[4];
  // for (int i = 0; i < n_snps; i++) {
  //   //cout << "\rsnp: " << i;
  //   gd.count_snp_alleles(i, counts);
  //   // cout << endl;
  //   // cout << counts[0] << " " << counts[1] << endl;
  //   // cout << counts[2] << " " << counts[3] << endl;
  //   //float crit = gd.chi_sq_val(counts);
  //   // cout << endl << "chi-sq: " << crit << endl;
  //   // int dof = 1; // 2x2 contingency table, hence (2-1)*(2-1) dof
  //   // boost::math::chi_squared distribution(1);
  //   // float p = 1 - boost::math::cdf(distribution, crit);
  //   // cout << "p-val: " << p << endl;
  // }

  ACO aco = ACO();
  aco.attach_gd_obj(gd);
  aco.init_from_gd();
  cout << "SNP count: " << aco.get_n_snps() << endl;
  cout << "Individual count: " << aco.get_n_individuals() << endl;
  aco.init_pheremone(1.0);
  aco.run(1000);
}
