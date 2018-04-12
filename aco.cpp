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
    ~ACO();
    int get_n_snps();
    void set_n_snps(int);
    int get_n_individuals();
    void set_n_individuals(int);
    double* get_pheremone_vals();
    void init_from_gd();
    void attach_gd_obj(GenomeData &obj);
    void init_pheremone(double);
    void evap_pheremone();
    void run(int);
};

ACO::ACO() {
  n_ants = 200;
  tournament_size = 100;
  initial_pheremone = 1;
  evap_rate = 0.99;
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

double* ACO::get_pheremone_vals() {
  return ACO::pheremone_vals;
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

void ACO::run(int max_iter) {
  for (int iter = 0; iter < max_iter; iter++) {
    cout << "\rIteration: " << iter;
    // Calculate sum of pheremone over all ant paths
    double pheremone_sum = 0;
    for (int i = 0; i < ACO::n_snps; i++) {
      pheremone_sum += ACO::pheremone_vals[i];
    }

    for (int ant = 0; ant < n_ants; ant++) {
      int snps[2] = {-1, -1};
      //cout << "getting random snps" << endl;
      //ACO::data->roulette_wheel_select(2, snps, ACO::pheremone_vals);
      ACO::data->tournament_select(2, snps, ACO::tournament_size, ACO::pheremone_vals);
      //cout << "SNPs: " << endl;
      //cout << "\t1 -> " << rand_snps[0] << endl;
      //cout << "\t2 -> " << rand_snps[1] << endl;

      int counts[32];
      for (int i = 0; i < 32; ++i) {
        counts[i] = 0;
      }
      ACO::data->count_pairwise_genotypes(snps[0], snps[1], counts);
      double crit = ACO::data->chi_sq_val_32(counts);
      if (crit > 60) {
        cout << endl << "SNPs " << snps[0] << " & " << snps[1] << " produced crit " << crit << endl;
      }
      //cout << "crit -> " << crit << endl;

      // Deposit Pheromone
      ACO::pheremone_vals[snps[0]] += crit;
      ACO::pheremone_vals[snps[1]] += crit;

      // Clamp pheromone to encourage exploration
      // if (ACO::pheremone_vals[rand_snps[0]] > 10) {
      //   cout << "clamping pheromone for snp " << rand_snps[0] << endl;
      // }

    }

    ACO::evap_pheremone();

  }

  int high_count = 0;
  float threshold = 1.0;
  for (int i = 0; i < ACO::n_snps; i++) {
    if (ACO::pheremone_vals[i] > threshold) {
      ++high_count;
    }
  }
  cout << endl << high_count << " SNPs above " << threshold << " threshold" << endl;
}

/*
This function should be fragmented up so that only the core algorithm is
implemented here, the fitness function should be called, not inline within
this method.
*/
// void ACO::run(int max_iter) {
//   for (int iter = 0; iter < max_iter; iter++) {
//     // Calculate sum of pheremone over all ant paths
//     double pheremone_sum = 0;
//     for (int i = 0; i < ACO::n_snps; i++) {
//       pheremone_sum += ACO::pheremone_vals[i];
//     }
//     cout << "Pheremone Sum: " << pheremone_sum << endl;
//
//     // Generate random number from uniform distribution to represent a SNP
//     std::random_device rd;
//     std::mt19937 rng(rd());
//     std::uniform_real_distribution<double> uni(0, pheremone_sum);
//
//     for (int ant = 0; ant < n_ants; ant++) {
//       double rand_doub = uni(rng);
//
//       for (int i = 0; i < ACO::n_snps; i++) {
//         if (rand_doub <= 0) {
//           int counts[4];
//           ACO::data->count_snp_alleles(i, counts);
//           double crit = ACO::data->chi_sq_val(counts) / 50;
//           // Deposit pheremone
//           ACO::pheremone_vals[i] += crit;
//           break;
//         } else {
//           rand_doub -= ACO::pheremone_vals[i];
//         }
//       }
//
//     }
//
//     ACO::evap_pheremone();
//   }
//   for (int i = 0; i < ACO::n_snps; i++) {
//     if (ACO::pheremone_vals[i] > 0.1) {
//       cout << "pheremone_vals[" << i << "] = " << pheremone_vals[i] << endl;
//     }
//   }
// }

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
  aco.run(20000);
}
