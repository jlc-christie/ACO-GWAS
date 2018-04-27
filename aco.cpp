#include <iostream>
#include <algorithm>
#include <iterator>
#include <random>
#include <string>
#include <vector>
#include <fstream>
#include <math.h>
#include <boost/math/distributions/chi_squared.hpp>
#include "GenomeData.h"

using namespace std;

class ACO {
    int n_ants;
    int n_snps;
    int n_individuals;
    int tournament_size;
    int total_evals = 0;
    double* pheromone_vals = nullptr;
    vector<double> best_fitnesses;
    vector<double> single_assoc_fits;
    int best_fit_snps[20];
    float initial_pheromone;
    float evap_rate;
    GenomeData* data = nullptr;

  public:
    ACO();
    ~ACO();
    int get_n_snps();
    void set_n_snps(int);
    int get_n_individuals();
    void set_n_individuals(int);
    double* get_pheromone_vals();
    void init_from_gd();
    void attach_gd_obj(GenomeData &obj);
    void init_pheromone(double);
    void evap_pheromone();
    void push_single_assoc(double ts);
    double next_generation_adjusted(int best_snps[2]);
    double next_generation(int best_snps[2]);
    int get_total_evals();
    void run(int);
    void run(int, double[]);
    void save_pheromone_matrix(string filename);
};

ACO::ACO() {
  n_ants = 200;
  tournament_size = 50;
  initial_pheromone = 1;
  evap_rate = 0.99;

  for (int i = 0; i < 10; ++i) {
    best_fitnesses.push_back(10.0);
  }
}

ACO::~ACO() {
  if (ACO::pheromone_vals != nullptr) {
    delete[] pheromone_vals;
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

double* ACO::get_pheromone_vals() {
  return ACO::pheromone_vals;
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

void ACO::init_pheromone(double d) {
  ACO::pheromone_vals = new double[ACO::n_snps];
  fill_n(ACO::pheromone_vals, ACO::n_snps, d);
}

void ACO::evap_pheromone() {
  for (int i = 0; i < ACO::n_snps; i++) {
    ACO::pheromone_vals[i] *= ACO::evap_rate;
  }
}

void ACO::push_single_assoc(double ts) {
  ACO::single_assoc_fits.push_back(ts);
}

double ACO::next_generation_adjusted(int best_snps[2]) {
  double best_fitness = 0.0;

  // Calculate the cumulative sum of pheromone for use with tournament selection
  double pheromone_sum = 0;
  for (int i = 0; i < ACO::n_snps; ++i) {
    pheromone_sum += ACO::pheromone_vals[i];
  }

  for (int ant = 0; ant < ACO::n_ants; ++ant) {
    int snps[2] = {-1, -1};
    ACO::data->tournament_select(2, snps, ACO::tournament_size, ACO::pheromone_vals);

    // Create contingency table and calulate critical value from it
    int counts[32];
    for (int i = 0; i < 32; ++i) { counts[i] = 0; }
    ACO::data->count_pairwise_genotypes(snps[0], snps[1], counts);
    ++ACO::total_evals;
    double crit = ACO::data->chi_sq_val_32(counts);

    double adjusted_crit = crit - ACO::single_assoc_fits[snps[0]] - ACO::single_assoc_fits[snps[1]];

    if (adjusted_crit > ACO::best_fitnesses.back()) {
      vector<double>::iterator lower, upper;
      lower = lower_bound(ACO::best_fitnesses.begin(), ACO::best_fitnesses.end(), adjusted_crit, greater<double>());

      ACO::best_fitnesses.insert(lower, adjusted_crit);
      ACO::best_fitnesses.pop_back();

      int pos = (lower - ACO::best_fitnesses.begin());
      ACO::best_fit_snps[pos*2] = snps[0];
      ACO::best_fit_snps[(pos*2)+1] = snps[1];


    }

    if (adjusted_crit >= best_fitness) {
      if (adjusted_crit > best_fitness) {
        best_fitness = adjusted_crit;
        for (int i = 0; i < 2; ++i) { best_snps[i] = snps[i]; }
      } else {
        // If theres a tie, do a coin-flip for the win
        static std::random_device rd;
        static std::mt19937 rng(rd());
        std::uniform_int_distribution<int> uni(0, 1);
        int choice = uni(rng);
        if (choice) {
          best_fitness = adjusted_crit;
          for (int i = 0; i < 2; ++i) { best_snps[i] = snps[i]; }
        }
      }
    }

    // Deposit Pheromone equally on SNPs
    for (int i = 0; i < 2; ++i) { ACO::pheromone_vals[snps[i]] += adjusted_crit; }
  }

  return best_fitness;
}

double ACO::next_generation(int best_snps[2]) {
  double best_fitness = 0.0;    // The fitness associated with this combination

  // Calculate the cumulative sum of pheromone for use with tournament selection
  double pheromone_sum = 0;
  for (int i = 0; i < ACO::n_snps; ++i) {
    pheromone_sum += ACO::pheromone_vals[i];
  }

  // Let each ant traverse the graph
  for (int ant = 0; ant < ACO::n_ants; ++ant) {
    // Select 2 SNPs using tournament selection
    int snps[2] = {-1, -1};
    ACO::data->tournament_select(2, snps, ACO::tournament_size, ACO::pheromone_vals);

    // Create contingency table and calulate critical value from it
    int counts[32];
    for (int i = 0; i < 32; ++i) { counts[i] = 0; }
    ACO::data->count_pairwise_genotypes(snps[0], snps[1], counts);
    ++ACO::total_evals;
    double crit = ACO::data->chi_sq_val_32(counts);

    if (crit > ACO::best_fitnesses.back()) {
      vector<double>::iterator lower, upper;
      lower = lower_bound(ACO::best_fitnesses.begin(), ACO::best_fitnesses.end(), crit, greater<double>());

      ACO::best_fitnesses.insert(lower, crit);
      ACO::best_fitnesses.pop_back();

      int pos = (lower - ACO::best_fitnesses.begin());
      ACO::best_fit_snps[pos*2] = snps[0];
      ACO::best_fit_snps[(pos*2)+1] = snps[1];


    }

    if (crit >= best_fitness) {
      if (crit > best_fitness) {
        best_fitness = crit;
        for (int i = 0; i < 2; ++i) { best_snps[i] = snps[i]; }
      } else {
        // If theres a tie, do a coin-flip for the win
        static std::random_device rd;
        static std::mt19937 rng(rd());
        std::uniform_int_distribution<int> uni(0, 1);
        int choice = uni(rng);
        if (choice) {
          best_fitness = crit;
          for (int i = 0; i < 2; ++i) { best_snps[i] = snps[i]; }
        }
      }
    }

    // Deposit Pheromone equally on SNPs
    for (int i = 0; i < 2; ++i) { ACO::pheromone_vals[snps[i]] += crit; }
  }

  return best_fitness;
}

int ACO::get_total_evals() {
  return ACO::total_evals;
}

void ACO::run(int max_iter) {
  // Stores best crit for each iteration
  double crit_progress[max_iter];

  for (int iter = 0; iter < max_iter; iter++) {
    int best_snps[2] = {-1, -1};
    double best_adj_crit = 0.0;
    double best_crit = 0.0;

    cout << "\rIteration: " << iter;

    // Calculate sum of pheromone over all ant paths
    double pheromone_sum = 0;
    for (int i = 0; i < ACO::n_snps; i++) {
      pheromone_sum += ACO::pheromone_vals[i];
    }

    for (int ant = 0; ant < n_ants; ant++) {
      // Select SNPs using a particular selection method
      int snps[2] = {-1, -1};
      //ACO::data->roulette_wheel_select(2, snps, ACO::pheromone_vals);
      ACO::data->tournament_select(2, snps, ACO::tournament_size, ACO::pheromone_vals);


      // Create contingency table and calculate critical value from it
      int counts[32];
      for (int i = 0; i < 32; ++i) {
        counts[i] = 0;
      }
      ACO::data->count_pairwise_genotypes(snps[0], snps[1], counts);
      double crit = ACO::data->chi_sq_val_32(counts);

      double ind_crits[2];
      int ind_counts[4];
      for (int i = 0; i < 2; ++i) {
        ind_counts[0] = 0;
        ind_counts[1] = 0;
        ind_counts[2] = 0;
        ind_counts[3] = 0;
        ACO::data->count_snp_alleles(snps[i], ind_counts);
        ind_crits[i] = ACO::data->chi_sq_val(ind_counts);
      }

      double adjusted_crit = crit - ind_crits[0] - ind_crits[1];

      // Keep track of best pair of SNPs for this iteration
      if (adjusted_crit > best_crit) {
        best_adj_crit = adjusted_crit;
        best_crit = crit;
        best_snps[0] = snps[0];
        best_snps[1] = snps[1];
      }

      if (adjusted_crit > 40) {
        cout << "best_crit = " << best_crit << " for " << best_snps[0] << " & " << best_snps[1] << " with " << ind_crits[0] << " " << ind_crits[1] << endl;
      }

      // Deposit Pheromone
      ACO::pheromone_vals[snps[0]] += adjusted_crit;
      ACO::pheromone_vals[snps[1]] += adjusted_crit;

    }
    // Add best crit for this iteration to array
    crit_progress[iter] = best_crit;

    ACO::evap_pheromone();
  }

  // Write progress for this iteration to file
  ofstream progress_file;
  ofstream prog_av;
  progress_file.open("results/progress_file_init1_tourn50_ants200_evap_099_1_adj.dat");
  prog_av.open("results/prog_av_init1_tourn50_ants200_evap_099_1_adj.dat");

  for (int i = 0; i < max_iter; ++i) {
    if (i % 100 == 0 && i != 0) {
      double tot = 0.0;
      for (int j = 0; j < 100; ++j) {
        tot += crit_progress[i - j];
      }
      prog_av << i << " " << tot / 100 << endl;
    }
    progress_file << i << " " << crit_progress[i] << endl;
  }
  progress_file.close();
  prog_av.close();

  int high_count = 0;
  float threshold = 1.0;
  for (int i = 0; i < ACO::n_snps; i++) {
    if (ACO::pheromone_vals[i] > threshold) {
      ++high_count;
    }
  }
  cout << endl << high_count << " SNPs above " << threshold << " threshold" << endl;
}

void ACO::run(int max_iter, double crit_progress[]) {
  for (int iter = 0; iter < max_iter; ++iter) {
    // Change tournament size based on progress through algorithm assuming
    // 1 million total comparisons before exit
    if (total_evals > 900000) {
      ACO::tournament_size = 100;
    } else if (total_evals > 500000) {
      ACO::tournament_size = 500;
    } else if (total_evals > 250000) {
      ACO::tournament_size = 1000;
    } else {
      ACO::tournament_size = 2000;
    }
    int best_snps[2];
    double best_crit = ACO::next_generation_adjusted(best_snps);
    crit_progress[iter] = best_crit;
    ACO::evap_pheromone();
  }

  // Output top 10 SNPs
  cout << "Top 10 SNP combinations: " << endl;
  for (int i = 0; i < 10; ++i) {
    cout << i << " -> " << ACO::best_fitnesses[i] << endl;
    cout << "\t" << ACO::best_fit_snps[i*2] << " (" << ACO::data->get_snp_id(ACO::best_fit_snps[i*2]) << ") CHISQ = " << ACO::single_assoc_fits[best_fit_snps[i*2]] << endl;
    cout << "\t" << ACO::best_fit_snps[(i*2)+1] << " (" << ACO::data->get_snp_id(ACO::best_fit_snps[(i*2)+1]) << ") CHISQ = " << ACO::single_assoc_fits[best_fit_snps[(i*2)+1]] << endl;
  }

}

/*
This function should be fragmented up so that only the core algorithm is
implemented here, the fitness function should be called, not inline within
this method.
*/
// void ACO::run(int max_iter) {
//   for (int iter = 0; iter < max_iter; iter++) {
//     // Calculate sum of pheromone over all ant paths
//     double pheromone_sum = 0;
//     for (int i = 0; i < ACO::n_snps; i++) {
//       pheromone_sum += ACO::pheromone_vals[i];
//     }
//     cout << "pheromone Sum: " << pheromone_sum << endl;
//
//     // Generate random number from uniform distribution to represent a SNP
//     std::random_device rd;
//     std::mt19937 rng(rd());
//     std::uniform_real_distribution<double> uni(0, pheromone_sum);
//
//     for (int ant = 0; ant < n_ants; ant++) {
//       double rand_doub = uni(rng);
//
//       for (int i = 0; i < ACO::n_snps; i++) {
//         if (rand_doub <= 0) {
//           int counts[4];
//           ACO::data->count_snp_alleles(i, counts);
//           double crit = ACO::data->chi_sq_val(counts) / 50;
//           // Deposit pheromone
//           ACO::pheromone_vals[i] += crit;
//           break;
//         } else {
//           rand_doub -= ACO::pheromone_vals[i];
//         }
//       }
//
//     }
//
//     ACO::evap_pheromone();
//   }
//   for (int i = 0; i < ACO::n_snps; i++) {
//     if (ACO::pheromone_vals[i] > 0.1) {
//       cout << "pheromone_vals[" << i << "] = " << pheromone_vals[i] << endl;
//     }
//   }
// }

void ACO::save_pheromone_matrix(string filename) {
  ofstream out_file;
  out_file.open(filename);
  int tot_snps = ACO::get_n_snps();
  int line_length = sqrt(tot_snps) + 1;

  for (int i = 0; i < line_length; ++i) {
    for (int j = 0; j < line_length; j++) {
      out_file << ACO::pheromone_vals[(i*line_length) + j] << " ";
    }
    out_file << endl;
  }
  out_file.close();
}

int main(void) {
  GenomeData gd = GenomeData();
  gd.init_individual_data("data/extend_directly_genotyped_binary.fam");
  gd.init_snp_data("data/extend_directly_genotyped_binary.bim");
  gd.init_binary_genotype_data("data/extend_directly_genotyped_binary.bed");
  gd.init_phenotype_file("data/phenotypes.txt");

  double* progresses[10];
  for (int i = 0; i < 10; ++i) {
    progresses[i] = new double[12000];
  }

  for (int run = 0; run < 10; ++run) {
    double* progress = progresses[run];

    ACO aco = ACO();
    aco.attach_gd_obj(gd);
    aco.init_from_gd();

    // Perform single association for adjusted crit
    int number_of_snps = gd.get_snp_count();
    for (int i = 0; i < number_of_snps; ++i) {
      int counts[4] = {0,0,0,0};
      gd.count_snp_alleles(i, counts);
      double ts = gd.chi_sq_val(counts);
      aco.push_single_assoc(ts);
    }

    aco.init_pheromone(25.0);
    aco.run(12000, progress);
    cout << "TOTAL EVALS -> " << aco.get_total_evals() << endl;

    ofstream progress_file;
    ofstream prog_av;
    progress_file.open("results/" + to_string(run) + ".dat");
    prog_av.open("results/" + to_string(run) + "_average.dat");
    for (int i = 0; i < 12000; ++i) {
      int n = 100;
      int window = 100;
      double sum = 0;
      // Moving window
      for (int j = 0; j < 100; ++j) {
        if (i - j < 0) {
          window -= 100 - j;
          break;
        }
        sum += progress[i-j];
      }
      prog_av << i << " " << sum / window << endl;

      // Keep tracking of only improving results
      if (progress[i] > 0 && progress[i-1] > progress[i]) {
        progress[i] = progress[i-1];
      }
      progress_file << i << " " << progress[i] << endl;
    }
    progress_file.close();
    prog_av.close();
  }

  ofstream prog_tot_av;
  prog_tot_av.open("results/total_av.dat");
  for (int i = 0; i < 12000; ++i) {
    int sum = 0;
    for (int j = 0; j < 10; ++j) {
      sum += progresses[j][i];
    }
    prog_tot_av << i << " " << sum << endl;
  }
  prog_tot_av.close();

  // Perform single association for adjusted crit
  // int number_of_snps = gd.get_snp_count();
  // for (int i = 0; i < number_of_snps; ++i) {
  //   int counts[4] = {0,0,0,0};
  //   gd.count_snp_alleles(i, counts);
  //   double ts = gd.chi_sq_val(counts);
  //   aco.push_single_assoc(ts);
  // }
  //
  // aco.init_pheromone(25.0);
  // aco.run(5200, progress);
  // cout << "TOTAL EVALS -> " << aco.get_total_evals() << endl;

  // int counts[4] = {0, 0, 0, 0};
  // gd.count_snp_alleles(54557, counts);
  // float crit = gd.chi_sq_val(counts);
  // cout << "crit 1 -> " << crit << endl;
  // counts[0] = 0;
  // counts[1] = 0;
  // counts[2] = 0;
  // counts[3] = 0;
  // gd.count_snp_alleles(298163, counts);
  // crit = gd.chi_sq_val(counts);
  // cout << "crit 2 -> " << crit << endl;
  //
  // int counts_2[32];
  // for (int i = 0; i < 32; ++i) {
  //   counts_2[i] = 0;
  // }
  // gd.count_pairwise_genotypes(54557, 298163, counts_2);
  // double crit_2 = gd.chi_sq_val_32(counts_2);
  //
  // cout << "combined crit -> " << crit_2 << endl;


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
  // double* progresses[10];
  // for (int i = 0; i < 10; ++i) {
  //   progresses[i] = new double[5200];
  // }
  //
  // for (int run = 0; run < 10; ++run) {
  //   double* progress = progresses[run];
  //   ACO aco = ACO();
  //   aco.attach_gd_obj(gd);
  //   aco.init_from_gd();
  //   cout << "SNP count: " << aco.get_n_snps() << endl;
  //   cout << "Individual count: " << aco.get_n_individuals() << endl;
  //   aco.init_pheromone(25.0);
  //   aco.run(5200, progress);
  //
  //   ofstream progress_file;
  //   ofstream prog_av;
  //   progress_file.open("results/200_25/newnew_" + to_string(run) + ".dat");
  //   prog_av.open("results/200_25/newnew_av_" + to_string(run) + ".dat");
  //   for (int i = 0; i < 5200; ++i) {
  //     if (i % 100 == 0 && i != 0) {
  //       double tot = 0.0;
  //       for (int j = 0; j < 100; ++j) {
  //         tot += progress[i - j];
  //       }
  //       prog_av << i << " " << tot / 100 << endl;
  //     }
  //     if (progress[i] > 0 && progress[i-1] > progress[i]) {
  //       progress[i] = progress[i-1];
  //     }
  //     progress_file << i << " " << progress[i] << endl;
  //   }
  //   progress_file.close();
  //   prog_av.close();
  // }
  //
  // ofstream prog_tot_av;
  // prog_tot_av.open("results/200_25/new_total_av.dat");
  // for (int i = 0; i < 5200; ++i) {
  //   int sum = 0;
  //   for (int j = 0; j < 10; ++j) {
  //     sum += progresses[j][i];
  //   }
  //   prog_tot_av << i << " " << sum << endl;
  // }
  // prog_tot_av.close();

  // 383256 or 258653 or 314033 or 419226 or 391440 (85) or 54557 (99) or 476876 (86) or 393965 (89) or 153425 (89)
  // 298163
  //  71.995

}
