#include <iostream>

using namespace std;

class ACO {
    int n_ants;
    int n_snps;
    int tournament_size;
    float initial_pheremone;
    float evap_rate;

  public:
    ACO(int);
    ACO(int, int, int, float, float);
    int getNSnps();
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

ACO::

int ACO::getNSnps() {
  return n_snps;
}

int main(void) {
  ACO aco = ACO(500000);
  cout << aco.getNSnps() << endl;
}
