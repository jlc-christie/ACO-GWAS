#include <iostream>
#include <fstream>
#include <bitset>
#include <iterator>
#include <map>
#include <string>
#include <vector>
#include <chrono>
#include <cmath>
#include <sstream>
#include <typeinfo>

#include "input.h"
#include "stats.h"

#define TEST 0

using namespace std;
using byte = unsigned char;

void process_case_controls(string filename, int individual_count, map<int,int> id_to_pos, int pheno_mask[]);

void process_byte(bitset<8> b, int n, int counts[8], int pat_no, int pheno_mask[]) {
  for (int i = 0; i < n; i += 2) {
    int pat = pat_no + (i/2);

    // Skip these two bits if patient has missing phenotype
    if (pheno_mask[pat] == 2) {
      continue;
    }

    int j = 0;
    if (pheno_mask[pat]) { j = 4; } // If patient is a control

    if (b[i] == 0) {
      if (b[i+1] == 0) {
        // Minor Homozygote
        counts[j] += 1;
      } else {
        // Heterozygote
        counts[j+1] += 1;
      }
    } else {
      if (b[i+1] == 1) {
        // Major Homozygote
        counts[j+2] += 1;
      } else {
        // Missing
        counts[j+3] += 1;
      }

    }
  }
}

int main() {
  vector<string> fam_range, bim_range;
  handle_input_files(&fam_range, &bim_range);

  // Set up map for individuals position within binary file to their id
  map<int, int> id_to_pos;
  int pos = 0;
  for (auto const &line : fam_range) {
    int id = stoi(line.substr(0, line.find(' ')));
    id_to_pos[id] = pos;
    pos++;
  }

  // Prints size of fam vector range (number of individuals in .fam file)
  unsigned int individual_count = fam_range.size();
  cout << individual_count << " individuals detected in .fam file" << endl;

  // Print size of bim vector range (number of SNPs in .bim file)
  unsigned int snp_count = bim_range.size();
  cout << snp_count << " SNPs detected in .bim file" << endl;

  // Open the .ped file to read in binary mode
  #if TEST
    ifstream bed_file ("data/test_data/binary.bed", ios::binary);
  #else
    ifstream bed_file ("data/extend_directly_genotyped_binary.bed", ios::binary);
  #endif
  if (!bed_file.is_open()) {
    cerr << "Unable to open bed file, exiting..." << endl;
    return -1;
  }

  int pheno_mask[individual_count];

  //map<int, string> phenotypes;
  #if TEST
    process_case_controls("data/test_data/phenotypes.txt", individual_count, id_to_pos, pheno_mask);
  #else
    process_case_controls("data/phenotypes.txt", individual_count, id_to_pos, pheno_mask);
  #endif

  // .bed format's 'magic' bytes
  // Third byte (0x01) means data is in SNP-major form, 0x00 means that the data
  // is in individual-major form and needs to be handled differently...
  bitset<8> magic_bytes[3] = {0x6c, 0x1b, 0x01};

  // Read first 3 bytes to ensure they match .bed format bytes
  for (unsigned int i = 0; i < 3; i++) {
    char buf[1];
    bed_file.read(buf, 1);
    bitset<8> b = buf[0];
    if (b != magic_bytes[i]) {
      cerr << "File not in bed format, exiting..." << endl;
      return -1;
    }
  }

  // SNP Loop
  for (int i = 0; i < snp_count; i++) {
    int counts[8] = {0, 0, 0, 0, 0, 0, 0, 0};

    // Hacky break in code
    // int xyz;
    // cin >> xyz;
    // cout << "\rSNP: " << i;

    int complete_bytes = individual_count / 4;
    int remaining_bits = 2 * (individual_count % 4);
    bool remainder_exists = remaining_bits > 0;

    for (int j = 0; j < complete_bytes; j++) {
      char buf[1];
      bed_file.read(buf, 1);
      bitset<8> b = buf[0];
      process_byte(b, 8, counts, j*4, pheno_mask);
    }

    if (remainder_exists) {
      int pat_no = individual_count - (remaining_bits / 2);
      char buf[1];
      bed_file.read(buf, 1);
      bitset<8> b = buf[0];
      process_byte(b, remaining_bits, counts, pat_no, pheno_mask);
    }

    // cout << endl;
    // cout << counts[0] << " " << counts[1] << " " << counts[2] << " " << counts[3] << endl;
    // cout << counts[4] << " " << counts[5] << " " << counts[6] << " " << counts[7] << endl;
    // //
    // cout << endl;
    // cout << 2*counts[0] + counts[1] << " " << counts[1] + 2*counts[2] << endl;
    // cout << 2*counts[4] + counts[5] << " " << counts[5] + 2*counts[6] << endl;
    //
    // Calculate Chi-Squared Value
    double chi_squared_val = chi_sq(counts);
    if (chi_squared_val > 20.0) {
      cout << "SNP " << i << " has CHISQ val of " << chi_squared_val << endl;
    }
    // cout << endl << chi_squared_val << endl;

  }

}

void process_case_controls(string filename, int individual_count, map<int,int> id_to_pos, int pheno_mask[]) {
  // Open phenotypes file to find case/controls
  ifstream phen_file (filename);
  if (!phen_file.is_open()) {
    cerr << "Unable to open phenotype file, exiting..." << endl;
    return;
  }

  cout << "Analysing phenotype file to find case and controls..." << endl;

  int controls = 0, cases = 0, missing = 0, individual_n = 0;

  string buf;
  getline(phen_file, buf, '\n');  // Skip header

  while (getline(phen_file, buf, '\n')) {
    string word_buf;
    stringstream ss(buf);
    vector<string> tokens;

    while (ss >> word_buf) {
      tokens.push_back(word_buf);
    }

    //tokens[4] = T2D presence, 1=no, 2=yes
    int phen_val;
    int pos;
    auto iter = id_to_pos.find(stoi(tokens[0]));
    if (iter != id_to_pos.end()) {
      pos = iter->second;
      //cout << pos << endl;
    }

    //phenotypes[individual_n] = tokens[4];

    try {
      phen_val = stoi(tokens[4]);
    } catch (exception& e) {
      pheno_mask[pos] = 2;
      missing++;
      individual_n++;
      continue;
    }

    if (phen_val == 2) {
      pheno_mask[pos] = 0;
      cases++;
    } else if (phen_val == 1) {
      pheno_mask[pos] = 1;
      controls++;
    }
    individual_n++;
  }

  cout << "Finished analysing phenotype file, results:" << endl;
  cout << "\tCases: " << cases << endl;
  cout << "\tControls: " << controls << endl;
  cout << "\tMissing: " << missing << endl;
}
