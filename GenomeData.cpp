#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <random>
#include <fstream>
#include <bitset>
#include <sstream>
#include <cmath>

#include <boost/math/distributions/chi_squared.hpp>

#include "GenomeData.h"

using namespace std;

GenomeData::GenomeData() {
  individual_count = -1;
  snp_count = -1;
}

GenomeData::~GenomeData() {
  // Close bed file
  if (GenomeData::bed_file->is_open()) {
    GenomeData::bed_file->close();
  }

  if (GenomeData::bin_genotype_data != NULL) {
    delete[] GenomeData::bin_genotype_data;
  }
}

bool GenomeData::init_individual_data(string filename) {
  /*
  In future method should get phenotype status from here ideally, rather than
  from a seperate phenotype file. Would be a lot easier at least. This means
  pheno_mask would be set up in here instead.
  */
  ifstream fam_file (filename);

  if (!fam_file.is_open()) {
    cerr << "Unable to open individual data (usually .fam) file" << endl;
    return false;
  }

  // Read data in to vector
  string buf;
  while (getline(fam_file, buf, '\n')) {
    GenomeData::individual_data.push_back(buf);
  }
  fam_file.close();

  GenomeData::individual_count = GenomeData::individual_data.size();

  cout << GenomeData::individual_count << " individuals detected in .fam file" << endl;

  int pos = 0;
  for (auto const &line : GenomeData::individual_data) {
    int id = stoi(line.substr(0, line.find(' ')));
    GenomeData::id_to_pos[id] = pos;
    pos++;
  }

  return true;
}

bool GenomeData::init_snp_data(string filename) {
  ifstream bim_file (filename);

  if (!bim_file.is_open()) {
    cerr << "Unable to open snp data (usually .bim) file" << endl;
    return false;
  }

  // Read data in to vector
  string buf;
  while (getline(bim_file, buf, '\n')) {
    GenomeData::snp_data.push_back(buf);
  }
  bim_file.close();

  GenomeData::snp_count = GenomeData::snp_data.size();

  cout << GenomeData::snp_count << " SNPs detected in .bed file" << endl;
  return true;
}

bool GenomeData::init_binary_genotype_data(string filename) {
  //ifstream temp_file (filename, ios::binary);
  GenomeData::bed_file = new ifstream(filename, ios::binary);
  if (!bed_file->is_open()) {
    cerr << "Unable to open binary data file (usually .bed)" << endl;
    return false;
  }

  // Find the size of the file in bytes according to seekg.
  size_t beginning, end;
  beginning = bed_file->tellg();
  bed_file->seekg(0, ios::end);
  end = bed_file->tellg();

  // Calculate total bytes in each SNP
  int complete_bytes = GenomeData::individual_count / 4;
  int remaining_bits = GenomeData::individual_count % 4;
  bool remainder_exists = remaining_bits != 0;
  int total_bytes = remainder_exists ? complete_bytes + 1 : complete_bytes;

  // Seekg can sometimes produce incorrect end position if file is buffered
  // (which it shouldn't be), so we check the number of bytes seekg finds to be
  // the same as we'd expect given the bim and fam files.
  if (3 + total_bytes*GenomeData::snp_count != end - beginning) {
    cerr << "Discrepancy in the number of bytes expected and found in" << endl;
    cerr << " the binary file, failed to initialize genotype data" << endl;
  }

  // We can now be somewhat sure the values obtained are correct
  GenomeData::remainder_exists = remainder_exists;
  GenomeData::total_bytes = total_bytes;

  // .bed format's 'magic' bytes
  // Third byte (0x01) means data is in SNP-major form, 0x00 means that the data
  // is in individual-major form and needs to be handled differently...
  bitset<8> magic_bytes[3] = {0x6c, 0x1b, 0x01};

  // Read first 3 bytes to ensure they match .bed format bytes
  bed_file->seekg(0);
  for (int i = 0; i < 3; i++) {
    char buf[1];
    GenomeData::bed_file->read(buf, 1);
    bitset<8> b = buf[0];
    if (b != magic_bytes[i]) {
      cerr << "File not in bed format" << endl;
      return false;
    }
  }
  bed_file->seekg(0);

  // If there is enough ram, load whole binary file in to ram
  // Rather not allocate this amount of memory in the heap, can it be done on
  //  the stack somehow?
  int total_file_bytes = 3 + total_bytes * GenomeData::snp_count;
  try {
    GenomeData::bin_genotype_data = new char[total_file_bytes];
    GenomeData::bed_file->read(GenomeData::bin_genotype_data, total_file_bytes);
  } catch (exception &e) {
    cerr << "Couldn't allocate memory for whole binary data, program will ";
    cerr << "continue to read data SNP-by-SNP instead" << endl;
  }

  GenomeData::bed_file->clear();
  GenomeData::bed_file->seekg(0);

  return true;
}

bool GenomeData::init_phenotype_file(string filename) {
  /*
  This is *very* file specific for the data I have, and init_individual_data
  should init pheno_mask if phenotype values are present within the fam file.
  */
  ifstream phen_file (filename);
  if (!phen_file.is_open()) {
    cerr << "Unable to open phenotype file" << endl;
    return false;
  }

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
    auto iter = GenomeData::id_to_pos.find(stoi(tokens[0]));
    if (iter != GenomeData::id_to_pos.end()) {
      pos = iter->second;
      //cout << pos << endl;
    }

    //phenotypes[individual_n] = tokens[4];
    try {
      phen_val = stoi(tokens[4]);
    } catch (exception& e) {
      GenomeData::pheno_mask[pos] = 2;
      missing++;
      individual_n++;
      continue;
    }

    if (phen_val == 2) {
      GenomeData::pheno_mask[pos] = 0;
      cases++;
    } else if (phen_val == 1) {
      GenomeData::pheno_mask[pos] = 1;
      controls++;
    }

    individual_n++;
  }

  cout << "Finished analysing phenotype file, results:" << endl;
  cout << "\tCases: " << cases << endl;
  cout << "\tControls: " << controls << endl;
  cout << "\tMissing: " << missing << endl;

  return true;
}

void GenomeData::count_pairwise_genotypes(int snp1, int snp2, int counts[32]) {
  //cout << "in pairwise func with snp1 -> " << snp1 << " and snp2 -> " << snp2 << endl;
  int snp1_offset = 3 + snp1*total_bytes;
  int snp2_offset = 3 + snp2*total_bytes;
  int complete_bytes = remainder_exists ? total_bytes - 1 : total_bytes;

  for (int i = 0; i < complete_bytes; i++) {
    // cout << "Inside for loop with offsets: " << endl;
    // cout << "\t1 -> " << snp1_offset << endl;
    // cout << "\t2 -> " << snp2_offset << endl;

    bitset<8> b1 = bin_genotype_data[snp1_offset + i];
    bitset<8> b2 = bin_genotype_data[snp2_offset + i];

    // cout << "i -> " << i << endl;
    // cout << "b1 -> " << b1 << endl;
    // cout << "b2 -> " << b2 << endl;

    GenomeData::process_byte_pair(b1, b2, 8, counts, i*4, GenomeData::pheno_mask);
  }
}

void GenomeData::count_snp_alleles(int snp_num, int allele_counts[4]) {
  /*
  TO DO:
  this function does not account for SNPs where there are remaining bits
  */
  int counts[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  int snp_offset = 3 + snp_num*GenomeData::total_bytes;
  int complete_bytes = GenomeData::remainder_exists ? GenomeData::total_bytes - 1 : GenomeData::total_bytes;

  // This allocation will not be needed if binary data did fit in ram, but very
  // small and fast operation so it's worth it to keep code clean and readable
  char snp_buf[GenomeData::total_bytes];

  // If binary data couldn't be loaded in to RAM
  if (GenomeData::bin_genotype_data == NULL) {
    // Seek to start of SNP
    GenomeData::bed_file->seekg(snp_offset);

    for (int i = 0; i < complete_bytes; i++) {
      bitset<8> b = snp_buf[i];
      GenomeData::process_byte(b, 8, counts, i*4, GenomeData::pheno_mask);
    }
  } else { // If binary data IS in ram
    for (int i = 0; i < complete_bytes; i++) {
      bitset<8> b = GenomeData::bin_genotype_data[snp_offset + i];
      GenomeData::process_byte(b, 8, counts, i*4, GenomeData::pheno_mask);
    }
  }

  allele_counts[0] = 2*counts[0] + counts[1];
  allele_counts[1] = counts[1] + 2*counts[2];
  allele_counts[2] = 2*counts[4] + counts[5];
  allele_counts[3] = counts[5] + 2*counts[6];

}

double GenomeData::chi_sq_val(int counts[4]) {
  // Create a new array of doubles so floating point calculations can be done
  double d_counts[4];
  for (int i = 0; i < 4; i++) {
    d_counts[i] = double(counts[i]);
  }

  // Sum up the rows and columns of contingency table
  double total_a1, total_a2, total_cases, total_controls, total;
  total_a1 = d_counts[0] + d_counts[2];
  total_a2 = d_counts[1] + d_counts[3];
  total_cases = d_counts[0] + d_counts[1];
  total_controls = d_counts[2] + d_counts[3];

  total = total_cases + total_controls;

  // Calculate expected values, follows same order as allele counts array
  double expected[4];
  expected[0] = (total_cases * total_a1) / total;
  expected[1] = (total_cases * total_a2) / total;
  expected[2] = (total_controls * total_a1) / total;
  expected[3] = (total_controls * total_a2) / total;

  double res = 0.0;
  res += pow(d_counts[0] - expected[0], 2) / expected[0];
  res += pow(d_counts[1] - expected[1], 2) / expected[1];
  res += pow(d_counts[2] - expected[2], 2) / expected[2];
  res += pow(d_counts[3] - expected[3], 2) / expected[3];

  return res;
}

double GenomeData::chi_sq_val_32(int counts[32]) {
  // cout << "counts:" << endl;
  // cout << counts[0] << " " << counts[1] << " " << counts[2] << " " << counts[3];
  // cout << " " << counts[4] << " " << counts[5] << " " << counts[6] << " " << counts[7] << endl;
  // cout << counts[8] << " " << counts[9] << " " << counts[10] << " " << counts[11];
  // cout << " " << counts[12] << " " << counts[13] << " " << counts[14] << " " << counts[15] << endl;
  // cout << counts[16] << " " << counts[17] << " " << counts[18] << " " << counts[19];
  // cout << " " << counts[20] << " " << counts[21] << " " << counts[22] << " " << counts[23] << endl;
  // cout << counts[24] << " " << counts[25] << " " << counts[26] << " " << counts[27];
  // cout << " " << counts[28] << " " << counts[29] << " " << counts[30] << " " << counts[31] << endl;

  // Remove uneeded missing values and cast needed values to doubles
  double d_counts[18];
  int counter = 0;
  for (int i = 0; i < 32; i++) {
    if ((i+1) % 4 == 0 || ((i > 11) && (i < 16)) || i > 27) {
      // One of the genotype values is missing
      continue;
    } else {
      d_counts[counter] = double(counts[i]);
      counter++;
    }
  }

  // cout << "d_counts:" << endl;
  // cout << "\t" << d_counts[0] << " " << d_counts[1] << " " << d_counts[2] << " " << d_counts[3] << " " << d_counts[4] << " " << d_counts[5] << " " << d_counts[6] << " " << d_counts[7] << " " << d_counts[8] << endl;
  // cout << "\t" << d_counts[9] << " " << d_counts[10] << " " << d_counts[11] << " " << d_counts[12] << " " << d_counts[13] << " " << d_counts[14] << " " << d_counts[15] << " " << d_counts[16] << " " << d_counts[17] << endl;

  double col_totals[9];
  double row_totals[2];
  double total;

  // Calculate column totals
  for (int i = 0; i < 9; i++) {
    col_totals[i] = d_counts[i] + d_counts[9+i];
  }

  // Calculate row totals
  for (int i = 0; i < 2; i++) {
    row_totals[i] = 0;
    for (int j = 0; j < 9; j++) {
      row_totals[i] += d_counts[(i*9)+j];
    }
  }

  // Calculate total
  total = row_totals[0] + row_totals[1];

  double expected[18];
  for (int i = 0; i < 18; ++i) {
    expected[i] = 0;
  }

  // Fill expected value table
  for (int i = 0; i < 18; i++) {
    int col = i % 9;
    int row = (i < 9) ? 0 : 1;

    expected[i] = (col_totals[col] * row_totals[row]) / total;
  }

  // for (int i = 0; i < 18; ++i) {
  //   if (i % 9 == 0) {
  //     cout << endl;
  //   }
  //   cout << expected[i] << " ";
  // }
  // cout << endl;

  double res = 0.0;
  for (int i = 0; i < 18; i++) {
    if (expected[i] == 0) {
      continue;
    }
    res += pow(d_counts[i] - expected[i], 2) / expected[i];
  }

  return res;
}

int GenomeData::get_individual_count(void) {
  return GenomeData::individual_count;
}

int GenomeData::get_snp_count(void) {
  return GenomeData::snp_count;
}

void GenomeData::roulette_wheel_select(int n, int snps[], double pheremone_vals[]) {
  // Calculate sum of pheremone over all ant paths
  double pheremone_sum = 0;
  for (int i = 0; i < GenomeData::snp_count; i++) {
    pheremone_sum += pheremone_vals[i];
  }

  //Generate random number from uniform distribution to represent a SNP
  static std::random_device rd;
  static std::mt19937 rng(rd());
  std::uniform_real_distribution<double> uni(0, pheremone_sum);

  for (int i = 0; i < n; i++) {
    int snp = -1;
    bool not_unique = true;

    while (not_unique) {
      double snp_doub = uni(rng);

      for (int j = 0; j < GenomeData::snp_count; j++) {
        if (snp_doub <= 0) {
          snp = j;
          not_unique = false;

          // Check it's not a duplicate
          for (int k = 0; k < i; k++) {
            if (snps[k] == snp) {
              not_unique = true;
              int z;
              cin >> z;
              break;
            }
          }

          if (!not_unique) {
            snps[i] = snp;
          }

          break;
        } else {
          snp_doub -= pheremone_vals[j];
        }
      }
    }
  }
}

void GenomeData::tournament_select(int n, int snps[], int tourn_size, double pheremone_vals[]) {
  static std::random_device rd;
  static std::mt19937 rng(rd());
  std::uniform_int_distribution<int> uni(0, GenomeData::snp_count);

  for (int i = 0; i < n; ++i) {
    double max_pheromone = 0;
    int max_pher_snp = 0;
    for (int j = 0; j < tourn_size; ++j) {
      int snp = uni(rng);
      if (pheremone_vals[snp] > max_pheromone) {
        max_pheromone = pheremone_vals[snp];
        max_pher_snp = snp;
      }
    }

    snps[i] = max_pher_snp;
  }

}

void GenomeData::process_byte(bitset<8> b, int n, int counts[8],
                              int pat_no, int pheno_mask[]) {
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

int GenomeData::genotype_val(int first_bit, int second_bit) {
  if (first_bit == 0) {
    if (second_bit == 0) {
      // minor homo
      return 0;
    } else {
      // hetero
      return 1;
    }
  } else {
    if (second_bit == 1) {
      // major homo
      return 2;
    } else {
      return 3;
    }
  }
}

void GenomeData::process_byte_pair(bitset<8> b1, bitset<8> b2, int n, int counts[32],
                       int pat_no, int pheno_mask[]) {

  for (int i = 0; i < n; i += 2) {
    int pat = pat_no + (i/2);

    // Skip these two bits if patient has missing phenotype
    if (pheno_mask[pat] == 2) {
      continue;
    }

    int j = 0;
    if (pheno_mask[pat]) { j = 16; } // If patient is a control

    j += 4 * genotype_val(b1[i], b1[i+1]);
    j += genotype_val(b2[i], b2[i+1]);

    counts[j] += 1;
  }
}
//
// int main(void) {
//   GenomeData gd = GenomeData();
//   gd.init_individual_data("data/extend_directly_genotyped_binary.fam");
//   gd.init_snp_data("data/extend_directly_genotyped_binary.bim");
//   gd.init_binary_genotype_data("data/extend_directly_genotyped_binary.bed");
//   gd.init_phenotype_file("data/phenotypes.txt");
//
//   int counts[4];
//   int n_snps = gd.get_snp_count();
//   for (int i = 0; i < n_snps; i++) {
//     cout << "\rsnp: " << i;
//     gd.count_snp_alleles(i, counts);
//     cout << endl;
//     cout << counts[0] << " " << counts[1] << endl;
//     cout << counts[2] << " " << counts[3] << endl;
//     float crit = gd.chi_sq_val(counts);
//     cout << endl << "chi-sq: " << crit << endl;
//     int dof = 1; // 2x2 contingency table, hence (2-1)*(2-1) dof
//     boost::math::chi_squared distribution(1);
//     float p = boost::math::cdf(distribution, crit);
//     cout << "p-val: " << p << endl;
//     int z;
//     cin >> z;
//   }
//   // cout << counts[0] << " " << counts[1] << endl;
//   // cout << counts[2] << " " << counts[3] << endl;
// }
