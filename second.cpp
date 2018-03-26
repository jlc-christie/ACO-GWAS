#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <bitset>
#include <sstream>

using namespace std;

void process_byte(bitset<8> b, int n, int counts[8], int pat_no, int pheno_mask[]);

class GenomeData {
    int individual_count;
    int snp_count;
    static const int MAX_INDIVIDUAL = 50000;
    int pheno_mask[MAX_INDIVIDUAL];
    map<int, int> id_to_pos;
    ifstream* bed_file;
    unsigned char* bin_genotype_data;
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
    int get_individual_count();
    int get_snp_count();
};

GenomeData::GenomeData() {
  individual_count = -1;
  snp_count = -1;
}

GenomeData::~GenomeData() {
  // Close bed file
  if (GenomeData::bed_file->is_open()) {
    GenomeData::bed_file->close();
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

  return true;
}

bool GenomeData::init_phenotype_file(string filename) {
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

void GenomeData::count_snp_alleles(int snp_num, int allele_counts[4]) {
  int counts[8] = {0, 0, 0, 0, 0, 0, 0, 0};

  // Seek to start of SNP
  int snp_offset = 3 + snp_num*GenomeData::total_bytes;
  GenomeData::bed_file->seekg(snp_offset);

  // Allocate memory for SNP data and then read data in
  char snp_buf[GenomeData::total_bytes];
  GenomeData::bed_file->read(snp_buf, GenomeData::total_bytes);

  // Can the GenomeData:: be removed for brevity?
  int complete_bytes = GenomeData::remainder_exists ? GenomeData::total_bytes - 1 : GenomeData::total_bytes;

  for (int i = 0; i < complete_bytes; i++) {
    bitset<8> b = snp_buf[i];
    process_byte(b, 8, counts, i*4, GenomeData::pheno_mask);
  }

  allele_counts[0] = 2*counts[0] + counts[1];
  allele_counts[1] = counts[1] + 3*counts[2];
  allele_counts[2] = 2*counts[4] + counts[5];
  allele_counts[3] = counts[5] + 2*counts[6];
}

int GenomeData::get_individual_count(void) {
  return GenomeData::individual_count;
}

int GenomeData::get_snp_count(void) {
  return GenomeData::snp_count;
}

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

int main(void) {
  GenomeData gd = GenomeData();
  gd.init_individual_data("data/extend_directly_genotyped_binary.fam");
  gd.init_snp_data("data/extend_directly_genotyped_binary.bim");
  gd.init_binary_genotype_data("data/extend_directly_genotyped_binary.bed");
  gd.init_phenotype_file("data/phenotypes.txt");

  int counts[4];
  int n_snps = gd.get_snp_count();
  for (int i = 0; i < n_snps; i++) {
    cout << "\rsnp: " << i;
    gd.count_snp_alleles(i, counts);
  }
  // cout << counts[0] << " " << counts[1] << endl;
  // cout << counts[2] << " " << counts[3] << endl;
}
