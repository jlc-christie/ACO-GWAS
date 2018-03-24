#include <iostream>
#include <fstream>
#include <bitset>
#include <iterator>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <vector>
#include <chrono>
#include <cmath>
#include <sstream>
#include <typeinfo>

using namespace std;
using byte = unsigned char;

int main() {
  // Open the .fam file and read in to vector
  ifstream fam_file ("data/extend_directly_genotyped_binary.fam");
  if (!fam_file.is_open()) {
    cerr << "Unable to open fam file, exiting..." << endl;
    return -1;
  }
  vector<string> fam_range;
  string buf;
  while (getline(fam_file, buf, '\n')) {
    fam_range.push_back(buf);
  }
  fam_file.close();

  // Print size of fam vector range (number of individuals in .fam file)
  unsigned int individual_count = fam_range.size();
  cout << individual_count << " individuals detected in .fam file" << endl;

  // Open the .bim file and read in to vector
  ifstream bim_file ("data/extend_directly_genotyped_binary.bim");
  if (!bim_file.is_open()) {
    cerr << "Unable to open bim file, exiting..." << endl;
    return -1;
  }
  vector<string> bim_range;
  while (getline(bim_file, buf, '\n')) {
    bim_range.push_back(buf);
  }
  bim_file.close();

  // Print size of bim vector range (number of SNPs in .bim file)
  unsigned int snp_count = bim_range.size();
  cout << snp_count << " SNPs detected in .bim file" << endl;

  // Open the .ped file to read in binary mode
  ifstream bed_file ("data/extend_directly_genotyped_binary.bed", ios::binary);
  if (!bed_file.is_open()) {
    cerr << "Unable to open bed file, exiting..." << endl;
    return -1;
  }

  // Initialise input stream iterator to read byte-by-byte
  istream_iterator<byte> btt (bed_file);

  // .bed format's 'magic' bytes
  // Third byte (0x01) means data is in SNP-major form, 0x00 means that the data
  // is in individual-major form and needs to be handled differently...
  bitset<8> magic_bytes[3] = {0x6c, 0x1b, 0x01};

  // Read first 3 bytes to ensure they match .bed format magic bytes
  for (unsigned int i = 0; i < 3; i++) {
    bitset<8> b = (byte) *btt;
    if (b != magic_bytes[i]) {
      cerr << "File not in bed format, exiting..." << endl;
      return -1;
    }
    btt++;
  }

  // This is slow, is it needed?
  // cout << "Reading genotype values in to vector..." << endl;
  // vector<byte> bed_range = {btt, istream_iterator<byte>()};
  //
  // cout << "Read bits for genotypes..." << endl;
  // int foo;
  // cin >> foo;

  // Read bits for genotypes
  // Format for remaining bytes:
  // 2 bits per genotype, 4 genotypes per byte, encoded right-to-left
  // bits are also backwards.
  // e.g. 01101100
  //            AB  00 -- homozygote (first)
  //          CD    11 -- other homozygote (second)
  //        EF      01 -- heterozygote (third)
  //      GH        10 -- missing genotype (fourth)
  // Note: must skip to start of new byte when at end of SNP (or individual if
  // in individual-major form)

  // for (bitset<8> b : bed_range) {
  //   cout << hex << b.to_string() << endl;
  //   cout << "First genotype values: " << b[0] << b[1] << endl;
  //   cout << "Second genotype values: " << b[2] << b[3] << endl;
  //   cout << "Third genotype values: " << b[4] << b[5] << endl;
  //   cout << "Fourth genotype values: " << b[6] << b[7] << endl;
  //   cin >> foo;
  // }

  // unsigned int major_homozygote_case = 0;
  // unsigned int heterozygote_case = 0;
  // unsigned int minor_homozygote_case = 0;
  //
  // unsigned int major_homozygote_control = 0;
  // unsigned int heterozygote_control = 0;
  // unsigned int minor_homozygote_control = 0;



    /*
    Store ID's for cases (as they'll be less of them) in a unordered_map
    */
    unordered_set<unsigned int> case_ids;

    // Open phenotypes file to find case/controls
    ifstream phen_file ("data/phenotypes.txt");
    if (!phen_file.is_open()) {
      cerr << "Unable to open phenotype file, exiting..." << endl;
      return -1;
    }

    unsigned int missing = 0;
    unsigned int controls = 0;
    getline(phen_file, buf, '\n');  // Skip header
    while (getline(phen_file, buf, '\n')) {
      string word_buf;
      stringstream ss(buf);
      vector<string> tokens;

      while (ss >> word_buf) {
        tokens.push_back(word_buf);
      }
      // tokens[4] = T2D presence, 1=no, 2=yes
      unsigned int t2d_val;
      try {
          t2d_val = stoi(tokens[4]);
      } catch (exception& e) {
        // Should perform more detailed catch for invalid argument
        // and remove from both cases and controls, adds to controls atm
        missing++;
      }
      if (t2d_val == 2) {
        case_ids.insert(stoi(tokens[0]));
        cout << t2d_val;
      } else if (t2d_val == 1) {
        controls++;
        cout << t2d_val;
      }
    }

    unsigned int major_allelle_case = 0;
    unsigned int minor_allelle_case = 0;

    unsigned int major_allelle_control = 0;
    unsigned int minor_allelle_control = 0;

    unsigned int complete_bytes = individual_count / 8; // Remainder truncated
    unsigned int individual_counter = 0;
    for (unsigned int i = 0; i < complete_bytes; i++) {
      cout << "Iteration " << i << "-->";
      cout << major_allelle_case << " " << minor_allelle_case << " ";
      cout << major_allelle_control << " " << minor_allelle_control << endl;
      bitset<8> b = (byte) *btt;
      for (unsigned int j = 8; j > 0; j -= 2) {
        string individual_id_str;
        stringstream ss(fam_range[individual_counter]);
        ss >> individual_id_str;
        unsigned int individual_id_int = stoi(individual_id_str);

        if (case_ids.find(individual_id_int) == case_ids.end()) {
          if (b[i] == 0) {
            major_allelle_control++;
          } else {
            minor_allelle_control++;
          }

          if (b[i-1] == 0) {
            major_allelle_control++;
          } else {
            minor_allelle_control++;
          }
        } else {
          if (b[i] == 0) {
            major_allelle_case++;
          } else {
            minor_allelle_case++;
          }

          if (b[i-1] == 0) {
            major_allelle_case++;
          } else {
            minor_allelle_case++;
          }
        }
        individual_counter++;
      }
      btt++;
    }

    // Handle remaining data in last byte
    unsigned int remaining_bits = individual_count % 8;
    bitset<8> b = (byte) *btt;
    for (unsigned int i = 8; i > remaining_bits; i -= 2) {
      string individual_id_str;
      stringstream ss(fam_range[individual_counter]);
      ss >> individual_id_str;
      unsigned int individual_id_int = stoi(individual_id_str);

      if (case_ids.find(individual_id_int) == case_ids.end()) {
        if (b[i] == 0) {
          major_allelle_control++;
        } else {
          minor_allelle_control++;
        }

        if (b[i-1] == 0) {
          major_allelle_control++;
        } else {
          minor_allelle_control++;
        }
      } else {
        if (b[i] == 0) {
          major_allelle_case++;
        } else {
          minor_allelle_case++;
        }

        if (b[i-1] == 0) {
          major_allelle_case++;
        } else {
          minor_allelle_case++;
        }
      }
      individual_counter++;
    }
    btt++;
    cout << endl << "Minor Allelle Case: " << minor_allelle_case << endl;
    cout << "Major Allelle Case: " << major_allelle_case << endl;
    cout << "Minor Allelle Control: " << minor_allelle_control << endl;
    cout << "Major Allelle Control: " << major_allelle_control << endl;

    // Doesn't account for missing bits
    float minor_e = (float(minor_allelle_control) / float(controls)) * float(case_ids.size());
    float minor_o = float(minor_allelle_case);

    float major_e = (float(major_allelle_control) / float(controls)) * float(case_ids.size()); //0.87697
    float major_o = float(major_allelle_case); //0.914326

    cout << "major_e = " << major_e << endl;
    cout << "major_o = " << major_o << endl;

    float chi_sq_minor = pow(minor_o - minor_e, 2) / minor_e;
    float chi_sq_major = pow(major_o - major_e, 2) / major_e;
    cout << "Minor Allelle Chi-Squared = " << chi_sq_minor << endl;
    cout << "Major Allelle Chi-Squared = " << chi_sq_major << endl;

    cout << "Cases: " << case_ids.size() << endl;
    cout << "Missing: " << missing << endl;
    cout << "Controls: " << controls << endl;

    // for (string line : fam_range) {
    //
    // }

  bed_file.close();











  // string header_string;
  // unordered_map<string, unordered_map<string, int>> patient_data_container;
  //
  // string line_buff;
  // ifstream i_file ("data/extend_directly_genotyped.txt.raw");
  //
  // vector<string> header_tokens;
  //
  // // Check file is open before doing anything else
  // if (!i_file.is_open()) {
  //   cout << "File not open, exiting..." << endl;
  //   return -1;
  // }
  //
  // // Read header line
  // if (getline(i_file, line_buff)) {
  //   stringstream ss(line_buff);
  //   string tok_buf;
  //   while (ss >> tok_buf) {
  //     header_tokens.push_back(tok_buf);
  //   }
  // } else {
  //   cout << "Error reading header, exiting..." << endl;
  //   return -1;
  // }
  //
  // int temp_counter = 0;
  // // Read reamining data from file
  // while (getline(i_file, line_buff)) {
  //   temp_counter++;
  //   cout << temp_counter << endl;
  //
  //   stringstream ss(line_buff);
  //   string tok_buf; // Buffer for holding tokens
  //   vector<string> tokens; // Will hold each word seperated by space
  //
  //   auto start = chrono::high_resolution_clock::now();
  //   while (ss >> tok_buf) {
  //     tokens.push_back(tok_buf);
  //   }
  //   auto stop = chrono::high_resolution_clock::now();
  //   auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
  //   cout << milliseconds.count() << endl;
  //
  //   // Setup initial variables
  //   // https://www.cog-genomics.org/plink2/formats#raw for info on encoding
  //   string FID, IID, PAT, MAT, SEX, PHENOTYPE;
  //   FID = tokens[0];
  //   IID = tokens[1];
  //   PAT = tokens[2];
  //   MAT = tokens[3];
  //   SEX = tokens[4];
  //   PHENOTYPE = tokens[5];
  //
  //   string UID = FID + IID; // Uniquely identifies an individual
  //
  //   start = chrono::high_resolution_clock::now();
  //   unordered_map<string, int> pat_data;
  //   for (unsigned int i = 6; i < tokens.size(); i++) {
  //     try {
  //       pat_data[header_tokens[i]] = stoi(tokens[i]);
  //     } catch (exception& e) {
  //       pat_data[header_tokens[i]] = -1;
  //     }
  //   }
  //   stop = chrono::high_resolution_clock::now();
  //   milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
  //   cout << milliseconds.count() << endl;
  //
  //   patient_data_container[UID] = pat_data;
  // }
  //
  // i_file.close();

  return 0;
}
