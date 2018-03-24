#include <iostream>
#include <fstream>
#include <vector>

#include "input.h"

#define TEST 0

using namespace std;

void handle_input_files(
    vector<string> *fam_range,
    vector<string> *bim_range) {

  // Open the .fam file and read in to vector
  #if TEST
    ifstream fam_file ("data/test_data/binary.fam");
  #else
    ifstream fam_file ("data/extend_directly_genotyped_binary.fam");
  #endif
  if (!fam_file.is_open()) {
    cerr << "Unable to open fam file, exiting..." << endl;
  }
  string buf;
  while (getline(fam_file, buf, '\n')) {
    fam_range->push_back(buf);
  }
  fam_file.close();


  // Open the .bim file and read in to vector
  #if TEST
    ifstream bim_file ("data/test_data/binary.bim");
  #else
    ifstream bim_file ("data/extend_directly_genotyped_binary.bim");
  #endif
  if (!bim_file.is_open()) {
    cerr << "Unable to open bim file, exiting..." << endl;
  }
  while (getline(bim_file, buf, '\n')) {
    bim_range->push_back(buf);
  }
  bim_file.close();

}
