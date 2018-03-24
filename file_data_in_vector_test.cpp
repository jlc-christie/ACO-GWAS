#include <iostream>
#include <fstream>
#include <string>
#include <vector>

int main() {
  std::string header_string;
  std::vector<std::string> patient_data_container;

  std::string line_buff;
  std::ifstream i_file ("data/extend_directly_genotyped.txt.raw");

  // Read data from file and save to vector
  if (i_file.is_open()) {
    // Read data from file
    while (getline(i_file, line_buff)) {
      std::string temp = line_buff;
      patient_data_container.push_back(temp);
    }

    i_file.close();
  }

  // Check size of contents of vector
  unsigned long int bytes_used = 0;
  for (unsigned int i = 0; i < patient_data_container.size(); i++) {
    unsigned int len = patient_data_container[i].size();
    std::cout << "Length of string[" << i << "] = " << len << std::endl;
    bytes_used += len;
  }

  std::cout << "Total bytes used = " << bytes_used << std::endl;
  std::cout << "Size of vector = ";
  std::cout << patient_data_container.size()*sizeof(patient_data_container[0]);
  std::cout << std::endl;
}
