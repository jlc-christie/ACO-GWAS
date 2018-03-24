#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <stdio.h>
#include <string.h>

int main() {
  std::vector<void*> pointer_container;

  int mem_size = 0;

  while (mem_size != -1) {
    std::cout << "How much memory? (MBs)" << std::endl;
    std::cin >> mem_size;
    std::cout << std::endl;

    if (mem_size > 0) {
      mem_size = int(mem_size*1000000); // Convert to megabytes
      void* ptr = ::operator new(mem_size);
      memset(ptr, '0', mem_size);
      pointer_container.push_back(ptr);
    } else if (mem_size == -1) {
      std::cout << "Exiting..." << std::endl;
      break;
    } else {
      std::cout << "Can't allocate negative memory!\nTry again..." << std::endl;
      continue;
    }
  }

  // Deallocate memory
  for (unsigned int i = 0; i < pointer_container.size(); i++) {
    std::cout << "Deleting dynamically alocated memory..." << std::endl;
    void* ptr = pointer_container[i];
    ::operator delete(ptr);
  }

  // std::cout << "Counting lines in file..." << std::endl;
  //
  // std::string line_buff;
  // std::ifstream i_file ("data/extend_directly_genotyped.txt.raw");
  //
  // if (i_file.is_open()) {
  //   while ( getline(i_file, line_buff) ) {
  //     std::cout << line_buff.substr(0,60) << std::endl;
  //   }
  // }
  //
  // i_file.close();
}
