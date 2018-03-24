#include <stdio.h>
#include <stdlib.h>

void read_byte(char byte, int n, int* h1, int* h2, int* het, int* miss) {
  
}

int main(void) {
  FILE* bed = fopen("../data/test_data/binary.bed", "rb");
  FILE* bim = fopen("../data/test_data/binary.bim", "r");
  FILE* fam = fopen("../data/test_data/binary.fam", "r");

  int individual_count = 0;
  int snp_count = 0;

  while (!feof(fam)) {
    if (fgetc(fam) == '\n') {
      individual_count++;
    }
  }

  while (!feof(bim)) {
    if (fgetc(bim) == '\n') {
      snp_count++;
    }
  }

  printf("Individual Count: %d\n", individual_count);
  printf("SNP Count: %d\n", snp_count);


}
