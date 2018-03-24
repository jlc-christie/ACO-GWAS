#include <cmath>
#include <iostream>
#include <iomanip>

double chi_sq(int counts[8]) {
  // Create a new array of doubles so floating point calculations can be done
  double d_counts[8];
  for (int i = 0; i < 8; i++) {
    d_counts[i] = double(counts[i]);
  }

  // Store array of allele counts rather than genotype counts
  // [0] -> A1 Case
  // [1] -> A2 Case
  // [2] -> A1 Control
  // [3] -> A2 Control
  double allele_counts[4];
  allele_counts[0] = 2*d_counts[0] + d_counts[1];
  allele_counts[1] = d_counts[1] + 2*d_counts[2];
  allele_counts[2] = 2*d_counts[4] + d_counts[5];
  allele_counts[3] = d_counts[5] + 2*d_counts[6];

  // Sum up the rows and columns of contingency table
  double total_a1, total_a2, total_cases, total_controls, total;
  total_a1 = allele_counts[0] + allele_counts[2];
  total_a2 = allele_counts[1] + allele_counts[3];
  total_cases = allele_counts[0] + allele_counts[1];
  total_controls = allele_counts[2] + allele_counts[3];

  total = total_cases + total_controls;

  // Calculate expected values, follows same order as allele counts array
  double expected[4];
  expected[0] = (total_cases * total_a1) / total;
  expected[1] = (total_cases * total_a2) / total;
  expected[2] = (total_controls * total_a1) / total;
  expected[3] = (total_controls * total_a2) / total;

  double res = 0.0;
  res += pow(allele_counts[0] - expected[0], 2) / expected[0];
  res += pow(allele_counts[1] - expected[1], 2) / expected[1];
  res += pow(allele_counts[2] - expected[2], 2) / expected[2];
  res += pow(allele_counts[3] - expected[3], 2) / expected[3];

  return res;

}
//
// double chi_sq(int a, int b, int c,
//             int d, int e, int f) {
//   /**
//    * Takes 6 ints to represent the counts for each different cell as
//    * represented below, fills res array with 2 different chi-squared values
//    * for the major and minor alleles. mcase and mcontrol represent the number
//    * of missing cases and controls respectively for that SNP.
//    *
//    *             GG    GT    TT
//    * Cases       a     b     c
//    * Controls    d     e     f
//    */
//
//    // Create doubles for a,b,c,d,e,f so calculations can be performed
//    double da, db, dc, dd, de, df;
//    da = double(a);
//    db = double(b);
//    dc = double(c);
//    dd = double(d);
//    de = double(e);
//    df = double(f);
//
//    // Create new variables to count alleles rather than genotypes
//    double a1_case, a2_case, a1_control, a2_control;
//    a1_case = (2*da) + db;
//    a2_case = db + (2*dc);
//    a1_control = (2*dd) + de;
//    a2_control = de + (2*df);
//
//    // Calculate column and row totals
//    double total_cases, total_controls, total_a1, total_a2, total;
//    total_cases = a1_case + a2_case;
//    total_controls = a1_control + a2_control;
//    total_a1 = a1_case + a1_control;
//    total_a2 = a2_case + a2_control;
//
//    // Equally this could be -> total_a1 + total_a2
//    total = total_cases + total_controls;
//
//    double exp_a1_case, exp_a2_case, exp_a1_control, exp_a2_control;
//    exp_a1_case = (total_cases * total_a1) / total;
//    exp_a2_case = (total_cases * total_a2) / total;
//    exp_a1_control = (total_controls * total_a1) / total;
//    exp_a2_control = (total_controls * total_a2) / total;
//
//    double res = 0.0;
//    res += pow(a1_case - exp_a1_case, 2) / exp_a1_case;
//    res += pow(a2_case - exp_a2_case, 2) / exp_a2_case;
//    res += pow(a1_control - exp_a1_control, 2) / exp_a1_control;
//    res += pow(a2_control - exp_a2_control, 2) / exp_a2_control;
//
//    return res;
// }
