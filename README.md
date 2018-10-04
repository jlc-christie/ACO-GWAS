# ACO-GWAS
Code used to carry our my final year project titled "*Using Ant Colony Optimization to explore gene-gene interactions associated with Type 2 Diabetes*" as part of my BSc Computer Science degree at the University of Exeter. 

This repo will not be maintained, but if anyone wants to turn this in to a polished piece of software feel free to submit pull requests. 

The issues I'm aware of have been added to issues, please read these and the warning below before using.

What it can do:
- Read in BED, BIM, FAM files (plink format) and perform various checks
- Carry out Chi-Squared test using cases and controls for each SNP, computed using allele counts (2 dof)
- Carry out SNP-SNP interaction test using Chi-Squared test (8 dof)
- Carry out a basic case/control GWAS

WARNING: 
This tool was created to carry out specific tests for **specific data**, rather than a generic tool to carry out epistatic interaction studies using ACO. The code has not been tested on other datasets, there are known bugs and it will likely not be maintained, **do not use this tool to produce results for published work**.
