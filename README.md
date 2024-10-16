This R script conducts genetic analysis using SNP (Single Nucleotide Polymorphism) data from a panel of 75 genotypes, aiming to estimate genetic variance components and narrow-sense heritability for multiple traits. 
The script begins by importing SNP genotype data in HamMap format and corresponding phenotypic measurements. It then converts the genotype data from character format to numeric codes and calculates genetic relationship matrices (GRM) for both additive and dominance effects using the rrBLUP package. 
The script employs mixed model analysis with the sommer package to estimate genetic variance components for each trait in a predefined list. 
Finally, it outputs the results, including additive and dominance variances, residual variances, and heritability estimates, to a data frame for easy reference, which can be saved as an Excel file. 
This comprehensive approach allows researchers to gain insights into the genetic architecture of the traits under study.
