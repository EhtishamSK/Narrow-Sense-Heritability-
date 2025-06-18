GWAS Analysis for Fruit Morphology and Yield
This repository provides a reproducible workflow for performing a Genome-Wide Association Study (GWAS) using Genotyping-by-Sequencing (GBS) data in HapMap format and phenotypic data to estimate additive and dominance variance components and narrow-sense heritability for multiple traits.
Overview
The R script processes GBS data, converts genotypes to numerical format, calculates Additive and Dominance Relationship Matrices (GRM and DRM), and fits mixed models to estimate variance components and heritability. The workflow is designed for a dataset with genotypes labeled V12, V13, ..., V136 and phenotypic data in CSV format.
Key steps include:

Loading and processing HapMap-formatted GBS data.
Converting genotypes to numerical values (AA=0, AG/GA=1, GG=2, missing=NA).
Computing GRM and DRM using the rrBLUP and sommer packages.
Loading phenotypic data and ensuring genotype name consistency.
Fitting mixed models to estimate additive (σ_A) and dominance (σ_D) variances.
Calculating narrow-sense heritability (h² = σ_A / (σ_A + σ_D + σ_E)).
Analyzing multiple traits and saving results to a CSV file.

This workflow is suitable for genetic variance partitioning, heritability estimation, and understanding the genetic architecture of complex traits in agricultural research.
Prerequisites
Install the required R packages:

dplyr
rrBLUP
sommer
Matrix
MASS
crayon

install.packages(c("dplyr", "rrBLUP", "sommer", "Matrix", "MASS", "crayon"))

Input Files

GBS Data: GWASD125_KNNimp_beagle_imputed.hmp.txt

HapMap format with SNP genotypes.
First 10 lines are metadata; genotype data starts at line 11.
Columns 1–11 contain marker metadata; genotypes start at column 12.


Phenotypic Data: mydata_yield.csv

CSV file with a geno column (matching GRM names: V12, V13, ..., V136).
Traits to analyze are in columns 2–5.



Important: Genotype names in mydata_yield.csv must match GRM names. Manually rename (e.g., GBS001 to V12) using Excel or a text editor if needed.
Usage

Place input files in the working directory.
Update the setwd() path in the script to your directory.
Verify genotype name consistency between mydata_yield.csv and GRM.
Run the script in R or RStudio.

Output

variance_components_heritability.csv:
Columns: Trait, Sigma_AA (additive variance), Sigma_AE (residual variance, additive model), Sigma_DD (dominance variance), Sigma_DE (residual variance, dominance model), Heritability.



Notes

The script assumes a trait named PODW for single-trait analysis. Update the trait name as needed.
For multiple traits, columns 2–5 are analyzed. Modify traits <- colnames(pheno)[2:5] for different columns.
Epistatic matrix calculation (E.mat) is commented out to save resources. Uncomment if required.
Separate additive and dominance models are used to avoid zero variance issues in joint models.

Troubleshooting

Genotype Mismatch: If all(geno_names_pheno == geno_names_grm) is FALSE, ensure geno column matches GRM names.
Unrecognized Genotypes: Non-standard genotype calls are set to NA. Check the HapMap file.
Package Issues: Confirm all packages are installed and loaded.


 
