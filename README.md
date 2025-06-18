# Variance Component and Heritability Estimation Using Additive and Dominance Relationship Matrices

## Overview

This R script performs the following tasks:

1. **Processes Genotyping-By-Sequencing (GBS) HapMap data** to create numeric genotype matrices.
2. **Calculates additive (A) and dominance (D) genomic relationship matrices**.
3. **Fits mixed linear models** using the `sommer` package to estimate:
   - Additive genetic variance
   - Dominance variance
   - Residual variance
   - Narrow-sense heritability (h²)
4. **Loops through multiple traits** in a phenotype dataset to summarize variance components and heritability estimates.
5. **Exports results** to a `.csv` file.

---
