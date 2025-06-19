# Set the working directory
setwd("C:/Users/ehtis/OneDrive - New Mexico State University/SUNNY/Research Projects/Fruit Morphology Projects/GWAS/NSH/yield")
# Verify the current working directory
getwd() 

# Load the HapMap formatted GBS data, skipping the first 10 lines that contain metadata
# This is because my GBS file is in HamMap format which includes non-genotypic header lines.
hapmap_data <- read.table("GWASD125_KNNimp_beagle_imputed.hmp.txt", header = FALSE, skip = 10)

# Preview the structure of the loaded data to ensure it's loaded correctly
head(hapmap_data)
str(hapmap_data)

# Load necessary library for data manipulation
library(dplyr)

# Extract the SNP genotype data starting from the 12th column (based on HapMap format)
# Columns before V12 contain marker metadata (chromosome, position, etc.)
genotype_data <- hapmap_data[, 12:ncol(hapmap_data)]
str(genotype_data)  # Check the structure of the extracted genotype data

head (genotype_data)
# Conversion of genotype calls to numerical format:
# AA = 0 (homozygous reference)
# AG/GA = 1 (heterozygous)
# GG = 2 (homozygous alternate)
# Missing data (NAs) are retained as NA

# Define the conversion function
convert_genotype <- function(genotype) {
  if (is.na(genotype)) {
    return(NA)  # Keep missing data as NA
  } else if (genotype == "AA") {
    return(0)
  } else if (genotype == "AG" | genotype == "GA") {
    return(1)
  } else if (genotype == "GG") {
    return(2)
  } else {
    return(NA)  # Handle unrecognized genotypes
  }
}

# Apply the conversion to the entire genotype dataset
numeric_genotypes <- apply(genotype_data, MARGIN = c(1, 2), FUN = convert_genotype)

# Transpose the numeric genotypes for GRM calculation (genotypes as rows, SNPs as columns)
numeric_genotypes <- t(numeric_genotypes)

# Load the necessary library for GRM calculation
library(rrBLUP)
library(sommer)

# Calculate the Additive Relationship Matrix (GRM)
A <- A.mat(numeric_genotypes, return.impute = TRUE)
str(A)  # Check the structure of the GRM
str(A$A)  # Check the structure of the matrix inside GRM

# Optionally, calculate Dominance Relationship Matrix (D) and Epistatic Matrix (E)
D <- D.mat(numeric_genotypes, return.imputed = TRUE)  # Dominance matrix
str(D)
str(D$D)

# Epistatic matrix optional. I am skipping it 
E <- E.mat(numeric_genotypes)  

# Extract genotype names from the GRM
geno_names_grm <- rownames(A$A)


# Read the phenotypic data (IMPORTANT: genotype names in this file must match GRM names: V12, V13, ..., V136)
# If your original file uses different genotype names like "GBS001", "GBS002", ..., you should manually rename them
# to match GRM names ("V12", "V13", ...) in the CSV before proceeding.
# This is straightforward if done manually using Excel or any text editor by replacing GBS IDs with V12–V136 in order.

pheno <- read.csv("mydata_yield.csv")  # Make sure 'geno' column contains V12, V13, ..., V136
pheno$geno <- as.factor(pheno$geno)    # Convert genotype column to factor

# Extract genotype names from the phenotypic data
geno_names_pheno <- pheno$geno

# Check if genotype names in GRM and phenotypic data match
all(geno_names_pheno == geno_names_grm)  # Should return TRUE if they match


# Load necessary libraries for mixed model analysis
library(Matrix)
library(MASS)
library(crayon)
library(sommer)

### NOTE: 
# When fitting mixed models for certain traits, either the additive or dominance variance 
# may be estimated as zero when modeled together, with all phenotypic variance attributed to residuals. 
# To address this, I fitted models with both additive and dominance effects combined, 
# as well as separate models for additive-only and dominance-only effects. 
# For robust results, compare model fit and choose the approach with the best performance


## Let's try modeling for additive and dominance variance in a single mixed model 

# I have observed that for some traits, either additive or dominance variance turns out to be zero when modeled together. 
# To avoid this, I opted for both approaches, model the additive and dominance together and also separately 
# You can use either approach or use the the model with best results 


# Fit a mix model for additive and dominance variance for a single trait
model.AD <- mmes(
  fixed = GRN ~ 1,
  random = ~ vsm(ism(geno), Gu = A$A) +
    vsm(ism(geno), Gu = D$D),
  rcov = ~ units,
  nIters = 2000,
  data = pheno,
  verbose = TRUE
)
(summary(model.AD)$varcomp)

# Extract variance components manually from model.AD
varcomp_model.AD <- summary(model.AD)$varcomp
print(varcomp_model.AD)

model.AD_sigma_A <- varcomp_model.AD["geno:A:mu:mu", "VarComp"]  # Additive variance
model.AD_sigma_D <- varcomp_model.AD["geno:D:mu:mu", "VarComp"]  # Dominance variance
model.AD_sigma_E <- varcomp_model.AD["units:mu:mu" , "VarComp"]  # Residual variance


# Calculate narrow-sense heritability for single trait 
h2 <- model.AD_sigma_A / (model.AD_sigma_A + model.AD_sigma_D + model.AD_sigma_E)
print(h2)

# you can also use inbuilt function from sommer to calculate narrow sense heritability, however, I opted formula mentioned before   
# V1=additive variance, V2=dominance variance, V3=residual variance  
ARA <- vpredict(model.AD, h2 ~ (V1) / ( V1+V3) ) 
print(ARA)

# Let's use a loop for multiple traits for model.AD 

# Initialize an empty data frame to store results
results_model.AD <- data.frame(Trait = character(),
                      Additive_Var = numeric(),
                      Dominance_Var = numeric(),
                      Residual_Var = numeric(),
                      Heritability = numeric(),
                      stringsAsFactors = FALSE)

# Define the traits to analyze from the imported pheno data
traits <- colnames(pheno)[2:5]

# Loop through each trait
for (trait in traits) {
  # Create formula for the current trait
  formula <- as.formula(paste(trait, "~ 1"))
  
  # Fit mixed model for additive and dominance variance
  model.AD <- mmes(
    fixed = formula,
    random = ~ vsm(ism(geno), Gu = A$A) +
      vsm(ism(geno), Gu = D$D),
    rcov = ~ units,
    nIters = 1000,
    data = pheno,
    verbose = TRUE
  )
  
  # Extract variance components
  varcomp_model.AD <- summary(model.AD)$varcomp
  sigma_A <- varcomp_model.AD["geno:A:mu:mu", "VarComp"]  # Additive variance
  sigma_D <- varcomp_model.AD["geno:D:mu:mu", "VarComp"]  # Dominance variance
  sigma_E <- varcomp_model.AD["units:mu:mu", "VarComp"]  # Residual variance
  
  # Calculate narrow-sense heritability
  h2 <- sigma_A / (sigma_A + sigma_D + sigma_E)
  
  # Store results in data frame
  results_model.AD <- rbind(results_model.AD, data.frame(
    Trait = trait,
    Additive_Var = sigma_A,
    Dominance_Var = sigma_D,
    Residual_Var = sigma_E,
    Heritability = h2
  ))
}

# Save results to CSV
write.csv(results_model.AD, "model.AD_heritability_results.csv", row.names = FALSE)

# Print results
print(results_model.AD)

## Now, let's model additive and dominance variance in a separate mix model 

# Fit a mixed model using the GRM for additive effects on a single trait
model_A <- mmes(fixed = GRN ~ 1,  # Fixed effect for the trait
                random = ~ vsm(ism(geno), Gu = A$A),  # Random effect using the GRM (Additive variance)
                rcov = ~ units,  # Residual variance
                data = pheno)  # Phenotypic data
summary(model_A)

# Extract variance components manually from model_A
varcomp_model_A <- summary(model_A)$varcomp
print(varcomp_model_A)

# Extract additive genetic variance (σA=AA) and residual variance (σe=AE) from model_A
sigma_AA <- varcomp_model_A["geno:A:mu:mu", "VarComp"]  # Additive genetic variance
print(sigma_AA)
sigma_AE <- varcomp_model_A["units:mu:mu", "VarComp"]    # Residual variance
print(sigma_AE)

# Fit a mixed model using the GRM for dominance effects on a single trait 
model_D <- mmes(fixed = PODW ~ 1,  # Fixed effect for the trait
                random = ~ vsm(ism(geno), Gu = D$D),  # Random effect using the GRM (Dominance variance)
                rcov = ~ units,  # Residual variance
                data = pheno)  # Phenotypic data
summary(model_D)

# Extract variance components from model_D
varcomp_model_D <- summary(model_D)$varcomp
print(varcomp_model_D)

# Extract dominance variance (sigma_DD) and residual variance (sigma_DE) from model_D
sigma_DD <- varcomp_model_D["geno:D:mu:mu", "VarComp"]   # Dominance genetic variance
print(sigma_DD)
sigma_DE <- varcomp_model_D["units:mu:mu", "VarComp"]    # Residual variance
print(sigma_DE)

# Calculate narrow-sense heritability for single trait 
h2_PODW <- sigma_AA / (sigma_AA + sigma_DD + sigma_DE)
print(h2_PODW)


# Load required library 
library(sommer)

# Define the traits to analyze from the imported pheno data by column numbers
traits <- colnames(pheno)[2:5]  

# Initialize a data frame to store results
results <- data.frame(
  Trait = character(),
  Sigma_AA = numeric(),
  Sigma_AE = numeric(),
  Sigma_DD = numeric(),
  Sigma_DE = numeric(),
  Heritability = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each trait
for (trait in traits) {
  # Fit mixed model for additive effects
  model_A <- mmes(
    fixed = as.formula(paste(trait, "~ 1")),
    random = ~ vsm(ism(geno), Gu = A$A),
    rcov = ~ units,
    data = pheno
  )
  
  # Extract variance components for additive model
  varcomp_A <- summary(model_A)$varcomp
  sigma_AA <- varcomp_A["geno:A:mu:mu", "VarComp"]  # Additive genetic variance
  sigma_AE <- varcomp_A["units:mu:mu", "VarComp"]   # Residual variance
  
  # Fit mixed model for dominance effects
  model_D <- mmes(
    fixed = as.formula(paste(trait, "~ 1")),
    random = ~ vsm(ism(geno), Gu = D$D),
    rcov = ~ units,
    data = pheno
  )
  
  # Extract variance components for dominance model
  varcomp_D <- summary(model_D)$varcomp
  sigma_DD <- varcomp_D["geno:D:mu:mu", "VarComp"]  # Dominance genetic variance
  sigma_DE <- varcomp_D["units:mu:mu", "VarComp"]   # Residual variance
  
  # Calculate narrow-sense heritability
  h2 <- sigma_AA / (sigma_AA + sigma_DD + sigma_DE)
  
  # Store results in data frame
  results <- rbind(results, data.frame(
    Trait = trait,
    Sigma_AA = sigma_AA,
    Sigma_AE = sigma_AE,
    Sigma_DD = sigma_DD,
    Sigma_DE = sigma_DE,
    Heritability = h2
  ))
}

# Save results to a CSV file
write.csv(results, "variance_components_heritability.csv", row.names = FALSE)

# Print results to console
print(results)


## Note on Variance Components for Mixed Models

# The mixed model for one of the traits (using mmes from sommer) estimates 
# zero additive and dominance variances, with all phenotypic variance 
# attributed to residuals. This is likely due to low genetic variation for the 
# trait, a small sample size (125 genotypes), and near-singular relationship matrices
# Simplifying to an additive-only model, regularizing matrices, and 
# log transformation did not resolve the issue.

## The model code is correct, but the data may lack sufficient genetic signal.
