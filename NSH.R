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

## Let's try modeling for additive and dominance variance in a single mixed model 
# I have observed that for some traits either additive or dominance variance turns out to be zero when modeled together 
# To avoid this, I opted for two different mix models, one each for additive and dominance variance  
# Residual variance of either model can be used to calculate heritability 

# Fit a mix model for additive and dominance variance for a single trait
model.AD <- mmes(
  fixed = PODW ~ 1,
  random = ~ vsm(ism(geno), Gu = A$A) +
    vsm(ism(geno), Gu = D$D),
  rcov = ~ units,
  nIters = 10,
  data = pheno,
  verbose = FALSE
)
(summary(model.AD)$varcomp)

#calculate narrow sense heritability. V1=additive variance, V2=dominance variance, V3=residual variance  
podw <- vpredict(model.AD, h2 ~ (V1) / ( V1+V3) ) 
print(podw)

## Now, let's model additive and dominance variance in a separate mix model 

# Fit a mixed model using the GRM for additive effects on a single trait
model_A <- mmes(fixed = PODW ~ 1,  # Fixed effect for the trait
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







# Loop over all traits in the phenotypic dataset and calculate heritability
trait_names <- colnames(pheno)[2:21]  # Select all traits except 'geno'
results <- data.frame(trait = character(),
                      additive_variance = numeric(),
                      residual_variance = numeric(),
                      heritability = numeric(), stringsAsFactors = FALSE)

for (trait in trait_names) {
  result <- calculate_heritability(trait, pheno, G)
  results <- rbind(results, as.data.frame(t(result)))
}

# Correct column names for the results
colnames(results) <- c("Trait", "Additive_Variance", "Residual_Variance", "Heritability")

# Convert numeric columns for proper formatting
results$Additive_Variance <- as.numeric(as.character(results$Additive_Variance))
results$Residual_Variance <- as.numeric(as.character(results$Residual_Variance))
results$Heritability <- as.numeric(as.character(results$Heritability))

# Save the results to a CSV file for future reference
write.csv(results, "heritability_results.csv", row.names = FALSE)

# Print the heritability results
print(results)


#Now, let's calculate the dominance variance. 
#Initially, I attempted to use the same model for extraction, but it proved challenging to obtain the variance directly. 
#Therefore, I decided to calculate it separately using the script below


# Load the necessary library for mixed models
library(sommer)

# Function to calculate dominance variance for a given trait
# The function takes in a trait name, the phenotypic data, and the dominance GRM (D).
calculate_dominance_variance <- function(trait, data, D) {
  # Fit the mixed model for the specified trait, with dominance effect using the GRM (D$D)
  model <- mmer(fixed = as.formula(paste(trait, "~ 1")),  # Fixed effect for the trait
                random = ~ vsr(geno, Gu = D$D),           # Random effect for dominance using the GRM
                rcov = ~ units,                           # Residual variance
                data = data)                              # The phenotypic data
  
  # Extract variance components from the fitted model
  varcomp <- summary(model)$varcomp
  
  # Extract dominance genetic variance (sigma_D) and residual variance (sigma_e) from the model
  sigma_D <- varcomp[paste0("u:geno.", trait, "-", trait), "VarComp"]  # Dominance variance
  sigma_e <- varcomp[paste0("units.", trait, "-", trait), "VarComp"]   # Residual variance
  
  # Return the trait name, dominance variance, and residual variance as a named vector
  return(c(trait = trait, dominance_variance = sigma_D, residual_variance = sigma_e))
}

# Create an empty data frame to store the dominance variance results
dominance_results <- data.frame(trait = character(),              # Column for trait names
                                dominance_variance = numeric(),   # Column for dominance variance values
                                residual_variance = numeric(),    # Column for residual variance values
                                stringsAsFactors = FALSE)

# List of all trait names in the phenotypic data (excluding the first column 'geno')
trait_names <- colnames(pheno)[2:21]

# Loop through all traits in the phenotypic data
for (trait in trait_names) {
  # Call the function for each trait to calculate dominance and residual variances
  result <- calculate_dominance_variance(trait, pheno, D)
  # Convert the result to a data frame and append it to the dominance_results data frame
  dominance_results <- rbind(dominance_results, as.data.frame(t(result)))
}

# Convert the columns for dominance variance and residual variance to numeric format
dominance_results$dominance_variance <- as.numeric(as.character(dominance_results$dominance_variance))
dominance_results$residual_variance <- as.numeric(as.character(dominance_results$residual_variance))

# Save the final results to a CSV file for later reference or sharing
write.csv(dominance_results, "dominance_variance_results.csv", row.names = FALSE)

# Print the results to check the output
print(dominance_results)












