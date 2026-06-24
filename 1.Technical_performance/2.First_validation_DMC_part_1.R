# 08/06/2026
# Alvaro Gutierrez Rodriguez
# With this script we calculate de common CpGs and the DMCs between the 15x and 30x datatset with both methodologies (STD & CTM)
# Min diff > 15 %, p value < 0.05
# Table Figure 2 E

library(methylKit)

load('norm_objs_min5.Rdata')

common_counts <- numeric(5)
dmc_counts    <- numeric(5)

### Effect of Coverage: Same methodology, different coverage.
 # 1. STD Method Coverage Effect Comparison: STD15x vs. STD30x

all_STD_samples <- c(STD15x_normalized, STD30x_normalized)
treatment_vector <- c(rep(0, 7), rep(1, 7))
combined_STD <- new("methylRawList", all_STD_samples, treatment = treatment_vector)
meth_unite_STD <- unite(combined_STD, destrand = T)

common_counts[1] <- nrow(meth_unite_STD)
print("Number of CpGs common STD15x vs. STD30x")
print(nrow(meth_unite_STD))

myDiff_STD=calculateDiffMeth(meth_unite_STD,mc.cores=4)
myDiff15p_STD=getMethylDiff(myDiff_STD,difference=15,qvalue=0.05)
dmc_counts[1] <- nrow(myDiff15p_STD)
print("Number of DMCs STD15x vs. STD30x")
print(nrow(myDiff15p_STD))

rm(all_STD_samples, combined_STD, meth_unite_STD, d, myDiff_STD, myDiff15p_STD); gc()

### Effect of Coverage: Same methodology, different coverage.
 # 2. CTM Method Coverage Effect Comparison: CTM15x vs. CTM30x

all_CTM_samples <- c(CTM15x_normalized, CTM30x_normalized)
treatment_vector <- c(rep(0, 7), rep(1, 7))
combined_CTM <- new("methylRawList", all_CTM_samples, treatment = treatment_vector)
meth_unite_CTM <- unite(combined_CTM, destrand = T)

common_counts[2] <- nrow(meth_unite_CTM)
print("Number of CpGs coomon CTM15x vs. CTM30x")
print(nrow(meth_unite_CTM))

myDiff_CTM=calculateDiffMeth(meth_unite_CTM,mc.cores=4)
myDiff15p_CTM=getMethylDiff(myDiff_CTM,difference=15,qvalue=0.05)
dmc_counts[2] <- nrow(myDiff15p_CTM)
print("Number of DMCs CTM15x vs. CTM30x")
print(nrow(myDiff15p_CTM))

rm(all_CTM_samples, combined_CTM, meth_unite_CTM, myDiff_CTM, myDiff15p_CTM); gc()

### Effect of Method: Same coverage, different methodology
# 3. 15x Datset Method Effect Comparison: STD15x vs. CTM15x

all_15x_samples <- c(STD15x_normalized, CTM15x_normalized)
treatment_vector <- c(rep(0, 7), rep(1, 7))
combined_15x <- new("methylRawList", all_15x_samples, treatment = treatment_vector)
meth_unite_15x <- unite(combined_15x, destrand = T)

common_counts[3] <- nrow(meth_unite_15x)
print("Number of CpGs common STD15x vs. CTM15x")
print(nrow(meth_unite_15x))

myDiff_15x=calculateDiffMeth(meth_unite_15x,mc.cores=4)
myDiff15p_15x=getMethylDiff(myDiff_15x,difference=15,qvalue=0.05)
dmc_counts[3] <- nrow(myDiff15p_15x)
print("Number of DMCs STD15x vs. CTM15x")
print(nrow(myDiff15p_15x))

rm(all_15x_samples, combined_15x, meth_unite_15x, myDiff_15x, myDiff15p_15x); gc()

### Effect of Method: Same coverage, different methodology
# 4. 30x Dataset Method Effect Comparison: STD30x vs. CTM30x

all_30x_samples <- c(STD30x_normalized, CTM30x_normalized)
treatment_vector <- c(rep(0, 7), rep(1, 7))
combined_30x <- new("methylRawList", all_30x_samples, treatment = treatment_vector)
meth_unite_30x <- unite(combined_30x, destrand = T)

common_counts[4] <- nrow(meth_unite_30x)
print("Number of CpGs common STD30x vs. CTM30x")
print(nrow(meth_unite_30x))

myDiff_30x=calculateDiffMeth(meth_unite_30x,mc.cores=1)
myDiff15p_30x=getMethylDiff(myDiff_30x,difference=15,qvalue=0.05)
dmc_counts[4] <- nrow(myDiff15p_30x)
print("Number of DMCs STD30x vs. CTM30x")
print(nrow(myDiff15p_30x))

rm(all_30x_samples, combined_30x, meth_unite_30x, myDiff_30x, myDiff15p_30x); gc()

# FINAL COMPARISON: CTM15x vs. STD30x. How similar are the CTM low cov dataset to the STD 30x

CTM15x_STD30x_samples <- c(CTM15x_normalized, STD30x_normalized)
treatment_vector <- c(rep(0, 7), rep(1, 7))
CTM15x_STD30x_combined <- new("methylRawList", CTM15x_STD30x_samples, treatment = treatment_vector)
CTM15x_STD30x_meth_unite <- unite(CTM15x_STD30x_combined, destrand = T)

common_counts[5] <- nrow(CTM15x_STD30x_meth_unite)
print("Number of CpGs common CTM15x vs. STD30x")
print(nrow(CTM15x_STD30x_meth_unite))

myDiff_CTM15x_STD30x=calculateDiffMeth(CTM15x_STD30x_meth_unite, mc.cores=1)
myDiff15p_CTM15x_STD30x=getMethylDiff(myDiff_CTM15x_STD30x, difference=15,qvalue=0.05)
dmc_counts[5] <- nrow(myDiff15p_CTM15x_STD30x)
print("Number of DMCs CTM15x vs. STD30x")
print(nrow(myDiff15p_CTM15x_STD30x))
rm(CTM15x_STD30x_samples, CTM15x_STD30x_combined, CTM15x_STD30x_meth_unite, myDiff_CTM15x_STD30x, myDiff15p_CTM15x_STD30x); gc()

# FINAL SUMMARY TABLE
summary_table <- data.frame(
  Comparison  = c("STD 15x vs 30x", "CTM 15x vs 30x", "15x: STD vs CTM", "30x: STD vs CTM", "CTM 15x vs STD 30x"),
  Common_CpGs = common_counts,
  DMCs_Found  = dmc_counts
)

print(summary_table)
