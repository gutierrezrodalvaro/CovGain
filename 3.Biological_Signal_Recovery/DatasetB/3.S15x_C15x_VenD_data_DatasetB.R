# Rscript to generate the counts for the eulerr diagrams
# 09/06/2026
# Alvaro Gutierrez Rodriguez

library(methylKit)

load('DMCs_F5.Rdata')

get_coords <- function(mDiff) {
  paste(mDiff$chr, mDiff$start, sep = "_")
}

# Define the common background, CpGs covered in common
coords_S15_bg <- get_coords(S15_diff)
coords_C15_bg <- get_coords(C15_diff)
common_background <- intersect(coords_S15_bg, coords_C15_bg)

# Obtain coords from the DMCs
all_dmc_S15 <- get_coords(S15_DMCs)
all_dmc_C15 <- get_coords(C15_DMCs)

# Filter common DMCs in the background
dmc_S15_in_bg <- all_dmc_S15[all_dmc_S15 %in% common_background]
dmc_C15_in_bg <- all_dmc_C15[all_dmc_C15 %in% common_background]

# Calculate the common DMCs and the unique per comparison ones
shared_DMCs <- length(intersect(dmc_S15_in_bg, dmc_C15_in_bg))
unique_S15  <- length(setdiff(dmc_S15_in_bg, dmc_C15_in_bg))
unique_C15  <- length(setdiff(dmc_C15_in_bg, dmc_S15_in_bg))

total_S15_in_bg <- length(dmc_S15_in_bg)
total_C15_in_bg <- length(dmc_C15_in_bg)

# Dataframe final
summary_table <- data.frame(
  Comparison = c("S15 DMCs (in Background)", "C15 DMCs (in Background)"),
  Total_DMCs = c(total_S15_in_bg, total_C15_in_bg),
  Shared     = c(shared_DMCs, shared_DMCs),
  Unique     = c(unique_S15, unique_C15)
)

print(summary_table, row.names = FALSE)


## DMCs unique of C15 in S30?

all_dmc_S30 <- get_coords(S30_DMCs)
unique_C15_coords <- setdiff(dmc_C15_in_bg, dmc_S15_in_bg)
dmcs_C15_unicos_en_S30 <- intersect(unique_C15_coords, all_dmc_S30)
num_coincidencias <- length(dmcs_C15_unicos_en_S30)

cat(length(unique_C15_coords), "DMCs unique of C15,", num_coincidencias, "are present in S30_DMCs.\n")
