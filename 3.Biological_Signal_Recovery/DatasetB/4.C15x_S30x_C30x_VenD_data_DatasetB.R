# This script caclulates the counts on unique and overlapping DMCs across the three main approaches
# Alvaro Gutierrez Rodriguez
# 09/06/2026

library(methylKit)

load('DMCs_F5.Rdata')

# Function to define chr cords
get_coords <- function(mDiff) {
  paste(mDiff$chr, mDiff$start, sep = "_")
}

# Triple intersection to determine the common background between datasets
coords_C15 <- get_coords(C15_diff)
coords_C30 <- get_coords(C30_diff)
coords_S30 <- get_coords(S30_diff)

# Extract coords of the DMCs
common_background <- intersect(intersect(coords_C15, coords_C30), coords_S30)
message(paste("Fondo común de CpGs:", length(common_background)))

# Extract coords of the DMCs
all_C15 <- get_coords(C15_DMCs)
all_C30 <- get_coords(C30_DMCs)
all_S30 <- get_coords(S30_DMCs)

set1 <- all_C15[all_C15 %in% common_background]
set2 <- all_C30[all_C30 %in% common_background]
set3 <- all_S30[all_S30 %in% common_background]

# Determine the regions for the eurlerr diagram
shared_all <- length(intersect(intersect(set1, set2), set3))

# Final Dataframe
shared_C15_C30 <- length(setdiff(intersect(set1, set2), set3))
shared_C15_S30 <- length(setdiff(intersect(set1, set3), set2))
shared_C30_S30 <- length(setdiff(intersect(set2, set3), set1))

unique_C15 <- length(setdiff(set1, union(set2, set3)))
unique_C30 <- length(setdiff(set2, union(set1, set3)))
unique_S30 <- length(setdiff(set3, union(set1, set2)))

total_C15 <- length(set1)
total_C30 <- length(set2)
total_S30 <- length(set3)

# Final Dataframe
summary_table <- data.frame(
  Dataset    = c("C15 DMCs", "C30 DMCs", "S30 DMCs"),
  Total_In_BG = c(total_C15, total_C30, total_S30),
  Strict_Unique = c(unique_C15, unique_C30, unique_S30),
  Shared_With_C15 = c("-", shared_C15_C30, shared_C15_S30),
  Shared_With_C30 = c(shared_C15_C30, "-", shared_C30_S30),
  Shared_With_S30 = c(shared_C15_S30, shared_C30_S30, "-"),
  Shared_By_All   = c(shared_all, shared_all, shared_all)
)

print(summary_table, row.names = FALSE)
