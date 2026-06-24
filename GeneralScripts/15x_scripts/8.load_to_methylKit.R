#!/usr/bin/env Rscript
# This script load the sorted and indexed 15x bam files to methylKit

library(methylKit)

# Loading the data
src <- "/data2/gonzalo2/alvaro_data2/dlabrax/15x/dedup_sorted_indexed_15x"
file.list <- list.files(path = src, pattern = "\\.bam$", full.names = TRUE)
sample.id <- sub("_.*$", "", basename(file.list))
file.list <- as.list(file.list)
sample.id <- as.list(sample.id)

# We create the first object in methylKit (myobj)
myobj <- processBismarkAln(
  location = file.list,                     
  sample.id = sample.id,               
  save.folder = "/data2/gonzalo2/alvaro_data2/dlabrax/15x/Rdata_15x",
  save.context = "CpG",
  treatment = c(1,3,3,3,1,1,1,2,2,2,2,3),
  assembly = "seabass_V1.0",           # Genome assembly (Tine et al., 2014. Nature Communications)
  read.context = "CpG",
  mincov = 5,                    # Minimum coverage to retain sites (5x)
  )

# Save the resulting methylation object to a file
save(myobj, file = "/data2/gonzalo2/alvaro_data2/dlabrax/15x/Rdata_15x/myobjdlabrax_15x.RData")
