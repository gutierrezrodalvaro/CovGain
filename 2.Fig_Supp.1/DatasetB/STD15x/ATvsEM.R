# EXTERNAL VALIDATION
# This script calculates de number of CpGs between AT and EM populations
# It uses de 15x datasets with the S method
# DATASET B

library(methylKit)

dir_path = "/data2/gonzalo2/alvaro_data2/EMSeq_cov/External_Validation/15x_CpG.reports/"

metadata <- read.table(paste0(dir_path, "Sample_codes_external_validation_coverage.txt"), 
                       header = TRUE, stringsAsFactors = FALSE)
file_list <- as.list(paste0(dir_path, metadata$SampleID, "_CpG_report.txt"))
sample_ids <- as.list(metadata$SampleID)
treatment_vector <- ifelse(metadata$Population == "AT", 0, 1)

myobj_STD15x <- methRead(file_list,
                  sample.id = sample_ids,
                  assembly = "seabass_V1.0",
                  pipeline = "bismarkCytosineReport",
                  treatment = treatment_vector,
                  context = "CpG",
   		  mincov = 5)

filtered.myobj_STD15x=filterByCoverage(myobj_STD15x,lo.count=5,lo.perc=NULL,
                                hi.count=NULL,hi.perc=99.9)

normalized.myobj_STD15x=normalizeCoverage(filtered.myobj_STD15x, method = "median")

rm(myobj_STD15x, filtered.myobj_STD15x); gc()

# DMCs
settings <- list(
#  "12/12 (0%)"   = NULL, 
#  "11/12 (91.6%)" = 11L,
#  "10/12 (83.3%)"= 10L,
  "9/12 (25%)"= 9L
)

final_table <- data.frame()

cat("Starting multi-parameter analysis...\n")

for (label in names(settings)) {
  cat("Processing:", label, "\n")

  meth_tmp <- unite(normalized.myobj_STD15x,
                    destrand = TRUE,
                    min.per.group = settings[[label]],
                    mc.cores = 4)

  diff_tmp <- calculateDiffMeth(meth_tmp, mc.cores = 4)

  dmc_all   <- getMethylDiff(diff_tmp, difference = 15, qvalue = 0.05)
  dmc_hyper <- getMethylDiff(diff_tmp, difference = 15, qvalue = 0.05, type = "hyper")
  dmc_hypo  <- getMethylDiff(diff_tmp, difference = 15, qvalue = 0.05, type = "hypo")

  row_data <- data.frame(
    Comparison   = "AT vs. EM",
    Dataset      = "STD15x",
    Missing_data = label,
    CpGs_common  = nrow(meth_tmp),
    DMCs         = nrow(dmc_all),
    Hyper        = nrow(dmc_hyper),
    Hypo         = nrow(dmc_hypo)
  )

  final_table <- rbind(final_table, row_data)

  rm(meth_tmp, diff_tmp, dmc_all, dmc_hyper, dmc_hypo)
  gc()
}

print(final_table)

cat("Finished\n\n")

