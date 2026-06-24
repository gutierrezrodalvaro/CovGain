# 26/03/2026
# AGR
# With this script we do the part one of the second validation
# Here we compare the three populations using the CTM30x dataset

library('methylKit')

source("/data2/gonzalo2/alvaro_data2/EMSeq_cov/30x/raw_methylation_data/readBismarkFiles.R")

# Load the files
file_path = "/data2/gonzalo2/alvaro_data2/EMSeq_cov/30x/raw_methylation_data/"

# 1. AT vs. EM
file_list_1 = list(
  paste0(file_path, "L2-AT.CpG_report.txt"),
  paste0(file_path, "L3-AT.CpG_report.txt"),
  paste0(file_path, "L4-AT.CpG_report.txt"),
  paste0(file_path, "L9-EM.CpG_report.txt"),
  paste0(file_path, "L10-EM.CpG_report.txt"),
  paste0(file_path, "L11-EM.CpG_report.txt"),
  paste0(file_path, "L12-EM.CpG_report.txt")
)

myobj_CTM30x_1 = methRead(file_list_1,
                                  sample.id = list("L2-AT-CTM30x", "L3-AT-CTM30x", "L4-AT-CTM30x",
                                                   "L9-EM-CTM30x", "L10-EM-CT30x", "L11-EM-CTM30x","L12-EM-CTM30x"),
                                  assembly = "seabass_V1.0",
                                  treatment = c(0,0,0,1,1,1,1),
                                  mincov=1,
				  pipeline='bismarkCytosineReport',
				  context = "CpG")

my_filtered_data_1 <- filterByCombinedCpG(myobj_CTM30x_1, min.combined.cov = 5)

CTM30x_normalized.obj_1=normalizeCoverage(my_filtered_data_1)

rm(myobj_CTM30x_1, my_filtered_data_1); gc()

######### DMC!!
CTM30x_meth_1=unite(CTM30x_normalized.obj_1, destrand=TRUE, mc.cores = 1)

cat("The number of common CpGs between the AT and the EM at 30x is:", nrow(CTM30x_meth_1), "\n")

myDiff_1=calculateDiffMeth(CTM30x_meth_1, mc.cores = 2)

myDiff15p_1=getMethylDiff(myDiff_1,difference=15,qvalue=0.05)
cat("AT vs. EM DMCs:", nrow(myDiff15p_1), "\n")
myDiff15p.hypo_1=getMethylDiff(myDiff_1,difference=15,qvalue=0.05,type="hypo")
cat("AT vs. EM Hypo DMCs:", nrow(myDiff15p.hypo_1), "\n")
myDiff15p.hyper_1=getMethylDiff(myDiff_1,difference=15,qvalue=0.05,type="hyper")
cat("AT vs. EM Hyper DMCs:", nrow(myDiff15p.hyper_1), "\n")
