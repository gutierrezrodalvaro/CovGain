# 08/06/2026
# AGR
# With this script we do the part one of the second validation
# Here we compare the two populations using the CTM15x dataset
# Data for figure supplementary 1

library('methylKit')

source("/data2/gonzalo2/alvaro_data2/EMSeq_cov/CovGain.R")

file_path = "/data2/gonzalo2/alvaro_data2/EMSeq_cov/15x/meth_extracted_15x/CpG_reports/"
file_list = list(
  paste0(file_path, "L2-AT.CpG_report.txt"),
  paste0(file_path, "L3-AT.CpG_report.txt"),
  paste0(file_path, "L4-AT.CpG_report.txt"),
  paste0(file_path, "L9-EM.CpG_report.txt"),
  paste0(file_path, "L10-EM.CpG_report.txt"),
  paste0(file_path, "L11-EM.CpG_report.txt"),
  paste0(file_path, "L12-EM.CpG_report.txt")
)

myobj = methRead(file_list,
                 sample.id = list("L2-AT-CTM15x", "L3-AT-CTM15x", "L4-AT-CTM15x",
                   		  "L9-EM-CTM15x", "L10-EM-CTM15x", "L11-EM-CTM15x","L12-EM-CTM15x"),
                 assembly = "seabass_V1.0",
                 treatment = c(0,0,0,1,1,1,1),
                 mincov=1,
                 pipeline = "bismarkCytosineReport",
                 context = "CpG")

my_filtered_data_1 <- filterByCombinedCpG(myobj, min.combined.cov = 5)

CTM15x_normalized.obj_1=normalizeCoverage(my_filtered_data_1)

######### DMC!!
CTM15x_meth_1=unite(CTM15x_normalized.obj_1, destrand=TRUE, mc.cores = 4)
rm(myobj_CTM15x_1, filtered.myobj_1, my_filtered_data_1); gc()
cat("The number of common CpGs between the AT and the EM at 15x is:", nrow(CTM15x_meth_1), "\n")

myDiff_1=calculateDiffMeth(CTM15x_meth_1, mc.cores = 4)

myDiff15p_1=getMethylDiff(myDiff_1,difference=15,qvalue=0.05)
cat("AT vs. EM DMCs:", nrow(myDiff15p_1), "\n")
myDiff15p.hypo_1=getMethylDiff(myDiff_1,difference=15,qvalue=0.05,type="hypo")
cat("AT vs. EM Hypo DMCs:", nrow(myDiff15p.hypo_1), "\n")
myDiff15p.hyper_1=getMethylDiff(myDiff_1,difference=15,qvalue=0.05,type="hyper")
cat("AT vs. EM Hyper DMCs:", nrow(myDiff15p.hyper_1), "\n")

