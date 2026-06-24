# 08/06/2026
# This script calculates de number of CpGs between AT and EM populations
# It uses de 15x datasets with the S method

library(methylKit)

file_path = "/data2/gonzalo2/alvaro_data2/EMSeq_cov/15x/meth_extracted_15x/CpG_reports/"
# AT vs EM 15x
file_list_1 = list(
  paste0(file_path, "L2-AT.CpG_report.txt"),
  paste0(file_path, "L3-AT.CpG_report.txt"),
  paste0(file_path, "L4-AT.CpG_report.txt"),
  paste0(file_path, "L9-EM.CpG_report.txt"),
  paste0(file_path, "L10-EM.CpG_report.txt"),
  paste0(file_path, "L11-EM.CpG_report.txt"),
  paste0(file_path, "L12-EM.CpG_report.txt")
)

myobj_STD15x_1 <- methRead(
    location   = as.list(file_list_1),
    sample.id  = list("L2-AT", "L3-AT", "L4-AT", "L9-EM", "L10-EM", "L11-EM", "L12-EM"),
    assembly   = "seabass_V1.0",
    treatment  = c(0,0,0,1,1,1,1),
    context    = "CpG",
    header     = FALSE,
    pipeline = "bismarkCytosineReport",
    mincov     = 5
)

filtered.myobj_1=filterByCoverage(myobj_STD15x_1,lo.count=5,lo.perc=NULL,
                                hi.count=NULL,hi.perc=99.9)

STD15x_normalized.obj_1=normalizeCoverage(filtered.myobj_1)

### DMCs!!

STD15x_meth_1=unite(STD15x_normalized.obj_1, destrand=TRUE, mc.cores = 4)
cat("The number of common CpGs between the AT and the EM at 15x is:", nrow(STD15x_meth_1), "\n")

# Diff meth DMCs
myDiff_1=calculateDiffMeth(STD15x_meth_1)

myDiff15p_1=getMethylDiff(myDiff_1,difference=15,qvalue=0.05)
cat("AT vs. EM DMCs:", nrow(myDiff15p_1), "\n")
myDiff15p.hypo_1=getMethylDiff(myDiff_1,difference=15,qvalue=0.05,type="hypo")
cat("AT vs. EM Hypo DMCs:", nrow(myDiff15p.hypo_1), "\n")
myDiff15p.hyper_1=getMethylDiff(myDiff_1,difference=15,qvalue=0.05,type="hyper")
cat("AT vs. EM Hyper DMCs:", nrow(myDiff15p.hyper_1), "\n")

message('Finnished')
