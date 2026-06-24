# 08/06/2026
# AGR
# 10x script methylKit Standard
# IMPORTANT: use methRead(..., mincov = 5)

library(methylKit)

source("/data2/gonzalo2/alvaro_data2/EMSeq_cov/CovGain.R")

file_path = "/data2/gonzalo2/alvaro_data2/EMSeq_cov/10x/meth_extracted_10x/CpG_reports/"
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
                 sample.id = list("L2-AT", "L3-AT", "L4-AT",
                                  "L9-EM", "L10-EM", "L11-EM","L12-EM"),
                 assembly = "seabass_V1.0",
                 treatment = c(0,0,0,1,1,1,1),
                 mincov=5,
                 pipeline = "bismarkCytosineReport",
                 context = "CpG")

filtered.myobj=filterByCoverage(myobj,lo.count=5,lo.perc=NULL,
                                hi.count=NULL,hi.perc=99.9)

normalized.myobj=normalizeCoverage(filtered.myobj, method = "median")

getFilterStats(normalized.myobj)
meth=unite(normalized.myobj, destrand=T)

meth
