# 08/06/2026
# AGR
# 15x script using CovGain
# IMPORTANT: use methRead(..., mincov = 1)

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
                 sample.id = list("L2-AT", "L3-AT", "L4-AT",
                                  "L9-EM", "L10-EM", "L11-EM","L12-EM"),
                 assembly = "seabass_V1.0",
                 treatment = c(0,0,0,1,1,1,1),
                 mincov=1,
                 pipeline = "bismarkCytosineReport",
                 context = "CpG")

my_filtered_data <- filterByCombinedCpG(myobj, min.combined.cov = 5, hi.perc = 99.9)

normalized=normalizeCoverage(my_filtered_data)

getFilterStats(normalized)

meth=unite(normalized, destrand=TRUE, mc.cores = 4)

meth
