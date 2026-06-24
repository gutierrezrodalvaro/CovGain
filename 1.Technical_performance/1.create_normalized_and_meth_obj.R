# 08/06/2026
# In this script we generate the normalized objects of the STD15x, STD30x, CTM15x and CTM30x datasets so as to do the Fisrt
# validation (crossed comparison)
# ONLY AT AND EM SAMPLES
# We also generate the meth objects, with no NA, to do the second part of the first validation

library('methylKit')
source("/data2/gonzalo2/alvaro_data2/EMSeq_cov/CovGain.R")

##### 1. CTM15x

file_path_CTM15x = "/data2/gonzalo2/alvaro_data2/EMSeq_cov/15x/meth_extracted_15x/CpG_reports/"

file_list_CTM15x = list(
  paste0(file_path_CTM15x, "L2-AT.CpG_report.txt"),
  paste0(file_path_CTM15x, "L3-AT.CpG_report.txt"),
  paste0(file_path_CTM15x, "L4-AT.CpG_report.txt"),
  paste0(file_path_CTM15x, "L9-EM.CpG_report.txt"),
  paste0(file_path_CTM15x, "L10-EM.CpG_report.txt"),
  paste0(file_path_CTM15x, "L11-EM.CpG_report.txt"),
  paste0(file_path_CTM15x, "L12-EM.CpG_report.txt")
)

myobj_CTM15x = methRead(file_list_CTM15x,
                 sample.id = list("L2-AT-CTM15x", "L3-AT-CTM15x", "L4-AT-CTM15x",
                                                   "L9-EM-CTM15x", "L10-EM-CTM15x", "L11-EM-CTM15x","L12-EM-CTM15x"),
                 assembly = "seabass_V1.0",
                 treatment = c(0,0,0,1,1,1,1),
                 mincov=1,
                 pipeline = "bismarkCytosineReport",
                 context = "CpG")
                  
my_filtered_data_CTM15x <- filterByCombinedCpG(myobj_CTM15x, min.combined.cov = 5,hi.perc=99.9)

CTM15x_normalized=normalizeCoverage(my_filtered_data_CTM15x)

rm(myobj_CTM15x, my_filtered_data_CTM15x)
gc()

CTM15x_meth=unite(CTM15x_normalized, destrand=TRUE, mc.cores = 2)

##### 2. STD30x

file_path_STD30x="/data2/gonzalo2/alvaro_data2/EMSeq_cov/30x/raw_methylation_data/"

file_list_STD30x = list(
  paste0(file_path_STD30x, "L2-AT.CpG_report.txt"),
  paste0(file_path_STD30x, "L3-AT.CpG_report.txt"),
  paste0(file_path_STD30x, "L4-AT.CpG_report.txt"),
  paste0(file_path_STD30x, "L9-EM.CpG_report.txt"),
  paste0(file_path_STD30x, "L10-EM.CpG_report.txt"),
  paste0(file_path_STD30x, "L11-EM.CpG_report.txt"),
  paste0(file_path_STD30x, "L12-EM.CpG_report.txt")
)

myobj_STD30x = methRead(file_list_STD30x,
                 sample.id = list("L2-AT-STD30x", "L3-AT-STD30x", "L4-AT-STD30x",
                                  "L9-EM-STD30x", "L10-EM-STD30x", "L11-EM-STD30x","L12-EM-STD30x"),
                 assembly = "seabass_V1.0",
                 treatment = c(0,0,0,1,1,1,1),
                 mincov=5,
                 pipeline = "bismarkCytosineReport",
                 context = "CpG")

filtered.myobj_STD30x=filterByCoverage(myobj_STD30x,lo.count=5,lo.perc=NULL, hi.count=NULL,hi.perc=99.9)
STD30x_normalized=normalizeCoverage(filtered.myobj_STD30x)

rm(myobj_STD30x, filtered.myobj_STD30x)
gc()

STD30x_meth=unite(STD30x_normalized, destrand=TRUE, mc.cores = 2)

##### 3. STD15x

file_path_STD15x = "/data2/gonzalo2/alvaro_data2/EMSeq_cov/15x/meth_extracted_15x/CpG_reports/"
file_list_STD15x = list(
  paste0(file_path_STD15x, "L2-AT.CpG_report.txt"),
  paste0(file_path_STD15x, "L3-AT.CpG_report.txt"),
  paste0(file_path_STD15x, "L4-AT.CpG_report.txt"),
  paste0(file_path_STD15x, "L9-EM.CpG_report.txt"),
  paste0(file_path_STD15x, "L10-EM.CpG_report.txt"),
  paste0(file_path_STD15x, "L11-EM.CpG_report.txt"),
  paste0(file_path_STD15x, "L12-EM.CpG_report.txt")
)

myobj_STD15x = methRead(file_list_STD15x,
                 sample.id = list("L2-AT-STD15x", "L3-AT-STD15x", "L4-AT-STD15x",
                                  "L9-EM-STD15x", "L10-EM-STD15x", "L11-EM-STD15x","L12-EM-STD15x"),
                 assembly = "seabass_V1.0",
                 treatment = c(0,0,0,1,1,1,1),
                 mincov=5,
                 pipeline = "bismarkCytosineReport",
                 context = "CpG")

filtered.myobj_STD15x=filterByCoverage(myobj_STD15x,lo.count=5,lo.perc=NULL, hi.count=NULL,hi.perc=99.9)
STD15x_normalized=normalizeCoverage(filtered.myobj_STD15x)

rm(myobj_STD15x, filtered.myobj_STD15x)
gc()

STD15x_meth=unite(STD15x_normalized, destrand=TRUE, mc.cores = 2)

##### 4. CTM30x

file_path_CTM30x = "/data2/gonzalo2/alvaro_data2/EMSeq_cov/30x/raw_methylation_data/"
file_list_CTM30x = list(
  paste0(file_path_CTM30x, "L2-AT.CpG_report.txt"),
  paste0(file_path_CTM30x, "L3-AT.CpG_report.txt"),
  paste0(file_path_CTM30x, "L4-AT.CpG_report.txt"),
  paste0(file_path_CTM30x, "L9-EM.CpG_report.txt"),
  paste0(file_path_CTM30x, "L10-EM.CpG_report.txt"),
  paste0(file_path_CTM30x, "L11-EM.CpG_report.txt"),
  paste0(file_path_CTM30x, "L12-EM.CpG_report.txt")
)

myobj_CTM30x = methRead(file_list_CTM30x,
                 sample.id = list("L2-AT-CTM30x", "L3-AT-CTM30x", "L4-AT-CTM30x",
                                  "L9-EM-CTM30x", "L10-EM-CTM30x", "L11-EM-CTM30x","L12-EM-CTM30x"),
                 assembly = "seabass_V1.0",
                 treatment = c(0,0,0,1,1,1,1),
                 mincov=1,
                 pipeline = "bismarkCytosineReport",
                 context = "CpG")

my_filtered_data_CTM30x <- filterByCombinedCpG(myobj_CTM30x, min.combined.cov = 5, hi.perc = 99.9)

CTM30x_normalized=normalizeCoverage(my_filtered_data_CTM30x)

rm(myobj_CTM30x, my_filtered_data_CTM30x)
gc()

CTM30x_meth=unite(CTM30x_normalized, destrand=TRUE, mc.cores = 2)

save(CTM15x_normalized, CTM30x_normalized, STD15x_normalized, STD30x_normalized, file = "norm_objs_min5.Rdata")
save(CTM15x_meth, STD30x_meth, STD15x_meth, CTM30x_meth, file = "meth_objs_min5.Rdata")

message("Finnished!")
