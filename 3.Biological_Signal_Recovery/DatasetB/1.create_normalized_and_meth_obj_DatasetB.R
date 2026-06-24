# 09/06/2026
# In this script we generate the normalized objects of the STD15x, STD30x, CTM15x and CTM30x datasets so as to do the Fisrt
# validation (crossed comparison)
# ONLY AT AND EM SAMPLES
# DATASET 2

library('methylKit')

##### 1. CTM15x

source("/data2/gonzalo2/alvaro_data2/EMSeq_cov/CovGain.R")
dir_path = "/data2/gonzalo2/alvaro_data2/EMSeq_cov/External_Validation/15x_CpG.reports/"

metadata <- read.table(paste0(dir_path, "Sample_codes_external_validation_coverage.txt"),
                       header = TRUE, stringsAsFactors = FALSE)
file_list <- as.list(paste0(dir_path, metadata$SampleID, "_CpG_report.txt"))
sample_ids <- as.list(metadata$SampleID)
treatment_vector <- ifelse(metadata$Population == "AT", 0, 1)

myobj_CTM15x = methRead(file_list,
                                  sample.id = sample_ids,
                                  assembly = "seabass_V1.0",
                                  treatment = treatment_vector,
                                  mincov=1,
				  pipeline='bismarkCytosineReport',
                                  context = "CpG")

my_filtered_data <- filterByCombinedCpG(myobj_CTM15x, min.combined.cov = 5)

cat("Normalizing\n\n")
CTM15x_normalized=normalizeCoverage(my_filtered_data, method = "median")
rm(myobj_CTM15x, my_filtered_data)
gc()

CTM15x_meth=unite(CTM15x_normalized, min.per.group=9L, destrand=TRUE, mc.cores = 2)

##### 2. STD30x

dir_path = "/data2/gonzalo2/alvaro_data2/EMSeq_cov/External_Validation/30x_CpG.reports/"
metadata <- read.table(paste0(dir_path, "Sample_codes_external_validation_coverage.txt"),
                       header = TRUE, stringsAsFactors = FALSE)
file_list <- as.list(paste0(dir_path, metadata$SampleID, "_CpG_report.txt"))
sample_ids <- as.list(metadata$SampleID)
treatment_vector <- ifelse(metadata$Population == "AT", 0, 1)

myobj_STD30x <- methRead(file_list,
                  sample.id = sample_ids,
                  assembly = "seabass_V1.0",
                  pipeline = "bismarkCytosineReport",
                  treatment = treatment_vector,
                  context = "CpG",
                  mincov = 5)

filtered.myobj_STD30x=filterByCoverage(myobj_STD30x,lo.count=5,lo.perc=NULL,
                                hi.count=NULL,hi.perc=99.9)


STD30x_normalized=normalizeCoverage(filtered.myobj_STD30x, method = "median")

STD30x_meth=unite(STD30x_normalized, min.per.group=9L, destrand=TRUE, mc.cores = 2)

rm(myobj_STD30x, filtered.myobj_STD30x); gc()

##### 3. STD15x

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

STD15x_normalized=normalizeCoverage(filtered.myobj_STD15x, method = "median")

rm(myobj_STD15x, filtered.myobj_STD15x)
gc()

STD15x_meth=unite(STD15x_normalized, min.per.group=9L, destrand=TRUE, mc.cores = 2)

##### 4. CTM30x

dir_path = "/data2/gonzalo2/alvaro_data2/EMSeq_cov/External_Validation/30x_CpG.reports/"
metadata <- read.table(paste0(dir_path, "Sample_codes_external_validation_coverage.txt"),
                       header = TRUE, stringsAsFactors = FALSE)
file_list <- as.list(paste0(dir_path, metadata$SampleID, "_CpG_report.txt"))
sample_ids <- as.list(metadata$SampleID)
treatment_vector <- ifelse(metadata$Population == "AT", 0, 1)

myobj_CTM30x = methRead(file_list,
                                  sample.id = sample_ids,
                                  assembly = "seabass_V1.0",
                                  treatment = treatment_vector,
                                  mincov=1,
				  pipeline='bismarkCytosineReport',
                                  context = "CpG")

my_filtered_data <- filterByCombinedCpG(myobj_CTM30x, min.combined.cov = 5)

CTM30x_normalized=normalizeCoverage(my_filtered_data, method = "median")
rm(myobj_CTM30x, my_filtered_data)
gc()

CTM30x_meth=unite(CTM30x_normalized, min.per.group=9L, destrand=TRUE, mc.cores = 2)

save(CTM15x_normalized, CTM30x_normalized, STD15x_normalized, STD30x_normalized, file = "norm_objs_min5.Rdata")
save(CTM15x_meth, STD30x_meth, STD15x_meth, CTM30x_meth, file = "meth_objs_min5.Rdata")

message("Finnished!")
