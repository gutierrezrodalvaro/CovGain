# Anotation DatasetB
# Genomic features and CpG features
# S5_15, C5_15 & C5_30

library(methylKit)
library(dplyr)
library(purrr)
library(txdbmaker)
library(rtracklayer)

load('DMCs_F5.Rdata')

datasets_list <- list(
  "CTM 15x" = C15_diff,
  "CTM 30x" = C30_diff,
  "STD 30x" = S30_diff
)

# Prepare data for genomic annotation
txdb <- makeTxDbFromGFF("/data2/gonzalo2/alvaro_data2/EMSeq_cov/genomic_ref_data/NEW_COMBINED_ANNOTATION_FUNCTION.gtf.txt", format = "gtf")
all_exons <- exons(txdb, columns = "exon_rank")
exon_ranks <- sapply(mcols(all_exons)$exon_rank, function(x) x[1])

gene.obj <- list(
  "Promoter"      = promoters(txdb, upstream = 2000, downstream = 200),
  "First Exon"    = all_exons[exon_ranks == 1],
  "Rest of Exons" = all_exons[exon_ranks > 1],
  "Intron"        = unlist(intronsByTranscript(txdb))
)

cpg_bed <- rtracklayer::import("/data2/gonzalo2/alvaro_data2/EMSeq_cov/genomic_ref_data/final_cpg_annotation.bed")
cpg_features <- list(
  "CpG Island" = cpg_bed[cpg_bed$name == "Island"],
  "CpG Shore"  = cpg_bed[cpg_bed$name == "Shore"],
  "CpG Shelf"  = cpg_bed[cpg_bed$name == "Shelf"],
  "Open Sea"   = cpg_bed[cpg_bed$name == "OpenSea"]
)

# Counts only the first category a position is in
count_annotations_final <- function(meth_obj, label, g_feat, c_feat) {
  gr <- as(meth_obj, "GRanges")
  total_sites <- length(gr)

  process_set <- function(feat_list, set_name) {
    counts_df <- map_df(names(feat_list), function(feat_name) {
      # Sitios que solapan con esta categoría
      hits <- countOverlaps(gr, feat_list[[feat_name]]) > 0
      data.frame(Region = feat_name, n = sum(hits))
    })

    sum_n <- sum(counts_df$n)
    gap <- total_sites - sum_n

    if(set_name == "Genomic") {
      counts_df <- bind_rows(counts_df, data.frame(Region = "Intergenic", n = gap))
    } else {
      if(gap > 0) counts_df <- bind_rows(counts_df, data.frame(Region = "Other/Unknown", n = gap))
    }

    counts_df$Set <- set_name
    counts_df$Dataset <- label
    return(counts_df)
  }

  bind_rows(
    process_set(g_feat, "Genomic"),
    process_set(c_feat, "CpG Context")
  )
}

# Results
final_counts_table <- imap_dfr(datasets_list, ~count_annotations_final(.x, .y, gene.obj, cpg_features))

final_summary <- final_counts_table %>%
  group_by(Dataset, Set) %>%
  mutate(
    Percentage = (n / sum(n)) * 100,
    Label_Full = paste0("n = ", format(n, big.mark = ","), " (", round(Percentage, 1), "%)")
  ) %>%
  ungroup()

final_annot_counts_datasetB <- as.data.frame(final_summary)
write.csv(final_annot_counts_datasetB, 
          file = "final_annot_counts_datasetB.csv", 
          row.names = FALSE)



print("Tabla de Resumen:")
print(final_summary)
