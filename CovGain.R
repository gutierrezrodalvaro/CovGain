#' @title Filter methylRawList by combined CpG strand coverage
#'
#' @description
#' Designed for methylKit objects. Aligns +/- strands, sums coverage,
#' and filters based on the combined physical CpG site strength.
#'
#' @param mList A methylRawList object (CpG context!).
#' @param min.combined.cov The threshold for the summed strands (default 5).
#' @param hi.perc A numeric percentile (0-100) used to filter out extremely
#'   high-coverage bases (default: 99.9).
#'
#' @details
#' This function assumes CpG symmetry. It aligns '-' strand cytosines
#' to the '+' strand by subtracting 1 from the start position.
#' This is not applicable to CHH or CHG contexts.
#' It applies a double filter:
#' a minimum absolute coverage and a high-coverage percentile threshold to
#' remove potential PCR bias or mapping artifacts.
#'
#' @note IMPORTANT: Ensure that destrand = TRUE is used in subsequent
#' methylKit::unite() calls to finalize the merging process.
#'
#' @return A \code{methylRawList} object with filtered sites, retaining the
#'   original treatment and assembly metadata.
#'
#' @author Alvaro Gutierrez-Rodriguez
#' @importFrom data.table as.data.table setkey
#' @importFrom stats quantile
#' @importFrom methods new
#' @import methylKit
#' @export

filterByCombinedCpG <- function(mList, min.combined.cov = 5, hi.perc = 99.9) {
  if (!inherits(mList, "methylRawList")) stop("Input must be a methylRawList")

  results <- lapply(mList, function(obj) {
    dt <- data.table::as.data.table(methylKit::getData(obj))

    dt[, site_id := ifelse(strand == "-", start - 1L, start)]
    site_stats <- dt[, .(combined_cov = sum(coverage),
                         combined_cs  = sum(numCs)),
                     by = .(chr, site_id)]

    if (is.null(hi.perc)) {
      hi.count <- Inf
    } else {
      hi.count <- stats::quantile(site_stats$combined_cov, hi.perc / 100)
    }

    passing <- site_stats[combined_cov >= min.combined.cov & combined_cov < hi.count]

    data.table::setkey(dt, chr, site_id)
    data.table::setkey(passing, chr, site_id)
    keep_dt <- dt[passing, nomatch = 0]

    list(data = obj[dt[passing, which=TRUE, nomatch=0], ])
  })

  return(new("methylRawList", lapply(results, `[[`, "data"), treatment = mList@treatment))
}

#' @title Generate CovGain Summary Statistics
#'
#' @description
#' Calculates summary statistics based on the physical CpG site (combined strands).
#' It aligns reciprocal strands to compute the number of unique CpG sites and the
#' real combined coverage, reflecting the logic of the CovGain pipeline.
#' This stats may be evaluated after \code{methylKit::normalizeCoverage}
#'
#' @param mList A \code{methylRawList} object (typically already filtered).
#'
#' @return A data.frame with:
#'   \itemize{
#'     \item Total_Rows: Number of individual cytosines (strands).
#'     \item Unique_CpGs: Number of physical CpG sites (merged strands).
#'     \item Mean_Combined_Cov: Average coverage per physical CpG site.
#'     \item Global_Meth: Global methylation percentage.
#'   }
#'
#' @author Alvaro Gutierrez Rodriguez
#' @importFrom data.table as.data.table
#' @import methylKit
#' @export

getFilterStats <- function(mList) {
  if (!inherits(mList, "methylRawList")) stop("Input must be a methylRawList")

  stats_list <- lapply(mList, function(obj) {
    dt <- data.table::as.data.table(methylKit::getData(obj))

    dt[, site_id := ifelse(strand == "-", start - 1L, start)]

    site_stats <- dt[, .(combined_cov = sum(coverage),
                         combined_cs  = sum(numCs)),
                     by = .(chr, site_id)]

    stats_row <- data.frame(
      Sample = methylKit::getSampleID(obj),
      Total_Rows = nrow(dt),
      Unique_CpGs = nrow(site_stats),
      Mean_Combined_Cov = round(mean(site_stats$combined_cov), 2),
      Global_Meth = round(sum(site_stats$combined_cs) / sum(site_stats$combined_cov) * 100, 2),
      stringsAsFactors = FALSE
    )
    return(stats_row)
  })

  df_stats <- do.call(rbind, stats_list)

  cat("\n--- CovGain Summary Statistics (Combined Strands) ---\n")
  
  return(df_stats)
}
