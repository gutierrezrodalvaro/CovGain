# This script does the Figure 3-Scatter plots

library(methylKit)
library(ggplot2)
library(patchwork)
library(dplyr)

load('DMCs_F5.Rdata') # "C15_diff" "C15_DMCs" "C30_diff" "C30_DMCs" "S30_diff" "S30_DMCs"

# Custom theme with precise text controls and vertical Y-axis title
theme_scatter <- function() {
  theme_minimal() %+replace%
    theme(panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black"),      
          text             = element_text(size = 14),     
	  axis.title.x     = element_text(size = 14, margin = margin(t = 8)), 
          axis.title.y     = element_text(size = 14, angle = 90, vjust = 0.5, margin = margin(r = 8)), 
          
          axis.text        = element_text(size = 12, color = "black"),        
          plot.title       = element_text(size = 14, face = "bold", hjust = 0.5),
          legend.text      = element_text(size = 14),                                  
          plot.margin = margin(t = 12, r = 15, b = 12, l = 15),          
          legend.position = "bottom",                     
          legend.background = element_blank(),
          title = element_blank())
}

# Helper function (Paired membership logic)
prep_comparison <- function(diff1, diff2, dmc1, dmc2, label_x, label_y) {
  df1 <- getData(diff1)
  df2 <- getData(diff2)

  df_merged <- merge(df1, df2, by = c("chr", "start"), suffixes = c(".x", ".y"))

  ids1 <- paste(getData(dmc1)$chr, getData(dmc1)$start, sep="_")
  ids2 <- paste(getData(dmc2)$chr, getData(dmc2)$start, sep="_")
  df_merged$id <- paste(df_merged$chr, df_merged$start, sep="_")

  df_merged <- df_merged %>%
    mutate(significance = case_when(
      id %in% ids1 & id %in% ids2 ~ "Both Sig.",
      id %in% ids1 & !(id %in% ids2) ~ paste("Sig.", label_x),
      !(id %in% ids1) & id %in% ids2 ~ paste("Sig.", label_y),
      TRUE ~ "Non-Sig."
    ))

  return(df_merged)
}

# Plotting function with strict Z-order and 2x2 layout
create_scatter <- function(df, title, xlab, ylab, colors) {

  mid_labels <- setdiff(names(colors), c("Non-Sig.", "Both Sig."))
  plot_order <- c("Non-Sig.", mid_labels, "Both Sig.")

  df_sorted <- df %>%
    mutate(significance = factor(significance, levels = plot_order)) %>%
    arrange(significance)

  ggplot(df_sorted, aes(x = meth.diff.x, y = meth.diff.y, color = significance)) +
    geom_point(alpha = 0.2, size = 0.4) +
    scale_color_manual(values = colors) +
    geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dashed") +
    labs(title = title, x = xlab, y = ylab, color = NULL) +
    theme_scatter() +
    guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 4, alpha = 1)))
}

# Comparison 1: S_min5_30x vs C_min5_30x
df_plot1 <- prep_comparison(S30_diff, C30_diff, S30_DMCs, C30_DMCs, "S_30x", "C_30x")
cols1 <- c(
  "Non-Sig."         = "gray80",
  "Sig. S_30x"  = "#d35400",
  "Sig. C_30x"  = "#2980b9",
  "Both Sig."        = "#6E8B3D"
)

# Comparison 2: C_min5_30x vs C_min5_15x
df_plot2 <- prep_comparison(C30_diff, C15_diff, C30_DMCs, C15_DMCs, "C_30x", "C_15x")
cols2 <- c(
  "Non-Sig."         = "gray80",
  "Sig. C_30x"  = "#d35400",
  "Sig. C_15x"  = "#2980b9",
  "Both Sig."        = "#6E8B3D"
)

# GENERATE PLOTS

p1 <- create_scatter(df_plot1, "C_30x vs S_30x", "S_30x Meth. Diff", "C_30x Meth. Diff", cols1)
p2 <- create_scatter(df_plot2, "C_15x vs C_30x", "C_30x Meth. Diff", "C_15x Meth. Diff", cols2)

combined_scatter <- (p1 | p2) +
  plot_layout(guides = "keep", widths = c(1, 1))

ggsave("Figure2_DMC_Scatters_Final_DatasetB.jpg", combined_scatter, width = 9.5, height = 5.2, dpi = 600)

# CALCULATE PEARSON CORRELATION

cor_p1 <- cor(df_plot1$meth.diff.x, df_plot1$meth.diff.y, method = "pearson")
cor_p2 <- cor(df_plot2$meth.diff.x, df_plot2$meth.diff.y, method = "pearson")

print(paste("S_min5_30x vs C_min5_30x correlation:", round(cor_p1, 4)))
print(paste("C_min5_30x vs C_min5_15x correlation:", round(cor_p2, 4)))
