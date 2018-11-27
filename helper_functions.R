require("ggplot2")
require("Rtsne")
require("gridExtra")
require("griph")
require("org.Hs.eg.db")
library("SingleCellExperiment")
library("scater")
library("monocle")
library("destiny")
library("ggbeeswarm")
library("scran")


#' This is the function for selection of overdispersed genes adapted from:
#' https://github.com/10XGenomics/single-cell-3prime-paper/
#'
#' @param m  a (protentially sparse) gene x cells count matrix
#' @param f  a number between 0 and 1 the fraction of overdispersed genes to keep
#' @return a vector of the indices of genes to keep.
select_variable_genes<-function(m,f) {
  df <- data.frame(mean = rowMeans(m + 1/ncol(m)), cv = apply(m,1,sd) / rowMeans(m + 1/ncol(m)), var = apply(m,1,var))
  df$dispersion <- with(df, var/mean)
  df$mean_bin <- with(df, cut(mean, breaks = c(-Inf, unique(quantile(mean, seq(0.1,1,0.02), na.rm = TRUE)), Inf)))
  var_by_bin <- data.frame(mean_bin = factor(levels(df$mean_bin), levels = levels(df$mean_bin)),
                           bin_median = as.numeric(tapply(df$dispersion, df$mean_bin, stats::median)),
                           bin_mad = as.numeric(tapply(df$dispersion, df$mean_bin, stats::mad)))[table(df$mean_bin) > 0,]
  df$bin_disp_median <- var_by_bin$bin_median[match(df$mean_bin, var_by_bin$mean_bin)]
  df$bin_disp_mad <- var_by_bin$bin_mad[match(df$mean_bin, var_by_bin$mean_bin)]
  df$dispersion_norm <- with(df, (dispersion - bin_disp_median)/(bin_disp_mad + 0.01) )

  n_genes_keep=ceiling(f*nrow(m) ) #In the end retain only the top 100*f% overdispersed genes
  disp_cut_off <- sort(df$dispersion_norm,decreasing=TRUE)[n_genes_keep]
  genes_keep <- which(df$dispersion_norm >= disp_cut_off)
  return(genes_keep)
}