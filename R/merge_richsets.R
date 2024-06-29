library(dplyr)
library(tidyverse)

#' Merge List of Enrichment Results
#'
#' This function merges multiple enrichment results ('richsets') into a single richset by
#' combining unique GeneID elements across each richset, and averaging Pvalue / Padj
#' values for each term across all richsets.
#'
#' @param richsets A list of richset dataframes containing columns c('Term', 'GeneID', 'Pvalue', 'Padj')
#'
#' @return A single merged richset dataframe with all original columns
#'         suffixed with the index of the richset, with new columns 'GeneID', 'Pvalue',
#'         'Padj' containing the merged values.
#'
#' @export
merge_richsets <- function(richsets) {
  # Note: Keep track of what index each richset has in the list

  COLNAME_SEP <- "_"

  # Preprocessing: Suffix all non 'Term' columns by their index in the list
  # Allows base::merge by 'Term'
  for (i in seq_along(richsets)) {
    rownames(richsets[[i]]) <- NULL # Prevents rownames from causing errors

    # Get renamed columns (suffixed by index)
    all_colnames <- colnames(richsets[[i]])
    nonterm_cols <- all_colnames[all_colnames != 'Term']
    all_colnames[all_colnames != 'Term'] <- paste(nonterm_cols, i, sep=COLNAME_SEP)

    colnames(richsets[[i]]) <- all_colnames
  }

  # Initialize merged_gs with first richset
  merged_gs <- richsets[[1]]
  if (length(richsets) == 1) {
    return(merged_gs) # Return the single richset if only one is provided
  }

  # Else, merge the rest of the richsets
  for (i in 2:length(richsets)) {
    merged_gs <- base::merge(merged_gs, richsets[[i]], by='Term', all=TRUE)
  }

  # For each row in merged_gs, combine unique GeneID elements
  geneid_cols <- paste("GeneID", seq_along(richsets), sep=COLNAME_SEP)
  merged_gs$GeneID <- apply(merged_gs[, geneid_cols], 1, function(x) {
    paste(unique(na.omit(x)), collapse = ',')
  })

  # Average the value columns across all richsets
  # Avg Pvalue
  pvalue_cols <- paste("Pvalue", seq_along(richsets), sep=COLNAME_SEP)
  merged_gs$Pvalue <- rowMeans(merged_gs[, pvalue_cols], na.rm=TRUE)

  # Avg Padj
  padj_cols <- paste("Padj", seq_along(richsets), sep=COLNAME_SEP)
  merged_gs$Padj <- rowMeans(merged_gs[, padj_cols], na.rm=TRUE)

  # Return the merged richset df
  return(merged_gs)

}
