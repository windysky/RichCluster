#' Filter MergedSeeds by Number of Terms
#'
#' This function filters clusters from the MergedSeeds result of RichCluster,
#' keeping only clusters with greater than or equal to min_terms # of terms.
#'
#' @param merged_seeds A dataframe containing the merged seeds with column named `ClusterIndices`.
#' @param min_terms An integer specifying the minimum number of terms required in a cluster.
#'
#' @return The filtered data frame with clusters filtered to include only those with at least `min_terms` terms.
#'
#' @export
filter_merged_seeds <- function(merged_seeds, min_terms)
{
  final_clusters <- merged_seeds %>%
    mutate(row_id = row_number()) %>%  # Add a row identifier
    separate_rows(TermIndices, sep = ", ") %>%  # Separate into individual rows
    group_by(row_id, Cluster) %>%  # Group by the original rows
    dplyr::filter(n() >= min_terms) %>%  # Filter groups with at least X terms
    summarise(TermIndices = paste(TermIndices, collapse = ", ")) %>%  # Collapse back to single strings
    ungroup() %>%  # Ungroup to finalize the data frame
    select(-row_id)  # Remove the temporary row identifier

  return(final_clusters)
}


#' Make Full Cluster Dataframe
#'
#' Use the final clusters and merged richsets to create a full cluster data frame
#' containing all Pvalue, Padj, and GeneID information for each term in each cluster.
#'
#' @param final_clusters The final clusters data frame containing a 'ClusterIndices' column with
#'                       comma-delimited indices corresponding to terms in the merged_richsets df.
#' @param merged_richsets The merged_richsets data frame containing the original
#'                        combined Pvalue, Padj, and GeneID information (from merge_richsets function)
#'
#' @return A reconstructed full_cluster dataframe, adds a 'Cluster' column to the merged_richsets
#'
#' @export
make_full_clusterdf <- function(final_clusters, merged_richsets) {
  # Initialize an empty data frame to store the results
  full_clusterdf <- data.frame()

  # Loop over each row in final_clusters
  for(i in seq_len(nrow(final_clusters))) {
    row <- final_clusters[i, ]
    TermIndices <- unlist(strsplit(row$TermIndices, ", "))

    # Loop over each term index in TermIndices
    for (termIndex in TermIndices) {
      R_termIndex <- as.integer(termIndex) + 1  # Convert termIndex to integer and adjust for 1-based indexing
      term_row <- merged_richsets[R_termIndex, ]  # Get the row corresponding to the termIndex

      # Create a new row with the cluster number and term row
      new_row <- c(Cluster = i, term_row)

      # Append the new row to the data frame
      full_clusterdf <- rbind(full_clusterdf, new_row)
    }
  }

  return(full_clusterdf)
}
