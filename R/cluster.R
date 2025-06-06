

#' Cluster Terms from Enrichment Results
#'
#' This function performs clustering on enrichment results by integrating
#' gene similarity scores and various clustering strategies.
#'
#' @param enrichment_results A list of dataframes, each containing enrichment results.
#'        Each dataframe should include at least the columns 'Term', 'GeneID', and 'Padj'.
#' @param df_names Optional, a character vector of names for the enrichment result dataframes. Must
#'        match the length of `enrichment_results`. Default is `NULL`.
#' @param min_terms Minimum number of terms each final cluster must include
#' @param min_value Minimum 'Pvalue' a term must have in order to be counted in final clustering
#' @param distance_metric A string specifying the distance metric to use (e.g., "kappa").
#' @param distance_cutoff A numeric value for the distance cutoff (0 < cutoff <= 1).
#' @param merge_strategy A string specifying the merge strategy to use (e.g., "DAVID").
#' @param membership_cutoff A numeric value between 0 and 1 for the membership cutoff.
#'
#' @return A named list containing:
#'         - `distance_matrix`: The distance matrix used in clustering.
#'         - `clusters`: The final clusters.
#'         - `df_list`: The original list of enrichment result dataframes.
#'         - `merged_df`: The merged dataframe containing combined results.
#'         - `cluster_options`: A list of clustering parameters used in the analysis.
#'         - `df_names` (optional): The names of the input dataframes if provided.
#'
#' @examples
#' enrichment_results <- list(
#'   df1 = data.frame(Term = c("A", "B"), GeneID = c("gene1", "gene2"), Padj = c(0.01, 0.02)),
#'   df2 = data.frame(Term = c("A", "C"), GeneID = c("gene3", "gene4"), Padj = c(0.03, 0.04))
#' )
#' result <- cluster(enrichment_results, distance_metric = "kappa", distance_cutoff = 0.5)
#' print(result$distance_matrix)
#' @export
cluster <- function(enrichment_results, df_names=NULL, min_terms=5, min_value=0.05,
                    distance_metric="kappa", distance_cutoff=0.5,
                    merge_strategy="DAVID", membership_cutoff=0.5) {

  if (is.null(df_names) || length(enrichment_results) != length(df_names)) {
    df_names <- as.character(seq_along(enrichment_results))
  }

  validate_inputs(enrichment_results, df_names, distance_metric, distance_cutoff,
                  merge_strategy, membership_cutoff)

  # accept a list of dataframes as input
  # call merge_enrichment_results
  merged_df <- merge_enrichment_results(enrichment_results)

  # call the cpp function
  term_vec <- merged_df$Term
  geneID_vec <- merged_df$GeneID
  padj_vec <- merged_df$Padj

  merged_df <- merged_df %>%
    filter(Pvalue < min_value) # as default, but user adjusts if they want

  # throw error if cluster options are invalid

  cluster_result <- RichCluster::RichCluster(
    distance_metric, distance_cutoff,
    merge_strategy, membership_cutoff,
    term_vec,
    geneID_vec,
    padj_vec
  )

  # add the original stuff to the cluster_result
  # (helps visualizations later)
  cluster_options <- list(
    distance_metric = distance_metric,
    distance_cutoff = distance_cutoff,
    merge_strategy = merge_strategy,
    membership_cutoff = membership_cutoff
  )

  cluster_result$df_list <- enrichment_results
  cluster_result$merged_df <- merged_df
  cluster_result$cluster_options <- cluster_options
  cluster_result$df_names <- df_names

  cluster_result$final_clusters <- filter_clusters(cluster_result$all_clusters, min_terms)
  cluster_result$cluster_df <- make_full_clusterdf(cluster_result$final_clusters, merged_df)

  return(cluster_result)
}


validate_inputs <- function(enrichment_results, df_names=NA_character_,
                            distance_metric="kappa", distance_cutoff=0.5,
                            merge_strategy="DAVID", membership_cutoff=0.5) {
  if (!is.list(enrichment_results)) {
    stop("enrichment_results must be a list of dataframes.")
  }
  if (any(!sapply(enrichment_results, is.data.frame))) {
    stop("Each element of enrichment_results must be a dataframe.")
  }
  if (distance_cutoff <= 0 || distance_cutoff > 1) {
    stop("distance_cutoff must be between 0 and 1.")
  }
  if (membership_cutoff <= 0 || membership_cutoff > 1) {
    stop("membership_cutoff must be between 0 and 1.")
  }
  if (distance_metric != "kappa" && distance_metric != "jaccard") {
    stop("Unsupported distance metric. Only 'kappa' and 'jaccard' are supported.")
  }
  if (merge_strategy != "DAVID") {
    stop("Unsupported merge_strategy. Only DAVID is supported.")
  }

}

#' Filter Clusters by Number of Terms
#'
#' Filters the full list of clusters by keeping only those with greater
#' than or equal to min_terms # of terms.
#'
#' @param full_clusters A dataframe containing the merged seeds with column named `ClusterIndices`.
#' @param min_terms An integer specifying the minimum number of terms required in a cluster.
#'
#' @return The filtered data frame with clusters filtered to include only those with at least `min_terms` terms.
#'
#' @export
filter_clusters <- function(all_clusters, min_terms)
{
  filtered_clusters <- all_clusters %>%
    mutate(row_id = row_number()) %>%  # Add a row identifier
    separate_rows(TermIndices, sep = ", ") %>%  # Separate into individual rows
    group_by(row_id, Cluster) %>%  # Group by the original rows
    dplyr::filter(n() >= min_terms) %>%  # Filter groups with at least X terms
    summarise(TermIndices = paste(TermIndices, collapse = ", ")) %>%  # Collapse back to single strings
    ungroup() %>%  # Ungroup to finalize the data frame
    select(-row_id)  # Remove the temporary row identifier

  return(filtered_clusters)
}


make_full_clusterdf <- function(final_clusters, merged_df) {
  # Initialize an empty data frame to store the results
  full_clusterdf <- data.frame()

  # Loop over each row in final_clusters
  for(i in seq_len(nrow(final_clusters))) {
    row <- final_clusters[i, ]
    TermIndices <- unlist(strsplit(row$TermIndices, ", "))

    # Loop over each term index in TermIndices
    for (termIndex in TermIndices) {
      R_termIndex <- as.integer(termIndex) + 1  # Convert termIndex to integer and adjust for 1-based indexing
      term_row <- merged_df[R_termIndex, ]  # Get the row corresponding to the termIndex

      # Create a new row with the cluster number and term row
      new_row <- c(Cluster = i, term_row)

      # Append the new row to the data frame
      full_clusterdf <- rbind(full_clusterdf, new_row)
    }
  }

  return(full_clusterdf)
}

