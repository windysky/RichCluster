#' Preprocess the data for heatmap plotting
#'
#' @param data The full cluster data frame containing clusterNumber and value columns
#' @param value_type The type of value to use for heatmap plotting (e.g., "Padj")
#' @param summary_metric The summary metric to use for summarizing the value columns (e.g., mean)
make_summary_data <- function(full_clusterdf, value_type="Padj") {

  # Replace NA, NaN, and Inf values with 0
  replace_na_nan_inf <- function(x) {
    x[is.na(x) | is.nan(x) | is.infinite(x)] <- 0
    return(x)
  }

  # Group by clusterNumber and summarize the value columns, replacing NA, NaN, and Inf with 0
  summary_data <- full_clusterdf %>%
    group_by(Cluster) %>%
    summarise(across(starts_with(value_type), mean, na.rm = TRUE)) %>%
    mutate(across(where(is.numeric), replace_na_nan_inf))

  return(summary_data)
}


#' Create a Heatmap for All Clusters
#'
#' This function generates a heatmap for all clusters based on the specified value type (e.g., Padj or Pvalue).
#'
#' @param full_clusterdf A dataframe containing the full cluster data.
#' @param value_type A string specifying the type of value to use for the heatmap (default is "Padj").
#' @param x_labels A vector of strings specifying the labels for the x-axis (default is NULL, which uses the original column names).
#' @param title A string specifying the title of the heatmap (default is NULL, which generates a default title).
#' @return An interactive heatmaply heatmap.
#' @export
all_clusters_hmap <- function(full_clusterdf, value_type = "Padj", x_labels = NULL, title = NULL) {
  # Preprocess the full_clusterdf
  summary_data <- make_summary_data(full_clusterdf, value_type)
  print(nrow(summary_data))

  # Select columns for heatmap (excluding the Cluster and the average column)
  value_cols <- grep(paste0("^", value_type, "_\\d+"), colnames(summary_data), value = TRUE)
  heatmap_data <- summary_data %>%
    select(Cluster, all_of(value_cols))

  # Find the term with the lowest average value_type for each cluster
  full_clusterdf <- full_clusterdf %>%
    group_by(Cluster) %>%
    mutate(avg_value = rowMeans(across(starts_with(value_type)), na.rm = TRUE)) %>%
    filter(avg_value == min(avg_value, na.rm = TRUE)) %>%
    slice(1)  # In case of ties, take the first one

  # Use the term with the lowest average value_type as the x-axis label
  cluster_labels <- full_clusterdf %>%
    select(Cluster, Term) %>%
    distinct() %>%
    arrange(Cluster) %>%
    pull(Term)

  # Convert to matrix
  heatmap_matrix <- heatmap_data %>%
    select(-Cluster) %>%
    as.matrix()
  rownames(heatmap_matrix) <- cluster_labels  # Set custom x-axis labels

  # Set custom x-axis labels if provided
  if (!is.null(x_labels) && length(x_labels) == ncol(heatmap_matrix)) {
    colnames(heatmap_matrix) <- x_labels
  }
  # Set y-axis labels based on term with lowest avg padj
  rownames(heatmap_matrix) <- cluster_labels

  # Apply log10 transformation if required
  heatmap_matrix <- -log10(heatmap_matrix)
  heatmap_matrix[is.infinite(heatmap_matrix)] <- 0  # Handle -log10(0) which results in -Inf

  # Set the title of the heatmap
  if (is.null(title)) {
    title <- paste0("-log10(", value_type, ") by Cluster")
  }

  # Create the heatmap
  s_hmap <- heatmaply::heatmaply(
    heatmap_matrix,
    xlab = "Richset",
    ylab = "Cluster",
    main = title,
    colors = c("grey", viridis::viridis(256)),
    na.value = "grey",  # Set the color for NA values to grey
    margins = c(50, 50, 50, 50),
    show_dendrogram = c(TRUE, FALSE),
    plot_method = "plotly",
    colorbar_title = paste0("-log10(", value_type, ")")
  )

  return(s_hmap)
}


# Example usage
# s <- all_clusters_hmap(full_clusterdf10, "Padj")
# s

