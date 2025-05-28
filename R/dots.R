
#' Create a Dot Plot of Cluster Enrichment
#'
#' Generates a dot plot visualizing enrichment scores (e.g., -log10(Padj)) for
#' different clusters. The size of the dots can represent the number of terms in each cluster.
#' Representative terms are used as labels for clusters.
#'
#' @param cluster_result The result object from `RichCluster::cluster()`.
#' @param clusters Optional. A character or numeric vector specifying which clusters to include.
#'                 If NULL, all clusters are included.
#' @param value_type A string indicating the type of value to plot (e.g., "Padj", "Pvalue").
#'                   Defaults to "Padj".
#' @param title Optional. A character string for the plot title. If NULL, a default
#'              title is generated.
#' @return A plotly object representing the dot plot.
#' @export
cluster_dot <- function(cluster_result, clusters=NULL, value_type="Padj", title=NULL) {
  cluster_df <- cluster_result$cluster_df
  df_names <- cluster_result$df_names

  if (is.null(clusters)) {
    clusters <- unique(cluster_df$Cluster)
  }

  # dot data
  dot_data <- cluster_df magrittr::%>%
    dplyr::group_by(Cluster) magrittr::%>%
    dplyr::mutate(n_terms = dplyr::n()) magrittr::%>% # get n_terms count
    dplyr::summarise(
      n_terms = dplyr::first(n_terms), # and assign each cluster with n_terms
      dplyr::across(dplyr::starts_with(paste0(value_type, "_")), function(x) mean(x, na.rm=TRUE))) magrittr::%>%
    dplyr::mutate(dplyr::across(dplyr::starts_with(paste0(value_type, "_")), function(x) -log10(x))) magrittr::%>%
    dplyr::mutate(dplyr::across(dplyr::starts_with(paste0(value_type, "_")), function(x) ifelse(is.infinite(x), 0, x)))

  # representative term making
  # use the Pvalue/Padj average column
  representative_terms <- cluster_df magrittr::%>%
    dplyr::group_by(Cluster) magrittr::%>%
    dplyr::filter(!!rlang::sym(value_type) == min(!!rlang::sym(value_type), na.rm=TRUE)) magrittr::%>%
    dplyr::slice(1) magrittr::%>%
    dplyr::ungroup() magrittr::%>%
    dplyr::pull(Term)

  dot_data$ClusterName <- representative_terms

  # generate default title if none supplied
  if (is.null(title)) {
    title <- paste0(value_type, ", all clusters")
  }

  # get all dot data across different dfs
  value_cols <- grep(paste0("^", value_type, "_"), names(dot_data), value = TRUE)

  # the plotly object
  dot <- plotly::plot_ly(
    dot_data,
    x = dot_data[[value_cols[1]]],
    y = ~ClusterName,
    type = 'scatter',
    mode = 'markers',
    marker = list(opacity = 0.5),
    size = ~n_terms,
    name = df_names[1],
    text = ~paste0(
      'Cluster ', Cluster, '<br>',
      ClusterName, '<br>',
      '-log10(',value_type,') = ', sprintf("%.2f", dot_data[[value_cols[1]]]), '<br>',
      '# terms = ', n_terms
    ),
    hovertemplate = "%{text}"
  )
  if (length(value_cols) > 1) {
    for (i in 2:length(value_cols)) {
      dot <- dot %>% plotly::add_trace(
        x = dot_data[[value_cols[i]]],
        type = 'scatter',
        marker = list(opacity = 0.5),
        size = ~n_terms,
        name = df_names[i],
        text = ~paste0(
          'Cluster ', Cluster, '<br>',
          ClusterName, '<br>',
          '-log10(',value_type,') = ', sprintf("%.2f", dot_data[[value_cols[i]]]), '<br>',
          '# terms = ', n_terms
        ),
        hovertemplate = "%{text}"
      )
    }
  }
  # add layout stuff
  dot <- dot%>% plotly::layout(
    title = title,
    margin = list(autoexpand=TRUE, t=50, b=-50),
    xaxis = list(title = "-log10(Padj)"),
    yaxis = list(title = "Cluster")
  )
  return(dot)
}
# cdot <- cluster_dot(cluster_result)
# cdot

#' Create a Dot Plot of Term Enrichment within a Cluster
#'
#' Generates a dot plot visualizing enrichment scores (e.g., -log10(Padj)) for
#' individual terms within a specified cluster. The size of the dots can represent
#' the number of genes associated with each term.
#'
#' @param cluster_result The result object from `RichCluster::cluster()`.
#' @param cluster The cluster to visualize. Can be a numeric cluster ID or a
#'                character string representing a term within the desired cluster.
#'                Defaults to cluster 1.
#' @param value_type A string indicating the type of value to plot (e.g., "Padj", "Pvalue").
#'                   Defaults to "Padj".
#' @param title Optional. A character string for the plot title. If NULL, a default
#'              title is generated based on the representative term of the cluster.
#' @return A plotly object representing the dot plot.
#' @export
term_dot <- function(cluster_result, cluster=1, value_type="Padj", title=NULL) {

  cluster_df <- cluster_result$cluster_df
  df_names <- cluster_result$df_names

  # get cluster # from the term
  if (is.character(cluster)) {
    cluster <- cluster_df[cluster_df$Term==cluster, ]$Cluster
  } else if (!is.numeric(cluster)) {
    stop("cluster must be numeric (a cluster number) or character (a term name)")
  }

  # dot data
  dot_data <- cluster_df magrittr::%>%
    dplyr::group_by(Cluster) magrittr::%>%
    dplyr::filter(Cluster==cluster) magrittr::%>%
    # dplyr::summarise(dplyr::across(dplyr::starts_with(paste0(value_type, "_")), function(x) mean(x, na.rm=TRUE))) magrittr::%>%
    dplyr::mutate(dplyr::across(dplyr::starts_with(paste0(value_type, "_")), function(x) -log10(x))) magrittr::%>%
    dplyr::mutate(dplyr::across(dplyr::starts_with(paste0(value_type, "_")), function(x) ifelse(is.infinite(x), 0, x))) magrittr::%>%
    dplyr::mutate(n_genes=sapply(strsplit(GeneID, ','), length))

  # representative term making
  # use the Pvalue/Padj average column
  representative_term <- dot_data magrittr::%>%
    dplyr::group_by(Cluster) magrittr::%>%
    dplyr::filter(!!rlang::sym(value_type) == min(!!rlang::sym(value_type), na.rm=TRUE)) magrittr::%>%
    dplyr::slice(1) magrittr::%>%
    dplyr::ungroup() magrittr::%>%
    dplyr::pull(Term)

  # use it in the default title (if none supplied)
  if (is.null(title)) {
    title <- paste0(value_type, ", ", representative_term, " (cluster ", cluster, ")")
  }

  # get all dot data across different dfs
  value_cols <- grep(paste0("^", value_type, "_"), names(dot_data), value = TRUE)

  # the plotly object
  dot <- plotly::plot_ly(
    dot_data,
    x = dot_data[[value_cols[1]]],
    y = ~Term,
    type = 'scatter',
    mode = 'markers',
    marker = list(opacity = 0.5),
    size = ~n_genes,
    name = df_names[1],
    text = ~paste0(
      Term, '<br>',
      '-log10(',value_type,') = ', sprintf("%.2f", dot_data[[value_cols[1]]]), '<br>',
      '# genes = ', n_genes
    ),
    hovertemplate = "%{text}"
  )
  if (length(value_cols) > 1) {
    for (i in 2:length(value_cols)) {
      dot <- dot %>% plotly::add_trace(
        x = dot_data[[value_cols[i]]],
        type = 'scatter',
        marker = list(opacity = 0.5),
        size = ~n_genes,
        name = df_names[i],
        text = ~paste0(
          Term, '<br>',
          '-log10(',value_type,') = ', sprintf("%.2f", dot_data[[value_cols[i]]]), '<br>',
          '# terms = ', n_genes
        ),
        hovertemplate = "%{text}"
      )
    }
  }
  # add layout stuff
  dot <- dot%>% plotly::layout(
    title = title,
    margin = list(autoexpand=TRUE, t=50, b=-50),
    xaxis = list(title = "-log10(Padj)"),
    yaxis = list(title = "Term")
  )
  return(dot)
}

# tdot <- term_dot(cluster_result, cluster=48)
# tdot
