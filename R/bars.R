
#' Create a Bar Plot of Cluster Enrichment
#'
#' Generates a bar plot visualizing enrichment scores (e.g., -log10(Padj)) for
#' different clusters. Representative terms are used as labels for clusters.
#'
#' @param cluster_result The result object from `RichCluster::cluster()`.
#' @param clusters Optional. A character or numeric vector specifying which clusters to include.
#'                 If NULL, all clusters are included.
#' @param value_type A string indicating the type of value to plot (e.g., "Padj", "Pvalue").
#'                   Defaults to "Padj".
#' @param title Optional. A character string for the plot title. If NULL, a default
#'              title is generated.
#' @return A plotly object representing the bar plot.
#' @export
cluster_bar <- function(cluster_result, clusters=NULL, value_type="Padj", title=NULL) {

  cluster_df <- cluster_result$cluster_df
  df_names <- cluster_result$df_names

  # bar data
  bar_data <- cluster_df magrittr::%>%
    dplyr::group_by(Cluster) magrittr::%>%
    dplyr::summarise(dplyr::across(dplyr::starts_with(paste0(value_type, "_")), function(x) mean(x, na.rm=TRUE))) magrittr::%>%
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

  bar_data$Cluster <- representative_terms

  # generate default title if none supplied
  if (is.null(title)) {
    title <- paste0(value_type, ", all clusters")
  }

  # get all bar data across different dfs
  value_cols <- grep(paste0("^", value_type, "_"), names(bar_data), value = TRUE)

  # the plotly object
  bar <- plotly::plot_ly(
    bar_data,
    x = bar_data[[value_cols[1]]],
    y = ~Cluster,
    type = 'bar',
    name = df_names[1]
  )
  if (length(value_cols) > 1) {
    for (i in 2:length(value_cols)) {
      bar <- bar %>% plotly::add_trace(
        x = bar_data[[value_cols[i]]],
        name = df_names[i]
      )
    }
  }
  # add layout stuff
  bar <- bar%>% plotly::layout(
    title = title,
    margin = list(autoexpand=TRUE, t=50, b=-50),
    xaxis = list(title = "-log10(Padj)"),
    yaxis = list(title = "Term")
  )
  return(bar)
}

# cbar <- cluster_bar(cluster_result)
# cbar

#' Create a Bar Plot of Term Enrichment within a Cluster
#'
#' Generates a bar plot visualizing enrichment scores (e.g., -log10(Padj)) for
#' individual terms within a specified cluster.
#'
#' @param cluster_result The result object from `RichCluster::cluster()`.
#' @param cluster The cluster to visualize. Can be a numeric cluster ID or a
#'                character string representing a term within the desired cluster.
#'                Defaults to cluster 1.
#' @param value_type A string indicating the type of value to plot (e.g., "Padj", "Pvalue").
#'                   Defaults to "Padj".
#' @param title Optional. A character string for the plot title. If NULL, a default
#'              title is generated based on the representative term of the cluster.
#' @return A plotly object representing the bar plot.
#' @export
term_bar <- function(cluster_result, cluster=1, value_type="Padj", title=NULL) {

  cluster_df <- cluster_result$cluster_df
  df_names <- cluster_result$df_names

  # get cluster # from the term
  if (is.character(cluster)) {
    cluster <- cluster_df[cluster_df$Term==cluster, ]$Cluster
  } else if (!is.numeric(cluster)) {
    stop("cluster must be numeric (a cluster number) or character (a term name)")
  }

  # bar data
  bar_data <- cluster_df magrittr::%>%
    dplyr::group_by(Cluster) magrittr::%>%
    dplyr::filter(Cluster==cluster) magrittr::%>%
    # dplyr::summarise(dplyr::across(dplyr::starts_with(paste0(value_type, "_")), function(x) mean(x, na.rm=TRUE))) magrittr::%>%
    dplyr::mutate(dplyr::across(dplyr::starts_with(paste0(value_type, "_")), function(x) -log10(x))) magrittr::%>%
    dplyr::mutate(dplyr::across(dplyr::starts_with(paste0(value_type, "_")), function(x) ifelse(is.infinite(x), 0, x)))

  # representative term making
  # use the Pvalue/Padj average column
  representative_term <- bar_data magrittr::%>%
    dplyr::group_by(Cluster) magrittr::%>%
    dplyr::filter(!!rlang::sym(value_type) == min(!!rlang::sym(value_type), na.rm=TRUE)) magrittr::%>%
    dplyr::slice(1) magrittr::%>%
    dplyr::ungroup() magrittr::%>%
    dplyr::pull(Term)

  # use it in the default title (if none supplied)
  if (is.null(title)) {
    title <- paste0(value_type, ", ", representative_term, " (cluster ", cluster, ")")
  }

  # get all bar data across different dfs
  value_cols <- grep(paste0("^", value_type, "_"), names(bar_data), value = TRUE)

  # the plotly object
  bar <- plotly::plot_ly(
    bar_data,
    x = bar_data[[value_cols[1]]],
    y = ~Term,
    type = 'bar',
    name = df_names[1]
  )
  if (length(value_cols) > 1) {
    for (i in 2:length(value_cols)) {
      bar <- bar %>% plotly::add_trace(
        x = bar_data[[value_cols[i]]],
        name = df_names[i]
      )
    }
  }
  # add layout stuff
  bar <- bar%>% plotly::layout(
    title = title,
    margin = list(autoexpand=TRUE, t=50, b=-50),
    xaxis = list(title = "-log10(Padj)"),
    yaxis = list(title = "Term")
  )
  return(bar)
}

# tbar <- term_bar(cluster_result, cluster=48)
# tbar
