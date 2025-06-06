% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hmaps.R
\name{cluster_hmap}
\alias{cluster_hmap}
\title{Create a Heatmap of Clustered Enrichment Results}
\usage{
cluster_hmap(
  cluster_result,
  clusters = NULL,
  value_type = "Padj",
  aggr_type = mean
)
}
\arguments{
\item{cluster_result}{A list containing a data frame (`cluster_df`) with clustering results.
The data frame must contain at least the columns `Cluster`, `Term`, and `value_type_*` values.}

\item{clusters}{Optional. A numeric or character vector specifying the clusters to include.
If NULL (default), all clusters are included.}

\item{value_type}{A character string specifying the column name prefix for values to display in hmap cells.
Defaults to `"Padj"`.}

\item{aggr_type}{A function used to aggregate values across clusters (e.g., `mean` or `median`).
Defaults to `mean`.}
}
\value{
An interactive heatmap object (`plotly`), displaying the -log10(Padj) values
  across clusters, with representative terms as row labels.
}
\description{
Generates an interactive heatmap from the given clustering results,
visualizing -log10(Padj) values for each cluster. The function aggregates
values per cluster and assigns representative terms as row names.
}
\details{
The function processes the given cluster data frame (`cluster_df`),
aggregating the `value_type_*` values per cluster using the specified `aggr_type` function.
The -log10 transformation is applied, and infinite values are replaced with 0.

Representative terms are selected by choosing the term with the lowest
`value_type` in each cluster.

The final heatmap is generated using `heatmaply::heatmaply()`, with
an interactive `plotly` visualization.
}
