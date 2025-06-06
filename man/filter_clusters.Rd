% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster.R
\name{filter_clusters}
\alias{filter_clusters}
\title{Filter Clusters by Number of Terms}
\usage{
filter_clusters(all_clusters, min_terms)
}
\arguments{
\item{min_terms}{An integer specifying the minimum number of terms required in a cluster.}

\item{full_clusters}{A dataframe containing the merged seeds with column named `ClusterIndices`.}
}
\value{
The filtered data frame with clusters filtered to include only those with at least `min_terms` terms.
}
\description{
Filters the full list of clusters by keeping only those with greater
than or equal to min_terms # of terms.
}
