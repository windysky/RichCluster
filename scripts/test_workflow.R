
load_cluster_data <- function(from_scratch=FALSE)
{
  if (from_scratch == FALSE)  {
    # load the cluster data from the files
    cluster_result <- readRDS('inst/extdata/cluster_result.rds')
  }
  else {
    # read and manually cluster
    rr1 <- read.delim(system.file("extdata", "go1.txt", package="RichCluster"))
    rr2 <- read.delim(system.file("extdata", "go2.txt", package="RichCluster"))

    enrichment_results <- list(rr1, rr2)
    rr_names <- c('7mo_DEG', '7mo_DMR')

    cluster_result <- RichCluster::cluster(
      enrichment_results, df_names=rr_names, min_terms=5,
      distance_metric="kappa", distance_cutoff=0.5,
      merge_strategy="DAVID", membership_cutoff=0.5
    )
    saveRDS(cluster_result, file = "inst/extdata/cluster_result.rds")
  }
  return(cluster_result)
}

# example of what workflow using the functions looks like
test_workflow <- function(cluster_result, min_terms=5) {

  # Testing: Choose one to comment out and test
  # result <- RichCluster::all_clusters_hmap(full_clusterdf, "Padj")
  # result <- RichCluster::cluster_correlation_hmap(final_clusters, distance_matrix, 3)
  result <- RichCluster::cluster_network(final_clusters, distance_matrix, 1)
  # result <- RichCluster::full_network(distance_matrix[1:30, 1:30])
  return(result)

}

cluster_result <- load_cluster_data(from_scratch=TRUE)
# cluster_result <- load_cluster_data(from_scratch=FALSE)

hmap <- test_workflow(cluster_result, 10)
hmap
