
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

cluster_result <- load_cluster_data(from_scratch=FALSE)

# ALL VISUALIZATION TESTS
# ---
c_hmap <- RichCluster::cluster_hmap(cluster_result)
c_hmap

# clusters 4, 6, 8 + terms from clusters x,y,z
clusters <- c("mating plug formation", "regulation of protein refolding", "regulation of plasma cell differentiation")
terms <- c("neuroblast proliferation", "regulation of tissue remodeling", "protein secretion")
t_hmap <- RichCluster::term_hmap(cluster_result, clusters, terms, value_type="Padj")
t_hmap

c_bar <- RichCluster::cluster_bar(cluster_result)
c_bar

t_bar <- RichCluster::term_bar(cluster_result, 48)
t_bar

c_dot <- RichCluster::cluster_dot(cluster_result)
c_dot

t_dot <- RichCluster::term_dot(cluster_result, 48)
t_dot

cluster_df <- RichCluster::export_df(cluster_result)
# write.csv(cluster_df, "~/Downloads/cluster_df.csv", row.names=FALSE)
