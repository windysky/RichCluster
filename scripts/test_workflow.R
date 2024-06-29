
load_cluster_data <- function(from_scratch=FALSE)
{
  if (from_scratch == FALSE)  {
    # Load the cluster data from CSV files
    distance_matrix_df <- read.csv("inst/extdata/distance_matrix.csv", row.names=1, check.names=FALSE)
    distance_matrix <- as.matrix(distance_matrix_df)
    merged_seeds <- read.csv("inst/extdata/merged_seeds.csv")
    # Make the fake cluster_result with the loaded data :^)
    cluster_result <- list(DistanceMatrix=distance_matrix, MergedSeeds=merged_seeds)
  }
  else {
    # Read the data files
    rr1 <- read.delim(system.file("extdata", "go1.txt", package="RichCluster"))
    rr2 <- read.delim(system.file("extdata", "go2.txt", package="RichCluster"))

    richsets <- list(rr1, rr2)
    richnames <- c('go1', 'go2')

    merged_richsets <- RichCluster::merge_richsets(richsets)

    term_vec <- merged_richsets$Term
    geneID_vec <- merged_richsets$GeneID
    padj_vec <- merged_richsets$Padj

    cluster_result <- RichCluster::RichCluster(
      "kappa", 0.5,
      "DAVID", 0.5,
      term_vec,
      geneID_vec,
      padj_vec
    )

    # Save cluster results to CSV for faster future testing
    dist_matrix_df <- as.data.frame(cluster_result$DistanceMatrix)
    write.csv(dist_matrix_df, "inst/extdata/distance_matrix.csv", row.names=TRUE)
    write.csv(cluster_result$MergedSeeds, "inst/extdata/merged_seeds.csv", row.names=FALSE)

  }
  return(cluster_result)
}

# Example of what workflow using the functions looks like
test_workflow <- function(cluster_result, min_terms=5) {
  # Unpack distance_matrix and merged_seeds from list
  distance_matrix <- cluster_result$DistanceMatrix
  merged_seeds <- cluster_result$MergedSeeds

  final_clusters <- RichCluster::filter_merged_seeds(merged_seeds, min_terms)
  full_clusterdf <- RichCluster::make_full_clusterdf(final_clusters, merged_richsets)

  if (FALSE) {
    # Sanity check: Does summary_data look okay
    summary_data <- RichCluster::make_summary_data(full_clusterdf, "Padj")
    print(nrow(summary_data))
  }

  # Testing: Choose one to comment out and test
  # all_hmap <- RichCluster::all_clusters_hmap(full_clusterdf, "Padj")
  # cluster_correlation_hmap <- RichCluster::cluster_correlation_hmap(final_clusters, distance_matrix, 1)
  # cluster_network <- RichCluster::cluster_network(final_clusters, distance_matrix, 1)
  full_network <- RichCluster::full_network(distance_matrix[1:30, 1:30])
  return(full_network)
}

# cluster_result <- load_cluster_data(from_scratch=TRUE)
cluster_result <- load_cluster_data(from_scratch=FALSE)

all_hmap_5 <- test_workflow(cluster_result, 5)
all_hmap_5

all_hmap_10 <- test_workflow(cluster_result, 10)
all_hmap_10




