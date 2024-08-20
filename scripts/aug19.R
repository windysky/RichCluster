

m7_deg <- read.csv(system.file("extdata", "7mo_DEG.csv", package="RichCluster"))
m7_dmr <- read.csv(system.file("extdata", "7mo_DMR.csv", package="RichCluster"))
m10_deg <- read.csv(system.file("extdata", "10mo_DEG.csv", package="RichCluster"))
m10_dmr <- read.csv(system.file("extdata", "10mo_DMR.csv", package="RichCluster"))

m7_richsets <- list(m7_deg, m7_dmr)
m10_richsets <- list(m10_deg, m10_dmr)

all_richsets <- list(m7_deg, m7_dmr, m10_deg, m10_dmr)

m7_names <- c('7mo_DEG', '7mo_DMR')
m10_names <- c('10mo_DEG', '10mo_DMR')
all_names <- c(m7_names, m10_names)

cluster_workflow <- function(merged_richsets) {
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
  return(cluster_result)
}

# For month 7
merged_m7_richsets <- RichCluster::merge_richsets(m7_richsets)
m7_result <- cluster_workflow(merged_m7_richsets)

m7_dist_matrix <- m7_result$DistanceMatrix
m7_merged_seeds <- m7_result$MergedSeeds

m7_final_clusters <- RichCluster::filter_merged_seeds(m7_merged_seeds, min_terms=3)
m7_full_clusterdf <- RichCluster::make_full_clusterdf(m7_final_clusters, merged_m7_richsets)

m7_hmap <- RichCluster::all_clusters_hmap(
  m7_full_clusterdf,
  value_type="Padj",
  x_labels=m7_names,
  title="Month 7 DEG vs. DMR"
)
m7_hmap # view the heatmap (best viewed in browser)

# For month 10
merged_m10_richsets <- RichCluster::merge_richsets(m10_richsets)
m10_result <- cluster_workflow(merged_m10_richsets)

m10_dist_matrix <- m10_result$DistanceMatrix
m10_merged_seeds <- m10_result$MergedSeeds

m10_final_clusters <- RichCluster::filter_merged_seeds(m10_merged_seeds, min_terms=3)
m10_full_clusterdf <- RichCluster::make_full_clusterdf(m10_final_clusters, merged_m10_richsets)

m10_hmap <- RichCluster::all_clusters_hmap(
  m10_full_clusterdf,
  value_type="Padj",
  x_labels=m10_names,
  title="Month 10 DEG vs. DMR"
)
m10_hmap # view the heatmap

# For all months
merged_all_richsets <- RichCluster::merge_richsets(all_richsets)
all_result <- cluster_workflow(merged_all_richsets)

all_dist_matrix <- all_result$DistanceMatrix
all_merged_seeds <- all_result$MergedSeeds

all_final_clusters <- RichCluster::filter_merged_seeds(all_merged_seeds, min_terms=3)
all_full_clusterdf <- RichCluster::make_full_clusterdf(all_final_clusters, merged_all_richsets)

all_hmap <- RichCluster::all_clusters_hmap(
  all_full_clusterdf,
  value_type="Padj",
  x_labels=all_names,
  title="All Months DEG vs. DMR"
)
all_hmap # view the heatmap
