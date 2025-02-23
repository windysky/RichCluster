
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

hmap <- test_workflow(cluster_result, 10)
hmap




# heatmap creation
library(dplyr)

# utils
na_to_zero <- function(x) {
  # returns the column vector but replaces all places where it's nan/na/inf
  # to be zero (for hmap)
  x[is.nan(x) | is.na(x) | is.infinite(x)] <- 0
  return(x)
}

# Create heatmap with all clusters
cluster_hmap <- function(cluster_result, clusters=NULL, value_type="Padj", aggr_type=mean){
  # the hmap processing flow
  # for full_hmap
  cluster_df <- cluster_result$cluster_df

  # hmap_matrix
  hmap_matrix <- cluster_df %>%
    group_by(Cluster) %>%
    summarise(across(starts_with("Padj_"), function(x) mean(x, na.rm=TRUE))) %>%
    mutate(across(where(is.numeric), na_to_zero)) %>%
    select(-Cluster) %>% # remove (don't wanna plot her)
    mutate(across(where(is.numeric), function(x) -log10(x))) %>%
    mutate(across(where(is.numeric), function(x) ifelse(is.infinite(x), 0, x))) %>%
    as.matrix()

  # representative term making
  # use the Pvalue/Padj average column
  representative_terms <- cluster_df %>%
    group_by(Cluster) %>%
    filter(value_type==min(value_type, na.rm=TRUE)) %>%
    slice(1) %>%
    ungroup %>%
    pull(Term)
  rownames(hmap_matrix) <- representative_terms
  colnames(hmap_matrix) <- cluster_result$df_names

  # the hmap object
  hmap <- heatmaply::heatmaply(
    hmap_matrix,
    xlab = "Enrichment Result",
    ylab = "Cluster",
    main = "Title here",
    colors = c("grey", viridis::viridis(256)),
    na.value = "grey",
    margins = c(50, 50, 50, 50),
    show_dendrogram = c(TRUE, FALSE),
    plot_method = "plotly",
    colorbar_title = paste0("-log10(", value_type, ")")
  )
  return(hmap)
}

hmap <- cluster_hmap(cluster_result)

# clusters 4, 6, 8
clusters <- c("mating plug formation", "regulation of protein refolding", "regulation of plasma cell differentiation")
terms <- c("neuroblast proliferation", "regulation of tissue remodeling", "protein secretion")

# Heatmap displaying all terms in the specified clusters
# optionally accepts explicit list of terms to visualize if specified
# and we display the union of the two clusters/terms vectors in final result
#
# --- Representative Term Handling ---
# Clusters can be specified by cluster #
# or by the name of any term in the cluster
term_hmap <- function(cluster_result, clusters, terms, value_type, aggr_type) {

  cluster_df <- cluster_result$cluster_df

  # --- representative term search ---
  # convert character cluster terms -> cluster #s if not numeric
  if (is.null(clusters)) {
    clusters <- unique(cluster_df$Cluster)
  } else if (is.numeric(clusters)) {
    cluster_terms <- cluster_df %>%
      filter(Cluster %in% clusters) %>%
      pull(Cluster, Term, starts_with("Padj_"))
  } else if (is.character(clusters)) {
    # get the cluster numbers
    clusters <- cluster_df %>%
      filter(Term %in% clusters) %>%
      pull(Cluster) %>%
      unique()

    # use cluster numbers to get all terms in specified clusters
    cluster_terms <- cluster_df %>%
      filter(Cluster %in% clusters) %>%
      select(Cluster, Term, starts_with("Padj_"))

  } else {
    stop("`clusters` must be either numeric (Cluster #s) or character (Term names).")
  }
  # search for specific terms if supplied
  if (is.null(terms)) {
    # skip the following checks
    specific_terms <- c()
  } else if (!is.character(terms)) {
    stop("`terms` must be character (Term names).")
  } else {
    specific_terms <- cluster_df %>%
      filter(Term %in% terms) %>%
      select(Cluster, Term, starts_with("Padj_"))
  }

  # get the UNION of all terms in specified clusters
  # and those specified by terms
  final_terms <- bind_rows(specific_terms, cluster_terms) %>%
    distinct()

  # create the hmap_matrix
  # performing final value updates
  hmap_matrix <- final_terms %>%
    group_by(Cluster) %>%
    mutate(across(where(is.numeric), na_to_zero)) %>%
    # select(-Cluster) %>% # remove (don't wanna plot her)
    mutate(across(where(is.numeric), function(x) -log10(x))) %>%
    mutate(across(where(is.numeric), function(x) ifelse(is.infinite(x), 0, x)))

  # keep these vars for labeling
  cluster_annots <- hmap_matrix$Cluster
  row_names <- hmap_matrix$Term

  hmap_matrix <- hmap_matrix %>%
    ungroup() %>%
    select(starts_with("Padj_")) %>%
    as.matrix()
  rownames(hmap_matrix) <- row_names
  colnames(hmap_matrix) <- cluster_result$df_names

  h <- heatmaply(
    hmap_matrix,
    xlab = "Enrichment Result",
    ylab = "Term",
    main = "-log10(Padj) Values",
    colors = viridis::viridis(256),
    row_side_colors = cluster_annots,
    row_text_angle = 0,
    margins = c(60, 120, 40, 10),
    plot_method = "plotly",
    colorbar_title = "-log10(Padj)",
    cluster_rows=FALSE, cluster_cols=FALSE,
    Rowv=FALSE,
    Colv=FALSE,
  )
  return(h)
}



s <- cluster_df %>%
  s(starts_with("Padj"))
