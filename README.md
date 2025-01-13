## RichCluster
RichCluster is a fast C++ agglomerative hierarchical clustering algorithm packaged into easily callable R functions, designed to help cluster biological 'terms' based on how similar of genes are expressed in their activation. 

Terms are clustered together based on how many genes are shared between them. We support two different types of similarity scores:
- Kappa score
- Jaccard index

As well as different linkage criteria for iteratively merging clusters together.
- Multiple linkage (from DAVID method)
- Ward linkage

## Installation
The package is currently under review in submission to CRAN, but for now users can install the package and try out the clustering and visualization features by installing through GitHub.

```shell
install.packages("devtools")
require(devtools)
install_github("hyuncat/RichCluster")
```

## Usage
The basic flow is as follows. A demo can be followed along through the attached vignettes.

### Geneset formatting
Load in your desired genesets as R dataframes. `merge_genesets()` accepts a named list of genesets, the names you input here will be displayed as labels in the visualizations. You can edit these later.

The script will try to combine the genesets so that shared 'Term' rows will have their GeneID's merged. If you are receiving errors, please try renaming your columns of interest to 'Term' and 'GeneID' across all genesets of interest, and make sure the formatting of the columns is consistent (eg, same spellings) across all genesets.

Note, since the goal of the visualizations is to show the expression level of these clusters, we need to specify the format of the 'value' to display. Hence, `merge_genesets()` accepts a second argument 'value_columns', a character vector to specify the names of all the columns containing numeric values of interest to plot later. 

```r
rr1 <- read.csv(system.file("extdata", "7mo_DEG.csv", package="RichCluster"))

rr2 <- read.csv(system.file("extdata", "7mo_DMR.csv", package="RichCluster"))

richsets <- list(rr1, rr2)
richnames <- c('7mo_DEG', '7mo_DMR')
```

The output of `merge_genesets()` is a MergeObject, which is a named list with the corresponding data to easily plug and chug into the rest of the workflow. If you want to extract any part of it and tinker around, the relevant attributes are
- MergeObject$df - The big merged dataframe of all geneset data
- MergeObject$gs_names - Character vector of names with indices corresponding to the original index of the geneset in the argument.

## Clustering
To cluster the merged geneset, run `rich_cluster()` on the MergeObject. (Note this is only for convenience, to preserve the naming scheme of the genesets. If you would like to skip this you can alternatively call the overloaded method acting only on a dataframe object.)

Note the clustering script automatically looks for 'Term' and 'GeneID' columns to determine group candidates for clustering. If you are facing issues, please make sure your columns are named correctly or raise a GitHub issue.

A default clustering option can be run with Kappa similarity, but RichCluster can be customized with the following other options.

### Similarity metric
Type
- Kappa
- Jaccard index

Cutoff
- Only consider clusters where its terms have above this cutoff value of similarity.

### Merge strategy
- DAVID multiple linkage membership
- Single
- Complete
- Average
- Ward

The output of the clustering is a ClusterResult which can be directly inputted into the visualizations or exported as a csv file with some additional options.

The name of each cluster is determined as the term in the cluster with the highest gene count.

## Visualizations
### Export as CSV
Users can simply export the final clustered data as a CSV file and save their results for later. Then users can load in the file directly and visualize the clusters again.

### Multi-geneset heatmap
RichCluster's advantage is in its ability for highly customizable heatmap generation and comparison of enrichment data across multiple genesets.

Value options
- Name of the column value to visualize (column must be numeric type, and present across all genesets in ClusterResult)
- Type of column aggregation (mean, median, min, max, mode)

Additional options:
- Plot title
- Edit geneset names

![all_clusters_hmap 1](https://github.com/user-attachments/assets/0a0bac15-f1a9-404f-985c-3aa154f36f52)


### Intra-cluster similarity

![cluster7_correlation](https://github.com/user-attachments/assets/18a1aeb2-b221-41f4-a79a-1dcc5f3b2522)

### Cluster Network

<img src="https://github.com/user-attachments/assets/22021356-d352-45ff-a627-7e0156413f86" width="300">


