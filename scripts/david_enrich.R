library(dplyr)
library(richR)

deg <- read.delim(system.file("extdata", "deg_mouse1.txt", package="RichCluster"))
filtered_deg <- deg %>% filter(padj < 0.05)
