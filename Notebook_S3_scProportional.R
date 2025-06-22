library("scProportionTest")
data <- readRDS("data_S_elsie_all.rds")
data_sub <- data[,!(data$seurat_clusters %in% c(9, 12, 14))]
data_sub$seurat_clusters <- factor(data_sub$seurat_clusters)
data_sub@meta.data$sex <- gsub(".*_", "", data_sub$orig.ident)
prop_test <- sc_utils(data_sub)
prop_test <- permutation_test(
  prop_test, cluster_identity = "seurat_clusters",
  sample_1 = "Females", sample_2 = "Males",
  sample_identity = "sex"
)
permutation_plot(prop_test)