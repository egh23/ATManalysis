## senmayo enrichment score
library(Seurat)
library(GSEABase)
library(escape)
library(AUCell)
library(tidyverse)
library(gridExtra)
library(patchwork)
library(readxl)
library(dittoSeq)
library(RColorBrewer)

data <- readRDS("data_S_elsie_all.rds")
data_sub <- data[,!(data$seurat_clusters %in% c(9, 12, 14))]
data_sub$seurat_clusters <- factor(data_sub$seurat_clusters)
data$sex <- data_sub$sex <- sub(".*_", "", data_sub$orig.ident)
data_male <- data_sub[, data_sub$sex=="Males"]
data_female <- data_sub[, data_sub$sex=="Females"]

gene_senmayo  <- read.csv("./SenMayo_mouse.csv",header=TRUE)
gene_senmayo <- unique(gene_senmayo$`Gene.murine.`)

gs1  <- GeneSet(gene_senmayo)
geneSets <- list(geneSet1=gs1)
geneSets$geneSet1@setName <- "SenMayo"
geneSets <- GeneSetCollection(geneSets)

ES_female <- enrichIt(obj = data_female, method = "UCell",gene.sets = geneSets, cores = 2)
data_female <- AddMetaData(data_female, ES_female)

ES_female <- enrichIt(obj = data_female, method = "UCell",gene.sets = geneSets, cores = 2)
data_female <- AddMetaData(data_female, ES_female)

## Pathway Enrichment
library(Seurat)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(DOSE)
library(pathview)
library(clusterProfiler)
library(AnnotationHub)
library(ensembldb)
library(tidyverse)
library(biomaRt)
library(msigdbr)
library(qvalue)
library(readxl)
library(gridExtra)

data <- readRDS("data_S_elsie_all.rds")
data_sub <- data[,!(data$seurat_clusters %in% c(9, 12, 14))]
data_sub$seurat_clusters <- factor(data_sub$seurat_clusters)

# helper function 
FindDEGs <- function(object, group.by = "all", compare.by = "orig.ident", assay = "RNA", features = NULL, min.cell = 25, 
                     min.pct = 0.01, logfc.threshold = log(2), pseudocount.use = 0.01, only.pos = T, ...){
  #object <- sub_data_S
  #compare.by  <- "combined_da_tmp"
  DefaultAssay(object) <- assay
  DEGs_list <- list()
  Idents(object) <- compare.by
  #print(Idents(object))
  if (group.by == "all"){
    object$tmp_all_for_DEG <- "all"
    group.by <- "tmp_all_for_DEG"
  }
  object_all <- object
  for (group_i in unique(na.omit(object_all[[group.by]][,1]))){
    #print(group_i)
    object <- object_all[, which(object_all[[group.by]] == group_i)]
    print(paste("Processing", group_i, "..."))
    if(min(table(object[[compare.by]])) < min.cell | length(table(object[[compare.by]])) == 1){
      print(table(object[[compare.by]]))
      print(paste("Cell number for some condition is below", min.cell))
      next
    }
    cluster_markers <- FindAllMarkers(
      object, assay = assay, features = features, only.pos = only.pos, min.pct = min.pct, logfc.threshold = logfc.threshold, 
      pseudocount.use = pseudocount.use, ...
    )
    cluster_markers$pct.diff <- cluster_markers$pct.1 - cluster_markers$pct.2
    #cluster_markers <- cbind(cluster_markers, AverageExpression(object, features = rownames(cluster_markers)), assays = assay)
    
    #for (compare_i in unique(object[[compare.by]][,1])){
    #sub_dat <- slot(object[, which(object[[compare.by]] == compare_i)]@assays[[assay]], name = "data")[rownames(cluster_markers),]
    #cluster_markers[, paste0("avg_", compare_i)] <- apply(sub_dat, 1, mean)
    #}
    DEGs_list[[group_i]] <- cluster_markers
  }
  DEGs_list
}

de_female_old_vs_young_per_cluster <- FindDEGs(data_female, group.by = "seurat_clusters", compare.by = "condition", features = rownames(data_female), logfc.threshold = 0, return.thresh = 1 )
de_male_old_vs_young_per_cluster <- FindDEGs(data_male, group.by = "seurat_clusters", compare.by = "condition", features = rownames(data_male), logfc.threshold = 0, return.thresh = 1 )

data_female$rest_vs_5 <- "5"
data_female$rest_vs_5[data_female$seurat_clusters != 5] <- "rest"
Idents(data_female) <- "rest_vs_5"
de_5_vs_rest <- FindMarkers(data_female, ident.1 = "5",ident.2 = "rest")
de_5_vs_rest$gene <- rownames(de_5_vs_rest)

ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
gene_ids <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                  filters = "external_gene_name",
                  values = de_5_vs_rest$gene,
                  mart = ensembl)
gene_ids <- distinct(gene_ids)
de_final <- merge(de_5_vs_rest, gene_ids, by.x = "gene", by.y = "external_gene_name")
de <- de_final[de_final$p_val_adj < 0.05,]

de <- de[de$avg_log2FC > 0,]

m_C5_BP <- msigdbr(species = "Mus musculus", category = "C5",subcategory = "GO:BP") %>% 
  dplyr::select(gs_name, ensembl_gene)
m_C5_CC <- msigdbr(species = "Mus musculus", category = "C5",subcategory = "GO:CC") %>% 
  dplyr::select(gs_name, ensembl_gene)
m_C5_MF <- msigdbr(species = "Mus musculus", category = "C5",subcategory = "GO:MF") %>% 
  dplyr::select(gs_name, ensembl_gene)
m_C5 <- bind_rows(m_C5_BP,m_C5_CC,m_C5_MF)

em_G <- enricher(de$ensembl_gene_id, TERM2GENE=m_C5, pAdjustMethod = "BH")