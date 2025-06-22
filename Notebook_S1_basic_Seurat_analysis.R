## Basic Analysis
library(Seurat)
library(dplyr)
library(Matrix)

dataFolder <- "/data/rihao/elsie_project_202307/data"

dataFolders <- list()
dataFolders[[1]] <- paste(dataFolder, "/Old_Females_MMT_cellranger/filtered_feature_bc_matrix", sep="")
dataFolders[[2]] <- paste(dataFolder, "/Old_Males_MMT_cellranger/filtered_feature_bc_matrix", sep="")
dataFolders[[3]] <- paste(dataFolder, "/Young_Females_MMT_cellranger/filtered_feature_bc_matrix", sep="")
dataFolders[[4]] <- paste(dataFolder, "/Young_Males_MMT_cellranger/filtered_feature_bc_matrix", sep="")

fdata <- list()
fdata[[1]] <- Read10X(dataFolders[[1]])
fdata[[2]] <- Read10X(dataFolders[[2]])
fdata[[3]] <- Read10X(dataFolders[[3]])
fdata[[4]] <- Read10X(dataFolders[[4]])


whole <- list()
whole[[1]] <- CreateSeuratObject(counts=fdata[[1]], project = "Old_Females")
whole[[2]] <- CreateSeuratObject(counts=fdata[[2]], project = "Old_Males")
whole[[3]] <- CreateSeuratObject(counts=fdata[[3]], project = "Young_Females")
whole[[4]] <- CreateSeuratObject(counts=fdata[[4]], project = "Young_Males")

#merge all samples
whole <- merge(x=whole[[1]], y=list(whole[[2]], whole[[3]], whole[[4]]), add.cell.ids = c("Old_Females", "Old_Males", "Young_Females", "Young_Males" ))
whole <- merge(x=whole[[2]], y=whole[[4]], add.cell.ids = c("Old_Males", "Young_Males" ))

## Number of cells before
cells.before <- length(colnames(x= whole))

## NORMALIZATION
mito.genes <- c(grep("^MT-", rownames(x = whole), value = T),
                grep("^mt-", rownames(x = whole), value = T))
percent.mito <- Matrix::colSums(x = GetAssayData(object = whole, slot = 'counts')[mito.genes, ]) / Matrix::colSums(x = GetAssayData(object = whole, slot = 'counts'))
whole[['percent.mito']] <- percent.mito
VlnPlot(whole, features = "percent.mito", pt.size = 0.1) + NoLegend()
whole <- subset(x = whole, subset = percent.mito <= 0.10 & nFeature_RNA >= 363)
table(whole$orig.ident)
whole <- NormalizeData(object = whole, normalization.method = "LogNormalize", scale.factor = 10000)
whole <- FindVariableFeatures(object = whole, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))
whole <- ScaleData(object = whole, features = VariableFeatures(object = whole), vars.to.regress = c("nCount_RNA", "percent.mito"))
gc()

## PCA and Dimensionality Reduction
whole <- RunPCA(object = whole,
                features =  VariableFeatures(object = whole),
                verbose=FALSE)

whole <- RunUMAP(object = whole, dims = 1:20)

## CLUSTERING
whole <- FindNeighbors(object = whole, dims = 1:20)
whole <- FindClusters(object = whole, resolution = 0.5)

## DEGs
whole.markers <- FindAllMarkers(object = whole,
                                only.pos = TRUE,
                                min.pct = 0.10,
                                logfc.threshold = 0.10)

## SAVING
save(whole, file = "whole_object.Robj")


# condition
data_S_elsie_all <- whole

DimPlot(data_S_elsie_all, group.by = "orig.ident")
table(data_S_elsie_all$orig.ident)

data_S_elsie_all$condition <- "Young"
data_S_elsie_all$condition[which(data_S_elsie_all$orig.ident %in% c("Old_Females","Old_Males"))] <- "Old"

DimPlot(data_S_elsie_all, group.by = "condition", shuffle = T)

data_S_elsie_all -> data_S
p <- ggplot(data_S@meta.data, aes(x=seurat_clusters, fill=condition)) + geom_bar() + #position = "fill"
  theme(text = element_text(size=15))  + theme(legend.title = element_blank()) 
p


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

data_S <- data_S_elsie_all
DEG_lists <- FindDEGs(data_S, group.by = "seurat_clusters", compare.by = "condition", features = rownames(data_S), logfc.threshold = 0, return.thresh = 1 )
names(DEG_lists)

DEG_lists_all_genes <- DEG_lists

setwd("/data/rihao/elsie_project_202307/DEGs_per_cluster_all_genes")
for (i in names(DEG_lists)){
  openxlsx::write.xlsx(DEG_lists[[i]], sprintf("cluster_%s_DEGs.xlsx",i))
}