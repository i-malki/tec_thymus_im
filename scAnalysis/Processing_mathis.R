
library(tidyverse)
library(Seurat)
library(patchwork)
library(RColorBrewer)
library(pheatmap)
library(cowplot)
library(ape)
library(scales)

setwd("/shared/projects/tec_thymus/tec_thymus/GSE194252_diane_mathis_scRNAseq/input/matrix/")

date=Sys.Date()
# For only one data
#####
#read data
gene_table <- Read10X(data.dir = "GSM5831744_adult_perinate_gex")
hash_table <- Read10X(data.dir = "GSM5831745_adult_perinate_hash" , gene.column = 1)

seurat_table <- CreateSeuratObject( counts = gene_table, project = "Diane_Mathis", min.cells=3, min.features = 0 )

#####
#integrate hashes
hashtag_table <- t(read.table(file = "metadata/GSM5831744_adult_perinate_metadata.txt", header = TRUE, row.names = 1))

length(intersect(colnames(seurat_table), colnames(hashtag_table))) / length(union(colnames(seurat_table), colnames(hashtag_table)))

keep <- colnames(seurat_table) %in% colnames(hashtag_table)
seurat_table <- seurat_table[,keep]

#####filter data

grep("^mt", rownames(seurat_table), value = TRUE)

seurat_table[["percent.mt"]] <- PercentageFeatureSet(seurat_table, pattern = "^mt-")


x <- ggplot( seurat_table@meta.data, aes( x=orig.ident, y=nFeature_RNA ) ) +
  geom_jitter( color="gray ", size=0.5 ) +
  geom_violin( fill="steelblue" ) +
  theme_bw() +
  xlab("")

y <- ggplot( seurat_table@meta.data, aes( x=orig.ident, y=nCount_RNA ) ) +
  geom_jitter( color="gray ", size=0.5 ) +
  geom_violin( fill="pink" ) +
  theme_bw() +
  xlab("")

z <- ggplot( seurat_table@meta.data, aes( x=orig.ident, y=percent.mt ) ) +
  geom_hline( yintercept=25, color="red" ) +
  geom_jitter( color="gray ", size=0.5 ) +
  geom_violin( fill="orchid4" ) +
  theme_bw() +
  xlab("")

pdf( paste0("singlecell/figures/mathis_processed/qc_violin_mathis_processed-",date,".pdf" ), height=3, width=5 )
plot_grid( x,y,z, ncol=3 )
dev.off()

VlnPlot(seurat_table, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

seurat_table <- subset(seurat_table, subset = percent.mt < 25)



#####normalize data and find variable features

seurat_table <- NormalizeData(seurat_table, normalization.method = "LogNormalize", scale.factor = 10000)

nfeatures = 2000
seurat_table <- FindVariableFeatures(seurat_table, selection.method = "vst", nfeatures = nfeatures)

top10 <- head(VariableFeatures(seurat_table), 10)

plot1 <- VariableFeaturePlot(seurat_table)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

pdf( paste0("singlecell/figures/mathis_processed/variable_feature-",date,".pdf" ), height=4, width=6 )
plot1
plot2
dev.off()

#####dimensionality reduction

all.genes <- rownames(seurat_table)
seurat_table <- ScaleData(seurat_table, features = all.genes)

seurat_table <- RunPCA(seurat_table, features = VariableFeatures(object = seurat_table))

print(seurat_table[["pca"]], dims = 1:30, nfeatures = 5)


VizDimLoadings(seurat_table, dims = 1:30, reduction = "pca")
dev.off()

DimPlot(seurat_table, reduction = "pca")

DimHeatmap(seurat_table, dims = 1:30, cells = 500, balanced = TRUE)

ElbowPlot(seurat_table, ndims=50)
dev.off()

#####cluster and umap

seurat_table <- FindNeighbors(seurat_table, dims = 1:40)
seurat_table <- FindClusters(seurat_table, resolution = 1.8)
head(Idents(seurat_table), 100)
seurat_table <- RunUMAP(seurat_table, dims = 1:40, seed.use = 12345)
DimPlot(seurat_table, reduction="umap", label=T, repel=T )
dev.off()

FeaturePlot(seurat_table , features = "Gabrp" , label = T)

### Find markers to define clusters 

cluster.markers <- FindAllMarkers(seurat_table , min.pct = 0.1 , logfc.threshold = 0.5 , only.pos = T )
cluster.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top10.markers


########### Pour certain cluster spécifiques type ptf1a ou ionocytes qui sont faiblement representé ############

# A faire pour chaque marqueur spécifique de cellule
expression_matrix <- seurat_table@assays$RNA$scale.data

# Définir un seuil d'expression (par exemple, un seuil de 1)
threshold <- 1

# Identifier les cellules où le gène est exprimé au-dessus du seuil
expressing_cells <- colnames(expression_matrix)[expression_matrix["Foxi1", ] > 1]

# Calculer le pourcentage de cellules exprimant le gène
percentage_expressing <- length(expressing_cells) / ncol(seurat_table) * 100

# Afficher le pourcentage
print(paste("Pourcentage de cellules exprimant le gène :", round(percentage_expressing, 2), "%"))


cell_clust <- as.character(seurat_table$seurat_clusters)
names(cell_clust) <- names(seurat_table$seurat_clusters)
cell_clust[expressing_cells3] <- "Goblets"
seurat_table$seurat_clusters <- cell_clust
seurat_table@active.ident <- as.factor(seurat_table$seurat_clusters)
DimPlot(seurat_table , reduction = "umap" , label = TRUE , repel = TRUE)

############################################################################################################
############ Pour les cellules s'exprimant dans plusieurs clusters d emaniere aspécifique ###########
# Créer une copie de seurat_clusters pour éviter la modification directe
cell_clust <- as.character(seurat_table$seurat_clusters)
names(cell_clust) <- names(seurat_table$seurat_clusters)

# Filtrer les cellules 
cells_in_cluster_ <- names(cell_clust)[cell_clust == "12"]

# Filtrer les cellules qui sont dans `expressing_cells` et dans le cluster 14
cells_to_modify <- intersect(cells_in_cluster_, expressing_cells2)

# Modifier `cell_clust` uniquement pour les cellules correspondantes
cell_clust[cells_to_modify] <- "Keratinocytes"

# Réassigner les nouveaux clusters à `seurat_clusters`
seurat_table$seurat_clusters <- cell_clust

# Mettre à jour l'identifiant actif
seurat_table@active.ident <- as.factor(seurat_table$seurat_clusters)


# Rename Idents avec les markers trouvés avant 

seurat_table <- RenameIdents(seurat_table, "0"="Tuft1")
seurat_table <- RenameIdents(seurat_table, "1"="Tuft1")
seurat_table <- RenameIdents(seurat_table, "9"="Tuft2")
seurat_table <- RenameIdents(seurat_table, "2"="Tuft2")
seurat_table <- RenameIdents(seurat_table, "15"="Immature mTEC")
seurat_table <- RenameIdents(seurat_table, "6"="Immature mTEC")
seurat_table <- RenameIdents(seurat_table, "5"="Immature mTEC")
seurat_table <- RenameIdents(seurat_table, "7"="Immature mTEC")
seurat_table <- RenameIdents(seurat_table, "3"="Immature mTEC")
seurat_table <- RenameIdents(seurat_table, "4"="Immature mTEC")
seurat_table <- RenameIdents(seurat_table, "17"="adult cTEC")
seurat_table <- RenameIdents(seurat_table, "18"="perinatal cTEC")
seurat_table <- RenameIdents(seurat_table, "16"="Transit Amplifying")
seurat_table <- RenameIdents(seurat_table, "11"="Microfold")
seurat_table <- RenameIdents(seurat_table, "20"="Ciliated")
seurat_table <- RenameIdents(seurat_table, "14"="Neuroendocrine")
seurat_table <- RenameIdents(seurat_table, "10"="Enterocytes/hepatocytes")
seurat_table <- RenameIdents(seurat_table, "13"="Basal(skin)")
seurat_table <- RenameIdents(seurat_table, "12"="Basal(lung)")
seurat_table <- RenameIdents(seurat_table, "19"="Muscle")
seurat_table <- RenameIdents(seurat_table, "8"="Aire expressing")





DimPlot(seurat_table, reduction="umap", label=T, repel=T   )


table(Idents(seurat_table))


### Save file ###
saveRDS(seurat_table , "../../../singlecelldata/processed_mathis_single_cell.rds")

#####
#build cluster tree
seurat_table <- BuildClusterTree(seurat_table, dims=T)
PlotClusterTree(seurat_table)

#