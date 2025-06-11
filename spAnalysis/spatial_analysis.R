library(Seurat)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(patchwork)
library(scigenex)

## Load basic 
setwd("spatial/Young_6wk_2021/")
st_data <- Load10X_Spatial("Visium_files")

st_data <- NormalizeData(st_data)
st_data <- ScaleData(st_data)
st_data <- FindVariableFeatures(st_data)
st_data <- RunPCA(st_data, verbose = FALSE)
st_data <- FindNeighbors(st_data, reduction = "pca", dims = 1:20)
st_data <- FindClusters(st_data, resolution = 0.9)
st_data <- RunUMAP(st_data, reduction = "pca", dims = 1:20)
SpatialDimPlot(st_data, pt.size.factor = 4)
SpatialDimPlot(st_data, pt.size.factor = 0)
SpatialFeaturePlot(st_data, features = "nCount_Spatial")

## Select spot to be kept
## based on HE (the tissue seems 
## to be folded in some regions)

coord <- scigenex::getFlippedTissueCoordinates(st_data, as_data_frame = TRUE)
plot(coord)
p1 <- list(x=950, y=1050)
p2 <- list(x=1030, y=400)
slop <- (p2$y - p1$y) / (p2$x - p1$x) 
intersect <- p2$y - slop * p2$x
abline(b=slop, a=intersect, col="green")
points(coord[coord$y > coord$x * slop + intersect,], pch=16, col="red")
coord <- coord[coord$y > coord$x * slop + intersect,]
abline(h=500, col="green")
plot(coord)
points(coord[coord$y > 500,], pch=16, col="red")
coord <- coord[coord$y > 500,]

spot_to_keep <- rownames(coord)

## Reload the dataset
## Elimate spot to be deleted
## Redo analysis

st_data <- Load10X_Spatial("Visium_files")
st_data <- subset(st_data, cells= spot_to_keep)
st_data <- NormalizeData(st_data)
st_data <- ScaleData(st_data)
st_data <- FindVariableFeatures(st_data)
st_data <- RunPCA(st_data, verbose = FALSE)
st_data <- FindNeighbors(st_data, reduction = "pca", dims = 1:20)
st_data <- FindClusters(st_data, resolution = 0.9)
st_data <- RunUMAP(st_data, reduction = "pca", dims = 1:20)
DimPlot(st_data)
SpatialDimPlot(st_data, pt.size.factor = 4)
SpatialDimPlot(st_data, pt.size.factor = 0)
SpatialFeaturePlot(st_data, features ="Epcam" ,  pt.size.factor = 4 , ncol = 10)


## Add some gene program to the object
## (would be better to read a file...)
clust <- f_d_mathis_scigenex@gene_clusters[c(1,2,3,4,8,9,10,11)]

st_data <- AddModuleScore(st_data, features = clust)
st_data
SpatialFeaturePlot(st_data ,clust, pt.size.factor = 4)

## Some additional commands of scigenex
## that may improve diagrams


colnames(st_data@meta.data)[colnames(st_data@meta.data) == "Cluster1"] <- "Module Immature mTEC"
colnames(st_data@meta.data)[colnames(st_data@meta.data) == "Cluster2"] <- "Module Microfold mTEC"
colnames(st_data@meta.data)[colnames(st_data@meta.data) == "Cluster3"] <- "Module Aire+ mTEC"
colnames(st_data@meta.data)[colnames(st_data@meta.data) == "Cluster4"] <- "Module Keratinocytes mTEC"
colnames(st_data@meta.data)[colnames(st_data@meta.data) == "Cluster5"] <- "Module Neuroendocrine mTEC"
colnames(st_data@meta.data)[colnames(st_data@meta.data) == "Cluster6"] <- "Module Tufts2"
colnames(st_data@meta.data)[colnames(st_data@meta.data) == "Cluster7"] <- "Module Tufts1"
colnames(st_data@meta.data)[colnames(st_data@meta.data) == "Cluster8"] <- "Module Muscle mTEC"


hull <- display_hull(st_data, 
                     ident=ifelse(Seurat::Idents(st_data) %in% 0, 1, 0),
                     color = "white", size_x=10, size_y=10, 
                     step_y = 8, step_x = 8, delta = 3, hull_type = "wall", size = 1)

scigenex::plot_spatial(st_data, gene_name = "Fezf1", pt_size = 2.8 , stroke = 0 , size_title = 20 , face_title = "bold",colours = viridis::turbo(5)) + coord_fixed() + ggtitle("Tufts")+hull
scigenex::plot_spatial(st_data, gene_name = "Ccl21a", pt_size = 2.8 , stroke = 0 , size_title = 20 , face_title = "bold",colours = viridis::inferno(5)) + coord_fixed() + ggtitle("Immature mTEC")+hull
scigenex::plot_spatial(st_data, gene_name = "Aire", pt_size = 2.8 , stroke = 0 , size_title = 20 , face_title = "bold",colours = viridis::inferno(5)) + coord_fixed() + ggtitle("Aire+ mTEC")+hull
scigenex::plot_spatial(st_data, gene_name = "Spib", pt_size = 2.8, stroke = 0 , size_title = 20 , face_title = "bold",colours = viridis::inferno(5)) + coord_fixed() + ggtitle("Microfold TEC")+hull
scigenex::plot_spatial(st_data, gene_name = "Insm1", pt_size = 2.8 , stroke = 0 , size_title = 20 , face_title = "bold",colours = viridis::inferno(5)) + coord_fixed() + ggtitle("Neuroendocrine TEC")+hull
scigenex::plot_spatial(st_data, gene_name = "Grhl1", pt_size = 2.8, stroke = 0 , size_title = 20 , face_title = "bold",colours = viridis::inferno(5)) + coord_fixed() + ggtitle("Keratinocytes TEC")+hull
scigenex::plot_spatial(st_data, gene_name = "Myog", pt_size = 2.8, stroke = 0 , size_title = 20 , face_title = "bold",colours = viridis::inferno(5)) + coord_fixed() + ggtitle("Muscle TEC")+hull

scigenex::plot_spatial(st_data, metadata =  "Module Immature mTEC", pt_size = 3 , stroke = 0 , size_title = 20 , face_title = "bold",colours = viridis::turbo(5)) + coord_fixed() + ggtitle("Module Immature mTEC")+hull
scigenex::plot_spatial(st_data, metadata =  "Module Aire+ mTEC", pt_size = 3 , stroke = 0 , size_title = 20 , face_title = "bold",colours = viridis::turbo(5)) + coord_fixed() + ggtitle("Module Aire+ mTEC")+hull
scigenex::plot_spatial(st_data, metadata =  "Module Keratinocytes mTEC", pt_size = 3 , stroke = 0 , size_title = 20 , face_title = "bold",colours = viridis::turbo(5)) + coord_fixed() + ggtitle("Module Keratinocytes mTEC")+hull
scigenex::plot_spatial(st_data, metadata =  "Module Microfold mTEC", pt_size = 3 , stroke = 0 , size_title = 20 , face_title = "bold",colours = viridis::turbo(5)) + coord_fixed() + ggtitle("Module Microfold mTEC")+hull
scigenex::plot_spatial(st_data, metadata =  "Module Tufts1", pt_size = 3 , stroke = 0 , size_title = 20 , face_title = "bold",colours = viridis::turbo(5)) + coord_fixed() + ggtitle("Module Tufts1")+hull
scigenex::plot_spatial(st_data, metadata =  "Module Tufts2", pt_size = 2.8 , stroke = 0 , size_title = 20 , face_title = "bold",colours = viridis::turbo(5)) + coord_fixed() + ggtitle("Module Tufts2")+hull
scigenex::plot_spatial(st_data, metadata =  "Module Neuroendocrine mTEC", pt_size = 3 , stroke = 0 , size_title = 20 , face_title = "bold",colours = viridis::turbo(5)) + coord_fixed() + ggtitle("Module Neuroendocrine mTEC")+hull
scigenex::plot_spatial(st_data, metadata =  "Module Muscle mTEC", pt_size = 3 , stroke = 0 , size_title = 20 , face_title = "bold",colours = viridis::turbo(5)) + coord_fixed() + ggtitle("Module Muscle mTEC")+hull

scigenex::plot_spatial(st_data, gene_name = "Krt19", pt_size = 3 , stroke = 0 , size_title = 20 , face_title = "bold",colours = viridis::turbo(5)) + coord_fixed() + ggtitle("Muscle mTEC")+hull

subset_st_data <- st_data

subset_st_data<-RenameIdents(subset_st_data,"0"="Medulla")
subset_st_data<-RenameIdents(subset_st_data,"5"="Medulla")
subset_st_data<-RenameIdents(subset_st_data,"1"="Non-medulla")
subset_st_data<-RenameIdents(subset_st_data,"2"="Non-medulla")
subset_st_data<-RenameIdents(subset_st_data,"3"="Non-medulla")
subset_st_data<-RenameIdents(subset_st_data,"4"="Non-medulla")
subset_st_data<-RenameIdents(subset_st_data,"6"="Non-medulla")

table(Idents(subset_st_data))

VlnPlot(subset_st_data , "Aire+ mTEC")
VlnPlot(subset_st_data , "Keratinocytes mTEC")
VlnPlot(subset_st_data , "Microfold mTEC")
VlnPlot(subset_st_data , "Neuroendocrine mTEC")
VlnPlot(subset_st_data , "Tufts1")
VlnPlot(subset_st_data , "Tufts2")
VlnPlot(subset_st_data , "Muscle mTEC")

VlnPlot(subset_st_data , "Cited2")


# Normalizing module scores
# to get a unique legend

nclust <- length(clust)

for(i in 1:nclust){ 
  tmp <- st_data[[paste0("Cluster", i, sep="")]] 
  max_tmp <- max(tmp)
  min_tmp <- min(tmp)
  st_data[[paste0("Cluster", i, sep="")]]  <- (tmp[,1] - min(tmp))/(max_tmp - min_tmp)
}






scigenex::plot_spatial_panel(st_data, metadata = c("Module Immature mTEC","Module Microfold mTEC","Module Aire+ mTEC","Module Keratinocytes mTEC","Module Neuroendocrine mTEC", "Module Tufts1","Module Tufts2","Module Muscle mTEC"),colours = viridis::inferno(5),
                             ncol_layout = 4, 
                             pt_size=1, 
                             guide='collect',
                             stroke = 0,
                             size_title = 10, 
                             face_title = 'plain', 
                             barwidth = 0.25, barheight = 1.5 , panel_names = NULL) &hull & coord_fixed()

scigenex::plot_spatial_panel(st_data, metadata = paste("Cluster", 1:8, sep=""),colours = viridis::turbo(5),
                             ncol_layout = 4, 
                             pt_size=2.5, 
                             guide='collect',
                             stroke = 0,
                             size_title = 10, 
                             face_title = 'plain', 
                             barwidth = 0.25, barheight = 1.5 , panel_names =NULL) &hull & coord_fixed()


## Adding a hull around spot classified
## as class 1 (medulla most probably)
## By seurat

scigenex::plot_spatial(st_2, metadata =  "seurat_clusters", pt_size = 3.2 ) + coord_fixed() + hull


scigenex::plot_spatial(st_data, gene_name  = "Rgs21", pt_size = 3.5 , colours = viridis::turbo(5)) + coord_fixed()+hull



p1 <- scigenex::plot_spatial(st_data, gene_name = "Ccl19", pt_size = 3.5 , colours = viridis::viridis(5)) + coord_fixed() + hull
p2 <- scigenex::plot_spatial(st_data, metadata ="seurat_clusters", pt_size = 3.5) + coord_fixed() + hull
p1 +p2
saveRDS(st_data,"../st_data.rds")

a<-scigenex::plot_spatial(st_data, metadata =  "Module Immature mTEC", pt_size = 3 , stroke = 0 , size_title = 20 , face_title = "bold",colours = viridis::turbo(5)) + coord_fixed() + ggtitle("Module Immature mTEC")+hull
b<-scigenex::plot_spatial(st_data, metadata =  "Module Aire+ mTEC", pt_size = 3 , stroke = 0 , size_title = 20 , face_title = "bold",colours = viridis::turbo(5)) + coord_fixed() + ggtitle("Module Aire+ mTEC")+hull
a+b
