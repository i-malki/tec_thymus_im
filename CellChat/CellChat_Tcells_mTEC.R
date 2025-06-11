setwd("/shared/projects/tec_thymus/tec_analysis/scripts/")
#devtools::install_github("sqjin/CellChat" )
library(CellChat)
library(Seurat)
library(reticulate) # from h5ad object
options(Seurat.object.assay.version = "v5")
library(circlize)
library(patchwork)
options(stringAsFactors = FALSE)


### ### ## 1 Create CellChat object ### ### ##

# From merged_seurat



seurat_merged <- readRDS("mTEC_T_cells_seurat_object.rds")
data.input_merged <- GetAssayData(seurat_merged , assay = "RNA" , layer = "data")

labels_merged <- Idents(seurat_merged)

meta_merged<- data.frame(group = labels_merged, row.names = names(labels_merged))

cell_chat_merged <- createCellChat(object= seurat_merged , meta = meta_merged , group.by = "group")


data.input_subset <- GetAssayData(subset_seurat_11 , assay = "RNA" , layer = "data")

labels_merged_subset <- Idents(subset_seurat_11)

meta_merged_subset<- data.frame(group = labels_merged_subset, row.names = names(labels_merged_subset))

cell_chat_merged_subset <- createCellChat(object= data.input_subset , meta = meta_merged_subset , group.by = "group")

# From 5had object ( voir technique conversion h5ad en seurat)


### ### ## 2 Set and update CellChatDB ### ### ##

# Load ligand-receptor DB
Cellchat.mouse <- CellChatDB.mouse

# Show the categories
showDatabaseCategory(Cellchat.mouse)

#Set the DB to CellChat 
CellChatDB.use <- Cellchat.mouse

# Or subset the categories we cant 

CellChatDB.use <- subsetDB(Cellchat.mouse , search = "Secreted Signaling")
CellChatDB.use <- subsetDB(Cellchat.mouse , search = "ECM-Receptor")
CellChatDB.use <- subsetDB(Cellchat.mouse , search = "Cell-Cell Contact")

# Update CellChatDB with our DB 
# Add DB to CellChat object 
cell_chat_merged@DB <- CellChatDB.mouse
cell_chat_merged_subset@DB <- CellChatDB.mouse

### ### ##  3 Inférence of Cell-cell network ## ## ## ## ## 

# Subset the cellchat data for saving computation cost 
cell_chat_merged <- subsetData(cell_chat_merged)
cell_chat_merged_subset <- subsetData(cell_chat_merged_subset)

# Pre-processing the expression data 

cell_chat_merged <- identifyOverExpressedGenes(cell_chat_merged)
cell_chat_merged <- identifyOverExpressedInteractions(cell_chat_merged)

cell_chat_merged_subset <- identifyOverExpressedGenes(cell_chat_merged_subset)
cell_chat_merged_subset <- identifyOverExpressedInteractions(cell_chat_merged_subset)

# Project gene expression in protein-protein interaction 
cell_chat_merged <- projectData(cell_chat_merged , PPI.mouse)
cell_chat_merged_subset <- projectData(cell_chat_merged_subset , PPI.mouse)

# Compute the communication probability and infer cellular network 
#cell_chat <- computeCommunProb(cell_chat) # No projected data 
cell_chat_merged <- computeCommunProb(cell_chat_merged, raw.use = FALSE) # Use the projected data
cell_chat_merged_subset <- computeCommunProb(cell_chat_merged_subset, raw.use = FALSE) # Use the projected data

# Filter the cell-cell communication datas , here we only take if min 10 cells interacts 
cell_chat_merged <- filterCommunication(cell_chat_merged , min.cells = 10)
cell_chat_merged_subset <- filterCommunication(cell_chat_merged_subset , min.cells = 10)

# Infer communication at a pathway 

cell_chat_merged <- computeCommunProbPathway(cell_chat_merged)
cell_chat_merged_subset <- computeCommunProbPathway(cell_chat_merged_subset)

# Calculate the aggregated cell-cell communication network 

cell_chat_merged <- aggregateNet(cell_chat_merged)
cell_chat_merged@net$count
cell_chat_merged@net$weight

cell_chat_merged_subset <- aggregateNet(cell_chat_merged_subset)
cell_chat_merged_subset@net$count
cell_chat_merged_subset@net$weight

# Visualize the aggregated cell-cell communication network 
groupSize_cell_chat_merged <- as.numeric(table(cell_chat_merged@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cell_chat_merged@net$count , vertex.weight = groupSize_cell_chat_merged ,
                 weight.scale = TRUE , label.edge = FALSE , title.name = "Number of Interactions")
netVisual_circle(cell_chat_merged@net$weight, vertex.weight = groupSize_cell_chat_merged ,
                 weight.scale = TRUE , label.edge = FALSE , title.name = "Interaction weight/strenght" )
dev.off()

groupSize_cell_chat_merged_subset <- as.numeric(table(cell_chat_merged_subset@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cell_chat_merged_subset@net$count , vertex.weight = groupSize_cell_chat_merged_subset ,
                 weight.scale = TRUE , label.edge = FALSE , title.name = "Number of Interactions")
netVisual_circle(cell_chat_merged_subset@net$weight, vertex.weight = groupSize_cell_chat_merged_subset ,
                 weight.scale = TRUE , label.edge = FALSE , title.name = "Interaction weight/strenght" )
dev.off()
# TO examine the signaling from each group separatly on the same plot
cell_chat_merged_mat_weight_subset <- cell_chat_merged_subset@net$weight
par(mfrow = c(2,8), xpd=TRUE)
for (i in 1:nrow(cell_chat_merged_mat_weight_subset)) {
  mat_2 <- matrix(0 , nrow = nrow(cell_chat_merged_mat_weight_subset) , ncol = ncol(cell_chat_merged_mat_weight_subset) , dimnames = dimnames(cell_chat_merged_mat_weight_subset))
  mat_2[i, ] <- cell_chat_merged_mat_weight_subset[i, ]
  netVisual_circle(mat_2, vertex.weight = groupSize_cell_chat_merged_subset ,
                   weight.scale = TRUE , label.edge = FALSE , title.name = "Interaction weight/strenght" )
}

dev.off()



cell_chat_merged_mat_weight <- cell_chat_merged@net$weight
par(mfrow = c(2,8), xpd=TRUE)
for (i in 1:nrow(cell_chat_merged_mat_weight)) {
  mat_2 <- matrix(0 , nrow = nrow(cell_chat_merged_mat_weight) , ncol = ncol(cell_chat_merged_mat_weight) , dimnames = dimnames(cell_chat_merged_mat_weight))
  mat_2[i, ] <- cell_chat_merged_mat_weight[i, ]
  netVisual_circle(mat_2, vertex.weight = groupSize_cell_chat_merged ,
                   weight.scale = TRUE , label.edge = FALSE , title.name = "Interaction weight/strenght" )
}

dev.off()
### ### ##  4 visualization of Cell-cell network ## ## ## ## ## 

extractEnrichedLR(cell_chat_merged , signaling = c(cell_chat_merged@netP[["pathways"]]),
                  geneLR.return = TRUE)

extractEnrichedLR(cell_chat_merged_subset , signaling = c(cell_chat_merged_subset@netP[["pathways"]]),
                  geneLR.return = TRUE)

## Visualize the contribution of each pair in the network 

cell_chat_merged@netP[["pathways"]]


netAnalysis_contribution(cell_chat_merged , signaling = c(cell_chat_merged@netP[["pathways"]]) ,
                         title = "Contribution of each LR pair" , thresh = 0.05)

cell_chat_merged_subset@netP[["pathways"]]


netAnalysis_contribution(cell_chat_merged_subset , signaling = c(cell_chat_merged_subset@netP[["pathways"]]) ,
                         title = "Contribution of each LR pair")



# Extract a pathway in particular ex:CMH-I
extractEnrichedLR(cell_chat_merged, signaling = "TENASCIN", geneLR.return = FALSE)

netAnalysis_contribution(cell_chat_merged, signaling = "LCK")


extractEnrichedLR(cell_chat_merged_subset, signaling = "MHC-I", geneLR.return = TRUE)

netAnalysis_contribution(cell_chat_merged_subset, signaling = "MHC-II")
## Circle plot 
par(mfrow = c(1,2), xpd=TRUE)
par()
netVisual_aggregate(cell_chat_merged, signaling = "LCK", layout = "circle")
netVisual_aggregate(cell_chat_merged, signaling = "MHC-II", layout = "circle")
dev.off()

par(mfrow = c(1,2), xpd=TRUE)
par()
netVisual_aggregate(cell_chat_merged_subset, signaling = "MHC-I", layout = "circle")
netVisual_aggregate(cell_chat_merged_subset, signaling = "MHC-II", layout = "circle")
dev.off()

par(mfrow = c(4,4) , xpd= TRUE)
netVisual_individual(cell_chat_merged, signaling = "APP", layout = "circle")
dev.off()

## Chord PLot
par(mfrow = c(1, 1), xpd=TRUE)
par(cex = 0.5)
netVisual_aggregate(cell_chat_merged, signaling = "LCK", layout = "chord" , vertex.label.cex = 2)
netVisual_chord_cell (cell_chat_merged, signaling = "LCK" , lab.cex = 1.5)


netVisual_chord_gene (cell_chat_merged, signaling = "TENASCIN" , lab.cex = 1.5)
dev.off()


par(mfrow = c(1, 2), xpd=TRUE)
par(cex = 0.5)
netVisual_aggregate(cell_chat_merged_subset, signaling = "MHC-I", layout = "chord" , vertex.label.cex = 2)
netVisual_chord_cell (cell_chat_merged_subset, signaling = "MHC-II" , lab.cex = 1.5)
netVisual_chord_gene (cell_chat_merged_subset, signaling = "MHC-II" , lab.cex = 1.5)
dev.off()

# for all interactions 
netVisual_chord_cell_internal(cell_chat_merged_mat_weight)

netVisual_chord_cell_internal(cell_chat_merged_mat_weight_subset)
dev.off()
## define source and targets 
# sources.use = la source de l'interaction , targets.use = les target que l'on veut afficher , to subset the visualization to specific cells for all the interactions between 
# here as exexmple microfolf againt intertypical and early progenitor 
netVisual_chord_gene(cell_chat_merged, sources.use = c(3,4,5,6), targets.use = c(1,2), 
                     lab.cex = 0.5,legend.pos.y = 30)
netVisual_chord_gene(cell_chat_merged, sources.use = 2, targets.use = 9, 
                     lab.cex = 0.5,legend.pos.y = 30)
netVisual_chord_gene(cell_chat_merged, sources.use = c(1,2), targets.use = 9,
                     lab.cex = 0.5, legend.pos.x = 15)


netVisual_chord_gene(cell_chat_merged_subset, sources.use = c(3,4,5,6), targets.use = c(1,2), 
                     lab.cex = 0.5,legend.pos.y = 30)
netVisual_chord_gene(cell_chat_merged_subset, sources.use = 2, targets.use = 9, 
                     lab.cex = 0.5,legend.pos.y = 30)
netVisual_chord_gene(cell_chat_merged_subset, sources.use = c(1,2), targets.use = 9,
                     lab.cex = 0.5, legend.pos.x = 15)

## Hierarchical plot
par(mfrow = c(1, 1), xpd=TRUE)

vertex.receiver_cell_chat_merged = seq(1,2) # define the left portion of cell groups  

netVisual_aggregate(cell_chat_merged, signaling = "LCK", 
                    vertex.receiver = vertex.receiver_cell_chat_merged, layout = "hierarchy")

netVisual_hierarchy1(cell_chat_merged_mat_weight , vertex.receiver =vertex.receiver_cell_chat_merged )


par(mfrow = c(1, 1), xpd=TRUE)

vertex.receiver_cell_chat_merged_subset = seq(1,2) # define the left portion of cell groups  

netVisual_aggregate(cell_chat_merged_subset, signaling = "MHC-II", 
                    vertex.receiver = vertex.receiver_cell_chat_merged_subset, layout = "hierarchy")

netVisual_hierarchy1(cell_chat_merged_mat_weight_subset , vertex.receiver =vertex.receiver_cell_chat_merged_subset )

## Heatmap

# For a pathway
netVisual_heatmap(cell_chat_merged, signaling = "RANKL", color.heatmap = "Reds")
netVisual_heatmap(cell_chat_merged, signaling = "MHC-II", color.heatmap = "Reds")


netVisual_heatmap(cell_chat_merged_subset, signaling = "RANKL", color.heatmap = "Reds")
netVisual_heatmap(cell_chat_merged_subset, signaling = "MHC-II", color.heatmap = "Reds",  font.size = 8)

#Heatmap global
netVisual_heatmap(cell_chat_merged, color.heatmap = "Reds")

netVisual_heatmap(cell_chat_merged_subset, color.heatmap = "Reds" ,  font.size = 8  )

## Violin Plot

plotGeneExpression(cell_chat_merged, signaling = "RANKL")
plotGeneExpression(cell_chat_merged_subset , signaling = "RANKL")




plotGeneExpression(cell_chat_merged, signaling = "MHC-II")
plotGeneExpression(cell_chat_merged_subset , signaling = "MHC-I")
## Bubbleplot

netVisual_bubble(cell_chat_merged, sources.use = 2, targets.use = c(1:6), 
                 remove.isolate = FALSE) 


netVisual_bubble(cell_chat_merged_subset, sources.use = 2, targets.use = c(1:6), 
                 remove.isolate = FALSE) 


### ### ##  5 Systematic analysis of cell-cell communication networks  ## ## ## ## ## 

library(NMF)
library(ggalluvial)

# 1. Compute the network centrality scores
cell_chat_merged <- netAnalysis_computeCentrality(cell_chat_merged, slot.name = "netP")

cell_chat_merged_subset <- netAnalysis_computeCentrality(cell_chat_merged_subset, slot.name = "netP")

# Scatter plot to visualize aggregated communication networks for each cell type
netAnalysis_signalingRole_scatter(cell_chat_merged) # all signaling pathways


netAnalysis_signalingRole_scatter(cell_chat_merged_subset) # all signaling pathways

# Scatter plot to Visualize selected communication networks
netAnalysis_signalingRole_scatter(cell_chat_merged, signaling = "MHC-II")

netAnalysis_signalingRole_scatter(cell_chat_merged_subset, signaling = "MHC-II")

# Heatmap to visualize dominant cell types for each signaling pathway
netAnalysis_signalingRole_heatmap(cell_chat_merged, pattern = "outgoing", height = 20)
netAnalysis_signalingRole_heatmap(cell_chat_merged, pattern = "incoming", height = 20)

netAnalysis_signalingRole_heatmap(cell_chat_merged_subset, pattern = "outgoing", height = 20 , width = 20)
netAnalysis_signalingRole_heatmap(cell_chat_merged_subset, pattern = "incoming", height = 20 , width = 20)
# Visualize selected outgoing/incoming signals and contributing cell types
netAnalysis_signalingRole_heatmap(cell_chat_merged, pattern = "outgoing",
                                  signaling = c("MHC-II", "MHC-I"))
netAnalysis_signalingRole_heatmap(cell_chat_merged, pattern = "incoming",
                                  signaling = c("MHC-II", "MHC-I"))


netAnalysis_signalingRole_heatmap(cell_chat_merged_subset, pattern = "outgoing",
                                  signaling = c("MHC-II", "MHC-I"))
netAnalysis_signalingRole_heatmap(cell_chat_merged_subset, pattern = "incoming",
                                  signaling = c("MHC-II", "MHC-I"))
# Heatmap to visualize major signaling roles of different cell groups
netAnalysis_signalingRole_network(cell_chat_merged, signaling = "MHC-II", width = 10, 
                                  height = 5, font.size = 10)


netAnalysis_signalingRole_network(cell_chat_merged_subset, signaling = "MHC-II", width = 10, 
                                  height = 5, font.size = 10)
# 2. Identify global communication patterns to explore how multiple cell types 
# and signaling pathways coordinate

# Identify and visualize outgoing communication pattern of secreting cells
selectK(cell_chat_merged, pattern = "outgoing") # infer the number of patterns, NMF
nPatterns = 3 # a suitable number of patterns is the one begin to drop suddenly.
cell_chat_merged <- identifyCommunicationPatterns(cell_chat_merged, pattern = "outgoing",
                                                k = nPatterns, width = 5, height = 9)

netAnalysis_river(cell_chat_merged, pattern = "outgoing") # river plot
netAnalysis_dot(cell_chat_park, pattern = "outgoing") # dot plot



selectK(cell_chat_merged_subset, pattern = "outgoing") # infer the number of patterns, NMF
nPatterns = 3 # a suitable number of patterns is the one begin to drop suddenly.
cell_chat_merged_subset <- identifyCommunicationPatterns(cell_chat_merged_subset, pattern = "outgoing",
                                                  k = nPatterns, width = 5, height = 9)

netAnalysis_river(cell_chat_merged_subset, pattern = "outgoing" , font.size = 6 , font.size.title = 15) # river plot
netAnalysis_dot(cell_chat_merged_subset, pattern = "outgoing") # dot plot

## Identify and visualize incoming communication pattern of target cells
selectK(cell_chat_park, pattern = "incoming")
nPatterns = 3
cell_chat_merged <- identifyCommunicationPatterns(cell_chat_merged,pattern = "incoming", 
                                                k = nPatterns, width = 5, height = 9)

netAnalysis_river(cell_chat_merged, pattern = "incoming") # river plot
netAnalysis_dot(cell_chat_merged, pattern = "incoming") # dot plot


selectK(cell_chat_merged_subset, pattern = "incoming")
nPatterns = 3
cell_chat_merged_subset <- identifyCommunicationPatterns(cell_chat_merged_subset,pattern = "incoming", 
                                                  k = nPatterns, width = 5, height = 9)

netAnalysis_river(cell_chat_merged_subset, pattern = "incoming", font.size = 6 , font.size.title = 15) # river plot
netAnalysis_dot(cell_chat_merged_subset, pattern = "incoming") # dot plot

# 3. Groups signaling pathways based on their functional/structural similarities
# Identify signaling groups based on functional similarity

cell_chat_merged <- computeNetSimilarity(cell_chat_merged, type = "functional")
cell_chat_merged <- netEmbedding(cell_chat_merged, type = "functional")
cell_chat_merged <- netClustering (cell_chat_merged, type = "functional", do.parallel = FALSE)


cell_chat_merged_subset <- computeNetSimilarity(cell_chat_merged_subset, type = "functional")
cell_chat_merged_subset <- netEmbedding(cell_chat_merged_subset, type = "functional")
cell_chat_merged_subset <- netClustering (cell_chat_merged_subset, type = "functional", do.parallel = FALSE)

# Visualization in 2D-space
netVisual_embedding(cell_chat_merged, type = "functional", label.size = 3.5)
netVisual_embeddingZoomIn(cell_chat_merged, type = "functional", nCol = 2)

netVisual_embedding(cell_chat_merged_subset, type = "functional", label.size = 3.5)
netVisual_embeddingZoomIn(cell_chat_merged_subset, type = "functional", nCol = 2)


# Identify signaling groups based on structure similarity
# multimeric ligand-receptor complexes, soluble agonists and antagonists, 
# stimulatory and inhibitory co-ligands and co-receptors
cell_chat_merged <- computeNetSimilarity(cell_chat_merged, type = "structural")
cell_chat_merged <- netEmbedding(cell_chat_merged, type = "structural")
cell_chat_merged <- netClustering(cell_chat_merged, type = "structural",do.parallel = FALSE)


cell_chat_merged_subset <- computeNetSimilarity(cell_chat_merged_subset, type = "structural")
cell_chat_merged_subset <- netEmbedding(cell_chat_merged_subset, type = "structural")
cell_chat_merged_subset <- netClustering(cell_chat_merged_subset, type = "structural",do.parallel = FALSE)

# Visualization in 2D-space
netVisual_embedding(cell_chat_merged, type = "structural", label.size = 3.5)
netVisual_embeddingZoomIn(cell_chat_merged, type = "structural", nCol = 2)

netVisual_embedding(cell_chat_merged_subset, type = "structural", label.size = 3.5)
netVisual_embeddingZoomIn(cell_chat_merged_subset, type = "structural", nCol = 2)
