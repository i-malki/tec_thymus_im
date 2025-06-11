library(dplyr)
library(Seurat)
options(Seurat.object.assay.version = "v5")
library(patchwork)
library(RColorBrewer)
library(zoo)
library(scigenex)
library(ggplot2)

#devtools::install_github("dputhier/scigenex@17ead4863b8765bd68da3889f00a16d57dba9122")
#reload_pac()
setwd("/shared/projects/tec_thymus/tec_thymus/")

############## Datas ############## 
## Apply scigenex to Dm dataset

# Read Diane Mathis dataset processed with scigenex

d_mathis <- readRDS("singlecell/singlecelldata/mathis_processed_with_scigenex.rds")

# UMAP of datas

dimplot2 <- DimPlot(d_mathis, reduction = "umap", pt.size = 0.5 ) + NoLegend()
 UMAP<-LabelClusters(dimplot2, id = "ident", size = 3, repel = T,  box.padding = 1, box = TRUE, max.overlaps=100 )


########### Cell counts ################

# Extrait les identifiants 
ident_table <- table(Idents(d_mathis))

# Convertit la table en data.frame pour ggplot
ident_df <- as.data.frame(ident_table)

colnames(ident_df) <- c("Cell type", "Number of cells")
custom_colors <- colorRampPalette(brewer.pal(12, "Paired"))(23)

# Crée un barplot 
ident_df$`Cell type` <- factor(ident_df$`Cell type`, 
                               levels = ident_df$`Cell type`[order(-ident_df$`Number of cells`)])

# Graphique ggplot2
ggplot(ident_df, aes(x = `Cell type`, y = `Number of cells`, fill = `Cell type`)) +
  geom_bar(stat = "identity", position = "dodge") +  # Barres côte à côte
  labs(title = "Cell distributions across cell types", 
       x = NULL,  # Supprime le nom de l'axe X
       y = "Number of cells") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),  # Supprime les étiquettes de l'axe X
        axis.ticks.x = element_blank(),
        panel.border = element_rect(fill=NA , size = 1),
        panel.grid = element_blank()) +  # Supprime les ticks de l'axe X
  scale_fill_manual(values = custom_colors_named, guide = guide_legend(title = "Cell type"))# Légende triée



# Extrait les identifiants 
ident_table <- table(Idents(obj_query))

# Convertit la table en data.frame pour ggplot
ident_df <- as.data.frame(ident_table)

colnames(ident_df) <- c("Cell type", "Number of cells")
custom_colors <- colorRampPalette(brewer.pal(12, "Paired"))(23)

# Crée un barplot 
ident_df$`Cell type` <- factor(ident_df$`Cell type`, 
                               levels = ident_df$`Cell type`[order(-ident_df$`Number of cells`)])

# Graphique ggplot2
ggplot(ident_df, aes(x = `Cell type`, y = `Number of cells`, fill = `Cell type`)) +
  geom_bar(stat = "identity", position = "dodge") +  # Barres côte à côte
  labs(title = "Cell distributions across cell types", 
       x = NULL,  # Supprime le nom de l'axe X
       y = "Number of cells") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),  # Supprime les étiquettes de l'axe X
        axis.ticks.x = element_blank(),
        panel.border = element_rect(fill=NA , size = 1),
        panel.grid = element_blank()) +  # Supprime les ticks de l'axe X
  scale_fill_manual(values = custom_colors, guide = guide_legend(title = "Cell type"))  # Légende triée

############## Scigenex ####################
# Select normalized data for which genes
# have  sufficient count in raw data

# my processed datas
row_sum_raw <- rowSums(d_mathis@assays$RNA$counts)

d_mathis_mat <- as.matrix(d_mathis@assays$RNA$data)

d_mathis_mat <- d_mathis_mat[row_sum_raw>100,]



# Select genes that are well correlated with others
# and probably part of a cluster

d_mathis_mat_scigenex <- select_genes(d_mathis_mat, k = 85, noise_level = 0.75, 
                                      distance_method = "pearson", 
                                      no_anti_cor = TRUE)


# Construct a graph and partition it with MCL.



d_mathis_mat_scigenex <- gene_clustering(d_mathis_mat_scigenex, 
                                         method="reciprocal_neighborhood", 
                                         inflation = 1.6, 
                                         threads = 4)



# Filter the gene module by size.

f_d_mathis_scigenex <- scigenex::filter_cluster_size(d_mathis_mat_scigenex, min_cluster_size = 15)
f_d_mathis_scigenex2 <- scigenex::filter_cluster_size(d_mathis_mat_scigenex2, min_cluster_size = 15)

# Filter the gene module by std dev 
# (trying to restrict to those having most dispersion/variance)

f_d_mathis_scigenex <- scigenex::filter_cluster_sd(f_d_mathis_scigenex, min_sd = 0.1)
f_d_mathis_scigenex2 <- scigenex::filter_cluster_sd(f_d_mathis_scigenex2, min_sd = 0.2)

# clust stat

pdf( paste0("scigenex/outputs/figures/cluster_stats_mathis_processed-",".pdf" ), height=10, width=15 )
plot_cluster_stats(cluster_stats(f_d_mathis_scigenex)) + 
  ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                 axis.text.x = element_text(angle=45, vjust = 0.5),
                 panel.grid = element_blank())
dev.off()

# Compute gene that are the most representative 
# of each cluster
f_d_mathis_scigenex <-  top_genes(f_d_mathis_scigenex, top = 30)
f_d_mathis_scigenex2 <-  top_genes(f_d_mathis_scigenex2, top = 30)

# Plotting the most representative genes (subselect 30 cells from each cell pop)

s_f_d_mathis_mat_scigenex <- subsample_by_ident(f_d_mathis_scigenex, ident=Idents(d_mathis), nbcell=11)
pdf( paste0("scigenex/outputs/figures/clusters_heatmap_mathis_processed-",".pdf" ), height=15, width=20 )
plot_ggheatmap(f_d_mathis_scigenex[9], use_top_genes = TRUE, ident = Idents(d_mathis), hide_gene_name = FALSE)

s_f_d_mathis_mat_scigenex <-  top_genes(s_f_d_mathis_mat_scigenex, top = 5)
plot_ggheatmap(s_f_d_mathis_mat_scigenex[c(5,6,7,8,9,10,11,12)], use_top_genes = TRUE, ident = Idents(d_mathis), hide_gene_name = FALSE)


plot_profiles(s_f_d_mathis_mat_scigenex, ident = Idents(d_mathis))


##Top by go FT##

FT <- top_by_go(f_d_mathis_scigenex , go_id = "GO:0003677", species = "mmusculus", gene_id = "mgi_symbol")
s_FT <- subsample_by_ident(FT, ident=Idents(d_mathis), nbcell=30)
plot_ggheatmap(s_FT[9], use_top_genes = TRUE, ident = Idents(d_mathis), hide_gene_name = FALSE)

### Topbygo cytokines 

cyto <- top_by_go(f_d_mathis_scigenex , go_id = "GO:0005125", species = "mmusculus", gene_id = "mgi_symbol")
s_cyto <- subsample_by_ident(cyto, ident=Idents(d_mathis), nbcell=30)
plot_ggheatmap(s_cyto[1], use_top_genes = TRUE, ident = Idents(d_mathis), hide_gene_name = FALSE)

############## Enrichment of clusters ####

library(clusterProfiler)
library(org.Mm.eg.db)

# Cluster 2 (supposed microfold)

edo_2 <- enrichGO(gene = f_d_mathis_scigenex@gene_clusters[[2]],
                  OrgDb = org.Mm.eg.db,
                  keyType = "SYMBOL",  
                  ont = "BP",
                  pvalueCutoff = 0.5,
                  qvalueCutoff = 0.5)
cnetplot(edo)
edo_2_filtered <- edo_2@result %>%
  filter(grepl("\\bCcl6\\b", geneID))
edo_2_filtered <- edo_2_filtered %>%
  top_n(10, Count) %>%
  ungroup()
a<-ggplot(edo_2_filtered, aes(x = Count, y = reorder(Description, Count))) + 
  geom_bar(stat = "identity" , fill = "red" , color = "black") +
  labs(fill = "Ontology" , y= NULL , title = "Cluster 2") +
  theme(
    strip.text = element_blank(),   # Retirer les titres des facettes
    strip.background = element_blank(),  # Retirer le fond des titres de facettes
    panel.background = element_blank(),  # Retirer le fond du panneau (graphique)
    plot.background = element_blank(),
    axis.ticks.y = element_blank(),  # Retirer les ticks de l'axe Y
    axis.text.y = element_text(margin = margin(r = -10)),  # Réduire l'espacement entre les étiquettes et l'axe Y
    panel.grid.major.x = element_line(color = "grey"),
    panel.grid.minor.x = element_line(color = "grey")
  )

pdf( paste0("scigenex/Enrichment/Cluster2_microfold-",".pdf" ), height=3, width=5 )
a
dev.off()

# Cluster 4 (supposed keratinocytes) 
edo_4 <- enrichGO(gene = fdmgs@gene_clusters[[4]],
                  OrgDb = org.Mm.eg.db,
                  keyType = "SYMBOL",  
                  ont = "CC",
                  pvalueCutoff = 0.5,
                  qvalueCutoff = 0.5)

edo_4 <-  edo_4[edo_4@result$p.adjust <0.05,]
edo_4 <- edo_4 %>%
  top_n(10, Count) %>%
  ungroup()
b<-ggplot(edo_4, aes(x = Count, y = reorder(Description, Count))) + 
  geom_bar(stat = "identity" , fill = "blue" , color = "black") +
  labs(fill = "Ontology" , y= NULL , title = "Cluster 4") +
  theme(
    strip.text = element_blank(),   # Retirer les titres des facettes
    strip.background = element_blank(),  # Retirer le fond des titres de facettes
    panel.background = element_blank(),  # Retirer le fond du panneau (graphique)
    plot.background = element_blank(),
    axis.ticks.y = element_blank(),  # Retirer les ticks de l'axe Y
    axis.text.y = element_text(margin = margin(r = -10)),  # Réduire l'espacement entre les étiquettes et l'axe Y
    panel.grid.major.x = element_line(color = "grey"),
    panel.grid.minor.x = element_line(color = "grey")
  )


pdf( paste0("scigenex/Enrichment/Cluster4_keratinocytes-",".pdf" ), height=3, width=5 )
b
dev.off()

# Cluster 8  (supposed neuroendocrine) 
edo_8 <- enrichGO(gene = fdmgs@gene_clusters[[8]],
                  OrgDb = org.Mm.eg.db,
                  keyType = "SYMBOL",  
                  ont = "CC",
                  pvalueCutoff = 0.5,
                  qvalueCutoff = 0.5)

edo_8 <-  edo_8[edo_8@result$p.adjust <0.05,]
edo_8 <- edo_8 %>%
  top_n(10, Count) %>%
  ungroup()

c<-ggplot(edo_8, aes(x = Count, y = reorder(Description, Count))) + 
  geom_bar(stat = "identity" , fill = "yellow" , color = "black") +
  labs(fill = "Ontology" , y= NULL , title = "Cluster 8") +
  theme(
    strip.text = element_blank(),   # Retirer les titres des facettes
    strip.background = element_blank(),  # Retirer le fond des titres de facettes
    panel.background = element_blank(),  # Retirer le fond du panneau (graphique)
    plot.background = element_blank(),
    axis.ticks.y = element_blank(),  # Retirer les ticks de l'axe Y
    axis.text.y = element_text(margin = margin(r = -10)),  # Réduire l'espacement entre les étiquettes et l'axe Y
    panel.grid.major.x = element_line(color = "grey"),
    panel.grid.minor.x = element_line(color = "grey")
  )
pdf( paste0("scigenex/Enrichment/Cluster8_neuroendocrine-",".pdf" ), height=3, width=5 )
c
dev.off()

# Cluster 11 (supposed muscle) 
edo_11 <- enrichGO(gene = fdmgs@gene_clusters[[11]],
                   OrgDb = org.Mm.eg.db,
                   keyType = "SYMBOL",  
                   ont = "CC",
                   pvalueCutoff = 0.5,
                   qvalueCutoff = 0.5)

edo_11 <-  edo_11[edo_11@result$p.adjust <0.05,]
edo_11 <- edo_11 %>%
  top_n(10, Count) %>%
  ungroup()

d<-ggplot(edo_11, aes(x = Count, y = reorder(Description, Count))) + 
  geom_bar(stat = "identity" , fill = "green" , color = "black") +
  labs(fill = "Ontology" , y= NULL , title = "Cluster 11") +
  theme(
    strip.text = element_blank(),   # Retirer les titres des facettes
    strip.background = element_blank(),  # Retirer le fond des titres de facettes
    panel.background = element_blank(),  # Retirer le fond du panneau (graphique)
    plot.background = element_blank(),
    axis.ticks.y = element_blank(),  # Retirer les ticks de l'axe Y
    axis.text.y = element_text(margin = margin(r = -10)),  # Réduire l'espacement entre les étiquettes et l'axe Y
    panel.grid.major.x = element_line(color = "grey"),
    panel.grid.minor.x = element_line(color = "grey")
  )

pdf( paste0("scigenex/Enrichment/Cluster11_muscle-",".pdf" ), height=3, width=5 )
d
dev.off()

# Cluster  13 (supposed goblets) 
edo_13 <- enrichGO(gene = fdmgs@gene_clusters[[13]],
                   OrgDb = org.Mm.eg.db,
                   keyType = "SYMBOL",  
                   ont = "CC",
                   pvalueCutoff = 0.5,
                   qvalueCutoff = 0.5)
edo_13 <-  edo_13[edo_13@result$p.adjust <0.05,]
edo_13 <- edo_13 %>%
  top_n(10, Count) %>%
  ungroup()

e<-ggplot(edo_13, aes(x = Count, y = reorder(Description, Count))) + 
  geom_bar(stat = "identity" , fill = "purple" , color = "black") +
  labs(fill = "Ontology" , y= NULL , title = "Cluster 13") +
  theme(
    strip.text = element_blank(),   # Retirer les titres des facettes
    strip.background = element_blank(),  # Retirer le fond des titres de facettes
    panel.background = element_blank(),  # Retirer le fond du panneau (graphique)
    plot.background = element_blank(),
    axis.ticks.y = element_blank(),  # Retirer les ticks de l'axe Y
    axis.text.y = element_text(margin = margin(r = -10)),  # Réduire l'espacement entre les étiquettes et l'axe Y
    panel.grid.major.x = element_line(color = "grey"),
    panel.grid.minor.x = element_line(color = "grey")
  )


pdf( paste0("scigenex/Enrichment/Cluster13_goblets-",".pdf" ), height=5, width=5 )
e
dev.off()

a+b+c+d+e + plot_layout(ncol = 2 , nrow=3 )


############## Enrichment stats for up/downregulated genes in clusters ###############
# Perform enrichment stat
## Different genes deregulated in the Comp scripts
# Ripmova vs OTII

cmp_to_list <- scigenex::cmp_to_a_list(f_d_mathis_scigenex, rownames(Ripmova_induced_OTII) , background = intersect(rownames(d_mathis),rownames(Ripmova_vs_OTII)))
cluster_ripmova_induced_v_OTII <- rownames(cmp_to_list[[2]][[3]][["data"]])[cmp_to_list[[2]][[3]][["data"]]$jaccard >1.3010]

pdf( paste0("scigenex/outputs/figures/RipmOVAvsOTII/cmp_to_list_upregulated_ripmova_vs_OTII-",".pdf" ), height=10, width=15 )
cmp_to_list
dev.off()

cmp_to_list <- scigenex::cmp_to_a_list(f_d_mathis_scigenex, rownames(OTII_induced_Ripmova) , background = intersect(rownames(d_mathis),rownames(Ripmova_vs_OTII)))
cluster_OTII_induced_v_Ripmova <- rownames(cmp_to_list[[2]][[3]][["data"]])[cmp_to_list[[2]][[3]][["data"]]$jaccard >1.3010]

pdf( paste0("scigenex/outputs/figures/RipmOVAvsOTII/cmp_to_list_upregulated_OTII_vs_Ripmova-",".pdf" ), height=10, width=15 )
cmp_to_list
dev.off()

# Ripmova vs WT

cmp_to_list <- scigenex::cmp_to_a_list(f_d_mathis_scigenex, rownames(RipmOVA_induced_WT) , background = intersect(rownames(d_mathis),rownames(WT_RipmOVA)))
cluster_ripmova_induced_v_WT <- rownames(cmp_to_list[[2]][[3]][["data"]])[cmp_to_list[[2]][[3]][["data"]]$jaccard >1.3010]

pdf( paste0("scigenex/outputs/figures/RipmOVAvsWT/cmp_to_list_upregulated_ripmova_v_WT-",".pdf" ), height=10, width=15 )
cmp_to_list
dev.off()


cmp_to_list <- scigenex::cmp_to_a_list(f_d_mathis_scigenex, rownames(WT_induced_Ripmova) , background = intersect(rownames(d_mathis),rownames(WT_RipmOVA)))
cluster_WT_induced_v_Ripmova <- rownames(cmp_to_list[[2]][[3]][["data"]])[cmp_to_list[[2]][[3]][["data"]]$jaccard >1.3010]

pdf( paste0("scigenex/outputs/figures/RipmOVAvsWT/cmp_to_list_upregulated_WT_v_Ripmova-",".pdf" ), height=10, width=15 )
cmp_to_list
dev.off()

# WT vs OTII

cmp_to_list <- scigenex::cmp_to_a_list(f_d_mathis_scigenex, rownames(WT_induced_OTII) , background = intersect(rownames(d_mathis),rownames(OTII_WT)))
cluster_WT_induced_v_OTII <- rownames(cmp_to_list[[2]][[3]][["data"]])[cmp_to_list[[2]][[3]][["data"]]$jaccard >1.3010]

pdf( paste0("scigenex/outputs/figures/OTIIvsWT/cmp_to_list_upregulated_by_WT-",".pdf" ), height=10, width=15 )
cmp_to_list
dev.off()


cmp_to_list <- scigenex::cmp_to_a_list(f_d_mathis_scigenex, rownames(OTII_induced_WT) , background = intersect(rownames(d_mathis),rownames(OTII_WT)))
cluster_OTII_induced_v_WT <- rownames(cmp_to_list[[2]][[3]][["data"]])[cmp_to_list[[2]][[3]][["data"]]$jaccard >1.3010]

pdf( paste0("scigenex/outputs/figures/OTIIvsWT/cmp_to_list_upregulated_by_OTII-",".pdf" ), height=10, width=15 )
cmp_to_list
dev.off()

## Addmodule score

# Ripmova vs OTII

d_mathis_test<- AddModuleScore(d_mathis,f_d_mathis_scigenex@gene_clusters)

colnames(d_mathis_test@meta.data)[colnames(d_mathis_test@meta.data) == "Cluster1"] <- "Module immature mTEC"
colnames(d_mathis_test@meta.data)[colnames(d_mathis_test@meta.data) == "Cluster2"] <- "Module Microfold"
colnames(d_mathis_test@meta.data)[colnames(d_mathis_test@meta.data) == "Cluster3"] <- "Module Aire+ mTEC"
colnames(d_mathis_test@meta.data)[colnames(d_mathis_test@meta.data) == "Cluster4"] <- "Module Keratinocytes"
colnames(d_mathis_test@meta.data)[colnames(d_mathis_test@meta.data) == "Cluster8"] <- "Module Neuroendocrine"
colnames(d_mathis_test@meta.data)[colnames(d_mathis_test@meta.data) == "Cluster9"] <- "Module Tuft2"
colnames(d_mathis_test@meta.data)[colnames(d_mathis_test@meta.data) == "Cluster10"] <- "Module Tuft1"
colnames(d_mathis_test@meta.data)[colnames(d_mathis_test@meta.data) == "Cluster11"] <- "Module Muscle"
colnames(d_mathis_test@meta.data)[colnames(d_mathis_test@meta.data) == "Cluster21"] <- "Module Ptfa1+ pancreatic"


pdf(paste0("scigenex/outputs/figures/RipmOVAvsOTII/mapping_cluster_enrichies_upregulated_by_RipmOVA-",".pdf" ), height=30, width=20)
FeaturePlot(d_mathis_test,c("Module immature mTEC","Module Microfold","Module Aire+ mTEC","Module Keratinocytes","Module Neuroendocrine", "Module Muscle","Module Tuft1","Module Tuft2"), pt.size = 0.8 , label =TRUE , label.color = "black" , repel = TRUE , label.size = 3 , ncol =2)  & scale_color_viridis_c(option = "turbo")
dev.off()

pdf(paste0("scigenex/outputs/figures/RipmOVAvsOTII/mapping_cluster_enrichies_upregulated_by_OTII-",".pdf" ), height=30, width=20)
FeaturePlot(d_mathis_test, paste0("Cluster" ,cluster_OTII_induced_v_Ripmova), pt.size = 0.8 , label =TRUE , label.color = "black" , repel = TRUE , label.size = 3 , ncol =2) & scale_color_viridis_c(option = "turbo")
dev.off()

# Ripmova vs WT

pdf(paste0("scigenex/outputs/figures/RipmOVAvsWT/mapping_cluster_enrichies_upregulated_by_Ripmova-",".pdf" ), height=30, width=20)
FeaturePlot(d_mathis_test, paste0("Cluster" ,cluster_ripmova_induced_v_WT), pt.size = 0.8 , label =TRUE , label.color = "black" , repel = TRUE , label.size = 3 , ncol =2)  & scale_color_viridis_c(option = "turbo")
dev.off()


pdf(paste0("scigenex/outputs/figures/RipmOVAvsWT/mapping_cluster_enrichies_upregulated_by_WT-",".pdf" ), height=30, width=20)
FeaturePlot(d_mathis_test, paste0("Cluster" ,cluster_WT_induced_v_Ripmova), pt.size = 0.8 , label =TRUE , label.color = "black" , repel = TRUE , label.size = 3 , ncol =2)  & scale_color_viridis_c(option = "turbo")
dev.off()
# WT vs OTII

pdf(paste0("scigenex/outputs/figures/OTIIvsWT/mapping_cluster_enrichies_upregulated_by_OTII-",".pdf" ), height=30, width=20)
FeaturePlot(d_mathis_test, paste0("Cluster" ,cluster_WT_induced_v_OTII), pt.size = 0.8 , label =TRUE , label.color = "black" , repel = TRUE , label.size = 3 , ncol =2)  & scale_color_viridis_c(option = "turbo")
dev.off()


pdf(paste0("scigenex/outputs/figures/OTIIvsWT/mapping_cluster_enrichies_upregulated_by_OTII-",".pdf" ), height=30, width=20)
FeaturePlot(d_mathis_test, paste0("Cluster" ,cluster_OTII_induced_v_WT), pt.size = 0.8 , label =TRUE , label.color = "black" , repel = TRUE , label.size = 3 , ncol =2)  & scale_color_viridis_c(option = "turbo")
dev.off()

FeaturePlot(d_mathis_test, "Cluster9", pt.size = 0.8 , label =TRUE , label.color = "black" , repel = TRUE , label.size = 3 , ncol =2)  & scale_color_viridis_c(option = "turbo")




# top by go inters mathis_2_isma
top_interse<-top_by_intersect(f_d_mathis_scigenex, rownames(Ripmova_induced_OTII) )

# cytokines
top_interse_cytokines <-top_by_go(top_interse, go_id = "GO:0005125", species = "mmusculus", gene_id = "mgi_symbol")
sub_top_interse_cytokines <- subsample_by_ident(top_interse_cytokines, ident=Idents(d_mathis), nbcell=23)
pdf( paste0("scigenex/outputs/figures/cytokines",".pdf" ), height=10, width=15 )
plot_ggheatmap(sub_top_interse_cytokines, use_top_genes = TRUE, ident = Idents(d_mathis), hide_gene_name = FALSE)+theme(axis.text.y = element_text(face = "bold"))
dev.off()

# FT in clusters
top_interse_FT <-top_by_go(top_interse, go_id = "GO:0003677", species = "mmusculus", gene_id = "mgi_symbol")
sub_top_interse_FT <- subsample_by_ident(top_interse_FT, ident=Idents(d_mathis), nbcell=23)
pdf( paste0("scigenex/outputs/figures/FT",".pdf" ), height=10, width=15 )
plot_ggheatmap(sub_top_interse_FT[10], use_top_genes = TRUE, ident = Idents(d_mathis), hide_gene_name = FALSE)+theme(axis.text.y = element_text(face = "bold"))
dev.off()

#Keratins
top_interse_keratins <-top_by_go(top_interse, go_id = "GO:0031424", species = "mmusculus", gene_id = "mgi_symbol")
sub_top_interse_keratins <- subsample_by_ident(top_interse_keratins, ident=Idents(d_mathis), nbcell=23)
pdf( paste0("scigenex/outputs/figures/keratins",".pdf" ), height=10, width=15 )
plot_ggheatmap(sub_top_interse_keratins, use_top_genes = TRUE, ident = Idents(d_mathis), hide_gene_name = FALSE)+theme(axis.text.y = element_text(face = "bold"))
dev.off()

#Cell_surface
top_interse_cell_surface <-top_by_go(top_interse, go_id = "GO:0007166", species = "mmusculus", gene_id = "mgi_symbol")
sub_top_interse_cell_surface <- subsample_by_ident(top_interse_cell_surface, ident=Idents(d_mathis), nbcell=23)
pdf( paste0("scigenex/outputs/figures/cell_surface_receptor",".pdf" ), height=10, width=15 )
plot_ggheatmap(sub_top_interse_cell_surface, use_top_genes = TRUE, ident = Idents(d_mathis), hide_gene_name = FALSE)+theme(axis.text.y = element_text(face = "bold"))
dev.off()

#Topbyggo for enriched cluster upregulated genes in cytokines
top_interse@gene_clusters <- top_interse@top_genes

#not empty clusters
no_empty <- c("1","2","3","4","5","6","8","9","10","11","13","14","16","18","19","20","21","22","24","28","29")

#cytokines
cytokine_inters <-top_by_go(top_interse[no_empty], go_id = "GO:0005125", species = "mmusculus", gene_id = "mgi_symbol")
scytokine_inters <- subsample_by_ident(cytokine_inters, ident=Idents(d_mathis), nbcell=23)
pdf( paste0("scigenex/outputs/figures/RipmOVAvsOTII/cytokines_upregulated_by_ripmova-",".pdf" ), height=10, width=15 )
plot_ggheatmap(scytokine_inters, use_top_genes = TRUE, ident = Idents(d_mathis), hide_gene_name = FALSE)+theme(axis.text.y = element_text(face = "bold"))
dev.off()


#Topbyggo for enriched cluster upregulated genes in FT
FT_inters <- top_by_go(top_interse[no_empty], go_id = "GO:0003677", species = "mmusculus", gene_id = "mgi_symbol", host="https://www.ensembl.org")
sFT_inters <- subsample_by_ident(FT_inters, ident=Idents(d_mathis), nbcell=23)
pdf( paste0("scigenex/outputs/figures/RipmOVAvsOTII/TF_upregulated_by_ripmova-",".pdf" ), height=10, width=15 )
plot_ggheatmap(sFT_inters, use_top_genes = TRUE, ident = Idents(d_mathis), hide_gene_name = FALSE)+theme(axis.text.y = element_text(face = "bold"))
dev.off()



#Topbyggo for enriched cluster upregulated genes in keratins
kera_inters <- top_by_go(top_interse[no_empty], go_id = "GO:0031424", species = "mmusculus", gene_id = "mgi_symbol", host="https://www.ensembl.org")
s_kera_inters <- subsample_by_ident(kera_inters, ident=Idents(d_mathis), nbcell=23)
pdf( paste0("scigenex/outputs/figures/RipmOVAvsOTII/keratins_upregulated_by_ripmova-",".pdf" ), height=10, width=15 )
plot_ggheatmap(s_kera_inters, use_top_genes = TRUE, ident = Idents(d_mathis), hide_gene_name = FALSE)+theme(axis.text.y = element_text(face = "bold"))
dev.off()

#Topbyggo for enriched cluster upregulated genes in cell_surface_recep
CR_inters <- top_by_go(top_interse[no_empty], go_id = "GO:0007166", species = "mmusculus", gene_id = "mgi_symbol", host="https://www.ensembl.org")
s_CR_inters <- subsample_by_ident(CR_inters, ident=Idents(mathis_processed_with_scigenex), nbcell=23)
pdf( paste0("scigenex/outputs/figures/RipmOVAvsOTII/cell_surface_receptor_upregulated_by_ripmova-",".pdf" ), height=10, width=15 )
plot_ggheatmap(s_CR_inters, use_top_genes = TRUE, ident = Idents(mathis_processed_with_scigenex), hide_gene_name = FALSE)+theme(axis.text.y = element_text(face = "bold"))
dev.off()



############## Save object #####
saveRDS(f_d_mathis_mat_scigenex , "scigenex/f_d_mathis_scigenex.rds")


################Tissus specificty #####################
tabula_muris.R needed to be executed before

################ Modules by cell_type ################"

clust_to_cell <- setNames(rep(NA,length(f_d_mathis_scigenex@gene_clusters)) , names(f_d_mathis_scigenex@gene_clusters))

for (clust in 1:length(f_d_mathis_scigenex@gene_clusters)) {
  tmp <- f_d_mathis_scigenex@gene_clusters[[clust]]
  tmp <- names(sort(table(tscore_cell_max[tmp]) , decreasing = TRUE , na.last =TRUE))[1]
  clust_to_cell[clust] <- tmp
}
print(clust_to_cell)
clust_to_cell <- as.data.frame(clust_to_cell)
clust_to_cell$cluster <- as.numeric(rownames(clust_to_cell))
clust_to_cell$cluster <- factor(clust_to_cell$cluster, levels = as.character(sort(as.numeric(unique(clust_to_cell$cluster)))))

pdf( paste0("TRA_analysis/figures/modules_by_cell_type",".pdf" ), height=10, width=15 )
ggplot(clust_to_cell , aes(x = clust_to_cell, fill = cluster ))+geom_bar(position = "stack" ) +theme_bw() + labs(x = NULL)+
  geom_text(
    aes(label = cluster), 
    stat = "count", 
    position = position_stack(vjust = 0.5), 
    size = 5, 
    color = "black"
  )+theme(axis.text.x = element_text(face = "bold", angle = 45, hjust = 1) , legend.position = "none") +coord_flip()
dev.off()

################ Modules by Tissue ################"

clust_to_tissue <- setNames(rep(NA,length(f_d_mathis_scigenex@gene_clusters)) , names(f_d_mathis_scigenex@gene_clusters))

for (clust in 1:length(f_d_mathis_scigenex@gene_clusters)) {
  tmp <- f_d_mathis_scigenex@gene_clusters[[clust]]
  tmp <- names(sort(table(tscore_cell_max[tmp]) , decreasing = TRUE , na.last =TRUE))[1]
  clust_to_tissue[clust] <- tmp
}
print(clust_to_tissue)
clust_to_tissue <- as.data.frame(clust_to_tissue)
clust_to_tissue$cluster <- as.numeric(rownames(clust_to_tissue))
clust_to_tissue$cluster <- factor(clust_to_tissue$cluster, levels = as.character(sort(as.numeric(unique(clust_to_tissue$cluster)))))

pdf( paste0("TRA_analysis/figures/modules__by_tissue",".pdf" ), height=10, width=15 )
ggplot(clust_to_tissue , aes(x = clust_to_tissue, fill = cluster ))+geom_bar(position = "stack" ) +theme_bw() + labs(x = NULL)+
  geom_text(
    aes(label = cluster), 
    stat = "count", 
    position = position_stack(vjust = 0.5), 
    size = 5, 
    color = "black"
  )+theme(axis.text.x = element_text(face = "bold", angle = 45, hjust = 1) , legend.position = "none") +coord_flip()
dev.off()


##### Tscore by modules / tissues ######

#Microfold vs Tissues

df_tscore <- data.frame(ts_max = ts_max_tissue , cluster = "All genes")# Cree data frame de base ou chaque gene avec son score max dans le contexte total du génome

df_tscore_microfold <- data.frame(ts_max = ts_max_tissue[names(ts_max_tissue) %in% f_d_mathis_scigenex@gene_clusters[[2]]] , cluster =tscore_tissue_max[names(tscore_tissue_max) %in% f_d_mathis_scigenex@gene_clusters[[2]]]) #creer dataframe pour chaque gene sa valeur max et dans quel cluster on le retrouve 

df_tscore_microfold <- rbind(df_tscore_microfold ,df_tscore)

clust_order <- sort(tapply(df_tscore_microfold$ts_max, df_tscore_microfold$cluster, median , na.rm=TRUE) )
df_tscore_microfold$cluster <- factor(df_tscore_microfold$cluster , levels = names(clust_order),ordered = TRUE)

pdf( paste0("TRA_analysis/figures/Ripmova_vs_OTII/induced_by_Ripmova_vs_OTII_FC_by_Cell_type",".pdf" ), height=10, width=15)

Q1 <- quantile(df_tscore_microfold$ts_max, 0.25, na.rm = TRUE)
Q3 <- quantile(df_tscore_microfold$ts_max, 0.75, na.rm = TRUE)
IQR <- Q3 - Q1

# Calcul des limites des outliers (valeurs en dehors de 1.5 * IQR)
lower_limit <- Q1 - 1.5 * IQR
upper_limit <- Q3 + 1.5 * IQR

# Trouver la valeur maximale parmi toutes les boîtes sans les outliers
max_value <- max(df_tscore_microfold$ts_max[df_tscore_microfold$ts_max <= upper_limit])

# Création du boxplot avec une échelle dynamique basée sur les valeurs sans outliers
ggplot(df_tscore_microfold, aes(x = cluster, y = ts_max, fill = cluster)) + 
  geom_boxplot(
    outlier.shape = NA,  # Taille des outliers
    width = 0.8  # Largeur des boîtes
  ) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Ligne horizontale à 0
  theme_minimal() +
  theme(
    axis.ticks.x = element_blank(),
    legend.position = "none"
  ) + 
  labs(x = NULL, y = "Tscores in microfold module") + 
  coord_flip() +
  ylim(c(-3, max_value)) 
# facet_grid(~facet ,space='free' , scales = "fixed" )
dev.off()

#Muscle
df_tscore <- data.frame(ts_max = ts_max_tissue , cluster = "All genes")# Cree data frame de base ou chaque gene avec son score max dans le contexte total du génome

df_tscore_muscle <- data.frame(ts_max = ts_max_tissue[names(ts_max_tissue) %in% f_d_mathis_scigenex@gene_clusters[[4]]] , cluster =tscore_tissue_max[names(tscore_tissue_max) %in% f_d_mathis_scigenex@gene_clusters[[4]]]) #creer dataframe pour chaque gene sa valeur max et dans quel cluster on le retrouve 

df_tscore_muscle <- rbind(df_tscore_muscle ,df_tscore)

clust_order <- tapply(df_tscore_muscle$ts_max, df_tscore_muscle$cluster, function(x) quantile(x, 0.25, na.rm = TRUE))
df_tscore_muscle$cluster <- factor(df_tscore_muscle$cluster , levels = names(sort(clust_order)),ordered = TRUE)

pdf( paste0("TRA_analysis/figures/Ripmova_vs_OTII/induced_by_Ripmova_vs_OTII_FC_by_Cell_type",".pdf" ), height=10, width=15)

Q1 <- quantile(df_tscore_muscle$ts_max, 0.25, na.rm = TRUE)
Q3 <- quantile(df_tscore_muscle$ts_max, 0.75, na.rm = TRUE)
IQR <- Q3 - Q1

# Calcul des limites des outliers (valeurs en dehors de 1.5 * IQR)
lower_limit <- Q1 - 1.5 * IQR
upper_limit <- Q3 + 1.5 * IQR

# Trouver la valeur maximale parmi toutes les boîtes sans les outliers
max_value <- max(df_tscore_muscle$ts_max[df_tscore_muscle$ts_max <= upper_limit])

# Création du boxplot avec une échelle dynamique basée sur les valeurs sans outliers
ggplot(df_tscore_muscle, aes(x = cluster, y = ts_max, fill = cluster)) + 
  geom_boxplot(
    outlier.shape = NA,  # Taille des outliers
    width = 0.8  # Largeur des boîtes
  ) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Ligne horizontale à 0
  theme_minimal() +
  theme(
    axis.ticks.x = element_blank(),
    legend.position = "none"
  ) + 
  labs(x = NULL, y = "Tscores in microfold module") + 
  coord_flip() +
  ylim(c(-3, max_value)) 
# facet_grid(~facet ,space='free' , scales = "fixed" )
dev.off()


################ score specificity ################"
####  cell_modules / tissues ###




#### By cell_types###

# All genes + modules
df_cscore <- data.frame(ts_max = ts_max_cell , cluster = "All genes")  # Cree data frame de base ou chaque gene avec son score max dans le contexte total du génome
df_tscore_per_cluster <- data.frame(ts_max = ts_max_cell[names(scigenex::gene_cluster(f_d_mathis_scigenex))] , cluster =gene_cluster(f_d_mathis_scigenex)) #creer dataframe pour chaque gene sa valeur max et dans quel cluster on le retrouve 

# df_tscore_hollander <- data.frame(ts_max = ts_max[TRE_hollander] , cluster = "hollander TRA")
# df_tscore_WT_up <- data.frame(ts_max = ts_max[WT_upregulated] , cluster = "WT up vs RipmOVA")
# df_tscore_Ripmova_up <- data.frame(ts_max = ts_max[RipmOVA_upregulated] , cluster = "RipmOVA up vs WT")
# df_tscore_OTII_up_vs_WT <- data.frame(ts_max = ts_max[RipmOVA_upregulated] , cluster = "OTII up vs WT")
# df_tscore_WT_up_vs_OTII <- data.frame(ts_max = ts_max[WT_upregulated_vsOTII] , cluster = "WT up vs OTII")

df_cs_merge <- rbind(df_cscore , df_tscore_per_cluster)

# Order clusters by value of the median from top to flop
clust_order <- sort(tapply(df_cs_merge$ts_max , df_cs_merge$cluster, median , na.rm=TRUE))
df_cs_merge$cluster <- factor(df_cs_merge$cluster , levels = names(clust_order),ordered = TRUE)

pdf( paste0("TRA_analysis/figures/cell_type_score_spécificité_by_modules",".pdf" ), height=10, width=15 )
ggplot(df_cs_merge , aes(x = cluster , y = log10(ts_max) , fill = cluster)) + 
  geom_boxplot(outliers = FALSE) +  # Ajouter des lignes horizontales à 0, 1, -1
  theme_minimal() +
  theme(
    axis.ticks.x = element_blank(), 
    axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
    legend.position = "none", 
    strip.text = element_text(size = 12) ,panel.spacing = unit(1, "lines") # Ajuster la taille du texte des facettes
  ) + 
  labs(x = NULL , y = "cell_type specificity score (log10)") + coord_flip() 
# facet_grid(~facet ,space='free' , scales = "free")
dev.off()


# All genes + genes not in modules ( Ripmova vs OTII)

not_in_cluster_induced_by_ripmova <- rownames(aire_reg_down)[rownames(aire_reg_down) %in% unlist(f_d_mathis_scigenex@gene_clusters , use.names = FALSE)]
not_in_cluster_induced_by_OTII<- rownames(aire_reg_up)[rownames(aire_reg_up) %in% unlist(f_d_mathis_scigenex@gene_clusters , use.names = FALSE)]


df_tscore_not_in_cluster_ripmova <- data.frame(ts_max = ts_max_cell[not_in_cluster_induced_by_ripmova] , cluster = "induced by RipmOVA vs OTII and not in clusters")
df_tscore_not_in_cluster_OTII<- data.frame(ts_max = ts_max_cell[not_in_cluster_induced_by_OTII] , cluster = "induced by OTII vs RipmOVA and not in clusters")

df_cs_merge <- rbind(df_cscore ,df_tscore_not_in_cluster_ripmova ,df_tscore_not_in_cluster_OTII)

# Order clusters by value of the median from top to flop
clust_order <- sort(tapply(df_cs_merge$ts_max , df_cs_merge$cluster, median , na.rm=TRUE))
df_cs_merge$cluster <- factor(df_cs_merge$cluster , levels = names(clust_order),ordered = TRUE)

pdf( paste0("TRA_analysis/figures/Ripmova_vs_OTII/cell_type_score_spécificité_not_inclusters",".pdf" ), height=10, width=15 )
ggplot(df_cs_merge , aes(x = cluster , y = log10(ts_max) , fill = cluster)) + 
  geom_boxplot(outliers = FALSE) +  # Ajouter des lignes horizontales à 0, 1, -1
  theme_minimal() +
  theme(
    axis.ticks.x = element_blank(), 
    axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
    legend.position = "none", 
    strip.text = element_text(size = 12) ,panel.spacing = unit(1, "lines") # Ajuster la taille du texte des facettes
  ) + 
  labs(x = NULL , y = "Cell_type specificity score (log10)") + coord_flip() 
# facet_grid(~facet ,space='free' , scales = "free")
dev.off()

# All genes + up/downregulated by Ripmova vs OTII
df_cscore <- data.frame(ts_max = ts_max_cell , cluster = "All genes")  # Cree data frame de base ou chaque gene avec son score max dans le contexte total du génome
df_tscore_aire_down <- data.frame(ts_max = ts_max_cell[rownames(aire_reg_down)] , cluster = "induced by RipmOVA vs OTII")
df_tscore_aire_up <- data.frame(ts_max = ts_max_cell[rownames(aire_reg_up)] , cluster = "induced by OTII vs RipmOVA")
df_cs_merge <- rbind(df_cscore ,df_tscore_aire_down,df_tscore_aire_up )

# Order clusters by value of the median from top to flop
clust_order <- sort(tapply(df_cs_merge$ts_max , df_cs_merge$cluster, median , na.rm=TRUE))
df_cs_merge$cluster <- factor(df_cs_merge$cluster , levels = names(clust_order),ordered = TRUE)

pdf( paste0("TRA_analysis/figures/Ripmova_vs_OTII/Ripmova_OTII_up_down_regulated_cell_type_score",".pdf" ), height=10, width=15 )
ggplot(df_cs_merge , aes(x = cluster , y = log10(ts_max) , fill = cluster)) + 
  geom_boxplot(outliers = FALSE) +  # Ajouter des lignes horizontales à 0, 1, -1
  theme_minimal() +
  theme(
    axis.ticks.x = element_blank(), 
    axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
    legend.position = "none", 
    strip.text = element_text(size = 12) ,panel.spacing = unit(1, "lines") # Ajuster la taille du texte des facettes
  ) + 
  labs(x = NULL , y = "Cell_type specificity score (log10)") + coord_flip() 
# facet_grid(~facet ,space='free' , scales = "free")
dev.off()


# All genes + up/downregulated by Ripmova vs WT
df_cscore <- data.frame(ts_max = ts_max_cell , cluster = "All genes")  # Cree data frame de base ou chaque gene avec son score max dans le contexte total du génome
df_tscore_ripm <- data.frame(ts_max = ts_max_cell[rownames(RipmOVA_upregulated)] , cluster = "induced by RipmOVA vs WT")
df_tscore_WT <- data.frame(ts_max = ts_max_cell[rownames(WT_upregulated)] , cluster = "induced by WT vs RipmOVA")
df_cs_merge <- rbind(df_cscore ,df_tscore_ripm,df_tscore_WT )

# Order clusters by value of the median from top to flop
clust_order <- sort(tapply(df_cs_merge$ts_max , df_cs_merge$cluster, median , na.rm=TRUE))
df_cs_merge$cluster <- factor(df_cs_merge$cluster , levels = names(clust_order),ordered = TRUE)

pdf( paste0("TRA_analysis/figures/Ripmova_vs_WT/Ripmova_WT_up_down_regulated_cell_type_score",".pdf" ), height=10, width=15 )
ggplot(df_cs_merge , aes(x = cluster , y = log10(ts_max) , fill = cluster)) + 
  geom_boxplot(outliers = FALSE) +  # Ajouter des lignes horizontales à 0, 1, -1
  theme_minimal() +
  theme(
    axis.ticks.x = element_blank(), 
    axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
    legend.position = "none", 
    strip.text = element_text(size = 12) ,panel.spacing = unit(1, "lines") # Ajuster la taille du texte des facettes
  ) + 
  labs(x = NULL , y = "Cell_type specificity score (log10)") + coord_flip() 
# facet_grid(~facet ,space='free' , scales = "free")
dev.off()

# All genes + genes not in modules ( Ripmova vs WT)

not_in_cluster_induced_by_ripmova2 <- rownames(RipmOVA_upregulated)[rownames(RipmOVA_upregulated) %in% unlist(f_d_mathis_scigenex@gene_clusters , use.names = FALSE)]
not_in_cluster_induced_by_WT<- rownames(Ripmova)[rownames(WT_upregulated) %in% unlist(f_d_mathis_scigenex@gene_clusters , use.names = FALSE)]


df_tscore_not_in_cluster_ripmova2 <- data.frame(ts_max = ts_max_cell[not_in_cluster_induced_by_ripmova2] , cluster = "induced by Ripmova vs WT and not in clusters")
df_tscore_not_in_cluster_WT<- data.frame(ts_max = ts_max_cell[not_in_cluster_induced_by_WT] , cluster = " induced by WT vs Ripmova and not in clusters")

df_cs_merge <- rbind(df_cscore ,df_tscore_not_in_cluster_ripmova2 ,df_tscore_not_in_cluster_WT)

# Order clusters by value of the median from top to flop
clust_order <- sort(tapply(df_cs_merge$ts_max , df_cs_merge$cluster, median , na.rm=TRUE))
df_cs_merge$cluster <- factor(df_cs_merge$cluster , levels = names(clust_order),ordered = TRUE)

pdf( paste0("TRA_analysis/figures/Ripmova_vs_WT/cell_type_score_spécificité_not_inclusters",".pdf" ), height=10, width=15 )
ggplot(df_cs_merge , aes(x = cluster , y = log10(ts_max) , fill = cluster)) + 
  geom_boxplot(outliers = FALSE) +  # Ajouter des lignes horizontales à 0, 1, -1
  theme_minimal() +
  theme(
    axis.ticks.x = element_blank(), 
    axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
    legend.position = "none", 
    strip.text = element_text(size = 12) ,panel.spacing = unit(1, "lines") # Ajuster la taille du texte des facettes
  ) + 
  labs(x = NULL , y = "Cell_type specificity score (log10)") + coord_flip() 
# facet_grid(~facet ,space='free' , scales = "free")
dev.off()



#### By tissues ###


# All genes + modules
df_tscore <- data.frame(ts_max = ts_max_tissue , cluster = "All genes")  # Cree data frame de base ou chaque gene avec son score max dans le contexte total du génome
df_tscore_per_cluster <- data.frame(ts_max = ts_max_tissue[names(scigenex::gene_cluster(f_d_mathis_scigenex))] , cluster =gene_cluster(f_d_mathis_scigenex)) #creer dataframe pour chaque gene sa valeur max et dans quel cluster on le retrouve 

df_ts_merge <- rbind(df_tscore , df_tscore_per_cluster)

# Order clusters by value of the median from top to flop
clust_order <- sort(tapply(df_ts_merge$ts_max , df_ts_merge$cluster, median , na.rm=TRUE))
df_ts_merge$cluster <- factor(df_ts_merge$cluster , levels = names(clust_order),ordered = TRUE)

# Modifier la colonne facet selon les conditions données
df_ts_merge$facet <- ifelse(df_ts_merge$cluster == "All genes", "All genes", "Modules scigenex")

pdf( paste0("TRA_analysis/figures/tissue_score_spécificité_by_modules",".pdf" ), height=10, width=15 )
ggplot(df_ts_merge , aes(x = cluster , y = log10(ts_max) , fill = cluster)) + 
  geom_boxplot(outliers = FALSE) +  # Ajouter des lignes horizontales à 0, 1, -1
  theme_minimal() +
  theme(
    axis.ticks.x = element_blank(), 
    axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
    legend.position = "none", 
    strip.text = element_text(size = 12) ,panel.spacing = unit(1, "lines") # Ajuster la taille du texte des facettes
  ) + 
  labs(x = NULL , y = "Tissus specificity score (log10)") + coord_flip()
# facet_grid(~facet ,space='free' , scales = "free")
dev.off()

# All genes +genes not in modules (RipmOVA v OTII)

not_in_cluster_induced_by_ripmova <- rownames(Ripmova_induced_OTII)[rownames(Ripmova_induced_OTII) %in% unlist(f_d_mathis_scigenex@gene_clusters , use.names = FALSE)]
not_in_cluster_induced_by_OTII<- rownames(OTII_induced_Ripmova)[rownames(OTII_induced_Ripmova) %in% unlist(f_d_mathis_scigenex@gene_clusters , use.names = FALSE)]


df_tscore_not_in_cluster_ripmova <- data.frame(ts_max = ts_max_tissue[not_in_cluster_induced_by_ripmova] , cluster = "induced by Ripmova vs OTII and not in clusters ")
df_tscore_not_in_cluster_OTII<- data.frame(ts_max = ts_max_tissue[not_in_cluster_induced_by_OTII] , cluster = "induced by OTII vs RipmOVA and not in clusters")

df_ts_merge <- rbind(df_tscore ,df_tscore_not_in_cluster_ripmova ,df_tscore_not_in_cluster_OTII)

# Order clusters by value of the median from top to flop
clust_order <- sort(tapply(df_ts_merge$ts_max , df_ts_merge$cluster, median , na.rm=TRUE))
df_ts_merge$cluster <- factor(df_ts_merge$cluster , levels = names(clust_order),ordered = TRUE)



pdf( paste0("TRA_analysis/figures/Ripmova_vs_OTII/tissue_score_spécificité_not_inclusters",".pdf" ), height=10, width=15 )
ggplot(df_ts_merge , aes(x = cluster , y = log10(ts_max) , fill = cluster)) + 
  geom_boxplot(outliers = FALSE) +  # Ajouter des lignes horizontales à 0, 1, -1
  theme_minimal() +
  theme(
    axis.ticks.x = element_blank(), 
    axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
    legend.position = "none", 
    strip.text = element_text(size = 12) ,panel.spacing = unit(1, "lines") # Ajuster la taille du texte des facettes
  ) + 
  labs(x = NULL , y = "Tissue specificity score (log10)") + coord_flip() 
# facet_grid(~facet ,space='free' , scales = "free")
dev.off()


# All genes +genes not in modules (RipmOVA v WT)

not_in_cluster_induced_by_ripmova2 <- rownames(RipmOVA_induced_WT)[rownames(RipmOVA_induced_WT) %in% unlist(f_d_mathis_scigenex@gene_clusters , use.names = FALSE)]
not_in_cluster_induced_by_WT<- rownames(WT_induced_Ripmova)[rownames(WT_induced_Ripmova) %in% unlist(f_d_mathis_scigenex@gene_clusters , use.names = FALSE)]

df_tscore_not_in_cluster_ripmova2 <- data.frame(ts_max = ts_max_tissue[not_in_cluster_induced_by_ripmova2] , cluster = "induced by Ripmova vs WT and not in clusters")
df_tscore_not_in_cluster_WT<- data.frame(ts_max = ts_max_tissue[not_in_cluster_induced_by_WT] , cluster = "induced by WT vs RipmOVA and not in clusters")

df_ts_merge <- rbind(df_tscore ,df_tscore_not_in_cluster_ripmova2 ,df_tscore_not_in_cluster_WT)

# Order clusters by value of the median from top to flop
clust_order <- sort(tapply(df_ts_merge$ts_max , df_ts_merge$cluster, median , na.rm=TRUE))
df_ts_merge$cluster <- factor(df_ts_merge$cluster , levels = names(clust_order),ordered = TRUE)

pdf( paste0("TRA_analysis/figures/Ripmova_vs_WT/tissue_score_spécificité_not_inclusters",".pdf" ), height=10, width=15 )
ggplot(df_ts_merge , aes(x = cluster , y = log10(ts_max) , fill = cluster)) + 
  geom_boxplot(outliers = FALSE) +  # Ajouter des lignes horizontales à 0, 1, -1
  theme_minimal() +
  theme(
    axis.ticks.x = element_blank(), 
    axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
    legend.position = "none", 
    strip.text = element_text(size = 12) ,panel.spacing = unit(1, "lines") # Ajuster la taille du texte des facettes
  ) + 
  labs(x = NULL , y = "Tissue specificity score (log10)") + coord_flip() 
# facet_grid(~facet ,space='free' , scales = "free")
dev.off()

# All genes + up/downregulated by Ripmova vs OTII
df_tscore <- data.frame(ts_max = ts_max_tissue , cluster = "All genes")  # Cree data frame de base ou chaque gene avec son score max dans le contexte total du génome
df_tscore_aire_down <- data.frame(ts_max = ts_max_tissue[rownames(aire_reg_down)] , cluster = "induced by RipmOVA vs OTII")
df_tscore_aire_up <- data.frame(ts_max = ts_max_tissue[rownames(aire_reg_up)] , cluster = "induced by OTII vs WT")
df_ts_merge <- rbind(df_tscore ,df_tscore_aire_down,df_tscore_aire_up )

# Order clusters by value of the median from top to flop
clust_order <- sort(tapply(df_ts_merge$ts_max , df_ts_merge$cluster, median , na.rm=TRUE))
df_ts_merge$cluster <- factor(df_ts_merge$cluster , levels = names(clust_order),ordered = TRUE)

pdf( paste0("TRA_analysis/figures/Ripmova_vs_OTII/Ripmova_OTII_up_down_regulated_tissue_type_score",".pdf" ), height=10, width=15 )
ggplot(df_ts_merge , aes(x = cluster , y = log10(ts_max) , fill = cluster)) + 
  geom_boxplot(outliers = FALSE) +  # Ajouter des lignes horizontales à 0, 1, -1
  theme_minimal() +
  theme(
    axis.ticks.x = element_blank(), 
    axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
    legend.position = "none", 
    strip.text = element_text(size = 12) ,panel.spacing = unit(1, "lines") # Ajuster la taille du texte des facettes
  ) + 
  labs(x = NULL , y = "Tissue specificity score (log10)") + coord_flip() 
# facet_grid(~facet ,space='free' , scales = "free")
dev.off()


# All genes + up/downregulated by Ripmova vs WT
df_tscore <- data.frame(ts_max = ts_max_tissue , cluster = "All genes")  # Cree data frame de base ou chaque gene avec son score max dans le contexte total du génome
df_tscore_ripm <- data.frame(ts_max = ts_max_tissue[rownames(RipmOVA_upregulated)] , cluster = "induced by RipmOVA vs WT")
df_tscore_WT <- data.frame(ts_max = ts_max_tissue[rownames(WT_upregulated)] , cluster = "induced by WT vs RipmOVA")
df_ts_merge <- rbind(df_tscore ,df_tscore_ripm,df_tscore_WT )

# Order clusters by value of the median from top to flop
clust_order <- sort(tapply(df_ts_merge$ts_max , df_ts_merge$cluster, median , na.rm=TRUE))
df_ts_merge$cluster <- factor(df_ts_merge$cluster , levels = names(clust_order),ordered = TRUE)

pdf( paste0("TRA_analysis/figures/Ripmova_vs_WT/Ripmova_WT_up_down_regulated_tissue_type_score",".pdf" ), height=10, width=15 )
ggplot(df_ts_merge , aes(x = cluster , y = log10(ts_max) , fill = cluster)) + 
  geom_boxplot(outliers = FALSE) +  # Ajouter des lignes horizontales à 0, 1, -1
  theme_minimal() +
  theme(
    axis.ticks.x = element_blank(), 
    axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
    legend.position = "none", 
    strip.text = element_text(size = 12) ,panel.spacing = unit(1, "lines") # Ajuster la taille du texte des facettes
  ) + 
  labs(x = NULL , y = "Tissue specificity score (log10)") + coord_flip() 
# facet_grid(~facet ,space='free' , scales = "free")
dev.off()


############  ############         FOLD CHANGES      ############  ############  

############ FC by modules #############

# All genes (Ripmova v OTII) +modules 
all_genes <- aire_reg[names(scigenex::gene_cluster(f_d_mathis_scigenex)),]

df_FC <- data.frame(FC = aire_reg$DESeq2_log2FoldChange , cluster = "All genes")  # Cree data frame de base ou chaque gene avec son score max dans le contexte total du génome
df_FC_per_clust <- data.frame(FC = all_genes$DESeq2_log2FoldChange , cluster =gene_cluster(f_d_mathis_scigenex)) #creer dataframe pour chaque gene sa valeur max et dans quel cluster on le retrouve 

#Boxplot by modules
df_FC_merge <- rbind(df_FC , df_FC_per_clust )
df_FC_merge$facet <- ifelse(df_FC_merge$cluster == "All genes", "All genes", "Modules")
clust_order <- sort(tapply(df_FC_merge$FC , df_FC_merge$cluster, median , na.rm=TRUE) )
df_FC_merge$cluster <- factor(df_FC_merge$cluster , levels = names(clust_order),ordered = TRUE)

pdf( paste0("TRA_analysis/figures/Ripmova_vs_OTII/Boxplot_FC_par_modules",".pdf" ), height=10, width=15 )
ggplot(df_FC_merge , aes(x = cluster , y = FC , fill = cluster)) + 
  geom_boxplot(outlier.shape = NA ) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Ajouter des lignes horizontales à 0, 1, -1
  theme_minimal() +
  theme(
    axis.ticks.x = element_blank(),legend.position = "none"
  ) + 
  labs(x = NULL , y = "log2 FC") + coord_flip()
# facet_grid(~facet ,space='free' , scales = "free" )
dev.off()


# All genes + in / not in clusters RipmOVA vs OTII
not_in_cluster_induced_by_ripmova <- rownames(aire_reg_down)[rownames(aire_reg_down) %in% unlist(f_d_mathis_scigenex@gene_clusters , use.names = FALSE)]
not_in_cluster_induced_by_OTII<- rownames(aire_reg_up)[rownames(aire_reg_up) %in% unlist(f_d_mathis_scigenex@gene_clusters , use.names = FALSE)]

# not in cluster Ripmova
df_FC_not_in_clusters_ripmova <- data.frame(FC = aire_reg$DESeq2_log2FoldChange[rownames(aire_reg) %in% not_in_cluster_induced_by_ripmova] , cluster = "induced by RipmOVA vs OTII and not in cluster ")

# not in cluster OTII
df_FC_not_in_clusters_otii <- data.frame(FC = aire_reg$DESeq2_log2FoldChange[rownames(aire_reg) %in% not_in_cluster_induced_by_OTII] , cluster = "induced by OTII vs RipmOVA and not in cluster")

df_FC <- data.frame(FC = aire_reg$DESeq2_log2FoldChange , cluster = "All genes")  # Cree data frame de base ou chaque gene avec son score max dans le contexte total du génome

df_FC_merge <- rbind(df_FC , df_FC_not_in_clusters_ripmova , df_FC_not_in_clusters_otii )
clust_order <- sort(tapply(df_FC_merge$FC , df_FC_merge$cluster, median , na.rm=TRUE) )
df_FC_merge$cluster <- factor(df_FC_merge$cluster , levels = names(clust_order),ordered = TRUE)

pdf( paste0("TRA_analysis/figures/Ripmova_vs_OTII/Boxplot_FC_not_inclusters_Ripmova_OTII",".pdf" ), height=10, width=15 )
ggplot(df_FC_merge , aes(x = cluster , y = FC , fill = cluster)) + 
  geom_boxplot(outlier.shape = NA ) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Ajouter des lignes horizontales à 0, 1, -1
  theme_minimal() +
  theme(
    axis.ticks.x = element_blank(),legend.position = "none"
  ) + 
  labs(x = NULL , y = "log2 FC") + coord_flip()
# facet_grid(~facet ,space='free' , scales = "free" )
dev.off()

# All genes + in / not in clusters RipmOVA vs WT

not_in_cluster_induced_by_ripmova2 <- rownames(RipmOVA_upregulated)[rownames(RipmOVA_upregulated) %in% unlist(f_d_mathis_scigenex@gene_clusters , use.names = FALSE)]
not_in_cluster_induced_by_WT<- rownames(WT_upregulated)[rownames(WT_upregulated) %in% unlist(f_d_mathis_scigenex@gene_clusters , use.names = FALSE)]

# not in cluster Ripmova
df_FC_not_in_clusters_ripmova2 <- data.frame(FC = WT_RipmOVA$DESeq2_log2FoldChange[rownames(WT_RipmOVA) %in% not_in_cluster_induced_by_ripmova2] , cluster = "induced by RipmOVA vs WT and not in cluster")

# not in cluster WT
df_FC_not_in_clusters_WT <- data.frame(FC = WT_RipmOVA$DESeq2_log2FoldChange[rownames(WT_RipmOVA) %in% not_in_cluster_induced_by_WT] , cluster = "induced by WT vs RipmOVA and not in cluster")

df_FC <- data.frame(FC = WT_RipmOVA$DESeq2_log2FoldChange , cluster = "All genes")  # Cree data frame de base ou chaque gene avec son score max dans le contexte total du génome

df_FC_merge <- rbind(df_FC , df_FC_not_in_clusters_ripmova2 , df_FC_not_in_clusters_WT )
clust_order <- sort(tapply(df_FC_merge$FC , df_FC_merge$cluster, median , na.rm=TRUE) )
df_FC_merge$cluster <- factor(df_FC_merge$cluster , levels = names(clust_order),ordered = TRUE)

pdf( paste0("TRA_analysis/figures/Ripmova_vs_WT/Boxplot_FC_not_inclusters_RipmOVA_WT",".pdf" ), height=10, width=15 )
ggplot(df_FC_merge , aes(x = cluster , y = FC , fill = cluster)) + 
  geom_boxplot(outlier.shape = NA ) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Ajouter des lignes horizontales à 0, 1, -1
  theme_minimal() +
  theme(
    axis.ticks.x = element_blank(),legend.position = "none"
  ) + 
  labs(x = NULL , y = "log2 FC") + coord_flip()
# facet_grid(~facet ,space='free' , scales = "free" )
dev.off()

# All genes (Ripmova v OTII) +modules 
all_genes <- WT_RipmOVA[names(scigenex::gene_cluster(f_d_mathis_scigenex)),]

df_FC <- data.frame(FC = WT_RipmOVA$DESeq2_log2FoldChange , cluster = "All genes")  # Cree data frame de base ou chaque gene avec son score max dans le contexte total du génome
df_FC_per_clust <- data.frame(FC = all_genes$DESeq2_log2FoldChange , cluster =gene_cluster(f_d_mathis_scigenex)) #creer dataframe pour chaque gene sa valeur max et dans quel cluster on le retrouve 

#Boxplot by modules
df_FC_merge <- rbind(df_FC , df_FC_per_clust )
df_FC_merge$facet <- ifelse(df_FC_merge$cluster == "All genes", "All genes", "Modules")
clust_order <- sort(tapply(df_FC_merge$FC , df_FC_merge$cluster, median , na.rm=TRUE) )
df_FC_merge$cluster <- factor(df_FC_merge$cluster , levels = names(clust_order),ordered = TRUE)

pdf( paste0("TRA_analysis/figures/Ripmova_vs_WT/Boxplot_FC_par_modules",".pdf" ), height=10, width=15 )
ggplot(df_FC_merge , aes(x = cluster , y = FC , fill = cluster)) + 
  geom_boxplot(outlier.shape = NA ) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Ajouter des lignes horizontales à 0, 1, -1
  theme_minimal() +
  theme(
    axis.ticks.x = element_blank(),legend.position = "none"
  ) + 
  labs(x = NULL , y = "log2 FC") + coord_flip()
# facet_grid(~facet ,space='free' , scales = "free" )
dev.off()


######################## FC by tissues ############################

# induced by RipmOVA vs OTII

Tissues_genes_ripmova <- aire_reg_down[rownames(aire_reg_down) %in% names(tscore_tissue_max),]
df_FC_by_tissue <- data.frame(FC = Tissues_genes_ripmova$DESeq2_log2FoldChange , cluster = tscore_tissue_max[rownames(Tissues_genes_ripmova)])
df_FC <- data.frame(FC = aire_reg$DESeq2_log2FoldChange , cluster = "All genes")  # Cree data frame de base ou chaque gene avec son score max dans le contexte total du génome
df_FC_merge_tissue <- rbind(df_FC , df_FC_by_tissue )

df_FC_merge_tissue$facet <- ifelse(df_FC_merge_tissue$cluster == "All genes", "All genes", "tissues")
clust_order <- sort(tapply(df_FC_merge_tissue$FC , df_FC_merge_tissue$cluster, median , na.rm=TRUE) )
df_FC_merge_tissue$cluster <- factor(df_FC_merge_tissue$cluster , levels = names(clust_order),ordered = TRUE)


pdf( paste0("TRA_analysis/figures/Ripmova_vs_OTII/induced_by_Ripmova_vs_OTII_FC_by_tissue",".pdf" ), height=10, width=15 )
ggplot(df_FC_merge_tissue , aes(x = cluster , y = FC , fill = cluster)) + 
  geom_boxplot(outlier.shape = NA ) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Ajouter des lignes horizontales à 0, 1, -1
  theme_minimal() +
  theme(
    axis.ticks.x = element_blank(),legend.position = "none"
  ) + 
  labs(x = NULL , y = "log2 FC") + coord_flip()
  # facet_grid(~facet ,space='free' , scales = "fixed" )
dev.off()



#Induced by OTII vs Ripmova

Tissues_genes_otii <- aire_reg_up[rownames(aire_reg_up) %in% names(tscore_tissue_max),]
df_FC_by_tissue <- data.frame(FC = Tissues_genes_otii$DESeq2_log2FoldChange , cluster = tscore_tissue_max[rownames(Tissues_genes_otii)])
df_FC <- data.frame(FC = aire_reg$DESeq2_log2FoldChange , cluster = "All genes")  # Cree data frame de base ou chaque gene avec son score max dans le contexte total du génome

df_FC_merge_tissue <- rbind(df_FC , df_FC_by_tissue )

clust_order <- sort(tapply(df_FC_merge_tissue$FC , df_FC_merge_tissue$cluster, median , na.rm=TRUE) )
df_FC_merge_tissue$cluster <- factor(df_FC_merge_tissue$cluster , levels = names(clust_order),ordered = TRUE)

pdf( paste0("TRA_analysis/figures/Ripmova_vs_OTII/induced_by_OTII_vs_Ripmova_FC_by_tissue",".pdf" ), height=10, width=15 )
ggplot(df_FC_merge_tissue , aes(x = cluster , y = FC , fill = cluster)) + 
  geom_boxplot(outlier.shape = NA ) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Ajouter des lignes horizontales à 0, 1, -1
  theme_minimal() +
  theme(
    axis.ticks.x = element_blank(),legend.position = "none"
  ) + 
  labs(x = NULL , y = "log2 FC") + coord_flip()
# facet_grid(~facet ,space='free' , scales = "fixed" )
dev.off()

#Induced by WT vs Ripmova

Tissues_genes <- WT_upregulated[rownames(WT_upregulated) %in% names(tscore_tissue_max),]
df_FC_by_tissue <- data.frame(FC = Tissues_genes$DESeq2_log2FoldChange , cluster = tscore_tissue_max[rownames(Tissues_genes)])
df_FC <- data.frame(FC = WT_RipmOVA$DESeq2_log2FoldChange , cluster = "All genes")  # Cree data frame de base ou chaque gene avec son score max dans le contexte total du génome


df_FC_merge_tissue <- rbind(df_FC , df_FC_by_tissue )

clust_order <- sort(tapply(df_FC_merge_tissue$FC , df_FC_merge_tissue$cluster, median , na.rm=TRUE) )
df_FC_merge_tissue$cluster <- factor(df_FC_merge_tissue$cluster , levels = names(clust_order),ordered = TRUE)

pdf( paste0("TRA_analysis/figures/Ripmova_vs_WT/induced_by_WT_vs_Ripmova_FC_by_tissue",".pdf" ), height=10, width=15 )
ggplot(df_FC_merge_tissue , aes(x = cluster , y = FC , fill = cluster)) + 
  geom_boxplot(outlier.shape = NA ) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Ajouter des lignes horizontales à 0, 1, -1
  theme_minimal() +
  theme(
    axis.ticks.x = element_blank(),legend.position = "none"
  ) + 
  labs(x = NULL , y = "log2 FC") + coord_flip()
# facet_grid(~facet ,space='free' , scales = "fixed" )
dev.off()

#Induced by Ripmova vs WT

Tissues_genes <- RipmOVA_upregulated[rownames(RipmOVA_upregulated) %in% names(tscore_tissue_max),]
df_FC_by_tissue <- data.frame(FC = Tissues_genes$DESeq2_log2FoldChange , cluster = tscore_tissue_max[rownames(Tissues_genes)])
df_FC <- data.frame(FC = WT_RipmOVA$DESeq2_log2FoldChange , cluster = "All genes")  # Cree data frame de base ou chaque gene avec son score max dans le contexte total du génome

df_FC_merge_tissue <- rbind(df_FC , df_FC_by_tissue )

clust_order <- sort(tapply(df_FC_merge_tissue$FC , df_FC_merge_tissue$cluster, median , na.rm=TRUE) )
df_FC_merge_tissue$cluster <- factor(df_FC_merge_tissue$cluster , levels = names(clust_order),ordered = TRUE)

pdf( paste0("TRA_analysis/figures/Ripmova_vs_WT/induced_by_Ripmova_vs_WT_FC_by_tissue",".pdf" ), height=10, width=15 )
ggplot(df_FC_merge_tissue , aes(x = cluster , y = FC , fill = cluster)) + 
  geom_boxplot(outlier.shape = NA ) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Ajouter des lignes horizontales à 0, 1, -1
  theme_minimal() +
  theme(
    axis.ticks.x = element_blank(),legend.position = "none"
  ) + 
  labs(x = NULL , y = "log2 FC") + coord_flip()
# facet_grid(~facet ,space='free' , scales = "fixed" )
dev.off()



######################## FC by cell_type ########################

# Induced by Ripmova vs OTII

cell_genes_ripmova <- aire_reg_down[rownames(aire_reg_down) %in% names(tscore_cell_max),]
df_FC_by_cell <- data.frame(FC = cell_genes_ripmova$DESeq2_log2FoldChange , cluster = tscore_cell_max[rownames(cell_genes_ripmova)])
df_FC <- data.frame(FC = aire_reg$DESeq2_log2FoldChange , cluster = "All genes")  # Cree data frame de base ou chaque gene avec son score max dans le contexte total du génome

df_FC_merge_cell_type <- rbind(df_FC , df_FC_by_cell )

clust_order <- sort(tapply(df_FC_merge_cell_type$FC , df_FC_merge_cell_type$cluster, median , na.rm=TRUE) )
df_FC_merge_cell_type$cluster <- factor(df_FC_merge_cell_type$cluster , levels = names(clust_order),ordered = TRUE)

pdf( paste0("TRA_analysis/figures/Ripmova_vs_OTII/induced_by_Ripmova_vs_OTII_FC_by_Cell_type",".pdf" ), height=10, width=15)
ggplot(df_FC_merge_cell_type , aes(x = cluster , y = FC , fill = cluster)) + 
  geom_boxplot(outlier.shape = NA ) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Ajouter des lignes horizontales à 0, 1, -1
  theme_minimal() +
  theme(
    axis.ticks.x = element_blank(),legend.position = "none"
  ) + 
  labs(x = NULL , y = "log2 FC") + coord_flip()
# facet_grid(~facet ,space='free' , scales = "fixed" )
dev.off()


# Induced by OTII vs Ripmova

cell_genes_otii <- aire_reg_up[rownames(aire_reg_up) %in% names(tscore_cell_max),]
df_FC_by_cell <- data.frame(FC = cell_genes_otii$DESeq2_log2FoldChange , cluster = tscore_cell_max[rownames(cell_genes_otii)])
df_FC <- data.frame(FC = aire_reg$DESeq2_log2FoldChange , cluster = "All genes")  # Cree data frame de base ou chaque gene avec son score max dans le contexte total du génome

df_FC_merge_cell_type <- rbind(df_FC , df_FC_by_cell )

clust_order <- sort(tapply(df_FC_merge_cell_type$FC , df_FC_merge_cell_type$cluster, median , na.rm=TRUE) )
df_FC_merge_cell_type$cluster <- factor(df_FC_merge_cell_type$cluster , levels = names(clust_order),ordered = TRUE)

pdf( paste0("TRA_analysis/figures/Ripmova_vs_OTII/induced_by_OTII_vs_Ripmova_FC_by_Cell_type",".pdf" ), height=10, width=15)
ggplot(df_FC_merge_cell_type , aes(x = cluster , y = FC , fill = cluster)) + 
  geom_boxplot(outlier.shape = NA ) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Ajouter des lignes horizontales à 0, 1, -1
  theme_minimal() +
  theme(
    axis.ticks.x = element_blank(),legend.position = "none"
  ) + 
  labs(x = NULL , y = "log2 FC") + coord_flip()
# facet_grid(~facet ,space='free' , scales = "fixed" )
dev.off()



# Induced by WT vs Ripmova

cell_genes <- WT_upregulated[rownames(WT_upregulated) %in% names(tscore_cell_max),]
df_FC_by_cell <- data.frame(FC = cell_genes$DESeq2_log2FoldChange , cluster = tscore_cell_max[rownames(cell_genes)])
df_FC <- data.frame(FC = WT_RipmOVA$DESeq2_log2FoldChange , cluster = "All genes")  # Cree data frame de base ou chaque gene avec son score max dans le contexte total du génome


df_FC_merge_cell_type <- rbind(df_FC , df_FC_by_cell )

clust_order <- sort(tapply(df_FC_merge_cell_type$FC , df_FC_merge_cell_type$cluster, median , na.rm=TRUE) )
df_FC_merge_cell_type$cluster <- factor(df_FC_merge_cell_type$cluster , levels = names(clust_order),ordered = TRUE)

pdf( paste0("TRA_analysis/figures/Ripmova_vs_WT//induced_by_WT_vs_Ripmova_FC_by_Cell_type",".pdf" ), height=10, width=15)
ggplot(df_FC_merge_cell_type , aes(x = cluster , y = FC , fill = cluster)) + 
  geom_boxplot(outlier.shape = NA ) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Ajouter des lignes horizontales à 0, 1, -1
  theme_minimal() +
  theme(
    axis.ticks.x = element_blank(),legend.position = "none"
  ) + 
  labs(x = NULL , y = "log2 FC") + coord_flip()
# facet_grid(~facet ,space='free' , scales = "fixed" )
dev.off()



# Induced by Ripmova vs WT

cell_genes <- RipmOVA_upregulated[rownames(RipmOVA_upregulated) %in% names(tscore_cell_max),]
df_FC_by_cell <- data.frame(FC = cell_genes$DESeq2_log2FoldChange , cluster = tscore_cell_max[rownames(cell_genes)])
df_FC <- data.frame(FC = WT_RipmOVA$DESeq2_log2FoldChange , cluster = "All genes")  # Cree data frame de base ou chaque gene avec son score max dans le contexte total du génome

df_FC_merge_cell_type <- rbind(df_FC , df_FC_by_cell )

clust_order <- sort(tapply(df_FC_merge_cell_type$FC , df_FC_merge_cell_type$cluster, median , na.rm=TRUE) )
df_FC_merge_cell_type$cluster <- factor(df_FC_merge_cell_type$cluster , levels = names(clust_order),ordered = TRUE)

pdf( paste0("TRA_analysis/figures/Ripmova_vs_WT//induced_by_Ripmova_vs_WT_FC_by_Cell_type",".pdf" ), height=10, width=15)
ggplot(df_FC_merge_cell_type , aes(x = cluster , y = FC , fill = cluster)) + 
  geom_boxplot(outlier.shape = NA ) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Ajouter des lignes horizontales à 0, 1, -1
  theme_minimal() +
  theme(
    axis.ticks.x = element_blank(),legend.position = "none"
  ) + 
  labs(x = NULL , y = "log2 FC") + coord_flip()
# facet_grid(~facet ,space='free' , scales = "fixed" )
dev.off()

########## FC by mimetic ########## 


# All genes (Ripmova v OTII) + modules 



all_genes <- Ripmova_vs_OTII[names(scigenex::gene_cluster(f_d_mathis_scigenex , cluster = c(1,2,3,4,8,9,10,11))),]
df_FC <- data.frame(FC = Ripmova_vs_OTII$DESeq2_log2FoldChange , cluster = "All genes")  # Cree data frame de base ou chaque gene avec son score max dans le contexte total du génome
df_FC_per_mimetic <- data.frame(FC = all_genes$DESeq2_log2FoldChange , cluster =gene_cluster(f_d_mathis_scigenex , cluster = c(1,2,3,4,8,9,10,11))) #creer dataframe pour chaque gene sa valeur max et dans quel cluster on le retrouve 
df_FC_per_mimetic$cluster[df_FC_per_mimetic$cluster == "1"] <- "Module Immature mTEC"
df_FC_per_mimetic$cluster[df_FC_per_mimetic$cluster == "2"] <- "Module Microfold mTEC"
df_FC_per_mimetic$cluster[df_FC_per_mimetic$cluster == "3"] <- "Module Aire+ mTEC"
df_FC_per_mimetic$cluster[df_FC_per_mimetic$cluster == "4"] <- "Module Keratinocytes mTEC"
df_FC_per_mimetic$cluster[df_FC_per_mimetic$cluster == "8"] <- "Module Neuroendocrine mTEC"
df_FC_per_mimetic$cluster[df_FC_per_mimetic$cluster == "9"] <- "Module Tuft 2 mTEC"
df_FC_per_mimetic$cluster[df_FC_per_mimetic$cluster == "10"] <- "Module Tuft 1 mTEC"
df_FC_per_mimetic$cluster[df_FC_per_mimetic$cluster == "11"] <- "Module Muscle mTEC"


df_FC_merge_mimetic <- rbind(df_FC , df_FC_per_mimetic )
clust_order <- sort(tapply(df_FC_merge_mimetic$FC , df_FC_merge_mimetic$cluster, median , na.rm=TRUE) )
df_FC_merge_mimetic$cluster <- factor(df_FC_merge_mimetic$cluster , levels = names(clust_order),ordered = TRUE)

pdf( paste0("TRA_analysis/figures/Ripmova_vs_OTII/Ripmova_OTII_FC_by_mimetics",".pdf" ), height=10, width=15)
ggplot(df_FC_merge_mimetic , aes(x = cluster , y = FC , fill = cluster)) + 
  geom_boxplot(outlier.shape = NA ) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Ajouter des lignes horizontales à 0, 1, -1
  theme_minimal() +
  theme(
    axis.ticks.x = element_blank(),legend.position = "none"
  ) + 
  labs(x = NULL , y = "log2 FC") + coord_flip() +ggtitle("RipmOVA_vs_OTII")
dev.off()



# All genes (Ripmova v WT) + modules 

WT_RipmOVA_genes <- WT_RipmOVA[names(scigenex::gene_cluster(f_d_mathis_scigenex , cluster = c(2,4,8,9,10,11,13,21))),]
df_FC <- data.frame(FC = WT_RipmOVA$DESeq2_log2FoldChange , cluster = "All genes")  # Cree data frame de base ou chaque gene avec son score max dans le contexte total du génome
df_FC_per_mimetic <- data.frame(FC = WT_RipmOVA_genes$DESeq2_log2FoldChange , cluster =gene_cluster(f_d_mathis_scigenex , cluster = c(2,4,8,9,10,11,13,21))) #creer dataframe pour chaque gene sa valeur max et dans quel cluster on le retrouve 
df_FC_per_mimetic$cluster[df_FC_per_mimetic$cluster == "2"] <- "Microfold"
df_FC_per_mimetic$cluster[df_FC_per_mimetic$cluster == "4"] <- "Keratinocytes / skin"
df_FC_per_mimetic$cluster[df_FC_per_mimetic$cluster == "8"] <- "Neuroendocrine"
df_FC_per_mimetic$cluster[df_FC_per_mimetic$cluster == "9"] <- "Tuft 2"
df_FC_per_mimetic$cluster[df_FC_per_mimetic$cluster == "10"] <- "Tuft1"
df_FC_per_mimetic$cluster[df_FC_per_mimetic$cluster == "11"] <- "Muscle"
df_FC_per_mimetic$cluster[df_FC_per_mimetic$cluster == "21"] <- "Ptf1a+ pancreatic"
df_FC_per_mimetic$cluster[df_FC_per_mimetic$cluster == "13"] <- "Goblets/secretory/basal lung"


df_FC_merge_mimetic <- rbind(df_FC , df_FC_per_mimetic )
clust_order <- sort(tapply(df_FC_merge_mimetic$FC , df_FC_merge_mimetic$cluster, median , na.rm=TRUE) )
df_FC_merge_mimetic$cluster <- factor(df_FC_merge_mimetic$cluster , levels = names(clust_order),ordered = TRUE)

pdf( paste0("TRA_analysis/figures/Ripmova_vs_WT/Ripmova_WT_FC_by_mimetics",".pdf" ), height=10, width=15)
ggplot(df_FC_merge_mimetic , aes(x = cluster , y = FC , fill = cluster)) + 
  geom_boxplot(outlier.shape = NA ) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Ajouter des lignes horizontales à 0, 1, -1
  theme_minimal() +
  theme(
    axis.ticks.x = element_blank(),legend.position = "none"
  ) + 
  labs(x = NULL , y = "log2 FC") + coord_flip()+ggtitle("RipmOVA_vs_WT")
dev.off()


# All genes + OTII v WT


WT_OTII_genes <- OTII_WT[names(scigenex::gene_cluster(f_d_mathis_scigenex , cluster = c(2,4,8,9,10,11,13,21))),]
df_FC <- data.frame(FC = OTII_WT$DESeq2_log2FoldChange , cluster = "All genes")  # Cree data frame de base ou chaque gene avec son score max dans le contexte total du génome
df_FC_per_mimetic <- data.frame(FC = WT_OTII_genes$DESeq2_log2FoldChange , cluster =gene_cluster(f_d_mathis_scigenex , cluster = c(2,4,8,9,10,11,13,21))) #creer dataframe pour chaque gene sa valeur max et dans quel cluster on le retrouve 
df_FC_per_mimetic$cluster[df_FC_per_mimetic$cluster == "2"] <- "Microfold"
df_FC_per_mimetic$cluster[df_FC_per_mimetic$cluster == "4"] <- "Keratinocytes / skin"
df_FC_per_mimetic$cluster[df_FC_per_mimetic$cluster == "8"] <- "Neuroendocrine"
df_FC_per_mimetic$cluster[df_FC_per_mimetic$cluster == "9"] <- "Tuft 2"
df_FC_per_mimetic$cluster[df_FC_per_mimetic$cluster == "10"] <- "Tuft1"
df_FC_per_mimetic$cluster[df_FC_per_mimetic$cluster == "11"] <- "Muscle"
df_FC_per_mimetic$cluster[df_FC_per_mimetic$cluster == "21"] <- "Ptf1a+ pancreatic"
df_FC_per_mimetic$cluster[df_FC_per_mimetic$cluster == "13"] <- "Goblets/secretory/basal lung"


df_FC_merge_mimetic <- rbind(df_FC , df_FC_per_mimetic )
clust_order <- sort(tapply(df_FC_merge_mimetic$FC , df_FC_merge_mimetic$cluster, median , na.rm=TRUE) )
df_FC_merge_mimetic$cluster <- factor(df_FC_merge_mimetic$cluster , levels = names(clust_order),ordered = TRUE)

pdf( paste0("TRA_analysis/figures/OTII_vs_WT/OTII_WT_FC_by_mimetics",".pdf" ), height=10, width=15)
ggplot(df_FC_merge_mimetic , aes(x = cluster , y = FC , fill = cluster)) + 
  geom_boxplot(outlier.shape = NA ) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Ajouter des lignes horizontales à 0, 1, -1
  theme_minimal() +
  theme(
    axis.ticks.x = element_blank(),legend.position = "none"
  ) + 
  labs(x = NULL , y = "log2 FC") + coord_flip()+ggtitle("WT_vs_OTII")
dev.off()

#####################

# 1. Générer le FeaturePlot
pdf( paste0("scigenex/outputs/figures/RipmOVAvsOTII/mapping_genes_not_in_cluster_induced_by_Ripmova",".pdf" ), height=10, width=15)
d_mathis_test<-AddModuleScore(d_mathis , not_in_cluster_induced_by_ripmova , name = "notinclusters_induced_by_Ripmova_vs_OTII")
FeaturePlot(d_mathis_test , "notinclusters_induced_by_Ripmova_vs_OTII1", pt.size = 0.2, repel = TRUE , label.size = 3 ,label = TRUE  )+ggtitle("Ripmova induced vs OTII not in clusters")  & scale_color_viridis_c(option = "turbo")
dev.off()

pdf( paste0("scigenex/outputs/figures/RipmOVAvsOTII/mapping_genes_not_in_cluster_induced_by_OTII",".pdf" ), height=10, width=15)
d_mathis_test<-AddModuleScore(d_mathis , not_in_cluster_induced_by_OTII , name = "notinclusters_induced_by_OTII_vsRipmova")
FeaturePlot(d_mathis_test , "notinclusters_induced_by_OTII_vsRipmova1", pt.size = 0.2 , label =TRUE , label.color = "black" , repel = TRUE , label.size = 3)+ggtitle("OTII induced vs RipmOVA not in clusters")  & scale_color_viridis_c(option = "turbo")
dev.off()

pdf( paste0("scigenex/outputs/figures/RipmOVAvsWT/mapping_genes_not_in_cluster_induced_by_WT",".pdf" ), height=10, width=15)
d_mathis_test<-AddModuleScore(d_mathis , not_in_cluster_induced_by_WT, name = "notinclusters_induced_by_WT_vs_Ripmova")
FeaturePlot(d_mathis_test , "notinclusters_induced_by_WT_vs_Ripmova1", pt.size = 0.2 , label =TRUE , label.color = "black" , repel = TRUE , label.size = 3)+ggtitle("WT induced vs RipmOVA not in clusters")  & scale_color_viridis_c(option = "turbo")
dev.off()

pdf( paste0("scigenex/outputs/figures/RipmOVAvsWT/mapping_genes_not_in_cluster_induced_by_Ripmova",".pdf" ), height=10, width=15)
d_mathis_test<-AddModuleScore(d_mathis , not_in_cluster_induced_by_ripmova2, name = "notinclusters_induced_by_Ripmova_vs_WT")
FeaturePlot(d_mathis_test , "notinclusters_induced_by_Ripmova_vs_WT1", pt.size = 0.2 , label =TRUE , label.color = "black" , repel = TRUE , label.size = 3)+ggtitle("Ripmova induced vs WT not in clusters")  & scale_color_viridis_c(option = "turbo")
dev.off()
