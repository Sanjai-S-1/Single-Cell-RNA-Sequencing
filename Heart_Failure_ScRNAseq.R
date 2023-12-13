#Heart failure Sc-RNAseq

#Importing Libraries

library(tidyverse)
library(Seurat)
library(Matrix)
library(celldex)
library(scater)
library(SingleR)
library(rafalib)  #mypar
library(cowplot) #plotgrid
library(gridExtra)
library(ggplot2)


#Loading Data and Making a Seurat Object

# get data location
dirs <- list.dirs(path = 'data/', recursive = F, full.names = F)

for(x in dirs){
  name <- x
  
  cts <- ReadMtx(mtx = paste0('data/',x,'/matrix.mtx.gz'),
                 features = paste0('data/',x,'/features.tsv.gz'),
                 cells = paste0('data/',x,'/barcodes.tsv.gz'))
  
  # create seurat objects
  assign(name, CreateSeuratObject(counts = cts))
}

# merge datasets

merged_seurat <- merge(CTL, HF,
                       add.cell.ids = c("Control", "Heart Failure" ),
                       project = 'Heart_Failure')


# QC & filtering -----------------------

View(merged_seurat@meta.data)

# create a sample column
merged_seurat$sample <- rownames(merged_seurat@meta.data)

# split sample column
merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'sample', into = c('Condition', 'Barcode'), 
                                    sep = '_')

# calculate mitochondrial percentage
merged_seurat$mitoPercent <- PercentageFeatureSet(merged_seurat, pattern='^MT-')

# explore QC

VlnPlot(merged_seurat, c("nCount_RNA", "nFeature_RNA", "mitoPercent"), pt.size = 0.1)

# filtering
merged_seurat_filtered <- subset(merged_seurat, subset = nCount_RNA > 1000 &
                                   nCount_RNA < 20000 &
                                   nFeature_RNA > 1000 &
                                   mitoPercent < 10)
remove("merged_seurat")

merged_seurat_filtered

VlnPlot(merged_seurat_filtered, c("nCount_RNA", "nFeature_RNA", "mitoPercent"), pt.size = 0.1)

merged_seurat

# perform standard workflow steps to figure out if we see any batch effects --------
merged_seurat_filtered <- NormalizeData(object = merged_seurat_filtered)
merged_seurat_filtered <- FindVariableFeatures(object = merged_seurat_filtered)
merged_seurat_filtered <- ScaleData(object = merged_seurat_filtered)
merged_seurat_filtered <- RunPCA(object = merged_seurat_filtered)
ElbowPlot(merged_seurat_filtered)
merged_seurat_filtered <- FindNeighbors(object = merged_seurat_filtered, dims = 1:20)

merged_seurat_filtered <- FindClusters(object = merged_seurat_filtered,
                                       resolution = c(0.1,0.3,0.5,0.7,1))

merged_seurat_filtered <- RunUMAP(merged_seurat_filtered, dims = 1:20)

DimPlot(merged_seurat_filtered,
        group.by = "RNA_snn_res.0.1",
        label = T)

Idents(merged_seurat_filtered) <- "RNA_snn_res.0.1"


DimPlot(merged_seurat_filtered, reduction = "umap", split.by = "Condition")


# perform integration to correct for batch effects ------
obj.list <- SplitObject(merged_seurat_filtered, split.by = 'Condition')

for(i in 1:length(obj.list)){
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
}


load("obj.list")

# select integration features
features <- SelectIntegrationFeatures(object.list = obj.list)

# find integration anchors (CCA)
anchors <- FindIntegrationAnchors(object.list = obj.list,
                                  anchor.features = features)

save(anchors, file = "anchors")

# integrate data

load("anchors")

seurat.integrated <- IntegrateData(anchorset = anchors)


seurat.integrated <- ScaleData(object = seurat.integrated)
seurat.integrated <- RunPCA(object = seurat.integrated)
ElbowPlot(seurat.integrated)
seurat.integrated <- RunUMAP(seurat.integrated, reduction = "pca", dims = 1:20)
seurat.integrated <- FindNeighbors(seurat.integrated, reduction = "pca", dims = 1:20)
seurat.integrated <- FindClusters(seurat.integrated, resolution = 0.1)
DimPlot(seurat.integrated, reduction = "umap", label = TRUE)
DimPlot(seurat.integrated, reduction = 'umap', split.by  = 'Condition', label = T,pt.size = 1,raster.dpi = c(800,800))
DimPlot(seurat.integrated, reduction = 'umap', label = T,pt.size = 1,raster.dpi = c(800,800))



#comparison of singleR and sctype for cell type annotation -----

#singleR -----

HPCA.data <- HumanPrimaryCellAtlasData(ensembl = FALSE)

# Leveraging cluster identity in your analysis
# now let's rerun our cluster identification using SingleR
seurat.integrated.sce <- as.SingleCellExperiment(seurat.integrated)
predictions <- SingleR(test=seurat.integrated.sce, assay.type.test=1,
                       ref=HPCA.data, labels=HPCA.data$label.main)

#now add back to singleCellExperiment object (or Seurat objects)
seurat.integrated.sce[["Cell_Type"]] <- predictions$labels
plotUMAP(seurat.integrated.sce, colour_by = "Cell_Type")

seurat_integrated2 <- as.Seurat(seurat.integrated.sce, counts = NULL)
DimPlot(seurat_integrated2, reduction = "UMAP",
        #split.by = "Condition", # this facets the plot
        group.by = "Cell_Type", # labels the cells with values from your group.by variable
        label = TRUE)


#sctype ----



# load libraries and functions

lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# get cell-type-specific gene sets from our in-built database (DB)
gs_list = gene_sets_prepare("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx", "Immune system") # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain

# assign cell types
scRNAseqData = readRDS(gzcon(url('https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/exampleData.RDS'))); #load example scRNA-seq matrix
es.max = sctype_score(scRNAseqData = scRNAseqData, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# View results, cell-type by cell matrix. See the complete example below
View(es.max)

# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")


# DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Heart" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)


# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = seurat.integrated[["integrated"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
# In case Seurat is used, it is either seurat.integrated[["RNA"]]@scale.data (default), seurat.integrated[["SCT"]]@scale.data, in case sctransform is used for normalization,
# or seurat.integrated[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(seurat.integrated@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(seurat.integrated@meta.data[seurat.integrated@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat.integrated@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])


seurat.integrated@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  seurat.integrated@meta.data$customclassif[seurat.integrated@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(seurat.integrated, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif', split.by = "Condition")        



# load libraries
lapply(c("ggraph","igraph","tidyverse", "data.tree"), library, character.only = T)

# prepare edges
cL_resutls=cL_resutls[order(cL_resutls$cluster),]; edges = cL_resutls; edges$type = paste0(edges$type,"_",edges$cluster); edges$cluster = paste0("cluster ", edges$cluster); edges = edges[,c("cluster", "type")]; colnames(edges) = c("from", "to"); rownames(edges) <- NULL

# prepare nodes
nodes_lvl1 = sctype_scores[,c("cluster", "ncells")]; nodes_lvl1$cluster = paste0("cluster ", nodes_lvl1$cluster); nodes_lvl1$Colour = "#f1f1ef"; nodes_lvl1$ord = 1; nodes_lvl1$realname = nodes_lvl1$cluster; nodes_lvl1 = as.data.frame(nodes_lvl1); nodes_lvl2 = c(); 
ccolss= c("#5f75ae","#92bbb8","#64a841","#e5486e","#de8e06","#eccf5a","#b5aa0f","#e4b680","#7ba39d","#b15928","#ffff99", "#6a3d9a","#cab2d6","#ff7f00","#fdbf6f","#e31a1c","#fb9a99","#33a02c","#b2df8a","#1f78b4","#a6cee3")
for (i in 1:length(unique(cL_resutls$cluster))){
  dt_tmp = cL_resutls[cL_resutls$cluster == unique(cL_resutls$cluster)[i], ]; nodes_lvl2 = rbind(nodes_lvl2, data.frame(cluster = paste0(dt_tmp$type,"_",dt_tmp$cluster), ncells = dt_tmp$scores, Colour = ccolss[i], ord = 2, realname = dt_tmp$type))
}
nodes = rbind(nodes_lvl1, nodes_lvl2); nodes$ncells[nodes$ncells<1] = 1;
files_db = openxlsx::read.xlsx(db_)[,c("cellName","shortName")]; files_db = unique(files_db); nodes = merge(nodes, files_db, all.x = T, all.y = F, by.x = "realname", by.y = "cellName", sort = F)
nodes$shortName[is.na(nodes$shortName)] = nodes$realname[is.na(nodes$shortName)]; nodes = nodes[,c("cluster", "ncells", "Colour", "ord", "shortName", "realname")]

mygraph <- graph_from_data_frame(edges, vertices=nodes)

# Make the graph
gggr<- ggraph(mygraph, layout = 'circlepack', weight=I(ncells)) + 
  geom_node_circle(aes(filter=ord==1,fill=I("#F5F5F5"), colour=I("#D3D3D3")), alpha=0.9) + geom_node_circle(aes(filter=ord==2,fill=I(Colour), colour=I("#D3D3D3")), alpha=0.9) +
  theme_void() + geom_node_text(aes(filter=ord==2, label=shortName, colour=I("#ffffff"), fill="white", repel = !1, parse = T, size = I(log(ncells,25)*1.5)))+ geom_node_label(aes(filter=ord==1,  label=shortName, colour=I("#000000"), size = I(3), fill="white", parse = T), repel = !0, segment.linetype="dotted")

library(scater)

scater::multiplot(DimPlot(seurat.integrated, reduction = "umap", label = TRUE, repel = TRUE, cols = ccolss), gggr, cols = 2)


# Differential Gene Expression ----

load("seurat.integrated") 


current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7)


new.cluster.ids <- c("Stromal cells", "Vascular endothelial cells", "Smooth muscle cells",
                     "Lymphoid cells", "Myeloid cells", "Schwann cells", "CLC_IL5RA positive cells", "7") 

names(new.cluster.ids) <- levels(seurat.integrated)
seurat.integrated <- RenameIdents(seurat.integrated, new.cluster.ids)


DimPlot(seurat.integrated, reduction = "umap", 
        split.by = "Condition",
        label = TRUE)


save(seurat.integrated, file = "seurat.integrated")

# finding diff exp genes in endothelial cells

seurat.integrated.Endothelial.cells <- subset(seurat.integrated, idents = "Vascular endothelial cells")


DimPlot(seurat.integrated.Endothelial.cells, reduction = "umap", label = TRUE)
Idents(seurat.integrated.Endothelial.cells)

# now we need to switch out 'Idents' to be treatment, rather than cluster
Idents(seurat.integrated.Endothelial.cells) <- seurat.integrated.Endothelial.cells$Condition
inf.vs.naive.markers <- FindMarkers(object = seurat.integrated.Endothelial.cells, 
                                    ident.1 = "Heart Failure", 
                                    ident.2 = "Control", 
                                    min.pct = 0)

inf.vs.naive.markers$pct.diff <- inf.vs.naive.markers$pct.1 - inf.vs.naive.markers$pct.2
inf.vs.naive.markers.df <- as_tibble(inf.vs.naive.markers, rownames = "geneID")
# Export DEGs for each cluster (ranked by avg_logFC > 0.5)
myTopHits <- inf.vs.naive.markers.df %>% arrange(desc(avg_log2FC))

FeaturePlot(seurat.integrated, 
            reduction = "umap", 
            features = "IGFBP5",
            pt.size = 0.4, 
            order = TRUE,
            split.by = "Condition",
            min.cutoff = 'q10',
            label = FALSE)



##############################################

seurat_integrated2 <- seurat.integrated
Idents(seurat_integrated2) <- seurat_integrated2$Condition
inf.vs.naive.markers <- FindMarkers(object = seurat_integrated2, 
                                    ident.1 = "Heart Failure", 
                                    ident.2 = "Control", 
                                    min.pct = 0)

inf.vs.naive.markers$pct.diff <- inf.vs.naive.markers$pct.1 - inf.vs.naive.markers$pct.2
inf.vs.naive.markers.df <- as_tibble(inf.vs.naive.markers, rownames = "geneID")
# Export DEGs for each cluster (ranked by avg_logFC > 0.5)
myTopHits <- inf.vs.naive.markers.df %>% arrange(desc(avg_log2FC))


FeaturePlot(seurat.integrated, 
            reduction = "umap", 
            features = "CKM",
            pt.size = 0.4, 
            order = TRUE,
            split.by = "Condition",
            min.cutoff = 'q10',
            label = FALSE)

my_fav_genes <- c("CXCL13", "IL17A")

FeaturePlot(seurat.integrated, 
            reduction = "umap", 
            features = my_fav_genes,
            pt.size = 0.4, 
            order = TRUE,
            split.by = "Condition",
            min.cutoff = 'q10',
            label = FALSE)



DotPlot(seurat.integrated, features = rev(as.character(unique(myTopHits$geneID[1:50]))), group.by = "customclassif",
        assay = "RNA", cols = c("Blue", "Red")) + coord_flip()


markers_genes <- FindAllMarkers(seurat.integrated, log2FC.threshold = 0.2, test.use = "wilcox",
                                min.pct = 0.1, min.diff.pct = 0.2, only.pos = TRUE, max.cells.per.ident = 50,
                                assay = "RNA")

seurat.integrated@meta.data$customclassif

mypar(2, 5, mar = c(4, 6, 3, 1))
for (i in unique(seurat.integrated@meta.data$customclassif)) {
  barplot(sort(setNames(myTopHits$avg_log2FC, myTopHits$geneID)[seurat.integrated@meta.data$customclassif == i], F),
          horiz = T, las = 1, main = paste0(i, " vs. rest"), border = "white", yaxs = "i")
  abline(v = c(0, 0.25), lty = c(1, 2))
}


save(seurat.integrated, file = "seurat.integrated2")


# Find all markers ------

load("seurat.integrated2")
load("markers_genes")

markers_genes <- FindAllMarkers(seurat.integrated, log2FC.threshold = 0.2, test.use = "wilcox",
                                min.pct = 0.1, min.diff.pct = 0.2, only.pos = TRUE, max.cells.per.ident = 50,
                                assay = "RNA")

markers_genes %>%
  group_by(cluster) %>%
  top_n(-25, p_val_adj) -> top25
top25

mypar(2, 4, mar = c(4, 6, 3, 1))
for (i in unique(top25$cluster)) {
  barplot(sort(setNames(top25$avg_log2FC, top25$gene)[top25$cluster == i], F),
          horiz = T, las = 1, main = paste0(i, " vs. rest"), border = "white", yaxs = "i")
  abline(v = c(0, 0.25), lty = c(1, 2))
}


markers_genes %>%
  group_by(cluster) %>%
  top_n(-5, p_val_adj) -> top5



# create a scale.data slot for the selected genes
seurat.integrated <- ScaleData(seurat.integrated, features = as.character(unique(top5$gene)), assay = "RNA")
DoHeatmap(seurat.integrated, features = as.character(unique(top5$gene)),
          assay = "RNA", angle = 90, size = 4)

DotPlot(seurat.integrated, features = rev(as.character(unique(top5$gene))),
        assay = "RNA", cols = c("Blue", "Red"), split.by = "Condition") + RotatedAxis()


top5 %>%
  group_by(cluster) %>%
  top_n(-1, p_val) -> top1


# set pt.size to zero if you do not want all the points to hide the violin
# shapes, or to a small value like 0.1


#all things not comming do one by one or loop

VlnPlot(object = seurat.integrated, features = c("FBLN1","CD69"), pt.size = 0)

gene_list <- as.character(unique(top1$gene))


mypar(7, 3, mar = c(4, 6, 3, 1))
plotlist <- list()
for (i in gene_list) {
  plotlist[[i]] <- VlnPlot(object = seurat.integrated, features = i , pt.size = 0)
}

grid.arrange(grobs = plotlist, ncol = 4)


# HF vs CTL DGE ----

up <- VlnPlot(seurat.integrated, features = as.character(unique(myTopHits$geneID[1:5])),
        ncol = 5, group.by = "Condition", assay = "RNA", pt.size = 0)

down <- VlnPlot(seurat.integrated, features = as.character(tail(myTopHits$geneID,5)),
        ncol = 5, group.by = "Condition", assay = "RNA", pt.size = 0)

plot_grid(up, down, ncol = 1)

VlnPlot(seurat.integrated, features = as.character(unique(myTopHits$geneID[1:5])), ncol = 5,
        split.by = "Condition", assay = "RNA", pt.size = 0)

topDEG = c("MYL2","TCAP","MB","TNNC1","ACTA1","TNNT2","TNNI3","ANKRD1",
           "DES","TPM1","TAGLN", "MT1X",  "MT1E",  "APOE",  "NRXN1", "ACTA2", "PDK4",  "MT1G",  "MYH11", "MT1A")

DotPlot(seurat.integrated, features = topDEG ,
        group.by = "Condition", assay = "RNA", cols = c("Blue", "Red")) + coord_flip()

FeaturePlot(object = seurat.integrated, features = c("MYL2", "APOE"), split.by = "Condition")


seurat.integrated <- ScaleData(seurat.integrated, features = as.character(unique(myTopHits$geneID)), assay = "RNA")

DoHeatmap(seurat.integrated, features = as.character(unique(myTopHits$geneID)),
          assay = "RNA", angle = 90, size = 4, group.by = "Condition") + scale_fill_gradientn(colors = c("blue", "white", "red"))



DoHeatmap(seurat.integrated, features = as.character(unique(myTopHits$geneID)),
          assay = "RNA", angle = 90, size = 4, group.by = "Condition", disp.min = -3, disp.max = 3,draw.lines = F) 


RidgePlot(seurat.integrated, features = c("C1R", "SLC9A3R2"))



DotPlot(seurat.integrated, features = top5$gene, dot.scale = 8, split.by = "Condition") +
  RotatedAxis()


library(SeuratDisk)


SaveH5Seurat(seurat.integrated, filename = "seurat.integrated.h5Seurat")
Convert("seurat.integrated.h5Seurat", dest = "h5ad")
write.csv(myTopHits, file = "myTopHits.csv")


seurat.subset <- subset(seurat.integrated, features = myTopHits$geneID)


nrow(seurat.subset@assays$RNA@data)

normalizedCounts <- as.matrix(seurat.subset@assays$RNA@data)

normalizedCounts2 <- as.data.frame(normalizedCounts)

head(normalizedCounts2)


cell_identities <- Idents(object = seurat.subset)


###################################################

library(gplots)
library(ComplexHeatmap)
library(pheatmap)

myheatcolors1 <- bluered(100)

clustRows <- hclust(as.dist(1-cor(t(diffGenes), method="pearson")), method="complete") 
clustColumns <- hclust(as.dist(1-cor(diffGenes, method="spearman")), method="complete")

plot(clustColumns)

colanno <- HeatmapAnnotation(group = targets$group, sex = targets$sex,
                             show_annotation_name = T, show_legend = T)
columnsplit <- as.vector(targets$group)
columnsplit <- gsub("_", " ", columnsplit)

scaled_mat = t(scale(t(normalizedCounts)))

Heatmap(scaled_mat, use_raster = T, show_column_names = F, show_row_dend = F, cluster_columns = T)

cluster_rows = clustRows, cluster_columns = F, row_names_side = "left", col = myheatcolors1,
              show_column_dend = F, column_split = columnsplit, row_split = 2, row_gap = unit(c(2, 1), "mm"),
              column_gap = unit(c(2, 1), "mm"), border = T,
              row_title = NULL, name = "Z-score", row_names_gp = gpar(fontsize = 2),
              show_row_names = F, show_row_dend = F, show_column_names = T)


log2_val <- as.data.frame(myTopHits[,1], row.names = rownames(myTopHits)) 

log2_val <-  as.data.frame(log2_val[rownames(diffGenes),], row.names = rownames(diffGenes)) 
colnames(log2_val)[1] <- "logFC"


colours <- bluered(100)

h2 <- Heatmap(as.matrix(log2_val) , 
              cluster_rows = clustRows, name="logFC", col = colours, show_row_names = F)

h = h1+h2
h

########### THE END ###########################