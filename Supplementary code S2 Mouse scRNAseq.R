#Top ----

#Analysis of the single cell RNA-seq data from two normal mouse lungs.
#Data generation: Nikita Joshi, Alexandra McQuattie-Pimentel, James Walter.
#Analysis: Paul Reyfman, Ziyou Ren, Kishore Anekalla, Scott Budinger, Alexander Misharin.
#Reference: Reyfman et al., AJRCCM, 2018.
#Analysis performed using Seurat 2.3.0.
#Code is based tutorials and vignettes from satijalab.org/seurat.  

#load libraries
library(Seurat)
library(dplyr)
library(Matrix)
library(cowplot)
library(ggplot2)

#Initial data processing, PCA and clustering ----
#SC01 - Naive mouse #1, SC02 - Naive mouse #2. Both mice are C57Bl/6 males, 3 months old. 
#Load data and create aggregated Seurat object
SC01.data <- Read10X(data.dir = "/SC01/mm10/")
SC01 <- CreateSeuratObject(raw.data = SC01.data, min.cells = 3, min.genes = 200, project = "SC01")
SC01@cell.names <- paste("SC01", SC01@cell.names, sep = "_")
colnames(x = SC01@raw.data) <- paste("SC01", colnames(x = SC01@raw.data), sep = "_")
rownames(x = SC01@meta.data) <- paste("SC01", rownames(x = SC01@meta.data), sep = "_")
head(SC01@cell.names)
tail(SC01@cell.names)
SC02.data <- Read10X(data.dir = "/SC02/mm10/")
lung <- AddSamples(object = SC01, new.data = SC02.data, add.cell.id = "SC02")
head(lung@cell.names)
tail(lung@cell.names)

#Data inspection and basic filtering
mito.genes <- grep(pattern = "^mt-", x = rownames(x = lung@data), value = TRUE)
percent.mito <- Matrix::colSums(lung@raw.data[mito.genes, ])/Matrix::colSums(lung@raw.data)
lung <- AddMetaData(object = lung, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = lung, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
GenePlot(object = lung, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = lung, gene1 = "nUMI", gene2 = "nGene")
#Let's filter out cells with percent.mito > 0.1 and nGene <300 and >4000
#This step is highly subjective, skipping filtering will increase number of ciliated and club cells downstream. 
lung <- FilterCells(object = lung, subset.names = c("nGene", "percent.mito"), 
                    low.thresholds = c(300, -Inf), high.thresholds = c(4000, 0.1))
#Normalize, find variable genes and scale
lung <- NormalizeData(object = lung, normalization.method = "LogNormalize", scale.factor = 10000)
lung <- FindVariableGenes(object = lung, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = lung@var.genes) #9217
lung <- ScaleData(object = lung, vars.to.regress = c("nUMI", "percent.mito"), do.par = T, num.cores = 2)
#Perform PCA
lung <- RunPCA(object = lung, pc.genes = lung@var.genes, do.print = TRUE, pcs.print = 1:5, 
               genes.print = 5, pcs.compute = 100)
PCAPlot(object = lung, dim.1 = 1, dim.2 = 2)
PCAPlot(object = lung, dim.1 = 1, dim.2 = 3)
PCHeatmap(object = lung, pc.use = 1:9, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = lung, pc.use = 10:18, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = lung, pc.use = 19:27, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = lung, pc.use = 28:36, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCElbowPlot(object = lung, num.pc = 50)
#Perform clustering using first 35 principal components. 
lung <- FindClusters(object = lung, reduction.type = "pca", dims.use = 1:35, 
                     resolution = 0.7, print.output = 0, save.SNN = TRUE, force.recalc = T)
PrintFindClustersParams(object = lung)
lung <- RunTSNE(object = lung, dims.use = 1:35, do.fast = TRUE, check_duplicates = FALSE)
TSNEPlot(object = lung, do.label = T)

lung.markers <- FindAllMarkers(object = lung, only.pos = TRUE, min.pct = 0.5, thresh.use = 0.5, max.cells.per.ident = 200)
write.csv(lung.markers, "initial_lung_markers_PC1_35_res0.7.csv")

FeaturePlot(lung, features.plot = c("Sftpc", "Ager","Foxj1", "Scgb3a1"))

#Initial list of clusters ----
#Cluster0 - B cells
#Cluster1 - Alveolar type 2 cells
#Cluster2 - T cells
#Cluster3 - Alveolar macrophages
#Cluster4 - Endothelial (alveolar capillaries)
#Cluster5 - NK cells
#Cluster6 - Neutrophils
#Cluster7 - T cells, NKT cells and ILC
#Cluster8 - Non-classical monocytes
#Cluster9 - Endothelial cells
#Cluster10 - Classical monocytes
#Cluster11 - Interstitial macrophages
#Cluster12 - Fibroblasts
#Cluster13 - Ciliated and club cells
#Cluster14 - DC2
#Cluster15 - Proliferating/cell cycle cells
#Cluster16 - Alveolar type 1 cells
#Cluster17 - DC1
#Cluster18 - Endothelial cells
#Cluster19 - B cell doublets
#Cluster20 - Smooth muscle cells
#Cluster21 - pDC
#Cluster22 - Alveolar macrophages doublets
#Cluster23 - Basophils
#Cluster24 - Ccr7 DC
#Cluster25 - Alveolar type 2 cells doublets

#Resolving megakaryocytes 
FeaturePlot(lung, features.plot = c("Ppbp","Pf4", "Gp9", "Gp1bb"))
TSNEPlot(object = lung, do.label = T, do.hover = T, data.hover = c("ident", "Ppbp"))
select.cells <- TSNEPlot(object = lung, do.identify = TRUE) #30 cells
lung <- SetIdent(object = lung, cells.use = select.cells, ident.use = "26")
#Cluster26 - Megakaryocytes
TSNEPlot(object = lung, do.label = T)

#Separating ciliated and club cells
FeaturePlot(lung, features.plot = c("Foxj1","Scgb1a1"))
TSNEPlot(object = lung, do.label = T, do.hover = T, data.hover = c("ident", "Foxj1"))
select.cells <- TSNEPlot(object = lung, do.identify = TRUE) #153 cells
lung <- SetIdent(object = lung, cells.use = select.cells, ident.use = "27")
#Cluster27 - Ciliated airway cells
TSNEPlot(object = lung, do.label = T)

#Resolve two subsets of interstitial macropahges
C11 <- SubsetData(lung, ident.use = c(11), do.clean = T) #subset on C11
mito.genes <- grep(pattern = "^mt-", x = rownames(x = C11@data), value = TRUE)
percent.mito <- Matrix::colSums(C11@raw.data[mito.genes, ])/Matrix::colSums(C11@raw.data)
C11 <- AddMetaData(object = C11, metadata = percent.mito, col.name = "percent.mito")
C11 <- NormalizeData(object = C11, normalization.method = "LogNormalize", scale.factor = 10000) #Log normalize
C11 <- FindVariableGenes(object = C11, mean.function = ExpMean, dispersion.function = LogVMR, 
                         x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.plot = F)
length(x = C11@var.genes)
C11 <- ScaleData(object = C11, vars.to.regress = c("nUMI"), do.par = T, num.cores = 2)
C11 <- RunPCA(object = C11, pc.genes = C11@var.genes, do.print = TRUE, pcs.print = 1:5, 
              genes.print = 5, pcs.compute = 50)
PCAPlot(object = C11, dim.1 = 1, dim.2 = 2)
PCAPlot(object = C11, dim.1 = 1, dim.2 = 3)
PCAPlot(object = C11, dim.1 = 2, dim.2 = 3, group.by = "orig.ident")
PCHeatmap(object = C11, pc.use = 1:9, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = C11, pc.use = 10:18, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = C11, pc.use = 19:27, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = C11, pc.use = 28:36, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCElbowPlot(object = C11, num.pc = 50)
C11 <- FindClusters(object = C11, reduction.type = "pca", dims.use = 1:3, 
                    resolution = 0.3, print.output = 0, save.SNN = TRUE, force.recalc = T)
PrintFindClustersParams(object = C11)
C11 <- RunTSNE(object = C11, dims.use = 1:3, do.fast = TRUE, check_duplicates = FALSE)
TSNEPlot(object = C11, do.label = T)
FeaturePlot(C11, features.plot = c("Lgals3","Ccr2", "Mrc1", "Lyve1"), cols.use = c("grey","blue"))
C11.markers <- FindAllMarkers(object = C11, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.5)
#Cluster 0 is characterized by Lgals3 and Ccr2, matching peribronchial interstitial macrophages. 
#Cluster 1 is characterized by Mrc1 and Lyve1, matching perivascular interstitial macrophages.
#Rename cells in the main object:
C28<-WhichCells(C11, ident = 0)
lung <- SetIdent(object = lung, cells.use = C28, ident.use = "28")
TSNEPlot(object = lung, do.label = T)

#For presenting data in the paper - remove doublets and damaged/proliferating cells. 
lung.clean <- SubsetData(lung, ident.remove = c(15,19,22,25))
TSNEPlot(object = lung.clean, do.label = T)

current.cluster.ids <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,17,18,20,21,23,24,26,27,28)
new.cluster.ids <- c("B cells", "AT2 cells", "Naive T cells", "Alveolar macrophages", 
                     "Endothelial (alveolar capillaries)", "NK cells", "Neutrophils",
                     "T cells, NKT cells and ILC", "Non-classical monocytes", "Endothelial cells",
                     "Classical monocytes", "Peribronchial interstitial macrophages", "Fibroblasts", 
                     "Club cells", "DC2", "AT1 cells", 
                     "DC1", "Endothelial cells", "Smooth muscle cells",
                     "pDC", "Basophils", "Ccr7+ DC", "Megakaryocytes", "Ciliated cells", "Perivascular interstitial macrophages")
lung.clean@ident <- plyr::mapvalues(x = lung.clean@ident, from = current.cluster.ids, to = new.cluster.ids)
TSNEPlot(object = lung.clean, do.label = T, no.legend = T, pt.size = 2) #Figure S2B

#To vizualize top 3 genes per cell type, #Figure S2C 
lung.markers <- FindAllMarkers(object = lung.clean, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.5, max.cells.per.ident = 200)
write.csv(lung.markers, "Supplementary_Table_2_Murine_SC_Marker_Genes.csv") #Supplementary_Table_2_Murine_SC_Marker_Genes
topgenes <- lung.markers %>% group_by(cluster) %>% top_n(3, avg_logFC)
cluster.averages <- AverageExpression(object = lung.clean, return.seurat = TRUE, show.progress = T)
DoHeatmap(object = cluster.averages, genes.use = topgenes$gene, group.label.rot = TRUE, group.cex = 0) #Figure S2C
save(lung.markers, file = "mouse_lung_markers.Robj")

#Generate barplot reflecting proportion of the individual libraries in each cell type , #Figure S2D 
barplotdata = lung.clean@meta.data
barplotdata$ident = as.character(lung.clean@ident)
toplotdata <- barplotdata %>% group_by(ident,orig.ident) %>% count()
ggplot(data=toplotdata,mapping = aes(x=ident,y=n,fill=orig.ident)) + geom_col(inherit.aes = T) 
sums <- toplotdata %>% group_by(ident) %>% summarise(Sum=sum(n))
sums <- rep(sums$Sum,each=2)
toplotdata$sums <- sums
toplotdata$proportion <- (toplotdata$n/toplotdata$sums)
ggplot(data=toplotdata,mapping = aes(x=ident,y=proportion,fill=orig.ident)) + geom_col(inherit.aes = T) +
  theme(axis.text = element_text(angle=90,hjust = 1)) #Figure S2D

save(lung.clean, file = "Figure_S2_mouse_lung.Robj")
load("Figure_S2_mouse_lung.Robj")

#Figure S4C, mouse Wnts in epithelial and mesenchymal cells ----
FigS4 <-SubsetData(lung, ident.use = c(1,12,13,16,20,27))
TSNEPlot(FigS4, do.label = T)
current.FigS4.ids <- c(1,12,13,16,20,27)
new.FigS4.ids <- c("AT2", "Fib","Club", "AT1", "SM", "Cil")
FigS4@ident <- plyr::mapvalues(x = FigS4@ident, from = current.FigS4.ids, to = new.FigS4.ids)
TSNEPlot(object = FigS4, do.label = T, pt.size = 0.5)
VlnPlot(FigS4, features.plot = c("Wnt2", "Wnt2b", "Wnt3a", "Wnt4", "Wnt5a", "Wnt5b", "Wnt6", "Wnt7a", "Wnt7b", "Wnt9a", "Wnt10b", "Wnt11", "Lgr4", "Lgr5", "Lgr6", "Axin2", "Porcn", "Wls", "Pdgfra"), same.y.lims = T, size.x.use = 0 ,size.title.use = 12, nCol = 4) #Figure S4C

#Figure S9 Mouse Wnts in AT2 ----
#Subset on AT2 cells
AT2 <- SubsetData(lung, ident.use = c(1), do.clean = T)
AT2 <- NormalizeData(object = AT2, normalization.method = "LogNormalize", scale.factor = 10000)
AT2 <- FindVariableGenes(object = AT2, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = AT2@var.genes)
AT2 <- ScaleData(object = AT2, vars.to.regress = c("nUMI", "percent.mito"))
AT2 <- RunPCA(object = AT2, pc.genes = AT2@var.genes, do.print = TRUE, pcs.print = 1:5, 
               genes.print = 5, pcs.compute = 20)
PCAPlot(object = AT2, dim.1 = 1, dim.2 = 2)
PCAPlot(object = AT2, dim.1 = 1, dim.2 = 3)
PCHeatmap(object = AT2, pc.use = 1:9, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = AT2, pc.use = 10:18, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = AT2, pc.use = 19:27, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = AT2, pc.use = 28:36, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCElbowPlot(object = AT2, num.pc = 20)
AT2 <- FindClusters(object = AT2, reduction.type = "pca", dims.use = 1:12, 
                     resolution = 0.1, print.output = 0, save.SNN = TRUE, force.recalc = T)
PrintFindClustersParams(object = AT2)
AT2 <- RunTSNE(object = AT2, dims.use = 1:12, do.fast = TRUE, check_duplicates = FALSE)
TSNEPlot(object = AT2, do.label = T)
TSNEPlot(object = AT2, do.label = T, group.by = "orig.ident")

FeaturePlot(AT2, features.plot = c("Sftpc", "Ager","Nkg7", "Cd79a"))
#Clusters 2 and 3 contain B and NK cell doublets, let's remove them. 
AT2 <- SubsetData(AT2, ident.use = c(0,1), do.clean = T)
#TSNEPlot(AT2, do.label = T)
AT2 <- NormalizeData(object = AT2, normalization.method = "LogNormalize", scale.factor = 10000)
AT2 <- FindVariableGenes(object = AT2, mean.function = ExpMean, dispersion.function = LogVMR, 
                         x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = AT2@var.genes)
AT2 <- ScaleData(object = AT2, vars.to.regress = c("nUMI", "percent.mito", "orig.ident"))
AT2 <- RunPCA(object = AT2, pc.genes = AT2@var.genes, do.print = TRUE, pcs.print = 1:5, 
              genes.print = 5, pcs.compute = 20)
PCAPlot(object = AT2, dim.1 = 1, dim.2 = 2)
PCAPlot(object = AT2, dim.1 = 1, dim.2 = 3)
PCHeatmap(object = AT2, pc.use = 1:9, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = AT2, pc.use = 10:18, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = AT2, pc.use = 19:27, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = AT2, pc.use = 28:36, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCElbowPlot(object = AT2, num.pc = 20)
AT2 <- FindClusters(object = AT2, reduction.type = "pca", dims.use = 1:6, 
                    resolution = 0.3, print.output = 0, save.SNN = TRUE, force.recalc = T)
PrintFindClustersParams(object = AT2)
AT2 <- RunTSNE(object = AT2, dims.use = 1:6, do.fast = TRUE, check_duplicates = FALSE)
TSNEPlot(object = AT2, do.label = T) #Figure S10A
TSNEPlot(object = AT2, do.label = F, group.by = "orig.ident") #Figure S10B

#To illustrate co-expression of Wnt ligands, Figure S9A and B
a=FeaturePlot(object = AT2, features.plot = c("Wnt3a", "Wnt7b"), cols.use = c("grey", "red", "blue", "green"), overlay = TRUE, no.legend = FALSE, do.return=T, pt.size = 0.5)
b=FeaturePlot(object = AT2, features.plot = c("Wnt3a", "Wnt4"), cols.use = c("grey", "red", "blue", "green"), overlay = TRUE, no.legend = FALSE, do.return=T, pt.size = 0.5)
c=FeaturePlot(object = AT2, features.plot = c("Wnt3a", "Wnt9a"), cols.use = c("grey", "red", "blue", "green"), overlay = TRUE, no.legend = FALSE, do.return=T, pt.size = 0.5)
d=FeaturePlot(object = AT2, features.plot = c("Wnt7b", "Wnt4"), cols.use = c("grey", "red", "blue", "green"), overlay = TRUE, no.legend = FALSE, do.return=T, pt.size = 0.5)
e=FeaturePlot(object = AT2, features.plot = c("Wnt7b", "Wnt9a"), cols.use = c("grey", "red", "blue", "green"), overlay = TRUE, no.legend = FALSE, do.return=T, pt.size = 0.5)
f=FeaturePlot(object = AT2, features.plot = c("Wnt4", "Wnt9a"), cols.use = c("grey", "red", "blue", "green"), overlay = TRUE, no.legend = FALSE, do.return=T, pt.size = 0.5)
g=FeaturePlot(object = AT2, features.plot = c("Wnt3a", "Axin2"), cols.use = c("grey", "red", "blue", "green"), overlay = TRUE, no.legend = FALSE, do.return=T, pt.size = 0.5)
h=FeaturePlot(object = AT2, features.plot = c("Wnt7a", "Axin2"), cols.use = c("grey", "red", "blue", "green"), overlay = TRUE, no.legend = FALSE, do.return=T, pt.size = 0.5)
i=FeaturePlot(object = AT2, features.plot = c("Wnt7b", "Axin2"), cols.use = c("grey", "red", "blue", "green"), overlay = TRUE, no.legend = FALSE, do.return=T, pt.size = 0.5)
cowplot::plot_grid(plotlist = c(a,b,c,d,e,f,g,h,i))

save(AT2, file = "AT2.Robj")

#Figure S10: Mouse AEPs. 
AT2.markers <- FindAllMarkers(object = AT2, only.pos = TRUE, min.pct = 0.50, thresh.use = 0.25, test.use = "roc")
write.csv(AT2.markers, "AT2.markers_PC1_6_res0.3roc.csv")
AT2.markers <- FindAllMarkers(object = AT2, only.pos = TRUE, min.pct = 0.50, thresh.use = 0.25)
write.csv(AT2.markers, "AT2.markers_PC1_6_res0.3.csv")
top10 <- AT2.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = AT2, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)

VlnPlot(AT2, features.plot = c("Scd1", "Fasn", "Sftpc", "Lyz2", "Ctsh", "Krt8", "Axin2", "Tm4sf1"), size.title.use = 14, size.x.use = 0, nCol = 4, point.size.use = 0.01) #Figure S10C
FeaturePlot(AT2, features.plot = "Tm4sf1")

#Let's test if any of these subclusters within AT2 cluster is enriched for AEP genes (Zacharias et al., 2018, Nature)
AEPgenes <- list(readLines(con = "AEPtop500.txt"))
AT2 <- AddModuleScore(AT2, genes.list = c(AEPgenes), enrich.name = "AEP", ctrl.size = 5)
head(x = AT2@meta.data)
FeaturePlot(AT2, features.plot = "AEP1")
VlnPlot(AT2, features.plot = "AEP1", size.x.use = 0, size.title.use = 0) #Figure S10D
RidgePlot(AT2, features.plot = "AEP1")
GenePlot(AT2, gene1 = "Tm4sf1", gene2 = "Axin2")
FeaturePlot(object = AT2, features.plot = c("Tm4sf1", "Axin2"), cols.use = c("grey", "red", "blue", "green"), overlay = TRUE, no.legend = FALSE)

rm(list=ls()) #clean environment
#End


