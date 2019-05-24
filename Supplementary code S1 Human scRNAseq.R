#Analysis of the single cell RNA-seq data from eight lung transplant donors, eight lung transplant recipients with pulmonary fibrosis,
#and one cryobiopsy specimen from a patient with pulmonary fibrosis early in the disease course.
#Data generation: Nikita Joshi, Alexandra McQuattie-Pimentel, James Walter.
#Analysis: Paul Reyfman, Ziyou Ren, Kishore Anekalla, Scott Budinger, Alexander Misharin.
#Reference: Reyfman et al., AJRCCM, 2018.
#Analysis performed using Seurat 2.3.0.
#Code is based tutorials and vignettes from satijalab.org/seurat.  

# Set R Environment ---------------------------------------------------------
# Load necessary packages
library(Seurat)
library(dplyr)
library(Matrix)
library('Cairo')
library(ggplot2)
library(tibble)
library(ggrepel)
library(cowplot)

#First we load and merge all of the donor samples

# donor01 analysis ---------------------------------------------------------

# Load the donor01 dataset
donor01.data <- Read10X(data.dir = "/donor01_tables/")

donor01 <- CreateSeuratObject(raw.data = donor01.data, min.cells = 3, min.genes = 200, 
                           project = "donor01")

# Identify mitochondrial genes
mito.genes <- grep(pattern = "^MT-", x = rownames(x = donor01@data), value = TRUE)

# Create meta data column containing percentage of mitochondrial genes for each cell
percent.mito <- Matrix::colSums(donor01@raw.data[mito.genes, ])/Matrix::colSums(donor01@raw.data)

# Add percent.mito to object meta data
donor01 <- AddMetaData(object = donor01, metadata = percent.mito, col.name = "percent.mito")

#Data inspection and basic filtering
VlnPlot(object = donor01, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
GenePlot(object = donor01, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = donor01, gene1 = "nUMI", gene2 = "nGene")

donor01 <- FilterCells(object = donor01, subset.names = c("nGene", "percent.mito"), 
                    low.thresholds = c(200, -Inf), high.thresholds = c(5000, 0.2))

#Normalize, find variable genes and scale
donor01 <- NormalizeData(object = donor01, normalization.method = "LogNormalize", 
                      scale.factor = 10000)
donor01 <- FindVariableGenes(object = donor01, mean.function = ExpMean, dispersion.function = LogVMR, 
                             x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = donor01@var.genes)
donor01 <- ScaleData(object = donor01, vars.to.regress = c("nUMI", "percent.mito"))


#Perform PCA
donor01 <- RunPCA(object = donor01, pc.genes = donor01@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5, pcs.compute = 40)
PrintPCA(object = donor01, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(object = donor01, pcs.use = 1:2)
PCAPlot(object = donor01, dim.1 = 1, dim.2 = 2)
donor01 <- ProjectPCA(object = donor01, do.print = FALSE)
PCHeatmap(object = donor01, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = donor01, pc.use = 13:20, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCElbowPlot(object = donor01, num.pc = 40)

#Perform clustering using first 13 principal components. 
donor01 <- FindClusters(object = donor01, reduction.type = "pca", dims.use = 1:13, 
                     resolution = 0.3, save.SNN = TRUE, print.output = FALSE, force.recalc = T)
PrintFindClustersParams(object = donor01)

donor01 <- RunTSNE(object = donor01, dims.use = 1:13, do.fast = TRUE)
TSNEPlot(object = donor01, do.label = T)

#Construct the cluster tree and merge nodes based on visual inspection
donor01 <- BuildClusterTree(object = donor01, do.reorder = TRUE, reorder.numeric = TRUE, do.plot = T)
donor01 <- MergeNode(donor01, node.use = 16, rebuild.tree = T)
TSNEPlot(object = donor01, do.label = T)

#Subset Cluster 4 (Monocytes and DCs)
donor01.cluster04 <- SubsetData(donor01, ident.use = c(4))
donor01.cluster04 <- ScaleData(object = donor01.cluster04, vars.to.regress = c("nUMI", "percent.mito"))
donor01.cluster04 <- FindVariableGenes(object = donor01.cluster04, mean.function = ExpMean, dispersion.function = LogVMR, 
                                    x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = donor01.cluster04@var.genes)

donor01.cluster04 <- RunPCA(object = donor01.cluster04, pc.genes = donor01.cluster04@var.genes, do.print = TRUE, pcs.print = 1:5, 
                         genes.print = 5)
donor01.cluster04 <- ProjectPCA(object = donor01.cluster04, do.print = FALSE)
PCHeatmap(object = donor01.cluster04, pc.use = 1:12, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)

PCElbowPlot(object = donor01.cluster04, num.pc = 12)

#Clustering of Subsetted Cells
donor01.cluster04 <- RunTSNE(object = donor01.cluster04, dims.use = 1:11, do.fast = TRUE)
donor01.cluster04 <- FindClusters(object = donor01.cluster04, reduction.type = "pca", dims.use = 1:11, 
                               resolution = 0.5, save.SNN = TRUE)
TSNEPlot(donor01.cluster04)

#Markers of Subsetted Cell Subclusters
donor01.cluster04.markers <- FindAllMarkers(object = donor01.cluster04, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
donor01.cluster04.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
FeaturePlot(object = donor01.cluster04, features.plot = c("HLA-DQA1", "CD1C", "CD14", "S100A8"), min.cutoff = "q9",
            cols.use = c("lightgrey", "blue"), pt.size = 1)

#Cell assignment for Subclusters
donor01.dcs <- WhichCells(object = donor01.cluster04, ident = c(0))
donor01.monos <- WhichCells(object = donor01.cluster04, ident = c(1))

#Rename clusters
donor01.ident <- c(1, 2, 3, 4, 5, 6, 7, 8)
donor01.new.ident <- c("Unassigned", "Mast Cells", "Endothelial Cells", "Monocytes and DCs", 
                     "Macrophages", "Unassigned", "AT1 Cells", "AT2 Cells")

donor01@ident <- plyr::mapvalues(x = donor01@ident, from = donor01.ident, to = donor01.new.ident)

donor01 <- SetIdent(object = donor01, cells.use = donor01.dcs, ident.use = "Dendritic Cells")
donor01 <- SetIdent(object = donor01, cells.use = donor01.monos, ident.use = "Monocytes")

donor01@ident <- ordered(donor01@ident,levels = c("Macrophages", "AT2 Cells", "Monocytes", "Dendritic Cells",
                                            "AT1 Cells", "Endothelial Cells", "Mast Cells", "Unassigned"))

#t-SNE Plot (Figure E5)
TSNEPlot(object = donor01, do.label = T)

#Find marker genes for clusters (Table E5)
donor01.markers <- FindAllMarkers(object = donor01, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.table(donor01.markers %>% group_by(cluster), file=paste(datadir, "/donor01.markers.txt",sep=""), row.names=FALSE, sep="\t")

#Add Cluster Identity Metadata to Separate Metadata Column
donor01 <- AddMetaData(object = donor01, metadata = donor01@ident, col.name = "indiv.cell.ident")
donor01@meta.data[, "diagnosis"] <- "donor"
donor01@meta.data[, "condition"] <- "donor"

donor01@meta.data$res.0.3 <- NULL
donor01@meta.data$tree.ident <- NULL

#Macrophage Heterogeneity Feature Plots (Figure E7E)
FeaturePlot(object = donor01, features.plot = c("MRC1", "FABP4", "SPP1", "CHI3L1"), min.cutoff = "q9",
              cols.use = c("lightgrey", "blue"), nCol = 1, pt.size = 0.5)

#Save object
save(donor01, file = "/donor01.Robj")
load(file = "/donor01.Robj")

# donor02 analysis ---------------------------------------------------------

# Load the donor02 dataset
donor02.data <- Read10X(data.dir = "/donor02_tables/")

donor02 <- CreateSeuratObject(raw.data = donor02.data, min.cells = 3, min.genes = 200, 
                              project = "donor02")

# Identify mitochondrial genes
mito.genes <- grep(pattern = "^MT-", x = rownames(x = donor02@data), value = TRUE)

# Create meta data column containing percentage of mitochondrial genes for each cell
percent.mito <- Matrix::colSums(donor02@raw.data[mito.genes, ])/Matrix::colSums(donor02@raw.data)

# Add percent.mito to object meta data
donor02 <- AddMetaData(object = donor02, metadata = percent.mito, col.name = "percent.mito")

#Data inspection and basic filtering
VlnPlot(object = donor02, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
GenePlot(object = donor02, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = donor02, gene1 = "nUMI", gene2 = "nGene")

donor02 <- FilterCells(object = donor02, subset.names = c("nGene", "percent.mito"), 
                       low.thresholds = c(200, -Inf), high.thresholds = c(6000, 0.15))

#Normalize, find variable genes and scale
donor02 <- NormalizeData(object = donor02, normalization.method = "LogNormalize", 
                         scale.factor = 10000)
donor02 <- FindVariableGenes(object = donor02, mean.function = ExpMean, dispersion.function = LogVMR, 
                             x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = donor02@var.genes)
donor02 <- ScaleData(object = donor02, vars.to.regress = c("nUMI", "percent.mito"))


#Perform PCA
donor02 <- RunPCA(object = donor02, pc.genes = donor02@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5, pcs.compute = 40)
PrintPCA(object = donor02, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(object = donor02, pcs.use = 1:2)
PCAPlot(object = donor02, dim.1 = 1, dim.2 = 2)
donor02 <- ProjectPCA(object = donor02, do.print = FALSE)
PCHeatmap(object = donor02, pc.use = 1:9, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = donor02, pc.use = 10:18, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = donor02, pc.use = 19:27, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCElbowPlot(object = donor02, num.pc = 40)

#Perform clustering using first 16 principal components. 
donor02 <- FindClusters(object = donor02, reduction.type = "pca", dims.use = 1:16, 
                     resolution = 0.3, print.output = 0, save.SNN = TRUE, force.recalc = T)
PrintFindClustersParams(object = donor02)
donor02 <- RunTSNE(object = donor02, dims.use = 1:16, do.fast = TRUE)
TSNEPlot(object = donor02, do.label = T)

#Construct the cluster tree and merge nodes based on visual inspection
donor02 <- BuildClusterTree(object = donor02, do.reorder = TRUE, reorder.numeric = TRUE, do.plot = T)
donor02 <- MergeNode(donor02, node.use = 19, rebuild.tree = T)
TSNEPlot(object = donor02, do.label = T)

#Subset Cluster 1 (Ciliated and Club Cells)
donor02.cluster01 <- SubsetData(donor02, ident.use = c(1))

donor02.cluster01 <- ScaleData(object = donor02.cluster01, vars.to.regress = c("nUMI", "percent.mito"))

donor02.cluster01 <- FindVariableGenes(object = donor02.cluster01, mean.function = ExpMean, dispersion.function = LogVMR, 
                                    x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

length(x = donor02.cluster01@var.genes)

donor02.cluster01 <- RunPCA(object = donor02.cluster01, pc.genes = donor02.cluster01@var.genes, do.print = TRUE, pcs.print = 1:5, 
                         genes.print = 5)

donor02.cluster01 <- ProjectPCA(object = donor02.cluster01, do.print = FALSE)

PCHeatmap(object = donor02.cluster01, pc.use = 1:12, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCElbowPlot(object = donor02.cluster01, num.pc = 12)

#Clustering of Subsetted Cells
donor02.cluster01 <- RunTSNE(object = donor02.cluster01, dims.use = 1:4, do.fast = TRUE)

donor02.cluster01 <- FindClusters(object = donor02.cluster01, reduction.type = "pca", dims.use = 1:4, 
                               resolution = 0.5, save.SNN = TRUE)
TSNEPlot(donor02.cluster01)

#Markers of Subsetted Cell Subclusters
donor02.cluster01.markers <- FindAllMarkers(object = donor02.cluster01, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

donor02.cluster01.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)

FeaturePlot(object = donor02.cluster01, features.plot = c("TPPP3", "SCGB3A2"), min.cutoff = "q9",
            cols.use = c("lightgrey", "blue"), pt.size = 1)

#Cell assignment for Subclusters
donor02.ciliated <- WhichCells(object = donor02.cluster01, ident = c(0))
donor02.clubs <- WhichCells(object = donor02.cluster01, ident = c(1))

#Rename clusters
donor02.ident <- c(1, 2, 3, 4, 5, 6, 7, 9, 10)

donor02.new.ident <- c("Ciliated and Club Cells", "T/NKT Cells", "Monocytes", "Macrophages", 
                     "Dendritic Cells", "AT1 Cells", "AT2 Cells", "Endothelial Cells", "Fibroblasts")

donor02@ident <- plyr::mapvalues(x = donor02@ident, from = donor02.ident, to = donor02.new.ident)

donor02 <- SetIdent(object = donor02, cells.use = donor02.ciliated, ident.use = "Ciliated Cells")
donor02 <- SetIdent(object = donor02, cells.use = donor02.clubs, ident.use = "Club Cells")

donor02@ident <- ordered(donor02@ident,levels = c("Macrophages", "AT2 Cells", "Monocytes", "Club Cells",
                                            "Dendritic Cells", "T/NKT Cells", "AT1 Cells", "Ciliated Cells",
                                            "Endothelial Cells", "Fibroblasts"))

#t-SNE Plot (Figure E5)
TSNEPlot(object = donor02, do.label = T)

#Find marker genes for clusters (Table E5)
donor02.markers <- FindAllMarkers(object = donor02, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.table(donor02.markers %>% group_by(cluster), file= "/donor02.markers.txt", row.names=FALSE, sep="\t")

#Add Cluster Identity Metadata to Separate Metadata Column
donor02 <- AddMetaData(object = donor02, metadata = donor02@ident, col.name = "indiv.cell.ident")

donor02@meta.data[, "diagnosis"] <- "donor"
donor02@meta.data[, "condition"] <- "donor"

donor02@meta.data$res.0.3 <- NULL
donor02@meta.data$tree.ident <- NULL

#Save object
save(donor02, file = "/donor02.Robj")
load(file = "/donor02.Robj")


# donor03 analysis ---------------------------------------------------------

# Load the donor03 dataset
donor03.data <- Read10X(data.dir = "/donor03_tables/")

donor03 <- CreateSeuratObject(raw.data = donor03.data, min.cells = 3, min.genes = 200, 
                              project = "donor03")

# Identify mitochondrial genes
mito.genes <- grep(pattern = "^MT-", x = rownames(x = donor03@data), value = TRUE)

# Create meta data column containing percentage of mitochondrial genes for each cell
percent.mito <- Matrix::colSums(donor03@raw.data[mito.genes, ])/Matrix::colSums(donor03@raw.data)

# Add percent.mito to object meta data
donor03 <- AddMetaData(object = donor03, metadata = percent.mito, col.name = "percent.mito")

#Data inspection and basic filtering
VlnPlot(object = donor03, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
GenePlot(object = donor03, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = donor03, gene1 = "nUMI", gene2 = "nGene")

donor03 <- FilterCells(object = donor03, subset.names = c("nGene", "percent.mito"), 
                       low.thresholds = c(200, -Inf), high.thresholds = c(5000, 0.1))

#Normalize, find variable genes and scale
donor03 <- NormalizeData(object = donor03, normalization.method = "LogNormalize", 
                         scale.factor = 10000)
donor03 <- FindVariableGenes(object = donor03, mean.function = ExpMean, dispersion.function = LogVMR, 
                             x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = donor03@var.genes)
donor03 <- ScaleData(object = donor03, vars.to.regress = c("nUMI", "percent.mito"))

#Perform PCA
donor03 <- RunPCA(object = donor03, pc.genes = donor03@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5, pcs.compute = 40)
PrintPCA(object = donor03, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(object = donor03, pcs.use = 1:2)
PCAPlot(object = donor03, dim.1 = 1, dim.2 = 2)
donor03 <- ProjectPCA(object = donor03, do.print = FALSE)
PCHeatmap(object = donor03, pc.use = 1:9, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = donor03, pc.use = 10:18, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = donor03, pc.use = 19:27, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCElbowPlot(object = donor03, num.pc = 40)

#Perform clustering using first 9 principal components. 
donor03 <- FindClusters(object = donor03, reduction.type = "pca", dims.use = 1:9, 
                        resolution = 0.3, print.output = 0, save.SNN = TRUE, force.recalc = T)
PrintFindClustersParams(object = donor03)
donor03 <- RunTSNE(object = donor03, dims.use = 1:9, do.fast = TRUE)
TSNEPlot(object = donor03, do.label = T)

#Construct the cluster tree and merge nodes based on visual inspection
donor03 <- BuildClusterTree(object = donor03, do.reorder = TRUE, reorder.numeric = TRUE, do.plot = T)
donor03 <- MergeNode(donor03, node.use = 19, rebuild.tree = T)
TSNEPlot(object = donor03, do.label = T)

#Subset Cluster 2 (Ciliated and Club Cells)
donor03.cluster02 <- SubsetData(donor03, ident.use = c(2))

donor03.cluster02 <- ScaleData(object = donor03.cluster02, vars.to.regress = c("nUMI", "percent.mito"))

donor03.cluster02 <- FindVariableGenes(object = donor03.cluster02, mean.function = ExpMean, dispersion.function = LogVMR, 
                                    x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

length(x = donor03.cluster02@var.genes)

donor03.cluster02 <- RunPCA(object = donor03.cluster02, pc.genes = donor03.cluster02@var.genes, do.print = TRUE, pcs.print = 1:5, 
                         genes.print = 5)

donor03.cluster02 <- ProjectPCA(object = donor03.cluster02, do.print = FALSE)

PCHeatmap(object = donor03.cluster02, pc.use = 1:12, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)

PCElbowPlot(object = donor03.cluster02, num.pc = 12)

#Clustering of Subsetted Cells
donor03.cluster02 <- RunTSNE(object = donor03.cluster02, dims.use = 1:2, do.fast = TRUE)

donor03.cluster02 <- FindClusters(object = donor03.cluster02, reduction.type = "pca", dims.use = 1:2, 
                               resolution = 0.2, save.SNN = TRUE)
TSNEPlot(donor03.cluster02)

#Markers of Subsetted Cell Subclusters
donor03.cluster02.markers <- FindAllMarkers(object = donor03.cluster02, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

donor03.cluster02.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)

FeaturePlot(object = donor03.cluster02, features.plot = c("TPPP3", "SCGB3A2"), min.cutoff = "q9",
            cols.use = c("lightgrey", "blue"), pt.size = 1)

#Cell assignment for Subclusters
donor03.ciliated <- WhichCells(object = donor03.cluster02, ident = c(0))
donor03.clubs <- WhichCells(object = donor03.cluster02, ident = c(1))

#Rename clusters
donor03.ident <- c(0, 1, 2, 3, 4, 5, 6, 7, 8)

donor03.new.ident <- c("Macrophages", "AT2 Cells", "Ciliated and Club Cells", "Monocytes", "Plasma Cells",
                     "AT1 Cells", "Endothelial Cells", "Unassigned", "Fibroblasts")

donor03@ident <- plyr::mapvalues(x = donor03@ident, from = donor03.ident, to = donor03.new.ident)

donor03 <- SetIdent(object = donor03, cells.use = donor03.clubs, ident.use = "Club Cells")
donor03 <- SetIdent(object = donor03, cells.use = donor03.ciliated, ident.use = "Ciliated Cells")

donor03@ident <- ordered(donor03@ident,levels = c("Macrophages", "AT2 Cells", "Monocytes", "Club Cells",
                                            "Plasma Cells", "AT1 Cells", "Ciliated Cells", "Endothelial Cells",  
                                            "Fibroblasts", "Unassigned"))

#t-SNE Plot (Figure E5)
TSNEPlot(object = donor03, do.label = T)

#Find marker genes for clusters (Table E5)
donor03.markers <- FindAllMarkers(object = donor03, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.table(donor03.markers %>% group_by(cluster), file= "/donor03.markers.txt", row.names=FALSE, sep="\t")

#Add Cluster Identity Metadata to Separate Metadata Column
donor03 <- AddMetaData(object = donor03, metadata = donor03@ident, col.name = "indiv.cell.ident")

donor03@meta.data[, "diagnosis"] <- "donor"
donor03@meta.data[, "condition"] <- "donor"

donor03@meta.data$res.0.3 <- NULL
donor03@meta.data$tree.ident <- NULL

#Save object
save(donor03, file = "/donor03.Robj")
load(file = "/donor03.Robj")

# donor04 analysis ---------------------------------------------------------

# Load the donor04 dataset
donor04.data <- Read10X(data.dir = "/donor04_tables/")

donor04 <- CreateSeuratObject(raw.data = donor04.data, min.cells = 3, min.genes = 200, 
                              project = "donor04")

# Identify mitochondrial genes
mito.genes <- grep(pattern = "^MT-", x = rownames(x = donor04@data), value = TRUE)

# Create meta data column containing percentage of mitochondrial genes for each cell
percent.mito <- Matrix::colSums(donor04@raw.data[mito.genes, ])/Matrix::colSums(donor04@raw.data)

# Add percent.mito to object meta data
donor04 <- AddMetaData(object = donor04, metadata = percent.mito, col.name = "percent.mito")

#Data inspection and basic filtering
VlnPlot(object = donor04, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
GenePlot(object = donor04, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = donor04, gene1 = "nUMI", gene2 = "nGene")

donor04 <- FilterCells(object = donor04, subset.names = c("nGene", "percent.mito"), 
                       low.thresholds = c(200, -Inf), high.thresholds = c(5000, 0.15))

#Normalize, find variable genes and scale
donor04 <- NormalizeData(object = donor04, normalization.method = "LogNormalize", 
                         scale.factor = 10000)
donor04 <- FindVariableGenes(object = donor04, mean.function = ExpMean, dispersion.function = LogVMR, 
                             x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = donor04@var.genes)
donor04 <- ScaleData(object = donor04, vars.to.regress = c("nUMI", "percent.mito"))

#Perform PCA
donor04 <- RunPCA(object = donor04, pc.genes = donor04@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5, pcs.compute = 40)
PrintPCA(object = donor04, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(object = donor04, pcs.use = 1:2)
PCAPlot(object = donor04, dim.1 = 1, dim.2 = 2)
donor04 <- ProjectPCA(object = donor04, do.print = FALSE)
PCHeatmap(object = donor04, pc.use = 1:9, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = donor04, pc.use = 10:18, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = donor04, pc.use = 19:27, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCElbowPlot(object = donor04, num.pc = 40)

#Perform clustering using first 9 principal components. 
donor04 <- FindClusters(object = donor04, reduction.type = "pca", dims.use = 1:9, 
                        resolution = 0.3, print.output = 0, save.SNN = TRUE, force.recalc = T)
PrintFindClustersParams(object = donor04)
donor04 <- RunTSNE(object = donor04, dims.use = 1:9, do.fast = TRUE)
TSNEPlot(object = donor04, do.label = T)

#Subset Cluster 4 (Mixed Epithelial Cells)
donor04.cluster04 <- SubsetData(donor04, ident.use = c(4))

donor04.cluster04 <- ScaleData(object = donor04.cluster04, vars.to.regress = c("nUMI", "percent.mito"))

donor04.cluster04 <- FindVariableGenes(object = donor04.cluster04, mean.function = ExpMean, dispersion.function = LogVMR, 
                                    x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

length(x = donor04.cluster04@var.genes)

donor04.cluster04 <- RunPCA(object = donor04.cluster04, pc.genes = donor04.cluster04@var.genes, do.print = TRUE, pcs.print = 1:5, 
                         genes.print = 5)

donor04.cluster04 <- ProjectPCA(object = donor04.cluster04, do.print = FALSE)

PCHeatmap(object = donor04.cluster04, pc.use = 1:12, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)

PCElbowPlot(object = donor04.cluster04, num.pc = 12)

#Clustering of Subsetted Cells
donor04.cluster04 <- RunTSNE(object = donor04.cluster04, dims.use = 1:4, do.fast = TRUE)

donor04.cluster04 <- FindClusters(object = donor04.cluster04, reduction.type = "pca", dims.use = 1:4, 
                               resolution = .5, save.SNN = TRUE)
TSNEPlot(donor04.cluster04)

#Markers of Subsetted Cell Subclusters
donor04.cluster04.markers <- FindAllMarkers(object = donor04.cluster04, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

donor04.cluster04.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)

FeaturePlot(object = donor04.cluster04, features.plot = c("SFTPC", "AGER", "SCGB3A2"), min.cutoff = "q9",
            cols.use = c("lightgrey", "blue"), pt.size = 1)

#Cell assignment for Subclusters
donor04.at2 <- WhichCells(object = donor04.cluster04, ident = c(0))
donor04.at1 <- WhichCells(object = donor04.cluster04, ident = c(1))
donor04.clubs <- WhichCells(object = donor04.cluster04, ident = c(2))

#Rename clusters
donor04.ident <- c(0, 1, 2, 3, 4, 5, 6, 7)

donor04.new.ident <- c("AT2 Cells", "Macrophages", "AT2 Cells", "Macrophages",
                     "AT1, AT2, and Club Cells", "Ciliated Cells", "Unassigned", "Plasma Cells")

donor04@ident <- plyr::mapvalues(x = donor04@ident, from = donor04.ident, to = donor04.new.ident)

donor04 <- SetIdent(object = donor04, cells.use = donor04.at2, ident.use = "AT2 Cells")
donor04 <- SetIdent(object = donor04, cells.use = donor04.at1, ident.use = "AT1 Cells")
donor04 <- SetIdent(object = donor04, cells.use = donor04.clubs, ident.use = "Club Cells")

donor04@ident <- ordered(donor04@ident,levels = c("Macrophages", "AT2 Cells", "Club Cells", "Plasma Cells",
                                            "AT1 Cells", "Ciliated Cells", "Unassigned"))

donor04.ident <- c("Macrophages", "AT2 Cells", "Club Cells", "Plasma Cells",
                 "AT1 Cells", "Ciliated Cells", "Unassigned")

donor04.new.ident <- c("Macrophages", "AT2 Cells", "Club Cells", "Plasma Cells",
                     "AT1 Cells", "Ciliated Cells", "Unassigned")

donor04@ident <- plyr::mapvalues(x = donor04@ident, from = donor04.ident, to = donor04.new.ident)

#t-SNE Plot (Figure E5)
TSNEPlot(object = donor04, do.label = T)

#Find marker genes for clusters (Table E5)
donor04.markers <- FindAllMarkers(object = donor04, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.table(donor04.markers %>% group_by(cluster), file= "/donor04.markers.txt", row.names=FALSE, sep="\t")

#Add Cluster Identity Metadata to Separate Metadata Column
donor04 <- AddMetaData(object = donor04, metadata = donor04@ident, col.name = "indiv.cell.ident")

donor04@meta.data[, "diagnosis"] <- "donor"
donor04@meta.data[, "condition"] <- "donor"

donor04@meta.data$res.0.3 <- NULL
donor04@meta.data$tree.ident <- NULL

#Save object
save(donor04, file = "/donor04.Robj")
load(file = "/donor04.Robj")

# donor05 analysis ---------------------------------------------------------

# Load the donor05 dataset
donor05.data <- Read10X(data.dir = "/donor05_tables/")

donor05 <- CreateSeuratObject(raw.data = donor05.data, min.cells = 3, min.genes = 200, 
                              project = "donor05")

# Identify mitochondrial genes
mito.genes <- grep(pattern = "^MT-", x = rownames(x = donor05@data), value = TRUE)

# Create meta data column containing percentage of mitochondrial genes for each cell
percent.mito <- Matrix::colSums(donor05@raw.data[mito.genes, ])/Matrix::colSums(donor05@raw.data)

# Add percent.mito to object meta data
donor05 <- AddMetaData(object = donor05, metadata = percent.mito, col.name = "percent.mito")

#Data inspection and basic filtering
VlnPlot(object = donor05, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
GenePlot(object = donor05, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = donor05, gene1 = "nUMI", gene2 = "nGene")

donor05 <- FilterCells(object = donor05, subset.names = c("nGene", "percent.mito"), 
                       low.thresholds = c(200, -Inf), high.thresholds = c(6000, 0.1))

#Normalize, find variable genes and scale
donor05 <- NormalizeData(object = donor05, normalization.method = "LogNormalize", 
                         scale.factor = 10000)
donor05 <- FindVariableGenes(object = donor05, mean.function = ExpMean, dispersion.function = LogVMR, 
                             x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = donor05@var.genes)
donor05 <- ScaleData(object = donor05, vars.to.regress = c("nUMI", "percent.mito"))

#Perform PCA
donor05 <- RunPCA(object = donor05, pc.genes = donor05@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5, pcs.compute = 40)
PrintPCA(object = donor05, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(object = donor05, pcs.use = 1:2)
PCAPlot(object = donor05, dim.1 = 1, dim.2 = 2)
donor05 <- ProjectPCA(object = donor05, do.print = FALSE)
PCHeatmap(object = donor05, pc.use = 1:9, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = donor05, pc.use = 10:18, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = donor05, pc.use = 19:27, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCElbowPlot(object = donor05, num.pc = 40)

#Perform clustering using first 27 principal components. 
donor05 <- FindClusters(object = donor05, reduction.type = "pca", dims.use = 1:27, 
                        resolution = 0.5, print.output = 0, save.SNN = TRUE, force.recalc = T)
PrintFindClustersParams(object = donor05)
donor05 <- RunTSNE(object = donor05, dims.use = 1:27, do.fast = TRUE)
TSNEPlot(object = donor05, do.label = T)

#Subset Cluster 3 (Epithelial Cells)
donor05.cluster03 <- SubsetData(donor05, ident.use = c(3))

donor05.cluster03 <- ScaleData(object = donor05.cluster03, vars.to.regress = c("nUMI", "percent.mito"))

donor05.cluster03 <- FindVariableGenes(object = donor05.cluster03, mean.function = ExpMean, dispersion.function = LogVMR, 
                                       x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

length(x = donor05.cluster03@var.genes)

donor05.cluster03 <- RunPCA(object = donor05.cluster03, pc.genes = donor05.cluster03@var.genes, do.print = TRUE, pcs.print = 1:5, 
                            genes.print = 5)

donor05.cluster03 <- ProjectPCA(object = donor05.cluster03, do.print = FALSE)

PCHeatmap(object = donor05.cluster03, pc.use = 1:12, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)

PCElbowPlot(object = donor05.cluster03, num.pc = 12)

#Clustering of Subsetted Cells
donor05.cluster03 <- RunTSNE(object = donor05.cluster03, dims.use = 1:7, do.fast = TRUE)

donor05.cluster03 <- FindClusters(object = donor05.cluster03, reduction.type = "pca", dims.use = 1:7, 
                                  resolution = .3, save.SNN = TRUE)
TSNEPlot(donor05.cluster03)

#Markers of Subsetted Cell Subclusters
donor05.cluster03.markers <- FindAllMarkers(object = donor05.cluster03, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

donor05.cluster03.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)

FeaturePlot(object = donor05.cluster03, features.plot = c("SFTPC", "AGER", "SCGB3A2", "CD68"), min.cutoff = "q9",
            cols.use = c("lightgrey", "blue"), pt.size = 1)

#Cell assignment for Subclusters
donor05.at2 <- WhichCells(object = donor05.cluster03, ident = c(0))
donor05.at1 <- WhichCells(object = donor05.cluster03, ident = c(1))
donor05.clubs <- WhichCells(object = donor05.cluster03, ident = c(2))
donor05.cluster03.unassigned <- WhichCells(object = donor05.cluster03, ident = c(3))

#Subset Cluster 5 (Monocytes and DCs)
donor05.cluster05 <- SubsetData(donor05, ident.use = c(5))

donor05.cluster05 <- ScaleData(object = donor05.cluster05, vars.to.regress = c("nUMI", "percent.mito"))

donor05.cluster05 <- FindVariableGenes(object = donor05.cluster05, mean.function = ExpMean, dispersion.function = LogVMR, 
                                       x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

length(x = donor05.cluster05@var.genes)

donor05.cluster05 <- RunPCA(object = donor05.cluster05, pc.genes = donor05.cluster05@var.genes, do.print = TRUE, pcs.print = 1:5, 
                            genes.print = 5)

donor05.cluster05 <- ProjectPCA(object = donor05.cluster05, do.print = FALSE)

PCHeatmap(object = donor05.cluster05, pc.use = 1:12, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)

PCElbowPlot(object = donor05.cluster05, num.pc = 12)

#Clustering of Subsetted Cells
donor05.cluster05 <- RunTSNE(object = donor05.cluster05, dims.use = 1:5, do.fast = TRUE)

donor05.cluster05 <- FindClusters(object = donor05.cluster05, reduction.type = "pca", dims.use = 1:5, 
                                  resolution = .3, save.SNN = TRUE)
TSNEPlot(donor05.cluster05)

#Markers of Subsetted Cell Subclusters
donor05.cluster05.markers <- FindAllMarkers(object = donor05.cluster05, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

donor05.cluster05.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)

FeaturePlot(object = donor05.cluster05, features.plot = c("HLA-DPB1", "CD14", "CD1D"), min.cutoff = "q9",
            cols.use = c("lightgrey", "blue"), pt.size = 1)

#Cell assignment for Subclusters
donor05.dcs <- WhichCells(object = donor05.cluster05, ident = c(0))
donor05.monos <- WhichCells(object = donor05.cluster05, ident = c(1))
donor05.cluster05.unassigned <- WhichCells(object = donor05.cluster05, ident = c(2))

#Rename clusters
donor05.ident <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)

donor05.new.ident <- c("Macrophages", "AT2 Cells", "AT2 Cells", "Epithelial Cells", 
                       "AT2 Cells", "Monocytes and DCs", "Unassigned", "Endothelial Cells",
                       "Fibroblasts", "Plasma Cells", "Ciliated Cells", "Mast Cells")

donor05@ident <- plyr::mapvalues(x = donor05@ident, from = donor05.ident, to = donor05.new.ident)

donor05 <- SetIdent(object = donor05, cells.use = donor05.at2, ident.use = "AT2 Cells")
donor05 <- SetIdent(object = donor05, cells.use = donor05.at1, ident.use = "AT1 Cells")
donor05 <- SetIdent(object = donor05, cells.use = donor05.clubs, ident.use = "Club Cells")
donor05 <- SetIdent(object = donor05, cells.use = donor05.dcs, ident.use = "Dendritic Cells")
donor05 <- SetIdent(object = donor05, cells.use = donor05.monos, ident.use = "Monocytes")
donor05 <- SetIdent(object = donor05, cells.use = donor05.cluster03.unassigned, ident.use = "Unassigned")
donor05 <- SetIdent(object = donor05, cells.use = donor05.cluster05.unassigned, ident.use = "Unassigned")

donor05@ident <- ordered(donor05@ident,levels = c("Macrophages", "AT2 Cells", "Monocytes", "Club Cells",
                                                  "Dendritic Cells", "Plasma Cells", "AT1 Cells", "Ciliated Cells", 
                                                  "Endothelial Cells", "Fibroblasts", "Mast Cells", "Unassigned"))

#t-SNE Plot (Figure E5)
TSNEPlot(object = donor05, do.label = T)

#Find marker genes for clusters (Table E5)
donor05.markers <- FindAllMarkers(object = donor05, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.table(donor05.markers %>% group_by(cluster), file= "/donor05.markers.txt", row.names=FALSE, sep="\t")

#Add Cluster Identity Metadata to Separate Metadata Column
donor05 <- AddMetaData(object = donor05, metadata = donor05@ident, col.name = "indiv.cell.ident")

donor05@meta.data[, "diagnosis"] <- "donor"
donor05@meta.data[, "condition"] <- "donor"

donor05@meta.data$res.0.3 <- NULL
donor05@meta.data$tree.ident <- NULL

#Save object
save(donor05, file = "/donor05.Robj")
load(file = "/donor05.Robj")

# donor06 analysis ---------------------------------------------------------

# Load the donor06 dataset
donor06.data <- Read10X(data.dir = "/donor06_tables/")

donor06 <- CreateSeuratObject(raw.data = donor06.data, min.cells = 3, min.genes = 200, 
                              project = "donor06")

# Identify mitochondrial genes
mito.genes <- grep(pattern = "^MT-", x = rownames(x = donor06@data), value = TRUE)

# Create meta data column containing percentage of mitochondrial genes for each cell
percent.mito <- Matrix::colSums(donor06@raw.data[mito.genes, ])/Matrix::colSums(donor06@raw.data)

# Add percent.mito to object meta data
donor06 <- AddMetaData(object = donor06, metadata = percent.mito, col.name = "percent.mito")

#Data inspection and basic filtering
VlnPlot(object = donor06, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
GenePlot(object = donor06, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = donor06, gene1 = "nUMI", gene2 = "nGene")

donor06 <- FilterCells(object = donor06, subset.names = c("nGene", "percent.mito"), 
                       low.thresholds = c(200, -Inf), high.thresholds = c(6000, 0.1))

#Normalize, find variable genes and scale
donor06 <- NormalizeData(object = donor06, normalization.method = "LogNormalize", 
                         scale.factor = 10000)
donor06 <- FindVariableGenes(object = donor06, mean.function = ExpMean, dispersion.function = LogVMR, 
                             x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = donor06@var.genes)
donor06 <- ScaleData(object = donor06, vars.to.regress = c("nUMI", "percent.mito"))

#Perform PCA
donor06 <- RunPCA(object = donor06, pc.genes = donor06@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5, pcs.compute = 40)
PrintPCA(object = donor06, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(object = donor06, pcs.use = 1:2)
PCAPlot(object = donor06, dim.1 = 1, dim.2 = 2)
donor06 <- ProjectPCA(object = donor06, do.print = FALSE)
PCHeatmap(object = donor06, pc.use = 1:9, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = donor06, pc.use = 10:18, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = donor06, pc.use = 19:27, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCElbowPlot(object = donor06, num.pc = 40)

#Perform clustering using first 20 principal components. 
donor06 <- FindClusters(object = donor06, reduction.type = "pca", dims.use = 1:20, 
                        resolution = 0.5, print.output = 0, save.SNN = TRUE, force.recalc = T)
PrintFindClustersParams(object = donor06)
donor06 <- RunTSNE(object = donor06, dims.use = 1:20, do.fast = TRUE)
TSNEPlot(object = donor06, do.label = T)

#Subset Cluster 2 (Monocytes and Macrophages)
donor06.cluster02 <- SubsetData(donor06, ident.use = c(2))

donor06.cluster02 <- ScaleData(object = donor06.cluster02, vars.to.regress = c("nUMI", "percent.mito"))

donor06.cluster02 <- FindVariableGenes(object = donor06.cluster02, mean.function = ExpMean, dispersion.function = LogVMR, 
                                       x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

length(x = donor06.cluster02@var.genes)

donor06.cluster02 <- RunPCA(object = donor06.cluster02, pc.genes = donor06.cluster02@var.genes, do.print = TRUE, pcs.print = 1:5, 
                            genes.print = 5)

donor06.cluster02 <- ProjectPCA(object = donor06.cluster02, do.print = FALSE)

PCHeatmap(object = donor06.cluster02, pc.use = 1:12, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)

PCElbowPlot(object = donor06.cluster02, num.pc = 12)

#Clustering of Subsetted Cells
donor06.cluster02 <- RunTSNE(object = donor06.cluster02, dims.use = 1:6, do.fast = TRUE)

donor06.cluster02 <- FindClusters(object = donor06.cluster02, reduction.type = "pca", dims.use = 1:6, 
                                  resolution = .2, save.SNN = TRUE)
TSNEPlot(donor06.cluster02)

#Markers of Subsetted Cell Subclusters
donor06.cluster02.markers <- FindAllMarkers(object = donor06.cluster02, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

donor06.cluster02.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)

FeaturePlot(object = donor06.cluster02, features.plot = c("FCN1", "MRC1"), min.cutoff = "q9",
            cols.use = c("lightgrey", "blue"), pt.size = 1)

#Cell assignment for Subclusters
donor06.monos <- WhichCells(object = donor06.cluster02, ident = c(0))
donor06.macrophages <- WhichCells(object = donor06.cluster02, ident = c(1))

#Subset Cluster 12 (Club and Ciliated Cells)
donor06.cluster12 <- SubsetData(donor06, ident.use = c(12))

donor06.cluster12 <- ScaleData(object = donor06.cluster12, vars.to.regress = c("nUMI", "percent.mito"))

donor06.cluster12 <- FindVariableGenes(object = donor06.cluster12, mean.function = ExpMean, dispersion.function = LogVMR, 
                                       x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

length(x = donor06.cluster12@var.genes)

donor06.cluster12 <- RunPCA(object = donor06.cluster12, pc.genes = donor06.cluster12@var.genes, do.print = TRUE, pcs.print = 1:5, 
                            genes.print = 5)

donor06.cluster12 <- ProjectPCA(object = donor06.cluster12, do.print = FALSE)

PCHeatmap(object = donor06.cluster12, pc.use = 1:12, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)

PCElbowPlot(object = donor06.cluster12, num.pc = 12)

#Clustering of Subsetted Cells
donor06.cluster12 <- RunTSNE(object = donor06.cluster12, dims.use = 1:8, do.fast = TRUE)

donor06.cluster12 <- FindClusters(object = donor06.cluster12, reduction.type = "pca", dims.use = 1:8, 
                                  resolution = 1, save.SNN = TRUE)
TSNEPlot(donor06.cluster12)

#Markers of Subsetted Cell Subclusters
donor06.cluster12.markers <- FindAllMarkers(object = donor06.cluster12, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

donor06.cluster12.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)

FeaturePlot(object = donor06.cluster12, features.plot = c("SCGB3A2", "TPPP3"), min.cutoff = "q9",
            cols.use = c("lightgrey", "blue"), pt.size = 1)

#Cell assignment for Subclusters
donor06.clubs <- WhichCells(object = donor06.cluster12, ident = c(0))
donor06.ciliated <- WhichCells(object = donor06.cluster12, ident = c(1))

#Rename clusters
donor06.ident <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)

donor06.new.ident <- c("Macrophages", "AT2 Cells", "Monocytes and Macrophages", "AT2 Cells",
                       "AT1 Cells", "AT2 Cells", "T/NKT Cells", "Macrophages",
                       "Dendritic Cells", "Endothelial Cells", "Unassigned", "Plasma Cells",
                       "Club and Ciliated Cells", "Fibroblasts")

donor06@ident <- plyr::mapvalues(x = donor06@ident, from = donor06.ident, to = donor06.new.ident)

donor06 <- SetIdent(object = donor06, cells.use = donor06.monos, ident.use = "Monocytes")
donor06 <- SetIdent(object = donor06, cells.use = donor06.macrophages, ident.use = "Macrophages")
donor06 <- SetIdent(object = donor06, cells.use = donor06.clubs, ident.use = "Club Cells")
donor06 <- SetIdent(object = donor06, cells.use = donor06.ciliated, ident.use = "Ciliated Cells")

donor06@ident <- ordered(donor06@ident,levels = c("Macrophages", "AT2 Cells", "Monocytes", "Club Cells",
                                                  "Dendritic Cells", "Plasma Cells", "T/NKT Cells", "AT1 Cells", 
                                                  "Ciliated Cells", "Endothelial Cells", "Fibroblasts", "Unassigned"))

#t-SNE Plot (Figure E5)
TSNEPlot(object = donor06, do.label = T)

#Find marker genes for clusters (Table E5)
donor06.markers <- FindAllMarkers(object = donor06, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.table(donor06.markers %>% group_by(cluster), file= "/donor06.markers.txt", row.names=FALSE, sep="\t")

#Monocyte Heterogeneity Feature Plots (Figure E6C)
FeaturePlot(object = donor06, features.plot = c("CD14", "FCGR3A", "MRC1"), min.cutoff = "q9",
            cols.use = c("lightgrey", "blue"), nCol = 1, pt.size = 0.5)

#Macrophage Heterogeneity Feature Plots (Figure E7F)
FeaturePlot(object = donor06, features.plot = c("MRC1", "FABP4", "CCL3", "CD14"), min.cutoff = "q9",
            cols.use = c("lightgrey", "blue"), nCol = 1, pt.size = 0.5)

#Add Cluster Identity Metadata to Separate Metadata Column
donor06 <- AddMetaData(object = donor06, metadata = donor06@ident, col.name = "indiv.cell.ident")

donor06@meta.data[, "diagnosis"] <- "donor"
donor06@meta.data[, "condition"] <- "donor"

donor06@meta.data$res.0.3 <- NULL
donor06@meta.data$tree.ident <- NULL

#Save object
save(donor06, file = "/donor06.Robj")
load(file = "/donor06.Robj")

# donor07 analysis ---------------------------------------------------------

# Load the donor07 dataset
donor07.data <- Read10X(data.dir = "/donor07_tables/")

donor07 <- CreateSeuratObject(raw.data = donor07.data, min.cells = 3, min.genes = 200, 
                              project = "donor07")

# Identify mitochondrial genes
mito.genes <- grep(pattern = "^MT-", x = rownames(x = donor07@data), value = TRUE)

# Create meta data column containing percentage of mitochondrial genes for each cell
percent.mito <- Matrix::colSums(donor07@raw.data[mito.genes, ])/Matrix::colSums(donor07@raw.data)

# Add percent.mito to object meta data
donor07 <- AddMetaData(object = donor07, metadata = percent.mito, col.name = "percent.mito")

#Data inspection and basic filtering
VlnPlot(object = donor07, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
GenePlot(object = donor07, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = donor07, gene1 = "nUMI", gene2 = "nGene")

donor07 <- FilterCells(object = donor07, subset.names = c("nGene", "percent.mito"), 
                       low.thresholds = c(200, -Inf), high.thresholds = c(5000, 0.1))

#Normalize, find variable genes and scale
donor07 <- NormalizeData(object = donor07, normalization.method = "LogNormalize", 
                         scale.factor = 10000)
donor07 <- FindVariableGenes(object = donor07, mean.function = ExpMean, dispersion.function = LogVMR, 
                             x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = donor07@var.genes)
donor07 <- ScaleData(object = donor07, vars.to.regress = c("nUMI", "percent.mito"))

#Perform PCA
donor07 <- RunPCA(object = donor07, pc.genes = donor07@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5, pcs.compute = 40)
PrintPCA(object = donor07, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(object = donor07, pcs.use = 1:2)
PCAPlot(object = donor07, dim.1 = 1, dim.2 = 2)
donor07 <- ProjectPCA(object = donor07, do.print = FALSE)
PCHeatmap(object = donor07, pc.use = 1:9, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = donor07, pc.use = 10:18, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = donor07, pc.use = 19:27, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCElbowPlot(object = donor07, num.pc = 40)

#Perform clustering using first 22 principal components. 
donor07 <- FindClusters(object = donor07, reduction.type = "pca", dims.use = 1:22, 
                        resolution = 0.3, print.output = 0, save.SNN = TRUE, force.recalc = T)
PrintFindClustersParams(object = donor07)
donor07 <- RunTSNE(object = donor07, dims.use = 1:20, do.fast = TRUE)
TSNEPlot(object = donor07, do.label = T)

#Subset Cluster 06 (Club and Airway Cells)
donor07.cluster06 <- SubsetData(donor07, ident.use = c(6))

donor07.cluster06 <- ScaleData(object = donor07.cluster06, vars.to.regress = c("nUMI", "percent.mito"))

donor07.cluster06 <- FindVariableGenes(object = donor07.cluster06, mean.function = ExpMean, dispersion.function = LogVMR, 
                                       x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

length(x = donor07.cluster06@var.genes)

donor07.cluster06 <- RunPCA(object = donor07.cluster06, pc.genes = donor07.cluster06@var.genes, do.print = TRUE, pcs.print = 1:5, 
                            genes.print = 5)

donor07.cluster06 <- ProjectPCA(object = donor07.cluster06, do.print = FALSE)

PCHeatmap(object = donor07.cluster06, pc.use = 1:12, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)

PCElbowPlot(object = donor07.cluster06, num.pc = 12)

#Clustering of Subsetted Cells
donor07.cluster06 <- RunTSNE(object = donor07.cluster06, dims.use = 1:5, do.fast = TRUE)

donor07.cluster06 <- FindClusters(object = donor07.cluster06, reduction.type = "pca", dims.use = 1:5, 
                                  resolution = 0.5, save.SNN = TRUE)
TSNEPlot(donor07.cluster06)

#Markers of Subsetted Cell Subclusters
donor07.cluster06.markers <- FindAllMarkers(object = donor07.cluster06, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

donor07.cluster06.markers %>% group_by(cluster) %>% top_n(3, avg_logFC)

FeaturePlot(object = donor07.cluster06, features.plot = c("SCGB3A2", "TPPP3"), min.cutoff = "q9",
            cols.use = c("lightgrey", "blue"), pt.size = 1)

#Cell assignment for Subclusters
donor07.ciliated <- WhichCells(object = donor07.cluster06, ident = c(0))
donor07.clubs <- WhichCells(object = donor07.cluster06, ident = c(1))

#Rename clusters
donor07.ident <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)

donor07.new.ident <- c("AT2 Cells", "Macrophages", "AT2 Cells", "Monocytes", 
                       "AT2 Cells", "Macrophages", "Club and Ciliated Cells", "T/NKT Cells",
                       "Endothelial Cells", "Unassigned", "Lymphatic Cells", "Fibroblasts",
                       "Mast Cells")

donor07@ident <- plyr::mapvalues(x = donor07@ident, from = donor07.ident, to = donor07.new.ident)

donor07 <- SetIdent(object = donor07, cells.use = donor07.clubs, ident.use = "Club Cells")
donor07 <- SetIdent(object = donor07, cells.use = donor07.ciliated, ident.use = "Ciliated Cells")

donor07@ident <- ordered(donor07@ident,levels = c("Macrophages", "AT2 Cells", "Monocytes", "Club Cells",
                                                  "T/NKT Cells", "Ciliated Cells", "Endothelial Cells", "Lymphatic Cells",
                                                  "Fibroblasts", "Mast Cells", "Unassigned"))

#t-SNE Plot (Figure E5)
TSNEPlot(object = donor07, do.label = T)

#Find marker genes for clusters (Table E5)
donor07.markers <- FindAllMarkers(object = donor07, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.table(donor07.markers %>% group_by(cluster), file= "/donor07.markers.txt", row.names=FALSE, sep="\t")

#Add Cluster Identity Metadata to Separate Metadata Column
donor07 <- AddMetaData(object = donor07, metadata = donor07@ident, col.name = "indiv.cell.ident")

donor07@meta.data[, "diagnosis"] <- "donor"
donor07@meta.data[, "condition"] <- "donor"

donor07@meta.data$res.0.3 <- NULL
donor07@meta.data$tree.ident <- NULL

#Save object
save(donor07, file = "/donor07.Robj")
load(file = "/donor07.Robj")

# donor08 analysis ---------------------------------------------------------

# Load the donor08 dataset
donor08.data <- Read10X(data.dir = "/donor08_tables/")

donor08 <- CreateSeuratObject(raw.data = donor08.data, min.cells = 3, min.genes = 200, 
                              project = "donor08")

# Identify mitochondrial genes
mito.genes <- grep(pattern = "^MT-", x = rownames(x = donor08@data), value = TRUE)

# Create meta data column containing percentage of mitochondrial genes for each cell
percent.mito <- Matrix::colSums(donor08@raw.data[mito.genes, ])/Matrix::colSums(donor08@raw.data)

# Add percent.mito to object meta data
donor08 <- AddMetaData(object = donor08, metadata = percent.mito, col.name = "percent.mito")

#Data inspection and basic filtering
VlnPlot(object = donor08, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
GenePlot(object = donor08, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = donor08, gene1 = "nUMI", gene2 = "nGene")

donor08 <- FilterCells(object = donor08, subset.names = c("nGene", "percent.mito"), 
                       low.thresholds = c(200, -Inf), high.thresholds = c(6000, 0.2))

#Normalize, find variable genes and scale
donor08 <- NormalizeData(object = donor08, normalization.method = "LogNormalize", 
                         scale.factor = 10000)
donor08 <- FindVariableGenes(object = donor08, mean.function = ExpMean, dispersion.function = LogVMR, 
                             x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = donor08@var.genes)
donor08 <- ScaleData(object = donor08, vars.to.regress = c("nUMI", "percent.mito"))

#Perform PCA
donor08 <- RunPCA(object = donor08, pc.genes = donor08@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5, pcs.compute = 40)
PrintPCA(object = donor08, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(object = donor08, pcs.use = 1:2)
PCAPlot(object = donor08, dim.1 = 1, dim.2 = 2)
donor08 <- ProjectPCA(object = donor08, do.print = FALSE)
PCHeatmap(object = donor08, pc.use = 1:9, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = donor08, pc.use = 10:18, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = donor08, pc.use = 19:27, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCElbowPlot(object = donor08, num.pc = 40)

#Perform clustering using first 12 principal components. 
donor08 <- FindClusters(object = donor08, reduction.type = "pca", dims.use = 1:12, 
                        resolution = 0.5, print.output = 0, save.SNN = TRUE, force.recalc = T)
PrintFindClustersParams(object = donor08)
donor08 <- RunTSNE(object = donor08, dims.use = 1:12, do.fast = TRUE)
TSNEPlot(object = donor08, do.label = T)

#Rename clusters
donor08.ident <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)

donor08.new.ident <- c("AT2 Cells", "Macrophages", "Macrophages", "AT2 Cells",
                       "AT2 Cells", "Macrophages", "Monocytes", "AT1 Cells",
                       "T/NKT Cells", "Endothelial Cells", "Unassigned", "Fibroblasts",
                       "Lymphatic Cells")

donor08@ident <- plyr::mapvalues(x = donor08@ident, from = donor08.ident, to = donor08.new.ident)

donor08@ident <- ordered(donor08@ident,levels = c("Macrophages", "AT2 Cells", "Monocytes", "T/NKT Cells",
                                                  "AT1 Cells", "Endothelial Cells", "Fibroblasts", "Lymphatic Cells",
                                                  "Unassigned"))

#t-SNE Plot (Figure E5)
TSNEPlot(object = donor08, do.label = T)

#Find marker genes for clusters (Table E5)
donor08.markers <- FindAllMarkers(object = donor08, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.table(donor08.markers %>% group_by(cluster), file= "/donor08.markers.txt", row.names=FALSE, sep="\t")

#Add Cluster Identity Metadata to Separate Metadata Column
donor08 <- AddMetaData(object = donor08, metadata = donor08@ident, col.name = "indiv.cell.ident")

donor08@meta.data[, "diagnosis"] <- "donor"
donor08@meta.data[, "condition"] <- "donor"

donor08@meta.data$res.0.3 <- NULL
donor08@meta.data$tree.ident <- NULL

#Save object
save(donor08, file = "/donor08.Robj")
load(file = "/donor08.Robj")


# ipf01 analysis ---------------------------------------------------------

# Load the ipf01 dataset
ipf01.data <- Read10X(data.dir = "/ipf01_tables/")

ipf01 <- CreateSeuratObject(raw.data = ipf01.data, min.cells = 3, min.genes = 200, 
                              project = "ipf01")

# Identify mitochondrial genes
mito.genes <- grep(pattern = "^MT-", x = rownames(x = ipf01@data), value = TRUE)

# Create meta data column containing percentage of mitochondrial genes for each cell
percent.mito <- Matrix::colSums(ipf01@raw.data[mito.genes, ])/Matrix::colSums(ipf01@raw.data)

# Add percent.mito to object meta data
ipf01 <- AddMetaData(object = ipf01, metadata = percent.mito, col.name = "percent.mito")

#Data inspection and basic filtering
VlnPlot(object = ipf01, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
GenePlot(object = ipf01, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = ipf01, gene1 = "nUMI", gene2 = "nGene")

ipf01 <- FilterCells(object = ipf01, subset.names = c("nGene", "percent.mito"), 
                       low.thresholds = c(200, -Inf), high.thresholds = c(6000, 0.2))

#Normalize, find variable genes and scale
ipf01 <- NormalizeData(object = ipf01, normalization.method = "LogNormalize", 
                         scale.factor = 10000)
ipf01 <- FindVariableGenes(object = ipf01, mean.function = ExpMean, dispersion.function = LogVMR, 
                             x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = ipf01@var.genes)
ipf01 <- ScaleData(object = ipf01, vars.to.regress = c("nUMI", "percent.mito"))

#Perform PCA
ipf01 <- RunPCA(object = ipf01, pc.genes = ipf01@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5, pcs.compute = 40)
PrintPCA(object = ipf01, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(object = ipf01, pcs.use = 1:2)
PCAPlot(object = ipf01, dim.1 = 1, dim.2 = 2)
ipf01 <- ProjectPCA(object = ipf01, do.print = FALSE)
PCHeatmap(object = ipf01, pc.use = 1:9, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = ipf01, pc.use = 10:18, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = ipf01, pc.use = 19:27, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCElbowPlot(object = ipf01, num.pc = 40)

#Perform clustering using first 26 principal components. 
ipf01 <- FindClusters(object = ipf01, reduction.type = "pca", dims.use = 1:26, 
                        resolution = 0.5, print.output = 0, save.SNN = TRUE, force.recalc = T)
PrintFindClustersParams(object = ipf01)
ipf01 <- RunTSNE(object = ipf01, dims.use = 1:26, do.fast = TRUE)
TSNEPlot(object = ipf01, do.label = T)

ipf01 <- BuildClusterTree(object = ipf01, do.reorder = TRUE, reorder.numeric = TRUE, do.plot = T)
TSNEPlot(object = ipf01, do.label = T)

#Subset Cluster 13 (Club Cells and Fibroblasts)
ipf01.cluster13 <- SubsetData(ipf01, ident.use = c(13))

ipf01.cluster13 <- ScaleData(object = ipf01.cluster13, vars.to.regress = c("nUMI", "percent.mito"))

ipf01.cluster13 <- FindVariableGenes(object = ipf01.cluster13, mean.function = ExpMean, dispersion.function = LogVMR, 
                                     x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

length(x = ipf01.cluster13@var.genes)

ipf01.cluster13 <- RunPCA(object = ipf01.cluster13, pc.genes = ipf01.cluster13@var.genes, do.print = TRUE, pcs.print = 1:5, 
                          genes.print = 5)

ipf01.cluster13 <- ProjectPCA(object = ipf01.cluster13, do.print = FALSE)

PCHeatmap(object = ipf01.cluster13, pc.use = 1:12, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)

PCElbowPlot(object = ipf01.cluster13, num.pc = 12)

#Clustering of Subsetted Cells
ipf01.cluster13 <- RunTSNE(object = ipf01.cluster13, dims.use = 1:7, do.fast = TRUE)

ipf01.cluster13 <- FindClusters(object = ipf01.cluster13, reduction.type = "pca", dims.use = 1:7, 
                                resolution = 1, save.SNN = TRUE)
TSNEPlot(ipf01.cluster13)

#Markers of Subsetted Cell Subclusters
ipf01.cluster13.markers <- FindAllMarkers(object = ipf01.cluster13, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

ipf01.cluster13.markers %>% group_by(cluster) %>% top_n(3, avg_logFC)

FeaturePlot(object = ipf01.cluster13, features.plot = c("KRT5", "COL3A1", "SCGB3A2"), min.cutoff = "q9",
            cols.use = c("lightgrey", "blue"), pt.size = 1)

#Cell assignment for Subclusters
ipf01.basals <- WhichCells(object = ipf01.cluster13, ident = c(0))
ipf01.clubs <- WhichCells(object = ipf01.cluster13, ident = c(1))
ipf01.fibroblasts <- WhichCells(object = ipf01.cluster13, ident = c(2))

#Rename clusters
ipf01.ident <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)

ipf01.new.ident <- c("Mast Cells", "AT2 Cells", "Plasma Cells", "Monocytes", 
                     "Macrophages", "Macrophages", "Dendritic Cells", "Endothelial Cells",
                     "Monocytes", "Unassigned", "B Cells", "T/NKT Cells",
                     "Basals, Clubs, and Fibroblasts", "Ciliated Cells")

ipf01@ident <- plyr::mapvalues(x = ipf01@ident, from = ipf01.ident, to = ipf01.new.ident)

ipf01 <- SetIdent(object = ipf01, cells.use = ipf01.basals, ident.use = "Basal Cells")
ipf01 <- SetIdent(object = ipf01, cells.use = ipf01.clubs, ident.use = "Club Cells")
ipf01 <- SetIdent(object = ipf01, cells.use = ipf01.fibroblasts, ident.use = "Fibroblasts")

ipf01@ident <- ordered(ipf01@ident,levels = c("Macrophages", "AT2 Cells", "Monocytes", "Club Cells",
                                              "Dendritic Cells", "Plasma Cells", "T/NKT Cells", "Ciliated Cells",
                                              "Endothelial Cells", "B Cells", "Fibroblasts", "Mast Cells",
                                              "Basal Cells", "Unassigned"))

#t-SNE Plot (Figure E5)
TSNEPlot(object = ipf01, do.label = T)

#Find marker genes for clusters (Table E5)
ipf01.markers <- FindAllMarkers(object = ipf01, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.table(ipf01.markers %>% group_by(cluster), file= "/ipf01.markers.txt", row.names=FALSE, sep="\t")

#Macrophage Heterogeneity Feature Plots (Figure E7E)
FeaturePlot(object = ipf01, features.plot = c("CD14", "FCGR3A", "MRC1"), min.cutoff = "q9",
            cols.use = c("lightgrey", "blue"), nCol = 1, pt.size = 0.5)

#Add Cluster Identity Metadata to Separate Metadata Column
ipf01 <- AddMetaData(object = ipf01, metadata = ipf01@ident, col.name = "indiv.cell.ident")

ipf01@meta.data[, "diagnosis"] <- "donor"
ipf01@meta.data[, "condition"] <- "donor"

ipf01@meta.data$res.0.3 <- NULL
ipf01@meta.data$tree.ident <- NULL

#Save object
save(ipf01, file = "/ipf01.Robj")
load(file = "/ipf01.Robj")

# ipf02 analysis ---------------------------------------------------------

# Load the ipf02 dataset
ipf02.data <- Read10X(data.dir = "/ipf02_tables/")

ipf02 <- CreateSeuratObject(raw.data = ipf02.data, min.cells = 3, min.genes = 200, 
                              project = "ipf02")

# Identify mitochondrial genes
mito.genes <- grep(pattern = "^MT-", x = rownames(x = ipf02@data), value = TRUE)

# Create meta data column containing percentage of mitochondrial genes for each cell
percent.mito <- Matrix::colSums(ipf02@raw.data[mito.genes, ])/Matrix::colSums(ipf02@raw.data)

# Add percent.mito to object meta data
ipf02 <- AddMetaData(object = ipf02, metadata = percent.mito, col.name = "percent.mito")

#Data inspection and basic filtering
VlnPlot(object = ipf02, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
GenePlot(object = ipf02, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = ipf02, gene1 = "nUMI", gene2 = "nGene")

ipf02 <- FilterCells(object = ipf02, subset.names = c("nGene", "percent.mito"), 
                       low.thresholds = c(200, -Inf), high.thresholds = c(6000, 0.2))

#Normalize, find variable genes and scale
ipf02 <- NormalizeData(object = ipf02, normalization.method = "LogNormalize", 
                         scale.factor = 10000)
ipf02 <- FindVariableGenes(object = ipf02, mean.function = ExpMean, dispersion.function = LogVMR, 
                             x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = ipf02@var.genes)
ipf02 <- ScaleData(object = ipf02, vars.to.regress = c("nUMI", "percent.mito"))

#Perform PCA
ipf02 <- RunPCA(object = ipf02, pc.genes = ipf02@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5, pcs.compute = 40)
PrintPCA(object = ipf02, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(object = ipf02, pcs.use = 1:2)
PCAPlot(object = ipf02, dim.1 = 1, dim.2 = 2)
ipf02 <- ProjectPCA(object = ipf02, do.print = FALSE)
PCHeatmap(object = ipf02, pc.use = 1:9, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = ipf02, pc.use = 10:18, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = ipf02, pc.use = 19:27, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCElbowPlot(object = ipf02, num.pc = 40)

#Perform clustering using first 26 principal components. 
ipf02 <- FindClusters(object = ipf02, reduction.type = "pca", dims.use = 1:26, 
                        resolution = 0.7, print.output = 0, save.SNN = TRUE, force.recalc = T)
PrintFindClustersParams(object = ipf02)
ipf02 <- RunTSNE(object = ipf02, dims.use = 1:26, do.fast = TRUE)
TSNEPlot(object = ipf02, do.label = T)

ipf02 <- BuildClusterTree(object = ipf02, do.reorder = TRUE, reorder.numeric = TRUE, do.plot = T)

#Rename clusters
ipf02.ident <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)

ipf02.new.ident <- c("Club Cells", "Plasma Cells", "Mast Cells", "Fibroblasts",
                     "Ciliated Cells", "Macrophages", "Unassigned", "Endothelial Cells",
                     "Basal Cells", "T/NKT Cells", "T/NKT Cells", "Dendritic Cells",
                     "B Cells")

ipf02@ident <- plyr::mapvalues(x = ipf02@ident, from = ipf02.ident, to = ipf02.new.ident)

ipf02@ident <- ordered(ipf02@ident,levels = c("Macrophages", "Club Cells", "Dendritic Cells", "Plasma Cells",
                                              "T/NKT Cells", "Ciliated Cells", "Endothelial Cells",
                                              "B Cells", "Fibroblasts", "Mast Cells", "Basal Cells",
                                              "Unassigned"))

#t-SNE Plot (Figure E5)
TSNEPlot(object = ipf02, do.label = T)

#Find marker genes for clusters (Table E5)
ipf02.markers <- FindAllMarkers(object = ipf02, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.table(ipf02.markers %>% group_by(cluster), file= "/ipf02.markers.txt", row.names=FALSE, sep="\t")

#Treg Feature Plots (Figure E6A)
FeaturePlot(object = ipf02, features.plot = c("CD4", "FOXP3", "CD8A", "ITGAE"), min.cutoff = "q9",
            cols.use = c("lightgrey", "blue"), nCol = 4, pt.size = 0.5)

#Monocyte Heterogeneity Feature Plots (Figure E6C)
FeaturePlot(object = ipf02, features.plot = c("CD14", "FCGR3A", "MRC1", "CHI3L1"), min.cutoff = "q9",
            cols.use = c("lightgrey", "blue"), nCol = 1, pt.size = 0.5)

#Macrophage Heterogeneity Feature Plots (Figure E7E)
FeaturePlot(object = ipf02, features.plot = c("MRC1", "FABP4", "SPP1", "CHI3L1"), min.cutoff = "q9",
            cols.use = c("lightgrey", "blue"), nCol = 1, pt.size = 0.5)

#Add Cluster Identity Metadata to Separate Metadata Column
ipf02 <- AddMetaData(object = ipf02, metadata = ipf02@ident, col.name = "indiv.cell.ident")

ipf02@meta.data[, "diagnosis"] <- "donor"
ipf02@meta.data[, "condition"] <- "donor"

ipf02@meta.data$res.0.3 <- NULL
ipf02@meta.data$tree.ident <- NULL

#Save object
save(ipf02, file = "/ipf02.Robj")
load(file = "/ipf02.Robj")

# ipf03 analysis ---------------------------------------------------------

# Load the ipf03 dataset
ipf03.data <- Read10X(data.dir = "/ipf03_tables/")

ipf03 <- CreateSeuratObject(raw.data = ipf03.data, min.cells = 3, min.genes = 200, 
                              project = "ipf03")

# Identify mitochondrial genes
mito.genes <- grep(pattern = "^MT-", x = rownames(x = ipf03@data), value = TRUE)

# Create meta data column containing percentage of mitochondrial genes for each cell
percent.mito <- Matrix::colSums(ipf03@raw.data[mito.genes, ])/Matrix::colSums(ipf03@raw.data)

# Add percent.mito to object meta data
ipf03 <- AddMetaData(object = ipf03, metadata = percent.mito, col.name = "percent.mito")

#Data inspection and basic filtering
VlnPlot(object = ipf03, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
GenePlot(object = ipf03, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = ipf03, gene1 = "nUMI", gene2 = "nGene")

ipf03 <- FilterCells(object = ipf03, subset.names = c("nGene", "percent.mito"), 
                       low.thresholds = c(200, -Inf), high.thresholds = c(5000, 0.2))

#Normalize, find variable genes and scale
ipf03 <- NormalizeData(object = ipf03, normalization.method = "LogNormalize", 
                         scale.factor = 10000)
ipf03 <- FindVariableGenes(object = ipf03, mean.function = ExpMean, dispersion.function = LogVMR, 
                             x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = ipf03@var.genes)
ipf03 <- ScaleData(object = ipf03, vars.to.regress = c("nUMI", "percent.mito"))

#Perform PCA
ipf03 <- RunPCA(object = ipf03, pc.genes = ipf03@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5, pcs.compute = 40)
PrintPCA(object = ipf03, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(object = ipf03, pcs.use = 1:2)
PCAPlot(object = ipf03, dim.1 = 1, dim.2 = 2)
ipf03 <- ProjectPCA(object = ipf03, do.print = FALSE)
PCHeatmap(object = ipf03, pc.use = 1:9, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = ipf03, pc.use = 10:18, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = ipf03, pc.use = 19:27, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCElbowPlot(object = ipf03, num.pc = 40)

#Perform clustering using first 36 principal components. 
ipf03 <- FindClusters(object = ipf03, reduction.type = "pca", dims.use = 1:36, 
                        resolution = 0.5, print.output = 0, save.SNN = TRUE, force.recalc = T)
PrintFindClustersParams(object = ipf03)
ipf03 <- RunTSNE(object = ipf03, dims.use = 1:36, do.fast = TRUE)
TSNEPlot(object = ipf03, do.label = T)

#Rename clusters
ipf03.ident <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)

ipf03.new.ident <- c("Macrophages", "Plasma Cells", "Macrophages", "AT2 Cells", 
                     "Unassigned", "Dendritic Cells", "Macrophages", "Unassigned", 
                     "Ciliated Cells", "Endothelial Cells", "Mast Cells", "Fibroblasts")

ipf03@ident <- plyr::mapvalues(x = ipf03@ident, from = ipf03.ident, to = ipf03.new.ident)

ipf03@ident <- ordered(ipf03@ident,levels = c("Macrophages", "AT2 Cells", "Dendritic Cells", "Plasma Cells",
                                              "Ciliated Cells", "Endothelial Cells", "Fibroblasts", "Mast Cells",
                                              "Unassigned"))

#t-SNE Plot (Figure E5)
TSNEPlot(object = ipf03, do.label = T)

#Find marker genes for clusters (Table E5)
ipf03.markers <- FindAllMarkers(object = ipf03, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.table(ipf03.markers %>% group_by(cluster), file= "/ipf03.markers.txt", row.names=FALSE, sep="\t")

#Add Cluster Identity Metadata to Separate Metadata Column
ipf03 <- AddMetaData(object = ipf03, metadata = ipf03@ident, col.name = "indiv.cell.ident")

ipf03@meta.data[, "diagnosis"] <- "donor"
ipf03@meta.data[, "condition"] <- "donor"

ipf03@meta.data$res.0.3 <- NULL
ipf03@meta.data$tree.ident <- NULL

#Save object
save(ipf03, file = "/ipf03.Robj")
load(file = "/ipf03.Robj")

# ipf04 analysis ---------------------------------------------------------

# Load the ipf04 dataset
ipf04.data <- Read10X(data.dir = "/ipf04_tables/")

ipf04 <- CreateSeuratObject(raw.data = ipf04.data, min.cells = 3, min.genes = 200, 
                              project = "ipf04")

# Identify mitochondrial genes
mito.genes <- grep(pattern = "^MT-", x = rownames(x = ipf04@data), value = TRUE)

# Create meta data column containing percentage of mitochondrial genes for each cell
percent.mito <- Matrix::colSums(ipf04@raw.data[mito.genes, ])/Matrix::colSums(ipf04@raw.data)

# Add percent.mito to object meta data
ipf04 <- AddMetaData(object = ipf04, metadata = percent.mito, col.name = "percent.mito")

#Data inspection and basic filtering
VlnPlot(object = ipf04, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
GenePlot(object = ipf04, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = ipf04, gene1 = "nUMI", gene2 = "nGene")

ipf04 <- FilterCells(object = ipf04, subset.names = c("nGene", "percent.mito"), 
                       low.thresholds = c(200, -Inf), high.thresholds = c(2500, 0.2))

#Normalize, find variable genes and scale
ipf04 <- NormalizeData(object = ipf04, normalization.method = "LogNormalize", 
                         scale.factor = 10000)
ipf04 <- FindVariableGenes(object = ipf04, mean.function = ExpMean, dispersion.function = LogVMR, 
                             x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = ipf04@var.genes)
ipf04 <- ScaleData(object = ipf04, vars.to.regress = c("nUMI", "percent.mito"))

#Perform PCA
ipf04 <- RunPCA(object = ipf04, pc.genes = ipf04@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5, pcs.compute = 40)
PrintPCA(object = ipf04, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(object = ipf04, pcs.use = 1:2)
PCAPlot(object = ipf04, dim.1 = 1, dim.2 = 2)
ipf04 <- ProjectPCA(object = ipf04, do.print = FALSE)
PCHeatmap(object = ipf04, pc.use = 1:9, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = ipf04, pc.use = 10:18, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = ipf04, pc.use = 19:27, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCElbowPlot(object = ipf04, num.pc = 40)

#Perform clustering using first 11 principal components. 
ipf04 <- FindClusters(object = ipf04, reduction.type = "pca", dims.use = 1:11, 
                        resolution = 0.5, print.output = 0, save.SNN = TRUE, force.recalc = T)
PrintFindClustersParams(object = ipf04)
ipf04 <- RunTSNE(object = ipf04, dims.use = 1:11, do.fast = TRUE)
TSNEPlot(object = ipf04, do.label = T)

#Subset Cluster 05 (AT1 and Club Cells)
ipf04.cluster05 <- SubsetData(ipf04, ident.use = c(05))

ipf04.cluster05 <- ScaleData(object = ipf04.cluster05, vars.to.regress = c("nUMI", "percent.mito"))

ipf04.cluster05 <- FindVariableGenes(object = ipf04.cluster05, mean.function = ExpMean, dispersion.function = LogVMR, 
                                     x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

length(x = ipf04.cluster05@var.genes)

ipf04.cluster05 <- RunPCA(object = ipf04.cluster05, pc.genes = ipf04.cluster05@var.genes, do.print = TRUE, pcs.print = 1:5, 
                          genes.print = 5)

ipf04.cluster05 <- ProjectPCA(object = ipf04.cluster05, do.print = FALSE)

PCHeatmap(object = ipf04.cluster05, pc.use = 1:12, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)

PCElbowPlot(object = ipf04.cluster05, num.pc = 12)

#Clustering of Subsetted Cells
ipf04.cluster05 <- RunTSNE(object = ipf04.cluster05, dims.use = 1:6, do.fast = TRUE)

ipf04.cluster05 <- FindClusters(object = ipf04.cluster05, reduction.type = "pca", dims.use = 1:6, 
                                resolution = 0.3, save.SNN = TRUE)
TSNEPlot(ipf04.cluster05)

#Markers of Subsetted Cell Subclusters
ipf04.cluster05.markers <- FindAllMarkers(object = ipf04.cluster05, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

ipf04.cluster05.markers %>% group_by(cluster) %>% top_n(3, avg_logFC)

FeaturePlot(object = ipf04.cluster05, features.plot = c("SCGB3A2", "AGER"), min.cutoff = "q9",
            cols.use = c("lightgrey", "blue"), pt.size = 1)

#Cell assignment for Subclusters
ipf04.at1 <- WhichCells(object = ipf04.cluster05, ident = c(0))
ipf04.clubs <- WhichCells(object = ipf04.cluster05, ident = c(1))

#Rename clusters
ipf04.ident <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)

ipf04.new.ident <- c("Macrophages", "AT2 Cells", "Macrophages", "Macrophages", 
                     "Dendritic Cells", "AT1 and Club Cells", "Ciliated Cells", "Unassigned",
                     "Endothelial Cells", "Plasma Cells")

ipf04@ident <- plyr::mapvalues(x = ipf04@ident, from = ipf04.ident, to = ipf04.new.ident)

ipf04 <- SetIdent(object = ipf04, cells.use = ipf04.at1, ident.use = "AT1 Cells")
ipf04 <- SetIdent(object = ipf04, cells.use = ipf04.clubs, ident.use = "Club Cells")

ipf04@ident <- ordered(ipf04@ident,levels = c("Macrophages", "AT2 Cells", "Club Cells", "Dendritic Cells",
                                              "Plasma Cells", "AT1 Cells", "Ciliated Cells", "Endothelial Cells",
                                              "Unassigned"))

#t-SNE Plot (Figure E5)
TSNEPlot(object = ipf04, do.label = T)

#Find marker genes for clusters (Table E5)
ipf04.markers <- FindAllMarkers(object = ipf04, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.table(ipf04.markers %>% group_by(cluster), file= "/ipf04.markers.txt", row.names=FALSE, sep="\t")

#Add Cluster Identity Metadata to Separate Metadata Column
ipf04 <- AddMetaData(object = ipf04, metadata = ipf04@ident, col.name = "indiv.cell.ident")

ipf04@meta.data[, "diagnosis"] <- "donor"
ipf04@meta.data[, "condition"] <- "donor"

ipf04@meta.data$res.0.3 <- NULL
ipf04@meta.data$tree.ident <- NULL

#Save object
save(ipf04, file = "/ipf04.Robj")
load(file = "/ipf04.Robj")

# hp01 analysis ---------------------------------------------------------

# Load the hp01 dataset
hp01.data <- Read10X(data.dir = "/hp01_tables/")

hp01 <- CreateSeuratObject(raw.data = hp01.data, min.cells = 3, min.genes = 200, 
                              project = "hp01")

# Identify mitochondrial genes
mito.genes <- grep(pattern = "^MT-", x = rownames(x = hp01@data), value = TRUE)

# Create meta data column containing percentage of mitochondrial genes for each cell
percent.mito <- Matrix::colSums(hp01@raw.data[mito.genes, ])/Matrix::colSums(hp01@raw.data)

# Add percent.mito to object meta data
hp01 <- AddMetaData(object = hp01, metadata = percent.mito, col.name = "percent.mito")

#Data inspection and basic filtering
VlnPlot(object = hp01, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
GenePlot(object = hp01, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = hp01, gene1 = "nUMI", gene2 = "nGene")

hp01 <- FilterCells(object = hp01, subset.names = c("nGene", "percent.mito"), 
                       low.thresholds = c(200, -Inf), high.thresholds = c(6000, 0.15))

#Normalize, find variable genes and scale
hp01 <- NormalizeData(object = hp01, normalization.method = "LogNormalize", 
                         scale.factor = 10000)
hp01 <- FindVariableGenes(object = hp01, mean.function = ExpMean, dispersion.function = LogVMR, 
                             x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = hp01@var.genes)
hp01 <- ScaleData(object = hp01, vars.to.regress = c("nUMI", "percent.mito"))

#Perform PCA
hp01 <- RunPCA(object = hp01, pc.genes = hp01@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5, pcs.compute = 40)
PrintPCA(object = hp01, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(object = hp01, pcs.use = 1:2)
PCAPlot(object = hp01, dim.1 = 1, dim.2 = 2)
hp01 <- ProjectPCA(object = hp01, do.print = FALSE)
PCHeatmap(object = hp01, pc.use = 1:9, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = hp01, pc.use = 10:18, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = hp01, pc.use = 19:27, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = hp01, pc.use = 28:36, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCElbowPlot(object = hp01, num.pc = 40)

#Perform clustering using first 32 principal components. 
hp01 <- FindClusters(object = hp01, reduction.type = "pca", dims.use = 1:32, 
                        resolution = 0.7, print.output = 0, save.SNN = TRUE, force.recalc = T)
PrintFindClustersParams(object = hp01)
hp01 <- RunTSNE(object = hp01, dims.use = 1:32, do.fast = TRUE)
TSNEPlot(object = hp01, do.label = T)

#Subset Cluster 02 (Macrophages and DCs)
hp01.cluster02 <- SubsetData(hp01, ident.use = c(02))

hp01.cluster02 <- ScaleData(object = hp01.cluster02, vars.to.regress = c("nUMI", "percent.mito"))

hp01.cluster02 <- FindVariableGenes(object = hp01.cluster02, mean.function = ExpMean, dispersion.function = LogVMR, 
                                    x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

length(x = hp01.cluster02@var.genes)

hp01.cluster02 <- RunPCA(object = hp01.cluster02, pc.genes = hp01.cluster02@var.genes, do.print = TRUE, pcs.print = 1:5, 
                         genes.print = 5)

hp01.cluster02 <- ProjectPCA(object = hp01.cluster02, do.print = FALSE)

PCHeatmap(object = hp01.cluster02, pc.use = 1:12, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)

PCElbowPlot(object = hp01.cluster02, num.pc = 12)

#Clustering of Subsetted Cells
hp01.cluster02 <- RunTSNE(object = hp01.cluster02, dims.use = 1:9, do.fast = TRUE)

hp01.cluster02 <- FindClusters(object = hp01.cluster02, reduction.type = "pca", dims.use = 1:9, 
                               resolution = 0.3, save.SNN = TRUE)
TSNEPlot(hp01.cluster02)

#Markers of Subsetted Cell Subclusters
hp01.cluster02.markers <- FindAllMarkers(object = hp01.cluster02, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

hp01.cluster02.markers %>% group_by(cluster) %>% top_n(3, avg_logFC)

FeaturePlot(object = hp01.cluster02, features.plot = c("CD68", "CLEC10A"), min.cutoff = "q9",
            cols.use = c("lightgrey", "blue"), pt.size = 1)

#Cell assignment for Subclusters
hp01.macrophages <- WhichCells(object = hp01.cluster02, ident = c(0))
hp01.dcs <- WhichCells(object = hp01.cluster02, ident = c(1))

#Rename clusters
hp01.ident <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16)

hp01.new.ident <- c("Macrophages", "AT2 Cells", "Macrophages and DCs", "Monocytes", 
                    "T/NKT Cells", "Macrophages", "Plasma Cells", "Unassigned",
                    "Ciliated Cells", "AT1 Cells", "T Cells", "B Cells",
                    "Endothelial Cells", "DC1", "pDC", "Fibroblasts", "Mast Cells")

hp01@ident <- plyr::mapvalues(x = hp01@ident, from = hp01.ident, to = hp01.new.ident)

hp01 <- SetIdent(object = hp01, cells.use = hp01.macrophages, ident.use = "Macrophages")
hp01 <- SetIdent(object = hp01, cells.use = hp01.dcs, ident.use = "Dendritic Cells")

hp01@ident <- ordered(hp01@ident,levels = c("Macrophages", "AT2 Cells", "Monocytes", "Dendritic Cells",
                                            "Plasma Cells", "T/NKT Cells", "AT1 Cells", "Ciliated Cells", 
                                            "Endothelial Cells", "B Cells", "Fibroblasts", "Mast Cells",
                                            "T Cells", "DC1", "pDC", "Unassigned"))

#t-SNE Plot (Figure E5)
TSNEPlot(object = hp01, do.label = T)

#Find marker genes for clusters (Table E5)
hp01.markers <- FindAllMarkers(object = hp01, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.table(hp01.markers %>% group_by(cluster), file= "/hp01.markers.txt", row.names=FALSE, sep="\t")

#Dendritic Cell Subset Feature Plots (Figure E6B)
FeaturePlot(object = hp01, features.plot = c("CLEC4C", "IRF8", "TLR7", "CLEC9A", "CD1C"), min.cutoff = "q9",
            cols.use = c("lightgrey", "blue"), nCol = 5, pt.size = 0.5)

#Monocyte Heterogeneity Feature Plots (Figure E6C)
FeaturePlot(object = hp01, features.plot = c("CD14", "FCGR3A", "MRC1"), min.cutoff = "q9",
            cols.use = c("lightgrey", "blue"), nCol = 1, pt.size = 0.5)

#Add Cluster Identity Metadata to Separate Metadata Column
hp01 <- AddMetaData(object = hp01, metadata = hp01@ident, col.name = "indiv.cell.ident")

hp01@meta.data[, "diagnosis"] <- "donor"
hp01@meta.data[, "condition"] <- "donor"

hp01@meta.data$res.0.3 <- NULL
hp01@meta.data$tree.ident <- NULL

#Save object
save(hp01, file = "/hp01.Robj")
load(file = "/hp01.Robj")

# ssc_ild01 analysis ---------------------------------------------------------

# Load the ssc_ild01 dataset
ssc_ild01.data <- Read10X(data.dir = "/ssc_ild01_tables/")

ssc_ild01 <- CreateSeuratObject(raw.data = ssc_ild01.data, min.cells = 3, min.genes = 200, 
                              project = "ssc_ild01")

# Identify mitochondrial genes
mito.genes <- grep(pattern = "^MT-", x = rownames(x = ssc_ild01@data), value = TRUE)

# Create meta data column containing percentage of mitochondrial genes for each cell
percent.mito <- Matrix::colSums(ssc_ild01@raw.data[mito.genes, ])/Matrix::colSums(ssc_ild01@raw.data)

# Add percent.mito to object meta data
ssc_ild01 <- AddMetaData(object = ssc_ild01, metadata = percent.mito, col.name = "percent.mito")

#Data inspection and basic filtering
VlnPlot(object = ssc_ild01, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
GenePlot(object = ssc_ild01, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = ssc_ild01, gene1 = "nUMI", gene2 = "nGene")

ssc_ild01 <- FilterCells(object = ssc_ild01, subset.names = c("nGene", "percent.mito"), 
                       low.thresholds = c(200, -Inf), high.thresholds = c(6000, 0.15))

#Normalize, find variable genes and scale
ssc_ild01 <- NormalizeData(object = ssc_ild01, normalization.method = "LogNormalize", 
                         scale.factor = 10000)
ssc_ild01 <- FindVariableGenes(object = ssc_ild01, mean.function = ExpMean, dispersion.function = LogVMR, 
                             x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = ssc_ild01@var.genes)
ssc_ild01 <- ScaleData(object = ssc_ild01, vars.to.regress = c("nUMI", "percent.mito"))

#Perform PCA
ssc_ild01 <- RunPCA(object = ssc_ild01, pc.genes = ssc_ild01@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5, pcs.compute = 40)
PrintPCA(object = ssc_ild01, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(object = ssc_ild01, pcs.use = 1:2)
PCAPlot(object = ssc_ild01, dim.1 = 1, dim.2 = 2)
ssc_ild01 <- ProjectPCA(object = ssc_ild01, do.print = FALSE)
PCHeatmap(object = ssc_ild01, pc.use = 1:9, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = ssc_ild01, pc.use = 10:18, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = ssc_ild01, pc.use = 19:27, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCElbowPlot(object = ssc_ild01, num.pc = 40)

#Perform clustering using first 19 principal components. 
ssc_ild01 <- FindClusters(object = ssc_ild01, reduction.type = "pca", dims.use = 1:19, 
                        resolution = 0.5, print.output = 0, save.SNN = TRUE, force.recalc = T)
PrintFindClustersParams(object = ssc_ild01)
ssc_ild01 <- RunTSNE(object = ssc_ild01, dims.use = 1:19, do.fast = TRUE)
TSNEPlot(object = ssc_ild01, do.label = T)

#Subset Cluster 10 (Lymphocytes)
ssc_ild01.cluster10 <- SubsetData(ssc_ild01, ident.use = c(10))

ssc_ild01.cluster10 <- ScaleData(object = ssc_ild01.cluster10, vars.to.regress = c("nUMI", "percent.mito"))

ssc_ild01.cluster10 <- FindVariableGenes(object = ssc_ild01.cluster10, mean.function = ExpMean, dispersion.function = LogVMR, 
                                         x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

length(x = ssc_ild01.cluster10@var.genes)

ssc_ild01.cluster10 <- RunPCA(object = ssc_ild01.cluster10, pc.genes = ssc_ild01.cluster10@var.genes, do.print = TRUE, pcs.print = 1:5, 
                              genes.print = 5)

ssc_ild01.cluster10 <- ProjectPCA(object = ssc_ild01.cluster10, do.print = FALSE)

PCHeatmap(object = ssc_ild01.cluster10, pc.use = 1:12, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)

PCElbowPlot(object = ssc_ild01.cluster10, num.pc = 12)

#Clustering of Subsetted Cells
ssc_ild01.cluster10 <- RunTSNE(object = ssc_ild01.cluster10, dims.use = 1:9, do.fast = TRUE)

ssc_ild01.cluster10 <- FindClusters(object = ssc_ild01.cluster10, reduction.type = "pca", dims.use = 1:9, 
                                    resolution = 0.5, save.SNN = TRUE)
TSNEPlot(ssc_ild01.cluster10)

#Markers of Subsetted Cell Subclusters
ssc_ild01.cluster10.markers <- FindAllMarkers(object = ssc_ild01.cluster10, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

ssc_ild01.cluster10.markers %>% group_by(cluster) %>% top_n(3, avg_logFC)

FeaturePlot(object = ssc_ild01.cluster10, features.plot = c("IGKC", "MRC1"), min.cutoff = "q9",
            cols.use = c("lightgrey", "blue"), pt.size = 1)

#Cell assignment for Subclusters
ssc_ild01.plasmas <- WhichCells(object = ssc_ild01.cluster10, ident = c(0))
ssc_ild01.macrophages <- WhichCells(object = ssc_ild01.cluster10, ident = c(1))

#Rename clusters
ssc_ild01.ident <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)

ssc_ild01.new.ident <- c("Macrophages", "AT2 Cells", "Macrophages", "Macrophages", 
                         "Macrophages", "Monocytes", "Ciliated Cells", "Club Cells",
                         "Dendritic Cells", "Unassigned", "Plasma Cells and Macrophages",
                         "Endothelial Cells", "AT1 Cells")

ssc_ild01@ident <- plyr::mapvalues(x = ssc_ild01@ident, from = ssc_ild01.ident, to = ssc_ild01.new.ident)

ssc_ild01 <- SetIdent(object = ssc_ild01, cells.use = ssc_ild01.plasmas, ident.use = "Plasma Cells")
ssc_ild01 <- SetIdent(object = ssc_ild01, cells.use = ssc_ild01.macrophages, ident.use = "Macrophages")

ssc_ild01@ident <- ordered(ssc_ild01@ident,levels = c("Macrophages", "AT2 Cells", "Monocytes", "Club Cells",
                                                      "Dendritic Cells", "Plasma Cells", "AT1 Cells", "Ciliated Cells",
                                                      "Endothelial Cells", "Unassigned"))

#t-SNE Plot (Figure E5)
TSNEPlot(object = ssc_ild01, do.label = T)

#Find marker genes for clusters (Table E5)
ssc_ild01.markers <- FindAllMarkers(object = ssc_ild01, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.table(ssc_ild01.markers %>% group_by(cluster), file= "/ssc_ild01.markers.txt", row.names=FALSE, sep="\t")

#Monocyte Heterogeneity Feature Plots (Figure E6C)
FeaturePlot(object = ssc_ild01, features.plot = c("CD14", "FCGR3A", "MRC1"), min.cutoff = "q9",
            cols.use = c("lightgrey", "blue"), nCol = 1, pt.size = 0.5)

#Macrophage Heterogeneity Feature Plots (Figure E7E)
FeaturePlot(object = ssc_ild01, features.plot = c("MRC1", "FABP4", "SPP1", "CHI3L1"), min.cutoff = "q9",
            cols.use = c("lightgrey", "blue"), nCol = 1, pt.size = 0.5)

#Add Cluster Identity Metadata to Separate Metadata Column
ssc_ild01 <- AddMetaData(object = ssc_ild01, metadata = ssc_ild01@ident, col.name = "indiv.cell.ident")

ssc_ild01@meta.data[, "diagnosis"] <- "donor"
ssc_ild01@meta.data[, "condition"] <- "donor"

ssc_ild01@meta.data$res.0.3 <- NULL
ssc_ild01@meta.data$tree.ident <- NULL

#Save object
save(ssc_ild01, file = "/ssc_ild01.Robj")
load(file = "/ssc_ild01.Robj")

# myositis_ild01 analysis ---------------------------------------------------------

# Load the myositis_ild01 dataset
myositis_ild01.data <- Read10X(data.dir = "/myositis_ild01_tables/")

myositis_ild01 <- CreateSeuratObject(raw.data = myositis_ild01.data, min.cells = 3, min.genes = 200, 
                              project = "myositis_ild01")

# Identify mitochondrial genes
mito.genes <- grep(pattern = "^MT-", x = rownames(x = myositis_ild01@data), value = TRUE)

# Create meta data column containing percentage of mitochondrial genes for each cell
percent.mito <- Matrix::colSums(myositis_ild01@raw.data[mito.genes, ])/Matrix::colSums(myositis_ild01@raw.data)

# Add percent.mito to object meta data
myositis_ild01 <- AddMetaData(object = myositis_ild01, metadata = percent.mito, col.name = "percent.mito")

#Data inspection and basic filtering
VlnPlot(object = myositis_ild01, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
GenePlot(object = myositis_ild01, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = myositis_ild01, gene1 = "nUMI", gene2 = "nGene")

myositis_ild01 <- FilterCells(object = myositis_ild01, subset.names = c("nGene", "percent.mito"), 
                       low.thresholds = c(200, -Inf), high.thresholds = c(4000, 0.3))

#Normalize, find variable genes and scale
myositis_ild01 <- NormalizeData(object = myositis_ild01, normalization.method = "LogNormalize", 
                         scale.factor = 10000)
myositis_ild01 <- FindVariableGenes(object = myositis_ild01, mean.function = ExpMean, dispersion.function = LogVMR, 
                             x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = myositis_ild01@var.genes)
myositis_ild01 <- ScaleData(object = myositis_ild01, vars.to.regress = c("nUMI", "percent.mito"))

#Perform PCA
myositis_ild01 <- RunPCA(object = myositis_ild01, pc.genes = myositis_ild01@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5, pcs.compute = 40)
PrintPCA(object = myositis_ild01, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(object = myositis_ild01, pcs.use = 1:2)
PCAPlot(object = myositis_ild01, dim.1 = 1, dim.2 = 2)
myositis_ild01 <- ProjectPCA(object = myositis_ild01, do.print = FALSE)
PCHeatmap(object = myositis_ild01, pc.use = 1:9, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = myositis_ild01, pc.use = 10:18, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = myositis_ild01, pc.use = 19:27, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCElbowPlot(object = myositis_ild01, num.pc = 40)

#Perform clustering using first 17 principal components. 
myositis_ild01 <- FindClusters(object = myositis_ild01, reduction.type = "pca", dims.use = 1:17, 
                        resolution = 0.7, print.output = 0, save.SNN = TRUE, force.recalc = T)
PrintFindClustersParams(object = myositis_ild01)
myositis_ild01 <- RunTSNE(object = myositis_ild01, dims.use = 1:17, do.fast = TRUE)
TSNEPlot(object = myositis_ild01, do.label = T)

#Rename clusters
myositis_ild01.ident <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)

myositis_ild01.new.ident <- c("Macrophages", "Macrophages", "Ciliated Cells", "AT2 Cells",
                              "Club Cells", "Macrophages", "Dendritic Cells", "KRT5 Positive Epithelial Cells",
                              "Macrophages", "AT2 Cells", "Mast Cells", "Unassigned",
                              "T/NKT Cells", "Monocytes", "Fibroblasts", "Endothelial Cells")

myositis_ild01@ident <- plyr::mapvalues(x = myositis_ild01@ident, from = myositis_ild01.ident, to = myositis_ild01.new.ident)

myositis_ild01@ident <- ordered(myositis_ild01@ident,levels = c("Macrophages", "AT2 Cells", "Monocytes", "Club Cells",
                                                                "Dendritic Cells", "T/NKT Cells", "Ciliated Cells", "Endothelial Cells",
                                                                "Fibroblasts", "Mast Cells", "KRT5 Positive Epithelial Cells", "Unassigned"))

#t-SNE Plot (Figure E5)
TSNEPlot(object = myositis_ild01, do.label = T)

#Find marker genes for clusters (Table E5)
myositis_ild01.markers <- FindAllMarkers(object = myositis_ild01, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.table(myositis_ild01.markers %>% group_by(cluster), file= "/myositis_ild01.markers.txt", row.names=FALSE, sep="\t")

#Add Cluster Identity Metadata to Separate Metadata Column
myositis_ild01 <- AddMetaData(object = myositis_ild01, metadata = myositis_ild01@ident, col.name = "indiv.cell.ident")

myositis_ild01@meta.data[, "diagnosis"] <- "donor"
myositis_ild01@meta.data[, "condition"] <- "donor"

myositis_ild01@meta.data$res.0.3 <- NULL
myositis_ild01@meta.data$tree.ident <- NULL

#Save object
save(myositis_ild01, file = "/myositis_ild01.Robj")
load(file = "/myositis_ild01.Robj")

# ssc_ild02 analysis ---------------------------------------------------------

# Load the ssc_ild02 dataset
ssc_ild02.data <- Read10X(data.dir = "/ssc_ild02_tables/")

ssc_ild02 <- CreateSeuratObject(raw.data = ssc_ild02.data, min.cells = 3, min.genes = 200, 
                              project = "ssc_ild02")

# Identify mitochondrial genes
mito.genes <- grep(pattern = "^MT-", x = rownames(x = ssc_ild02@data), value = TRUE)

# Create meta data column containing percentage of mitochondrial genes for each cell
percent.mito <- Matrix::colSums(ssc_ild02@raw.data[mito.genes, ])/Matrix::colSums(ssc_ild02@raw.data)

# Add percent.mito to object meta data
ssc_ild02 <- AddMetaData(object = ssc_ild02, metadata = percent.mito, col.name = "percent.mito")

#Data inspection and basic filtering
VlnPlot(object = ssc_ild02, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
GenePlot(object = ssc_ild02, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = ssc_ild02, gene1 = "nUMI", gene2 = "nGene")

ssc_ild02 <- FilterCells(object = ssc_ild02, subset.names = c("nGene", "percent.mito"), 
                       low.thresholds = c(200, -Inf), high.thresholds = c(6000, 0.3))

#Normalize, find variable genes and scale
ssc_ild02 <- NormalizeData(object = ssc_ild02, normalization.method = "LogNormalize", 
                         scale.factor = 10000)
ssc_ild02 <- FindVariableGenes(object = ssc_ild02, mean.function = ExpMean, dispersion.function = LogVMR, 
                             x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = ssc_ild02@var.genes)
ssc_ild02 <- ScaleData(object = ssc_ild02, vars.to.regress = c("nUMI", "percent.mito"))

#Perform PCA
ssc_ild02 <- RunPCA(object = ssc_ild02, pc.genes = ssc_ild02@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5, pcs.compute = 40)
PrintPCA(object = ssc_ild02, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(object = ssc_ild02, pcs.use = 1:2)
PCAPlot(object = ssc_ild02, dim.1 = 1, dim.2 = 2)
ssc_ild02 <- ProjectPCA(object = ssc_ild02, do.print = FALSE)
PCHeatmap(object = ssc_ild02, pc.use = 1:9, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = ssc_ild02, pc.use = 10:18, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = ssc_ild02, pc.use = 19:27, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCElbowPlot(object = ssc_ild02, num.pc = 40)

#Perform clustering using first 17 principal components. 
ssc_ild02 <- FindClusters(object = ssc_ild02, reduction.type = "pca", dims.use = 1:17, 
                        resolution = 0.5, print.output = 0, save.SNN = TRUE, force.recalc = T)
PrintFindClustersParams(object = ssc_ild02)
ssc_ild02 <- RunTSNE(object = ssc_ild02, dims.use = 1:17, do.fast = TRUE)
TSNEPlot(object = ssc_ild02, do.label = T)

#Rename clusters
ssc_ild02.ident <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)

ssc_ild02.new.ident <- c("AT2 Cells", "Macrophages", "Macrophages", "AT2 Cells",
                         "Macrophages", "Monocytes", "T/NKT Cells", "Unassigned",
                         "AT1 Cells", "Endothelial Cells", "Ciliated Cells", "Plasma Cells")

ssc_ild02@ident <- plyr::mapvalues(x = ssc_ild02@ident, from = ssc_ild02.ident, to = ssc_ild02.new.ident)

ssc_ild02@ident <- ordered(ssc_ild02@ident,levels = c("Macrophages", "AT2 Cells", "Monocytes", "Plasma Cells",
                                                      "T/NKT Cells", "AT1 Cells", "Ciliated Cells", "Endothelial Cells",
                                                      "Unassigned"))
#t-SNE Plot (Figure E5)
TSNEPlot(object = ssc_ild02, do.label = T)

#Find marker genes for clusters (Table E5)
ssc_ild02.markers <- FindAllMarkers(object = ssc_ild02, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.table(ssc_ild02.markers %>% group_by(cluster), file= "/ssc_ild02.markers.txt", row.names=FALSE, sep="\t")

#Add Cluster Identity Metadata to Separate Metadata Column
ssc_ild02 <- AddMetaData(object = ssc_ild02, metadata = ssc_ild02@ident, col.name = "indiv.cell.ident")

ssc_ild02@meta.data[, "diagnosis"] <- "donor"
ssc_ild02@meta.data[, "condition"] <- "donor"

ssc_ild02@meta.data$res.0.3 <- NULL
ssc_ild02@meta.data$tree.ident <- NULL

#Save object
save(ssc_ild02, file = "/ssc_ild02.Robj")
load(file = "/ssc_ild02.Robj")

# cryo analysis ---------------------------------------------------------

# Load the cryo dataset
cryo.data <- Read10X(data.dir = "/cryo_tables/")

cryo <- CreateSeuratObject(raw.data = cryo.data, min.cells = 3, min.genes = 200, 
                              project = "cryo")

# Identify mitochondrial genes
mito.genes <- grep(pattern = "^MT-", x = rownames(x = cryo@data), value = TRUE)

# Create meta data column containing percentage of mitochondrial genes for each cell
percent.mito <- Matrix::colSums(cryo@raw.data[mito.genes, ])/Matrix::colSums(cryo@raw.data)

# Add percent.mito to object meta data
cryo <- AddMetaData(object = cryo, metadata = percent.mito, col.name = "percent.mito")

#Data inspection and basic filtering
VlnPlot(object = cryo, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
GenePlot(object = cryo, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = cryo, gene1 = "nUMI", gene2 = "nGene")

cryo <- FilterCells(object = cryo, subset.names = c("nGene", "percent.mito"), 
                       low.thresholds = c(200, -Inf), high.thresholds = c(5000, 0.2))

#Normalize, find variable genes and scale
cryo <- NormalizeData(object = cryo, normalization.method = "LogNormalize", 
                         scale.factor = 10000)
cryo <- FindVariableGenes(object = cryo, mean.function = ExpMean, dispersion.function = LogVMR, 
                             x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = cryo@var.genes)
cryo <- ScaleData(object = cryo, vars.to.regress = c("nUMI", "percent.mito"))

#Perform PCA
cryo <- RunPCA(object = cryo, pc.genes = cryo@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5, pcs.compute = 40)
PrintPCA(object = cryo, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(object = cryo, pcs.use = 1:2)
PCAPlot(object = cryo, dim.1 = 1, dim.2 = 2)
cryo <- ProjectPCA(object = cryo, do.print = FALSE)
PCHeatmap(object = cryo, pc.use = 1:9, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = cryo, pc.use = 10:18, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = cryo, pc.use = 19:27, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCElbowPlot(object = cryo, num.pc = 40)

#Perform clustering using first 23 principal components. 
cryo <- FindClusters(object = cryo, reduction.type = "pca", dims.use = 1:23, 
                        resolution = 0.5, print.output = 0, save.SNN = TRUE, force.recalc = T)
PrintFindClustersParams(object = cryo)
cryo <- RunTSNE(object = cryo, dims.use = 1:23, do.fast = TRUE)
TSNEPlot(object = cryo, do.label = T)

#Subset Cluster 05 (Club and Airway Cells)
cryo.cluster05 <- SubsetData(cryo, ident.use = c(5))

cryo.cluster05 <- ScaleData(object = cryo.cluster05, vars.to.regress = c("nUMI", "percent.mito"))

cryo.cluster05 <- FindVariableGenes(object = cryo.cluster05, mean.function = ExpMean, dispersion.function = LogVMR, 
                                    x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

length(x = cryo.cluster05@var.genes)

cryo.cluster05 <- RunPCA(object = cryo.cluster05, pc.genes = cryo.cluster05@var.genes, do.print = TRUE, pcs.print = 1:5, 
                         genes.print = 5)

cryo.cluster05 <- ProjectPCA(object = cryo.cluster05, do.print = FALSE)

PCHeatmap(object = cryo.cluster05, pc.use = 1:12, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)

PCElbowPlot(object = cryo.cluster05, num.pc = 12)

#Clustering of Subsetted Cells
cryo.cluster05 <- RunTSNE(object = cryo.cluster05, dims.use = 1:10, do.fast = TRUE)

cryo.cluster05 <- FindClusters(object = cryo.cluster05, reduction.type = "pca", dims.use = 1:10, 
                               resolution = 0.5, save.SNN = TRUE)
TSNEPlot(cryo.cluster05)

#Markers of Subsetted Cell Subclusters
cryo.cluster05.markers <- FindAllMarkers(object = cryo.cluster05, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

cryo.cluster05.markers %>% group_by(cluster) %>% top_n(3, avg_logFC)

FeaturePlot(object = cryo.cluster05, features.plot = c("SCGB3A2", "TPPP3"), min.cutoff = "q9",
            cols.use = c("lightgrey", "blue"), pt.size = 1)

#Cell assignment for Subclusters
cryo.clubs <- WhichCells(object = cryo.cluster05, ident = c(0))
cryo.ciliated <- WhichCells(object = cryo.cluster05, ident = c(1))

#Rename clusters
cryo.ident <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)

cryo.new.ident <- c("T/NKT Cells", "Macrophages", "AT2 Cells", "Dendritic Cells", 
                    "T/NKT Cells", "Club and Ciliated Cells", "Endothelial Cells", "Lymphatic Cells",
                    "Fibroblasts", "Mast Cells", "AT1 Cells", "Endothelial Cells",
                    "Lymphatic Progenitors", "Smooth Muscle Cells")

cryo@ident <- plyr::mapvalues(x = cryo@ident, from = cryo.ident, to = cryo.new.ident)

cryo <- SetIdent(object = cryo, cells.use = cryo.clubs, ident.use = "Club Cells")
cryo <- SetIdent(object = cryo, cells.use = cryo.ciliated, ident.use = "Ciliated Cells")

cryo@ident <- ordered(cryo@ident,levels = c("Macrophages", "AT2 Cells", "Club Cells", "Dendritic Cells",
                                            "T/NKT Cells", "AT1 Cells", "Ciliated Cells",
                                            "Endothelial Cells", "Lymphatic Cells", "Lymphatic Progenitors",
                                            "Fibroblasts", "Mast Cells", "Smooth Muscle Cells"))

#t-SNE Plot (Figure 10C)
TSNEPlot(object = cryo, do.label = T)

#Find marker genes for clusters (Table E5)
cryo.markers <- FindAllMarkers(object = cryo, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.table(cryo.markers %>% group_by(cluster), file= "/cryo.markers.txt", row.names=FALSE, sep="\t")

#Endothelial Cell Subset Feature Plots (Figure 10D)
FeaturePlot(object = cryo, features.plot = c("VWF", "SOX17", "PTGS1", "PROX1"), min.cutoff = "q9",
            cols.use = c("lightgrey", "blue"), nCol = 2, pt.size = 0.5)

#Add Cluster Identity Metadata to Separate Metadata Column
cryo <- AddMetaData(object = cryo, metadata = cryo@ident, col.name = "indiv.cell.ident")

cryo@meta.data[, "diagnosis"] <- "donor"
cryo@meta.data[, "condition"] <- "donor"

cryo@meta.data$res.0.3 <- NULL
cryo@meta.data$tree.ident <- NULL

#Save object
save(cryo, file = "/cryo.Robj")
load(file = "/cryo.Robj")

# Combining 16 individual samples-------------------------------

load(file = paste(datadir, "/donor01.Robj",sep=""))
load(file = paste(datadir, "/donor02.Robj",sep=""))
load(file = paste(datadir, "/donor03.Robj",sep=""))
load(file = paste(datadir, "/donor04.Robj",sep=""))
load(file = paste(datadir, "/donor05.Robj",sep=""))
load(file = paste(datadir, "/donor06.Robj",sep=""))
load(file = paste(datadir, "/donor07.Robj",sep=""))
load(file = paste(datadir, "/donor08.Robj",sep=""))

load(file = paste(datadir, "/ipf01.Robj",sep=""))
load(file = paste(datadir, "/ipf02.Robj",sep=""))
load(file = paste(datadir, "/ipf03.Robj",sep=""))
load(file = paste(datadir, "/ipf04.Robj",sep=""))
load(file = paste(datadir, "/hp01.Robj",sep=""))
load(file = paste(datadir, "/ssc_ild01.Robj",sep=""))
load(file = paste(datadir, "/myositis_ild01.Robj",sep=""))
load(file = paste(datadir, "/ssc_ild02.Robj",sep=""))

donor01 <- RenameCells(donor01, add.cell.id = "donor01")
donor02 <- RenameCells(donor02, add.cell.id = "donor02")
donor03 <- RenameCells(donor03, add.cell.id = "donor03")
donor04 <- RenameCells(donor04, add.cell.id = "donor04")
donor05 <- RenameCells(donor05, add.cell.id = "donor05")
donor06 <- RenameCells(donor06, add.cell.id = "donor06")
donor07 <- RenameCells(donor07, add.cell.id = "donor07")
donor08 <- RenameCells(donor08, add.cell.id = "donor08")

ipf01 <- RenameCells(ipf01, add.cell.id = "ipf01")
ipf02 <- RenameCells(ipf02, add.cell.id = "ipf02")
ipf03 <- RenameCells(ipf03, add.cell.id = "ipf03")
ipf04 <- RenameCells(ipf04, add.cell.id = "ipf04")
hp01 <- RenameCells(hp01, add.cell.id = "hp01")
ssc_ild01 <- RenameCells(ssc_ild01, add.cell.id = "ssc_ild01")
myositis_ild01 <- RenameCells(myositis_ild01, add.cell.id = "myositis_ild01")
ssc_ild02 <- RenameCells(ssc_ild02, add.cell.id = "ssc_ild02")

#Merge each individual donor into the combined object
hl.all <- MergeSeurat(donor01, donor02, project = "all")

hl.all <- MergeSeurat(hl.all, donor03)
hl.all <- MergeSeurat(hl.all, donor04)
hl.all <- MergeSeurat(hl.all, donor05)
hl.all <- MergeSeurat(hl.all, donor06)
hl.all <- MergeSeurat(hl.all, donor07)
hl.all <- MergeSeurat(hl.all, donor08)

hl.all <- MergeSeurat(hl.all, ipf01)
hl.all <- MergeSeurat(hl.all, ipf02)
hl.all <- MergeSeurat(hl.all, ipf03)
hl.all <- MergeSeurat(hl.all, ipf04)
hl.all <- MergeSeurat(hl.all, hp01)
hl.all <- MergeSeurat(hl.all, ssc_ild01)
hl.all <- MergeSeurat(hl.all, myositis_ild01)
hl.all <- MergeSeurat(hl.all, ssc_ild02)

hl.all <- AddMetaData(object = hl.all, metadata = hl.all@ident, col.name = "subject.id")

#Save object
save(hl.all, file = "/hl.all.Robj")
load(file = "/hl.all.Robj")

# Merging Donors for CCA Object ---------------------------------------------------

#Load donor objects

load(file = paste(datadir, "/donor01.Robj",sep=""))
load(file = paste(datadir, "/donor02.Robj",sep=""))
load(file = paste(datadir, "/donor03.Robj",sep=""))
load(file = paste(datadir, "/donor04.Robj",sep=""))
load(file = paste(datadir, "/donor05.Robj",sep=""))
load(file = paste(datadir, "/donor06.Robj",sep=""))
load(file = paste(datadir, "/donor07.Robj",sep=""))
load(file = paste(datadir, "/donor08.Robj",sep=""))

#Rename donor cells

donor01 <- RenameCells(donor01, add.cell.id = "donor01")
donor02 <- RenameCells(donor02, add.cell.id = "donor02")
donor03 <- RenameCells(donor03, add.cell.id = "donor03")
donor04 <- RenameCells(donor04, add.cell.id = "donor04")
donor05 <- RenameCells(donor05, add.cell.id = "donor05")
donor06 <- RenameCells(donor06, add.cell.id = "donor06")
donor07 <- RenameCells(donor07, add.cell.id = "donor07")
donor08 <- RenameCells(donor08, add.cell.id = "donor08")


#Merge each individual donor into the combined object
donor <- MergeSeurat(donor01, donor02, project = "donor")

donor <- MergeSeurat(donor, donor03)
donor <- MergeSeurat(donor, donor04)
donor <- MergeSeurat(donor, donor05)
donor <- MergeSeurat(donor, donor06)
donor <- MergeSeurat(donor, donor07)
donor <- MergeSeurat(donor, donor08)

save(donor, "/donor.hl.Robj")
load("/donor.hl.Robj", verbose=TRUE)

# Merging fibrosis samples ---------------------------------------------------

#Load fibrosis objects

load(file = paste(datadir, "/ipf01.Robj",sep=""))
load(file = paste(datadir, "/ipf02.Robj",sep=""))
load(file = paste(datadir, "/ipf03.Robj",sep=""))
load(file = paste(datadir, "/ipf04.Robj",sep=""))
load(file = paste(datadir, "/hp01.Robj",sep=""))
load(file = paste(datadir, "/ssc_ild01.Robj",sep=""))
load(file = paste(datadir, "/myositis_ild01.Robj",sep=""))
load(file = paste(datadir, "/ssc_ild02.Robj",sep=""))

#Rename fibrosis cells

ipf01 <- RenameCells(ipf01, add.cell.id = "ipf01")
ipf02 <- RenameCells(ipf02, add.cell.id = "ipf02")
ipf03 <- RenameCells(ipf03, add.cell.id = "ipf03")
ipf04 <- RenameCells(ipf04, add.cell.id = "ipf04")
hp01 <- RenameCells(hp01, add.cell.id = "hp01")
ssc_ild01 <- RenameCells(ssc_ild01, add.cell.id = "ssc_ild01")
myositis_ild01 <- RenameCells(myositis_ild01, add.cell.id = "myositis_ild01")
ssc_ild02 <- RenameCells(ssc_ild02, add.cell.id = "ssc_ild02")


fibrosis <- MergeSeurat(ipf01, ipf02, project = "fibrosis")

#fibrosis <- MergeSeurat(fibrosis, sc13, add.cell.id2 = "sc13")
fibrosis <- MergeSeurat(fibrosis, ipf03)
fibrosis <- MergeSeurat(fibrosis, ipf04)
fibrosis <- MergeSeurat(fibrosis, hp01)
fibrosis <- MergeSeurat(fibrosis, ssc_ild01)
fibrosis <- MergeSeurat(fibrosis, myositis_ild01)
fibrosis <- MergeSeurat(fibrosis, ssc_ild02)

save(fibrosis, "/fibrosis.hl.Robj")
load("/fibrosis.hl.Robj", verbose=TRUE)

# Prepare donor and fibrosis data for alignment ------------------------------------------------

# take the union of the top 1500 variable genes in each dataset for alignment
hvg.donor <- rownames(x = head(x = donor@hvg.info, n = 1500))
hvg.fibrosis <- rownames(x = head(x = fibrosis@hvg.info, n = 1500))
hvg.union <- union(x = hvg.donor, y = hvg.fibrosis)

# set the 'condition' in each dataset for easy identification
donor@meta.data[, "condition"] <- "donor"
fibrosis@meta.data[, "condition"] <- "fibrosis"

hl <- RunCCA(object = donor, object2 = fibrosis, genes.use = hvg.union, num.cc = 40)

# visualize results of CCA plot CC1 versus CC2 and look at a violin plot
p1 <- DimPlot(object = hl, reduction.use = "cca", group.by = "condition", pt.size = 0.5, 
              do.return = TRUE)
p2 <- VlnPlot(object = hl, features.plot = "CC1", group.by = "condition", do.return = TRUE)
plot_grid(p1, p2)
dev.off()

p3 <- MetageneBicorPlot(hl, grouping.var = "condition", dims.eval = 1:40, 
                        display.progress = FALSE)

PrintDim(object = hl, reduction.type = "cca", dims.print = 1:2, genes.print = 10)

DimHeatmap(object = hl, reduction.type = "cca", cells.use = 500, dim.use = 1:9, 
           do.balanced = TRUE)
DimHeatmap(object = hl, reduction.type = "cca", cells.use = 500, dim.use = 10:18, 
           do.balanced = TRUE)
DimHeatmap(object = hl, reduction.type = "cca", cells.use = 500, dim.use = 19:27, 
           do.balanced = TRUE)
DimHeatmap(object = hl, reduction.type = "cca", cells.use = 500, dim.use = 28:36, 
           do.balanced = TRUE)
DimHeatmap(object = hl, reduction.type = "cca", cells.use = 500, dim.use = 37:40, 
           do.balanced = TRUE)

# Alignment ------------------------------------------------

hl <- AlignSubspace(object = hl, reduction.type = "cca", grouping.var = "condition", 
                      dims.align = 1:31)

# violin plot for cc1 and cc2
p1 <- VlnPlot(object = hl, features.plot = "ACC1", group.by = "condition", 
              do.return = TRUE)
p2 <- VlnPlot(object = hl, features.plot = "ACC2", group.by = "condition", 
              do.return = TRUE)
plot_grid(p1, p2)
dev.off()

# Integrated Seurat analysis of HL ------------------------------------------------

hl <- FindClusters(object = hl, reduction.type = "cca.aligned", dims.use = 1:31, resolution = 0.5)
hl <- RunTSNE(object = hl, reduction.use = "cca.aligned", dims.use = 1:31, 
              do.fast = TRUE)


TSNEPlot(object = hl, do.label = TRUE, do.return = TRUE, pt.size = 0.01)

# Subset out Clusters 10 and 19  ------------------------------------------------
hl <- SubsetData(hl, ident.use = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 20, 21, 22))

hl <- FindClusters(object = hl, reduction.type = "cca.aligned", dims.use = 1:31, resolution = 0.2, force.recalc = TRUE)
hl <- RunTSNE(object = hl, reduction.use = "cca.aligned", dims.use = 1:31, 
              do.fast = TRUE)

# Explore cluster tree and merge nodes ------------------------------------------------

hl <- BuildClusterTree(object = hl, do.reorder = TRUE, reorder.numeric = TRUE)

node.scores <- AssessNodes(object = hl)
node.scores[order(node.scores$oobe, decreasing = TRUE), ] -> node.scores

write.table(node.scores, file="/subset.cluster.tree.node.scores.txt", row.names=FALSE, sep="\t")

# choose the nodes to merge
selected.nodes <- c(27, 42, 29)
nodes.merge <- node.scores[node.scores$node %in% selected.nodes,]
nodes.to.merge <- sort(x = nodes.merge$node)

for (n in nodes.to.merge) {
  hl <- MergeNode(object = hl, node.use = n)
}

TSNEPlot(object = hl, do.label = TRUE, do.return = TRUE, pt.size = 0.01)

# Rename and reorder clusters------------------------------

current.ident <- c(1, 2, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23)

new.ident <- c("Plasma Cells", "Macrophages", "AT2 Cells", "Club Cells", "AT2 Cells",
               "B Cells", "T/NKT Cells", "Monocytes", "Dendritic Cells", "Basal Cells",
               "Mast Cells", "T/NKT Cells", "Endothelial/Lymphatic Cells", "Fibroblasts", "AT1 Cells",
               "Ciliated Cells")

hl@ident <- plyr::mapvalues(x = hl.merged@ident, from = current.ident, to = new.ident)

hl@ident <- ordered(hl@ident,levels = c("AT1 Cells", "AT2 Cells", "Club Cells", "Ciliated Cells", "Basal Cells",
                                        "Macrophages", "Monocytes", "Dendritic Cells", "T/NKT Cells",
                                        "B Cells", "Plasma Cells", "Mast Cells",
                                        "Endothelial/Lymphatic Cells", "Fibroblasts"))

current.org.ident <- unique(hl@meta.data$orig.ident)

paper.id <- c("Donor 1", "Donor 2", "Donor 3", "Donor 4", "Donor 5", "Donor 6", "Donor 7", "Donor 8",
              "IPF 1", "IPF 2", "IPF 3", "IPF 4", "HP", "SSc-ILD 1", "PM-ILD", "SSc-ILD 2")

hl@ident <- plyr::mapvalues(x = hl@ident, from = current.org.ident, to = paper.id)

# Cluster markers ------------------------------------------------
hl.markers <- FindAllMarkers(object = hl)

#Cluster marker table (Table E1)
write.table(hl.markers %>% group_by(cluster), file="/cluster.markers.txt", row.names=FALSE, sep="\t")

# tSNE analysis ------------------------------------------------

#t-SNE plot with labeled clusters (Figure 1A)
TSNEPlot(object = hl, do.label = TRUE, do.return = TRUE, pt.size = 0.01)

#t-SNE plot with labeled condition (Figure 1B)
TSNEPlot(object = hl, group.by = "condition")

#t-SNE plot with labeled individual (Figure 1C)
TSNEPlot(object = hl, group.by = "orig.ident")

#Cluster marker feature plots (Figure 1D)
FeaturePlot(object = hl, features.plot = c("AGER", "SFTPC", "SCGB3A2", "TPPP3", "KRT5",
                                           "CD68", "FCN1", "CLEC10A", "CD3D", "IGHG4",
                                           "MS4A1", "TPSB2", "VWF", "DCN"), min.cutoff = "q9",
            cols.use = c("lightgrey", "blue"), nCol = 5, pt.size = 0.5)

# Top marker heatmap ------------------------------------------------

hl.markers <- read.table("/cluster.markers.txt", header=T, sep="\t")

heatmap.order <- factor(c("Macrophages", "AT2 Cells", "Club Cells", "Monocytes", "Dendritic Cells",
                          "T/NKT Cells", "Ciliated Cells", "Plasma Cells", "AT1 Cells", "Endothelial/Lymphatic Cells",
                          "Basal Cells", "Mast Cells", "Fibroblasts", "B Cells"), 
                        levels = (c("Macrophages", "AT2 Cells", "Club Cells", "Monocytes", "Dendritic Cells",
                                    "T/NKT Cells", "Ciliated Cells", "Plasma Cells", "AT1 Cells", "Endothelial/Lymphatic Cells",
                                    "Basal Cells", "Mast Cells", "Fibroblasts", "B Cells")))

hl.markers$cluster <- factor(hl.markers$cluster, 
                                 levels = (c("Macrophages", "AT2 Cells", "Club Cells", "Monocytes", "Dendritic Cells",
                                             "T/NKT Cells", "Ciliated Cells", "Plasma Cells", "AT1 Cells", "Endothelial/Lymphatic Cells",
                                             "Basal Cells", "Mast Cells", "Fibroblasts", "B Cells")))

hl.markers <- hl.markers[order(hl.markers$cluster),]

top5.markers <- hl.markers %>% group_by(cluster) %>% top_n(5, avg_logFC)

#Top marker heatmap (Figure E2A)
DoHeatmap(hl, genes.use = top5.merged.markers$gene, 
          group.order = heatmap.order, slim.col.label = T, group.label.rot = T, group.spacing=0.3, 
          group.cex=26, cex.row=16,
          remove.key = T)  + 
  theme(axis.text.x=element_text(angle=45,hjust=1))

# Differential Expression Analysis for Macrophages, AT2 cells, and Fibroblasts ------------------------------------------------

#Output Background List
gene.background.list <- hl@raw.data@Dimnames[[1]]
gene.background.list = paste0('"', gene.background.list, '"')
write.table(gene.background.list, file="/gene.background.txt", row.names=FALSE, sep="\t")

# Macrophage Differential Expression--------------------------------------

Macrophages <- SubsetData(hl, ident.use = c("Macrophages"))
Macrophages <- SetAllIdent(Macrophages, id = "condition")

Macrophage.DE <- FindMarkers(object = Macrophages, ident.1 = "fibrosis", logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1)

Macrophage.DE$adj.p <- p.adjust(Macrophage.DE$p_val, method = "fdr", n = length(Macrophage.DE$p_val))

up <- Macrophage.DE$adj.p < 0.05 & Macrophage.DE$avg_logFC > 0
flag <- rep(0, nrow(Macrophage.DE))
flag[up] <- 1
Macrophage.DE$up <- flag

down <- Macrophage.DE$adj.p < 0.05 & Macrophage.DE$avg_logFC < 0
flag <- rep(0, nrow(Macrophage.DE))
flag[down] <- 1
Macrophage.DE$down <- flag

Macrophage.DE <- Macrophage.DE[order(-Macrophage.DE$up, -Macrophage.DE$down, -Macrophage.DE$avg_logFC),]

Macrophage.DE <- cbind(gene = row.names(Macrophage.DE), Macrophage.DE)

Macrophage.DE$gene <- paste0("'", Macrophage.DE$gene, "'")

Macrophage.DE <- rbind(Macrophage.DE[(Macrophage.DE$up == 1),], Macrophage.DE[(Macrophage.DE$down == 1),])

write.table(Macrophage.DE, file="/Macrophage_DE.txt", row.names=FALSE, sep="\t")

#Macrophage DE Heatmap (Figure 2A)
DoHeatmap(Macrophages, genes.use = Macrophage.DE$gene, group.by = "condition", 
          slim.col.label = T, remove.key = T,
          cex.row=32, group.cex=14,
          title = "Macrophages")

save(Macrophages, file = "/macrophages.Robj")
load("/macrophages.Robj", verbose=TRUE)

# AT2 Cell Differential Expression--------------------------------------

AT2s <- SubsetData(hl, ident.use = c("AT2 Cells"))
AT2s <- SetAllIdent(AT2s, id = "condition")

#Unadjusted Counts DE

AT2.DE <- FindMarkers(object = AT2s, ident.1 = "fibrosis", logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1)

AT2.DE$adj.p <- p.adjust(AT2.DE$p_val, method = "fdr", n = length(AT2.DE$p_val))

up <- AT2.DE$adj.p < 0.05 & AT2.DE$avg_logFC > 0
flag <- rep(0, nrow(AT2.DE))
flag[up] <- 1
AT2.DE$up <- flag

down <- AT2.DE$adj.p < 0.05 & AT2.DE$avg_logFC < 0
flag <- rep(0, nrow(AT2.DE))
flag[down] <- 1
AT2.DE$down <- flag

AT2.DE <- AT2.DE[order(-AT2.DE$up, -AT2.DE$down, -AT2.DE$avg_logFC),]

AT2.DE <- cbind(gene = row.names(AT2.DE), AT2.DE)

AT2.DE$gene <- paste0("'", AT2.DE$gene, "'")

AT2.DE <- rbind(AT2.DE[(AT2.DE$up == 1),], AT2.DE[(AT2.DE$down == 1),])

write.table(AT2.DE, file="/AT2_DE.txt", row.names=FALSE, sep="\t")

#AT2 DE heatmap (Figure 2B)
DoHeatmap(AT2s, genes.use = AT2.DE$gene, group.by = "condition", 
          slim.col.label = T, remove.key = T,
          cex.row=10, group.cex=14,
          title = "AT2s")

save(AT2s, file = "/at2.Robj")
load("/at2.Robj", verbose=TRUE)

# Fibroblast Differential Expression--------------------------------------

Fibroblasts <- SubsetData(hl, ident.use = c("Fibroblasts"))
Fibroblasts <- SetAllIdent(Fibroblasts, id = "condition")

#Unadjusted Counts DE

Fibroblast.DE <- FindMarkers(object = Fibroblasts, ident.1 = "fibrosis", logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1)

Fibroblast.DE$adj.p <- p.adjust(Fibroblast.DE$p_val, method = "fdr", n = length(Fibroblast.DE$p_val))

up <- Fibroblast.DE$adj.p < 0.05 & Fibroblast.DE$avg_logFC > 0
flag <- rep(0, nrow(Fibroblast.DE))
flag[up] <- 1
Fibroblast.DE$up <- flag

down <- Fibroblast.DE$adj.p < 0.05 & Fibroblast.DE$avg_logFC < 0
flag <- rep(0, nrow(Fibroblast.DE))
flag[down] <- 1
Fibroblast.DE$down <- flag

Fibroblast.DE <- Fibroblast.DE[order(-Fibroblast.DE$up, -Fibroblast.DE$down, -Fibroblast.DE$avg_logFC),]

Fibroblast.DE <- cbind(gene = row.names(Fibroblast.DE), Fibroblast.DE)

Fibroblast.DE$gene <- paste0("'", Fibroblast.DE$gene, "'")

Fibroblast.DE <- rbind(Fibroblast.DE[(Fibroblast.DE$up == 1),], Fibroblast.DE[(Fibroblast.DE$down == 1),])

write.table(Fibroblast.DE, file="/Fibroblast_DE.txt", row.names=FALSE, sep="\t")

#Fibroblast DE Heatmap (Figure 2C)
DoHeatmap(Fibroblasts, genes.use = Fibroblast.DE$gene, group.by = "condition", 
          slim.col.label = T, remove.key = T,
          cex.row=32, group.cex=14,
          title = "Fibroblasts")

save(Fibroblasts, file = "/fibroblasts.Robj")
load("/fibroblasts.Robj", verbose=TRUE)

# Donor versus fibrosis Violin Plots------------------------------------------------

#AM DE Violin Plots (Figure 2J)
VlnPlot(object = Macrophages, features.plot = c("CHI3L1"), x.lab.rot=T,
        do.return = T, cols.use = c("#00BFC4", "#F8766D"))+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())
VlnPlot(object = Macrophages, features.plot = c("FN1"), x.lab.rot=T, point.size.use=0.001,
        do.return = T, cols.use = c("#00BFC4", "#F8766D"))+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())
VlnPlot(object = Macrophages, features.plot = c("SPP1"), x.lab.rot=T, point.size.use=0.001,
        do.return = T, cols.use = c("#00BFC4", "#F8766D"))+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())
VlnPlot(object = Macrophages, features.plot = c("APOE"), x.lab.rot=T, point.size.use=0.001,
        do.return = T, cols.use = c("#00BFC4", "#F8766D"))+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())

#AT2 DE Violin Plots (Figure 2K)
VlnPlot(object = AT2s, features.plot = c("FN1"), x.lab.rot=T, point.size.use=0.001,
        do.return = T, cols.use = c("#00BFC4", "#F8766D"))+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())
VlnPlot(object = AT2s, features.plot = c("CHI3L1"), x.lab.rot=T, point.size.use=0.001,
        do.return = T, cols.use = c("#00BFC4", "#F8766D"))+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())
VlnPlot(object = AT2s, features.plot = c("HIF1A"), x.lab.rot=T, point.size.use=0.001,
        do.return = T, cols.use = c("#00BFC4", "#F8766D"))+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())
VlnPlot(object = AT2s, features.plot = c("MMP7"), x.lab.rot=T, point.size.use=0.001,
        do.return = T, cols.use = c("#00BFC4", "#F8766D"))+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())

#Fibroblast DE Violin Plots (Figure 2L)
VlnPlot(object = Fibroblasts, features.plot = c("COL1A1"), x.lab.rot=T, point.size.use=0.001,
        do.return = T, cols.use = c("#00BFC4", "#F8766D"))+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())
VlnPlot(object = Fibroblasts, features.plot = c("ACTA2"), x.lab.rot=T, point.size.use=0.001,
        do.return = T, cols.use = c("#00BFC4", "#F8766D"))+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())
VlnPlot(object = Fibroblasts, features.plot = c("POSTN"), x.lab.rot=T, point.size.use=0.001,
        do.return = T, cols.use = c("#00BFC4", "#F8766D"))+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())
VlnPlot(object = Fibroblasts, features.plot = c("TIMP1"), x.lab.rot=T, point.size.use=0.001,
        do.return = T, cols.use = c("#00BFC4", "#F8766D"))+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())

# Single-cell versus bulk-cell DE----------------------------

AT2.DE.all <- FindMarkers(object = AT2s, ident.1 = "fibrosis", logfc.threshold = -Inf, test.use = "wilcox", min.pct = -Inf)

AT2.bulk.DE <- read.table(file = "/at2_donor_fibrosis_DESeq.txt",
  header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")

AT2.DE.all$bulkFC <- AT2.bulk.DE$log2FoldChange[match(AT2.DE.all$gene, AT2.bulk.DE$gene)]

AT2.DE.all <- AT2.DE.all[!(is.na(AT2.DE.all$bulkFC)),]

at2_callout_genes <- c("SERPINA1", "CHI3L1", "NAMPT", "FN1", "HLA-DQA1", "HIF1A", "HES1", "FGG", "DMBT1")

AT2.DE.all <- AT2.DE.all%>%mutate(
  threshold = ifelse(gene %in% at2_callout_genes, "A", ifelse((abs(avg_logFC) < 0.05) | (abs(bulkFC) <0.05), "B", "C")))

#AT2 bulk versus single-cell DE (Figure E4A)
ggplot(AT2.DE.all, aes(x=avg_logFC, y=bulkFC))  +
  geom_point(aes(colour = threshold), size=1.5) +
  scale_colour_manual(values = c("A" = "red", "B"= "grey", "C"="black")) +
  geom_vline(xintercept=0, linetype="dashed", color = "red") +
  geom_hline(yintercept=0, linetype="dashed", color = "red") + 
  geom_abline(intercept = 0, slope = 1, color = "green", size = 1) +
  labs(title="AT2 Cell DE", x ="Single-cell logFC", y = "Bulk Cell logFC") +
  geom_label_repel(data=subset(AT2.DE.all, gene %in% at2_callout_genes),
                   aes(x=avg_logFC, y=bulkFC,label=gene),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'red') +
  theme(legend.position="none")

AT2.DE.all.subset <- subset(AT2.DE.all, abs(avg_logFC) >= 0.05 & abs(bulkFC) >= 0.05)

cor.test(AT2.DE.all.subset$avg_logFC, AT2.DE.all.subset$bulkFC)

Macrophage.DE.all <- FindMarkers(object = Macrophages, ident.1 = "fibrosis", logfc.threshold = -Inf, test.use = "wilcox", min.pct = -Inf)

Macrophage.bulk.DE <- read.table(file = "\am_donor_fibrosis_DESeq.txt",
  header = TRUE, sep = "\t", stringsAsFactors = FALSE)

Macrophage.DE.all$bulkFC <- Macrophage.bulk.DE$log2FoldChange[match(Macrophage.DE.all$gene, Macrophage.bulk.DE$gene)]

Macrophage.DE.all <- Macrophage.DE.all[!(is.na(Macrophage.DE.all$bulkFC)),]

macrophage_callout_genes <- c("CHI3L1", "MARCKS", "IL1RN", "PLA2G7", "MMP7", "SPP1")

Macrophage.DE.all <- Macrophage.DE.all%>%mutate(
  threshold = ifelse(gene %in% macrophage_callout_genes, "A", ifelse((abs(avg_logFC) < 0.05) | (abs(bulkFC) <0.05), "B", "C")))

#AT2 bulk versus single-cell DE (Figure E4B)
ggplot(Macrophage.DE.all, aes(x=avg_logFC, y=bulkFC)) +
  geom_point(aes(colour = threshold), size=1.5) +
  scale_colour_manual(values = c("A" = "red", "B"= "grey", "C"="black")) +
  geom_vline(xintercept=0, linetype="dashed", color = "red") +
  geom_hline(yintercept=0, linetype="dashed", color = "red") + 
  geom_abline(intercept = 0, slope = 1, color = "green", size = 1) +
  labs(title="Macrophage Cell DE", x ="Single-cell logFC", y = "Bulk Cell logFC") +
  geom_label_repel(data=subset(Macrophage.DE.all, gene %in% macrophage_callout_genes),
                   aes(x=avg_logFC, y=bulkFC,label=gene),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'red') +
  theme(legend.position="none")

Macrophage.DE.all.subset <- subset(Macrophage.DE.all, abs(avg_logFC) >= 0.05 & abs(bulkFC) >= 0.05)

cor.test(Macrophage.DE.all.subset$avg_logFC, Macrophage.DE.all.subset$bulkFC)

# Contribution of different cell populations to pulmonary fibrosis pathways ------------------------------------------------

cell.idents <- c("AT2 Cells", "Basal Cells", "Endothelial/Lymphatic Cells", "Ciliated Cells", "AT1 Cells",
                 "Club Cells", "Fibroblasts", "Macrophages")

dotplot.hl <- SubsetData(hl, ident.use = cell.idents)

dotplot.hl@ident <- factor(dotplot.hl@ident, levels = (c("Macrophages", "Basal Cells", "Club Cells", "Ciliated Cells", 
                                                         "AT2 Cells", "AT1 Cells", "Endothelial/Lymphatic Cells", "Fibroblasts")))

dotplot.genes <- c("WNT2", "WNT3A", "WNT7A", "WNT7B", "WNT5A",
                   "WNT9A", "RSPO1", "RSPO2", "RSPO3",
                   "LGR4", "LGR6", "LRP5", "LRP6", "AXIN2", "ZNRF3", "RNF43", 
                   "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", "HES1", "HEY1", "JAG1", "JAG2",
                   "SHH", "GLI1", "GLI2", "GLI3", 
                   "YAP1", "TAZ",
                   "SNAI1", "SNAI2", "SNAI3",
                   "TWIST1", "TWIST2",
                   "ZEB1", "ZEB2",
                   "HIF1A",
                   "PDGFRA", "PDGFRB",
                   "KRT5", "KRT14",
                   "ITGA6", "ITGB4",
                   "CDH1")

genegroups <- c(rep(c("EMT-Related Genes"), times = 21), rep(c("Notch Pathway Genes"), times = 8), rep(c("WNT Pathway Genes"), times = 16))

genegroups <- factor(genegroups, levels = (c("WNT Pathway Genes", "Notch Pathway Genes", "EMT-Related Genes")))

#Dot Plots (Figure 4A)
SplitDotPlotGG(dotplot.hl, genes.plot = rev(dotplot.genes), gene.groups = rev(genegroups), cols.use = c("blue", "red"), 
               x.lab.rot = T, plot.legend = T, dot.scale = 8, do.return = T, grouping.var = "condition")

#Violin Plots (Figure 4B)

dotplot.hl@ident <- factor(dotplot.hl@ident, levels = rev(c("Macrophages", "Basal Cells", "Club Cells", "Ciliated Cells", 
                                                            "AT2 Cells", "AT1 Cells", "Endothelial/Lymphatic Cells", "Fibroblasts")))

VlnPlot(object = dotplot.hl, features.plot = c("WNT2"), do.return = T)+
  theme(axis.title.x=element_blank(), 
        axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1),
        plot.title = element_text(face="italic"))
VlnPlot(object = dotplot.hl, features.plot = c("WNT7B"), do.return = T)+
  theme(axis.title.x=element_blank(), 
        axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1),
        plot.title = element_text(face="italic"))
VlnPlot(object = dotplot.hl, features.plot = c("RSPO3"), do.return = T)+
  theme(axis.title.x=element_blank(), 
        axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1),
        plot.title = element_text(face="italic"))
VlnPlot(object = dotplot.hl, features.plot = c("AXIN2"), do.return = T)+
  theme(axis.title.x=element_blank(), 
        axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1),
        plot.title = element_text(face="italic"))

# Subset Macrophages---------------------------------------------------

All.Macrophages <- SubsetData(hl.all, ident.use = c("Macrophages"))
All.Macrophages <- ScaleData(object = All.Macrophages, vars.to.regress = c("nUMI", "percent.mito", "orig.ident"))
All.Macrophages <- FindVariableGenes(object = All.Macrophages, mean.function = ExpMean, dispersion.function = LogVMR, 
                                 x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = All.Macrophages@var.genes)

All.Macrophages <- RunPCA(object = All.Macrophages, pc.genes = All.Macrophages@var.genes, do.print = TRUE, pcs.print = 1:5, 
                      genes.print = 5, pcs.compute = 40)
All.Macrophages <- ProjectPCA(object = All.Macrophages, do.print = FALSE)
PCHeatmap(object = All.Macrophages, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = All.Macrophages, pc.use = 13:24, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = All.Macrophages, pc.use = 25:36, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCElbowPlot(object = All.Macrophages, num.pc = 40)

All.Macrophages <- RunTSNE(object = All.Macrophages, dims.use = 1:27, do.fast = TRUE)
All.Macrophages <- FindClusters(object = All.Macrophages, reduction.type = "pca", dims.use = 1:27, 
                            resolution = 0.1, save.SNN = TRUE, force.recalc = T)

All.Macrophages <- SubsetData(All.Macrophages, ident.use = c(0, 1, 2, 3, 5))
All.Macrophages <- ScaleData(object = All.Macrophages, vars.to.regress = c("nUMI", "percent.mito", "orig.ident"))
All.Macrophages <- FindVariableGenes(object = All.Macrophages, mean.function = ExpMean, dispersion.function = LogVMR, 
                                 x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = All.Macrophages@var.genes)

All.Macrophages <- RunPCA(object = All.Macrophages, pc.genes = All.Macrophages@var.genes, do.print = TRUE, pcs.print = 1:5, 
                      genes.print = 5, pcs.compute = 40)
All.Macrophages <- ProjectPCA(object = All.Macrophages, do.print = FALSE)
PCHeatmap(object = All.Macrophages, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = All.Macrophages, pc.use = 13:24, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = All.Macrophages, pc.use = 25:36, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCElbowPlot(object = All.Macrophages, num.pc = 40)

All.Macrophages <- RunTSNE(object = All.Macrophages, dims.use = 1:23, do.fast = TRUE)
All.Macrophages <- FindClusters(object = All.Macrophages, reduction.type = "pca", dims.use = 1:23, 
                            resolution = 0.06, save.SNN = TRUE)

#t-SNE of all macrophages by cluster (Figure 5A)
TSNEPlot(object = All.Macrophages, pt.size = 0.5)

#t-SNE of all macrophages by condition (Figure 5B)
TSNEPlot(object = All.Macrophages, group.by = "condition",
         colors.use = c("#00BFC4", "#F8766D"),
         pt.size = 0.5)

#t-sne of all macrophages by diagnosis (Figure E7A)
TSNEPlot(object = All.Macrophages, group.by = "diagnosis", do.return = TRUE, pt.size = 0.5,
         colors.use = c("#4292c6", "#fb6a4a", "#99000d", "#ef3b2c", "#fcbba1"))

#All macrophage feature plots (Figure 5D)
FeaturePlot(object = All.Macrophages, features.plot = c("PPARG", "APOE", "MARCO",
                                                    "MRC1", "MAFB", "SIGLEC1",
                                                    "IL1RN", "MMP9", "CHI3L1",
                                                    "SPP1", "MARCKS", "PLA2G7"), min.cutoff = "q9",
            cols.use = c("lightgrey", "blue"), pt.size = 0.5, nCol = 3)


macrophage_freq_table <- data.frame(table(All.Macrophages@ident, All.Macrophages@meta.data[, "condition"]), 2)

#Bar plot by condition (Figure 5C)
ggplot(data=macrophage_freq_table, aes(x=Var1, y=Freq, fill = Subject)) +
  geom_bar(stat="identity", color="black", position = 'fill') + 
  labs(title = "Cell Distribution per Condition\nby Macrophage Cluster", x="Cluster", y="Proportion") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

macrophage_freq_df_orig <- data.frame(table(All.Macrophages@ident, All.Macrophages@meta.data[, "orig.ident"]), 2)
colnames(macrophage_freq_df_orig)[colnames(macrophage_freq_df_orig) == 'Var2'] <- 'Subject'

#Bar plot by individaul (Figure E7B)
ggplot(data=macrophage_freq_df_orig, aes(x=Var1, y=Freq, fill = Subject)) +
  geom_bar(stat="identity", color="black", position = 'fill') + 
  labs(title = "Cell Distribution per Subject\nby Macrophage Cluster", x="Cluster", y="Proportion") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

All.Macrophages <- SetAllIdent(object = All.Macrophages, id = 'orig.ident')

#Violin plots of CHI3L1, MMP9, and SPP1 (Figure 5F)

VlnPlot(object = All.Macrophages, features.plot = c("CHI3L1"), x.lab.rot = T, y.max = 5, point.size.use = 0.1,
        ident.include = c("donor01", "donor02", "donor03", "donor04",
                          "donor05", "donor06", "donor07", "donor08"))
VlnPlot(object = All.Macrophages, features.plot = c("CHI3L1"), x.lab.rot = T, y.max = 5, point.size.use = 0.1,
        ident.include = c("ipf01", "ipf02", "ipf03", "ipf04", "hp01", "ssc_ild01", "myositis_ild01", "ssc_ild02"))

VlnPlot(object = All.Macrophages, features.plot = c("MMP9"), x.lab.rot = T, y.max = 5, point.size.use = 0.1,
        ident.include = c("donor01", "donor02", "donor03", "donor04",
                          "donor05", "donor06", "donor07", "donor08"))
VlnPlot(object = All.Macrophages, features.plot = c("MMP9"), x.lab.rot = T, y.max = 5, point.size.use = 0.1,
        ident.include = c("ipf01", "ipf02", "ipf03", "ipf04", "hp01", "ssc_ild01", "myositis_ild01", "ssc_ild02"))

VlnPlot(object = All.Macrophages, features.plot = c("SPP1"), x.lab.rot = T, y.max = 5, point.size.use = 0.1,
        ident.include = c("donor01", "donor02", "donor03", "donor04",
                          "donor05", "donor06", "donor07", "donor08"))
VlnPlot(object = All.Macrophages, features.plot = c("SPP1"), x.lab.rot = T, y.max = 5, point.size.use = 0.1,
        ident.include = c("ipf01", "ipf02", "ipf03", "ipf04", "hp01", "ssc_ild01", "myositis_ild01", "ssc_ild02"))

#Find marker genes for clusters
Macrophages.markers <- FindAllMarkers(object = All.Macrophages, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

#Macrophage subcluster markers (Table E6)
write.table(Macrophages.markers %>% group_by(cluster), file="?macrophage_markers.txt", row.names=FALSE, sep="\t")

save(Macrophages, file = "/all.macrophages.Robj")
load(file = "/all.macrophages.Robj")

# Subset Epithelials---------------------------------------------------

hl <- SetAllIdent(hl, id="indiv.cell.ident")

Epithelials <- SubsetData(hl, ident.use = c("AT2 Cells", "AT1 Cells", "Ciliated Cells", "Club Cells",
                                            "Basal Cells", "KRT5 Positive Epithelial Cells"))

Epithelials <- ScaleData(object = Epithelials, vars.to.regress = c("nUMI", "percent.mito", "orig.ident"))
Epithelials <- FindVariableGenes(object = Epithelials, mean.function = ExpMean, dispersion.function = LogVMR, 
                                 x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = Epithelials@var.genes)

Epithelials <- RunPCA(object = Epithelials, pc.genes = Epithelials@var.genes, pcs.compute = 40)
Epithelials <- ProjectPCA(object = Epithelials, do.print = FALSE)
PCHeatmap(object = Epithelials, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = Epithelials, pc.use = 13:24, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = Epithelials, pc.use = 25:36, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCElbowPlot(object = Epithelials, num.pc = 40)

Epithelials <- RunTSNE(object = Epithelials, dims.use = 1:24, do.fast = TRUE)
Epithelials <- FindClusters(object = Epithelials, reduction.type = "pca", dims.use = 1:24, 
                            resolution = 0.1, save.SNN = TRUE)

Epithelials <- SubsetData(Epithelials, ident.use = c(0, 1, 2, 3, 4, 5))
Epithelials <- ScaleData(object = Epithelials, vars.to.regress = c("nUMI", "percent.mito", "orig.ident"))
Epithelials <- FindVariableGenes(object = Epithelials, mean.function = ExpMean, dispersion.function = LogVMR, 
                                 x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = Epithelials@var.genes)

Epithelials <- RunPCA(object = Epithelials, pc.genes = Epithelials@var.genes, pcs.compute = 40)
Epithelials <- ProjectPCA(object = Epithelials, do.print = FALSE)
PCHeatmap(object = Epithelials, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = Epithelials, pc.use = 13:24, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = Epithelials, pc.use = 25:36, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCElbowPlot(object = Epithelials, num.pc = 40)

Epithelials <- RunTSNE(object = Epithelials, dims.use = 1:20, do.fast = TRUE)
Epithelials <- FindClusters(object = Epithelials, reduction.type = "pca", dims.use = 1:20, 
                            resolution = 0.07, save.SNN = TRUE)

#t-SNE plot of reclustered epithelial cells by cluster (Figure 6A)
TSNEPlot(object = Epithelials, pt.size = 0.5)

#t-SNE plot of reclustered epithelial cells by condition (Figure 6B)
TSNEPlot(object = Epithelials, group.by = "condition",
         colors.use = c( "#EE82EE", "#00BFC4", "#F8766D"),
         pt.size = 0.5)

#t-SNE plot of reclustered epithelial cells by diagnosis (Figure E7C)
TSNEPlot(object = Epithelials, group.by = "diagnosis", do.return = TRUE, pt.size = 0.5,
         colors.use = c("#4292c6", "#fb6a4a", "#99000d", "#ef3b2c", "#fcbba1"))

epithelial_freq_table <- data.frame(table(Epithelials@ident, Epithelials@meta.data[, "condition"]), 2)

#bar plot of cluster frequencies by condition (Figure 6C)
ggplot(data=epithelial_freq_table, aes(x=Var1, y=Freq, fill = Var2)) +
  geom_bar(stat="identity", color="black", position = 'fill', show.legend = FALSE) + 
  labs(x="Cluster", y="Proportion") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = c("#00BFC4", "#F8766D"))

epithelial_freq_df_orig <- data.frame(table(Epithelials@ident, Epithelials@meta.data[, "orig.ident"]), 2)
colnames(epithelial_freq_df_orig)[colnames(epithelial_freq_df_orig) == 'Var2'] <- 'Subject'

#bar plot of cluster frequencies by individual (Figure E7D)
ggplot(data=epithelial_freq_df_orig, aes(x=Var1, y=Freq, fill = Subject)) +
  geom_bar(stat="identity", color="black", position = 'fill') + 
  labs(title = "Cell Distribution per Subject\nby Epithelial Cluster", x="Cluster", y="Proportion") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#Epithelial cell feature plots (Figure 6D)
FeaturePlot(object = Epithelials, features.plot = c("SFTPC", "AGER", "SCGB1A1",
                                                    "MUC5B", "FOXJ1", "RFX2", 
                                                    "HIF1A", "CHI3L1", "NKX2-1",
                                                    "HHIP", "FASN", "HES1"), min.cutoff = "q9",
            cols.use = c("lightgrey", "blue"), nCol = 4, pt.size = 0.5)


#Epithelial cell stem cell features (Figure 9A)
FeaturePlot(object = Epithelials, features.plot = c("KRT5", "TP63", "SOX2",
                                                    "ITGA6", "NGFR", "MYC"), min.cutoff = "q9",
            cols.use = c("lightgrey", "blue"), nCol = 3, pt.size = 0.5)

Epithelials <- SetAllIdent(object = Epithelials, id = 'orig.ident')

#Violin plots of SERPINA1 and CHI3L1 (Figure 6E)

VlnPlot(object = Epithelials, features.plot = c("SERPINA1"), x.lab.rot = T, y.max = 5, point.size.use = 0.1,
        ident.include = c("donor01", "donor02", "donor03", "donor04",
                          "donor05", "donor06", "donor07", "donor08"))
VlnPlot(object = Epithelials, features.plot = c("SERPINA1"), x.lab.rot = T, y.max = 5, point.size.use = 0.1,
        ident.include = c("ipf01", "ipf02", "ipf03", "ipf04", "hp01", "ssc_ild01", "myositis_ild01", "ssc_ild02"))

VlnPlot(object = Epithelials, features.plot = c("CHI3L1"), x.lab.rot = T, y.max = 5, point.size.use = 0.1,
        ident.include = c("donor01", "donor02", "donor03", "donor04",
                          "donor05", "donor06", "donor07", "donor08"))
VlnPlot(object = Epithelials, features.plot = c("CHI3L1"), x.lab.rot = T, y.max = 5, point.size.use = 0.1,
        ident.include = c("ipf01", "ipf02", "ipf03", "ipf04", "hp01", "ssc_ild01", "myositis_ild01", "ssc_ild02"))


#Epithelial cell Wnt pathway violin plots (Figure E8E)
VlnPlot(object = Epithelials.wnt, features.plot = c("WNT5A", "WNT7A", "WNT7B", "AXIN2"),
        group.by = "orig.ident", remove.legend = T)

sen_table <- read.table(file="senescence_table.txt", header=F, sep="\t")

sen_genes <- list(levels(sen_table[1,]))

Epithelials <- AddModuleScore(Epithelials, genes.list = sen_genes, ctrl.size = 5, enrich.name = "Senescence")

#Epithelial cell senescence score feature plot (Figure 9B)
FeaturePlot(Epithelials, features.plot = "Senescence1", cols.use = c("lightgrey", "blue"), min.cutoff = 'q7')

#Epithelial cell senescence score ridge plot by condition (Figure 9C)
RidgePlot(Epithelials, features.plot = "Senescence1", group.by = "condition", cols.use = c("#00BFC4", "#F8766D"), do.return = TRUE)

#Epithelial cell senescence score ridge plot by subcluster (Figure 9D)
RidgePlot(Epithelials, features.plot = "Senescence1", do.return = TRUE)

#Find marker genes for clusters
Epithelials.markers <- FindAllMarkers(object = Epithelials, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

#Epithelial subcluster markers (Table E8)
write.table(Epithelials.markers %>% group_by(cluster), file=paste(datadir, "tables/subcluster/20180718_epithelial_markers.txt",sep=""), row.names=FALSE, sep="\t")
Epithelials.average <- Epithelials

Epithelials.average@meta.data$Epithelials <- rep("Epithelials",nrow(Epithelials.average@meta.data))

Epithelials.average <- SetAllIdent(object = Epithelials.average, id = 'Epithelials')

epithelial_background <- AverageExpression(Epithelials.average)
epithelial_background$gene <- row.names(epithelial_background)
names(epithelial_background) <- c("expression", "gene")
epithelial_background <- epithelial_background[,c("gene", "expression")]
epithelial_background$gene <- paste0("'", epithelial_background$gene, "',")
epithelial_background <- epithelial_background[order(-epithelial_background$expression),]

write.table(epithelial_background, file=paste(datadir, "tables/subcluster/20180330_epithelial_background.txt",sep=""), row.names=FALSE, sep="\t")

write.table(Epithelials@meta.data, file=paste(datadir, "tables/subcluster/20180330_epithelials_sensecence_meta.txt",sep=""), row.names=FALSE, sep="\t")

#Epithelial cell senescence gene feature heatmaps (Figure E10G)
FeatureHeatmap(Epithelials, features.plot = c("CDKN2A", "GLB1", "SERPINE1", "IL6"), group.by = "condition",
               pt.size = 0.25, key.position = "top", 
               max.exp = 4)

save(Epithelials, file = "/Epithelials.Robj")
load(file = "/Epithelials.Robj")

# Cryo Macrophage and Epithelial Analysis ------------------------

cryo <- RenameCells(cryo, add.cell.id = "cryo")

cryo.macs <- SubsetData(cryo, ident.use = c("Macrophages"))
cryo_cells <- row.names(cryo@meta.data)

Macrophages.All.Cryo <- MergeSeurat(Macrophages.All, cryo.macs)

Macrophages.All.Cryo <- ScaleData(object = Macrophages.All.Cryo, vars.to.regress = c("nUMI", "percent.mito", "orig.ident"))
Macrophages.All.Cryo <- FindVariableGenes(object = Macrophages.All.Cryo, mean.function = ExpMean, dispersion.function = LogVMR, 
                                 x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = Macrophages.All.Cryo@var.genes)

Macrophages.All.Cryo <- RunPCA(object = Macrophages.All.Cryo, pc.genes = Macrophages.All.Cryo@var.genes, do.print = TRUE, pcs.print = 1:5, 
                      genes.print = 5, pcs.compute = 40)
Macrophages.All.Cryo <- ProjectPCA(object = Macrophages.All.Cryo, do.print = FALSE)
PCHeatmap(object = Macrophages.All.Cryo, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = Macrophages.All.Cryo, pc.use = 13:24, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = Macrophages.All.Cryo, pc.use = 25:36, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCElbowPlot(object = Macrophages.All.Cryo, num.pc = 40)

Macrophages.All.Cryo <- RunTSNE(object = Macrophages.All.Cryo, dims.use = 1:31, do.fast = TRUE)
Macrophages.All.Cryo <- FindClusters(object = Macrophages.All.Cryo, reduction.type = "pca", dims.use = 1:31, 
                            resolution = 0.05, save.SNN = TRUE, force.recalc = T)

#t-SNE of all macrophages including cryo by subcluster (Figure E11B)
TSNEPlot(object = Macrophages.All.Cryo, pt.size = 1)

#t-SNE of all macrophages including cryo by condition (Figure E11C)
TSNEPlot(object = Macrophages.All.Cryo, group.by = "condition", colors.use = c( "#EE82EE", "#00BFC4", "#F8766D"),
         pt.size = 0.5)

cryo_macrophage_freq_df <- subset(data.frame(table(Macrophages.All.Cryo@ident, Macrophages.All.Cryo@meta.data[, "condition"])), Var2 == "cryo_fibrosis")

#Macrophage pie chart by subcluster (Figure E11D)
bp<- ggplot(cryo_macrophage_freq_df, aes(x="", y=Freq/sum(Freq)*100, fill=Var1))+
  geom_bar(width = 1, stat = "identity")
pie <- bp + coord_polar("y", start=0)
pie + labs(title = "Cryobiopsy Cell Distribution by \nMacrophage Cell Cluster")

cryo.epi <- SubsetData(cryo, ident.use = c("AT2 Cells", "AT1 Cells", "Club Cells", "Ciliated Cells"))

Epithelials.All <- MergeSeurat(Epithelials, cryo.epi)
Epithelials.All <- ScaleData(object = Epithelials.All, vars.to.regress = c("nUMI", "percent.mito", "orig.ident"))
Epithelials.All <- FindVariableGenes(object = Epithelials.All, mean.function = ExpMean, dispersion.function = LogVMR, 
                                     x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = Epithelials.All@var.genes)

Epithelials.All <- RunPCA(object = Epithelials.All, pc.genes = Epithelials.All@var.genes, do.print = TRUE, pcs.print = 1:5, 
                          genes.print = 5, pcs.compute = 40)
Epithelials.All <- ProjectPCA(object = Epithelials.All, do.print = FALSE)
PCHeatmap(object = Epithelials.All, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = Epithelials.All, pc.use = 13:24, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCHeatmap(object = Epithelials.All, pc.use = 25:36, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCElbowPlot(object = Epithelials.All, num.pc = 40)

Epithelials.All <- RunTSNE(object = Epithelials.All, dims.use = 1:24, do.fast = TRUE)
Epithelials.All <- FindClusters(object = Epithelials.All, reduction.type = "pca", dims.use = 1:24, 
                                resolution = 0.05, save.SNN = TRUE, force.recalc = T)

#t-SNE of all epithelial cells by subcluster (Figure E11E)
TSNEPlot(object = Epithelials.All, pt.size = 1)

#t-SNE of all epithelial cells by condition (Figure E11F)
TSNEPlot(object = Epithelials.All, group.by = "condition", colors.use = c( "#EE82EE", "#00BFC4", "#F8766D"),
         pt.size = 0.5)

cryo_epithelial_freq_df <- subset(data.frame(table(Epithelials.All@ident, Epithelials.All@meta.data[, "condition"])), Var2 == "cryo_fibrosis")

# Pie Chart of all epithelial cell distrubtion by subcluster (Figure E11G)
bp<- ggplot(cryo_epithelial_freq_df, aes(x="", y=Freq/sum(Freq)*100, fill=Var1))+
  geom_bar(width = 1, stat = "identity")
pie <- bp + coord_polar("y", start=0)
pie + labs(title = "Cryobiopsy Cell Distribution by \nEpithelial Cell Cluster")

# Cell cycle analyisis ------------------------------------
cc.genes <- readLines(con = "/regev_lab_cell_cycle_genes.txt")

s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]

hl <- CellCycleScoring(object = hl, s.genes = s.genes, g2m.genes = g2m.genes)

#G2M score ridge plots (Figure E3A)
p <- list()
for (i in levels(hl@ident)){
  i_rep <- gsub("/","\\/" , i)
  plot <- RidgePlot(hl, ident.include = i, features.plot = "G2M.Score", group.by = "condition",
                    cols.use = c("#00BFC4", "#F8766D"), do.return = TRUE) + 
    labs(title = i_rep) + theme(plot.title = element_text(size=10), axis.title.y=element_blank())
  p[[i]] <- plot
}
do.call(grid.arrange,p)

#S score ridge plots (Figure E3B)
p <- list()
for (i in levels(hl@ident)){
  i_rep <- gsub("/","\\/" , i)
  plot <- RidgePlot(hl, ident.include = i, features.plot = "S.Score", group.by = "condition",
                    cols.use = c("#00BFC4", "#F8766D"), do.return = TRUE) + 
    labs(title = i_rep) + theme(plot.title = element_text(size=10), axis.title.y=element_blank())
  p[[i]] <- plot
}
do.call(grid.arrange,p)

