#This is the basic work flow for single cell RNA seq analysis

#1 read the data
library(Seurat)
setwd("F:/data science/R/1_PBMC_Seurat_Analysis")
hg19_data <- Read10X(data.dir = "data/")
seurat_obj <- CreateSeuratObject(counts=hg19_data, project="PBMChg19", min.cells = 3, min.features = 200)

#explore the dataset
dim(seurat_obj@assays$RNA@layers$counts) #show how many genes and samples were measured
seurat_obj@assays$RNA@layers$counts[1:32,1]
max(seurat_obj@assays$RNA@layers$counts)
min(seurat_obj@assays$RNA@layers$counts)
countmatrix <- GetAssayData(seurat_obj, slot = "counts")
count <- as.data.frame(countmatrix)
# calculate the maximal value of each gene
max_per_gene <- apply(count, 1, max)
min_per_gene <- apply(count, 1, min)
median_per_gene <- apply(count,1,median)
# Apply over rows (genes), count number of cells where expression > 0
gene_frequencies <- apply(count, 1, function(x) sum(x > 0))
gene_frequencies[1:10]
max_per_gene[1:10]
count_meta <- data.frame(genes=rownames(count),max=max_per_gene,median=median_per_gene, min=min_per_gene,frequency=gene_frequencies)
hist(median_per_gene,breaks = 100)
plot(density(median_per_gene))
# Show top 10 genes and their max values
top10_values <- sort(max_per_gene, decreasing = TRUE)[(length(max_per_gene)-9):length(max_per_gene)]
seurat_obj@assays$RNA@layers$data[1:32,1] # after normalization check the normalized data
seurat_obj@assays$RNA@layers$scale.data[1:32,1] # after scale check the normalized data.
rownames(seurat_obj)[1:5] #have a look at the rownames/genes
colnames(seurat_obj)[1:5] #have a look at the colnames/samples


# 2 QC 
#mitochondrial genes were used to check the contamination of mito/cyto
#add one column named percent.mt to meta data
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj,pattern = "^MT-")
View(seurat_obj@meta.data)
#calculate how many RNA molecules are expressed per gene on average
seurat_obj[["averageRNAmolecule"]] <-seurat_obj@meta.data$nCount_RNA/seurat_obj@meta.data$nFeature_RNA

#check the distribution of number of RNA or unique gene per cell
VlnPlot(seurat_obj,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)

#check the ratio between number of RNA and number of unique gene. 
FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# 3 Filtering
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt <5)

# 4 Normalize data by default method
seurat_obj <- NormalizeData(seurat_obj)

# 5 Identify highly variable features and visualization
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000) #default is 2000
top10 <- head(VariableFeatures(seurat_obj), 10)
plot1 <- VariableFeaturePlot(seurat_obj)
LabelPoints(plot=plot1, points = top10)

# 6 Scaling only work on the 2000 most variable genes
seurat_obj <- ScaleData(seurat_obj)

# 7 perform linear dimension reduction PCA
seurat_obj <- RunPCA(seurat_obj, features=VariableFeatures(object=seurat_obj))
print(seurat_obj[["pca"]], dims=1:5, nfeatures=5)
DimHeatmap(seurat_obj, dims=1, cells=500, balanced=TRUE)
DimHeatmap(seurat_obj, dims=2, cells=500, balanced=TRUE)
ElbowPlot(seurat_obj)
seurat_obj <- FindNeighbors(seurat_obj,dims = 1:9,k.param = 30)

#when using 2 more resolution, the last resolution results is taken as defaults cluster
seurat_obj <- FindClusters(seurat_obj,algorithm=4,resolution = c(0.2,0.4,1,2)) #pick one resolution that make sense on the graph

View(seurat_obj@meta.data) # the results of different resolution were stored in meta data

DimPlot(seurat_obj, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "pca")

seurat_obj <- RunUMAP(seurat_obj, dims = 1:9)

DimPlot(seurat_obj, label = TRUE, reduction = "umap")

# find markers
library(presto)
library(dplyr)
de <- wilcoxauc(seurat_obj, 'RNA_snn_res.0.2') # by presto method
dim(de)

top_markers <- de %>%
  group_by(group) %>%
  filter(auc > 0.6, logFC > 1, pct_in > 85, pct_out < 25, padj < 0.05) %>%
  slice_max(order_by = auc, n = 10)

#rename clusters
Idents(seurat_obj) <- "RNA_snn_res.0.2"
new_cluster_id <- c("CD4 T", "Monocytes", "NK and CD8 T", "B cells")

names(new_cluster_id) <- levels(seurat_obj)
levels((seurat_obj))
seurat_obj <- RenameIdents(seurat_obj, new_cluster_id)

VlnPlot(seurat_obj, features = c("CD3D", "CST3", "NKG7", "CD79A"))

FeaturePlot(seurat_obj,label=TRUE, features = c("CD3D", "CST3", "NKG7", "CD79A"))

RidgePlot(seurat_obj,features = c("CD3D", "CST3", "NKG7", "CD79A"))

CellScatter(seurat_obj,cell1 = "CATTTGACCACACA-1", cell2 = "CATTTGTGGGATCT-1")

DotPlot(seurat_obj,features = c("CD3D", "CST3", "NKG7", "CD79A"))
