# hPBMC_scRNA-seq-classical-analysis-with-seurat
This is the example of Seurat analysis, include all the explanations of each step and plot of variable graphs.

Steps of seurat analysis:
1. input data and creat seurat object
2. explore the data, such as dim, max, min expression levels, distribution of counts, etc
3. QC by mitochondrial percentage and gene numbers per cell to exclude low quanlity cell and doubledrops
4. stand procedure including normalization, variables selection, scale, dimension reduction, findneightbors, findclusters, runumap
5. visualization by a variety of methods including dimplot, featureplot, vlnplot,etc
