# alDFC
DFC extraction with Adaptive Lasso

## What’s this for ?




This R package `alDFC` is dedicated for the purpose, and extract a small gene subset characterizing target cell populations.

Please see also [*Discriminative feature of cells characterizes cell populations of interest by a small subset of genes.*](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009579)

## Installation

`install.packages("devtools")` if you haven't installed it.

Then run,

``` r
devtools::install_github("tfwis/alDFC")
```

## The workflow

* Step1: Prepare single cell data, especially single cell RNA-seq
* Step2: Perform standard single cell analysis procedure and set target cell cluster
* Step3: Characterize target cell cluster by discrimination
* Step4: Classify DFC set into 3 groups; *Strong*, *Weak* or *Niche*

## Tutorial using `Seurat`

### The purpose

Characterize target cluster by small subset of genes.

### 1. Set target clusters

Perform preprocesssing and standard analysis following [Seurat tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html). The data was also provided in tutorial.

```r
library(Seurat)
pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc <- NormalizeData(pbmc)
pbmc <- ScaleData(pbmc)
pbmc <- FindVariableFeatures(pbmc,nfeatures = 2000)
pbmc <- RunPCA(pbmc, features=VariableFeatures(pbmc))
pbmc <- FindNeighbors(object = pbmc)
pbmc <- FindClusters(pbmc)
pbmc <- RunUMAP(pbmc,dims = 1:10)
DimPlot(pbmc,reduction = 'umap')
```

![pbmc_umap](man/pbmc_umap.png)

Here, set cluster8 as target cluster.

### 2. DFC extraction

Using `Seurat` object and target number, `dfc()` function extract DFC subset. `target_cluster` can be given a vector; `c(3,8)`

```r
library(alDFC)
dfc_res <- dfc(pbmc, target_clusters = 8, return_Model = TRUE)
```

The solution path and the cross varidation result of each model are checked as follows.

```r
layout(matrix(1:4, ncol=2))

par(mar=c(3.5, 4, 3.5, 2))
plot(dfc_mod[['Ridge']], xlab="")
mtext(side=3, text = "Ridge", line =2.3)
mtext(side=1, text = "log(Lambda)", line =2)

par(mar=c(3.5, 4, 3.5, 2))
plot(dfc_mod[['Ridge']]$glmnet.fit,xvar = "lambda", xlab="")
mtext(side=3, text = "Ridge", line =2.3)
mtext(side=1, text = "log(Lambda)", line =2)

par(mar=c(3.5, 4, 3.5, 2))
plot(dfc_mod[['AdaLasso']], xlab="")
mtext(side=3, text = "Adaptive Lasso", line =2.3)
mtext(side=1, text = "log(Lambda)", line =2)

par(mar=c(3.5, 4, 3.5, 2))
plot(dfc_mod[['AdaLasso']]$glmnet.fit,xvar = "lambda", xlab="")
mtext(side=3, text = "Adaptive Lasso", line =2.3)
mtext(side=1, text = "log(Lambda)", line =2)
```

![regress_plot](man/regression_results.png)

### 3. DFC classification

Features in DFC subset are classified into about three groups; *Strong*, *Weak* or *Niche* feature.
`dfc_classify()` function is for the classification.
```r
dfc_class <- dfc_classify(dfc_mod$weights,pbmc)
```
