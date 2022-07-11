library(Seurat)
library(dplyr)
library(SeuratData)
library(Matrix)
library(CelliD)
library(igraph)
library(umap)
library(ggplot2)
library(msigdbr)
library(RANN)
library(graphlayouts)
library(ggraph)
library(node2vec)

clean_glist <- function(glist){
  idx <- grep("^RPS|XIST|RP11|^MT|^LINC|\\.|\\-", glist)
  glist <- glist[-idx]
  return(glist)
}

data("pbmc3k")
pbmc3k <- pbmc3k %>% 
  NormalizeData() %>% 
  ScaleData(features = rownames(pbmc3k))

pbmc3k <- RunMCA(pbmc3k)

# get cluster information 'seurat clusters'
pbmc3k <- pbmc3k %>% 
  FindVariableFeatures() %>%
  RunPCA() %>%
  FindNeighbors() %>%
  FindClusters()

# only consider protein coding gene for downstream analysis
glist <- clean_glist(rownames(pbmc3k))
# plan 1: build gene-gene network; learn patterns of GO nodes (scattered? aggregated? bi-modal?)

coembed <- rbind(pbmc3k@reductions$mca@feature.loadings[glist, ], pbmc3k@reductions$mca@cell.embeddings)
ElbowPlot(pbmc3k, ndims = 20, reduction = "mca") # use top 10 dimensions

# 
el_nn_comp <- function(dat, set1, set2, k){
  n1 <- nn2(dat[set1, ], query = dat[set2, ], k = k)
  el1 <- el_nn_search(n1)
  #print(head(el1))
  #print(set1[el1[, 2]])
  #el1[, 1] <- set2[el1[, 1]]
  #el1[, 2] <- set1[el1[, 2]]
  el_1 <- cbind(set2[el1[, 1]], set1[el1[, 2]])
  #print(head(el1))
  n2 <- nn2(dat[set2, ], query = dat[set1, ], k = k)
  el2 <- el_nn_search(n2)
  #el2[, 1] <- set1[el2[, 1]]
  #el2[, 2] <- set2[el2[, 2]]
  el_2 <- cbind(set1[el2[, 1]], set2[el2[, 2]])
  el <- rbind(el_1, el_2)
  return(el)
}

el_nn_search <- function(nn2_out){
  n.df <- nn2_out$nn.idx
  n.df <- cbind(1:nrow(n.df), n.df)
  el <- cbind(n.df[, 1], c(n.df[, -1]))
  return(el)
}
genes <- glist
cells <- colnames(pbmc3k)
el <- el_nn_comp(dat = coembed, set1 = genes, set2 = cells, k = 50)
