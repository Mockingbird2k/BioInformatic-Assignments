setwd("~/Code/bioinformatics/Project/phase1/")
library(Biobase)
library(GEOquery)
library(limma)
library(pheatmap)
library(umap)
library(ggplot2)
library(plyr)
library(Rtsne) # for t-SNE dimensionality reduction

series <- "GSE48558"
platform <- "GPL6244"

gset <- getGEO(series, GSEMatrix=TRUE, AnnotGPL=TRUE, destdir="./data/")

if (length(gset) > 1) idx <- grep(platform, attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

gsms <- paste0(
  "1111111111111XXXXXXXXXXXXXXXXXXXXXXXXXXX0XXX0XXXXX",
  "XXXXXXXXXXXXXXXXXX0X0XXX0X0000X0XX00XX00X0X0X0X0X0",
  "XXX0XXX0XXXXXXXXXXXXXXXXXXXXXXXXXXXXX0000000110111",
  "00000000000000000000")
sml <- strsplit(gsms, split = "")[[1]]
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[, sel]
gr <- c()

for (i in 1:length(sml)) {
  if (sml[i] == "1") gr[i] <- "healthy" else gr[i] <- "aml"
}

ex <- exprs(gset)

## Test data quality
pdf("results/boxplot.pdf", width=length(sml))
boxplot(ex)
dev.off()

pdf("results/cor_heatmap.pdf", width=15, height=15)
pheatmap(cor(ex), labels_row=gr, labels_col=gr)
dev.off()

## Dimensionality Reduction
# PCA
ex.scale <- t(scale(t(ex), scale=FALSE))
pc <- prcomp(ex.scale)
pcr <- data.frame(pc$rotation[, 1:3], Group=gr)
pdf("results/pca.pdf")
ggplot(pcr, aes(PC1, PC2, color=Group)) + geom_point(size=3) + theme_bw()
dev.off()

# MDS
ex.t <- t(ex)
d <- as.matrix(dist(ex.t))
fit <- cmdscale(d, eig=TRUE, k=2)
x <- fit$points[, 1]
y <- fit$points[, 2]
mds <- data.frame(x, y, Group=gr)
pdf("results/mds.pdf")
ggplot(mds, aes(x, y, color = Group)) + geom_point(size=3) + theme_bw()
dev.off()

# t-SNE
tsne_results <- Rtsne(ex.t, perplexity=10, check_duplicates=FALSE)
tsne <- data.frame(tsne_results$Y, Group=gr)
pdf("results/tSNE.pdf")
ggplot(tsne, aes(tsne_results$Y[, 1], tsne_results$Y[, 2], color=Group)) +
  geom_point(size=3) + theme_bw()
dev.off()

## checking source name correlation
source_names = c(
  rep("aml", 13),
  rep("granulocytes", 2),
  "Bcells",
  "Tcells",
  rep("granulocytes", 2),
  rep("monocytes", 2),
  "Bcells",
  rep("Tcells", 5),
  rep(c("Bcells", "Tcells"), 2),
  rep("CD34", 3),
  rep("granulocytes", 7),
  rep("aml", 2),
  "Tcells",
  rep("aml", 3),
  rep("Bcells", 7),
  "Tcells",
  rep("monocytes", 4),
  "granulocytes",
  rep("Tcells", 7)
)
pdf("results/group_cor_heatmap.pdf", width=15, height=15)
pheatmap(cor(ex), labels_row=source_names, labels_col=source_names)
dev.off()
