##setwd("~\phase 2")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")
library(GEOquery)
library(pheatmap)
library(ggplot2)
library(Rtsne)
library(limma)

## load the data
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

## display group correlation heatmap
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
pheatmap(cor(ex), labels_col=source_names, labels_row=source_names)
dev.off()

## Perform dimensionality reduction using T-SNE
tsne_res <- Rtsne(t(ex), perplexity = 10, check_duplicates = FALSE)
tsne <- data.frame(tsne_res$Y, Group=source_names)

# store results as scatter plot
pdf("results/tSNE_results.pdf", width=15, height=15)
ggplot(tsne, aes(tsne_res$Y[, 1], tsne_res$Y[, 2], color = Group)) + geom_point(size = 3) + theme_bw()
dev.off()

## According to the results above, either CD34 or Monocytes are the candidate
source_names.small = c(
  rep("aml", 13),
  rep(0, 2),
  0,
  0,
  rep(0, 2),
  rep("monocytes", 2),
  0,
  rep(0, 5),
  rep(c(0, 0), 2),
  rep("CD34", 3),
  rep(0, 7),
  rep("aml", 2),
  0,
  rep("aml", 3),
  rep(0, 7),
  0,
  rep("monocytes", 4),
  0,
  rep(0, 7)
)

ex.small <- ex[, which(source_names.small != 0)]
source_names.small <- source_names.small[which(source_names.small != 0)]

# display correlation heatmap for CD34 and Monocyte groups
pdf("results/monocyte_cd34_corr_heatmap.pdf", width=15, height=15)
pheatmap(cor(ex.small), labels_col = source_names.small, labels_row = source_names.small)
dev.off()

# display t-SNE scatter plot for only AML, CD34, and Monocytes
ex.scale <- t(scale(t(ex.small), scale = FALSE))
pc <- prcomp(ex.scale)
pcr <- data.frame(pc$rotation[, 1:3], Group = source_names.small)
pdf("results/tSNE_CD34_Monocyte.pdf", width=15, height=15)
ggplot(pcr, aes(PC1, PC2, color = Group)) + geom_point(size = 3) + theme_bw()
dev.off()

## According to the scatter plot using visual inspection we find that
## CD34 is the best candidate for AML genes
gr <- factor(source_names)
gset$group <- gr
design <- model.matrix(~group + 0, gset)

fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(contrasts='groupaml-groupCD34', levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
tT <- subset(tT, select=c("Gene.symbol", "Gene.ID", "adj.P.Val","logFC"))

aml.up <- subset(tT, logFC > 1 & adj.P.Val < 0.05)
aml.up.genes <- unique(as.character(strsplit2(aml.up$Gene.symbol, "///")))
write.table(aml.up.genes, file="results/aml_up_genes.csv", row.names=FALSE,
    col.names = FALSE, quote = FALSE, sep=",")

aml.down <- subset(tT, logFC <- 1 & adj.P.Val < 0.05)
aml.down.genes <- unique(as.character(strsplit2(aml.down$Gene.symbol, "///")))
write.table(aml.down.genes, file="results/aml_down_genes.csv", row.names=FALSE,
    col.names = FALSE, quote = FALSE, sep=",")
