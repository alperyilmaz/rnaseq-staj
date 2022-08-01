library(DESeq2)
library(gplots)
library("RColorBrewer")
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100) ## hmcol <- heat.colors
dds <- readRDS("deseq_results")
res <- results(dds)
resOrdered <- res[order(res$padj),]

sig <- resOrdered[!is.na(resOrdered$padj) &  resOrdered$padj<0.01 &  abs(resOrdered$log2FoldChange)>=2,]
selected <- rownames(sig);selected
png(file = "heatmap2.png" )
heatmap.2(counts(dds,normalized=TRUE)[rownames(dds) %in% selected,], col = hmcol, scale="row", Rowv = TRUE, Colv = FALSE, dendrogram="row", trace="none", margin=c(4,6), cexRow=0.5, cexCol=1, keysize=1 )

 dev.off()
