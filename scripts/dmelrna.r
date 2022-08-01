library(DESeq2)

count_data <- read.csv("merged_counts",sep = "\t")
rownames(count_data) <- count_data$gene
count_data <- count_data[-c(1)]
count_data <- as.matrix(count_data) #??

coldata <- read.csv("annotation",sep = "\t", row.names = 1)
coldata$condition <- factor(coldata$condition)
coldata$type <- factor(coldata$type)

dds <- DESeqDataSetFromMatrix(countData = count_data,
                                  colData = coldata,
                                  design = ~ condition)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds <- DESeq(dds)
res <- results(dds)

write.csv(res,"deseq2_results.csv")
saveRDS(dds,"deseq_results")
