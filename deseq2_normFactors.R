library("DESeq2")

counts = read.table("intron_counts", row.names=1, header=T)
norm_counts = read.table("junction_counts", row.names=1, header=T)
trimLastChar  =  function (x) substr(x, 1, nchar(x) - 1)
colData = data.frame(condition=strsplit(commandArgs(T)[1])[[1]])
dds = DESeqDataSetFromMatrix(as.matrix(counts), colData, ~condition)

# dispersion on not normalized
dds_raw = estimateSizeFactors(dds)
dds_raw = estimateDispersions(dds_raw)

# normalization counts
nGrows = nrow(norm_counts)
gene_mat = as.matrix(norm_counts)
gene_mat[which(gene_mat==0)] = 0.1
normFactors = gene_mat/exp(rowMeans(log(gene_mat)))
normalizationFactors(dds) = normFactors

# deseqqing
dds = estimateDispersions(dds)
dispersions(dds) = dispersions(dds_raw)
dds = nbinomWaldTest(dds)
res = results(dds, independentFiltering=F)
res = res[order(res$pval), ]
rld = rlog(dds)

# output
svg("MAplot.svg")
plotMA(res)
dev.off()
svg("PCAplot.svg",width=10,height=4)
plotPCA(rld, intgroup="condition")
dev.off()
svg("Dplot.svg")
plotDispEsts(dds)
dev.off()
normCounts = counts(dds,normalized=T)
colnames(normCounts)  =  names(counts)
write.table(normCounts, file="normCounts", col.names=T, row.names=T, quote=F)
write.table(res,file="results",col.names=F,row.names=T,quote=F)

save.image()
