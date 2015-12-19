library("DESeq2")
library("sva")

#load(".RData")
counts=read.table("intron_counts",row.names=1)
norm_counts=read.table("junction_counts",row.names=1)
names(counts)=c("A4_1","A4_2","A4_3","A4_4","A5_1","A5_2","A5_3","A5_4","wt_1","wt_2","wt_3","wt_4")
names(norm_counts)=names(counts)
trimLastChar = function (x) substr(x, 1, nchar(x) - 1)
colData=data.frame(condition=c(rep("KO",8),rep("wt",4)), samples=sapply(names(counts), trimLastChar))
dds=DESeqDataSetFromMatrix(as.matrix(counts),colData,~condition)

# dispersion on not normalized
dds_raw=estimateSizeFactors(dds)
dds_raw=estimateDispersions(dds_raw)

# normalization counts
nGrows=nrow(norm_counts)
gene_mat=as.matrix(norm_counts)
gene_mat[which(gene_mat==0)]=0.1
normFactors=gene_mat/exp(rowMeans(log(gene_mat)))
normalizationFactors(dds)=normFactors

# sva
dat = counts(dds, normalized=TRUE)
idx = rowMeans(dat) > 5
dat = dat[idx,]
mod = model.matrix(~condition,colData(dds))
mod0 = model.matrix(~1,colData(dds))
svseq = svaseq(dat, mod, mod0, n.sv=1)
dds$sv = svseq$sv
design(dds) = ~ sv + condition

# deseqqing
dds=estimateDispersions(dds)
dispersions(dds)=dispersions(dds_raw)
dds=nbinomWaldTest(dds)
res=results(dds,independentFiltering=F)#,lfcThreshold=0.3785)
res=res[order(res$pval),]
rld=rlog(dds)

# output
svg("MAplot.svg")
plotMA(res)
dev.off()
svg("PCAplot.svg",width=10,height=4)
plotPCA(rld, intgroup="samples")
dev.off()
svg("Dplot.svg")
plotDispEsts(dds)
dev.off()
normCounts=counts(dds,normalized=T)
colnames(normCounts) = names(counts)
write.table(normCounts, file="normCounts", col.names=T, row.names=T, quote=F)
write.table(res,file="results_noFC",col.names=F,row.names=T,quote=F)

save.image()
