library("DESeq2")

counts=read.table("count.txt",row.names=1)
norm_counts=read.table("norm.txt",row.names=1)
samples = c("A1_1","A1_2","A1_3", "A1_4", "A4_1","A4_2","A4_3", "A4_4", "A5_1","A5_2","A5_3", "A5_4", "wt_1","wt_2","wt_3", "wt_4")

#----------------------------------------------------------------------------------
# threshhold for low values
#----------------------------------------------------------------------------------


threshhold = function(row){

  return(sum(norm_counts[row, c(1,2,3,4)]) > 20 && 
         sum(norm_counts[row, c(5,6,7,8)]) > 20 && 
         sum(norm_counts[row, c(9,10,11,12)]) > 20 && 
         sum(norm_counts[row, c(13,14,15,16)]) > 20 && 
         sum(counts[row, c(1,2,3,4)]) > 20 &&
         sum(counts[row, c(5,6,7,8)]) > 20 && 
         sum(counts[row, c(9,10,11,12)]) > 20 && 
         sum(counts[row, c(13,14,15,16)]) > 20)
    
}


filter = sapply(seq(nrow(counts)), threshhold)

counts = counts[filter,]
norm_counts = norm_counts[filter,]

names(counts) = samples
names(norm_counts) = samples
colData=data.frame(row.names= samples, condition=c(rep("FUSK_KO", 12), rep("control", 4)))
dds_norm=DESeqDataSetFromMatrix(countData=as.matrix(norm_counts),colData=colData,design=~condition)

# normalization sizeFacs+counts
nGrows=nrow(counts)
gene_mat=as.matrix(counts)
gene_mat[which(gene_mat==0)]=0.1
normFactors=gene_mat/exp(rowMeans(log(gene_mat)))
normalizationFactors(dds_norm)=normFactors

dds_norm=DESeq(dds_norm)
res=results(dds_norm,independentFiltering=F)
res=res[order(res$pval),]
rld=rlog(dds_norm)

normCounts=counts(dds_norm,normalized=T)
write.table(res,file="results",col.names=T,row.names=T,quote=F)

pdf('DESeq2_plots.pdf')
plotMA(res, main = "MA-plot")


plotPCA(rld)

dds = estimateSizeFactors(dds_norm)
dds = estimateDispersions(dds)
plotDispEsts(dds)

dev.off()

#----------------------------------------------------------------------------------
# get reads of top results
#----------------------------------------------------------------------------------


top20 = head(rownames(res),20)
top_res = data.frame(row.names = top20)

for (col in seq(16)){
  
  count_name = paste(samples[col], "_count", sep = "")
  norm_name = paste(samples[col], "_norm" , sep = "")
  
  top_res[count_name] = counts[top20, col]
  top_res[norm_name] = norm_counts[top20, col]
    
}

write.table(top_res,file="top_results", col.names=T,row.names=T,quote=F)

#----------------------------------------------------------------------------------
# plots for overview of data
#----------------------------------------------------------------------------------

pdf("hist_count.pdf")
par(mfrow = c(4,4))

for (col in seq(16)){
  
  hist(log10(counts[,col]), freq = F, main = sampels[col])
}

dev.off()


pdf("hist_norm.pdf")
par(mfrow = c(4,4))

for (col in seq(16)){
  
  hist(log10(norm_counts[,col]), freq = F, main = sampels[col])
}

dev.off()

pdf("scatter_plots.pdf")
par(mfrow = c(4,4))

for (col in seq(16)){
  
  plot(counts[,col], norm_counts[,col] ,xlim = c(0, 2000), ylim = c(0,2000), pch = 20, cex = 0.5)
}

dev.off()

plot(log10(counts[,1]), log10(norm_counts[,1]), pch = 20, cex = 0.5)


