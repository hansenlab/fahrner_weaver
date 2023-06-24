library(TxDb.Mmusculus.UCSC.mm10.ensGene)
library(tximport)
library(DESeq2)
library(sva)
library(biomaRt)
library(ggplot2)
library(ggrepel)
library(tximeta)
library(SummarizedExperiment)
library(edgeR)
library(limma)
library(statmod)
library(tidyverse)


### BEGIN IMPORT
samples <- read.csv('RNAseq-samples.csv')

dir <- '~/Desktop/temp/Fahrner lab/Data/Weaver/fahrner_weaver/RNAseq/GSE138980/quants_final'
files <- file.path(dir, paste(samples$Sample.ID, '_quant', sep=''), 'quant.sf')

file.exists(files)
coldata <- data.frame(files, names=samples$Sample.ID, condition=samples$Genotype, stringsAsFactors = FALSE)

# import transcript abundances
se <- tximeta(coldata)

# summarize abundances to gene level
gse <- summarizeToGene(se)

# load SummarizedExperiment into DESeqDataSet
dds <- DESeqDataSet(gse, design=~condition)
### END IMPORT


### BEGIN ANALYSIS 
#dds$condition <- relevel(dds$condition, ref="Ctrl")
plot <- plotPCA(vst(dds))
setEPS()
postscript('GSE138980_PCA-uncollapsed.eps')
plot + geom_label_repel(aes(label=gsub("-.*", "", colnames(dds))), show.legend = FALSE)
dev.off()

# # collapse technical replicates
# dds$sample <- factor(gsub("-.*", "", colnames(dds)))
# dds.coll <- collapseReplicates(dds, dds$sample)
# stopifnot(all(rowSums(counts(dds[,which(dds$sample=="705")])) == counts(dds.coll[,1])))
# dds.coll.counts <- counts(dds.coll)

# # filter genes with low median counts
# keep <- rowMedians(counts(dds.coll)) > 10
# dds.coll.filtered <- dds.coll[keep,]
keep <- rowMedians(counts(dds)) > 10
dds.filtered <- dds[keep,]

plot <- plotPCA(vst(dds.filtered))

setEPS()
postscript('GSE138980_PCA-filtered.eps')
plot + geom_label_repel(aes(label=gsub("-.*", "", colnames(dds.filtered))), show.legend = FALSE)
dev.off()

# # find surrogate variables for batch effects using sva
# dds.coll.filtered$condition <- relevel(dds.coll.filtered$condition, ref="WT")
# dds.coll.filtered.counts <- counts(dds.coll.filtered)
# mod <- model.matrix(~condition, colData(dds.coll.filtered))
# mod0 <- model.matrix(~1, colData(dds.coll.filtered))
# svseq <- svaseq(dds.coll.filtered.counts, mod, mod0)
dds.filtered$condition <- relevel(dds.filtered$condition, ref="Ctrl")
dds.filtered.counts <- counts(dds.filtered)
mod <- model.matrix(~condition, colData(dds.filtered))
mod0 <- model.matrix(~1, colData(dds.filtered))
svseq <- svaseq(dds.filtered.counts, mod, mod0)

# append surrogate variables
ddssva <- dds.filtered
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
# ddssva$SV3 <- svseq$sv[,3]
# design(ddssva) <- formula(~SV1 + SV2 + SV3 + condition)
design(ddssva) <- formula(~SV1 + SV2 + condition)

# DESeq2 call
ddssva <- DESeq(ddssva)
dds <- ddssva
save(dds, file='dds_GSE138980.rda')
res <- results(dds)

diffexp.subset <- as.data.frame(res[which(res$padj <0.1),])
### END ANALYSIS

# gene annotation
ensembl.id <- rownames(diffexp.subset)

ensembl102.mmusculus <- useEnsembl(biomart='genes', dataset='mmusculus_gene_ensembl', version=102)
mgi <- getBM(attributes=c('ensembl_gene_id','mgi_symbol','mgi_description'),
             filters='ensembl_gene_id',
             values=ensembl.id,
             mart=ensembl102.mmusculus)
diffexp.subset$ensembl_gene_id <- ensembl.id
diffexp.subset <- merge(diffexp.subset, mgi, by='ensembl_gene_id')
diffexp.subset <- diffexp.subset[order(diffexp.subset$log2FoldChange, decreasing=TRUE),]
write.csv(diffexp.subset, file="diffexp_GSE138980.csv")

hist(res$pvalue, freq=FALSE, breaks=50, 
     main='GSK126 vs Ctrl',
     xlab='p-values',
     ylim=c(0,8))

setwd("~/Desktop/temp/Fahrner lab/Data/Weaver/fahrner_weaver/RNAseq/final_analysis/BigStuff")
load(file='dds_final.rda')
dds_final <- dds
setwd("~/Desktop/temp/Fahrner lab/Data/Weaver/fahrner_weaver/RNAseq/GSE138980/BigStuff")
load(file='dds_GSE138980.rda')
dds_GSE138980 <- dds
rm(dds)
res_final <- results(dds_final)
res_GSE138980 <- results(dds_GSE138980)
gene_ids_final_less <- rownames(res_final[which(res_final$padj<0.1),])
gene_ids_final_more <- rownames(res_final[which(res_final$padj>0.1),])
gene_ids_GSE138980_less <- rownames(res_GSE138980[which(res_GSE138980$padj<0.1),])
gene_ids_GSE138980_more <- rownames(res_GSE138980[which(res_GSE138980$padj>0.1),])

setEPS()
postscript(file='../density_cond-D14.eps', width=2.36, height=1.77)
par(mar=c(1.35,1.1,0,0), mgp=c(0.5,0.075,0))
plot(density(res_final[which(rownames(res_final) %in% gene_ids_GSE138980_less),]$pvalue, 
             from=0.01,to=0.97, bw=0.025), 
     asp=NA, lwd=2, col='orange', main='', bty='n', yaxt='n', xaxt='n', ylab='', xlab='')
lines(density(res_final[which(rownames(res_final) %in% gene_ids_GSE138980_more),]$pvalue, 
             from=0.01,to=0.97, bw=0.025),
     lwd=2, col='black')
axis(1, at=c(0,1), tck=-0.02, cex.axis=0.7)
axis(2, at=c(.6,2.1), tck=-0.02, cex.axis=0.7, pos=-0.02)
mtext('p-values, D14', side=1, line=0.3, cex=0.7)
mtext('Density', side=2, line=0.3, cex=0.7)
legend('topright', legend=c('p-adj > 0.1, GSK-126', 'p-adj < 0.1, GSK-126'),
       col=c("black", "orange"), lty=c(1,1), lwd=2, cex=0.7,
       box.lty=0)
dev.off()

setEPS()
postscript(file='../density_cond-GSE138980.eps', width=2.36, height=1.77)
par(mar=c(1.35,1.1,0,0), mgp=c(0.5,0.075,0))
plot(density(res_GSE138980[which(rownames(res_GSE138980) %in% gene_ids_final_less),]$pvalue, 
             from=0.01,to=0.98, bw=0.025), 
     asp=NA, lwd=2, col='orange', main='', bty='n', yaxt='n', xaxt='n', ylab='', xlab='')
lines(density(res_GSE138980[which(rownames(res_GSE138980) %in% gene_ids_final_more),]$pvalue, 
              from=0.01,to=0.98, bw=0.025),
      lwd=2, col='black')
axis(1, at=c(0,1), tck=-0.02, cex.axis=0.7)
axis(2, at=c(.5,4.2), tck=-0.02, cex.axis=0.7, pos=-0.02)
mtext('p-values, GSK-126', side=1, line=0.3, cex=0.7)
mtext('Density', side=2, line=0.3, cex=0.7)
legend('topright', legend=c('p-adj > 0.1, D14', 'p-adj < 0.1, D14'),
       col=c("black", "orange"), lty=c(1,1), lwd=2, cex=0.7,
       box.lty=0)
dev.off()