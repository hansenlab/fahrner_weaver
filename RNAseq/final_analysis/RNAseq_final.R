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

dir <- '~/Desktop/temp/Fahrner lab/Data/Weaver/fahrner_weaver/RNAseq/final_analysis/quants_final'
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
#dds$condition <- relevel(dds$condition, ref="WT")
plot <- plotPCA(vst(dds))
setEPS()
postscript('221016_PCA-uncollapsed.eps')
plot + geom_label_repel(aes(label=gsub("-.*", "", colnames(dds))), show.legend = FALSE)
dev.off()

# collapse technical replicates
dds$sample <- factor(gsub("-.*", "", colnames(dds)))
dds.coll <- collapseReplicates(dds, dds$sample)
stopifnot(all(rowSums(counts(dds[,which(dds$sample=="705")])) == counts(dds.coll[,1])))
dds.coll.counts <- counts(dds.coll)

# filter genes with low median counts
keep <- rowMedians(counts(dds.coll)) > 10
dds.coll.filtered <- dds.coll[keep,]

plot <- plotPCA(vst(dds.coll.filtered))

setEPS()
postscript('221016_PCA-collapsed.eps')
plot + geom_label_repel(aes(label=gsub("-.*", "", colnames(dds.coll.filtered))), show.legend = FALSE)
dev.off()

# find surrogate variables for batch effects using sva
dds.coll.filtered$condition <- relevel(dds.coll.filtered$condition, ref="WT")
dds.coll.filtered.counts <- counts(dds.coll.filtered)
mod <- model.matrix(~condition, colData(dds.coll.filtered))
mod0 <- model.matrix(~1, colData(dds.coll.filtered))
svseq <- svaseq(dds.coll.filtered.counts, mod, mod0)

# append surrogate variables
ddssva <- dds.coll.filtered
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
# ddssva$SV3 <- svseq$sv[,3]
# design(ddssva) <- formula(~SV1 + SV2 + SV3 + condition)
design(ddssva) <- formula(~SV1 + SV2 + condition)

# DESeq2 call
ddssva <- DESeq(ddssva)
dds <- ddssva
save(dds, file='dds_final.rda')
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
write.csv(diffexp.subset, file="diffexp_final.csv")

### PLOTS
# p-value density line
setEPS()
postscript('221016_p-val-density.eps')
plot(density(res$pvalue, from =0, to=1, adjust=1), 
     col='cornflowerblue', lwd=3, xlab='Unadjusted p-values', main='')
dev.off()

# p-value density histogram
setEPS()
postscript('221016_p-val-hist.eps')
hist(res$pvalue, 
     freq=FALSE,
     breaks = 40, 
     ylim = c(0, 3), 
     col = 'cornflowerblue', 
     xlab = "p-values", 
     main = "R684C/+ vs. +/+", 
     font.main = 1)
dev.off()

# Wilcoxon rank-sum, osteogenesis genes
genes_ost <- read_csv('~/Desktop/temp/Fahrner lab/Data/Weaver/fahrner_weaver/RNAseq/final_analysis/jill_genes.csv', col_names = TRUE)
genes_ost <- genes_ost$`Gene ID`
genes_ost <- toupper(genes_ost)
genes_ost <- genes_ost[-which(duplicated(genes_ost))]
gene_ids <- genes_ost

rank_WS <- wilcox.test(res$pvalue[which(rownames(res) %in% gene_ids)], 
                       res$pvalue[-which(rownames(res) %in% gene_ids)])$statistic

length1 <- length(which(rownames(res) %in% gene_ids))

permutation_rank_ost <- replicate(10000, {
  indices <- sample(1:length(rownames(res)), length1)
  wilcox.test(res$pvalue[indices], res$pvalue[-indices])$statistic
  })

save(permutation_rank_ost, file='permutation_rank_ost.rda')
prop.table(table(permutation_rank_ost<rank_WS))

setEPS()
postscript('221016_osteo.eps')
hist(permutation_rank_ost, col = "cornflowerblue", lty = 0, 
     breaks = 50, freq = FALSE, xlab = "Wilcoxon rank-sum test stat", cex.lab = 1.2, yaxt = 'n',
     main = "osteogenesis genes", cex.main = 1.2, font.main = 1, xlim = c(min(permutation_rank_ost)-0.05, max(permutation_rank_ost)+0.05), xaxt = 'n')
axis(1, at = round(quantile(permutation_rank_ost, c(0.01, 0.99))), cex.axis = 1.1)
axis(2, at = c(0, 0.000012), cex.axis = 1.1)
abline(v = rank_WS, col = "red", lwd = 2.5)
legend <- legend("topright", legend = c("random", "observed"), 
                 col = c("cornflowerblue", "red"), bty = 'n',
                 cex = 1, lty = "solid", lwd = 2.5)
dev.off()

# Wilcoxon rank-sum, Bmp2 pathway
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
edb_mouse <- EnsDb.Mmusculus.v79
proms_mouse <- promoters(edb_mouse, filter = TxBiotypeFilter("protein_coding"), upstream = 2000, downstream = 2000, columns = c("gene_name", "tx_id", "gene_id"))
genome(seqinfo(proms_mouse)) <- "mm10"
seqlevelsStyle(proms_mouse) <- "ucsc"
proms_mouse <- proms_mouse[which(seqnames(proms_mouse) %in% seqnames(Mmusculus)[1:21])]
proms_mouse$gene_name <- toupper(proms_mouse$gene_name)

genes_bmp2 <- read_csv('~/Desktop/temp/Fahrner lab/Data/Weaver/fahrner_weaver/RNAseq/final_analysis/bailey_genes_CWG.csv', col_names = FALSE)
genes_bmp2 <- genes_bmp2$X1
genes_bmp2 <- toupper(genes_bmp2)
genes_bmp2 <- genes_bmp2[-which(duplicated(genes_bmp2))]
gene_ids <- unique(proms_mouse$gene_id[which(proms_mouse$gene_name %in% genes_bmp2)])

rank_WS <- wilcox.test(res$pvalue[which(rownames(res) %in% gene_ids)], 
                       res$pvalue[-which(rownames(res) %in% gene_ids)])$statistic

length1 <- length(which(rownames(res) %in% gene_ids))

permutation_rank_Bmp2 <- replicate(10000, {
  indices <- sample(1:length(rownames(res)), length1)
  wilcox.test(res$pvalue[indices], res$pvalue[-indices])$statistic
  })

save(permutation_rank_Bmp2, file='permutation_rank_Bmp2.rda')
prop.table(table(permutation_rank_Bmp2<rank_WS))

setEPS()
postscript('221016_Bmp2.eps')
hist(permutation_rank_Bmp2, col = "cornflowerblue", lty = 0, 
     breaks = 50, freq = FALSE, xlab = "Wilcoxon rank-sum test stat", cex.lab = 1.2, yaxt = 'n',
     main = "BMP pathway genes", cex.main = 1.2, font.main = 1, xlim = c(rank_WS-50000, max(permutation_rank_Bmp2)+0.05), xaxt = 'n')
axis(1, at = c(round(quantile(permutation_rank_Bmp2, c(0.01, 0.99))),rank_WS), cex.axis = 1.1)
axis(2, at = c(0, 0.000008), cex.axis = 1.1)
abline(v = rank_WS, col = "red", lwd = 2.5)
legend <- legend("topright", legend = c("random", "observed"), 
                 col = c("cornflowerblue", "red"), bty = 'n', 
                 cex = 1, lty = "solid", lwd = 2.5)
dev.off()

gene_ids <- rownames(res[which(res$padj<0.1),])
makeVolcanoPlot <- function(DEgenes_df, lfc_cutoff, gene_ids, color){
  tab = data.frame(logFC = DEgenes_df$log2FoldChange, negLogPval = -log10(DEgenes_df$pvalue))
  rownames(tab) <- rownames(DEgenes_df)
  par(mar=c(5,5,5,5))
  plot(tab[-which(rownames(tab) %in% gene_ids), ], pch = 16, cex = 0.75, 
       xlab = expression(log[2]~fold~change),
       ylab = expression(-log[10]~pvalue), 
       cex.axis=1.5,
       cex.lab=1.5,
       col = "gray60", bty = 'l', 
       xlim = c(min(tab$logFC), max(tab$logFC)), 
       ylim = c(0, max(-log10(res$pvalue)[which(is.na(-log10(res$pvalue))==FALSE)])))
  points(tab[gene_ids, ], pch = 19, cex = 0.75, col = color)
  # abline(h = -log10(0.00139149), col= "cornflowerblue", lty=2, lwd=2)
  legend <- legend('topleft', legend=c('FDR < 0.1','other'),
                   col=c('red','gray60'), lty=c(NA, NA), pch=c(19,19), lwd=c(NA,NA), cex=1.25, bty='n')
  # legend <- legend('topright', legend=c('FDR = 0.1'),
  #                  col=c('cornflowerblue'), lty=c(2), pch=c(NA), lwd=c(2), cex=1.25, bty='n')
  if (lfc_cutoff != "none"){
    abline(v = c(-lfc_cutoff, lfc_cutoff), col = rgb(0,0,0,0.75), lty = "longdash", lwd = 1) 
  }
}

setEPS()
postscript('221019_WS_volcano.eps', width=7, height=6)
makeVolcanoPlot(res, "none", gene_ids, "red") #gene_ids have been set to the chondrogenesis gene ids
dev.off()

fdr <- res[which(res$padj>0.1 & is.na(res$padj)==FALSE),]
fdr[which(fdr$padj==min(fdr$padj)),]
