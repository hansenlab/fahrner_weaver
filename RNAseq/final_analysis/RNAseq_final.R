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
postscript('fig4b.eps', width = 7, height = 5)
hist(res$pvalue, 
     freq=FALSE,
     breaks = 40, 
     ylim = c(0, 3), 
     col = 'gray60', 
     xlab = "p-values", 
     main = "R684C/+ vs. +/+", 
     font.main = 1,
     cex.axis=1.5,
     cex.lab=1.5)
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
save(gene_ids, file='gene_ids_bmp.rda')

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
postscript('fig4d.eps', width=7, height=6)
hist(permutation_rank_Bmp2, col = "gray60", lty = 0, 
     breaks = 50, freq = FALSE, xlab = "Wilcoxon rank-sum test stat", cex.lab = 1.2, yaxt = 'n',
     main = "BMP pathway genes", cex.main = 1.2, font.main = 1, xlim = c(rank_WS-50000, max(permutation_rank_Bmp2)+0.05), xaxt = 'n')
axis(1, at = c(round(quantile(permutation_rank_Bmp2, c(0.01, 0.99))),rank_WS), cex.axis = 1.1)
axis(2, at = c(0, 0.000008), cex.axis = 1.1)
abline(v = rank_WS, col = "magenta", lwd = 3)
text(rank_WS+10000, 0.000007,"p < 9.9e-06", adj=0, col = "magenta")
legend <- legend("topright", legend = c("random", "observed"), 
                 col = c("gray60", "magenta"), bty = 'n', 
                 cex = 1, lty = "solid", lwd = 3)
dev.off()

makeVolcanoPlot <- function(DEgenes_df, lfc_cutoff, gene_ids, circ_ids, circ_names, color){
  tab = data.frame(logFC = DEgenes_df$log2FoldChange, negLogPval = -log10(DEgenes_df$pvalue))
  rownames(tab) <- rownames(DEgenes_df)
  par(mar=c(5,5,5,5))
  plot(tab[-which(rownames(tab) %in% gene_ids), ], pch = 19, cex = 0.75, 
       xlab = expression(log[2]~fold~change),
       ylab = expression(-log[10]~pvalue), 
       cex.axis=1.5,
       cex.lab=1.5,
       col = "gray60", bty = 'l', 
       xlim = c(min(tab$logFC), max(tab$logFC)), 
       ylim = c(0, 10))#max(-log10(res$pvalue)[which(is.na(-log10(res$pvalue))==FALSE)])))
  points(tab[gene_ids, ], pch = 19, cex = 0.75, col = color)
  circles(tab[circ_ids, ]$logFC, tab[circ_ids, ]$negLogPval, r = rep(0.12, length(circ_ids)), col='black')
  text(tab[circ_ids, ]$logFC+0.15, tab[circ_ids, ]$negLogPval+0.25, labels=circ_names, adj=0, col='black')
  abline(h = -log10(0.00139149), col= "red", lty=2, lwd=2)
  legend <- legend('topleft', legend=c('BMP pathway','other'),
                   col=c(color,'gray60'), lty=c(NA, NA), pch=c(19,19), lwd=c(NA,NA), cex=1, bty='n')
  # legend <- legend('topright', legend=c('FDR = 0.1'),
  #                  col=c('cornflowerblue'), lty=c(2), pch=c(NA), lwd=c(2), cex=1.25, bty='n')
  if (lfc_cutoff != "none"){
    abline(v = c(-lfc_cutoff, lfc_cutoff), col = rgb(0,0,0,0.75), lty = "longdash", lwd = 1) 
  }
}

setEPS()
postscript('fig4c.eps', width=7, height=7)
gene_ids <- rownames(res[which(res$padj<0.1),])
makeVolcanoPlot(res, "none", gene_ids, circ_ids='', circ_names='', "orange") 
dev.off()

setEPS()
postscript('fig4e.eps', width=7, height=7) #set gene_ids to bmp pathway, code above
circ_ids <- c('ENSMUSG00000051279', 
              'ENSMUSG00000020427',
              'ENSMUSG00000021540',
              'ENSMUSG00000039153',
              'ENSMUSG00000060284')
circ_names <- c('Gdf6',
                'Igfbp3',
                'Smad5',
                'Runx2',
                'Sp7')
makeVolcanoPlot(res, "none", gene_ids, circ_ids, circ_names, "magenta")
text(-6,-log10(0.00139149)+0.5, "FDR = 0.1", col="red")
dev.off()

fdr <- res[which(res$padj>0.1 & is.na(res$padj)==FALSE),]
fdr[which(fdr$padj==min(fdr$padj)),]

plot <- plotPCA(vst(dds.coll.filtered), returnData=TRUE)
plot$group <- relevel(plot$group, ref='WT')
group.colors <- c('#0F80FF','#FB0207')
setEPS()
postscript('Fig4a.eps', width=9, height=5)
ggplot(data = plot, aes(x = PC1, y = PC2, color=group)) + 
  geom_point(size = 4) + 
  scale_color_manual(labels=c('+/+','R684C/+'), values=group.colors) + 
  xlab('PC1: 75% variance') +
  ylab('PC2: 14% variance') +
  theme_bw(base_size=18) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()
