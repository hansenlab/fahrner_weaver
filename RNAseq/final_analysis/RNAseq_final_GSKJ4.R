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
samples <- read.csv('RNAseq-samples_GSKJ4.csv')

dir <- '~/Desktop/temp/Fahrner lab/Data/Weaver/fahrner_weaver/RNAseq/final_analysis/quants_final_GSKJ4'
files <- file.path(dir, paste(samples$Sample.ID, '_quant', sep=''), 'quant.sf')

file.exists(files)
coldata <- data.frame(files, 
                      names=samples$Sample.ID, 
                      genotype=as.factor(samples$Genotype), 
                      mouse=as.factor(samples$Mouse),
                      tx=as.factor(samples$Tx),
                      alizarin=as.factor(samples$Alizarin),
                      stringsAsFactors = FALSE)

# import transcript abundances
se <- tximeta(coldata)

# summarize abundances to gene level
gse <- summarizeToGene(se)

# load SummarizedExperiment into DESeqDataSet
dds <- DESeqDataSet(gse, design = ~genotype*tx)
### END IMPORT


### BEGIN ANALYSIS 
# #dds$condition <- relevel(dds$condition, ref="WT")
# plot <- plotPCA(vst(dds), intgroup='tx')
# setEPS()
# #postscript('221016_PCA-uncollapsed.eps')
# plot + geom_label_repel(aes(label=gsub("_S.*", "", colnames(dds))), show.legend = FALSE)
# dev.off()

# collapse technical replicates
dds$sample <- factor(gsub("_S.*", "", colnames(dds)))
dds.coll <- collapseReplicates(dds, dds$sample)
stopifnot(all(rowSums(counts(dds[,which(dds$sample=="614")])) == counts(dds.coll[,1])))
dds.coll.counts <- counts(dds.coll)

# filter genes with low median counts
keep <- rowMedians(counts(dds.coll)) > 10
dds.coll.filtered <- dds.coll[keep,]
save(dds.coll.filtered, file='BigStuff/dds.coll.filtered_GSKJ4.rda')

plot <- plotPCA(vst(dds.coll.filtered), intgroup=c('genotype','tx'), returnData=TRUE)

# setEPS()
# postscript('plots_final/pca_GSKJ4.eps', width=9, height=5)
plot$group <- factor(plot$group, levels=c('WT:DMSO','WT:GSKJ4','Mut:DMSO','Mut:GSKJ4'))
group.colors <- c('deepskyblue','deepskyblue','red','red')
ggplot(data = plot, aes(x = PC1, y = PC2, color=group, shape=group)) + 
  geom_point(size = 4) + 
  scale_shape_manual(values=c(16,17,16,17)) +
  scale_color_manual(values=group.colors) + 
  xlab('PC1: 77% variance') +
  ylab('PC2: 6% variance') +
  theme_bw(base_size=18) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
# dev.off()

# setEPS()
#postscript('221016_PCA-collapsed.eps')
plot + geom_label_repel(aes(label=gsub("_S.*", "", colnames(dds.coll.filtered))), show.legend = FALSE)
# dev.off()

# find surrogate variables for batch effects using sva
dds.coll.filtered$genotype <- relevel(dds.coll.filtered$genotype, ref="WT")
dds.coll.filtered.counts <- counts(dds.coll.filtered)
mod <- model.matrix(~genotype*tx, colData(dds.coll.filtered))
mod0 <- model.matrix(~1, colData(dds.coll.filtered))
svseq <- svaseq(dds.coll.filtered.counts, mod, mod0)

# append surrogate variables
ddssva <- dds.coll.filtered
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
ddssva$SV3 <- svseq$sv[,3]
ddssva$SV4 <- svseq$sv[,4]
design(ddssva) <- formula(~SV1 + SV2 + SV3 + SV4 + genotype*tx)

# DESeq2 call
ddssva <- DESeq(ddssva)
dds <- ddssva
save(dds, file='BigStuff/dds_final_GSKJ4.rda')

mod_mat <- model.matrix(design(dds), colData(dds))
res <- results(dds, contrast = list(c('tx_GSKJ4_vs_DMSO',
                                      'genotypeMut.txGSKJ4'))) # Mut_GSKJ4 vs Mut_DMSO
res <- results(dds, contrast = list(c('genotype_Mut_vs_WT'))) # Mut_DMSO vs WT_DMSO

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

# write.csv(diffexp.subset, file="diffexp_MutGSKJ4-MutDMSO.csv") # Mut_GSKJ4 vs Mut_DMSO
write.csv(diffexp.subset, file="diffexp_MutDMSO-WtDMSO.csv") # Mut_DMSO vs WT_DMSO


aggregate(counts(dds)['ENSMUSG00000116408',]~dds$genotype+dds$tx, data=counts(dds), mean)


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
       ylim = c(0, 75))#max(-log10(res$pvalue)[which(is.na(-log10(res$pvalue))==FALSE)])))
  points(tab[gene_ids, ], pch = 19, cex = 0.75, col = color)
  circles(tab[circ_ids, ]$logFC, tab[circ_ids, ]$negLogPval, r = rep(0.12, length(circ_ids)), col='black')
  text(tab[circ_ids, ]$logFC+0.15, tab[circ_ids, ]$negLogPval+0.25, labels=circ_names, adj=0, col='black')
  abline(h = -log10(0.0222893), col= "red", lty=2, lwd=2)
  legend <- legend('topleft', legend=c('BMP pathway','other'),
                   col=c(color,'gray60'), lty=c(NA, NA), pch=c(19,19), lwd=c(NA,NA), cex=1, bty='n')
  # legend <- legend('topright', legend=c('FDR = 0.1'),
  #                  col=c('cornflowerblue'), lty=c(2), pch=c(NA), lwd=c(2), cex=1.25, bty='n')
  if (lfc_cutoff != "none"){
    abline(v = c(-lfc_cutoff, lfc_cutoff), col = rgb(0,0,0,0.75), lty = "longdash", lwd = 1) 
  }
}

fdr <- res.MutGSKJ4.MutDMSO[which(res.MutGSKJ4.MutDMSO$padj>0.1 
                                  & is.na(res.MutGSKJ4.MutDMSO$padj)==FALSE),]
fdr[which(fdr$padj==min(fdr$padj)),]

#setEPS()
#postscript('fig4c.eps', width=7, height=7)
gene_ids <- rownames(res.MutGSKJ4.MutDMSO[which(res.MutGSKJ4.MutDMSO$padj<0.1),])
makeVolcanoPlot(res.MutGSKJ4.MutDMSO, "none", gene_ids, circ_ids='', circ_names='', "orange") 
#dev.off()