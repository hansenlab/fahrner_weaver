###
genes <- read_csv('~/Downloads/jill_genes.csv', col_names = TRUE)
genes <- genes$`Gene ID`
genes <- toupper(genes)
genes <- genes[-which(duplicated(genes))]
gene_ids <- genes
#gene_ids <- unique(proms_mouse$gene_id[which(proms_mouse$gene_id %in% genes)])

rank_WS <- wilcox.test(res$pvalue[which(rownames(res) %in% gene_ids)], 
                       res$pvalue[-which(rownames(res) %in% gene_ids)])$statistic

length1 <- length(which(rownames(res) %in% gene_ids))

permutation_rank_KS1 <- replicate(10000, {
  indices <- sample(1:length(rownames(res)), length1)
  wilcox.test(res$pvalue[indices], res$pvalue[-indices])$statistic
})

quartz(file = "osteogenesis_ranks.pdf", width = 2.2, height = 2.2, pointsize = 8, type = "pdf")
hist(permutation_rank_KS1, col = "cornflowerblue", lty = 0, 
     breaks = 50, freq = FALSE, xlab = "Wilcoxon rank-sum test stat", cex.lab = 1.2, yaxt = 'n',
     main = "osteogenesis genes", cex.main = 1.2, font.main = 1, xlim = c(min(permutation_rank_KS1)-0.05, max(permutation_rank_KS1)+0.05), xaxt = 'n')
axis(1, at = quantile(permutation_rank_KS1, c(0.01, 0.99)), cex.axis = 1.1)
axis(2, at = c(0, 0.000012), cex.axis = 1.1)
abline(v = rank_WS, col = alpha("red", 0.6), lwd = 2.5)
#legend <- legend("topright", legend = c("random", "observed"), col = c("cornflowerblue", alpha("red", 0.6)), bty = 'n', 
#                 cex = 1.25, lty = "solid", lwd = 2.5)
dev.off()

###
proms_mouse$gene_name <- toupper(proms_mouse$gene_name)
genes <- read_csv('~/Downloads/bailey_genes.csv', col_names = FALSE)
genes <- genes$X1
genes <- toupper(genes)
genes <- genes[-which(duplicated(genes))]
gene_ids <- unique(proms_mouse$gene_id[which(proms_mouse$gene_name %in% genes)])

rank_WS <- wilcox.test(res$pvalue[which(rownames(res) %in% gene_ids)], 
                       res$pvalue[-which(rownames(res) %in% gene_ids)])$statistic

length1 <- length(which(rownames(res) %in% gene_ids))

permutation_rank_KS1 <- replicate(10000, {
  indices <- sample(1:length(rownames(res)), length1)
  wilcox.test(res$pvalue[indices], res$pvalue[-indices])$statistic
})

quartz(file = "bmp_ranks.pdf", width = 2.2, height = 2.2, pointsize = 8, type = "pdf")
hist(permutation_rank_KS1, col = "cornflowerblue", lty = 0, 
     breaks = 50, freq = FALSE, xlab = "Wilcoxon rank-sum test stat", cex.lab = 1.2, yaxt = 'n',
     main = "BMP pathway genes", cex.main = 1.2, font.main = 1, xlim = c(rank_WS-0.05, max(permutation_rank_KS1)+0.05), xaxt = 'n')
axis(1, at = c(900000, 1200000), cex.axis = 1.1)
axis(2, at = c(0, 0.000007), cex.axis = 1.1)
abline(v = rank_WS, col = alpha("red", 0.6), lwd = 2.5)
#legend <- legend("topright", legend = c("random", "observed"), col = c("cornflowerblue", alpha("red", 0.6)), bty = 'n', 
#                 cex = 1.25, lty = "solid", lwd = 2.5)
dev.off()


makeVolcanoPlot <- function(DEgenes_df, lfc_cutoff, gene_ids, color){
  tab = data.frame(logFC = DEgenes_df$log2FoldChange, negLogPval = -log10(DEgenes_df$pvalue))
  rownames(tab) <- rownames(DEgenes_df)
  #head(tab)
  #par(mar = c(5, 4, 4, 4))
  plot(tab[-which(rownames(tab) %in% gene_ids), ], pch = 16, cex = 0.5, xlab = expression(log[2]~fold~change), 
       ylab = expression(-log[10]~pvalue), col = "gray60", bty = 'l', 
       xlim = c(min(tab$logFC), max(tab$logFC)), ylim = c(0, max(tab$negLogPval)))
  points(tab[gene_ids, ], pch = 19, cex = 0.5, col = color)
  if (lfc_cutoff != "none"){
    abline(v = c(-lfc_cutoff, lfc_cutoff), col = rgb(0,0,0,0.75), lty = "longdash", lwd = 1) 
  }
}

quartz(file = "weaver_rnaseq_volcano_plots.pdf", height = 2.2, width = 6.5, pointsize = 8, type = "pdf")
par(mfrow = c(1, 3))
makeVolcanoPlot(res, "none", gene_ids, "orange") #gene_ids have been set to the BMP pathway gene ids
legend <- legend("topright", legend = c("BMP pathway genes", "other"), bty = 'n', 
                 col = c("orange", "gray60"), pch = 19, cex = 0.9)

makeVolcanoPlot(res, "none", rownames(res)[which(res$padj < 0.05)], "red")
legend <- legend("topright", legend = c("FDR < 0.05", "other"), bty = 'n', 
                 col = c("red", "gray60"), pch = 19, cex = 0.9)

makeVolcanoPlot(res, 1, rownames(res)[which(res$padj < 0.05)], "red")
legend <- legend("topright", legend = c("FDR < 0.05", "other"), bty = 'n', 
                 col = c("red", "gray60"), pch = 19, cex = 0.9)
dev.off()



