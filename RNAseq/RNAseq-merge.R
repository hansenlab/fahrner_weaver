library(data.table)

LqLc <- data.frame(fread("Leandros-quant-Leandros-code.csv",
              select=c("mgi_symbol", "log2FoldChange"),
              col.names=c("mgi_symbol", "LqLc")))
CqLc <- data.frame(fread("Christine-quant-Leandros-code.csv",
              select=c("mgi_symbol", "log2FoldChange"),
              col.names=c("mgi_symbol", "CqLc")))
CqCc <- data.frame(fread("Christine-quant-Christine-code.csv",
              select=c("mgi_symbol", "log2FoldChange"),
              col.names=c("mgi_symbol", "CqCc")))
LqCc <- data.frame(fread("Leandros-quant-Christine-code.csv",
              select=c("mgi_symbol", "log2FoldChange"),
              col.names=c("mgi_symbol", "LqCc")))

df.list <- list(LqLc, CqLc, CqCc, LqCc)
df.combined <- Reduce(function(x,y) merge(x, y, by='mgi_symbol', all=TRUE), df.list)

write.csv(df.combined, file="RNAseq-merge-log2FC.csv")

# import and rename dds files
load(file='dds_Christine-quant-Christine-code.rda')
dds_CqCc <- dds
rm(dds)
load(file='dds_Leandros-quant-Christine-code.rda')
dds_LqCc <- dds
rm(dds)
load(file='dds_Leandros-quant-Leandros-code.rda')
dds_LqLc <- dds
rm(dds)
load(file='dds_Christine-quant-Leandros-code.rda')
dds_CqLc <- dds
rm(dds)
load(file='dds_Frankenstein.rda')
dds_Frankenstein <- dds
load(file="dds_L-import.rda")
dds_L.import <- dds
load(file="dds_CWG-import.rda")
dds_CWG.import <- dds
load(file="dds_L-import_Cq.rda")
dds_L.import_Cq <- dds
load(file="dds_CWG-import_Cq.rda")
dds_CWG.import_Cq <- dds
load(file="dds_L-import_Cq_EnsDb102.rda")
dds_L.import_Cq_EnsDb102 <- dds
rm(dds)
res_CqCc <- results(dds_CqCc)
res_LqCc <- results(dds_LqCc)
res_LqLc <- results(dds_LqLc)
res_CqLc <- results(dds_CqLc)
res_Frankenstein <- results(dds_Frankenstein)
res_L.import <- results(dds_L.import)
res_CWG.import <- results(dds_CWG.import)
res_L.import_Cq <- results(dds_L.import_Cq)
res_CWG.import_Cq <- results(dds_CWG.import_Cq)
res_L.import_Cq_EnsDb102 <- results(dds_L.import_Cq_EnsDb102)

# log2FC scatter plots for LqLc vs CqLc 
png(filename = 'LqLc_CqLc_all_col.png', width=700, height=700)

LqLc_CqLc <- merge(as.data.frame(res_LqLc), as.data.frame(res_CqLc), by="row.names", suffixes=c('.LqLc','.CqLc'))

# rbPal <- colorRampPalette(c('red','yellow','green','blue'))
# LqLc_CqLc$col <- rbPal(20)[as.numeric(cut(log(LqLc_CqLc$baseMean.LqLc), breaks=30))]

plot(LqLc_CqLc$log2FoldChange.LqLc, LqLc_CqLc$log2FoldChange.CqLc,
     pch=19, cex=0.5, col=LqLc_CqLc$col,
     xlim=c(-8.5,3), ylim=c(-8.5,3),
     xlab='log2FC, LqLc', ylab='log2FC, CqLc', main='LqLc vs CqLc, all genes')
abline(a=0, b=1, col='gray', lty=2)
abline(h=0)
abline(v=0)
dev.off()

png(filename = 'LqLc_CqLc_padj_all.png', width=700, height=700)
plot(LqLc_CqLc$padj.LqLc, LqLc_CqLc$padj.CqLc,
     pch=19, cex=0.5, col='black',
     xlim=c(0,1), ylim=c(0,1),
     xlab='padj, LqLc', ylab='padj, CqLc', main='LqLc vs CqLc, all genes')
abline(a=0, b=1, col='red', lty=2)
abline(h=0)
abline(v=0)
dev.off()

png(filename = 'LqLc_CqLc_sig_or.png', width=700, height=700)
LqLc_CqLc_sig_or <- LqLc_CqLc[which(LqLc_CqLc$padj.LqLc <0.1 | LqLc_CqLc$padj.CqLc <0.1),]
plot(LqLc_CqLc_sig_or$log2FoldChange.LqLc, LqLc_CqLc_sig_or$log2FoldChange.CqLc,
     pch=19, cex=0.5, col='black',
     xlim=c(-8.5,3), ylim=c(-8.5,3),
     xlab='log2FC, LqLc', ylab='log2FC, CqLc', main='LqLc vs CqLc, sig genes or')
abline(a=0, b=1, col='red', lty=2)
abline(h=0)
abline(v=0)
dev.off()

png(filename = 'LqLc_CqLc_fil.png', width=700, height=700)
LqLc_CqLc_fil <- LqLc_CqLc[which(LqLc_CqLc$baseMean.LqLc >100 & LqLc_CqLc$baseMean.CqLc >100),]
plot(LqLc_CqLc_fil$log2FoldChange.LqLc, LqLc_CqLc_fil$log2FoldChange.CqLc,
     pch=19, cex=0.5, col=LqLc_CqLc$col,
     xlim=c(-5,3), ylim=c(-5,3),
     xlab='log2FC, LqLc', ylab='log2FC, CqLc', main='LqLc vs CqLc, exp filtered > 100')
abline(a=0, b=1, col='gray', lty=2)
abline(h=0)
abline(v=0)
dev.off()

# log2FC scatter plots for LqLc vs LqCc
png(filename = 'LqLc_LqCc_all.png', width=700, height=700)
LqLc_LqCc <- merge(as.data.frame(res_LqLc), as.data.frame(res_LqCc), by="row.names", suffixes=c('.LqLc','.LqCc'))
plot(LqLc_LqCc$log2FoldChange.LqLc, LqLc_LqCc$log2FoldChange.LqCc,
     pch=19, cex=0.5,
     xlim=c(-6,3), ylim=c(-6,3),
     xlab='log2FC, LqLc', ylab='log2FC, LqCc', main='LqLc vs LqCc, all genes')
abline(a=0, b=1, col='red', lty=2)
abline(h=0)
abline(v=0)
dev.off()

LqLc_LqCc[which(LqLc_LqCc$log2FoldChange.LqLc<(-4)),]

png(filename = 'LqLc_LqCc_padj_all.png', width=700, height=700)
plot(LqLc_LqCc$padj.LqLc, LqLc_LqCc$padj.LqCc,
     pch=19, cex=0.5, col='black',
     xlim=c(0,1), ylim=c(0,1),
     xlab='padj, LqLc', ylab='padj, LqCc', main='LqLc vs LqCc, all genes')
abline(a=0, b=1, col='red', lty=2)
abline(h=0)
abline(v=0)
dev.off()

png(filename = 'LqLc_LqCc_sig.png', width=700, height=700)
LqLc_LqCc_sig <- LqLc_LqCc[which(LqLc_LqCc$padj.LqLc <0.1 & LqLc_LqCc$padj.LqCc <0.1),]
plot(LqLc_LqCc_sig$log2FoldChange.LqLc, LqLc_LqCc_sig$log2FoldChange.LqCc,
     pch=19, cex=0.5,
     xlim=c(-6,3), ylim=c(-6,3),
     xlab='log2FC, LqLc', ylab='log2FC, LqCc', main='LqLc vs LqCc, sig genes')
abline(a=0, b=1, col='red', lty=2)
abline(h=0)
abline(v=0)
dev.off()

png(filename = 'LqLc_LqCc_fil.png', width=700, height=700)
LqLc_LqCc_fil <- LqLc_LqCc[which(LqLc_LqCc$baseMean.LqLc >100 & LqLc_LqCc$baseMean.LqCc >100),]
plot(LqLc_LqCc_fil$log2FoldChange.LqLc, LqLc_LqCc_fil$log2FoldChange.LqCc,
     pch=19, cex=0.5,
     xlim=c(-5,3), ylim=c(-5,3),
     xlab='log2FC, LqLc', ylab='log2FC, LqCc', main='LqLc vs LqCc, exp filtered > 100')
abline(a=0, b=1, col='gray', lty=2)
abline(h=0)
abline(v=0)
dev.off()

# log2FC scatter plots for CqLc vs CqCc
png(filename = 'CqLc_CqCc_all.png', width=700, height=700)
CqLc_CqCc <- merge(as.data.frame(res_CqLc), as.data.frame(res_CqCc), by="row.names", suffixes=c('.CqLc','.CqCc'))
plot(CqLc_CqCc$log2FoldChange.CqLc, CqLc_CqCc$log2FoldChange.CqCc,
     pch=19, cex=0.5,
     xlim=c(-8,3), ylim=c(-8,3),
     xlab='log2FC, CqLc', ylab='log2FC, CqCc', main='CqLc vs CqCc, all genes')
abline(a=0, b=1, col='red', lty=2)
abline(h=0)
abline(v=0)
dev.off()

png(filename = 'CqLc_CqCc_padj_all.png', width=700, height=700)
plot(CqLc_CqCc$padj.CqLc, CqLc_CqCc$padj.CqCc,
     pch=19, cex=0.5, col='black',
     xlim=c(0,1), ylim=c(0,1),
     xlab='padj, CqLc', ylab='padj, CqCc', main='CqLc vs CqCc, all genes')
abline(a=0, b=1, col='red', lty=2)
abline(h=0)
abline(v=0)
dev.off()

png(filename = 'CqLc_CqCc_sig.png', width=700, height=700)
CqLc_CqCc_sig <- CqLc_CqCc[which(CqLc_CqCc$padj.CqLc <0.1 & CqLc_CqCc$padj.CqCc <0.1),]
plot(CqLc_CqCc_sig$log2FoldChange.CqLc, CqLc_CqCc_sig$log2FoldChange.CqCc,
     pch=19, cex=0.5,
     xlim=c(-8,3), ylim=c(-8,3),
     xlab='log2FC, CqLc', ylab='log2FC, CqCc', main='CqLc vs CqCc, sig genes')
abline(a=0, b=1, col='red', lty=2)
abline(h=0)
abline(v=0)
dev.off()

png(filename = 'CqLc_CqCc_fil.png', width=700, height=700)
CqLc_CqCc_fil <- CqLc_CqCc[which(CqLc_CqCc$baseMean.CqLc >100 & CqLc_CqCc$baseMean.CqCc >100),]
plot(CqLc_CqCc_fil$log2FoldChange.CqLc, CqLc_CqCc_fil$log2FoldChange.CqCc,
     pch=19, cex=0.5,
     xlim=c(-8,3), ylim=c(-8,3),
     xlab='log2FC, CqLc', ylab='log2FC, CqCc', main='CqLc vs CqCc, exp filtered > 100')
abline(a=0, b=1, col='gray', lty=2)
abline(h=0)
abline(v=0)
dev.off()

# log2FC scatter plots for LqCc vs CqCc
png(filename = 'LqCc_CqCc_all.png', width=700, height=700)
LqCc_CqCc <- merge(as.data.frame(res_LqCc), as.data.frame(res_CqCc), by="row.names", suffixes=c('.LqCc','.CqCc'))
plot(LqCc_CqCc$log2FoldChange.LqCc, LqCc_CqCc$log2FoldChange.CqCc,
     pch=19, cex=0.5,
     xlim=c(-8,3), ylim=c(-8,3),
     xlab='log2FC, LqCc', ylab='log2FC, CqCc', main='LqCc vs CqCc, all genes')
abline(a=0, b=1, col='red', lty=2)
abline(h=0)
abline(v=0)
dev.off()

png(filename = 'LqCc_CqCc_padj_all.png', width=700, height=700)
plot(LqCc_CqCc$padj.LqCc, LqCc_CqCc$padj.CqCc,
     pch=19, cex=0.5, col='black',
     xlim=c(0,1), ylim=c(0,1),
     xlab='padj, LqCc', ylab='padj, CqCc', main='LqCc vs CqCc, all genes')
abline(a=0, b=1, col='red', lty=2)
abline(h=0)
abline(v=0)
dev.off()

png(filename = 'LqCc_CqCc_sig.png', width=700, height=700)
LqCc_CqCc_sig <- LqCc_CqCc[which(LqCc_CqCc$padj.LqCc <0.1 & LqCc_CqCc$padj.CqCc <0.1),]
plot(LqCc_CqCc_sig$log2FoldChange.LqCc, LqCc_CqCc_sig$log2FoldChange.CqCc,
     pch=19, cex=0.5,
     xlim=c(-8,3), ylim=c(-8,3),
     xlab='log2FC, LqCc', ylab='log2FC, CqCc', main='LqCc vs CqCc, sig genes')
abline(a=0, b=1, col='red', lty=2)
abline(h=0)
abline(v=0)
dev.off()

png(filename = 'LqCc_CqCc_fil.png', width=700, height=700)
LqCc_CqCc_fil <- LqCc_CqCc[which(LqCc_CqCc$baseMean.LqCc >100 & LqCc_CqCc$baseMean.CqCc >100),]
plot(LqCc_CqCc_fil$log2FoldChange.LqCc, LqCc_CqCc_fil$log2FoldChange.CqCc,
     pch=19, cex=0.5,
     xlim=c(-8,3), ylim=c(-8,3),
     xlab='log2FC, LqCc', ylab='log2FC, CqCc', main='LqCc vs CqCc, exp filtered > 100')
abline(a=0, b=1, col='gray', lty=2)
abline(h=0)
abline(v=0)
dev.off()

# log2FC scatter plots for LqLc vs Frankenstein
png(filename = 'LqLc_Frankenstein_all.png', width=700, height=700)
LqLc_Frankenstein <- merge(as.data.frame(res_LqLc), as.data.frame(res_Frankenstein), by="row.names", suffixes=c('.LqLc','.Frankenstein'))
plot(LqLc_Frankenstein$log2FoldChange.LqLc, LqLc_Frankenstein$log2FoldChange.Frankenstein,
     pch=19, cex=0.5,
     xlim=c(-6,3), ylim=c(-6,3),
     xlab='log2FC, LqLc', ylab='log2FC, Frankenstein', main='LqLc vs Frankenstein, all genes')
abline(a=0, b=1, col='red', lty=2)
abline(h=0)
abline(v=0)
dev.off()

# log2FC scatter plots for testing import methods
png(filename = 'import_Frankenstein_all.png', width=700, height=700)
import_Frankenstein <- merge(as.data.frame(res_L.import), as.data.frame(res_CWG.import), by="row.names", suffixes=c('.L.import','.CWG.import'))
plot(import_Frankenstein$log2FoldChange.L.import, import_Frankenstein$log2FoldChange.CWG.import,
     pch=19, cex=0.5,
     xlim=c(-6,3), ylim=c(-6,3),
     xlab='log2FC, L.import', ylab='log2FC, CWG.import', main='import Frankenstein, all genes')
abline(a=0, b=1, col='red', lty=2)
abline(h=0)
abline(v=0)
dev.off()

png(filename = 'import_Frankenstein_all_Cq.png', width=700, height=700)
import_Frankenstein_Cq <- merge(as.data.frame(res_L.import_Cq), as.data.frame(res_CWG.import_Cq), by="row.names", suffixes=c('.L.import_Cq','.CWG.import_Cq'))
plot(import_Frankenstein_Cq$log2FoldChange.L.import_Cq, import_Frankenstein_Cq$log2FoldChange.CWG.import_Cq,
     pch=19, cex=0.5,
     xlim=c(-8,3), ylim=c(-8,3),
     xlab='log2FC, L.import_Cq', ylab='log2FC, CWG.import_Cq', main='import Frankenstein, Cq, all genes')
abline(a=0, b=1, col='red', lty=2)
abline(h=0)
abline(v=0)
dev.off()

png(filename = 'import_Frankenstein_all_Cq_EnsDb102.png', width=700, height=700)
import_Frankenstein_Cq_EnsDb102 <- merge(as.data.frame(res_CWG.import_Cq), as.data.frame(res_L.import_Cq_EnsDb102), by="row.names", suffixes=c('.CWG.import_Cq','.L.import_Cq_EnsDb102'))
plot(import_Frankenstein_Cq_EnsDb102$log2FoldChange.CWG.import_Cq, import_Frankenstein_Cq_EnsDb102$log2FoldChange.L.import_Cq_EnsDb102,
     pch=19, cex=0.5, 
     xlim=c(-8,3), ylim=c(-8,3),
     xlab='log2FC, CWG.import_Cq', ylab='log2FC, L.import_Cq_EnsDb102', main='import Frankenstein, Cq, EnsDb102, all genes')
abline(a=0, b=1, col='red', lty=2)
abline(h=0)
abline(v=0)
dev.off()
