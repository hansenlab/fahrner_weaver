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
