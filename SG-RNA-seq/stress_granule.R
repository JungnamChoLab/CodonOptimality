# Use the raw read count got from featurecounts to do DEG analysis and got SG-enriched and depleted transcripts
g <-read.delim("stress_granule_raw_read_count")
colnames(g)[7:10] <- c("Total1","Total2","SG1","SG2")
gg <- g[,c(1,7:10)]
g <- merge(t_type,gg,by.x="V1",by.y="Geneid")
table(g$V2)
library("DESeq2")
gg <- data.frame(g[,c(3:6)],row.names=g$V1)
countData <- as.matrix(gg)
colData <- data.frame("condition" = rep(c("Total","SG"),each=2))
rownames(colData) <- colnames(countData)
colData
all(rownames(colData) %in% colnames(countData))
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~condition)
dds$condition
dds$condition <- factor(dds$condition, levels = c("Total","SG"))
dds$condition
dds <- DESeq(dds)
res <- results(dds)
res
hist( res$pvalue, breaks=20, col="grey" )
plotDispEsts( dds , ylim = c(1e-6, 1e1))
summary(res)
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds)

#use ggmaplot draw the MA figure https://rpkgs.datanovia.com/ggpubr/reference/ggmaplot.html
install.packages("ggpubr")
library("ggpubr")
ggmaplot(res, fdr = 0.05, fc = 4, genenames = NULL,
                               detection_call = NULL, size = 0.4, font.label = NULL, label.rectangle = FALSE, palette = c("#B31B21", "#1465AC",
                            "darkgray"), top = 0, select.top.method = c("padj", "fc"),
                               main = NULL, xlab = "Log2 Mean expression",
                               ylab = "Log2 Fold change (SG/Total)", ggtheme = theme_classic()+theme(legend.position='none'))
#export img in eps format w*h 600*400

sum(res$padj < 0.05, na.rm=TRUE)
re <- as.data.frame(res)
re$Significant  <- "No"
re[is.na(res$padj),]$Significant  <- "Not test"
#re[re$padj < 0.05 &abs(re$log2FoldChange)>=1,]$sign <-"yes"
re[re$padj <= 0.05 &re$log2FoldChange>=1 &!is.na(res$padj),]$Significant <-"Up"
re[re$padj <= 0.05 &re$log2FoldChange<=-1&!is.na(res$padj),]$Significant  <-"Down"
re$Significant  <- factor(re$Significant,levels = c("Up","Down","No"))
table(re[,c(7)])
v <- ggplot(re,aes(log2FoldChange,-log10(padj)))
v + geom_point(aes(colour=Significant))+geom_vline(xintercept=c(-1,1),linetype=4,colour="grey")+geom_hline(yintercept=-log10(0.05),linetype=4,colour="grey")+xlim(-17, 17)
write.table(re,"allDEGs.txt",sep="\t",quote=F)
