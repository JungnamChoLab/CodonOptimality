###calculation of codon usage of TEs###
setwd("~/gccon")
library(ade4) #install.packages("ade4")
library(seqinr) #install.packages("seqinr")
library(ggplot2)

#Fly
fly_g <- read.delim("fly_gene.cds_longest_len",header = F)
fly_t <- read.delim("fly_transdecoder.longest.cds_len",header=F)
geneorf <- read.fasta("Drosophila_melanogaster.BDGP6.22.cds.all.fa.gz")
orf <- read.fasta("fly_dm6_ucsc_.fa.transdecoder.cds", set.attributes = F)
geneorf_df <- GetCodonFrequency(geneorf)
geneorf_df_t <- data.frame(t(geneorf_df))
geneorf_df_t$sum <- apply(geneorf_df_t[,1:64],1,sum)
geneorf_df_t$gc1<-apply(geneorf_df_t[,17:48],1,sum)
geneorf_df_t$GC1<-geneorf_df_t$gc1/geneorf_df_t$sum*100
geneorf_df_t$gc2<-apply(geneorf_df_t[,c(5:12,21:28,37:44,53:60)],1,sum)
geneorf_df_t$GC2<-geneorf_df_t$gc2/geneorf_df_t$sum*100
geneorf_df_t$gc3<-apply(geneorf_df_t[,c(2:3,6:7,10:11,14:15,18:19,22:23,26:27,30:31,34:35,38:39,42:43,46:47,50:51,54:55,58:59,62:63)],1,sum)
geneorf_df_t$GC3<-geneorf_df_t$gc3/geneorf_df_t$sum*100

geneorf_df_t_l <- geneorf_df_t[rownames(geneorf_df_t) %in% fly_g$V2,]
ggc <- stack(geneorf_df_t[,c(67,69,71)])
ggc <- stack(geneorf_df_t_l[,c(67,69,71)])
ggc$type <- "Gene"
ggplot(data=ggc, aes(x=ind,y=values))+geom_violin(aes(fill=ind))+labs(y="GC content (%)",x="")+theme( axis.text.x =element_text(size=12))

orf_df <- GetCodonFrequency(orf)
orf_df_t <- data.frame(t(orf_df))
orf_df_t$sum <- apply(orf_df_t[,1:64],1,sum)
orf_df_t$gc1<-apply(orf_df_t[,17:48],1,sum)
orf_df_t$GC1<-orf_df_t$gc1/orf_df_t$sum*100
orf_df_t$gc2<-apply(orf_df_t[,c(5:12,21:28,37:44,53:60)],1,sum)
orf_df_t$GC2<-orf_df_t$gc2/orf_df_t$sum*100
orf_df_t$gc3<-apply(orf_df_t[,c(2:3,6:7,10:11,14:15,18:19,22:23,26:27,30:31,34:35,38:39,42:43,46:47,50:51,54:55,58:59,62:63)],1,sum)
orf_df_t$GC3<-orf_df_t$gc3/orf_df_t$sum*100

orf_df_t_l <- orf_df_t[rownames(orf_df_t) %in% fly_t$V2,]

#gc <- stack(orf_df_t[,c(67,69,71)]) #not filter the longest transposon ORF
gc <- stack(orf_df_t_l[,c(67,69,71)])
gc$type <- "Transposon"

g <- rbind(gc,ggc)
ggplot(data=g, aes(x=ind,y=values))+geom_violin(draw_quantiles = TRUE,aes(fill=ind))+labs(y="GC content (%)",x="D. melanogaster")+
  theme(strip.background = element_rect(fill = "white"),strip.text = element_text(colour ="black",size=12),
        axis.text =element_text(colour ="black",size=12),
        axis.title = element_text(size=12),axis.title.x = element_text(face = "italic"))+facet_wrap(~type)+
  theme(legend.position='none')



#zebrafish
fish_g <- read.delim("zebrafish_gene.cds_longest_len",header = F)
fish_t <- read.delim("zebrafish_transdecoder.longest.cds_len",header=F)
geneorf <- read.fasta("Danio_rerio.GRCz11.cds.all.fa")
orf <- read.fasta("zebrafish_danRer11_ucsc_.fa.transdecoder.cds", set.attributes = F)
geneorf_df <- GetCodonFrequency(geneorf)
igeneorf_df_t <- data.frame(t(geneorf_df))
igeneorf_df_t$sum <- apply(igeneorf_df_t[,1:64],1,sum)
igeneorf_df_t$gc1<-apply(igeneorf_df_t[,17:48],1,sum)
igeneorf_df_t$GC1<-igeneorf_df_t$gc1/igeneorf_df_t$sum*100
igeneorf_df_t$gc2<-apply(igeneorf_df_t[,c(5:12,21:28,37:44,53:60)],1,sum)
igeneorf_df_t$GC2<-igeneorf_df_t$gc2/igeneorf_df_t$sum*100
igeneorf_df_t$gc3<-apply(igeneorf_df_t[,c(2:3,6:7,10:11,14:15,18:19,22:23,26:27,30:31,34:35,38:39,42:43,46:47,50:51,54:55,58:59,62:63)],1,sum)
igeneorf_df_t$GC3<-igeneorf_df_t$gc3/igeneorf_df_t$sum*100

geneorf_df_t_l <- igeneorf_df_t[rownames(igeneorf_df_t) %in% fish_g$V2,]
ggc <- stack(igeneorf_df_t[,c(67,69,71)])
ggc <- stack(geneorf_df_t_l[,c(67,69,71)])
ggc$type <- "Gene"
ggplot(data=ggc, aes(x=ind,y=values))+geom_violin(aes(fill=ind))+labs(y="GC content (%)",x="")+theme( axis.text.x =element_text(size=12))

orf_df <- GetCodonFrequency(orf)
iorf_df_t <- data.frame(t(orf_df))
iorf_df_t$sum <- apply(iorf_df_t[,1:64],1,sum)
iorf_df_t$gc1<-apply(iorf_df_t[,17:48],1,sum)
iorf_df_t$GC1<-iorf_df_t$gc1/iorf_df_t$sum*100
iorf_df_t$gc2<-apply(iorf_df_t[,c(5:12,21:28,37:44,53:60)],1,sum)
iorf_df_t$GC2<-iorf_df_t$gc2/iorf_df_t$sum*100
iorf_df_t$gc3<-apply(iorf_df_t[,c(2:3,6:7,10:11,14:15,18:19,22:23,26:27,30:31,34:35,38:39,42:43,46:47,50:51,54:55,58:59,62:63)],1,sum)
iorf_df_t$GC3<-iorf_df_t$gc3/iorf_df_t$sum*100

orf_df_t_l <- iorf_df_t[rownames(iorf_df_t) %in% fish_t$V2,]

#gc <- stack(iorf_df_t[,c(67,69,71)]) #not filter the longest transposon ORF
gc <- stack(orf_df_t_l[,c(67,69,71)])
gc$type <- "Transposon"

g <- rbind(gc,ggc)
ggplot(data=g, aes(x=ind,y=values))+geom_violin(draw_quantiles = TRUE,aes(fill=ind))+labs(y="GC content (%)",x="Danio rerio")+
  theme(strip.background = element_rect(fill = "white"),strip.text = element_text(colour ="black",size=12),
        axis.text =element_text(colour ="black",size=12),
        axis.title = element_text(size=12),axis.title.x = element_text(face = "italic"))+facet_wrap(~type)+
  theme(legend.position='none')


#human
human_g <- read.delim("human_hg38_gene.cds_longest_len",header = F)
human_t <- read.delim("human_transdecoder.longest.cds_len",header=F)
geneorf <- read.fasta("Homo_sapiens.GRCh38.cds.all.fa")
orf <- read.fasta("human_hg38_ucsc_.fa.transdecoder.cds", set.attributes = F)
geneorf_df <- GetCodonFrequency(geneorf)
hgeneorf_df_t <- data.frame(t(geneorf_df))
hgeneorf_df_t$sum <- apply(hgeneorf_df_t[,1:64],1,sum)
hgeneorf_df_t$gc1<-apply(hgeneorf_df_t[,17:48],1,sum)
hgeneorf_df_t$GC1<-hgeneorf_df_t$gc1/hgeneorf_df_t$sum*100
hgeneorf_df_t$gc2<-apply(hgeneorf_df_t[,c(5:12,21:28,37:44,53:60)],1,sum)
hgeneorf_df_t$GC2<-hgeneorf_df_t$gc2/hgeneorf_df_t$sum*100
hgeneorf_df_t$gc3<-apply(hgeneorf_df_t[,c(2:3,6:7,10:11,14:15,18:19,22:23,26:27,30:31,34:35,38:39,42:43,46:47,50:51,54:55,58:59,62:63)],1,sum)
hgeneorf_df_t$GC3<-hgeneorf_df_t$gc3/hgeneorf_df_t$sum*100

geneorf_df_t_l <- hgeneorf_df_t[rownames(hgeneorf_df_t) %in% human_g$V2,]
ggc <- stack(mgeneorf_df_t[,c(67,69,71)])
ggc <- stack(geneorf_df_t_l[,c(67,69,71)])
ggc$type <- "Gene"
ggplot(data=ggc, aes(x=ind,y=values))+geom_violin(aes(fill=ind))+labs(y="GC content (%)",x="")+theme( axis.text.x =element_text(size=12))

geneorf_df_t_l <- hgeneorf_df_t[rownames(hgeneorf_df_t) %in% human_g$V2,]
geneorf_df_t_l$cdslen <- geneorf_df_t_l$sum*3
cor.test(geneorf_df_t_l$cdslen,geneorf_df_t_l$GC3)
plot(geneorf_df_t_l$cdslen,geneorf_df_t_l$GC3,pch=20,ylab="GC3 content",xlab="CDS length",main="Homo_sapiens",xlim=c(0,15000))

orf_df <- GetCodonFrequency(orf)
horf_df_t <- data.frame(t(orf_df))
horf_df_t$sum <- apply(horf_df_t[,1:64],1,sum)
horf_df_t$gc1<-apply(horf_df_t[,17:48],1,sum)
horf_df_t$GC1<-horf_df_t$gc1/horf_df_t$sum*100
horf_df_t$gc2<-apply(horf_df_t[,c(5:12,21:28,37:44,53:60)],1,sum)
horf_df_t$GC2<-horf_df_t$gc2/horf_df_t$sum*100
horf_df_t$gc3<-apply(horf_df_t[,c(2:3,6:7,10:11,14:15,18:19,22:23,26:27,30:31,34:35,38:39,42:43,46:47,50:51,54:55,58:59,62:63)],1,sum)
horf_df_t$GC3<-horf_df_t$gc3/horf_df_t$sum*100

orf_df_t_l <- horf_df_t[rownames(horf_df_t) %in% human_t$V2,]

#gc <- stack(horf_df_t[,c(67,69,71)]) #not filter the longest transposon ORF
gc <- stack(orf_df_t_l[,c(67,69,71)])
gc$type <- "Transposon"

g <- rbind(gc,ggc)
ggplot(data=g, aes(x=ind,y=values))+geom_violin(draw_quantiles = TRUE,aes(fill=ind))+labs(y="GC content (%)",x="Homo sapiens")+
  theme(strip.background = element_rect(fill = "white"),strip.text = element_text(colour ="black",size=12),
        axis.text =element_text(colour ="black",size=12),
        axis.title = element_text(size=12),axis.title.x = element_text(face = "italic"))+facet_wrap(~type)+
  theme(legend.position='none')

#mouse
mouse_g <- read.delim("mouse_mm10_gene.cds_longest_len",header = F)
mouse_t <- read.delim("mouse_transdecoder.longest.cds_len",header=F)
geneorf <- read.fasta("Mus_musculus.GRCm38.cds.all.fa")
orf <- read.fasta("mouse_mm10_ucsc_.fa.transdecoder.cds", set.attributes = F)
geneorf_df <- GetCodonFrequency(geneorf)
mgeneorf_df_t <- data.frame(t(geneorf_df))
mgeneorf_df_t$sum <- apply(mgeneorf_df_t[,1:64],1,sum)
mgeneorf_df_t$gc1<-apply(mgeneorf_df_t[,17:48],1,sum)
mgeneorf_df_t$GC1<-mgeneorf_df_t$gc1/mgeneorf_df_t$sum*100
mgeneorf_df_t$gc2<-apply(mgeneorf_df_t[,c(5:12,21:28,37:44,53:60)],1,sum)
mgeneorf_df_t$GC2<-mgeneorf_df_t$gc2/mgeneorf_df_t$sum*100
mgeneorf_df_t$gc3<-apply(mgeneorf_df_t[,c(2:3,6:7,10:11,14:15,18:19,22:23,26:27,30:31,34:35,38:39,42:43,46:47,50:51,54:55,58:59,62:63)],1,sum)
mgeneorf_df_t$GC3<-mgeneorf_df_t$gc3/mgeneorf_df_t$sum*100

mgeneorf_df_t_l <- mgeneorf_df_t[rownames(mgeneorf_df_t) %in% mouse_g$V2,]
mgeneorf_df_t_l$cdslen <- mouse_g[mouse_g$V2 %in% rownames(mgeneorf_df_t),]$V3
cor.test(mgeneorf_df_t_l$cdslen,mgeneorf_df_t_l$GC3)
plot(mgeneorf_df_t_l$cdslen,mgeneorf_df_t_l$GC3,pch=20,xlim=c(0,20000))

geneorf_df_t_l <- mgeneorf_df_t[rownames(mgeneorf_df_t) %in% mouse_g$V2,]

#ggc <- stack(mgeneorf_df_t[,c(67,69,71)])
ggc <- stack(geneorf_df_t_l[,c(67,69,71)])
ggc$type <- "Gene"
ggplot(data=ggc, aes(x=ind,y=values))+geom_violin(aes(fill=ind))+labs(y="GC content (%)",x="")+theme( axis.text.x =element_text(size=12))

orf_df <- GetCodonFrequency(orf)
morf_df_t <- data.frame(t(orf_df))
morf_df_t$sum <- apply(morf_df_t[,1:64],1,sum)
morf_df_t$gc1<-apply(morf_df_t[,17:48],1,sum)
morf_df_t$GC1<-morf_df_t$gc1/morf_df_t$sum*100
morf_df_t$gc2<-apply(morf_df_t[,c(5:12,21:28,37:44,53:60)],1,sum)
morf_df_t$GC2<-morf_df_t$gc2/morf_df_t$sum*100
morf_df_t$gc3<-apply(morf_df_t[,c(2:3,6:7,10:11,14:15,18:19,22:23,26:27,30:31,34:35,38:39,42:43,46:47,50:51,54:55,58:59,62:63)],1,sum)
morf_df_t$GC3<-morf_df_t$gc3/morf_df_t$sum*100

orf_df_t_l <- morf_df_t[rownames(morf_df_t) %in% mouse_t$V2,]

#gc <- stack(morf_df_t[,c(67,69,71)]) #not filter the longest transposon ORF
gc <- stack(orf_df_t_l[,c(67,69,71)])
gc$type <- "Transposon"

g <- rbind(gc,ggc)
ggplot(data=g, aes(x=ind,y=values))+geom_violin(draw_quantiles = TRUE,aes(fill=ind))+labs(y="GC content (%)",x="Mus musculus")+
  theme(strip.background = element_rect(fill = "white"),strip.text = element_text(colour ="black",size=12),
        axis.text =element_text(colour ="black",size=12),
        axis.title = element_text(size=12),axis.title.x = element_text(face = "italic"))+facet_wrap(~type)+
  theme(legend.position='none')





