rice_type<- read.delim("MSU7_gene_trans_type",header = F)
first600 <- read.delim("first600_step4codon",header = F)
last600 <- read.delim("last600_step4codon",header = F)
f <- merge(first600,rice_type,by.x="V1",by.y="V1")
l <- merge(last600,rice_type,by.x="V1",by.y="V1")
library("reshape2")
start <- melt(f[,-48],id=c("V1","V3.y"))
start$variable <- rep(1:46,each=23630)
library(Rmisc)
co <- summarySE(start, measurevar="value", groupvars=c("variable","V3.y"))
colnames(co)[1:2] <- c("pos","type")
co$cds <- "start"
end <- melt(l[,-48],id=c("V1","V3.y"))
end$variable <- rep(1:46,each=23630)
coe <- summarySE(end, measurevar="value", groupvars=c("variable","V3.y"))
colnames(coe)[1:2] <- c("pos","type")
coe$cds <- "end"
codon <- rbind(co,coe)
codon$cds <- factor(codon$cds,levels=c("start","end"))

library("ggplot2")
ggplot(data=codon, aes(x=pos,y=value,fill=type)) +
geom_line(aes(col=type))+geom_ribbon(aes(ymin=value-ci, ymax=value+ci),alpha=0.4)+
facet_wrap(~cds)+
scale_color_manual(values=c("black","red")) +
scale_fill_manual(values=c("black","red")) +
labs(y="Suboptimal codon frequency (%)",x="Codon position relative to CDS")+
theme( axis.text.x =element_text(size=12))+
theme(legend.position='none')+
scale_x_continuous(breaks=seq(1,46,by=11.25),labels=c("Start","50","100","150","200"))+
theme(strip.background = element_rect(fill = "white"),strip.text = element_text(colour ="white",size=12),
axis.text =element_text(colour ="black",size=12),
axis.title = element_text(size=12))
