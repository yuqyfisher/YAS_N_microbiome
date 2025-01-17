####Fig.3
#R4.2.3
library(ggplot2)
library("phyloseq")
#install.packages("vegan")
library("ggplot2")
library("scales")
library("grid")
theme_set(theme_bw())
library("DESeq2")
library("ape")
library("vegan")
library("data.table")
library("RColorBrewer")
library(colorRamps)
library("svglite")
library(VennDiagram)


#Fig.3a-f
#1
data<-read.table("ACI.txt",header=TRUE,sep="\t")
combined =data##用于计算显著性

AN_UP <- data[,c(1,2,4)]
ratio <- data.frame(cbind(AN_UP,AN_UP$AN_up/AN_UP$AN_total))
dat <- ratio[,c(1,4)]
colnames(dat) <-c("Annotation","AN_UP_ratio")
merAN<- dat


AU_UP <- data[,c(1,3,7)]
ratio2 <- data.frame(cbind(AU_UP,AU_UP$AU_up/AU_UP$AU_total))
dat2 <- ratio2[,c(1,4)]
colnames(dat2) <- c("Annotation","AU_UP_ratio")
merAU <- dat2

rm <- merge(merAN,merAU,by="Annotation")
rm.m <- melt(rm)
rm.m$variable <- factor(rm.m$variable, levels=c("AN_UP_ratio","AU_UP_ratio"))
dat <- rm.m[order(-rm.m[,3]),]
dat$Annotation <- factor(dat$Annotation, levels=c(
  "TSL",
  "GBM (S)",
  "CM (Y)",
  "BC (A)",
  "CC (A)",
  "TSP",
  "AAM (Y)",
  "OAAM (Y)",
  "CGD",
  "SM (A)",
  "OSM (S)",
  "NM (A)",
  "FSD",
  "ACD",
  "RR (S)",
  "MT (A)",
  "QS",
  "MMDE",
  "EB (S)",
  "Carbon M (Y)",
  "EM (A)",
  "AAB (Y)",
  "BSM (S)",
  "LM (Y)",
  "CFPP (A)",
  "TPM",
  "CB",
  "CVM",
  "Nucleotide M (Y)",
  "ST",
  "NSB (Y)",
  "BF (S)",
  "FAB (Y)",
  "FAM (Y)",
  "LB (S)"
))

samples_new <- sapply(as.character(dat$variable), function(x) strsplit(x, "[_]"))

Timepoint <- {}
for(i in 1:70){
  T <- samples_new[[i]][1]
  Timepoint <- append(Timepoint,T)
}

SampleType <- {}
for(i in 1:70){
  T <- samples_new[[i]][2]
  SampleType <- append(SampleType,T)
}

dat$SampleType = SampleType
dat$Timepoint = Timepoint


combined <- read.table("combined.ACI.txt",header=TRUE,sep="\t")#用于计算显著性
combined$AN_UP_pvalue <- phyper(combined$AN_up-1,combined$AN_total,sum(combined$AN_total)-combined$AN_total,sum(combined$AN_up),
                                lower.tail = FALSE, log.p = FALSE)

combined$AU_UP_pvalue <- phyper(combined$AU_up-1,combined$AU_total,sum(combined$AU_total)-combined$AU_total,sum(combined$AU_up),
                                lower.tail = FALSE, log.p = FALSE)

combined[which(combined$AN_UP_pvalue < 0.05),'sig_AN'] <- 'Significant'
combined[which(combined$AU_UP_pvalue < 0.05),'sig_AU'] <- 'Significant'


dat[which(dat$value <= 1 ),'sig'] <- 'Non-significant'
dat[which(dat$value > 1),'sig'] <- 'Significant'
dat1 <- subset(dat,Timepoint=="AN50")
dat2 <- subset(dat,Timepoint=="AU50")

ggplot(dat, aes(x=value, y=Annotation,colour=Timepoint,shape=sig)) +
  scale_colour_manual(name="",values = c("royalblue1","#F24F50","#9B4F0F","#6FB98F"))+
  scale_shape_manual(name="",values=c(1,1,19,19))+facet_wrap(~SampleType,scales="free_x",ncol=2)+
  scale_shape_manual(values  = c(21,21), limits = c('Non-significant', 'Significant')) +
  geom_point(size = 4,stroke = 1)+theme_bw()+
  geom_vline(aes(xintercept=1),colour="Black",size=1,linetype="dashed")+xlab("Enrichment ratio")+
  theme(axis.text.x=element_text(size=11,color="black",face="bold",angle=90), 
        axis.text.y=element_text(size=11,color="black",face="bold"), 
        axis.title=element_text(size=11,face="bold"),text=element_text(size=10))+
  xlim(0, 4.6)  +
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) +
  theme(strip.text = element_text(size=11,face="bold"))

#2
data<-read.table("ACTINO.txt",header=TRUE,sep="\t")
combined =data##用于计算显著性

AN_UP <- data[,c(1,2,4)]
ratio <- data.frame(cbind(AN_UP,AN_UP$AN_up/AN_UP$AN_total))
dat <- ratio[,c(1,4)]
colnames(dat) <-c("Annotation","AN_UP_ratio")
merAN<- dat


AU_UP <- data[,c(1,3,7)]
ratio2 <- data.frame(cbind(AU_UP,AU_UP$AU_up/AU_UP$AU_total))
dat2 <- ratio2[,c(1,4)]
colnames(dat2) <- c("Annotation","AU_UP_ratio")
merAU <- dat2



rm <- merge(merAN,merAU,by="Annotation")
rm.m <- melt(rm)
rm.m$variable <- factor(rm.m$variable, levels=c("AN_UP_ratio","AU_UP_ratio"))
dat <- rm.m[order(-rm.m[,3]),]
dat$Annotation <- factor(dat$Annotation, levels=c(
  "MMDE",
  "TPM",
  "OSM (S)",
  "MT (A)",
  "AAM (Y)",
  "CFPP (A)",
  "CM (Y)",
  "CC (A)",
  "BC (A)",
  "LB (S)",
  "CGD",
  "BF (S)",
  "ST",
  "NSB (Y)",
  "EB (S)",
  "SM (A)",
  "FAB (Y)",
  "GBM (S)",
  "TSL",
  "TSP",
  "Nucleotide M (Y)",
  "EM (A)",
  "ACD",
  "QS",
  "CVM",
  "LM (Y)",
  "FAM (Y)",
  "CB",
  "OAAM (Y)",
  "BSM (S)",
  "RR (S)",
  "Carbon M (Y)",
  "NM (A)",
  "AAB (Y)",
  "FSD" 
))

samples_new <- sapply(as.character(dat$variable), function(x) strsplit(x, "[_]"))

Timepoint <- {}
for(i in 1:70){
  T <- samples_new[[i]][1]
  Timepoint <- append(Timepoint,T)
}

SampleType <- {}
for(i in 1:70){
  T <- samples_new[[i]][2]
  SampleType <- append(SampleType,T)
}

dat$SampleType = SampleType
dat$Timepoint = Timepoint


combined <- read.table("combined.ACTINO.txt",header=TRUE,sep="\t")#用于计算显著性
combined$AN_UP_pvalue <- phyper(combined$AN_up-1,combined$AN_total,sum(combined$AN_total)-combined$AN_total,sum(combined$AN_up),
                                lower.tail = FALSE, log.p = FALSE)

combined$AU_UP_pvalue <- phyper(combined$AU_up-1,combined$AU_total,sum(combined$AU_total)-combined$AU_total,sum(combined$AU_up),
                                lower.tail = FALSE, log.p = FALSE)

combined[which(combined$AN_UP_pvalue < 0.05),'sig_AN'] <- 'Significant'
combined[which(combined$AU_UP_pvalue < 0.05),'sig_AU'] <- 'Significant'


dat[which(dat$value <= 1 ),'sig'] <- 'Non-significant'
dat[which(dat$value > 1),'sig'] <- 'Significant'
dat1 <- subset(dat,Timepoint=="AN50")
dat2 <- subset(dat,Timepoint=="AU50")



ggplot(dat, aes(x=value, y=Annotation,colour=Timepoint,shape=sig)) +
  scale_colour_manual(name="",values = c("royalblue1","#F24F50","#9B4F0F","#6FB98F"))+
  scale_shape_manual(name="",values=c(1,1,19,19))+facet_wrap(~SampleType,scales="free_x",ncol=2)+
  scale_shape_manual(values  = c(21,21), limits = c('Non-significant', 'Significant')) +
  geom_point(size = 4,stroke = 1)+theme_bw()+
  geom_vline(aes(xintercept=1),colour="Black",size=1,linetype="dashed")+xlab("Enrichment ratio")+
  theme(axis.text.x=element_text(size=11,color="black",face="bold",angle=90), 
        axis.text.y=element_text(size=11,color="black",face="bold"), 
        axis.title=element_text(size=11,face="bold"),text=element_text(size=10))+
  xlim(0, 4.6)  +
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) +
  theme(strip.text = element_text(size=11,face="bold"))



#3
data<-read.table("BAC.txt",header=TRUE,sep="\t")
combined =data##用于计算显著性

AN_UP <- data[,c(1,2,4)]
ratio <- data.frame(cbind(AN_UP,AN_UP$AN_up/AN_UP$AN_total))
dat <- ratio[,c(1,4)]
colnames(dat) <-c("Annotation","AN_UP_ratio")
merAN<- dat


AU_UP <- data[,c(1,3,7)]
ratio2 <- data.frame(cbind(AU_UP,AU_UP$AU_up/AU_UP$AU_total))
dat2 <- ratio2[,c(1,4)]
colnames(dat2) <- c("Annotation","AU_UP_ratio")
merAU <- dat2


rm <- merge(merAN,merAU,by="Annotation")
rm.m <- melt(rm)
rm.m$variable <- factor(rm.m$variable, levels=c("AN_UP_ratio","AU_UP_ratio"))
dat <- rm.m[order(-rm.m[,3]),]
dat$Annotation <- factor(dat$Annotation, levels=c(
  "CVM",
  "BSM (S)",
  "AAB (Y)",
  "CB",
  "TSL",
  "ST",
  "ACD",
  "FSD",
  "MT (A)",
  "MMDE",
  "BF (S)",
  "EB (S)",
  "CC (A)",
  "TPM",
  "EM (A)",
  "OAAM (Y)",
  "Carbon M (Y)",
  "CFPP (A)",
  "FAM (Y)",
  "AAM (Y)",
  "LM (Y)",
  "QS",
  "OSM (S)",
  "NM (A)",
  "SM (A)",
  "RR (S)",
  "TSP",
  "CGD",
  "FAB (Y)",
  "BC (A)",
  "CM (Y)",
  "GBM (S)",
  "Nucleotide M (Y)",
  "NSB (Y)",
  "LB (S)"
))

samples_new <- sapply(as.character(dat$variable), function(x) strsplit(x, "[_]"))

Timepoint <- {}
for(i in 1:70){
  T <- samples_new[[i]][1]
  Timepoint <- append(Timepoint,T)
}

SampleType <- {}
for(i in 1:70){
  T <- samples_new[[i]][2]
  SampleType <- append(SampleType,T)
}

dat$SampleType = SampleType
dat$Timepoint = Timepoint

combined <- read.table("combined.BAC.txt",header=TRUE,sep="\t")#用于计算显著性

combined$AN_UP_pvalue <- phyper(combined$AN_up-1,combined$AN_total,sum(combined$AN_total)-combined$AN_total,sum(combined$AN_up),
                                lower.tail = FALSE, log.p = FALSE)

combined$AU_UP_pvalue <- phyper(combined$AU_up-1,combined$AU_total,sum(combined$AU_total)-combined$AU_total,sum(combined$AU_up),
                                lower.tail = FALSE, log.p = FALSE)

combined[which(combined$AN_UP_pvalue < 0.05),'sig_AN'] <- 'Significant'
combined[which(combined$AU_UP_pvalue < 0.05),'sig_AU'] <- 'Significant'



dat[which(dat$value <= 1 ),'sig'] <- 'Non-significant'
dat[which(dat$value > 1),'sig'] <- 'Significant'
dat1 <- subset(dat,Timepoint=="AN50")
dat2 <- subset(dat,Timepoint=="AU50")
dat3 <- subset(dat,SampleType=="UP")



ggplot(dat, aes(x=value, y=Annotation,colour=Timepoint,shape=sig)) +
  scale_colour_manual(name="",values = c("royalblue1","#F24F50","#9B4F0F","#6FB98F"))+
  scale_shape_manual(name="",values=c(1,1,19,19))+facet_wrap(~SampleType,scales="free_x",ncol=2)+
  scale_shape_manual(values  = c(21,21), limits = c('Non-significant', 'Significant')) +
  geom_point(size = 4,stroke = 1)+theme_bw()+
  geom_vline(aes(xintercept=1),colour="Black",size=1,linetype="dashed")+xlab("Enrichment ratio")+
  theme(axis.text.x=element_text(size=11,color="black",face="bold",angle=90), 
        axis.text.y=element_text(size=11,color="black",face="bold"), 
        axis.title=element_text(size=11,face="bold"),text=element_text(size=10))+
  xlim(0, 4.6)  +
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) +
  theme(strip.text = element_text(size=11,face="bold"))


#4
data<-read.table("GEMM.txt",header=TRUE,sep="\t")
combined =data##用于计算显著性

AN_UP <- data[,c(1,2,4)]
ratio <- data.frame(cbind(AN_UP,AN_UP$AN_up/AN_UP$AN_total))
dat <- ratio[,c(1,4)]
colnames(dat) <-c("Annotation","AN_UP_ratio")
merAN<- dat


AU_UP <- data[,c(1,3,7)]
ratio2 <- data.frame(cbind(AU_UP,AU_UP$AU_up/AU_UP$AU_total))
dat2 <- ratio2[,c(1,4)]
colnames(dat2) <- c("Annotation","AU_UP_ratio")
merAU <- dat2


rm <- merge(merAN,merAU,by="Annotation")
rm.m <- melt(rm)
rm.m$variable <- factor(rm.m$variable, levels=c("AN_UP_ratio","AU_UP_ratio"))
dat <- rm.m[order(-rm.m[,3]),]
dat$Annotation <- factor(dat$Annotation, levels=c(
  "ST",
  "CB",
  "Nucleotide M (Y)",
  "FSD",
  "FAM (Y)",
  "TSP",
  "NSB (Y)",
  "CC (A)",
  "ACD",
  "TPM",
  "TSL",
  "QS",
  "RR (S)",
  "MMDE",
  "EB (S)",
  "SM (A)",
  "MT (A)",
  "OSM (S)",
  "CM (Y)",
  "BSM (S)",
  "EM (A)",
  "LM (Y)",
  "Carbon M (Y)",
  "CVM",
  "OAAM (Y)",
  "AAB (Y)",
  "NM (A)",
  "CFPP (A)",
  "BF (S)",
  "FAB (Y)",
  "CGD",
  "AAM (Y)",
  "GBM (S)",
  "LB (S)",
  "BC (A)"
))


samples_new <- sapply(as.character(dat$variable), function(x) strsplit(x, "[_]"))

Timepoint <- {}
for(i in 1:70){
  T <- samples_new[[i]][1]
  Timepoint <- append(Timepoint,T)
}

SampleType <- {}
for(i in 1:70){
  T <- samples_new[[i]][2]
  SampleType <- append(SampleType,T)
}

dat$SampleType = SampleType
dat$Timepoint = Timepoint

combined <- read.table("combined.GEMM.txt",header=TRUE,sep="\t")#用于计算显著性

combined$AN_UP_pvalue <- phyper(combined$AN_up-1,combined$AN_total,sum(combined$AN_total)-combined$AN_total,sum(combined$AN_up),
                                lower.tail = FALSE, log.p = FALSE)


combined$AU_UP_pvalue <- phyper(combined$AU_up-1,combined$AU_total,sum(combined$AU_total)-combined$AU_total,sum(combined$AU_up),
                                lower.tail = FALSE, log.p = FALSE)

combined[which(combined$AN_UP_pvalue < 0.05),'sig_AN'] <- 'Significant'
combined[which(combined$AU_UP_pvalue < 0.05),'sig_AU'] <- 'Significant'


dat[which(dat$value <= 1 ),'sig'] <- 'Non-significant'
dat[which(dat$value > 1),'sig'] <- 'Significant'
dat1 <- subset(dat,Timepoint=="AN50")
dat2 <- subset(dat,Timepoint=="AU50")
dat3 <- subset(dat,SampleType=="UP")


ggplot(dat, aes(x=value, y=Annotation,colour=Timepoint,shape=sig)) +
  scale_colour_manual(name="",values = c("royalblue1","#F24F50","#9B4F0F","#6FB98F"))+
  scale_shape_manual(name="",values=c(1,1,19,19))+facet_wrap(~SampleType,scales="free_x",ncol=2)+
  scale_shape_manual(values  = c(21,21), limits = c('Non-significant', 'Significant')) +
  geom_point(size = 4,stroke = 1)+theme_bw()+
  geom_vline(aes(xintercept=1),colour="Black",size=1,linetype="dashed")+xlab("Enrichment ratio")+
  theme(axis.text.x=element_text(size=11,color="black",face="bold",angle=90), 
        axis.text.y=element_text(size=11,color="black",face="bold"), 
        axis.title=element_text(size=11,face="bold"),text=element_text(size=10))+
  xlim(0, 4.6)  +
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) +
  theme(strip.text = element_text(size=11,face="bold"))



#5
data<-read.table("PLAN.txt",header=TRUE,sep="\t")
combined =data##用于计算显著性

AN_UP <- data[,c(1,2,4)]
ratio <- data.frame(cbind(AN_UP,AN_UP$AN_up/AN_UP$AN_total))
dat <- ratio[,c(1,4)]
colnames(dat) <-c("Annotation","AN_UP_ratio")
merAN<- dat


AU_UP <- data[,c(1,3,7)]
ratio2 <- data.frame(cbind(AU_UP,AU_UP$AU_up/AU_UP$AU_total))
dat2 <- ratio2[,c(1,4)]
colnames(dat2) <- c("Annotation","AU_UP_ratio")
merAU <- dat2


rm <- merge(merAN,merAU,by="Annotation")
rm.m <- melt(rm)
rm.m$variable <- factor(rm.m$variable, levels=c("AN_UP_ratio","AU_UP_ratio"))
dat <- rm.m[order(-rm.m[,3]),]
dat$Annotation <- factor(dat$Annotation, levels=c(
  "CM (Y)",
  "SM (A)",
  "LB (S)",
  "TSL",
  "ACD",
  "RR (S)",
  "EB (S)",
  "QS",
  "Nucleotide M (Y)",
  "OAAM (Y)",
  "TPM",
  "MMDE",
  "AAM (Y)",
  "Carbon M (Y)",
  "BSM (S)",
  "FAB (Y)",
  "FSD",
  "CVM",
  "CFPP (A)",
  "AAB (Y)",
  "CB",
  "MT (A)",
  "CGD",
  "EM (A)",
  "NM (A)",
  "FAM (Y)",
  "OSM (S)",
  "GBM (S)",
  "LM (Y)",
  "NSB (Y)",
  "TSP",
  "ST",
  "BF (S)",
  "BC (A)",
  "CC (A)"
))


samples_new <- sapply(as.character(dat$variable), function(x) strsplit(x, "[_]"))

Timepoint <- {}
for(i in 1:70){
  T <- samples_new[[i]][1]
  Timepoint <- append(Timepoint,T)
}

SampleType <- {}
for(i in 1:70){
  T <- samples_new[[i]][2]
  SampleType <- append(SampleType,T)
}

dat$SampleType = SampleType
dat$Timepoint = Timepoint

combined <- read.table("combined.PLAN.txt",header=TRUE,sep="\t")#用于计算显著性

combined$AN_UP_pvalue <- phyper(combined$AN_up-1,combined$AN_total,sum(combined$AN_total)-combined$AN_total,sum(combined$AN_up),
                                lower.tail = FALSE, log.p = FALSE)


combined$AU_UP_pvalue <- phyper(combined$AU_up-1,combined$AU_total,sum(combined$AU_total)-combined$AU_total,sum(combined$AU_up),
                                lower.tail = FALSE, log.p = FALSE)

combined[which(combined$AN_UP_pvalue < 0.05),'sig_AN'] <- 'Significant'
combined[which(combined$AU_UP_pvalue < 0.05),'sig_AU'] <- 'Significant'


dat[which(dat$value <= 1 ),'sig'] <- 'Non-significant'
dat[which(dat$value > 1),'sig'] <- 'Significant'
dat1 <- subset(dat,Timepoint=="AN50")
dat2 <- subset(dat,Timepoint=="AU50")
dat3 <- subset(dat,SampleType=="UP")


ggplot(dat, aes(x=value, y=Annotation,colour=Timepoint,shape=sig)) +
  scale_colour_manual(name="",values = c("royalblue1","#F24F50","#9B4F0F","#6FB98F"))+
  scale_shape_manual(name="",values=c(1,1,19,19))+facet_wrap(~SampleType,scales="free_x",ncol=2)+
  scale_shape_manual(values  = c(21,21), limits = c('Non-significant', 'Significant')) +
  geom_point(size = 4,stroke = 1)+theme_bw()+
  geom_vline(aes(xintercept=1),colour="Black",size=1,linetype="dashed")+xlab("Enrichment ratio")+
  theme(axis.text.x=element_text(size=11,color="black",face="bold",angle=90), 
        axis.text.y=element_text(size=11,color="black",face="bold"), 
        axis.title=element_text(size=11,face="bold"),text=element_text(size=10))+
  xlim(0, 4.6)  +
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) +
  theme(strip.text = element_text(size=11,face="bold"))



#6
data<-read.table("PRO.txt",header=TRUE,sep="\t")
combined =data##用于计算显著性

AN_UP <- data[,c(1,2,4)]
ratio <- data.frame(cbind(AN_UP,AN_UP$AN_up/AN_UP$AN_total))
dat <- ratio[,c(1,4)]
colnames(dat) <-c("Annotation","AN_UP_ratio")
merAN<- dat


AU_UP <- data[,c(1,3,7)]
ratio2 <- data.frame(cbind(AU_UP,AU_UP$AU_up/AU_UP$AU_total))
dat2 <- ratio2[,c(1,4)]
colnames(dat2) <- c("Annotation","AU_UP_ratio")
merAU <- dat2


rm <- merge(merAN,merAU,by="Annotation")
rm.m <- melt(rm)
rm.m$variable <- factor(rm.m$variable, levels=c("AN_UP_ratio","AU_UP_ratio"))
dat <- rm.m[order(-rm.m[,3]),]
dat$Annotation <- factor(dat$Annotation, levels=c(
  "TSP",
  "MMDE",
  "ACD",
  "FSD",
  "EM (A)",
  "GBM (S)",
  "TSL",
  "AAB (Y)",
  "TPM",
  "RR (S)",
  "Nucleotide M (Y)",
  "NM (A)",
  "OSM (S)",
  "BSM (S)",
  "NSB (Y)",
  "CM (Y)",
  "CFPP (A)",
  "AAM (Y)",
  "Carbon M (Y)",
  "CB",
  "CVM",
  "OAAM (Y)",
  "QS",
  "LM (Y)",
  "FAM (Y)",
  "FAB (Y)",
  "MT (A)",
  "SM (A)",
  "BF (S)",
  "ST",
  "EB (S)",
  "LB (S)",
  "CGD",
  "BC (A)",
  "CC (A)" 
))

samples_new <- sapply(as.character(dat$variable), function(x) strsplit(x, "[_]"))

Timepoint <- {}
for(i in 1:70){
  T <- samples_new[[i]][1]
  Timepoint <- append(Timepoint,T)
}

SampleType <- {}
for(i in 1:70){
  T <- samples_new[[i]][2]
  SampleType <- append(SampleType,T)
}

dat$SampleType = SampleType
dat$Timepoint = Timepoint

combined <- read.table("combined.PRO.txt",header=TRUE,sep="\t")#用于计算显著性

combined$AN_UP_pvalue <- phyper(combined$AN_up-1,combined$AN_total,sum(combined$AN_total)-combined$AN_total,sum(combined$AN_up),
                                lower.tail = FALSE, log.p = FALSE)


combined$AU_UP_pvalue <- phyper(combined$AU_up-1,combined$AU_total,sum(combined$AU_total)-combined$AU_total,sum(combined$AU_up),
                                lower.tail = FALSE, log.p = FALSE)

combined[which(combined$AN_UP_pvalue < 0.05),'sig_AN'] <- 'Significant'
combined[which(combined$AU_UP_pvalue < 0.05),'sig_AU'] <- 'Significant'




dat[which(dat$value <= 1 ),'sig'] <- 'Non-significant'
dat[which(dat$value > 1),'sig'] <- 'Significant'
dat1 <- subset(dat,Timepoint=="AN50")
dat2 <- subset(dat,Timepoint=="AU50")
dat3 <- subset(dat,SampleType=="UP")



ggplot(dat, aes(x=value, y=Annotation,colour=Timepoint,shape=sig)) +
  scale_colour_manual(name="",values = c("royalblue1","#F24F50","#9B4F0F","#6FB98F"))+
  scale_shape_manual(name="",values=c(1,1,19,19))+facet_wrap(~SampleType,scales="free_x",ncol=2)+
  scale_shape_manual(values  = c(21,21), limits = c('Non-significant', 'Significant')) +
  geom_point(size = 4,stroke = 1)+theme_bw()+
  geom_vline(aes(xintercept=1),colour="Black",size=1,linetype="dashed")+xlab("Enrichment ratio")+
  theme(axis.text.x=element_text(size=11,color="black",face="bold",angle=90), 
        axis.text.y=element_text(size=11,color="black",face="bold"), 
        axis.title=element_text(size=11,face="bold"),text=element_text(size=10))+
  xlim(0, 4.6)  +
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) +
  theme(strip.text = element_text(size=11,face="bold"))

######计算组间显著性#######################
#1
combined <- read.table("combined.ACI2.txt",header=TRUE,sep="\t")#用于计算显著性

combined$AN_UP_pvalue <- phyper(combined$AN_up-1,combined$ANAU_up_total,sum(combined$ANAU_up_total)-combined$ANAU_up_total,sum(combined$AN_up),
                                lower.tail = FALSE, log.p = FALSE)

combined$AU_UP_pvalue <- phyper(combined$AU_up-1,combined$ANAU_up_total,sum(combined$ANAU_up_total)-combined$ANAU_up_total,sum(combined$AU_up),
                                lower.tail = FALSE, log.p = FALSE)

combined[which(combined$AN_UP_pvalue < 0.05),'sig_AN'] <- 'Significant'
combined[which(combined$AU_UP_pvalue < 0.05),'sig_AU'] <- 'Significant'



#2
combined <- read.table("combined.ACTINO2.txt",header=TRUE,sep="\t")#用于计算显著性

combined$AN_UP_pvalue <- phyper(combined$AN_up-1,combined$ANAU_up_total,sum(combined$ANAU_up_total)-combined$ANAU_up_total,sum(combined$AN_up),
                                lower.tail = FALSE, log.p = FALSE)

combined$AU_UP_pvalue <- phyper(combined$AU_up-1,combined$ANAU_up_total,sum(combined$ANAU_up_total)-combined$ANAU_up_total,sum(combined$AU_up),
                                lower.tail = FALSE, log.p = FALSE)

combined[which(combined$AN_UP_pvalue < 0.05),'sig_AN'] <- 'Significant'
combined[which(combined$AU_UP_pvalue < 0.05),'sig_AU'] <- 'Significant'



#3
combined <- read.table("combined.BAC2.txt",header=TRUE,sep="\t")#用于计算显著性

combined$AN_UP_pvalue <- phyper(combined$AN_up-1,combined$ANAU_up_total,sum(combined$ANAU_up_total)-combined$ANAU_up_total,sum(combined$AN_up),
                                lower.tail = FALSE, log.p = FALSE)

combined$AU_UP_pvalue <- phyper(combined$AU_up-1,combined$ANAU_up_total,sum(combined$ANAU_up_total)-combined$ANAU_up_total,sum(combined$AU_up),
                                lower.tail = FALSE, log.p = FALSE)

combined[which(combined$AN_UP_pvalue < 0.05),'sig_AN'] <- 'Significant'
combined[which(combined$AU_UP_pvalue < 0.05),'sig_AU'] <- 'Significant'



#4
combined <- read.table("combined.GEMM2.txt",header=TRUE,sep="\t")#用于计算显著性

combined$AN_UP_pvalue <- phyper(combined$AN_up-1,combined$ANAU_up_total,sum(combined$ANAU_up_total)-combined$ANAU_up_total,sum(combined$AN_up),
                                lower.tail = FALSE, log.p = FALSE)

combined$AU_UP_pvalue <- phyper(combined$AU_up-1,combined$ANAU_up_total,sum(combined$ANAU_up_total)-combined$ANAU_up_total,sum(combined$AU_up),
                                lower.tail = FALSE, log.p = FALSE)

combined[which(combined$AN_UP_pvalue < 0.05),'sig_AN'] <- 'Significant'
combined[which(combined$AU_UP_pvalue < 0.05),'sig_AU'] <- 'Significant'


#5
combined <- read.table("combined.PLAN2.txt",header=TRUE,sep="\t")#用于计算显著性

combined$AN_UP_pvalue <- phyper(combined$AN_up-1,combined$ANAU_up_total,sum(combined$ANAU_up_total)-combined$ANAU_up_total,sum(combined$AN_up),
                                lower.tail = FALSE, log.p = FALSE)

combined$AU_UP_pvalue <- phyper(combined$AU_up-1,combined$ANAU_up_total,sum(combined$ANAU_up_total)-combined$ANAU_up_total,sum(combined$AU_up),
                                lower.tail = FALSE, log.p = FALSE)

combined[which(combined$AN_UP_pvalue < 0.05),'sig_AN'] <- 'Significant'
combined[which(combined$AU_UP_pvalue < 0.05),'sig_AU'] <- 'Significant'


#6
combined <- read.table("combined.PRO2.txt",header=TRUE,sep="\t")#用于计算显著性

combined$AN_UP_pvalue <- phyper(combined$AN_up-1,combined$ANAU_up_total,sum(combined$ANAU_up_total)-combined$ANAU_up_total,sum(combined$AN_up),
                                lower.tail = FALSE, log.p = FALSE)

combined$AU_UP_pvalue <- phyper(combined$AU_up-1,combined$ANAU_up_total,sum(combined$ANAU_up_total)-combined$ANAU_up_total,sum(combined$AU_up),
                                lower.tail = FALSE, log.p = FALSE)

combined[which(combined$AN_UP_pvalue < 0.05),'sig_AN'] <- 'Significant'
combined[which(combined$AU_UP_pvalue < 0.05),'sig_AU'] <- 'Significant'



#Fig.3g-k
data<-read.table("enrichment.all.txt",header=TRUE,sep="\t")
combined =data##用于计算显著性


##Armatimonadota
Armatimonadota_UP <- data[,c(1,3,2)]
ratio <- data.frame(cbind(Armatimonadota_UP,Armatimonadota_UP$Armatimonadota_up/Armatimonadota_UP$total))
dat <- ratio[,c(1,4)]
colnames(dat) <-c("Annotation","Armatimonadota_UP_ratio")

merArmatimonadota<- dat



##Chloroflexota
Chloroflexota_UP <- data[,c(1,5,2)]
ratio1 <- data.frame(cbind(Chloroflexota_UP,Chloroflexota_UP$Chloroflexota_up/Chloroflexota_UP$total))
dat1 <- ratio1[,c(1,4)]
colnames(dat1) <-c("Annotation","Chloroflexota_UP_ratio")

merChloroflexota<- dat1


##Eremiobacterota
Eremiobacterota_UP <- data[,c(1,7,2)]
ratio2 <- data.frame(cbind(Eremiobacterota_UP,Eremiobacterota_UP$Eremiobacterota_up/Eremiobacterota_UP$total))
dat2 <- ratio2[,c(1,4)]
colnames(dat2) <-c("Annotation","Eremiobacterota_UP_ratio")

merEremiobacterota<- dat2

##Patescibacteria
Patescibacteria_UP <- data[,c(1,11,2)]
ratio4 <- data.frame(cbind(Patescibacteria_UP,Patescibacteria_UP$Patescibacteria_up/Patescibacteria_UP$total))
dat4 <- ratio4[,c(1,4)]
colnames(dat4) <-c("Annotation","Patescibacteria_UP_ratio")

merPatescibacteria<- dat4

##Verrucomicrobiota
Verrucomicrobiota_UP <- data[,c(1,15,2)]
ratio6 <- data.frame(cbind(Verrucomicrobiota_UP,Verrucomicrobiota_UP$Verrucomicrobiota_up/Verrucomicrobiota_UP$total))
dat6 <- ratio6[,c(1,4)]
colnames(dat6) <-c("Annotation","Verrucomicrobiota_UP_ratio")

merVerrucomicrobiota<- dat6




#strat:
#Armatimonadota
rm <- merArmatimonadota
rm.m <- melt(rm)
rm.m$variable <- factor(rm.m$variable, levels=c("Armatimonadota_UP_ratio"))
dat <- rm.m[order(-rm.m[,3]),]
#dat <- arrange(dat, desc(value))
dat$Annotation <- factor(dat$Annotation, levels=c(
  "TSP",
  "EB (S)",
  "ACD",
  "TSL",
  "FAM (Y)",
  "OAAM (Y)",
  "CFPP (A)",
  "CGD",
  "QS",
  "RR (S)",
  "AAB (Y)",
  "AAM (Y)",
  "MMDE",
  "Carbon M (Y)",
  "FAB (Y)",
  "BSM (S)",
  "FSD",
  "TPM",
  "ST",
  "EM (A)",
  "MT (A)",
  "LM (Y)",
  "CVM",
  "Nucleotide M (Y)",
  "NM (A)",
  "CB",
  "CM (Y)",
  "NSB (Y)",
  "SM (A)",
  "OSM (S)",
  "BC (A)",
  "BF (S)",
  "GBM (S)",
  "LB (S)",
  "CC (A)"
))


samples_new <- sapply(as.character(dat$variable), function(x) strsplit(x, "[_]"))

Timepoint <- {}
for(i in 1:35){
  T <- samples_new[[i]][1]
  Timepoint <- append(Timepoint,T)
}

SampleType <- {}
for(i in 1:35){
  T <- samples_new[[i]][2]
  SampleType <- append(SampleType,T)
}

dat$SampleType = SampleType
dat$Timepoint = Timepoint

combined<-read.table("ARM.combined.txt",header=TRUE,sep="\t")

combined$AN_UP_pvalue <- phyper(combined$Armatimonadota_up-1,combined$total,sum(combined$total)-combined$total,sum(combined$Armatimonadota_up),
                                lower.tail = FALSE, log.p = FALSE)
combined[which(combined$AN_UP_pvalue < 0.05),'sig'] <- 'Significant'

dat[which(dat$value <= 1 ),'sig'] <- 'Non-significant'
dat[which(dat$value > 1),'sig'] <- 'Significant'

ggplot(dat, aes(x=value, y=Annotation,colour=Timepoint,shape=sig)) +
  scale_colour_manual(name="",values = c("#F24F50","#F24F50","royalblue1"))+
  scale_shape_manual(name="",values=c(1,1,19,19))+facet_wrap(~SampleType,scales="free_x",ncol=2)+
  scale_shape_manual(values  = c(21,21), limits = c('Non-significant', 'Significant')) +
  geom_point(size = 4,stroke = 1)+theme_bw()+
  geom_vline(aes(xintercept=1),colour="Black",size=1,linetype="dashed")+xlab("Enrichment ratio")+
  theme(axis.text.x=element_text(size=11,color="black",face="bold",angle=90), 
        axis.text.y=element_text(size=11,color="black",face="bold"), 
        axis.title=element_text(size=11,face="bold"),text=element_text(size=10))+
  xlim(0, 4.6)  +
  theme(legend.text = element_text(colour="black", size = 11, face = "bold"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) +
  theme(strip.text = element_text(size=11,face="bold"))



##Chloroflexota
rm <- merChloroflexota
rm.m <- melt(rm)
rm.m$variable <- factor(rm.m$variable, levels=c("Chloroflexota_UP_ratio"))
dat1 <- rm.m[order(-rm.m[,3]),]
#dat <- arrange(dat, desc(value))
dat1$Annotation <- factor(dat1$Annotation, levels=c(
  "EB (S)",
  "CGD",
  "BF (S)",
  "CC (A)",
  "TSL",
  "NM (A)",
  "OAAM (Y)",
  "ST",
  "FSD",
  "RR (S)",
  "LM (Y)",
  "TPM",
  "FAM (Y)",
  "CB",
  "NSB (Y)",
  "ACD",
  "LB (S)",
  "BSM (S)",
  "BC (A)",
  "SM (A)",
  "CVM",
  "Nucleotide M (Y)",
  "CM (Y)",
  "MMDE",
  "TSP",
  "GBM (S)",
  "EM (A)",
  "OSM (S)",
  "Carbon M (Y)",
  "AAB (Y)",
  "QS",
  "AAM (Y)",
  "CFPP (A)",
  "MT (A)",
  "FAB (Y)"
))


samples_new <- sapply(as.character(dat1$variable), function(x) strsplit(x, "[_]"))

Timepoint <- {}
for(i in 1:35){
  T <- samples_new[[i]][1]
  Timepoint <- append(Timepoint,T)
}

SampleType <- {}
for(i in 1:35){
  T <- samples_new[[i]][2]
  SampleType <- append(SampleType,T)
}

dat1$SampleType = SampleType
dat1$Timepoint = Timepoint

combined<-read.table("CHLO.combined.txt",header=TRUE,sep="\t")
combined$AN_UP_pvalue <- phyper(combined$Chloroflexota_up-1,combined$total,sum(combined$total)-combined$total,sum(combined$Chloroflexota_up),
                                lower.tail = FALSE, log.p = FALSE)

combined[which(combined$AN_UP_pvalue < 0.05),'sig'] <- 'Significant'

dat1[which(dat1$value <= 1 ),'sig'] <- 'Non-significant'
dat1[which(dat1$value > 1),'sig'] <- 'Significant'

ggplot(dat1, aes(x=value, y=Annotation,colour=Timepoint,shape=sig)) +
  scale_colour_manual(name="",values = c("royalblue1","#F24F50","#9B4F0F","#6FB98F"))+
  scale_shape_manual(name="",values=c(1,1,19,19))+facet_wrap(~SampleType,scales="free_x",ncol=2)+
  scale_shape_manual(values  = c(21,21), limits = c('Non-significant', 'Significant')) +
  geom_point(size = 4,stroke = 1)+theme_bw()+
  geom_vline(aes(xintercept=1),colour="Black",size=1,linetype="dashed")+xlab("Enrichment ratio")+
  theme(axis.text.x=element_text(size=11,color="black",face="bold",angle=90), 
        axis.text.y=element_text(size=11,color="black",face="bold"), 
        axis.title=element_text(size=11,face="bold"),text=element_text(size=10))+
  xlim(0, 4.6)  +
  theme(legend.text = element_text(colour="black", size = 11, face = "bold"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) +
  theme(strip.text = element_text(size=11,face="bold"))




#Eremiobacterota
rm <- merEremiobacterota
rm.m <- melt(rm)
rm.m$variable <- factor(rm.m$variable, levels=c("Eremiobacterota_UP_ratio"))
dat2 <- rm.m[order(-rm.m[,3]),]
#dat <- arrange(dat, desc(value))
dat2$Annotation <- factor(dat2$Annotation, levels=c(
  "NM (A)",
  "FSD",
  "RR (S)",
  "CGD",
  "TSP",
  "MT (A)",
  "CB",
  "EB (S)",
  "QS",
  "CM (Y)",
  "EM (A)",
  "BF (S)",
  "CFPP (A)",
  "LM (Y)",
  "Carbon M (Y)",
  "AAB (Y)",
  "BSM (S)",
  "ST",
  "CVM",
  "Nucleotide M (Y)",
  "MMDE",
  "TSL",
  "AAM (Y)",
  "GBM (S)",
  "FAM (Y)",
  "FAB (Y)",
  "NSB (Y)",
  "SM (A)",
  "TPM",
  "OAAM (Y)",
  "OSM (S)",
  "LB (S)",
  "ACD",
  "BC (A)",
  "CC (A)"
))


samples_new <- sapply(as.character(dat2$variable), function(x) strsplit(x, "[_]"))

Timepoint <- {}
for(i in 1:35){
  T <- samples_new[[i]][1]
  Timepoint <- append(Timepoint,T)
}

SampleType <- {}
for(i in 1:35){
  T <- samples_new[[i]][2]
  SampleType <- append(SampleType,T)
}

dat2$SampleType = SampleType
dat2$Timepoint = Timepoint


combined<-read.table("ERM.combined.txt",header=TRUE,sep="\t")
combined$AN_UP_pvalue <- phyper(combined$Eremiobacterota_up-1,combined$total,sum(combined$total)-combined$total,sum(combined$Eremiobacterota_up),
                                lower.tail = FALSE, log.p = FALSE)


combined[which(combined$AN_UP_pvalue < 0.05),'sig'] <- 'Significant'

dat2[which(dat2$value <= 1 ),'sig'] <- 'Non-significant'
dat2[which(dat2$value > 1),'sig'] <- 'Significant'

ggplot(dat2, aes(x=value, y=Annotation,colour=Timepoint,shape=sig)) +
  scale_colour_manual(name="",values = c("royalblue1","#F24F50","#9B4F0F","#6FB98F"))+
  scale_shape_manual(name="",values=c(1,1,19,19))+facet_wrap(~SampleType,scales="free_x",ncol=2)+
  scale_shape_manual(values  = c(21,21), limits = c('Non-significant', 'Significant')) +
  geom_point(size = 4,stroke = 1)+theme_bw()+
  geom_vline(aes(xintercept=1),colour="Black",size=1,linetype="dashed")+xlab("Enrichment ratio")+
  theme(axis.text.x=element_text(size=11,color="black",face="bold",angle=90), 
        axis.text.y=element_text(size=11,color="black",face="bold"), 
        axis.title=element_text(size=11,face="bold"),text=element_text(size=10))+
  xlim(0, 4.6)  +
  theme(legend.text = element_text(colour="black", size = 11, face = "bold"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) +
  theme(strip.text = element_text(size=11,face="bold"))




#Patescibacteria
rm <- merPatescibacteria
rm.m <- melt(rm)
rm.m$variable <- factor(rm.m$variable, levels=c("Patescibacteria_UP_ratio"))
dat4 <- rm.m[order(-rm.m[,3]),]
#dat <- arrange(dat, desc(value))
dat4$Annotation <- factor(dat4$Annotation, levels=c(
  "CC (A)",
  "BC (A)",
  "OSM (S)",
  "CB",
  "AAB (Y)",
  "CGD",
  "CFPP (A)",
  "FAM (Y)",
  "TSL",
  "ACD",
  "TSP",
  "Carbon M (Y)",
  "CVM",
  "TPM",
  "QS",
  "OAAM (Y)",
  "MMDE",
  "AAM (Y)",
  "NSB (Y)",
  "ST",
  "BSM (S)",
  "EM (A)",
  "CM (Y)",
  "MT (A)",
  "LM (Y)",
  "FSD",
  "LB (S)",
  "NM (A)",
  "BF (S)",
  "SM (A)",
  "GBM (S)",
  "Nucleotide M (Y)",
  "RR (S)",
  "FAB (Y)",
  "EB (S)"
))


samples_new <- sapply(as.character(dat4$variable), function(x) strsplit(x, "[_]"))

Timepoint <- {}
for(i in 1:35){
  T <- samples_new[[i]][1]
  Timepoint <- append(Timepoint,T)
}

SampleType <- {}
for(i in 1:35){
  T <- samples_new[[i]][2]
  SampleType <- append(SampleType,T)
}

dat4$SampleType = SampleType
dat4$Timepoint = Timepoint


combined<-read.table("PASE.combined.txt",header=TRUE,sep="\t")
combined$AN_UP_pvalue <- phyper(combined$Patescibacteria_up-1,combined$total,sum(combined$total)-combined$total,sum(combined$Patescibacteria_up),
                                lower.tail = FALSE, log.p = FALSE)



combined[which(combined$AN_UP_pvalue < 0.05),'sig'] <- 'Significant'

dat4[which(dat4$value <= 1 ),'sig'] <- 'Non-significant'
dat4[which(dat4$value > 1),'sig'] <- 'Significant'

ggplot(dat4, aes(x=value, y=Annotation,colour=Timepoint,shape=sig)) +
  scale_colour_manual(name="",values = c("#F24F50","#9B4F0F","#6FB98F"))+
  scale_shape_manual(name="",values=c(1,1,19,19))+facet_wrap(~SampleType,scales="free_x",ncol=2)+
  scale_shape_manual(values  = c(21,21), limits = c('Non-significant', 'Significant')) +
  geom_point(size = 4,stroke = 1)+theme_bw()+
  geom_vline(aes(xintercept=1),colour="Black",size=1,linetype="dashed")+xlab("Enrichment ratio")+
  theme(axis.text.x=element_text(size=11,color="black",face="bold",angle=90), 
        axis.text.y=element_text(size=11,color="black",face="bold"), 
        axis.title=element_text(size=11,face="bold"),text=element_text(size=10))+
  xlim(0, 4.6)  +
  theme(legend.text = element_text(colour="black", size = 11, face = "bold"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) +
  theme(strip.text = element_text(size=11,face="bold"))



#Verrucomicrobiota
rm <- merVerrucomicrobiota
rm.m <- melt(rm)
rm.m$variable <- factor(rm.m$variable, levels=c("Verrucomicrobiota_UP_ratio"))
dat6 <- rm.m[order(-rm.m[,3]),]
#dat <- arrange(dat, desc(value))
dat6$Annotation <- factor(dat6$Annotation, levels=c(
  "BC (A)",
  "CFPP (A)",
  "CC (A)",
  "ACD",
  "FAM (Y)",
  "NM (A)",
  "QS",
  "Carbon M (Y)",
  "MMDE",
  "MT (A)",
  "CGD",
  "ST",
  "RR (S)",
  "LM (Y)",
  "CM (Y)",
  "FSD",
  "AAB (Y)",
  "EM (A)",
  "OSM (S)",
  "OAAM (Y)",
  "BSM (S)",
  "BF (S)",
  "AAM (Y)",
  "CVM",
  "FAB (Y)",
  "Nucleotide M (Y)",
  "EB (S)",
  "CB",
  "TPM",
  "SM (A)",
  "GBM (S)",
  "NSB (Y)",
  "TSL",
  "TSP",
  "LB (S)"
  
))


samples_new <- sapply(as.character(dat6$variable), function(x) strsplit(x, "[_]"))

Timepoint <- {}
for(i in 1:35){
  T <- samples_new[[i]][1]
  Timepoint <- append(Timepoint,T)
}

SampleType <- {}
for(i in 1:35){
  T <- samples_new[[i]][2]
  SampleType <- append(SampleType,T)
}

dat6$SampleType = SampleType
dat6$Timepoint = Timepoint


combined<-read.table("VERRU.combined.txt",header=TRUE,sep="\t")
combined$AN_UP_pvalue <- phyper(combined$Verrucomicrobiota_up-1,combined$total,sum(combined$total)-combined$total,sum(combined$Verrucomicrobiota_up),
                                lower.tail = FALSE, log.p = FALSE)



combined[which(combined$AN_UP_pvalue < 0.05),'sig'] <- 'Significant'

dat6[which(dat6$value <= 1 ),'sig'] <- 'Non-significant'
dat6[which(dat6$value > 1),'sig'] <- 'Significant'

ggplot(dat6, aes(x=value, y=Annotation,colour=Timepoint,shape=sig)) +
  scale_colour_manual(name="",values = c("royalblue1","#F24F50","#9B4F0F","#6FB98F"))+
  scale_shape_manual(name="",values=c(1,1,19,19))+facet_wrap(~SampleType,scales="free_x",ncol=2)+
  scale_shape_manual(values  = c(21,21), limits = c('Non-significant', 'Significant')) +
  geom_point(size = 4,stroke = 1)+theme_bw()+
  geom_vline(aes(xintercept=1),colour="Black",size=1,linetype="dashed")+xlab("Enrichment ratio")+
  theme(axis.text.x=element_text(size=11,color="black",face="bold",angle=90), 
        axis.text.y=element_text(size=11,color="black",face="bold"), 
        axis.title=element_text(size=11,face="bold"),text=element_text(size=10))+
  xlim(0, 4.6)  +
  theme(legend.text = element_text(colour="black", size = 11, face = "bold"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) +
  theme(strip.text = element_text(size=11,face="bold"))





