#Fig.S15-S21

#R3.6.2
library(vegan)
library(ggplot2)

library(SoDA)
library(picante)
library(vegan)
library(ecodist)
library(MASS)
library(hier.part)
library(car)
library(MuMIn)


library(vegan)
library(PMCMRplus)
library(MASS)
library(permute)
library(lattice)
library(car)
library(nlme)
library(mgcv)
library(labdsv)
library(sciplot)
library(spdep)
library(agricolae)

#####Fig.S15
#Fig.S15a
rm(list = ls())
otu.16S <- read.delim('data/amplicon.OTU_NMDS.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
ko.16S <- read.delim('data/shotgun.KO_NMDS.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)

# Procrustes analysis，details ?procrustes
proc <- procrustes(Y = otu.16S, X = ko.16S, symmetric = TRUE)
summary(proc)

#rotation plot
plot(proc, kind = 1, type = 'text')
#main result
names(proc)
head(proc$Yrot)  #Procrustes  Y 
head(proc$X)  #Procrustes  X 
proc$ss  # M2 
proc$rotation  #rotation coordinate axis
#residual error plot
plot(proc, kind = 2)
residuals(proc)  

#PROTEST test，datails ?protest
set.seed(123)
prot <- protest(Y = otu.16S, X = ko.16S, permutations = how(nperm = 999))
prot

#main result
names(prot)
prot$signif  #p value
prot$ss  #M2 

library(ggplot2)
#Procrustes coordinate axis
Y <- cbind(data.frame(proc$Yrot), data.frame(proc$X))
X <- data.frame(proc$rotation)

#group
group <- read.delim('data/group_procrustes.txt', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
Y$samples <- rownames(Y)
Y <- merge(Y, group, by = 'samples')


Y2 <- read.csv("data/Y2.amplicon.shotgun.csv",head=T,row.names=1)
##border
df_points <- read.csv("data/Y2.amplicon.shotgun.csv",head=T,row.names=1)
library(plyr)
df <-df_points[,c("X", "Y", "Type")]
find_hull <- function(df_points) df_points[chull(df_points$X, df_points$Y), ]
hulls <- ddply(df, "Type", find_hull )###border

#plot
p1<-ggplot(data=Y2)+ 
    geom_point(aes(x=X, y=Y, color = treat,shape=Category), size = 3.5,alpha=0.7) +
    scale_color_manual(values = c("#C7E2FF","lightblue1","#95E9FF","#5EB5FF","#339CFF","blue",
                                  "#FFE1FF","#FFC8FF","#FFAFFF","#FF63FF","#FF6371","red"), 
                      limits = c('AN00','AN02','AN05','AN10','AN20','AN50',
                                 'AU00','AU02','AU05','AU10','AU20','AU50')) +
  scale_shape_manual(values=c(19, 1))+
  
  geom_polygon(data=hulls,aes(x=X,y=Y,fill = Type), alpha=.1)+
  scale_fill_manual(values=c("blue","red"))+
  geom_segment(data=Y, aes(x = X1, y = X2, xend = x, yend = y,color = group), 
             arrow = arrow(length = unit(0.3, 'cm')),size = 0.6) +
  theme_bw()+
  guides(color=guide_legend(title= "Treatment"))+
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        legend.key = element_rect(fill = 'transparent')) +
  theme(legend.position = 'left',legend.direction = 'vertical',
        legend.title = element_text(colour="black", size=13, face="bold"),
        legend.text = element_text(colour="black", size=12),
        axis.text=element_text(size=15,face="bold",colour="black"),
        axis.title=element_text(size=17,face="bold"),
        plot.title=element_text(face = "bold", colour = "black", size = 15))+
  labs(x = 'NMDS1', y = 'NMDS2', color = '') +
  geom_vline(xintercept = 0, color = 'gray34', linetype = 2, linewidth = 1) +
  geom_hline(yintercept = 0, color = 'gray34', linetype = 2, linewidth = 1) +
  geom_abline(intercept = 0, slope = X[1,2]/X[1,1], linewidth = 1.3) +
  geom_abline(intercept = 0, slope = X[2,2]/X[2,1], linewidth = 1.3) +
  ggtitle(expression(M ^2~'= 0.113***'))

p1
ggsave("Fig.S15a.pdf", p1,height=5.6,width =6.5,limitsize = FALSE )


#Fig.15b
#R3.6.2
library(EnhancedVolcano)
library(patchwork)
library(edgeR)
library(statmod)

rm(list = ls())
gene5<-read.csv("data/2DP_AN50.AU50.csv",row.names=1,head=T)

gene5[which(gene5$FDR < 0.05 & gene5$AU <= -1),'sig'] <- 'Down'
gene5[which(gene5$FDR < 0.05 & gene5$AU >=1),'sig'] <- 'Up'
gene5[which(gene5$FDR >= 0.05 | abs(gene5$AU) < 1),'sig'] <- 'None'

sum(gene5$sig=='Down')
sum(gene5$sig=="Up")
sum(gene5$sig=="None")

p2<-ggplot(data=gene5, aes(x= AU , y= FDR,group=sig))+ 
  geom_point(aes(colour=sig),size=4,alpha=0.5)+ 
  labs(x = "2DP", y = "P value (FDR)")+
  scale_colour_manual(values=c("red","blue","gray70","#5EB5FF","#339CFF","blue",
                               "#FFE1FF","#FFC8FF","#FFAFFF","#FF63FF","#FF6371"),
                      breaks=c("Up","Down", "None"),
                      labels=c("AU preference", "AN preference","None"))+
  
  theme_bw()+
  ggtitle(expression('AU versus AN'),subtitle = "KOs in MAGs: Total = 6715 variables")+
  theme(legend.position = 'left',legend.direction = 'vertical',
        legend.title = element_text(colour="black", size=13, face="bold"),
        legend.text = element_text(colour="black", size=12),
        axis.text=element_text(size=15,face="bold",colour="black"),
        axis.title=element_text(size=17,face="bold"))+
  geom_vline(xintercept = 1, color = 'gray34', linetype = 2, linewidth = 1.3) +
  geom_vline(xintercept = -1, color = 'gray34', linetype = 2, linewidth = 1.3) +
  geom_hline(yintercept = 0.05, color = 'gray34', linetype = 2, linewidth = 1.3) 

p2
ggsave("Fig.S15b", p2,height=4.5,width =6,limitsize = FALSE )


#Fig.S15c
#R4.2.3
library(ggplot2)
library("phyloseq")
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

rm(list = ls())
data<-read.table("data/2DP_ANAU50.enrichment.YAS.function.txt",header=TRUE,sep="\t")

#AN50
AN_UP <- data[,c(1,2,4)]
ratio <- data.frame(cbind(AN_UP,AN_UP$AN_up/AN_UP$total))
dat <- ratio[,c(1,4)]
colnames(dat) <-c("Annotation","AN_UP_ratio")

merAN<- dat

#AU50
AU_UP <- data[,c(1,3,4)]
ratio2 <- data.frame(cbind(AU_UP,AU_UP$AU_up/AU_UP$total))
dat2 <- ratio2[,c(1,4)]
colnames(dat2) <- c("Annotation","AU_UP_ratio")

merAU <- dat2


rm <- merge(merAN,merAU,by="Annotation")
rm.m <- melt(rm)
rm.m$variable <- factor(rm.m$variable, levels=c("AN_UP_ratio","AU_UP_ratio"))
dat <- rm.m[order(-rm.m[,3]),]
dat$Annotation <- factor(dat$Annotation, levels=c(
  "Fatty acid metabolism(Y)",
  "Amino acid metabolism(Y)",
  "Biosynthesis of cofactors",
  "Nucleotide metabolism(Y)",
  "Degradation of aromatic compounds",
  "Sulfur metabolism(A)",
  "Metabolism of cofactors and vitamins",
  "Bacterial chemotaxis(A)",
  "Lipid metabolism(Y)",
  "Cell motility(A)",
  "Microbial metabolism in diverse environments",
  "Energy metabolism(A)",
  "Metabolism of other amino acids(Y)",
  "Carbon fixation pathways in prokaryotes(A)",
  "Glycan biosynthesis and metabolism(S)",
  "Transcription",
  "Biosynthesis of other secondary metabolites(S)",
  "Quorum sensing",
  "Metabolism of terpenoids and polyketides",
  "Biosynthesis of nucleotide sugars(Y)",
  "Signal transduction",
  "Biofilm formation(S)",
  "Nitrogen metabolism(A)",
  "Lipopolysaccharide biosynthesis(S)",
  "Carbohydrate metabolism(Y)",
  "Membrane transport(A)",
  "Biosynthesis of amino acids(Y)",
  "Biosynthesis of secondary metabolites(S)",
  "Carbon metabolism(Y)",
  "Fatty acid biosynthesis(Y)",
  "Translation",
  "Exopolysaccharide biosynthesis(S)",
  "Replication and repair(S)",
  "Folding, sorting and degradation"
))

samples_new <- sapply(as.character(dat$variable), function(x) strsplit(x, "[_]"))

Timepoint <- {}
for(i in 1:68){
  T <- samples_new[[i]][1]
  Timepoint <- append(Timepoint,T)
}

SampleType <- {}
for(i in 1:68){
  T <- samples_new[[i]][2]
  SampleType <- append(SampleType,T)
}

dat$SampleType = SampleType
dat$Timepoint = Timepoint

combined <- read.table("data/2DP_AnAU50.enrichment.YAS.function.combined.txt",header=TRUE,sep="\t")#用于计算显著性##用于计算显著性
combined$AN_UP_pvalue <- phyper(combined$AN_up-1,combined$total,sum(combined$total)-combined$total,sum(combined$AN_up),
                                lower.tail = FALSE, log.p = FALSE)
combined$AU_UP_pvalue <- phyper(combined$AU_up-1,combined$total,sum(combined$total)-combined$total,sum(combined$AU_up),
                                lower.tail = FALSE, log.p = FALSE)


combined[which(combined$AN_UP_pvalue < 0.05),'sigAN50'] <- 'Significant'
combined[which(combined$AU_UP_pvalue < 0.05),'sigAU50'] <- 'Significant'


dat[which(dat$value <= 1 ),'sig'] <- 'Non-entiched'
dat[which(dat$value > 1),'sig'] <- 'enriched'


ggplot(dat, aes(x=value, y=Annotation,colour=Timepoint,shape=sig)) +
  scale_colour_manual(name="",values = c("royalblue1","#F24F50","#9B4F0F","#6FB98F"))+
  scale_shape_manual(name="",values=c(1,1,19,19))+facet_wrap(~SampleType,scales="free_x",ncol=2)+
  scale_shape_manual(values  = c(21,21), limits = c('Non-entiched', 'enriched')) +
  geom_point(size = 4,stroke = 1)+theme_bw()+
  geom_vline(aes(xintercept=1),colour="Black",size=1,linetype="dashed")+xlab("Enrichment ratio")+
  theme(axis.text.x=element_text(size=11,color="black",face="bold",angle=90), 
        axis.text.y=element_text(size=11,color="black",face="bold"), 
        axis.title=element_text(size=11,face="bold"),text=element_text(size=10))+
  xlim(0, 3)  +
  theme(legend.text = element_text(colour="black", size = 11, face = "bold"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) +
  theme(strip.text = element_text(size=11,face="bold"))


####Fig.S16
#R.3.6.2
library(vegan)
library(ggplot2)

library(SoDA)
library(picante)
library(vegan)
library(ecodist)
library(MASS)
library(hier.part)
library(car)
library(MuMIn)

library(vegan)
library(PMCMRplus)
library(MASS)
library(permute)
library(lattice)
library(car)
library(nlme)
library(mgcv)
library(labdsv)
library(sciplot)
library(spdep)
library(agricolae)

rm(list = ls())
otu.16S <- read.delim('data/shotgun.OTU_NMDS.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
ko.16S <- read.delim('data/shotgun.KO_NMDS.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)

#do Procrustes analysis，details ?procrustes
proc <- procrustes(Y = otu.16S, X = ko.16S, symmetric = TRUE)
summary(proc)

#rotation plot
plot(proc, kind = 1, type = 'text')

names(proc)
head(proc$Yrot)  #Procrustes  Y 
head(proc$X)  #Procrustes  X 
proc$ss  # M2 
proc$rotation  #rotation coordinate axis
#residual error plot
plot(proc, kind = 2)
residuals(proc)  

#PROTEST test，datails ?protest
set.seed(123)
prot <- protest(Y = otu.16S, X = ko.16S, permutations = how(nperm = 999))
prot

#main result
names(prot)
prot$signif  #p value
prot$ss  # M2 

#Procrustes coordinate axis
Y <- cbind(data.frame(proc$Yrot), data.frame(proc$X))
X <- data.frame(proc$rotation)

#group
group <- read.delim('data/group_procrustes.txt', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)

Y$samples <- rownames(Y)
Y <- merge(Y, group, by = 'samples')
Y2 <- read.csv("data/Y2.shotgun.shotgun.csv",head=T,row.names=1)

##border
df_points <- read.csv("data/Y2.shotgun.shotgun.csv",head=T,row.names=1)
library(plyr)
df <-df_points[,c("X", "Y", "Type")]
find_hull <- function(df_points) df_points[chull(df_points$X, df_points$Y), ]
hulls <- ddply(df, "Type", find_hull )###寻找边界

#####plot
p3<-ggplot(data=Y2)+ 
  geom_point(aes(x=X, y=Y, color = treat,shape=Category), size = 3.5,alpha=0.7) +
  scale_color_manual(values = c("#C7E2FF","lightblue1","#95E9FF","#5EB5FF","#339CFF","blue",
                                "#FFE1FF","#FFC8FF","#FFAFFF","#FF63FF","#FF6371","red"), 
                     limits = c('AN00','AN02','AN05','AN10','AN20','AN50',
                                'AU00','AU02','AU05','AU10','AU20','AU50')) +
  scale_shape_manual(values=c(19, 1))+
  
  geom_polygon(data=hulls,aes(x=X,y=Y,fill = Type), alpha=.1)+
  scale_fill_manual(values=c("blue","red","yellow"))+
  geom_segment(data=Y, aes(x = X1, y = X2, xend = x, yend = y,color = group), 
               arrow = arrow(length = unit(0.3, 'cm')),size = 0.6) +
  theme_bw()+
  guides(color=guide_legend(title= "Treatment"))+
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        legend.key = element_rect(fill = 'transparent')) +
  theme(legend.position = 'left',legend.direction = 'vertical',
        legend.title = element_text(colour="black", size=13, face="bold"),
        legend.text = element_text(colour="black", size=12),
        axis.text=element_text(size=15,face="bold",colour="black"),
        axis.title=element_text(size=17,face="bold"),
        plot.title=element_text(face = "bold", colour = "black", size = 15))+
  labs(x = 'NMDS1', y = 'NMDS2', color = '') +
  geom_vline(xintercept = 0, color = 'gray34', linetype = 2, linewidth = 1) +
  geom_hline(yintercept = 0, color = 'gray34', linetype = 2, linewidth = 1) +
  geom_abline(intercept = 0, slope = X[1,2]/X[1,1], linewidth = 1.3) +
  geom_abline(intercept = 0, slope = X[2,2]/X[2,1], linewidth = 1.3) +
  ggtitle(expression(M ^2~'= 0.127***'))

p3
ggsave("Fig.S16.pdf", p3,height=5.6,width =6.5,limitsize = FALSE )



####Fig.S17
library(ggplot2)
mag1 <- read.csv("data/all.N.dreped.MAGs.csv")

test_rect <- data.frame(xmin = 80,
                        xmax = 100,
                        ymin = 0,
                        ymax = 5,
                        grp = "High quality (n = 134)")
test_rect

p1 <- ggplot(mag1) + 
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent'),
        legend.title = element_text(face = "bold", colour = "black", size = 11)) +
  labs(x = 'Completeness %', y = 'Contamination %')+
  
  geom_point(aes(Completeness,  Contamination,colour=Phylum),size = 2.5) +
  geom_rect(data = test_rect,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = grp), alpha = 0.15) +
  guides(fill = guide_legend(" ",order=0),
         colour = guide_legend("Phylum",order=1)) +
  
  ggtitle("Non-redundant MAGs from N treatment; n = 366")+
  geom_vline(xintercept=80,size=1,linetype=2,color="gray30")+
  
  geom_hline(yintercept=5, size=1,linetype=2,color="gray30")+
  
  theme(panel.grid.major = element_line(color = 'gray', size = 0.2),
        axis.title = element_text(face = "bold", size = 12,color = 'black'),
        axis.text.x =element_text(color="black",size = 11),
        axis.text.y =element_text(color="black",size = 11))+
  scale_color_manual(
    
    #Phylum
    breaks = c("Actinobacteriota (47.0%)", "Proteobacteria (23.5%)", "Patescibacteria (7.1%)","Bacteroidota (6.0%)",
               "Gemmatimonadota (4.9%)","Acidobacteriota (2.7%)","Verrucomicrobiota (2.2%)","Chloroflexota (1.6%)","Thermoproteota (1.6%)",
               "Eremiobacterota (1.4%)","Armatimonadota (0.8%)","Planctomycetota (0.8%)","Myxococcota (0.3%)"),
    values = c("lightgoldenrod1","paleturquoise1","slateblue2","violet","seagreen1","lightskyblue","wheat1","mediumpurple1",'gray43',
               "darkseagreen3",'sienna1', "lightpink2",'deepskyblue'))    
p1
ggsave("Fig.S17", p1,height=5,width =6.6,limitsize = FALSE )


####Fig.S18
library(vegan)
library(ggplot2)

library(SoDA)
library(picante)
library(PCNM)
library(vegan)
library(ecodist)
library(MASS)
library(hier.part)
library(car)
library(MuMIn)


library(vegan)
library(PMCMRplus)
library(MASS)
library(permute)
library(lattice)
library(car)
library(nlme)
library(mgcv)
library(labdsv)
library(sciplot)
library(spdep)
library(agricolae)


#Fig.S20a
rm(list = ls())
otu.16S <- read.delim('data/ANAU.16S.OTU.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
ko.16S <- read.delim('data/ANAU.traits.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)

#执行 Procrustes 分析，详情 ?procrustes
#以对称分析为例（symmetric = TRUE）
proc <- procrustes(Y = otu.16S, X = ko.16S, symmetric = TRUE)
summary(proc)
#???procrustes
#旋转图
plot(proc, kind = 1, type = 'text')
#一些重要的结果提取
names(proc)

head(proc$Yrot)  #Procrustes 分析后 Y 的坐标
head(proc$X)  #Procrustes 分析后 X 的坐标
proc$ss  #偏差平方和 M2 统计量
proc$rotation  #通过该值可获得旋转轴的坐标位置
#残差图
plot(proc, kind = 2)
residuals(proc)  #残差值
#PROTEST 检验，详情 ?protest
#以 999 次置换为例
#注：protest() 中执行的是对称 Procrustes 分析，X 和 Y 的分配调换不影响 M2 统计量的计算
set.seed(123)
prot <- protest(Y = otu.16S, X = ko.16S, permutations = how(nperm = 999))
prot

#重要统计量的提取
names(prot)
prot$signif  #p 值
prot$ss  #偏差平方和 M2 统计量



library(ggplot2)

#提取 Procrustes 分析的坐标
Y <- cbind(data.frame(proc$Yrot), data.frame(proc$X))
X <- data.frame(proc$rotation)

#添加分组信息
group <- read.delim('data/ANAU.group.txt', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)

Y$samples <- rownames(Y)
Y <- merge(Y, group, by = 'samples')
#write.csv(Y,"Y.16S.KO.csv")##整理成txt格式,加入type信息
####################################################
Y2 <- read.csv("data/Y2.16S.OTU_traits.csv",head=T,row.names=1)
##边界
df_points <- read.csv("data/Y2.16S.OTU_traits.csv",head=T,row.names=1)
library(plyr)
df <-df_points[,c("X", "Y", "Type")]
find_hull <- function(df_points) df_points[chull(df_points$X, df_points$Y), ]
hulls <- ddply(df, "Type", find_hull )###寻找边界


p1<-ggplot(data=Y2)+ 
  geom_point(aes(x=X, y=Y, color = treat,shape=Category), size = 3.5,alpha=0.7) +
  scale_color_manual(values = c("#C7E2FF","lightblue1","#95E9FF","#5EB5FF","#339CFF","blue",
                                "#FFE1FF","#FFC8FF","#FFAFFF","#FF63FF","#FF6371","red"), 
                     limits = c('AN00','AN02','AN05','AN10','AN20','AN50',
                                'AU00','AU02','AU05','AU10','AU20','AU50')) +
  scale_shape_manual(values=c(19, 1))+
  
  geom_polygon(data=hulls,aes(x=X,y=Y,fill = Type), alpha=.1)+
  #geom_polygon(data=hulls2, alpha=.1,aes(fill = Type))+
  scale_fill_manual(values=c("blue","red"))+
  
  
  geom_segment(data=Y, aes(x = X1, y = X2, xend = x, yend = y,color = group), 
               arrow = arrow(length = unit(0.3, 'cm')),size = 0.6) +
  
  theme_bw()+
  guides(color=guide_legend(title= "Treatment"))+
  
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        legend.key = element_rect(fill = 'transparent')) +
  
  theme(legend.position = 'left',legend.direction = 'vertical',
        legend.title = element_text(colour="black", size=13, face="bold"),
        legend.text = element_text(colour="black", size=12),
        axis.text=element_text(size=15,face="bold",colour="black"),
        axis.title=element_text(size=17,face="bold"),
        plot.title=element_text(face = "bold", colour = "black", size = 15))+
  labs(x = 'NMDS1', y = 'NMDS2', color = '') +
  geom_vline(xintercept = 0, color = 'gray34', linetype = 2, linewidth = 1) +
  geom_hline(yintercept = 0, color = 'gray34', linetype = 2, linewidth = 1) +
  geom_abline(intercept = 0, slope = X[1,2]/X[1,1], linewidth = 1.3) +
  geom_abline(intercept = 0, slope = X[2,2]/X[2,1], linewidth = 1.3) +
  
  ggtitle(expression(M ^2~'= 0.315***'))
p1
ggsave("Fig.S20a.pdf", p1,height=5.6,width =6.5,limitsize = FALSE )





#Fig.S20b
rm(list = ls())
otu.16S <- read.delim('data/ANAU.shotgun.OTU.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
ko.16S <- read.delim('data/ANAU.traits.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)

#执行 Procrustes 分析，详情 ?procrustes
#以对称分析为例（symmetric = TRUE）
proc <- procrustes(Y = otu.16S, X = ko.16S, symmetric = TRUE)
summary(proc)

#旋转图
plot(proc, kind = 1, type = 'text')
#一些重要的结果提取
names(proc)

head(proc$Yrot)  #Procrustes 分析后 Y 的坐标
head(proc$X)  #Procrustes 分析后 X 的坐标
proc$ss  #偏差平方和 M2 统计量
proc$rotation  #通过该值可获得旋转轴的坐标位置
#残差图
plot(proc, kind = 2)
residuals(proc)  #残差值
#PROTEST 检验，详情 ?protest
#以 999 次置换为例
#注：protest() 中执行的是对称 Procrustes 分析，X 和 Y 的分配调换不影响 M2 统计量的计算
set.seed(123)
prot <- protest(Y = otu.16S, X = ko.16S, permutations = how(nperm = 999))
prot

#重要统计量的提取
names(prot)
prot$signif  #p 值
prot$ss  #偏差平方和 M2 统计量



library(ggplot2)

#提取 Procrustes 分析的坐标
Y <- cbind(data.frame(proc$Yrot), data.frame(proc$X))
X <- data.frame(proc$rotation)

#添加分组信息
group <- read.delim('data/ANAU.group.txt', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)

Y$samples <- rownames(Y)
Y <- merge(Y, group, by = 'samples')
#write.csv(Y,"Y.shotgun.16S.Traits.csv")##整理成txt格式,加入type信息


Y2 <- read.csv("data/Y2.shotgun.16S.Traits.csv",head=T,row.names=1)
##边界
df_points <- read.csv("data/Y2.shotgun.16S.Traits.csv",head=T,row.names=1)
library(plyr)
df <-df_points[,c("X", "Y", "Type")]
find_hull <- function(df_points) df_points[chull(df_points$X, df_points$Y), ]
hulls <- ddply(df, "Type", find_hull )###寻找边界



p1<-ggplot(data=Y2)+ 
  geom_point(aes(x=X, y=Y, color = treat,shape=Category), size = 3.5,alpha=0.7) +
  scale_color_manual(values = c("#C7E2FF","lightblue1","#95E9FF","#5EB5FF","#339CFF","blue",
                                "#FFE1FF","#FFC8FF","#FFAFFF","#FF63FF","#FF6371","red"), 
                     limits = c('AN00','AN02','AN05','AN10','AN20','AN50',
                                'AU00','AU02','AU05','AU10','AU20','AU50')) +
  scale_shape_manual(values=c(19, 1))+
  
  geom_polygon(data=hulls,aes(x=X,y=Y,fill = Type), alpha=.1)+
  #geom_polygon(data=hulls2, alpha=.1,aes(fill = Type))+
  scale_fill_manual(values=c("blue","red","yellow"))+
  
  
  geom_segment(data=Y, aes(x = X1, y = X2, xend = x, yend = y,color = group), 
               arrow = arrow(length = unit(0.3, 'cm')),size = 0.6) +
  
  theme_bw()+
  guides(color=guide_legend(title= "Treatment"))+
  
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        legend.key = element_rect(fill = 'transparent')) +
  
  theme(legend.position = 'left',legend.direction = 'vertical',
        legend.title = element_text(colour="black", size=13, face="bold"),
        legend.text = element_text(colour="black", size=12),
        axis.text=element_text(size=15,face="bold",colour="black"),
        axis.title=element_text(size=17,face="bold"),
        plot.title=element_text(face = "bold", colour = "black", size = 15))+
  labs(x = 'NMDS1', y = 'NMDS2', color = '') +
  geom_vline(xintercept = 0, color = 'gray34', linetype = 2, linewidth = 1) +
  geom_hline(yintercept = 0, color = 'gray34', linetype = 2, linewidth = 1) +
  geom_abline(intercept = 0, slope = X[1,2]/X[1,1], linewidth = 1.3) +
  geom_abline(intercept = 0, slope = X[2,2]/X[2,1], linewidth = 1.3) +
  
  ggtitle(expression(M ^2~'= 0.332***'))
p1
ggsave("Fig.S20b.pdf", p1,height=5.6,width =6.5,limitsize = FALSE )



####Fig.S21
#R.3.6.2#
library(ggplot2)
library(vegan)
library(colorRamps)
library(psych)
library(tidyr)
library(nlme)
library(MuMIn)
library(ggbreak)

library(tidyverse)
library(ggh4x)
library(ggsci)
library(patchwork)
library(readxl)
library(writexl)
library(agricolae)
library(broom)



rm(list=ls())
env.amp1 <- read.csv("data/weighed_traits.csv",head=T,row.names=1)

env.amp.N <- env.amp1[1:96,]
env.amp.AN <- env.amp1[1:48,]
env.amp.AU <- env.amp1[49:96,]


###S.KO
p1 <-ggplot()+
  geom_jitter(data=env.amp.N, aes(y=S.KO, x= log(level3+1), color=type1),size=10,alpha=0.8)+
  geom_smooth( data=env.amp.N,aes(y=S.KO, x= log(level3+1),color=type1),method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('Ln(N addition rate ('*g~N~m ^-2~yr^-1*'))')) +
  ylab(bquote('S.KO (CWMs)'))+
  scale_color_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+ 

  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R)["AN"] ^2~'= 0.540***'~~~~~italic(R)["AU"] ^2~'= 0.440***'))

p1
summary(lm(env.amp.AN$S.KO~log(env.amp.AN$level3+1)))
summary(lm(env.amp.AU$S.KO~log(env.amp.AU$level3+1)))


#####genome_size
p2 <-ggplot()+
  geom_jitter(data=env.amp.N, aes(y=genome_size, x= log(level3+1), color=type1),size=10,alpha=0.8)+
  geom_smooth(data=subset(env.amp.N,env.amp.N$type1=="AN"),aes(x= log(level3+1), y=genome_size,color=type1), method ="lm", span=0.9,linewidth=2)+
  geom_smooth(data=subset(env.amp.N,env.amp.N$type1=="AU"),aes(x= log(level3+1), y=genome_size,color=type1), method ="lm", span=0.9,linewidth=2)+
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('Ln(N addition rate ('*g~N~m ^-2~yr^-1*'))')) +
  ylab(bquote('Genome Size (CWMs)'))+
  scale_color_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+ 
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R)["AN"] ^2~'= 0.299***'~~~~~italic(R)["AU"] ^2~'= 0.286***'))

p2
summary(lm(env.amp.AN$genome_size~log(env.amp.AN$level3+1)))
summary(lm(env.amp.AU$genome_size~log(env.amp.AU$level3+1)))



####S.ARG
p3 <-ggplot()+
  geom_jitter(data=env.amp.N, aes(y=S.ARG, x= log(level3+1), color=type1),size=10,alpha=0.8)+
  geom_smooth(data=subset(env.amp.N,env.amp.N$type1=="AN"),aes(x= log(level3+1), y=S.ARG,color=type1), method ="lm", span=0.9,linewidth=2)+
  geom_smooth(data=subset(env.amp.N,env.amp.N$type1=="AU"),aes(x= log(level3+1), y=S.ARG,color=type1), method ="lm", span=0.9,linewidth=2)+
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('Ln(N addition rate ('*g~N~m ^-2~yr^-1*'))')) +
  ylab(bquote('S.ARG (CWMs)'))+
  scale_color_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+ 
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R)["AN"] ^2~'= 0.515***'~~~~~italic(R)["AU"] ^2~'= 0.214***'))

p3
summary(lm(env.amp.AN$S.ARG~log(env.amp.AN$level3+1)))
summary(lm(env.amp.AU$S.ARG~log(env.amp.AU$level3+1)))


###d
p4 <-ggplot()+
  geom_jitter(data=env.amp.N, aes(y=d, x= log(level3+1), color=type1),size=10,alpha=0.8)+
  geom_smooth( data=env.amp.N,aes(y=d, x= log(level3+1),color=type1),method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('Ln(N addition rate ('*g~N~m ^-2~yr^-1*'))')) +
  ylab(bquote('Minimal Doubling Time (CWMs)'))+
  scale_color_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+ 
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R)["AN"] ^2~'= 0.673***'~~~~~italic(R)["AU"] ^2~'= 0.585***'))

p4
summary(lm(env.amp.AN$d~log(env.amp.AN$level3+1)))
summary(lm(env.amp.AU$d~log(env.amp.AU$level3+1)))


###GC
p5 <-ggplot()+
  geom_jitter(data=env.amp.N, aes(y=GC, x= log(level3+1), color=type1),size=10,alpha=0.8)+
  geom_smooth( data=env.amp.N,aes(y=GC, x= log(level3+1),color=type1),method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('Ln(N addition rate ('*g~N~m ^-2~yr^-1*'))')) +
  ylab(bquote('GC Content (% CWMs)'))+
  scale_color_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_y_continuous(labels = scales::percent)+
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+ 
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R)["AN"] ^2~'= 0.395***'~~~~~italic(R)["AU"] ^2~'= 0.325***'))
p5
summary(lm(env.amp.AN$GC~log(env.amp.AN$level3+1)))
summary(lm(env.amp.AU$GC~log(env.amp.AU$level3+1)))


p <- 
  p1 + p3 + p4 + p5 + p2 +
  plot_layout(guides = "collect",ncol=5,nrow=1) +
  theme(legend.position='right')

p
ggsave("Fig.S21.pdf", p,height=7,width =32,limitsize = FALSE )
