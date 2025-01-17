#Fig.S9-S11


####Fig.S9
#R.3.6.2#
library(vegan)

#Fig.S9a
otu <- read.csv("data/shotgun.function.table.csv")

bac.nmds<-metaMDS(vegdist(otu[1:96,4:15433],method = 'bray'),trymax=100)
bac.nmds
#####permanova
permanova<-adonis(otu[1:96,4:15433] ~ otu$treat)###perm anova### 
permanova

x=bac.nmds$points[,1]
y=bac.nmds$points[,2]
max(x)
min(x)
max(y)
min(y)

library(ggplot2)
df_points <- as.data.frame(bac.nmds$points)
df_points$samples <- otu$treat
df_points$Type <-c("AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN",
                   "AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN",
                   "AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU",
                   "AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU")
library(plyr)
df <-df_points[,c("MDS1", "MDS2", "Type")]
find_hull <- function(df_points) df_points[chull(df_points$MDS1, df_points$MDS2), ]
hulls <- ddply(df, "Type", find_hull )###border

#准备构建旋转后的坐标
#write.csv(df_points,"functional_potential.NMDS.csv")

####读入新的旋转后的坐标
df_points2 <- read.csv("data/functional_potential.NMDS.csv",header = T,row.names = 1)
df <-df_points2[,c("MDS1", "MDS2", "Type")]
find_hull <- function(df_points2) df_points2[chull(df_points2$MDS1, df_points2$MDS2), ]

hulls <- ddply(df, "Type", find_hull )###border
                   

#plot
p1<-ggplot(data=df_points2, aes(x= MDS1 , y= MDS2,group=Type))+ 
  geom_point(aes(colour=samples),size=6,alpha=0.7)+ 
  geom_polygon(data=hulls, alpha=.1,aes(fill = Type))+
  labs(x = "NMDS1", y = "NMDS2")+
  scale_colour_manual(values=c("#C7E2FF","lightblue1","#95E9FF","#5EB5FF","#339CFF","blue",
                                        "#FFE1FF","#FFC8FF","#FFAFFF","#FF63FF","#FF6371","red"))+
  scale_fill_manual(values=c("blue","red"))+
  ggtitle(expression('Stress= 0.037'~~~~'Treatment:'~italic(R)^2~'= 0.773***'))+
  theme_bw()+
  guides(color=guide_legend(title= "Treatment"))+
  theme(legend.position = 'left',legend.direction = 'vertical',
        legend.title = element_text(colour="black", size=13, face="bold"),
        legend.text = element_text(colour="black", size=12),
        axis.text=element_text(size=15,face="bold",colour="black"),
        axis.title=element_text(size=17,face="bold"),
        plot.title=element_text(face = "bold", colour = "black", size = 15))+
  geom_vline(xintercept = 0, color = 'gray34', linetype = 1, linewidth = 1.5) +
  geom_hline(yintercept = 0, color = 'gray34', linetype = 1, linewidth = 1.5) 

p1
ggsave("Fig.S9a.pdf", p1,height=5,width =6,limitsize = FALSE )


#Fig.S9b
library(EnhancedVolcano)
library(patchwork)

gene5<-read.csv("data/LOGFC.kegg.ANAU50.csv",row.names=1,head=T)

gene5[which(gene5$FDR < 0.05 & gene5$logFC <= -0.5),'sig'] <- 'Down'
gene5[which(gene5$FDR < 0.05 & gene5$logFC >=0.5),'sig'] <- 'Up'
gene5[which(gene5$FDR >= 0.05 | abs(gene5$logFC) < 0.5),'sig'] <- 'None'
sum(gene5$sig=="Down")
sum(gene5$sig=="Up")
sum(gene5$sig=="None")

keyvals5 <- ifelse(gene5$logFC < -0.5 & gene5$FDR <0.05,"royalblue1",
                   ifelse(gene5$logFC > 0.5 & gene5$FDR <0.05,"#F24F50","gray70"))

names(keyvals5)[keyvals5=="#F24F50"] <- "AU50 up"
names(keyvals5)[keyvals5=="royalblue1"] <- "AN50 up"
names(keyvals5)[keyvals5=="gray70"] <- "None"

p5 <- EnhancedVolcano(gene5,
                      lab = NA,
                      x = 'logFC',
                      y = 'FDR',
                      title = 'AU50 versus AN50',
                      subtitle = "KEGG: Total = 15554 variables",
                      caption = bquote("AN50 up = 2722 variables; AU50 up = 1311 variables"),
                      pCutoff = 0.05,
                      FCcutoff = 0.5,
                      pointSize = 2.0,
                      labSize = 6,
                      col=c('gray70', 'gray70', 'gray70', 'red3'),
                      colAlpha = 0.5,
                      ylim = c(0, 55),
                      xlim = c(-5,8),
                      legendPosition = 'right',
                      selectLab = row.names(gene5)[which(names(keyvals5) %in% c("Up","Down","None"))],colCustom = keyvals5)
p5
ggsave("Fig.S9b.pdf", p5,height=6,width =7,limitsize = FALSE )


#Fig.S9cd
#R4.0.4
library(rfPermute)
library(randomForest)
library(rfUtilities)

fun_env<-read.csv("data/random_forest.csv",head=T,row.names=1)

###NMDS1
fun.rf <- randomForest(fun_env$NMDS1_function~IN+pH+plantrichness+NH3g,
                       data=fun_env,num.cores = 4, nrep = 1000)##ntree=500 mtry = 6 default
rf.significance(fun.rf,fun_env[,c(7:10)])##change the combined table with different columns,计算显著性 the same result P=0.001,R2=0.399
fun.rf##the details of model

fun.rfP1 <- rfPermute( fun_env$NMDS1_function~IN+pH+plantrichness+NH3g,
                       data=fun_env,num.cores = 4, nrep = 1000)
importance(fun.rfP1)
plotImportance(fun.rfP1,imp.type = NULL,plot = TRUE,alpha = 0.05,ranks = TRUE,ylab = "Increase in MSE (%)",scale = TRUE,
               sig.only = FALSE)


###NMDS2
fun.rf2 <- randomForest(fun_env$NMDS2_function~IN+pH+plantrichness+NH3g,
                        data=fun_env,num.cores = 4, nrep = 1000)##ntree=500 mtry = 6 default
rf.significance(fun.rf2,fun_env[,c(7:10)])##change the combined table with different columns,计算显著性 the same result P=0.001,R2=0.399
fun.rf2##the details of model

fun.rfP12 <- rfPermute( fun_env$NMDS2_function~IN+pH+plantrichness+NH3g,
                        data=fun_env,num.cores = 4, nrep = 1000)
importance(fun.rfP12)
plotImportance(fun.rfP12,imp.type = NULL,plot = TRUE,alpha = 0.05,ranks = TRUE,ylab = "Increase in MSE (%)",scale = TRUE,
               sig.only = FALSE)


#Fig.S9e
#better R4.2.3
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

data<-read.table("data/ANAU50.enrichment.YAS.function.txt",header=TRUE,sep="\t")
#AN50 input
AN50_UP <- data[,c(1,2,4)]
ratio <- data.frame(cbind(AN50_UP,AN50_UP$AN50_up/AN50_UP$total))
dat <- ratio[,c(1,4)]
colnames(dat) <-c("Annotation","AN50_UP_ratio")
merAN<- dat

#AU50 input
AU50_UP <- data[,c(1,3,4)]
ratio2 <- data.frame(cbind(AU50_UP,AU50_UP$AU50_up/AU50_UP$total))
dat2 <- ratio2[,c(1,4)]
colnames(dat2) <- c("Annotation","AU50_UP_ratio")

merAU <- dat2

rm <- merge(merAN,merAU,by="Annotation")
rm.m <- melt(rm)
rm.m$variable <- factor(rm.m$variable, levels=c("AN50_UP_ratio","AU50_UP_ratio"))
dat <- rm.m[order(-rm.m[,3]),]
dat$Annotation <- factor(dat$Annotation, levels=c(
  "Bacterial chemotaxis(A)",
  "Carbon fixation pathways in prokaryotes(A)",
  "Metabolism of terpenoids and polyketides",
  "Degradation of aromatic compounds",
  "Carbon metabolism(Y)",
  "Carbohydrate metabolism(Y)",
  "Biosynthesis of secondary metabolites(S)",
  "Biosynthesis of amino acids(Y)",
  "Amino acid metabolism(Y)",
  "Nucleotide metabolism(Y)",
  "Glycan biosynthesis and metabolism(S)",
  "Biosynthesis of cofactors",
  "Metabolism of cofactors and vitamins",
  "Quorum sensing",
  "Biosynthesis of nucleotide sugars(Y)",
  "Biosynthesis of other secondary metabolites(S)",
  "Metabolism of other amino acids(Y)",
  "Energy metabolism(A)",
  "Signal transduction",
  "Microbial metabolism in diverse environments",
  "Fatty acid biosynthesis(Y)",
  "Sulfur metabolism(A)",
  "Lipopolysaccharide biosynthesis(S)",
  "Cell motility(A)",
  "Biofilm formation(S)",
  "Exopolysaccharide biosynthesis(S)",
  "Membrane transport(A)",
  "Nitrogen metabolism(A)",
  "Fatty acid metabolism(Y)",
  "Lipid metabolism(Y)",
  "Translation",
  "Cell growth and death",
  "Replication and repair(S)",
  "Transport and catabolism",
  "Folding, sorting and degradation",
  "Transcription"
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

combined <- read.table("data/ANAU50.enrichment.YAS.function.combined.txt",header=TRUE,sep="\t")
combined$AN50_UP_pvalue <- phyper(combined$AN50_up,combined$total,sum(combined$total)-combined$total,sum(combined$AN50_up),
                                  lower.tail = FALSE, log.p = FALSE)
combined$AU50_UP_pvalue <- phyper(combined$AU50_up,combined$total,sum(combined$total)-combined$total,sum(combined$AU50_up),
                                  lower.tail = FALSE, log.p = FALSE)
combined[which(combined$AN50_UP_pvalue < 0.05),'sigAN50'] <- 'Significant'
combined[which(combined$AU50_UP_pvalue < 0.05),'sigAU50'] <- 'Significant'


dat[which(dat$value <= 1 ),'sig'] <- 'Non-enriched'
dat[which(dat$value > 1),'sig'] <- 'enriched'

ggplot(dat, aes(x=value, y=Annotation,colour=Timepoint,shape=sig)) +
  scale_colour_manual(name="",values = c("royalblue1","#F24F50","#9B4F0F","#6FB98F"))+
  scale_shape_manual(name="",values=c(1,1,19,19))+facet_wrap(~SampleType,scales="free_x",ncol=2)+
  scale_shape_manual(values  = c(21,21), limits = c('Non-enriched', 'enriched')) +
  geom_point(size = 4,stroke = 1)+theme_bw()+
  geom_vline(aes(xintercept=1),colour="Black",size=1,linetype="dashed")+xlab("Enrichment ratio")+
  theme(axis.text.x=element_text(size=11,color="black",face="bold",angle=90), 
        axis.text.y=element_text(size=11,color="black",face="bold"), 
        axis.title=element_text(size=11,face="bold"),text=element_text(size=10))+
  xlim(0, 4)  +
  theme(legend.text = element_text(colour="black", size = 11, face = "bold"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) +
  theme(strip.text = element_text(size=11,face="bold"))


#Fig.S9f &Fig. S11
#R3.6.2
library(readxl)
library(psych)
library(tidyr)
library(ggplot2)
library(vegan)
library(colorRamps)
library(ggrepel)
library(grid)

ko.sel0 <-read.csv("data/function.table3_NMDS.csv",header=T,row.names=1)
ko.sel.ID0 <-  ko.sel0[, c(1:3)]
ko.sel <-ko.sel0[, c(-1:-3)]
ko.sel2 <- log(ko.sel,10)
ko.sel3 <- scale(ko.sel2)
ko.sel4 <- sqrt(ko.sel)
env.amp1 <-read.csv("data/random_forest.csv",header=T,row.names=1)
colnames(ko.sel) == env.amp1$name

env.tmp1<-env.amp1[,c("NMDS1_function","NMDS2_function")]

#Fig.S11
#env.tmp1<-env.amp1[,c("plantrichness","NH3g","pH","IN")]

env.tmp2 <- data.frame(t(ko.sel2))
spman.d12 = corr.test(env.tmp1, env.tmp2,use="pairwise",method="spearman",adjust="fdr",alpha=.05,ci=FALSE)

r<-data.frame(spman.d12$r)
r[is.na(r)] <-0
r$X<-row.names(r)

r.long <- gather(r, Y, r, 1 : ncol(r)-1,  factor_key=TRUE)
p<-data.frame(spman.d12$p)
p$X<-row.names(p)
p.long <- gather(p, Y, p, 1 : ncol(p)-1, factor_key=TRUE)
cor.out<-cbind(r.long,p.long$p)
cor.out$r <- round(as.numeric(cor.out$r), 2)
str(cor.out$r)

cor.out$X<- factor(cor.out$X, 
                   levels = c("NMDS2_function","NMDS1_function"))

cor.out$X<- factor(cor.out$X, 
                   labels = c("NMDS2", "NMDS1"))

#Fig.S11
#cor.out$X<- factor(cor.out$X, 
#levels = c("plantrichness","NH3g","pH","IN"))

#cor.out$X<- factor(cor.out$X, 
#labels = c("Plant richness","NH3 stress","pH","IN"))

merged_df <- merge(cor.out, ko.sel.ID0, by = "Y")
cor.out<-merged_df

library(splitstackshape)
cor.out$category <- factor(cor.out$category, levels = c("Amino acid metabolism(Y)","Metabolism of other amino acids(Y)","Nucleotide metabolism(Y)","Fatty acid biosynthesis(Y)", "Lipid metabolism(Y)","Cell growth and death(Y)",
                                                        "Biofilm formation(S)","Lipopolysaccharide biosynthesis(S)","Exopolysaccharide biosynthesis(S)","Glycan biosynthesis and metabolism(S)" ,"Biosynthesis of other secondary metabolites(S)",
                                                        "Cell motility(A)","Membrane transport(A)","Nitrogen metabolism(A)","Sulfur metabolism(A)","Carbon fixation pathways in prokaryotes(A)",
                                                        "Replication and repair(YS)","Carbohydrate metabolism(YA)","Polycyclic aromatic hydrocarbon degradation(AS)","Metabolism of terpenoids and polyketides(AS)"
                                                        ))

cor.out$category <- factor(cor.out$category,labels = c("Amino acid metabolism", "OAAM","Nucleotide M","FB" ,"Lipid metabolism","Cell growth and death",
                                                       "BF","LB","EB","GBM","OSM", 
                                                       "CM","Membrane transport","NM","SM","CP",
                                                       "Replication and repair(YS)","Carbohydrate metabolism(YA)","AD(AS)","MTP(AS)"))
p<- ggplot(cor.out, aes(Y,X)) +
  facet_grid(.~ category,  scales = "free", space = "free" )+
  geom_tile(aes(fill = r), size=0)+
  #geom_hline(yintercept = c(8.5,9.5), color = "yellow")+
  #scale_fill_gradient(guide = "legend", high='green', low='blue',name="rho")+
  scale_fill_gradient2(guide = "legend", high='red3',mid="white", low='dodgerblue3',name="rho")+
  theme(strip.text = element_text(size = 12,face="bold", color = c("white")),
        strip.background = element_rect(fill = "red"),
        panel.spacing = unit(0.2, "lines"),
        axis.title= element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y=element_text(colour="black", size=12, face="bold"))
p

g <- ggplot_gtable(ggplot_build(p))
strips <- which(grepl('strip-', g$layout$name))
fills <- c("green4","green4","green4","green4","green4","green4",
           "orange","orange","orange","orange","orange",
           "purple","purple","purple","purple","purple",
           "gray","gray","gray","gray")
for (i in seq_along(strips)) {
  k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  l <- which(grepl('titleGrob', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- fills[i]
}
plot(g)


ggsave("Fig.S9f.pdf", g,height=2,width =33,limitsize = FALSE )


#Fig.S10
#R3.6.2
library(ecodist)
library(pheatmap)
library(vegan)

data2 <- read.csv("data/heatmap.YAS.table2.csv",header = T,row.names = 1)

apply(data2,1,sd)

mycolors <- colorRampPalette(c("#2c7bb6","gray99","red3"))(50)
group_sample <- read.csv("data/heatmap.YAS.group.csv",header = T,row.names = 1)
gene_group <-read.csv("data/heatmap.YAS.group_gene2.csv",header = T,row.names = 1)
ann_colors = list(level=c(AN00="#C7E2FF",AN02="lightblue1",AN05="#95E9FF",AN10="#5EB5FF",AN20="#339CFF",AN50="blue",
                          AU00="#FFE1FF",AU02="#FFC8FF",AU05="#FFAFFF",AU10="#FF63FF",AU20="#FF6371",AU50="red"),
                  type=c(AN="royalblue1",AU="red4"),
                  
                  category =c("Amino acid metabolism(Y)"="green4","Metabolism of other amino acids(Y)"="green4",
                               "Fatty acid biosynthesis(Y)"="green4","Lipid metabolism(Y)"="green4",
                              "Nucleotide metabolism(Y)"="green4","Cell growth and death(Y)"="green4",
                              
                              "Biofilm formation(S)"="orange","Biosynthesis of other secondary metabolites(S)"="orange",
                              "Exopolysaccharide biosynthesis(S)"="orange","Lipopolysaccharide biosynthesis(S)"="orange","Glycan biosynthesis and metabolism(S)"="orange", 
                              
                              "Cell motility(A)"="purple","Membrane transport(A)"="purple","Nitrogen metabolism(A)"="purple","Sulfur metabolism(A)"="purple","Carbon fixation pathways in prokaryotes(A)"="purple",
                              "Carbohydrate metabolism(YA)"="gray", "Replication and repair(YS)"="gray", "Polycyclic aromatic hydrocarbon degradation(AS)"="gray","Metabolism of terpenoids and polyketides(AS)"="gray"))

pdf('mine/plot/Fig.S10.pdf', height =12, width = 28)
pheatmap(scale="column",as.matrix(data2),fontsize=10,fontsize_row=10,cex=1,cellwidth=1,cellheight =3,col=mycolors,
         annotation_col = gene_group, 
         annotation_row=group_sample,
         annotation_colors = ann_colors,
         border=FALSE,show_rownames=F,cluster_cols=F,cluster_rows=F,show_colnames = F,
         gaps_col =c(176,229,400,495,548,558,612,631,751,762,773,782,800,957,975,995,1126,1348,1397,1404))

dev.off()




