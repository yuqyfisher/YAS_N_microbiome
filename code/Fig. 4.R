####Fig.4
#R4.2.2
library(circlize)  
library(ComplexHeatmap)  
library(grid)  
library(graphics)

#Fig.4ab
pdf('Fig.4ab.pdf', height =18, width = 6)
par(mfrow = c(6,2),cex = 1.1,pin=c(0.5,0.5))

###AN00
rm(list=ls())
mat<-read.csv("data/AN00.HGT.phylum.csv",row.names = 1,header = T)
library(circlize)
circos.clear()

grid.col = c(A = "lightgoldenrod1", B = "violet", C = "mediumpurple1",D = "springgreen3", E = "paleturquoise1", 
             a = "lightgoldenrod1", b = "violet", c = "mediumpurple1",d = "purple4",e = "paleturquoise1")

circos.par(gap.after = c("A" = 6.07, "B" = 6.07, "C" = 6.07, "D" = 6.07, "E" = 6.07,
                         "a" = 6.07,"b" = 6.07, "c" = 6.07, "d" = 6.07, "e" = 6.07),clock.wise = FALSE)

chordDiagram(mat, order = c("A", "B", "C", "D", "E", 
                            "a", "b", "c", "d","e"),
             grid.col = grid.col,
             annotationTrackHeight = c(0.15, 0.05),
             transparency = 0.4,
             link.lwd = 2,
             link.lty =2)
title("AN00: a total of 354 HGT events",line = 0,cex.main=0.9)
circos.clear()


###AU00
rm(list=ls())
mat<-read.csv("data/AU00.HGT.phylum.csv",row.names = 1,header = T)
library(circlize)
circos.clear()

grid.col = c(A = "lightgoldenrod1", B = "violet", C = "mediumpurple1",D = "springgreen3",E = "lightpink2", "F" = "paleturquoise1", G="wheat1",
             a = "lightgoldenrod1", b = "mediumpurple1",c = "springgreen3",d = "lightpink2",e = "paleturquoise1")

circos.par(gap.after = c("A" = 5.88, "B" = 5.88, "C" = 5.88, "D" = 5.88, "E" = 5.88,"F"=5.88,"G"=5.88,
                         "a" = 5.88,"b" = 5.88, "c" = 5.88, "d" = 5.88, "e" = 5.88),clock.wise = FALSE)

chordDiagram(mat, order = c("A", "B", "C", "D", "E","F","G", 
                            "a", "b", "c", "d","e"),
             grid.col = grid.col,
             annotationTrackHeight = c(0.15, 0.05),
             transparency = 0.4,
             link.lwd = 2, link.lty = 2)
title("AU00: a total of 362 HGT events",line = 0,cex.main=0.9)
circos.clear()


###AN02
rm(list=ls())
mat<-read.csv("data/AN02.HGT.phylum.csv",row.names = 1,header = T)
library(circlize)
circos.clear()

grid.col = c(A = "lightgoldenrod1", B = "violet", C = "mediumpurple1",D = "lightpink2", E = "paleturquoise1", 
             a = "lightgoldenrod1", b = "violet", c = "mediumpurple1",d = "lightpink2",e = "paleturquoise1")

circos.par(gap.after = c("A" = 7.08, "B" = 7.08, "C" = 7.08, "D" = 7.08, "E" = 7.08,
                         "a" = 7.08,"b" = 7.08, "c" = 7.08, "d" = 7.08, "e" = 7.08),clock.wise = FALSE)

chordDiagram(mat, order = c("A", "B", "C", "D", "E", 
                            "a", "b", "c", "d","e"),
             grid.col = grid.col,
             annotationTrackHeight = c(0.15, 0.05),
             transparency = 0.4,
             link.lwd = 2, link.lty = 2)
title("AN02: a total of 340 HGT events",line = 0,cex.main=0.9)
circos.clear()


###AU02
rm(list=ls())
mat<-read.csv("data/AU02.HGT.phylum.csv",row.names = 1,header = T)
library(circlize)
circos.clear()

grid.col = c(A = "lightgoldenrod1", B = "violet", C = "mediumpurple1",D = "springgreen3", E = "paleturquoise1", 
             a = "lightgoldenrod1", b = "violet", c = "mediumpurple1",d = "springgreen3",e = "lightpink2",f = "paleturquoise1", g="wheat1")

circos.par(gap.after = c("A" = 7.04, "B" = 7.04, "C" = 7.04, "D" = 7.04, "E" = 7.04,
                         "a" = 7.04,"b" = 7.04, "c" = 7.04, "d" = 7.04, "e" = 7.04,"f"=7.04, "g"=7.04 ),clock.wise = FALSE)

chordDiagram(mat, order = c("A", "B", "C", "D", "E", 
                            "a", "b", "c", "d","e","f","g"),
             grid.col = grid.col,
             annotationTrackHeight = c(0.15, 0.05),
             transparency = 0.4,
             link.lwd = 2, link.lty = 2)
title("AU02: a total of 347 HGT events",line = 0,cex.main=0.9)
circos.clear()


###AN05
rm(list=ls())
mat<-read.csv("data/AN05.HGT.phylum.csv",row.names = 1,header = T)
library(circlize)
circos.clear()

grid.col = c(A = "lightgoldenrod1", B = "violet", C = "mediumpurple1",D= "purple4",E = "lightpink2", "F" = "paleturquoise1", G="wheat1" ,
             a = "lightgoldenrod1", b = "violet", c = "springgreen3",d = "lightpink2",e = "paleturquoise1")

circos.par(gap.after = c("A" = 4.21, "B" = 4.21, "C" = 4.21, "D" = 4.21, "E" = 4.21," F"=4.21, "G"=4.21,
                         "a" = 4.21,"b" = 4.21, "c" = 4.21, "d" = 4.21, "e" = 4.21),clock.wise = FALSE)

chordDiagram(mat, order = c("A", "B", "C", "D", "E", "F","G", 
                            "a", "b", "c", "d","e"),
             grid.col = grid.col,
             annotationTrackHeight = c(0.15, 0.05),#名称离圆弧的距离，以及圆弧的宽度
             transparency = 0.4,
             link.lwd = 2, link.lty = 2)
title("AN05: a total of 430 HGT events",line = 0,cex.main=0.9)
circos.clear()


###AU05
rm(list=ls())
mat<-read.csv("data/AU05.HGT.phylum.csv",row.names = 1,header = T)
library(circlize)
circos.clear()

grid.col = c(A= "lightskyblue",B = "lightgoldenrod1", C = "violet",D= "purple4",E = "springgreen3", "F" ="seagreen1", G= "paleturquoise1", H="wheat1" ,
             a= "lightskyblue", b = "lightgoldenrod1", c = "violet",d = "mediumpurple1", e= "purple4",f = "springgreen3",g = "paleturquoise1",h="wheat1")

circos.par(gap.after = c("A" = 8, "B" = 8, "C" = 8, "D" = 8, "E" = 8," F"=8, "G"=8,"H"=8,
                         "a" = 8,"b" = 8, "c" = 8, "d" = 8, "e" = 8, "f"=8 ,"g"=8,"h"=8),clock.wise = FALSE)

chordDiagram(mat, order = c("A", "B", "C", "D", "E", "F","G","H", 
                            "a", "b", "c", "d","e","f","g","h"),
             grid.col = grid.col,
             annotationTrackHeight = c(0.15, 0.05),
             transparency = 0.4,
             link.lwd = 2, link.lty = 2)
title("AU05: a total of 338 HGT events",line = 0,cex.main=0.9)
circos.clear()


###AN10
rm(list=ls())
mat<-read.csv("data/AN10.HGT.phylum.csv",row.names = 1,header = T)
library(circlize)
circos.clear()

grid.col = c(A= "lightskyblue",B = "lightgoldenrod1", C = "violet", D = "mediumpurple1",E= "springgreen3","F" = "lightpink2", G = "paleturquoise1", H="wheat1" ,
             a = "lightgoldenrod1", b = "violet",c = "mediumpurple1", d = "springgreen3",e = "seagreen1",f = "lightpink2",g = "paleturquoise1", h="wheat1")

circos.par(gap.after = c("A" = 3.62, "B" = 3.62, "C" = 3.62, "D" = 3.62, "E" = 3.62," F"=3.62, "G"=3.62,"H"=3.62,
                         "a" = 3.62,"b" = 3.62, "c" = 3.62, "d" = 3.62, "e" = 3.62,"f"=3.62, "g"=3.62,"h"=3.62),clock.wise = FALSE)

chordDiagram(mat, order = c("A", "B", "C", "D", "E", "F","G","H", 
                            "a", "b", "c", "d","e","f","g","h"),
             grid.col = grid.col,
             annotationTrackHeight = c(0.15, 0.05),
             transparency = 0.4,
             link.lwd = 2, link.lty = 2,)
title("AN10: a total of 588 HGT events",line = 0,cex.main=0.9)
circos.clear()


###AU10
rm(list=ls())
mat<-read.csv("data/AU10.HGT.phylum.csv",row.names = 1,header = T)
library(circlize)
circos.clear()

grid.col = c(A= "lightskyblue",B = "lightgoldenrod1", C = "mediumpurple1",D= "springgreen3",E = "lightpink2", "F" = "paleturquoise1",
             a = "lightgoldenrod1", b = "mediumpurple1", c = "springgreen3",d = "seagreen1",e = "lightpink2",f = "paleturquoise1")

circos.par(gap.after = c("A" = 4.30, "B" = 4.30, "C" =4.30, "D" = 4.30, "E" = 4.30," F"=4.30, 
                         "a" = 4.30,"b" = 4.30, "c" = 4.30, "d" = 4.30, "e" = 4.30,"f"=4.30),clock.wise = FALSE)

chordDiagram(mat, order = c("A", "B", "C", "D", "E", "F", 
                            "a", "b", "c", "d","e","f"),
             grid.col = grid.col,
             annotationTrackHeight = c(0.15, 0.05),
             transparency = 0.4,
             link.lwd = 2, link.lty = 2)
title("AU10: a total of 421 HGT events",line = 0,cex.main=0.9)
circos.clear()


###AN20
rm(list=ls())
mat<-read.csv("data/AN20.HGT.phylum.csv",row.names = 1,header = T)
library(circlize)
circos.clear()

grid.col = c(A= "lightskyblue",B = "lightgoldenrod1", C = "violet", D = "mediumpurple1",E= "springgreen3","F" ="paleturquoise1" ,G="wheat1",
             a= "lightskyblue",b = "lightgoldenrod1", c = "violet",d = "mediumpurple1", e = "springgreen3",f = "paleturquoise1", g="wheat1")

circos.par(gap.after = c("A" = 3.54, "B" = 3.54, "C" = 3.54, "D" = 3.54, "E" = 3.54,"F"=3.54, "G"=3.54,
                         "a" = 3.54,"b" = 3.54, "c" = 3.54, "d" = 3.54, "e" = 3.54,"f"=3.54, "g"=3.54),clock.wise = FALSE)

chordDiagram(mat, order = c("A", "B", "C", "D", "E", "F","G", 
                            "a", "b", "c", "d","e","f","g"),
             grid.col = grid.col,
             annotationTrackHeight = c(0.15, 0.05),#名称离圆弧的距离，以及圆弧的宽度
             transparency = 0.4,
             link.lwd = 2, link.lty = 2)
title("AN20: a total of 607 HGT events",line = 0,cex.main=0.9)
circos.clear()


###AU20
rm(list=ls())
mat<-read.csv("data/AU20.HGT.phylum.csv",row.names = 1,header = T)
library(circlize)
circos.clear()

grid.col = c(A= "lightskyblue",B = "lightgoldenrod1", C = "violet", D= "springgreen3",E="seagreen1","F" ="lightpink2",G="paleturquoise1" ,H="wheat1",
             a= "lightskyblue",b = "lightgoldenrod1", c = "violet",d = "mediumpurple1", e = "purple4",f= "springgreen3",g ="lightpink2",h = "paleturquoise1", i="magenta")

circos.par(gap.after = c("A" = 3.4, "B" = 3.4, "C" = 3.4, "D" = 3.4, "E" = 3.4,"F"=3.4, "G"=3.4,"H"=3.4,
                         "a" = 3.4,"b" = 3.4, "c" = 3.4, "d" = 3.4, "e" = 3.4,"f"=3.4, "g"=3.4, "h"=3.4,"i"=3.4),clock.wise = FALSE)

chordDiagram(mat, order = c("A", "B", "C", "D", "E", "F","G","H" ,
                            "a", "b", "c", "d","e","f","g","h","i"),
             grid.col = grid.col,
             annotationTrackHeight = c(0.15, 0.05),
             transparency = 0.4,
             link.lwd = 2, link.lty = 2)
title("AU20: a total of 708 HGT events",line = 0,cex.main=0.9)
circos.clear()


###AN50
rm(list=ls())
mat<-read.csv("data/AN50.HGT.phylum.csv",row.names = 1,header = T)
library(circlize)
circos.clear()

grid.col = c(A= "lightskyblue", B = "lightgoldenrod1", C = "springgreen3", D = "seagreen1",E = "paleturquoise1", "F" = "wheat1", 
             a= "lightskyblue", b = "lightgoldenrod1", c = "violet", d = "purple4",e = "springgreen3",f = "seagreen1",g = "paleturquoise1",h="mediumpurple1")

circos.par(gap.after = c("A" = 3.47, "B" = 3.47, "C" = 3.47, "D" = 3.47, "E" = 3.47,"F" = 3.47,
                         "a" = 3.47,"b" = 3.47, "c" = 3.47, "d" = 3.47, "e" = 3.47, "f"=3.47, "g"=3.47,"h"=3.47),clock.wise = FALSE)

chordDiagram(mat, order = c("A", "B", "C", "D", "E", "F",
                            "a", "b", "c", "d","e","f","g","h"),
             grid.col = grid.col,
             annotationTrackHeight = c(0.15, 0.05),
             transparency = 0.4,
             link.lwd = 2, link.lty = 2)
title("AN50: a total of 704 HGT events",line = 0,cex.main=0.9)
circos.clear()


###AU50
rm(list=ls())
mat<-read.csv("data/AU50.HGT.phylum.csv",row.names = 1,header = T)
library(circlize)
circos.clear()

grid.col = c(A= "lightskyblue", B = "lightgoldenrod1", C = "violet",D="mediumpurple1",E="red", "F" = "seagreen1",G = "paleturquoise1", H = "wheat1", 
             a= "lightskyblue", b = "lightgoldenrod1", c = "violet", d = "mediumpurple1",e = "springgreen3",f = "seagreen1",g = "paleturquoise1",h="mediumpurple1")

circos.par(gap.after = c("A" = 3.38, "B" = 3.38, "C" = 3.38, "D" = 3.38, "E" = 3.38,"F" = 3.38,"G"=3.38,"H"=3.38,
                         "a" = 3.38,"b" = 3.38, "c" = 3.38, "d" = 3.38, "e" = 3.38, "f"=3.38, "g"=3.38,"h"=3.38),clock.wise = FALSE)

chordDiagram(mat, order = c("A", "B", "C", "D", "E", "F","G","H",
                            "a", "b", "c", "d","e","f","g","h"),
             grid.col = grid.col,
             annotationTrackHeight = c(0.15, 0.05),#名称离圆弧的距离，以及圆弧的宽度
             transparency = 0.4,
             link.lwd = 2, link.lty = 2)
title("AU50: a total of 1297 HGT events",line = 0,cex.main=0.9)
circos.clear()
dev.off()


###############Figure legend###########################
##donor
pdf('Fig.4ab.legend.pdf', height =18, width = 6)

legend11 <- Legend(
  at = c(1, 2, 3, 4, 5, 6, 7, 8,9,10,11,12),
  labels = c('Acidobacteriota', 'Actinobacteriota',"Bacteroidota","Chloroflexota","Chrysiogenetes", "Deinococcus_Thermus","Firmicutes",
             "Gemmatimonadota","Planctomycetota","Proteobacteria","Spirochaetes", "I_Verrucomicrobiota"),
  labels_gp = gpar(fontsize = 11),  title_gp = gpar(fontsize = 13),title = 'Phylum',#title_position = "leftcenter",
  grid_height = unit(0.5, 'cm'), grid_width = unit(0.5, 'cm'), type = 'points',background = NA,
  legend_gp = gpar(col= c("lightskyblue", "lightgoldenrod1","violet","mediumpurple1","red","purple4","springgreen3",
                          "seagreen1","lightpink2","paleturquoise1","magenta","wheat1")),
  ncol = 1,
  gap = unit(0.45, "cm"),
  pch = c(15, 15, 15, 15, 15, 15, 15, 15,15,15,15,15) ) 


pushViewport(viewport(x = 0.5, y = 0.5))
lgd_list_vertical2 <- packLegend(legend11)
grid.draw(lgd_list_vertical2)

upViewport()

dev.off()
circos.clear()


#Fig.4c
#R3.6.2
library(ecodist)
library(pheatmap)
library(vegan)

data2 <- read.csv("data/HGT.KO.csv",header = T,row.names = 1)

apply(data2,1,sd)

mycolors <- colorRampPalette(c("yellow","gray99","purple"))(50)

group <- read.csv("data/heatmap.YAS.group.csv",header = T,row.names = 1)

ann_colors = list(level=c(AN00="#C7E2FF",AN02="lightblue1",AN05="#95E9FF",AN10="#5EB5FF",AN20="#339CFF",AN50="blue",
                          AU00="#FFE1FF",AU02="#FFC8FF",AU05="#FFAFFF",AU10="#FF63FF",AU20="#FF6371",AU50="red"),type=c(AN="royalblue1",AU="red4"))

pdf('Fig.4c.pdf', height =17, width = 10)

pheatmap(scale="row",as.matrix(data2),fontsize=10,fontsize_row=10,cex=1,cellwidth=5,cellheight =0.5,col=mycolors,
         annotation_col = group,annotation_colors = ann_colors,
         border=FALSE,show_rownames=F,cluster_cols=F,cluster_rows=F)

dev.off()


#Fig.4d
#R.3.6.2#
library(ggplot2)
library(vegan)
library(colorRamps)
library(psych)
library(tidyr)
library(nlme)
library(MuMIn)
library(ggbreak)

data<-read.csv("data/HGT.KO.richness.csv",row.names=1,head=T)
data.AN <- data[1:48,]
data.AU <- data[49:96,]

p4 <-ggplot()+
  geom_jitter(data=data, aes(y=Richness, x= log(level3+1), color=type1),size=5,alpha=0.8)+
  geom_smooth(data=subset(data,data$type=="AN"),aes(x= log(level3+1), y=Richness,color=type1), method ="lm", span=0.9)+
  geom_smooth(data=subset(data,data$type=="AU"),aes(x= log(level3+1), y=Richness,color=type1), method ="lm", span=0.9)+
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('ln(N addition rate ('*g~N~m ^-2~yr^-1*'))')) +
  ylab(bquote('KO richness of HGT'))+
  scale_color_manual(values = c( "royalblue1","#F24F50", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=16,face="bold"),
        axis.text.x =element_text(color="black",size = 14),
        axis.text.y =element_text(color="black",size = 14),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 12),
        plot.title=element_text(face = "bold", colour = "black", size = 18),
        legend.title = element_text(face = "bold", colour = "black", size = 15))+ 
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R)["AN"] ^2~'= 0.410***'~~~~~italic(R)["AU"] ^2~'= 0.503***'))

p4
summary(lm(data.AN$Richness~log(data.AN$level3+1)))
summary(lm(data.AU$Richness~log(data.AU$level3+1)))

ggsave("Fig.4d.pdf", p4,height=5,width =5.8,limitsize = FALSE )


#Fig.4e
#R3.6.2
library(EnhancedVolcano)
library(patchwork)
library(edgeR)
library(statmod)

gene5<-read.csv("data/AN50AU50.HGT.preference.csv",row.names=1,head=T)

gene5[which(gene5$FDR < 0.05 & gene5$AU50 <= -1),'sig'] <- 'Down'
gene5[which(gene5$FDR < 0.05 & gene5$AU50 >=1),'sig'] <- 'Up'
gene5[which(gene5$FDR >= 0.05 | abs(gene5$AU50) < 1),'sig'] <- 'None'

sum(gene5$sig=='Down')
sum(gene5$sig=="Up")
sum(gene5$sig=="None")

p1<-ggplot(data=gene5, aes(x= AU50 , y= FDR,group=sig))+ 
  geom_point(aes(colour=sig),size=4,alpha=0.5)+ 
  labs(x = "2DP", y = "P value (FDR)")+
  scale_colour_manual(values=c("red","blue","gray70","#5EB5FF","#339CFF","blue",
                               "#FFE1FF","#FFC8FF","#FFAFFF","#FF63FF","#FF6371","red"),
                      breaks=c("Up","Down", "None"),
                      labels=c("AU50 preference", "AN50 preference","None"))+
  
  theme_bw()+
  ggtitle(expression('AU50 versus AN50'),subtitle = "HGT genes (KO): Total = 1219 variables")+
  theme(legend.position = 'left',legend.direction = 'vertical',
        legend.title = element_text(colour="black", size=13, face="bold"),
        legend.text = element_text(colour="black", size=12),
        axis.text=element_text(size=15,face="bold",colour="black"),
        axis.title=element_text(size=17,face="bold"))+
  geom_vline(xintercept = 1, color = 'gray34', linetype = 2, linewidth = 1.3) +
  geom_vline(xintercept = -1, color = 'gray34', linetype = 2, linewidth = 1.3) +
  geom_hline(yintercept = 0.05, color = 'gray34', linetype = 2, linewidth = 1.3) 

p1
ggsave("Fig.4e.pdf", p1,height=4.5,width =6,limitsize = FALSE )



#Fig.4f
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

data<-read.table("data/HGT.preference.enrichment.AN50AU50.txt",header=TRUE,sep="\t")

#AN50
AN50_UP <- data[,c(1,2,4)]
ratio <- data.frame(cbind(AN50_UP,AN50_UP$AN50_up/AN50_UP$total))
dat <- ratio[,c(1,4)]
colnames(dat) <-c("Annotation","AN50_UP_ratio")
merAN<- dat

#AU50 
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
  "Exopolysaccharide biosynthesis(S)",
  "Glycan biosynthesis and metabolism(S)",
  "Biosynthesis of nucleotide sugars(Y)",
  "Cell growth and death",
  "Amino acid metabolism(Y)",
  "Biofilm formation(S)",
  "Translation",
  "Sulfur metabolism(A)",
  "Fatty acid biosynthesis(Y)",
  "Biosynthesis of other secondary metabolites(S)",
  "Lipopolysaccharide biosynthesis(S)",
  "Bacterial chemotaxis(A)",
  "Degradation of aromatic compounds",
  "Metabolism of terpenoids and polyketides",
  "Folding, sorting and degradation",
  "Membrane transport(A)",
  "Biosynthesis of secondary metabolites(S)",
  "Carbohydrate metabolism(Y)",
  "Microbial metabolism in diverse environments",
  "Signal transduction",
  "Carbon metabolism(Y)",
  "Biosynthesis of amino acids(Y)",
  "Carbon fixation pathways in prokaryotes(A)",
  "Energy metabolism(A)",
  "Lipid metabolism(Y)",
  "Quorum sensing",
  "Metabolism of other amino acids(Y)",
  "Fatty acid metabolism(Y)",
  "Nitrogen metabolism(A)",
  "Transcription",
  "Biosynthesis of cofactors",
  "Metabolism of cofactors and vitamins",
  "Nucleotide metabolism(Y)",
  "Replication and repair(S)",
  "Cell motility(A)"))

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

combined <- read.table("data/HGT.preference.enrichment.AN50AU50.combined.txt",header=TRUE,sep="\t")
combined$AN50_UP_pvalue <- phyper(combined$AN50_up-1,combined$total,sum(combined$total)-combined$total,sum(combined$AN50_up),
                                  lower.tail = FALSE, log.p = FALSE)

combined$AU50_UP_pvalue <- phyper(combined$AU50_up-1,combined$total,sum(combined$total)-combined$total,sum(combined$AU50_up),
                                  lower.tail = FALSE, log.p = FALSE)

combined[which(combined$AN50_UP_pvalue < 0.05),'sig_AN'] <- 'Significant'
combined[which(combined$AU50_UP_pvalue < 0.05),'sig_AU'] <- 'Significant'


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


