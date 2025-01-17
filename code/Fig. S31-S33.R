#Fig.S31-S33

####Fig.S31
#R.3.6.2#
library(ggplot2)
library(vegan)
library(colorRamps)
library(psych)
library(tidyr)
library(nlme)
library(MuMIn)
library(ggbreak)
library(patchwork)

rm(list=ls())

#######Genome size
#ËáÏ¸¾ú
rm(list=ls())
env.amp.N <- read.csv("Genome_size_Aci.csv",head=T,row.names = 1)

p1 <-ggplot(data=env.amp.N,aes(y=corrected_size, x= Type,fill=Type))+
  geom_boxplot(outlier.shape = NA,size=0.8, width = 0.5,alpha=0.2,aes(fill=Type),position = position_dodge(0.5),) +
  
  geom_jitter(aes(Type, corrected_size,color = Type),position=position_jitter(width=0.01,height=0),
              alpha=0.6,stroke=2,shape=1, size = 3)+
  
  theme_bw()+
  xlab(bquote('Nitrogen type ')) +
  ylab(bquote('Genome size (Mbp)'))+
  scale_fill_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_color_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_y_continuous(limits = c(2, 9))+
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 18),
        axis.text.y =element_text(color="black",size = 18),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 12),
        plot.title=element_text(face = "bold", colour = "black", size = 16),
        legend.title = element_text(face = "bold", colour = "black", size = 14))+ 
  
  labs(color="Nitrogen type")+
  ggtitle(expression(~Acidobacteriota~~~italic(P)~'= 0.800'))
p1

#·ÅÏß¾ú
env.amp.N <- read.csv("Genome_size_Actino.csv",head=T,row.names = 1)

p2 <-ggplot(data=env.amp.N,aes(y=corrected_size, x= Type,fill=Type))+
  geom_boxplot(outlier.shape = NA,size=0.8, width = 0.5,alpha=0.2,aes(fill=Type),position = position_dodge(0.5),) +
  
  geom_jitter(aes(Type, corrected_size,color = Type),position=position_jitter(width=0.01,height=0),
              alpha=0.6,stroke=2,shape=1, size = 3)+
  
  theme_bw()+
  xlab(bquote('Nitrogen type ')) +
  ylab(bquote('Genome size (Mbp)'))+
  scale_fill_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_color_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_y_continuous(limits = c(2, 9))+
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 18),
        axis.text.y =element_text(color="black",size = 18),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 12),
        plot.title=element_text(face = "bold", colour = "black", size = 16),
        legend.title = element_text(face = "bold", colour = "black", size = 14))+ 
  
  labs(color="Nitrogen type")+
  ggtitle(expression(~Actinobacteriota~~~italic(P)~'= 0.879'))
p2


#Äâ¸Ë¾ú
env.amp.N <- read.csv("Genome_size_Bact.csv",head=T,row.names = 1)

p3 <-ggplot(data=env.amp.N,aes(y=corrected_size, x= Type,fill=Type))+
  geom_boxplot(outlier.shape = NA,size=0.8, width = 0.5,alpha=0.2,aes(fill=Type),position = position_dodge(0.5),) +
  
  geom_jitter(aes(Type, corrected_size,color = Type),position=position_jitter(width=0.01,height=0),
              alpha=0.6,stroke=2,shape=1, size = 3)+
  
  theme_bw()+
  xlab(bquote('Nitrogen type ')) +
  ylab(bquote('Genome size (Mbp)'))+
  scale_fill_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_color_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_y_continuous(limits = c(2, 9))+
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 18),
        axis.text.y =element_text(color="black",size = 18),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 12),
        plot.title=element_text(face = "bold", colour = "black", size = 16),
        legend.title = element_text(face = "bold", colour = "black", size = 14))+ 
  
  labs(color="Nitrogen type")+
  ggtitle(expression(~Bacteroidota~~~italic(P)~'= 0.063'))
p3



#GEMM
env.amp.N <- read.csv("Genome_size_GEMM.csv",head=T,row.names = 1)


p4 <-ggplot(data=env.amp.N,aes(y=corrected_size, x= Type,fill=Type))+
  geom_boxplot(outlier.shape = NA,size=0.8, width = 0.5,alpha=0.2,aes(fill=Type),position = position_dodge(0.5),) +
  
  geom_jitter(aes(Type, corrected_size,color = Type),position=position_jitter(width=0.01,height=0),
              alpha=0.6,stroke=2,shape=1, size = 3)+
  
  theme_bw()+
  xlab(bquote('Nitrogen type ')) +
  ylab(bquote('Genome size (Mbp)'))+
  scale_fill_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_color_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_y_continuous(limits = c(2, 9))+
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 18),
        axis.text.y =element_text(color="black",size = 18),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 12),
        plot.title=element_text(face = "bold", colour = "black", size = 16),
        legend.title = element_text(face = "bold", colour = "black", size = 14))+ 
  
  labs(color="Nitrogen type")+
  ggtitle(expression(~Gemmatimonadota~~~italic(P)~'= 0.905'))
p4


#Plancto

env.amp.N <- read.csv("Genome_size_Plancto.csv",head=T,row.names = 1)

p5 <-ggplot(data=env.amp.N,aes(y=corrected_size, x= Type,fill=Type))+
  geom_boxplot(outlier.shape = NA,size=0.8, width = 0.5,alpha=0.2,aes(fill=Type),position = position_dodge(0.5),) +
  
  geom_jitter(aes(Type, corrected_size,color = Type),position=position_jitter(width=0.01,height=0),
              alpha=0.6,stroke=2,shape=1, size = 3)+
  
  theme_bw()+
  xlab(bquote('Nitrogen type ')) +
  ylab(bquote('Genome size (Mbp)'))+
  scale_fill_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_color_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_y_continuous(limits = c(2, 9))+
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 18),
        axis.text.y =element_text(color="black",size = 18),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 12),
        plot.title=element_text(face = "bold", colour = "black", size = 16),
        legend.title = element_text(face = "bold", colour = "black", size = 14))+ 
  
  labs(color="Nitrogen type")+
  ggtitle(expression(~Planctomycetota~~~italic(P)~'= 1.000'))
p5


#±äÐÎ¾ú

env.amp.N <- read.csv("Genome_size_Pro.csv",head=T,row.names = 1)


p6 <-ggplot(data=env.amp.N,aes(y=corrected_size, x= Type,fill=Type))+
  geom_boxplot(outlier.shape = NA,size=0.8, width = 0.5,alpha=0.2,aes(fill=Type),position = position_dodge(0.5),) +
  
  geom_jitter(aes(Type, corrected_size,color = Type),position=position_jitter(width=0.01,height=0),
              alpha=0.6,stroke=2,shape=1, size = 3)+
  
  theme_bw()+
  xlab(bquote('Nitrogen type ')) +
  ylab(bquote('Genome size (Mbp)'))+
  scale_fill_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_color_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_y_continuous(limits = c(2, 9))+
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 18),
        axis.text.y =element_text(color="black",size = 18),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 12),
        plot.title=element_text(face = "bold", colour = "black", size = 16),
        legend.title = element_text(face = "bold", colour = "black", size = 14))+ 
  
  labs(color="Nitrogen type")+
  ggtitle(expression(~Proteobacteria~~~italic(P)~'= 0.799'))
p6



p11 <- 
  p1 + p2 + p3 + p4 + p5 + p6+  
  plot_layout(guides = "collect",ncol=6,nrow=1) +
  theme(legend.position='right')

p11




####GCº¬Á¿
#ËáÏ¸¾ú
rm(list=ls())
env.amp.N <- read.csv("Genome_size_Aci.csv",head=T,row.names = 1)


p1 <-ggplot(data=env.amp.N,aes(y=GC_percent, x= Type,fill=Type))+
  geom_boxplot(outlier.shape = NA,size=0.8, width = 0.5,alpha=0.2,aes(fill=Type),position = position_dodge(0.5),) +
  
  geom_jitter(aes(Type, GC_percent,color = Type),position=position_jitter(width=0.01,height=0),
              alpha=0.6,stroke=2,shape=1, size = 3)+
  
  theme_bw()+
  xlab(bquote('Nitrogen type ')) +
  ylab(bquote('GC content (%)'))+
  scale_fill_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_color_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_y_continuous(labels = scales::percent,limits = c(0.30, 0.75))+
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 18),
        axis.text.y =element_text(color="black",size = 18),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 12),
        plot.title=element_text(face = "bold", colour = "black", size = 16),
        legend.title = element_text(face = "bold", colour = "black", size = 14))+ 
  
  labs(color="Nitrogen type")+
  ggtitle(expression(~Acidobacteriota~~~italic(P)~'= 0.400'))
p1


#·ÅÏß¾ú
env.amp.N <- read.csv("Genome_size_Actino.csv",head=T,row.names = 1)

p2 <-ggplot(data=env.amp.N,aes(y=GC_percent, x= Type,fill=Type))+
  geom_boxplot(outlier.shape = NA,size=0.8, width = 0.5,alpha=0.2,aes(fill=Type),position = position_dodge(0.5),) +
  
  geom_jitter(aes(Type, GC_percent,color = Type),position=position_jitter(width=0.01,height=0),
              alpha=0.6,stroke=2,shape=1, size = 3)+
  
  theme_bw()+
  xlab(bquote('Nitrogen type ')) +
  ylab(bquote('GC content (%)'))+
  scale_fill_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_color_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_y_continuous(labels = scales::percent,limits = c(0.30, 0.75))+
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 18),
        axis.text.y =element_text(color="black",size = 18),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 12),
        plot.title=element_text(face = "bold", colour = "black", size = 16),
        legend.title = element_text(face = "bold", colour = "black", size = 14))+ 
  
  labs(color="Nitrogen type")+
  ggtitle(expression(~Actinobacteriota~~~italic(P)~'= 0.371'))
p2



#Äâ¸Ë¾ú
env.amp.N <- read.csv("Genome_size_Bact.csv",head=T,row.names = 1)

p3 <-ggplot(data=env.amp.N,aes(y=GC_percent, x= Type,fill=Type))+
  geom_boxplot(outlier.shape = NA,size=0.8, width = 0.5,alpha=0.2,aes(fill=Type),position = position_dodge(0.5),) +
  
  geom_jitter(aes(Type, GC_percent,color = Type),position=position_jitter(width=0.01,height=0),
              alpha=0.6,stroke=2,shape=1, size = 3)+
  
  theme_bw()+
  xlab(bquote('Nitrogen type ')) +
  ylab(bquote('GC content (%)'))+
  scale_fill_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_color_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_y_continuous(labels = scales::percent,limits = c(0.30, 0.75))+
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 18),
        axis.text.y =element_text(color="black",size = 18),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 12),
        plot.title=element_text(face = "bold", colour = "black", size = 16),
        legend.title = element_text(face = "bold", colour = "black", size = 14))+ 
  
  labs(color="Nitrogen type")+
  ggtitle(expression(~Bacteroidota~~~italic(P)~'= 0.659'))
p3



#GEMM
env.amp.N <- read.csv("Genome_size_GEMM.csv",head=T,row.names = 1)

p4 <-ggplot(data=env.amp.N,aes(y=GC_percent, x= Type,fill=Type))+
  geom_boxplot(outlier.shape = NA,size=0.8, width = 0.5,alpha=0.2,aes(fill=Type),position = position_dodge(0.5),) +
  
  geom_jitter(aes(Type, GC_percent,color = Type),position=position_jitter(width=0.01,height=0),
              alpha=0.6,stroke=2,shape=1, size = 3)+
  
  theme_bw()+
  xlab(bquote('Nitrogen type ')) +
  ylab(bquote('GC content (%)'))+
  scale_fill_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_color_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_y_continuous(labels = scales::percent,limits = c(0.30, 0.75))+
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 18),
        axis.text.y =element_text(color="black",size = 18),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 12),
        plot.title=element_text(face = "bold", colour = "black", size = 16),
        legend.title = element_text(face = "bold", colour = "black", size = 14))+ 
  
  labs(color="Nitrogen type")+
  ggtitle(expression(~Gemmatimonadota~~~italic(P)~'= 0.225'))
p4


#Plancto
env.amp.N <- read.csv("Genome_size_Plancto.csv",head=T,row.names = 1)


p5 <-ggplot(data=env.amp.N,aes(y=GC_percent, x= Type,fill=Type))+
  geom_boxplot(outlier.shape = NA,size=0.8, width = 0.5,alpha=0.2,aes(fill=Type),position = position_dodge(0.5),) +
  
  geom_jitter(aes(Type, GC_percent,color = Type),position=position_jitter(width=0.01,height=0),
              alpha=0.6,stroke=2,shape=1, size = 3)+
  
  theme_bw()+
  xlab(bquote('Nitrogen type ')) +
  ylab(bquote('GC content (%)'))+
  scale_fill_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_color_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_y_continuous(labels = scales::percent,limits = c(0.30, 0.75))+
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 18),
        axis.text.y =element_text(color="black",size = 18),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 12),
        plot.title=element_text(face = "bold", colour = "black", size = 16),
        legend.title = element_text(face = "bold", colour = "black", size = 14))+ 
  
  labs(color="Nitrogen type")+
  ggtitle(expression(~Planctomycetota~~~italic(P)~'= 1.000'))
p5


#±äÐÎ¾ú
env.amp.N <- read.csv("Genome_size_Pro.csv",head=T,row.names = 1)


p6 <-ggplot(data=env.amp.N,aes(y=GC_percent, x= Type,fill=Type))+
  geom_boxplot(outlier.shape = NA,size=0.8, width = 0.5,alpha=0.2,aes(fill=Type),position = position_dodge(0.5),) +
  
  geom_jitter(aes(Type, GC_percent,color = Type),position=position_jitter(width=0.01,height=0),
              alpha=0.6,stroke=2,shape=1, size = 3)+
  
  theme_bw()+
  xlab(bquote('Nitrogen type ')) +
  ylab(bquote('GC content (%)'))+
  scale_fill_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_color_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_y_continuous(labels = scales::percent,limits = c(0.30, 0.75))+
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 18),
        axis.text.y =element_text(color="black",size = 18),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 12),
        plot.title=element_text(face = "bold", colour = "black", size = 16),
        legend.title = element_text(face = "bold", colour = "black", size = 14))+ 
  
  labs(color="Nitrogen type")+
  ggtitle(expression(~Proteobacteria~~~italic(P)~'= 0.846'))
p6



p12 <- 
  p1 + p2 + p3 + p4 + p5 + p6+  
  plot_layout(guides = "collect",ncol=6,nrow=1) +
  theme(legend.position='right')
p12



####S.ARG
#ËáÏ¸¾ú
rm(list=ls())
env.amp.N <- read.csv("S.ARG.Aci.csv",head=T,row.names = 1)


p1 <-ggplot(data=env.amp.N,aes(y=S.ARG, x= Type,fill=Type))+
  geom_boxplot(outlier.shape = NA,size=0.8, width = 0.5,alpha=0.2,aes(fill=Type),position = position_dodge(0.5),) +
  
  geom_jitter(aes(Type, S.ARG,color = Type),position=position_jitter(width=0.01,height=0),
              alpha=0.6,stroke=2,shape=1, size = 3)+
  
  theme_bw()+
  xlab(bquote('Nitrogen type ')) +
  ylab(bquote('S.ARG'))+
  scale_fill_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_color_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_y_continuous(limits = c(3, 6))+
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 18),
        axis.text.y =element_text(color="black",size = 18),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 12),
        plot.title=element_text(face = "bold", colour = "black", size = 16),
        legend.title = element_text(face = "bold", colour = "black", size = 14))+ 
  
  labs(color="Nitrogen type")+
  ggtitle(expression(~Acidobacteriota~~~italic(P)~'= 0.500'))
p1



#·ÅÏß¾ú
env.amp.N <- read.csv("S.ARG.Actino.csv",head=T,row.names = 1)


p2 <-ggplot(data=env.amp.N,aes(y=S.ARG, x= Type,fill=Type))+
  geom_boxplot(outlier.shape = NA,size=0.8, width = 0.5,alpha=0.2,aes(fill=Type),position = position_dodge(0.5),) +
  
  geom_jitter(aes(Type, S.ARG,color = Type),position=position_jitter(width=0.01,height=0),
              alpha=0.6,stroke=2,shape=1, size = 3)+
  
  theme_bw()+
  xlab(bquote('Nitrogen type ')) +
  ylab(bquote('S.ARG'))+
  scale_fill_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_color_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_y_continuous(limits = c(3, 6))+
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 18),
        axis.text.y =element_text(color="black",size = 18),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 12),
        plot.title=element_text(face = "bold", colour = "black", size = 16),
        legend.title = element_text(face = "bold", colour = "black", size = 14))+ 
  
  labs(color="Nitrogen type")+
  ggtitle(expression(~Actinobacteriota~~~italic(P)~'= 0.295'))
p2


#Äâ¸Ë¾ú
env.amp.N <- read.csv("S.ARG.Bact.csv",head=T,row.names = 1)

p3 <-ggplot(data=env.amp.N,aes(y=S.ARG, x= Type,fill=Type))+
  geom_boxplot(outlier.shape = NA,size=0.8, width = 0.5,alpha=0.2,aes(fill=Type),position = position_dodge(0.5),) +
  
  geom_jitter(aes(Type, S.ARG,color = Type),position=position_jitter(width=0.01,height=0),
              alpha=0.6,stroke=2,shape=1, size = 3)+
  
  theme_bw()+
  xlab(bquote('Nitrogen type ')) +
  ylab(bquote('S.ARG'))+
  scale_fill_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_color_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_y_continuous(limits = c(3, 6))+
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 18),
        axis.text.y =element_text(color="black",size = 18),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 12),
        plot.title=element_text(face = "bold", colour = "black", size = 16),
        legend.title = element_text(face = "bold", colour = "black", size = 14))+ 
  
  labs(color="Nitrogen type")+
  ggtitle(expression(~Bacteroidota~~~italic(P)~'= 0.331'))
p3




#GEMM
env.amp.N <- read.csv("S.ARG.GEMM.csv",head=T,row.names = 1)


p4 <-ggplot(data=env.amp.N,aes(y=S.ARG, x= Type,fill=Type))+
  geom_boxplot(outlier.shape = NA,size=0.8, width = 0.5,alpha=0.2,aes(fill=Type),position = position_dodge(0.5),) +
  
  geom_jitter(aes(Type, S.ARG,color = Type),position=position_jitter(width=0.01,height=0),
              alpha=0.6,stroke=2,shape=1, size = 3)+
  
  theme_bw()+
  xlab(bquote('Nitrogen type ')) +
  ylab(bquote('S.ARG'))+
  scale_fill_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_color_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_y_continuous(limits = c(3, 6))+
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 18),
        axis.text.y =element_text(color="black",size = 18),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 12),
        plot.title=element_text(face = "bold", colour = "black", size = 16),
        legend.title = element_text(face = "bold", colour = "black", size = 14))+ 
  
  labs(color="Nitrogen type")+
  ggtitle(expression(~Gemmatimonadota~~~italic(P)~'= 0.925'))
p4



#Plancto
env.amp.N <- read.csv("S.ARG.Plancto.csv",head=T,row.names = 1)


p5 <-ggplot(data=env.amp.N,aes(y=S.ARG, x= Type,fill=Type))+
  geom_boxplot(outlier.shape = NA,size=0.8, width = 0.5,alpha=0.2,aes(fill=Type),position = position_dodge(0.5),) +
  
  geom_jitter(aes(Type, S.ARG,color = Type),position=position_jitter(width=0.01,height=0),
              alpha=0.6,stroke=2,shape=1, size = 3)+
  
  theme_bw()+
  xlab(bquote('Nitrogen type ')) +
  ylab(bquote('S.ARG'))+
  scale_fill_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_color_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_y_continuous(limits = c(3, 6))+
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 18),
        axis.text.y =element_text(color="black",size = 18),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 12),
        plot.title=element_text(face = "bold", colour = "black", size = 16),
        legend.title = element_text(face = "bold", colour = "black", size = 14))+ 
  
  labs(color="Nitrogen type")+
  ggtitle(expression(~Planctomycetota~~~italic(P)~'= 1.000'))
p5


#±äÐÎ¾ú
env.amp.N <- read.csv("S.ARG.Pro.csv",head=T,row.names = 1)


p6 <-ggplot(data=env.amp.N,aes(y=S.ARG, x= Type,fill=Type))+
  geom_boxplot(outlier.shape = NA,size=0.8, width = 0.5,alpha=0.2,aes(fill=Type),position = position_dodge(0.5),) +
  
  geom_jitter(aes(Type, S.ARG,color = Type),position=position_jitter(width=0.01,height=0),
              alpha=0.6,stroke=2,shape=1, size = 3)+
  
  theme_bw()+
  xlab(bquote('Nitrogen type ')) +
  ylab(bquote('S.ARG'))+
  scale_fill_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_color_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_y_continuous(limits = c(3, 6))+
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 18),
        axis.text.y =element_text(color="black",size = 18),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 12),
        plot.title=element_text(face = "bold", colour = "black", size = 16),
        legend.title = element_text(face = "bold", colour = "black", size = 14))+ 
  
  labs(color="Nitrogen type")+
  ggtitle(expression(~Proteobacteria~~~italic(P)~'= 0.868'))
p6


p13 <- 
  p1 + p2 + p3 + p4 + p5 + p6+  
  plot_layout(guides = "collect",ncol=6,nrow=1) +
  theme(legend.position='right')

p13




####S.KO
#ËáÏ¸¾ú
rm(list=ls())
env.amp.N <- read.csv("S.KO.Aci.csv",head=T,row.names = 1)


p1 <-ggplot(data=env.amp.N,aes(y=S.KO, x= Type,fill=Type))+
  geom_boxplot(outlier.shape = NA,size=0.8, width = 0.5,alpha=0.2,aes(fill=Type),position = position_dodge(0.5),) +
  
  geom_jitter(aes(Type, S.KO,color = Type),position=position_jitter(width=0.01,height=0),
              alpha=0.6,stroke=2,shape=1, size = 3)+
  
  theme_bw()+
  xlab(bquote('Nitrogen type ')) +
  ylab(bquote('S.KO'))+
  scale_fill_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_color_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_y_continuous(limits = c(9.5, 11))+
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 18),
        axis.text.y =element_text(color="black",size = 18),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 12),
        plot.title=element_text(face = "bold", colour = "black", size = 16),
        legend.title = element_text(face = "bold", colour = "black", size = 14))+ 
  
  labs(color="Nitrogen type")+
  ggtitle(expression(~Acidobacteriota~~~italic(P)~'= 0.400'))
p1




#·ÅÏß¾ú
env.amp.N <- read.csv("S.KO.Actino.csv",head=T,row.names = 1)


p2 <-ggplot(data=env.amp.N,aes(y=S.KO, x= Type,fill=Type))+
  geom_boxplot(outlier.shape = NA,size=0.8, width = 0.5,alpha=0.2,aes(fill=Type),position = position_dodge(0.5),) +
  
  geom_jitter(aes(Type, S.KO,color = Type),position=position_jitter(width=0.01,height=0),
              alpha=0.6,stroke=2,shape=1, size = 3)+
  
  theme_bw()+
  xlab(bquote('Nitrogen type ')) +
  ylab(bquote('S.KO'))+
  scale_fill_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_color_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_y_continuous(limits = c(9.5, 11))+
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 18),
        axis.text.y =element_text(color="black",size = 18),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 12),
        plot.title=element_text(face = "bold", colour = "black", size = 16),
        legend.title = element_text(face = "bold", colour = "black", size = 14))+ 
  
  labs(color="Nitrogen type")+
  ggtitle(expression(~Actinobacteriota~~~italic(P)~'= 0.315'))
p2



#Äâ¸Ë¾ú
env.amp.N <- read.csv("S.KO.Bact.csv",head=T,row.names = 1)


p3 <-ggplot(data=env.amp.N,aes(y=S.KO, x= Type,fill=Type))+
  geom_boxplot(outlier.shape = NA,size=0.8, width = 0.5,alpha=0.2,aes(fill=Type),position = position_dodge(0.5),) +
  
  geom_jitter(aes(Type, S.KO,color = Type),position=position_jitter(width=0.01,height=0),
              alpha=0.6,stroke=2,shape=1, size = 3)+
  
  theme_bw()+
  xlab(bquote('Nitrogen type ')) +
  ylab(bquote('S.KO'))+
  scale_fill_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_color_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_y_continuous(limits = c(9.5, 11))+
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 18),
        axis.text.y =element_text(color="black",size = 18),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 12),
        plot.title=element_text(face = "bold", colour = "black", size = 16),
        legend.title = element_text(face = "bold", colour = "black", size = 14))+ 
  
  labs(color="Nitrogen type")+
  ggtitle(expression(~Bacteroidota~~~italic(P)~'= 0.549'))
p3



#GEMM
env.amp.N <- read.csv("S.KO.GEMM.csv",head=T,row.names = 1)


p4 <-ggplot(data=env.amp.N,aes(y=S.KO, x= Type,fill=Type))+
  geom_boxplot(outlier.shape = NA,size=0.8, width = 0.5,alpha=0.2,aes(fill=Type),position = position_dodge(0.5),) +
  
  geom_jitter(aes(Type, S.KO,color = Type),position=position_jitter(width=0.01,height=0),
              alpha=0.6,stroke=2,shape=1, size = 3)+
  
  theme_bw()+
  xlab(bquote('Nitrogen type ')) +
  ylab(bquote('S.KO'))+
  scale_fill_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_color_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_y_continuous(limits = c(9.5, 11))+
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 18),
        axis.text.y =element_text(color="black",size = 18),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 12),
        plot.title=element_text(face = "bold", colour = "black", size = 16),
        legend.title = element_text(face = "bold", colour = "black", size = 14))+ 
  
  labs(color="Nitrogen type")+
  ggtitle(expression(~Gemmatimonadota~~~italic(P)~'= 0.414'))
p4


#Plancto
env.amp.N <- read.csv("S.KO.Plancto.csv",head=T,row.names = 1)


p5 <-ggplot(data=env.amp.N,aes(y=S.KO, x= Type,fill=Type))+
  geom_boxplot(outlier.shape = NA,size=0.8, width = 0.5,alpha=0.2,aes(fill=Type),position = position_dodge(0.5),) +
  
  geom_jitter(aes(Type, S.KO,color = Type),position=position_jitter(width=0.01,height=0),
              alpha=0.6,stroke=2,shape=1, size = 3)+
  
  theme_bw()+
  xlab(bquote('Nitrogen type ')) +
  ylab(bquote('S.KO'))+
  scale_fill_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_color_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_y_continuous(limits = c(9.5, 11))+
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 18),
        axis.text.y =element_text(color="black",size = 18),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 12),
        plot.title=element_text(face = "bold", colour = "black", size = 16),
        legend.title = element_text(face = "bold", colour = "black", size = 14))+ 
  
  labs(color="Nitrogen type")+
  ggtitle(expression(~Planctomycetota~~~italic(P)~'= 1.000'))
p5




#±äÐÎ¾ú
env.amp.N <- read.csv("S.KO.Pro.csv",head=T,row.names = 1)


p6 <-ggplot(data=env.amp.N,aes(y=S.KO, x= Type,fill=Type))+
  geom_boxplot(outlier.shape = NA,size=0.8, width = 0.5,alpha=0.2,aes(fill=Type),position = position_dodge(0.5),) +
  
  geom_jitter(aes(Type, S.KO,color = Type),position=position_jitter(width=0.01,height=0),
              alpha=0.6,stroke=2,shape=1, size = 3)+
  
  theme_bw()+
  xlab(bquote('Nitrogen type ')) +
  ylab(bquote('S.KO'))+
  scale_fill_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_color_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_y_continuous(limits = c(9.5, 11))+
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 18),
        axis.text.y =element_text(color="black",size = 18),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 12),
        plot.title=element_text(face = "bold", colour = "black", size = 16),
        legend.title = element_text(face = "bold", colour = "black", size = 14))+ 
  
  labs(color="Nitrogen type")+
  ggtitle(expression(~Proteobacteria~~~italic(P)~'= 0.372'))
p6




p14 <- 
  p1 + p2 + p3 + p4 + p5 + p6+  
  plot_layout(guides = "collect",ncol=6,nrow=1) +
  theme(legend.position='right')


p14




####Éú³¤ËÙÂÊ
#ËáÏ¸¾ú
rm(list=ls())
env.amp.N <- read.csv("speed.Aci.csv",head=T,row.names = 1)


p1 <-ggplot(data=env.amp.N,aes(y=d, x= Type,fill=Type))+
  geom_boxplot(outlier.shape = NA,size=0.8, width = 0.5,alpha=0.2,aes(fill=Type),position = position_dodge(0.5),) +
  
  geom_jitter(aes(Type, d,color = Type),position=position_jitter(width=0.01,height=0),
              alpha=0.6,stroke=2,shape=1, size = 3)+
  
  theme_bw()+
  xlab(bquote('Nitrogen type ')) +
  ylab(bquote('Minimal doubling time (hours)'))+
  scale_fill_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_color_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_y_continuous(limits = c(0, 16))+
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 18),
        axis.text.y =element_text(color="black",size = 18),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 12),
        plot.title=element_text(face = "bold", colour = "black", size = 16),
        legend.title = element_text(face = "bold", colour = "black", size = 14))+ 
  
  labs(color="Nitrogen type")+
  ggtitle(expression(~Acidobacteriota~~~italic(P)~'= 0.400'))
p1



#·ÅÏß¾ú
env.amp.N <- read.csv("speed.Actino.csv",head=T,row.names = 1)


p2 <-ggplot(data=env.amp.N,aes(y=d,x= Type,fill=Type))+
  geom_boxplot(outlier.shape = NA,size=0.8, width = 0.5,alpha=0.2,aes(fill=Type),position = position_dodge(0.5),) +
  
  geom_jitter(aes(Type, d,color = Type),position=position_jitter(width=0.01,height=0),
              alpha=0.6,stroke=2,shape=1, size = 3)+
  
  theme_bw()+
  xlab(bquote('Nitrogen type ')) +
  ylab(bquote('Minimal doubling time (hours)'))+
  scale_fill_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_color_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_y_continuous(limits = c(0, 16))+
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 18),
        axis.text.y =element_text(color="black",size = 18),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 12),
        plot.title=element_text(face = "bold", colour = "black", size = 16),
        legend.title = element_text(face = "bold", colour = "black", size = 14))+ 
  
  labs(color="Nitrogen type")+
  ggtitle(expression(~Actinobacteriota~~~italic(P)~'= 0.062'))
p2


#Äâ¸Ë¾ú
env.amp.N <- read.csv("speed.Bact.csv",head=T,row.names = 1)


p3 <-ggplot(data=env.amp.N,aes(y=d, x= Type,fill=Type))+
  geom_boxplot(outlier.shape = NA,size=0.8, width = 0.5,alpha=0.2,aes(fill=Type),position = position_dodge(0.5),) +
  
  geom_jitter(aes(Type, d,color = Type),position=position_jitter(width=0.01,height=0),
              alpha=0.6,stroke=2,shape=1, size = 3)+
  
  theme_bw()+
  xlab(bquote('Nitrogen type ')) +
  ylab(bquote('Minimal doubling time (hours)'))+
  scale_fill_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_color_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_y_continuous(limits = c(0, 16))+
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 18),
        axis.text.y =element_text(color="black",size = 18),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 12),
        plot.title=element_text(face = "bold", colour = "black", size = 16),
        legend.title = element_text(face = "bold", colour = "black", size = 14))+ 
  
  labs(color="Nitrogen type")+
  ggtitle(expression(~Bacteroidota~~~italic(P)~'= 0.937'))
p3



#GEMM
env.amp.N <- read.csv("speed.GEMM.csv",head=T,row.names = 1)


p4 <-ggplot(data=env.amp.N,aes(y=d, x= Type,fill=Type))+
  geom_boxplot(outlier.shape = NA,size=0.8, width = 0.5,alpha=0.2,aes(fill=Type),position = position_dodge(0.5),) +
  
  geom_jitter(aes(Type, d,color = Type),position=position_jitter(width=0.01,height=0),
              alpha=0.6,stroke=2,shape=1, size = 3)+
  
  theme_bw()+
  xlab(bquote('Nitrogen type ')) +
  ylab(bquote('Minimal doubling time (hours)'))+
  scale_fill_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_color_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_y_continuous(limits = c(0, 16))+
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 18),
        axis.text.y =element_text(color="black",size = 18),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 12),
        plot.title=element_text(face = "bold", colour = "black", size = 16),
        legend.title = element_text(face = "bold", colour = "black", size = 14))+ 
  
  labs(color="Nitrogen type")+
  ggtitle(expression(~Gemmatimonadota~~~italic(P)~'= 0.165'))
p4


#Plancto
env.amp.N <- read.csv("speed.Plancto.csv",head=T,row.names = 1)


p5 <-ggplot(data=env.amp.N,aes(y=d, x= Type,fill=Type))+
  geom_boxplot(outlier.shape = NA,size=0.8, width = 0.5,alpha=0.2,aes(fill=Type),position = position_dodge(0.5),) +
  
  geom_jitter(aes(Type, d,color = Type),position=position_jitter(width=0.01,height=0),
              alpha=0.6,stroke=2,shape=1, size = 3)+
  
  theme_bw()+
  xlab(bquote('Nitrogen type ')) +
  ylab(bquote('Minimal doubling time (hours)'))+
  scale_fill_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_color_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_y_continuous(limits = c(0, 16))+
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 18),
        axis.text.y =element_text(color="black",size = 18),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 12),
        plot.title=element_text(face = "bold", colour = "black", size = 16),
        legend.title = element_text(face = "bold", colour = "black", size = 14))+ 
  
  labs(color="Nitrogen type")+
  ggtitle(expression(~Planctomycetota~~~italic(P)~'= 1.000'))
p5



#±äÐÎ¾ú
env.amp.N <- read.csv("speed.Pro.csv",head=T,row.names = 1)


p6 <-ggplot(data=env.amp.N,aes(y=d, x= Type,fill=Type))+
  geom_boxplot(outlier.shape = NA,size=0.8, width = 0.5,alpha=0.2,aes(fill=Type),position = position_dodge(0.5),) +
  
  geom_jitter(aes(Type, d,color = Type),position=position_jitter(width=0.01,height=0),
              alpha=0.6,stroke=2,shape=1, size = 3)+
  
  theme_bw()+
  xlab(bquote('Nitrogen type ')) +
  ylab(bquote('Minimal doubling time (hours)'))+
  scale_fill_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_color_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_y_continuous(limits = c(0, 16))+
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x =element_text(color="black",size = 18),
        axis.text.y =element_text(color="black",size = 18),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 12),
        plot.title=element_text(face = "bold", colour = "black", size = 16),
        legend.title = element_text(face = "bold", colour = "black", size = 14))+ 
  
  labs(color="Nitrogen type")+
  ggtitle(expression(~Proteobacteria~~~italic(P)~'= 0.514'))
p6


p15 <- 
  p1 + p2 + p3 + p4 + p5 + p6+  
  plot_layout(guides = "collect",ncol=6,nrow=1) +
  theme(legend.position='right')

p15




####Fig.S32
library(EnhancedVolcano)
library(patchwork)
library(edgeR)
library(statmod)

#Acidobacteriota
gene5<-read.csv("data/2DP.ACI.YAS.AN.csv",row.names=1,head=T)
gene5[which(gene5$FDR < 0.05 & gene5$Acidobacteriota <= -1),'sig'] <- 'None'
gene5[which(gene5$FDR < 0.05 & gene5$Acidobacteriota >=1),'sig'] <- 'Up'
gene5[which(gene5$FDR >= 0.05 | abs(gene5$Acidobacteriota) < 1),'sig'] <- 'None'

sum(gene5$sig=='Down')
sum(gene5$sig=="Up")
sum(gene5$sig=="None")

p1<-ggplot(data=gene5, aes(x= Acidobacteriota , y= FDR,group=sig))+ 
  geom_point(aes(colour=sig),size=4,alpha=0.5)+ 
  labs(x = "2DP (AN)", y = "P value (FDR)")+
  scale_colour_manual(values=c("lightskyblue","blue","gray70","#5EB5FF","#339CFF","blue",
                               "#FFE1FF","#FFC8FF","#FFAFFF","#FF63FF","#FF6371"),
                      breaks=c("Up","Down", "None"),
                      labels=c("Acidobacteriota preference", "AN preference","None"))+
  theme_bw()+
  ggtitle(expression('Acidobacteriota'),subtitle = "Preferenced KOs = 649 variables")+
  theme(legend.position = 'left',legend.direction = 'vertical',
        legend.title = element_text(colour="black", size=13, face="bold"),
        legend.text = element_text(colour="black", size=12),
        axis.text=element_text(size=15,face="bold",colour="black"),
        axis.title=element_text(size=17,face="bold"))+
  geom_vline(xintercept = 1, color = 'gray34', linetype = 2, linewidth = 1.3) +
  geom_hline(yintercept = 0.05, color = 'gray34', linetype = 2, linewidth = 1.3)

p1


#Actinobacteriota
gene6<-read.csv("data/2DP.ACTINO.YAS.AN.csv",row.names=1,head=T)
gene6[which(gene6$FDR < 0.05 & gene6$Actinobacteriota <= -1),'sig'] <- 'None'
gene6[which(gene6$FDR < 0.05 & gene6$Actinobacteriota >=1),'sig'] <- 'Up'
gene6[which(gene6$FDR >= 0.05 | abs(gene6$Actinobacteriota) < 1),'sig'] <- 'None'

sum(gene6$sig=='Down')
sum(gene6$sig=="Up")
sum(gene6$sig=="None")

p2<-ggplot(data=gene6, aes(x= Actinobacteriota , y= FDR,group=sig))+ 
  geom_point(aes(colour=sig),size=4,alpha=0.5)+ 
  labs(x = "2DP (AN)", y = "P value (FDR)")+
  scale_colour_manual(values=c("lightgoldenrod1","#FFC8FF","gray70","#5EB5FF","#339CFF","blue",
                               "#FFC8FF","#FFAFFF","#FF63FF","#FF6371","red"),
                      breaks=c("Up","Down", "None"),
                      labels=c("Actinobacteriota preference", "AU00 preference","None"))+
  theme_bw()+
  ggtitle(expression('Actinobacteriota'),subtitle = "Preferenced KOs = 1654 variables")+
  theme(legend.position = 'left',legend.direction = 'vertical',
        legend.title = element_text(colour="black", size=13, face="bold"),
        legend.text = element_text(colour="black", size=12),
        axis.text=element_text(size=15,face="bold",colour="black"),
        axis.title=element_text(size=17,face="bold"))+
  geom_vline(xintercept = 1, color = 'gray34', linetype = 2, linewidth = 1.3) +
  geom_hline(yintercept = 0.05, color = 'gray34', linetype = 2, linewidth = 1.3) 

p2


#Bacteroidota
gene7<-read.csv("data/2DP.BAC.YAS.AN.csv",row.names=1,head=T)
gene7[which(gene7$FDR < 0.05 & gene7$Bacteroidota <= -1),'sig'] <- 'None'
gene7[which(gene7$FDR < 0.05 & gene7$Bacteroidota >=1),'sig'] <- 'Up'
gene7[which(gene7$FDR >= 0.05 | abs(gene7$Bacteroidota) < 1),'sig'] <- 'None'

sum(gene7$sig=='Down')
sum(gene7$sig=="Up")
sum(gene7$sig=="None")

p3<-ggplot(data=gene7, aes(x= Bacteroidota , y= FDR,group=sig))+ 
  geom_point(aes(colour=sig),size=4,alpha=0.5)+ 
  labs(x = "2DP (AN)", y = "P value (FDR)")+
  scale_colour_manual(values=c("violet","#FFC8FF","gray70","#5EB5FF","#339CFF","blue",
                               "#FFC8FF","#FFAFFF","#FF63FF","#FF6371","red"),
                      breaks=c("Up","Down", "None"),
                      labels=c("Bacteroidota preference", "AU00 preference","None"))+
  theme_bw()+
  ggtitle(expression('Bacteroidota'),subtitle = "Preferenced KOs = 1101 variables")+
  theme(legend.position = 'left',legend.direction = 'vertical',
        legend.title = element_text(colour="black", size=13, face="bold"),
        legend.text = element_text(colour="black", size=12),
        axis.text=element_text(size=15,face="bold",colour="black"),
        axis.title=element_text(size=17,face="bold"))+
  geom_vline(xintercept = 1, color = 'gray34', linetype = 2, linewidth = 1.3) +
  geom_hline(yintercept = 0.05, color = 'gray34', linetype = 2, linewidth = 1.3)

p3


#Gemmatimonadota
gene8<-read.csv("data/2DP.GEMM.YAS.AN.csv",row.names=1,head=T)
gene8[which(gene8$FDR < 0.05 & gene8$Gemmatimonadota <= -1),'sig'] <- 'None'
gene8[which(gene8$FDR < 0.05 & gene8$Gemmatimonadota >=1),'sig'] <- 'Up'
gene8[which(gene8$FDR >= 0.05 | abs(gene8$Gemmatimonadota) < 1),'sig'] <- 'None'

sum(gene8$sig=='Down')
sum(gene8$sig=="Up")
sum(gene8$sig=="None")

p4<-ggplot(data=gene8, aes(x= Gemmatimonadota , y= FDR,group=sig))+ 
  geom_point(aes(colour=sig),size=4,alpha=0.5)+ 
  labs(x = "2DP (AN)", y = "P value (FDR)")+
  scale_colour_manual(values=c("seagreen1","#FFC8FF","gray70","#5EB5FF","#339CFF","blue",
                               "#FFC8FF","#FFAFFF","#FF63FF","#FF6371","red"),
                      breaks=c("Up","Down", "None"),
                      labels=c("Gemmatimonadota preference", "AU00 preference","None"))+
  theme_bw()+
  ggtitle(expression('Gemmatimonadota'),subtitle = "Preferenced KOs = 983 variables")+
  theme(legend.position = 'left',legend.direction = 'vertical',
        legend.title = element_text(colour="black", size=13, face="bold"),
        legend.text = element_text(colour="black", size=12),
        axis.text=element_text(size=15,face="bold",colour="black"),
        axis.title=element_text(size=17,face="bold"))+
  geom_vline(xintercept = 1, color = 'gray34', linetype = 2, linewidth = 1.3) +
  geom_hline(yintercept = 0.05, color = 'gray34', linetype = 2, linewidth = 1.3) 

p4


#Planctomycetota
gene9<-read.csv("data/2DP.PLAN.YAS.AN.csv",row.names=1,head=T)
gene9[which(gene9$FDR < 0.05 & gene9$Planctomycetota <= -1),'sig'] <- 'None'
gene9[which(gene9$FDR < 0.05 & gene9$Planctomycetota >=1),'sig'] <- 'Up'
gene9[which(gene9$FDR >= 0.05 | abs(gene9$Planctomycetota) < 1),'sig'] <- 'None'

sum(gene9$sig=='Down')
sum(gene9$sig=="Up")
sum(gene9$sig=="None")

p5<-ggplot(data=gene9, aes(x= Planctomycetota , y= FDR,group=sig))+ 
  geom_point(aes(colour=sig),size=4,alpha=0.5)+ 
  labs(x = "2DP (AN)", y = "P value (FDR)")+
  scale_colour_manual(values=c("lightpink2","#FFC8FF","gray70","#5EB5FF","#339CFF","blue",
                               "#FFC8FF","#FFAFFF","#FF63FF","#FF6371","red"),
                      breaks=c("Up","Down", "None"),
                      labels=c("Planctomycetota preference", "AU00 preference","None"))+
  theme_bw()+
  ggtitle(expression('Planctomycetota'),subtitle = "Preferenced KOs = 832 variables")+
  theme(legend.position = 'left',legend.direction = 'vertical',
        legend.title = element_text(colour="black", size=13, face="bold"),
        legend.text = element_text(colour="black", size=12),
        axis.text=element_text(size=15,face="bold",colour="black"),
        axis.title=element_text(size=17,face="bold"))+
  geom_vline(xintercept = 1, color = 'gray34', linetype = 2, linewidth = 1.3) +
  geom_hline(yintercept = 0.05, color = 'gray34', linetype = 2, linewidth = 1.3) 

p5


#Proteobacteria
gene10<-read.csv("data/2DP.PRO.YAS.AN.csv",row.names=1,head=T)
gene10[which(gene10$FDR < 0.05 & gene10$Proteobacteria <= -1),'sig'] <- 'None'
gene10[which(gene10$FDR < 0.05 & gene10$Proteobacteria >=1),'sig'] <- 'Up'
gene10[which(gene10$FDR >= 0.05 | abs(gene10$Proteobacteria) < 1),'sig'] <- 'None'

sum(gene10$sig=='Down')
sum(gene10$sig=="Up")
sum(gene10$sig=="None")

p6<-ggplot(data=gene10, aes(x= Proteobacteria , y= FDR,group=sig))+ 
  geom_point(aes(colour=sig),size=4,alpha=0.5)+ 
  labs(x = "2DP (AN)", y = "P value (FDR)")+
  scale_colour_manual(values=c("paleturquoise1","#FFC8FF","gray70","#5EB5FF","#339CFF","blue",
                               "#FFC8FF","#FFAFFF","#FF63FF","#FF6371","red"),
                      breaks=c("Up","Down", "None"),
                      labels=c("Proteobacteria preference", "AU00 preference","None"))+
  theme_bw()+
  ggtitle(expression('Proteobacteria'),subtitle = "Preferenced KOs = 1810 variables")+
  theme(legend.position = 'left',legend.direction = 'vertical',
        legend.title = element_text(colour="black", size=13, face="bold"),
        legend.text = element_text(colour="black", size=12),
        axis.text=element_text(size=15,face="bold",colour="black"),
        axis.title=element_text(size=17,face="bold"))+
  geom_vline(xintercept = 1, color = 'gray34', linetype = 2, linewidth = 1.3) +
  geom_hline(yintercept = 0.05, color = 'gray34', linetype = 2, linewidth = 1.3)

p6


p <- 
  p1+ p2 +p3+p4+p5+p6+
  plot_layout(guides = "collect",ncol=3,nrow=2) +
  theme(legend.position='right')

p
ggsave("Fig.S32.pdf", p,height=9,width =15,limitsize = FALSE )


####Fig.S33
#Acidobacteriota
gene5<-read.csv("data/2DP.ACI.YAS.AU.csv",row.names=1,head=T)
gene5[which(gene5$FDR < 0.05 & gene5$Acidobacteriota <= -1),'sig'] <- 'None'
gene5[which(gene5$FDR < 0.05 & gene5$Acidobacteriota >=1),'sig'] <- 'Up'
gene5[which(gene5$FDR >= 0.05 | abs(gene5$Acidobacteriota) < 1),'sig'] <- 'None'

sum(gene5$sig=='Down')
sum(gene5$sig=="Up")
sum(gene5$sig=="None")

p1<-ggplot(data=gene5, aes(x= Acidobacteriota , y= FDR,group=sig))+ 
  geom_point(aes(colour=sig),size=4,alpha=0.5)+ 
  labs(x = "2DP (AU)", y = "P value (FDR)")+
  scale_colour_manual(values=c("lightskyblue","blue","gray70","#5EB5FF","#339CFF","blue",
                               "#FFE1FF","#FFC8FF","#FFAFFF","#FF63FF","#FF6371"),
                      breaks=c("Up","Down", "None"),
                      labels=c("Acidobacteriota preference", "AN preference","None"))+
  theme_bw()+
  ggtitle(expression('Acidobacteriota'),subtitle = "Preferenced KOs = 1194 variables")+
  theme(legend.position = 'left',legend.direction = 'vertical',
        legend.title = element_text(colour="black", size=13, face="bold"),
        legend.text = element_text(colour="black", size=12),
        axis.text=element_text(size=15,face="bold",colour="black"),
        axis.title=element_text(size=17,face="bold"))+
  geom_vline(xintercept = 1, color = 'gray34', linetype = 2, linewidth = 1.3) +
  geom_hline(yintercept = 0.05, color = 'gray34', linetype = 2, linewidth = 1.3) 

p1


#Actinobacteriota
gene6<-read.csv("data/2DP.ACTINO.YAS.AU.csv",row.names=1,head=T)
gene6[which(gene6$FDR < 0.05 & gene6$Actinobacteriota <= -1),'sig'] <- 'None'
gene6[which(gene6$FDR < 0.05 & gene6$Actinobacteriota >=1),'sig'] <- 'Up'
gene6[which(gene6$FDR >= 0.05 | abs(gene6$Actinobacteriota) < 1),'sig'] <- 'None'

sum(gene6$sig=='Down')
sum(gene6$sig=="Up")
sum(gene6$sig=="None")

p2<-ggplot(data=gene6, aes(x= Actinobacteriota , y= FDR,group=sig))+ 
  geom_point(aes(colour=sig),size=4,alpha=0.5)+ 
  labs(x = "2DP (AU)", y = "P value (FDR)")+
  scale_colour_manual(values=c("lightgoldenrod1","#FFC8FF","gray70","#5EB5FF","#339CFF","blue",
                               "#FFC8FF","#FFAFFF","#FF63FF","#FF6371","red"),
                      breaks=c("Up","Down", "None"),
                      labels=c("Actinobacteriota preference", "AU00 preference","None"))+
  theme_bw()+
  ggtitle(expression('Actinobacteriota'),subtitle = "Preferenced KOs = 2058 variables")+
  theme(legend.position = 'left',legend.direction = 'vertical',
        legend.title = element_text(colour="black", size=13, face="bold"),
        legend.text = element_text(colour="black", size=12),
        axis.text=element_text(size=15,face="bold",colour="black"),
        axis.title=element_text(size=17,face="bold"))+
  geom_vline(xintercept = 1, color = 'gray34', linetype = 2, linewidth = 1.3) +
  geom_hline(yintercept = 0.05, color = 'gray34', linetype = 2, linewidth = 1.3) 

p2


#Bacteroidota
gene7<-read.csv("data/2DP.BAC.YAS.AU.csv",row.names=1,head=T)
gene7[which(gene7$FDR < 0.05 & gene7$Bacteroidota <= -1),'sig'] <- 'None'
gene7[which(gene7$FDR < 0.05 & gene7$Bacteroidota >=1),'sig'] <- 'Up'
gene7[which(gene7$FDR >= 0.05 | abs(gene7$Bacteroidota) < 1),'sig'] <- 'None'

sum(gene7$sig=='Down')
sum(gene7$sig=="Up")
sum(gene7$sig=="None")

p3<-ggplot(data=gene7, aes(x= Bacteroidota , y= FDR,group=sig))+ 
  geom_point(aes(colour=sig),size=4,alpha=0.5)+ 
  labs(x = "2DP (AU)", y = "P value (FDR)")+
  scale_colour_manual(values=c("violet","#FFC8FF","gray70","#5EB5FF","#339CFF","blue",
                               "#FFC8FF","#FFAFFF","#FF63FF","#FF6371","red"),
                      breaks=c("Up","Down", "None"),
                      labels=c("Bacteroidota preference", "AU00 preference","None"))+
  theme_bw()+
  ggtitle(expression('Bacteroidota'),subtitle = "Preferenced KOs = 1375 variables")+
  theme(legend.position = 'left',legend.direction = 'vertical',
        legend.title = element_text(colour="black", size=13, face="bold"),
        legend.text = element_text(colour="black", size=12),
        axis.text=element_text(size=15,face="bold",colour="black"),
        axis.title=element_text(size=17,face="bold"))+
  geom_vline(xintercept = 1, color = 'gray34', linetype = 2, linewidth = 1.3) +
  geom_hline(yintercept = 0.05, color = 'gray34', linetype = 2, linewidth = 1.3) 

p3


#Gemmatimonadota
gene8<-read.csv("data/2DP.GEMM.YAS.AU.csv",row.names=1,head=T)
gene8[which(gene8$FDR < 0.05 & gene8$Gemmatimonadota <= -1),'sig'] <- 'None'
gene8[which(gene8$FDR < 0.05 & gene8$Gemmatimonadota >=1),'sig'] <- 'Up'
gene8[which(gene8$FDR >= 0.05 | abs(gene8$Gemmatimonadota) < 1),'sig'] <- 'None'

sum(gene8$sig=='Down')
sum(gene8$sig=="Up")
sum(gene8$sig=="None")

p4<-ggplot(data=gene8, aes(x= Gemmatimonadota , y= FDR,group=sig))+ 
  geom_point(aes(colour=sig),size=4,alpha=0.5)+ 
  labs(x = "2DP (AU)", y = "P value (FDR)")+
  scale_colour_manual(values=c("seagreen1","#FFC8FF","gray70","#5EB5FF","#339CFF","blue",
                               "#FFC8FF","#FFAFFF","#FF63FF","#FF6371","red"),
                      breaks=c("Up","Down", "None"),
                      labels=c("Gemmatimonadota preference", "AU00 preference","None"))+
  theme_bw()+
  ggtitle(expression('Gemmatimonadota'),subtitle = "Preferenced KOs = 997 variables")+
  theme(legend.position = 'left',legend.direction = 'vertical',
        legend.title = element_text(colour="black", size=13, face="bold"),
        legend.text = element_text(colour="black", size=12),
        axis.text=element_text(size=15,face="bold",colour="black"),
        axis.title=element_text(size=17,face="bold"))+
  geom_vline(xintercept = 1, color = 'gray34', linetype = 2, linewidth = 1.3) +
  geom_hline(yintercept = 0.05, color = 'gray34', linetype = 2, linewidth = 1.3) 

p4


#Planctomycetota
gene9<-read.csv("data/2DP.PLAN.YAS.AU.csv",row.names=1,head=T)
gene9[which(gene9$FDR < 0.05 & gene9$Planctomycetota <= -1),'sig'] <- 'None'
gene9[which(gene9$FDR < 0.05 & gene9$Planctomycetota >=1),'sig'] <- 'Up'
gene9[which(gene9$FDR >= 0.05 | abs(gene9$Planctomycetota) < 1),'sig'] <- 'None'

sum(gene9$sig=='Down')
sum(gene9$sig=="Up")
sum(gene9$sig=="None")

p5<-ggplot(data=gene9, aes(x= Planctomycetota , y= FDR,group=sig))+ 
  geom_point(aes(colour=sig),size=4,alpha=0.5)+ 
  labs(x = "2DP (AU)", y = "P value (FDR)")+
  scale_colour_manual(values=c("lightpink2","#FFC8FF","gray70","#5EB5FF","#339CFF","blue",
                               "#FFC8FF","#FFAFFF","#FF63FF","#FF6371","red"),
                      breaks=c("Up","Down", "None"),
                      labels=c("Planctomycetota preference", "AU00 preference","None"))+
  theme_bw()+
  ggtitle(expression('Planctomycetota'),subtitle = "Preferenced KOs = 839 variables")+
  theme(legend.position = 'left',legend.direction = 'vertical',
        legend.title = element_text(colour="black", size=13, face="bold"),
        legend.text = element_text(colour="black", size=12),
        axis.text=element_text(size=15,face="bold",colour="black"),
        axis.title=element_text(size=17,face="bold"))+
  geom_vline(xintercept = 1, color = 'gray34', linetype = 2, linewidth = 1.3) +
  geom_hline(yintercept = 0.05, color = 'gray34', linetype = 2, linewidth = 1.3)

p5


#Proteobacteria
gene10<-read.csv("data/2DP.PRO.YAS.AU.csv",row.names=1,head=T)
gene10[which(gene10$FDR < 0.05 & gene10$Proteobacteria <= -1),'sig'] <- 'None'
gene10[which(gene10$FDR < 0.05 & gene10$Proteobacteria >=1),'sig'] <- 'Up'
gene10[which(gene10$FDR >= 0.05 | abs(gene10$Proteobacteria) < 1),'sig'] <- 'None'

sum(gene10$sig=='Down')
sum(gene10$sig=="Up")
sum(gene10$sig=="None")

p6<-ggplot(data=gene10, aes(x= Proteobacteria , y= FDR,group=sig))+ 
  geom_point(aes(colour=sig),size=4,alpha=0.5)+ 
  labs(x = "2DP (AU)", y = "P value (FDR)")+
  scale_colour_manual(values=c("paleturquoise1","#FFC8FF","gray70","#5EB5FF","#339CFF","blue",
                               "#FFC8FF","#FFAFFF","#FF63FF","#FF6371","red"),
                      breaks=c("Up","Down", "None"),
                      labels=c("Proteobacteria preference", "AU00 preference","None"))+
  theme_bw()+
  ggtitle(expression('Proteobacteria'),subtitle = "Preferenced KOs = 1827 variables")+
  theme(legend.position = 'left',legend.direction = 'vertical',
        legend.title = element_text(colour="black", size=13, face="bold"),
        legend.text = element_text(colour="black", size=12),
        axis.text=element_text(size=15,face="bold",colour="black"),
        axis.title=element_text(size=17,face="bold"))+
  geom_vline(xintercept = 1, color = 'gray34', linetype = 2, linewidth = 1.3) +
  geom_hline(yintercept = 0.05, color = 'gray34', linetype = 2, linewidth = 1.3)

p6

p2 <- 
  p1+ p2 +p3+p4+p5+p6+
  plot_layout(guides = "collect",ncol=3,nrow=2) +
  theme(legend.position='right')

p2
ggsave("Fig.S33.pdf", p2,height=9,width =15,limitsize = FALSE )

