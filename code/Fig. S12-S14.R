#Fig.S12-S14
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

#Fig.S12a-j
rm(list=ls())

load("data/metanaper.data.prep.2022.10.25.RData")
env.amp.N <- env.amp1[1:96,]
env.amp.AN <- env.amp1[1:48,]
env.amp.AU <- env.amp1[49:96,]
env.amp.N <- env.amp.N[!is.na(env.amp.N$pH),]
env.amp.AN <- env.amp.AN[!is.na(env.amp.AN$pH),]
env.amp.AU <- env.amp.AU[!is.na(env.amp.AU$pH),]

#S.Bacteria
p1 <-ggplot()+
  geom_jitter(data=env.amp.N, aes(y=bac_shannon, x= log(level3+1), color=type1),size=10,alpha=0.8)+
  geom_smooth( data=env.amp.N,aes(y=bac_shannon, x= log(level3+1),color=type1),method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('Ln(N addition rate ('*g~N~m ^-2~yr^-1*'))')) +
  ylab(bquote('S.Bacteria'))+
  scale_color_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        plot.margin = unit(c(0,1,0,1),"cm"),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+ 
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R)["AN"] ^2~'= 0.462***'~~~~~italic(R)["AU"] ^2~'= 0.489***'))

p1
summary(lm(env.amp.AN$bac_shannon~log(env.amp.AN$level3+1)))
summary(lm(env.amp.AU$bac_shannon~log(env.amp.AU$level3+1)))

library(simba)
#To compare the regression coefficient. p<0.05 means significant. Details: ?diffslope
set.seed(123)
diffslope(env.amp.AN$bac_shannon, (env.amp.AN$level3), env.amp.AU$bac_shannon, (env.amp.AU$level3))


#S.fun
p2 <-ggplot()+
  geom_jitter(data=env.amp.N, aes(y=fun_shannon, x= log(level3+1), color=type1),size=10,alpha=0.8)+
  geom_smooth(data=subset(env.amp.N,env.amp.N$type1=="AN"),aes(x= log(level3+1), y=fun_shannon,color=type1), method ="lm", span=0.9,linewidth=2)+
  geom_smooth(data=subset(env.amp.N,env.amp.N$type1=="AU"),aes(x= log(level3+1), y=fun_shannon,color=type1), method ="lm", span=0.9,linewidth=2)+
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('Ln(N addition rate ('*g~N~m ^-2~yr^-1*'))')) +
  ylab(bquote('S.Fungi'))+
  scale_color_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        plot.margin = unit(c(0,1,0,1),"cm"),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+ 
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R)["AN"] ^2~'= 0.652***'~~~~~italic(R)["AU"] ^2~'= 0.204***'))

p2
summary(lm(env.amp.AN$fun_shannon~log(env.amp.AN$level3+1)))
summary(lm(env.amp.AU$fun_shannon~log(env.amp.AU$level3+1)))

library(simba)
#To compare the regression coefficient. p<0.05 means significant. Details: ?diffslope
set.seed(123)
diffslope(env.amp.AN$fun_shannon, (env.amp.AN$level3), env.amp.AU$fun_shannon, (env.amp.AU$level3))


#HGT
data3<-read.csv("data/HGT_number.csv",row.names=1,head=T)
data3.AN <- data3[1:48,]
data3.AU <- data3[49:96,]

p3 <-ggplot()+
  geom_jitter(data=data3, aes(y=HGT_per_1K_assembled_genes, x= log(level3+1), color=type1),size=10,alpha=0.8)+
  geom_smooth(data=subset(data3,data3$type=="AN"),aes(x= log(level3+1), y=HGT_per_1K_assembled_genes,color=type1), method ="lm", span=0.9,linewidth=2)+
  geom_smooth(data=subset(data3,data3$type=="AU"),aes(x= log(level3+1), y=HGT_per_1K_assembled_genes,color=type1), method ="lm", span=0.9,linewidth=2)+
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('Ln(N addition rate ('*g~N~m ^-2~yr^-1*'))')) +
  ylab(bquote('HGT Events \n per 1K Assembled Genes'))+
  scale_color_manual(values = c( "royalblue1","#F24F50", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        plot.margin = unit(c(0,1,0,1),"cm"),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+ 
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R)["AN"] ^2~'= 0.483***'~~~~~italic(R)["AU"] ^2~'= 0.576***'))

p3
summary(lm(data3.AN$HGT_per_1K_assembled_genes~log(data3.AN$level3+1)))
summary(lm(data3.AU$HGT_per_1K_assembled_genes~log(data3.AU$level3+1)))

library(simba)
#To compare the regression coefficient. p<0.05 means significant. Details: ?diffslope
set.seed(123)
diffslope(data3.AN$HGT_per_1K_assembled_genes, (data3.AN$level3+1), data3.AU$HGT_per_1K_assembled_genes, (data3.AU$level3+1))


#S.Phage
p4 <-ggplot()+
  geom_jitter(data=env.amp.N, aes(y=vir_shannon_mgm, x= log(level3+1), color=type1),size=10,alpha=0.8)+
  geom_smooth(data=subset(env.amp.N,env.amp.N$type1=="AN"),aes(x= log(level3+1), y=vir_shannon_mgm,color=type1), method ="lm", span=0.9,linewidth=2)+
  geom_smooth(data=subset(env.amp.N,env.amp.N$type1=="AU"),aes(x= log(level3+1), y=vir_shannon_mgm,color=type1), method ="lm", span=0.9,linewidth=2)+
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('Ln(N addition rate ('*g~N~m ^-2~yr^-1*'))')) +
  ylab(bquote('S.Phage'))+
  scale_color_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        plot.margin = unit(c(0,1,0,1),"cm"),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+ 
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R)["AN"] ^2~'= 0.665***'~~~~~italic(R)["AU"] ^2~'= 0.702***'))

p4
summary(lm(env.amp.AN$vir_shannon_mgm~log(env.amp.AN$level3+1)))
summary(lm(env.amp.AU$vir_shannon_mgm~log(env.amp.AU$level3+1)))

library(simba)
#To compare the regression coefficient. p<0.05 means significant. Details: ?diffslope
set.seed(123)
diffslope(env.amp.AN$vir_shannon_mgm, (env.amp.AN$level3+1), env.amp.AU$vir_shannon_mgm, (env.amp.AU$level3+1))


#Minimal Doubling Time (hours)
data5 <- read.csv("data/d.csv",row.names=1,header=T) 
data5.AN <- data5[1:48,]
data5.AU <- data5[49:96,]

p5 <-ggplot()+
  geom_jitter(data=data5, aes(y=d, x= log(level3+1), color=type1),size=10,alpha=0.8)+
  geom_smooth( data=data5,aes(y=d, x= log(level3+1),color=type1),method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('Ln(N addition rate ('*g~N~m ^-2~yr^-1*'))')) +
  ylab(bquote('Minimal Doubling Time (hours)'))+
  scale_color_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        plot.margin = unit(c(0,1,0,1),"cm"),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+ 
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R)["AN"] ^2~'= 0.521***'~~~~~italic(R)["AU"] ^2~'= 0.440***'))

p5
summary(lm(data5.AN$d~log(data5.AN$level3+1)))
summary(lm(data5.AU$d~log(data5.AU$level3+1)))

library(simba)
#To compare the regression coefficient. p<0.05 means significant. Details: ?diffslope
set.seed(123)
diffslope(data5.AN$d, (data5.AN$level3+1), data5.AU$d, (data5.AU$level3+1))


###Average 16S Copy Number
data6<-read.csv("data/average 16S rRNA copy number.csv",row.names=1,head=T)
data6.AN <- data6[1:48,]
data6.AU <- data6[49:96,]

p6 <-ggplot()+
  geom_jitter(data=data6, aes(y=rrn, x= log(level3+1), color=type1),size=10,alpha=0.8)+
  geom_smooth(data=subset(data6,data6$type=="AN"),aes(x= log(level3+1), y=rrn,color=type1), method ="lm", span=0.9,linetype=2,se=F,linewidth=2)+
  geom_smooth(data=subset(data6,data6$type=="AU"),aes(x= log(level3+1), y=rrn,color=type1), method ="lm", span=0.9,linetype=1,linewidth=2)+
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('Ln(N addition rate ('*g~N~m ^-2~yr^-1*'))')) +
  
  ylab(bquote('Average 16S Copy Number'))+
  scale_color_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        plot.margin = unit(c(0,1,0,1),"cm"),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+ 
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R)["AN"] ^2~'= 0.049'^ns~~~~~italic(R)["AU"] ^2~'= 0.405***'))

p6
summary(lm(data6.AN$rrn~log(data6.AN$level3+1)))
summary(lm(data6.AU$rrn~log(data6.AU$level3+1)))

library(simba)
#To compare the regression coefficient. p<0.05 means significant. Details: ?diffslope
set.seed(123)
diffslope(data6.AN$rrn, (data6.AN$level3+1), data6.AU$rrn, (data6.AU$level3+1))


#ARG
data7 <- read.csv("data/S.ARG.csv",row.names=1,header=T) 
data7.AN <- data7[1:48,]
data7.AU <- data7[49:96,]

p7 <-ggplot()+
  geom_jitter(data=data7, aes(y=Resfam, x= log(level3+1), color=type1),size=10,alpha=0.8)+
  geom_smooth(data=subset(data7,data7$type=="AN"),aes(x= log(level3+1), y=Resfam,color=type1), method ="lm", span=0.9,linetype=2,se=F,linewidth=2)+
  geom_smooth(data=subset(data7,data7$type=="AU"),aes(x= log(level3+1), y=Resfam,color=type1), method ="lm", span=0.9,linetype=1,linewidth=2)+
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('Ln(N addition rate ('*g~N~m ^-2~yr^-1*'))')) +
  ylab(bquote('S.ARG'))+
  scale_color_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        plot.margin = unit(c(0,1,0,1),"cm"),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+ 
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R)["AN"] ^2~'= 0.006'^ns~~~~~italic(R)["AU"] ^2~'= 0.536***'))

p7
summary(lm(data7.AN$Resfam~log(data7.AN$level3+1)))
summary(lm(data7.AU$Resfam~log(data7.AU$level3+1)))

library(simba)
#To compare the regression coefficient. p<0.05 means significant. Details: ?diffslope
set.seed(123)
diffslope(data7.AN$Resfam, (data7.AN$level3+1), data7.AU$Resfam, (data7.AU$level3+1))


#average genome size
data8<-read.csv("data/genomesize.csv",row.names=1,head=T)
data8.AN <- data8[1:48,]
data8.AU <- data8[49:96,]

p8 <-ggplot()+
  geom_jitter(data=data8, aes(y=Ags, x= log(level3+1), color=type1),size=10,alpha=0.8)+
  geom_smooth(data=subset(data8,data8$type=="AN"),aes(x= log(level3+1), y=Ags,color=type1), method ="lm", span=0.9,linewidth=2)+
  geom_smooth(data=subset(data8,data8$type=="AU"),aes(x= log(level3+1), y=Ags,color=type1), method ="lm",se=F, linetype=2,linewidth=2)+
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('Ln(N addition rate ('*g~N~m ^-2~yr^-1*'))')) +
  ylab(bquote('Average Genome Size (Mbp)'))+
  scale_color_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        plot.margin = unit(c(0,1,0,1),"cm"),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+ 
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R)["AN"] ^2~'= 0.228***'~~~~~italic(R)["AU"] ^2~'= -0.022'^ns))

p8
summary(lm(data8.AN$Ags~log(data8.AN$level3+1)))
summary(lm(data8.AU$Ags~log(data8.AU$level3+1)))

library(simba)
#To compare the regression coefficient. p<0.05 means significant. Details: ?diffslope
set.seed(123)
diffslope(data8.AN$Ags, (data8.AN$level3+1), data8.AU$Ags, (data8.AU$level3+1))


#S.KO
p9 <-ggplot()+
  geom_jitter(data=env.amp.N, aes(y=ko_shannon_mgm, x= log(level3+1), color=type1),size=10,alpha=0.8)+
  geom_smooth( data=env.amp.N,aes(y=ko_shannon_mgm, x= log(level3+1),color=type1),method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('Ln(N addition rate ('*g~N~m ^-2~yr^-1*'))')) +
  ylab(bquote('S.KO'))+
  scale_color_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        plot.margin = unit(c(0,1,0,1),"cm"),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+ 
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R)["AN"] ^2~'= 0.498***'~~~~~italic(R)["AU"] ^2~'= 0.652***'))

p9
summary(lm(env.amp.AN$ko_shannon_mgm~log(env.amp.AN$level3+1)))
summary(lm(env.amp.AU$ko_shannon_mgm~log(env.amp.AU$level3+1)))

library(simba)
#To compare the regression coefficient. p<0.05 means significant. Details: ?diffslope
set.seed(123)
diffslope(env.amp.AN$ko_shannon_mgm, (env.amp.AN$level3+1), env.amp.AU$ko_shannon_mgm, (env.amp.AU$level3+1))

#GC
data10<-read.csv("data/GC_content.csv",row.names=1,head=T)
data10.AN <- data10[1:48,]
data10.AU <- data10[49:96,]

p10 <-ggplot()+
  geom_jitter(data=data10, aes(y=GC, x= log(level3+1), color=type1),size=10,alpha=0.8)+
  geom_smooth( data=data10,aes(y=GC, x= log(level3+1),color=type1),method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('Ln(N addition rate ('*g~N~m ^-2~yr^-1*'))')) +
  ylab(bquote('GC Content (%)'))+
  scale_color_manual(values = c("royalblue1","#F24F50", "gold", "#00ff00"))+ 
  scale_shape_manual(values=c(19, 19, 19, 1))+
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        plot.margin = unit(c(0,1,0,1),"cm"),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+ 
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R)["AN"] ^2~'= 0.464***'~~~~~italic(R)["AU"] ^2~'= 0.382***'))

p10
summary(lm(data10.AN$GC~log(data10.AN$level3+1)))
summary(lm(data10.AU$GC~log(data10.AU$level3+1)))

library(simba)
#To compare the regression coefficient. p<0.05 means significant. Details: ?diffslope
set.seed(123)
diffslope(data10.AN$GC, (data10.AN$level3+1), data10.AU$GC, (data10.AU$level3+1))

p <- 
  p1 + p2 + p4 + p3 + p9 + p8+ p7+ p5 + p6 + p10 +
  plot_layout(guides = "collect",ncol=5,nrow=2) +
  theme(legend.position='right')

p
ggsave("Fig.S12a-j.pdf", p,height=13,width =35,limitsize = FALSE )





#Fig.S12k-t
#R4.0.4
library(rfPermute)
library(randomForest)
library(rfUtilities)

fun_env<-read.csv("data/microbial traits_random forest.csv",head=T,row.names=1)

#S.Bacteria
fun.rf <- randomForest(fun_env$bac_shannon~IN+pH+plantrichness+NH3g,
                       data=fun_env,num.cores = 4, nrep = 1000)##ntree=500 mtry = 6 default
rf.significance(fun.rf,fun_env[,c(14:17)])##change the combined table with different columns
fun.rf##the details of model

fun.rfP1 <- rfPermute( fun_env$bac_shannon~IN+pH+plantrichness+NH3g,
                       data=fun_env,num.cores = 4, nrep = 1000)
importance(fun.rfP1)
plotImportance(fun.rfP1,imp.type = NULL,plot = TRUE,alpha = 0.05,ranks = TRUE,ylab = "Increase in MSE (%)",scale = TRUE,
               sig.only = FALSE)

#S.KO
fun_env<-read.csv("data/microbial traits_random forest.csv",head=T,row.names=1)

fun.rf <- randomForest(fun_env$KO_shannon~IN+pH+plantrichness+NH3g,
                       data=fun_env,num.cores = 4, nrep = 1000)##ntree=500 mtry = 6 default
rf.significance(fun.rf,fun_env[,c(14:17)])##change the combined table with different columns
fun.rf##the details of model

fun.rfP1 <- rfPermute( fun_env$KO_shannon~IN+pH+plantrichness+NH3g,
                       data=fun_env,num.cores = 4, nrep = 1000)
importance(fun.rfP1)
plotImportance(fun.rfP1,imp.type = NULL,plot = TRUE,alpha = 0.05,ranks = TRUE,ylab = "Increase in MSE (%)",scale = TRUE,
               sig.only = FALSE)


#S.Phage
fun_env<-read.csv("data/microbial traits_random forest.csv",head=T,row.names=1)

fun.rf <- randomForest(fun_env$vir_shannon_mgm~IN+pH+plantrichness+NH3g,
                       data=fun_env,num.cores = 4, nrep = 1000)##ntree=500 mtry = 6 default
rf.significance(fun.rf,fun_env[,c(14:17)])##change the combined table with different columns
fun.rf##the details of model

fun.rfP1 <- rfPermute( fun_env$vir_shannon_mgm~IN+pH+plantrichness+NH3g,
                       data=fun_env,num.cores = 4, nrep = 1000)
importance(fun.rfP1)
plotImportance(fun.rfP1,imp.type = NULL,plot = TRUE,alpha = 0.05,ranks = TRUE,ylab = "Increase in MSE (%)",scale = TRUE,
               sig.only = FALSE)


###S.FUngi
fun_env<-read.csv("data/microbial traits_random forest.csv",head=T,row.names=1)

fun.rf <- randomForest(fun_env$fun_shannon~IN+pH+plantrichness+NH3g,
                       data=fun_env,num.cores = 4, nrep = 1000)##ntree=500 mtry = 6 default
rf.significance(fun.rf,fun_env[,c(14:17)])##change the combined table with different columns
fun.rf##the details of model

fun.rfP1 <- rfPermute( fun_env$fun_shannon~IN+pH+plantrichness+NH3g,
                       data=fun_env,num.cores = 4, nrep = 1000)
importance(fun.rfP1)
plotImportance(fun.rfP1,imp.type = NULL,plot = TRUE,alpha = 0.05,ranks = TRUE,ylab = "Increase in MSE (%)",scale = TRUE,
               sig.only = FALSE)


#HGT
fun_env<-read.csv("data/microbial traits_random forest.csv",head=T,row.names=1)

fun.rf <- randomForest(fun_env$HGT_per_1K_assembled_genes~IN+pH+plantrichness+NH3g,
                       data=fun_env,num.cores = 4, nrep = 1000)##ntree=500 mtry = 6 default
rf.significance(fun.rf,fun_env[,c(14:17)])##change the combined table with different columns
fun.rf##the details of model

fun.rfP1 <- rfPermute( fun_env$HGT_per_1K_assembled_genes~IN+pH+plantrichness+NH3g,
                       data=fun_env,num.cores = 4, nrep = 1000)
importance(fun.rfP1)
plotImportance(fun.rfP1,imp.type = NULL,plot = TRUE,alpha = 0.05,ranks = TRUE,ylab = "Increase in MSE (%)",scale = TRUE,
               sig.only = FALSE)


#growth_rate
fun_env<-read.csv("data/microbial traits_random forest.csv",head=T,row.names=1)

fun.rf <- randomForest(fun_env$growth_rate~IN+pH+plantrichness+NH3g,
                       data=fun_env,num.cores = 4, nrep = 1000)##ntree=500 mtry = 6 default
rf.significance(fun.rf,fun_env[,c(14:17)])##change the combined table with different columns
fun.rf##the details of model

fun.rfP1 <- rfPermute( fun_env$growth_rate~IN+pH+plantrichness+NH3g,
                       data=fun_env,num.cores = 4, nrep = 1000)
importance(fun.rfP1)
plotImportance(fun.rfP1,imp.type = NULL,plot = TRUE,alpha = 0.05,ranks = TRUE,ylab = "Increase in MSE (%)",scale = TRUE,
               sig.only = FALSE)


#Ags
fun_env<-read.csv("data/microbial traits_random forest.csv",head=T,row.names=1)

fun.rf <- randomForest(fun_env$Ags~IN+pH+plantrichness+NH3g+lysogeny_phage,
                       data=fun_env,num.cores = 4, nrep = 1000)##ntree=500 mtry = 6 default
rf.significance(fun.rf,fun_env[,c(14:17,19)])##change the combined table with different columns,¼ÆËãÏÔÖøÐÔ the same result P=0.001,R2=0.399
fun.rf##the details of model

fun.rfP1 <- rfPermute( fun_env$Ags~IN+pH+plantrichness+NH3g+lysogeny_phage,
                       data=fun_env,num.cores = 4, nrep = 1000)
importance(fun.rfP1)
plotImportance(fun.rfP1,imp.type = NULL,plot = TRUE,alpha = 0.05,ranks = TRUE,ylab = "Increase in MSE (%)",scale = TRUE,
               sig.only = FALSE)


###rrn
fun_env<-read.csv("data/microbial traits_random forest.csv",head=T,row.names=1)

fun.rf <- randomForest(fun_env$rrn~IN+pH+plantrichness+NH3g,
                       data=fun_env,num.cores = 4, nrep = 1000)##ntree=500 mtry = 6 default
rf.significance(fun.rf,fun_env[,c(14:17)])##change the combined table with different columns
fun.rf##the details of model

fun.rfP1 <- rfPermute( fun_env$rrn~IN+pH+plantrichness+NH3g,
                       data=fun_env,num.cores = 4, nrep = 1000)
importance(fun.rfP1)
plotImportance(fun.rfP1,imp.type = NULL,plot = TRUE,alpha = 0.05,ranks = TRUE,ylab = "Increase in MSE (%)",scale = TRUE,
               sig.only = FALSE)


#GC
fun_env<-read.csv("data/microbial traits_random forest.csv",head=T,row.names=1)

fun.rf <- randomForest(fun_env$GC~IN+pH+plantrichness+NH3g,
                       data=fun_env,num.cores = 4, nrep = 1000)##ntree=500 mtry = 6 default
rf.significance(fun.rf,fun_env[,c(14:17)])##change the combined table with different columns
fun.rf##the details of model

fun.rfP1 <- rfPermute( fun_env$GC~IN+pH+plantrichness+NH3g,
                       data=fun_env,num.cores = 4, nrep = 1000)
importance(fun.rfP1)
plotImportance(fun.rfP1,imp.type = NULL,plot = TRUE,alpha = 0.05,ranks = TRUE,ylab = "Increase in MSE (%)",scale = TRUE,
               sig.only = FALSE)


###Resfam
fun_env<-read.csv("data/microbial traits_random forest.csv",head=T,row.names=1)

fun.rf <- randomForest(fun_env$Resfam~IN+pH+plantrichness+NH3g,
                       data=fun_env,num.cores = 4, nrep = 1000)##ntree=500 mtry = 6 default
rf.significance(fun.rf,fun_env[,c(14:17)])##change the combined table with different columns
fun.rf##the details of model

fun.rfP1 <- rfPermute( fun_env$Resfam~IN+pH+plantrichness+NH3g,
                       data=fun_env,num.cores = 4, nrep = 1000)
importance(fun.rfP1)
plotImportance(fun.rfP1,imp.type = NULL,plot = TRUE,alpha = 0.05,ranks = TRUE,ylab = "Increase in MSE (%)",scale = TRUE,
               sig.only = FALSE)



####Fig.S13
dat2<-read.csv("data/Shotgun.16SrRNA.diversity.csv",head=T,row.names=1)
dat2.AN <- dat2[1:48,]
dat2.AU <- dat2[49:96,]

p2 <-ggplot()+
  geom_jitter(data=dat2, aes(y=Shannon, x= log(level3+1), color=type1),size=10,alpha=0.8)+
  geom_smooth( data=dat2,aes(y=Shannon, x= log(level3+1),color=type1),method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('Ln(N addition rate ('*g~N~m ^-2~yr^-1*'))')) +
  ylab(bquote('S.Bacteria from Metagenomic Data'))+
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
  ggtitle(expression(~italic(R)["AN"] ^2~'= 0.598***'~~~~~italic(R)["AU"] ^2~'= 0.291***'))
p2
#summary(lm(dat2.AN$Shannon~(dat2.AN$level3)))
#summary(lm(dat2.AU$Shannon~(dat2.AU$level3)))


ggsave("Fig.S13.pdf", p2,height=7,width =8,limitsize = FALSE )

library(simba)

set.seed(123)
diffslope(dat2.AN$Shannon, log(dat2.AN$level3+1), dat2.AU$Shannon, log(dat2.AU$level3+1))




####Fig.S14
#Fig.S14a-c
data<-read.csv("data/lysogenic phage.csv",row.names=1,head=T)
data.AN <- data[1:48,]
data.AU <- data[49:96,]


p4 <-ggplot()+
  geom_jitter(data=data, aes(y=LPRA, x= log(level3+1), color=type1),size=10,alpha=0.8)+
  geom_smooth(data=subset(data,data$type=="AN"),aes(x= log(level3+1), y=LPRA,color=type1), method ="lm", se=F,linetype=2,linewidth=2)+
  geom_smooth(data=subset(data,data$type=="AU"),aes(x= log(level3+1), y=LPRA,color=type1), method ="lm", span=0.9,linetype=1,linewidth=2)+
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('Ln(N addition rate ('*g~N~m ^-2~yr^-1*'))')) +
  ylab(bquote('Relative Abunbdance of Lysogenic Phage (%)'))+
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
  ggtitle(expression(~italic(R)["AN"] ^2~'= -0.021'^ns~~~~~italic(R)["AU"] ^2~'= 0.300***'))

p4
summary(lm(data.AN$LPRA~log(data.AN$level3+1)))
summary(lm(data.AU$LPRA~log(data.AU$level3+1)))

library(simba)
set.seed(123)

diffslope(data.AN$LPRA, log(data.AN$level3+1), data.AU$LPRA, log(data.AU$level3+1))



p5 <-ggplot()+
  geom_jitter(data=data, aes(y=lysogeny_phage_contigs, x= log(level3+1), color=type1),size=10,alpha=0.8)+
  geom_smooth(data=subset(data,data$type=="AN"),aes(x= log(level3+1), y=lysogeny_phage_contigs,color=type1), method ="lm", se=T,linetype=1,linewidth=2)+
  geom_smooth(data=subset(data,data$type=="AU"),aes(x= log(level3+1), y=lysogeny_phage_contigs,color=type1), method ="lm", span=0.9,linetype=1,linewidth=2)+
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('Ln(N addition rate ('*g~N~m ^-2~yr^-1*'))')) +
  ylab(bquote('The Number of Lysogenic Phage contigs'))+
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
  ggtitle(expression(~italic(R)["AN"] ^2~'= 0.565***'~~~~~italic(R)["AU"] ^2~'= 0.526***'))

p5
summary(lm(data.AN$lysogeny_phage_contigs~log(data.AN$level3+1)))
summary(lm(data.AU$lysogeny_phage_contigs~log(data.AU$level3+1)))

library(simba)
set.seed(123)

diffslope(data2.AN$lysogeny_phage_contigs, log(data2.AN$level3+1), data2.AU$lysogeny_phage_contigs, log(data2.AU$level3+1))



data3<-read.csv("data/number of phage.csv",row.names=1,head=T)
data3.AN <- data3[1:48,]
data3.AU <- data3[49:96,]

p6 <-ggplot()+
  geom_jitter(data=data3, aes(y=npc, x= log(level3+1), color=type1),size=10,alpha=0.8)+
  geom_smooth(data=subset(data3,data3$type=="AN"),aes(x= log(level3+1), y=npc,color=type1), method ="lm", span=0.9,linetype=1,linewidth=2)+
  geom_smooth(data=subset(data3,data3$type=="AU"),aes(x= log(level3+1), y=npc,color=type1), method ="lm", span=0.9,linetype=1,linewidth=2)+
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('Ln(N addition rate ('*g~N~m ^-2~yr^-1*'))')) +
  ylab(bquote('The Number of Phage Contigs'))+
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
  ggtitle(expression(~italic(R)["AN"] ^2~'= 0.628***'~~~~~italic(R)["AU"] ^2~'= 0.540***'))

p6
summary(lm(data3.AN$npc~log(data3.AN$level3+1)))
summary(lm(data3.AU$npc~log(data3.AU$level3+1)))

library(simba)
set.seed(123)

diffslope(data3.AN$npc, log(data3.AN$level3+1), data3.AU$npc, log(data3.AU$level3+1))


p <- 
  p6 + p5 + p4 + 
  plot_layout(guides = "collect",ncol=4,nrow=1) +
  theme(legend.position='right')

p
ggsave("Fig.S14a-c.pdf", p,height=8,width =27,limitsize = FALSE )



#Fig.S14d-e
#number of lysogeny_phage
fun_env<-read.csv("data/microbial traits_random forest.csv",head=T,row.names=1)

fun.rf <- randomForest(fun_env$lysogeny_phage~IN+pH+plantrichness+NH3g,
                       data=fun_env,num.cores = 4, nrep = 1000)##ntree=500 mtry = 6 default
rf.significance(fun.rf,fun_env[,c(14:17)])##change the combined table with different columns
fun.rf##the details of model

fun.rfP1 <- rfPermute( fun_env$lysogeny_phage~IN+pH+plantrichness+NH3g,
                       data=fun_env,num.cores = 4, nrep = 1000)
importance(fun.rfP1)
plotImportance(fun.rfP1,imp.type = NULL,plot = TRUE,alpha = 0.05,ranks = TRUE,ylab = "Increase in MSE (%)",scale = TRUE,
               sig.only = FALSE)


#number pf phage
fun_env<-read.csv("data/microbial traits_random forest.csv",head=T,row.names=1)

fun.rf <- randomForest(fun_env$phage~IN+pH+plantrichness+NH3g,
                       data=fun_env,num.cores = 4, nrep = 1000)##ntree=500 mtry = 6 default
rf.significance(fun.rf,fun_env[,c(14:17)])##change the combined table with different columns
fun.rf##the details of model

fun.rfP1 <- rfPermute( fun_env$phage~IN+pH+plantrichness+NH3g,
                       data=fun_env,num.cores = 4, nrep = 1000)
importance(fun.rfP1)
plotImportance(fun.rfP1,imp.type = NULL,plot = TRUE,alpha = 0.05,ranks = TRUE,ylab = "Increase in MSE (%)",scale = TRUE,
               sig.only = FALSE)




