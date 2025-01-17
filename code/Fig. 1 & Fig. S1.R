#Fig.1 & Fig. S1
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



env.amp.N<-read.csv("data/env_data.csv",head=T,row.names = 1)
env.amp.AN <- env.amp.N[1:48,]
env.amp.AU <- env.amp.N[49:96,]

#####Fig.1
#IN
p1 <-ggplot()+
  geom_jitter(data=env.amp.N, aes(y=log(IN), x= log(level3+1), color=type1),size=10,alpha=0.8)+
  geom_smooth(data=subset(env.amp.N,env.amp.N$type1=="AN"),aes(x= log(level3+1), y=log(IN),color=type1), method ="lm", span=0.9,linewidth=2)+
  geom_smooth(data=subset(env.amp.N,env.amp.N$type1=="AU"),aes(x= log(level3+1), y=log(IN),color=type1), method ="lm",linetype=1,linewidth=2)+
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('Ln(N addition rate ('*g~N~m ^-2~yr^-1*'))')) +
  ylab(bquote(Ln*(IN)~~(mg~kg ^-1)))+
  scale_color_manual(values = c("royalblue1","#F24F50"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+ 
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R)["AN"] ^2~'= 0.580***'~~~~~italic(R)["AU"] ^2~'= 0.420***'))

p1

summary(lm(env.amp.AN$IN~log(env.amp.AN$level3+1)))
summary(lm(env.amp.AU$IN~log(env.amp.AU$level3+1)))

library(simba)

#To compare the regression coefficient£¬p<0.05 means significant. Details: ?diffslope
set.seed(123)
diffslope(env.amp.AN$IN, (env.amp.AN$level3+1), env.amp.AU$IN, (env.amp.AU$level3+1))


#pH
p2 <-ggplot()+
  geom_jitter(data=env.amp.N, aes(y=pH, x= log(level3+1), color=type1),size=10,alpha=0.8)+
  geom_smooth( data=env.amp.N,aes(y=pH, x= log(level3+1),color=type1),method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('Ln(N addition rate ('*g~N~m ^-2~yr^-1*'))')) +
  ylab(bquote('Soil pH'))+
  scale_color_manual(values = c("royalblue1","#F24F50"))+ 
  scale_y_continuous(limits = c(5, 8))+
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+ 
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R)["AN"] ^2~'= 0.604***'~~~~~italic(R)["AU"] ^2~'= 0.539***'))

p2

summary(lm(env.amp.AN$pH~log(env.amp.AN$level3+1)))
summary(lm(env.amp.AU$pH~log(env.amp.AU$level3+1)))

set.seed(123)
diffslope(env.amp.AN$pH, (env.amp.AN$level3+1), env.amp.AU$pH, (env.amp.AU$level3+1))


#NH3
p3 <-ggplot()+
  geom_jitter(data=env.amp.N, aes(y=log(NH3g+1), x= log(level3+1), color=type1),size=10,alpha=0.8)+
  geom_smooth(data=subset(env.amp.N,env.amp.N$type1=="AN"),aes(x= log(level3+1), y=log(NH3g+1),color=type1), method ="lm", linetype=1,linewidth=2)+
  geom_smooth(data=subset(env.amp.N,env.amp.N$type1=="AU"),aes(x= log(level3+1), y=log(NH3g+1),color=type1), method ="lm",linetype=1,linewidth=2)+
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('Ln(N addition rate ('*g~N~m ^-2~yr^-1*'))')) +
  ylab(bquote(Ln*(NH[3])~~(g~N~m ^-2)))+
  scale_color_manual(values = c("royalblue1","#F24F50"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+ 
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R)["AN"] ^2~'= 0.267***'~~~~~italic(R)["AU"] ^2~'= 0.721***'))

p3

summary(lm(env.amp.AN$NH3g~log(env.amp.AN$level3+1)))
summary(lm(env.amp.AU$NH3g~log(env.amp.AU$level3+1)))

set.seed(123)
diffslope(env.amp.AN$NH3g, (env.amp.AN$level+1), env.amp.AU$NH3g, (env.amp.AU$level+1))


#plant richness
p4 <-ggplot()+
  geom_jitter(data=env.amp.N, aes(y=sr, x= log(level3+1), color=type1),size=10,alpha=0.8)+
  geom_smooth(data=subset(env.amp.N,env.amp.N$type1=="AN"),aes(x= log(level3+1), y=sr,color=type1), method ="lm", linewidth=2,span=0.9)+
  geom_smooth(data=subset(env.amp.N,env.amp.N$type1=="AU"),aes(x= log(level3+1), y=sr,color=type1), method ="lm",linetype=1,linewidth=2)+
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('Ln(N addition rate ('*g~N~m ^-2~yr^-1*'))')) +
  ylab(bquote('Species Richness of Plant'))+
  scale_color_manual(values = c("royalblue1","#F24F50"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+ 
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R)["AN"] ^2~'= 0.378***'~~~~~italic(R)["AU"] ^2~'= 0.302***'))

p4
summary(lm(env.amp.AN$sr~log(env.amp.AN$level3+1)))
summary(lm(env.amp.AU$sr~log(env.amp.AU$level3+1)))

set.seed(123)
diffslope(env.amp.AN$sr, (env.amp.AN$level3+1), env.amp.AU$sr, (env.amp.AU$level3+1))

#all
p <- 
  p1 +p2 + p3 + p4 +
  plot_layout(guides = "collect",ncol=2,nrow=2) +
  theme(legend.position='right',plot.margin = unit(c(0,1,0,2),"cm"))
  
p

ggsave("Fig.1.pdf", p, height=13,width =14,limitsize = FALSE )


#####Fig.S1
NH3<-read.csv("data/NH3-meta.csv",head=T,row.names = 1)
NH3.AN <- NH3[1:30,]
NH3.AU <- NH3[31:60,]

p5 <-ggplot()+
  geom_jitter(data=NH3, aes(y=log(NH3g+1), x= log(level+1), color=type),size=10,alpha=0.8)+
  geom_smooth(data=subset(NH3,NH3$type=="AN"),aes(x= log(level+1), y=log(NH3g+1),color=type), method ="lm", linetype=1,linewidth=2)+
  geom_smooth(data=subset(NH3,NH3$type=="AU"),aes(x= log(level+1), y=log(NH3g+1),color=type), method ="lm",linetype=1,linewidth=2)+
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('Ln(N addition rate ('*g~N~m ^-2~yr^-1*'))')) +
  ylab(bquote(Ln*(NH[3])~~(g~N~m ^-2)))+
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
  ggtitle(expression(~italic(R)["AN"] ^2~'= 0.135*'~~~~~italic(R)["AU"] ^2~'= 0.673***'))

p5
summary(lm(NH3.AN$NH3g~log(NH3.AN$level+1)))
summary(lm(NH3.AU$NH3g~log(NH3.AU$level+1)))

ggsave("Fig.S1.pdf", p5,height=7,width =8,limitsize = FALSE )


library(simba)
set.seed(123)
diffslope(NH3.AN$NH3g, (NH3.AN$level+1), NH3.AU$NH3g, (NH3.AU$level+1))