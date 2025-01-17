#Fig.2 & 
#Fig. S2-S8
#Fig,S22-S30



####Fig.2
#Fig.2a
#R3.6.2
ibrary(vegan)
library(omicade4)


functional.traits<-read.csv("data/functional traits_MCOA.csv",head=T,row.names=1)
functional<-read.csv("data/functional potential_MCOA.csv",head=T,row.names=1)
bacteria<-read.csv("data/bacterial taxa_MCOA.csv",head=T,row.names=1)

#log10
functional.traits2 <- log(functional.traits,10)
functional2<-log(functional,10)
bacteria2<-log(bacteria,10)
#z-score
functional.traits3 <- scale(functional.traits2)
functional3 <-scale(functional2)
bacteria3 <- scale(bacteria2)

#construct list
lista <- list(functional.traits3,functional3,bacteria3)
names(lista)<-c("functional.traits","functional","bacteria")

head(lista$functional.traits[1:6]) #Agilent platform.
head(lista$functional[1:6])        #H-GU133 platform
head(lista$bacteria[1:6])    #H-GU133 plus 2.0 platform

# MCoIA analysis，details ?mcia
mcoin <- mcia(lista, cia.nf = 3, cia.scan = FALSE, nsc = TRUE)
mcoin
summary(mcoin)

#main results
names(mcoin$mcoa)

#pseudo-eigenvalue of MCoIA
mcoin$mcoa$pseudoeig

#RV value
mcoin$mcoa$RV

mcoin$mcoa
sample_type <- colnames(lista$functional.traits)
sample_type <- sapply(strsplit(sample_type, split="\\."), function(x) x[1])
sample_type


##plot
library(vegan)
library(ggplot2)
df_points <- as.data.frame(mcoin$mcoa$SynVar)
df_points$samples <- sample_type
df_points$Type <-c("AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN",
                   "AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN",
                   "AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU",
                   "AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU")
library(plyr)
df <-df_points[,c("SynVar1", "SynVar2", "Type")]
find_hull <- function(df_points) df_points[chull(df_points$SynVar1, df_points$SynVar2), ]
hulls <- ddply(df, "Type", find_hull )###border

#准备构建旋转后的坐标
#write.csv(df_points,"MCOA_amplicon.NMDS.csv")


####读入新的旋转后的坐标
df_points2 <- read.csv("data/MCOA_amplicon.NMDS.csv",header = T,row.names = 1)
df <-df_points2[,c("SynVar1", "SynVar2", "Type")]
find_hull <- function(df_points2) df_points2[chull(df_points2$SynVar1, df_points2$SynVar2), ]

hulls <- ddply(df, "Type", find_hull )###border




#读取环境数据
env <- read.csv("data/env.csv",head=T,row.names=1)

mcoa_env<-envfit(df_points2[,c(6,7)]~.,data=env,perm=999,choices=c(1,2),display='sites')


mcoa_env_adj<-mcoa_env
mcoa_env_adj$vectors$pvals<-p.adjust(mcoa_env_adj$vectors$pvals,method='bonferroni')
mcoa_env_adj

mcoa_env2<-data.frame(cbind(mcoa_env_adj$vectors$arrows,mcoa_env_adj$vectors$r,mcoa_env_adj$vectors$pvals))
names(mcoa_env2)<-c('SynVar1','SynVar2','r2','p.adj')
#write.csv(mcoa_env2,'mcoa_env.csv',quote=FALSE)

mcoa_env3<- read.csv("mcoa_env.csv",head=T,row.names=1)



p1<-ggplot(data=df_points2, aes(x= SynVar1 , y= SynVar2))+ 
  
  geom_point(aes(colour=samples),size=6,alpha=0.7)+ 
  geom_polygon(data=hulls, alpha=.1,aes(fill = Type))+
  labs(x = "MCOA1", y = "MCOA2")+
  scale_colour_manual(values=c("#C7E2FF","lightblue1","#95E9FF","#5EB5FF","#339CFF","blue",
                               "#FFE1FF","#FFC8FF","#FFAFFF","#FF63FF","#FF6371","red"))+
  scale_fill_manual(values=c("blue","red"))+
  geom_segment(data=mcoa_env3,aes(x=0,y=0,xend=SynVar1*2,yend=SynVar2*2),arrow=arrow(length=unit(0.3,'cm')),size=0.9,color='black')+
  geom_text(data=mcoa_env3,aes(SynVar1*2.5,SynVar2*2.5,label=group),color='black',size=3)+
  scale_x_continuous(limits = c(-3, 5))+
  theme_bw()+
  guides(color=guide_legend(title= "Treatment"))+
  theme(legend.position = 'left',legend.direction = 'vertical',
        legend.title = element_text(colour="black", size=13, face="bold"),
        legend.text = element_text(colour="black", size=12),
        axis.text=element_text(size=15,face="bold",colour="black"),
        axis.title=element_text(size=17,face="bold"))+
  geom_vline(xintercept = 0, color = 'gray34', linetype = 1, linewidth = 1.1) +
  geom_hline(yintercept = 0, color = 'gray34', linetype = 1, linewidth = 1.1) 


p1
ggsave("Fig.2a.pdf", p1,height=5,width =6,limitsize = FALSE )


####Fig.S22
#Fig.S22a
eigen_value <- as.data.frame(mcoin$mcoa$pseudoeig)
eigen_value$number <-c(1:20)
names(eigen_value)[1:2] <- c('eigenvalues', 'dimension')
#write.csv(eigen_value,"variation and dimension_amplicon.csv")

eigen_value2<- read.csv("data/variation and dimension_amplicon.csv",head=T,row.names=1)

p2 <- ggplot(data = eigen_value2, aes(x =dimension, y = variation100))+
  geom_bar(stat = "identity", width = 0.8, colour = "black", linewidth = 0.8,fill = "thistle3", alpha = 1) +
  ylim(0, 80) + # y_scale
  scale_x_continuous(breaks = seq(0, 20, by = 1),minor_breaks=NULL)+
  labs(y = "% of variation", x = "MCOA dimension")+
  theme_bw()+
  theme(axis.title = element_text(size = 17, face = "plain", color = "black"), 
        axis.text = element_text(size = 12, face = "plain", color = "black"))

p2
ggsave("Fig.S22a.pdf", p2,height=5,width =6,limitsize = FALSE )

#Fig.S22d
df_cov <- as.data.frame(mcoin$mcoa$cov2)
df_cov$Dataset <- c("Microbial Traits","Functional Potential","Bacterial Taxa")

#write.csv(df_cov,"contribution_amplicon.csv")


con2 <- read.csv("data/contribution_amplicon.csv",head=T,row.names=1)

#MCOA1
p1 <- ggplot(con2,aes(Dataset,mcoa1))+
  geom_col(aes(fill=Dataset))+
  scale_fill_manual(values=c("#7552cc","royalblue1","#00c16e","#7552cc"),
                    breaks=c("Bacterial.taxa","Functional.potential","Microbiome.parameters"),
                    labels=c("Bacterial taxa","Functional.potential",  "Microbiome parameters"))+
  labs(y = "Contribution to MCOA1", x = "Dateset Category")+
  theme_bw()+
  guides(fill=guide_legend(title=NULL,reverse=TRUE))+
  theme(
    panel.background = element_rect(fill = "gray99" ),
    legend.title = element_text(colour="black", size=13, face="bold"),
    legend.text = element_text(colour="black", size=12),
    axis.text=element_text(size=13,face="bold",colour="black"),
    axis.title=element_text(size=16,face="bold"),
    legend.position="top")+
  coord_flip()

p1

#MCOA2
p2 <- ggplot(con2,aes(Dataset,mcoa2))+
  geom_col(aes(fill=Dataset))+
  scale_fill_manual(values=c("#7552cc","royalblue1","#00c16e","#7552cc"),
                    breaks=c("Bacterial.taxa","Functional.potential","Microbiome.parameters"),
                    labels=c("Bacterial taxa","Functional.potential",  "Microbiome parameters"))+
  labs(y = "Contribution to MCOA2", x = "Dateset Category")+
  theme_bw()+
  guides(fill=guide_legend(title=NULL,reverse=TRUE))+
  theme(
    panel.background = element_rect(fill = "gray99" ),
    legend.title = element_text(colour="black", size=13, face="bold"),
    legend.text = element_text(colour="black", size=12),
    axis.text=element_text(size=13,face="bold",colour="black"),
    axis.title=element_text(size=16,face="bold"),
    legend.position="top")+
  coord_flip()

p2

#Fig.S22bc
library(rfPermute)
library(randomForest)
library(rfUtilities)

fun_env<-read.csv("data/random_forest.csv",head=T,row.names=1)

###MCOA1
fun.rf <- randomForest(fun_env$MCOA1_amplicon~IN+pH+plantrichness+NH3g,
                       data=fun_env,num.cores = 4, nrep = 1000)
rf.significance(fun.rf,fun_env[,c(7:10)])##change the combined table with different columns
fun.rf##the details of model

fun.rfP1 <- rfPermute( fun_env$MCOA1_amplicon~IN+pH+plantrichness+NH3g,
                       data=fun_env,num.cores = 4, nrep = 1000)
importance(fun.rfP1)
plotImportance(fun.rfP1,imp.type = NULL,plot = TRUE,alpha = 0.05,ranks = TRUE,ylab = "Increase in MSE (%)",scale = TRUE,
               sig.only = FALSE)


###MCOA2
fun.rf2 <- randomForest(fun_env$MCOA2_amplicon~IN+pH+plantrichness+NH3g,
                        data=fun_env,num.cores = 4, nrep = 1000)##ntree=500 mtry = 6 default
rf.significance(fun.rf2,fun_env[,c(7:10)])##change the combined table with different columns
fun.rf2##the details of model

fun.rfP12 <- rfPermute( fun_env$MCOA2_amplicon~IN+pH+plantrichness+NH3g,
                        data=fun_env,num.cores = 4, nrep = 1000)
importance(fun.rfP12)
plotImportance(fun.rfP12,imp.type = NULL,plot = TRUE,alpha = 0.05,ranks = TRUE,ylab = "Increase in MSE (%)",scale = TRUE,
               sig.only = FALSE)


####Fig.2b-d
#R3.6.2
library(ggplot2)

####Fig.2b: bacterial taxa
#MCOA1
dat<-read.csv("data/bacteria taxa_contribution_amplicon.csv",head=T,row.names=1)

#####1轴筛选变量
#拟合各环境变量与一轴的一元回归，并提取各自 R2 和 p 值,

dat_plot <- reshape2::melt(dat[,c(-3,-4)], id = 'SynVar1')

p <- ggplot(dat_plot, aes( SynVar1,value)) +
  geom_point() +
  facet_wrap(~variable, ncol = 3, scale = 'free') +
  geom_smooth(method = 'lm')

p


env <- colnames(dat[,5:21])
R2_adj <- c()
p_value <- c()

for (i in env) {
  fit_stat <- summary(lm(dat[['SynVar1']]~dat[[i]]))  #lm
  R2_adj <- c(R2_adj, fit_stat$adj.r.squared)  #R2
  p_value <- c(p_value, fit_stat$coefficients[2,4])  #p value
}

env_stat <- data.frame(row.names = env, R2_adj, p_value)
env_stat

fit1 <- lm(SynVar1~Acidobacteria+Actinobacteria+Bacteroidetes+Chloroflexi+Armatimonadetes
           +Gemmatimonadetes+Planctomycetes+Verrucomicrobia, data = dat)
summary(fit1)

library(relaimpo)
crf <- calc.relimp(fit1, rela = T, 
                   type = c('lmg', 'last', 'first', 'betasq', 'pratt', 'genizi', 'car'))
crf
plot(crf)

#lmg
con<-read.csv("data/bacteria taxa to MCOA1_amplicon.csv",head=T,row.names=1)
con1<-within(con,{
  type<-factor(type,levels=c("Acidobacteria","Gemmatimonadetes", "Verrucomicrobia","Planctomycetes","Armatimonadetes",
                             "Actinobacteria","Bacteroidetes","Chloroflexi"))
})

ggplot(con1,aes(type,lmg))+
  geom_col(aes(fill=group))+
  scale_fill_manual(values=c("steelblue1","darksalmon","#00c16e","#7552cc"),
                    breaks=c("n", "p"),
                    labels=c("Trait negatively associated with MCOA1", "Trait positively associated with MCOA1"))+
  labs(y = "Contribution to MCOA1", x = "Bacteria taxa")+
  theme_bw()+
  guides(fill=guide_legend(title=NULL,reverse=TRUE))+
  theme(
    panel.background = element_rect(fill = "gray99" ),
    legend.title = element_text(colour="black", size=13, face="bold"),
    legend.text = element_text(colour="black", size=12),
    axis.text=element_text(size=13,face="bold",colour="black"),
    axis.title=element_text(size=16,face="bold"),
    legend.position="top")+
  coord_flip()


#MCOA2
dat<-read.csv("data/bacteria taxa_contribution_amplicon.csv",head=T,row.names=1)

env <- colnames(dat[,5:21])
R2_adj <- c()
p_value <- c()

for (i in env) {
  fit_stat <- summary(lm(dat[['SynVar2']]~dat[[i]]))  #lm
  R2_adj <- c(R2_adj, fit_stat$adj.r.squared)  #R2
  p_value <- c(p_value, fit_stat$coefficients[2,4])  #p value
}

env_stat <- data.frame(row.names = env, R2_adj, p_value)
env_stat  

fit1 <- lm(SynVar2~Acidobacteria+Actinobacteria+Bacteroidetes+Armatimonadetes+Gemmatimonadetes+Planctomycetes
           +Proteobacteria, data = dat)
summary(fit1)

library(relaimpo)
crf <- calc.relimp(fit1, rela = T, 
                   type = c('lmg', 'last', 'first', 'betasq', 'pratt', 'genizi', 'car'))
crf
plot(crf)

##lmg
con<-read.csv("data/bacteria taxa to MCOA2_amplicon.csv",head=T,row.names=1)

con1<-within(con,{
  type<-factor(type,levels=c("Actinobacteria","Acidobacteria", "Armatimonadetes","Planctomycetes","Gemmatimonadetes",  "Proteobacteria","Bacteroidetes"
  ))
})

ggplot(con1,aes(type,lmg))+
  geom_col(aes(fill=group))+
  scale_fill_manual(values=c("steelblue1","darksalmon","#00c16e","#7552cc"),
                    breaks=c("n", "p"),
                    labels=c("Trait negatively associated with MCOA2", "Trait positively associated with MCOA2"))+
  labs(y = "Contribution to MCOA2", x = "Bacteria taxa")+
  theme_bw()+
  guides(fill=guide_legend(title=NULL,reverse=TRUE))+
  theme(
    panel.background = element_rect(fill = "gray99" ),
    legend.title = element_text(colour="black", size=13, face="bold"),
    legend.text = element_text(colour="black", size=12),
    axis.text=element_text(size=13,face="bold",colour="black"),
    axis.title=element_text(size=16,face="bold"),
    legend.position="top")+
  coord_flip()


####Fig.2c: funtcional potential
dat<-read.csv("data/1functional potential contribution_amplicon.csv",head=T,row.names=1)

#MCOA1
env <- colnames(dat[,4:30])
R2_adj <- c()
p_value <- c()

for (i in env) {
  fit_stat <- summary(lm(dat[['SynVar1']]~dat[[i]]))  #lm
  R2_adj <- c(R2_adj, fit_stat$adj.r.squared)  #R2
  p_value <- c(p_value, fit_stat$coefficients[2,4])  #p value
}

env_stat <- data.frame(row.names = env, R2_adj, p_value)
env_stat  


fit1 <- lm(SynVar1~YAS5+YAS6+YAS7+YAS9+YAS11+YAS12+YAS14+YAS16+YAS17+YAS20+YAS22+YAS23+YAS24+YAS25+YAS26+YAS27,data = dat)
summary(fit1)

library(relaimpo)
crf <- calc.relimp(fit1, rela = T, type = 'lmg')
crf
plot(crf)


#lmg方法
con<-read.csv("data/functional potential to MCOA1_amplicon.csv",head=T,row.names=1)

con1<-within(con,{
  type<-factor(category,levels=c("Metabolism of cofactors and vitamins","Sulfur metabolism(A)","Nitrogen metabolism(A)","Membrane transport(A)",
                                 "Signal transduction","Carbohydrate metabolism(YA)","Energy metabolism(A)","Carbon fixation pathways in prokaryotes(A)",
                                 "Quorum sensing","Fatty acid biosynthesis(Y)","Transcription","Cell growth and death(Y)","Lipid metabolism(Y)","Folding, sorting and degradation",
                                 "Translation","Replication and repair(YS)"
  ))
})


library(ggplot2)
ggplot(con1,aes(type,lmg))+
  geom_col(aes(fill=group))+
  scale_fill_manual(values=c("steelblue1","darksalmon","#00c16e","#7552cc"),
                    breaks=c("n", "p"),
                    labels=c("Trait negatively associated with MCOA1", "Trait positively associated with MCOA1"))+
  labs(y = "Contribution to MCOA1", x = "Microbiome parameters")+
  theme_bw()+
  guides(fill=guide_legend(title=NULL,reverse=TRUE))+
  #scale_y_continuous( limits=c(-0.09, 0.1))+
  theme(
    panel.background = element_rect(fill = "gray99" ),
    legend.title = element_text(colour="black", size=13, face="bold"),
    legend.text = element_text(colour="black", size=12),
    axis.text=element_text(size=13,face="bold",colour="black"),
    axis.title=element_text(size=16,face="bold"),
    legend.position="top")+
  
  coord_flip()

#MCOA2
dat<-read.csv("data/2functional potential contribution_amplicon.csv",head=T,row.names=1)

env <- colnames(dat[,4:30])
R2_adj <- c()
p_value <- c()

for (i in env) {
  fit_stat <- summary(lm(dat[['SynVar2']]~dat[[i]]))  #lm
  R2_adj <- c(R2_adj, fit_stat$adj.r.squared)  #R2
  p_value <- c(p_value, fit_stat$coefficients[2,4])  #p value
}

env_stat <- data.frame(row.names = env, R2_adj, p_value)
env_stat  


fit1 <- lm(SynVar2~YAS2+YAS3+YAS4+YAS6+YAS8+YAS9+YAS10+YAS12+YAS13+YAS15+YAS16+YAS17+YAS20+YAS22+YAS24,data = dat)
summary(fit1)


library(relaimpo)
crf <- calc.relimp(fit1, rela = T, type = 'lmg')
crf
plot(crf)

#lmg
con<-read.csv("data/functional potential to MCOA2_amplicon.csv",head=T,row.names=1)
con1<-within(con,{
  type<-factor(category,levels=c("Folding, sorting and degradation","Nitrogen metabolism(A)","Metabolism of cofactors and vitamins","Membrane transport(A)",
                                 "Glycan biosynthesis and metabolism(S)","Signal transduction","Carbon fixation pathways in prokaryotes(A)",
                                 "Bacterial chemotaxis(A)","Energy metabolism(A)",
                                 "Biofilm formation(S)","Biosynthesis of other secondary metabolites(S)","Lipopolysaccharide biosynthesis(S)",
                                 "Cell motility(A)",
                                 "Quorum sensing","Exopolysaccharide biosynthesis(S)"))})
library(ggplot2)
ggplot(con1,aes(type,lmg))+
  geom_col(aes(fill=group))+
  scale_fill_manual(values=c("steelblue1","darksalmon","#00c16e","#7552cc"),
                    breaks=c("n", "p"),
                    labels=c("Trait negatively associated with MCOA2", "Trait positively associated with MCOA2"))+
  labs(y = "Contribution to MCOA2", x = "Microbiome parameters")+
  theme_bw()+
  guides(fill=guide_legend(title=NULL,reverse=TRUE))+
  theme(
    panel.background = element_rect(fill = "gray99" ),
    legend.title = element_text(colour="black", size=13, face="bold"),
    legend.text = element_text(colour="black", size=12),
    axis.text=element_text(size=13,face="bold",colour="black"),
    axis.title=element_text(size=16,face="bold"),
    legend.position="top")+
  coord_flip()


####Fig.2d: microbial traits
dat<-read.csv("data/1microbial traits contribution_amplicon.csv",head=T,row.names=1)

#MCOA1
dat_plot <- reshape2::melt(dat, id = 'SynVar1')
p <- ggplot(dat_plot, aes(value, SynVar1)) +
  geom_point() +
  facet_wrap(~variable, ncol = 3, scale = 'free') +
  geom_smooth(method = 'lm')
p

env <- c('bac_shannon', 'KO_shannon', 'vir_shannon_mgm', 'fun_shannon', 'HGT_per_assembled_genes', 
         'growth_rate',"Ags","rrn","Resfam","GC")
R2_adj <- c()
p_value <- c()

for (i in env) {
  fit_stat <- summary(lm(dat[['SynVar1']]~dat[[i]]))  #lm
  R2_adj <- c(R2_adj, fit_stat$adj.r.squared)  #R2
  p_value <- c(p_value, fit_stat$coefficients[2,4])  #p value
}

env_stat <- data.frame(row.names = env, R2_adj, p_value)
env_stat  

fit1 <- lm(SynVar1~bac_shannon+ KO_shannon+vir_shannon_mgm+fun_shannon+HGT_per_assembled_genes+growth_rate+Ags
           +rrn+Resfam+GC,data = dat)
summary(fit1) 

# calc.relimp() for contribution，details ?calc.relimp
library(relaimpo)
crf <- calc.relimp(fit1, rela = T, 
                   type = c('lmg', 'last', 'first', 'betasq', 'pratt', 'genizi', 'car'))
crf
plot(crf)

#lmg
con<-read.csv("data/microbial traits to MCOA1_amplicon.csv",head=T,row.names=1)

con1<-within(con,{
  type<-factor(type,levels=c("GC", "S.Bacteria","S.Fungi","rrn" ,"Average.Genome.Size","S.ARG", "HGT events","S.Phage" ,"S.KO","Growth rate"))
})

p1 <- ggplot(con1,aes(type,lmg))+
  geom_col(aes(fill=group))+
  scale_fill_manual(values=c("steelblue1","darksalmon","#00c16e","#7552cc"),
                    breaks=c("n", "p"),
                    labels=c("Trait negatively associated with MCOA1", "Trait positively associated with MCOA1"))+
  labs(y = "Contribution to MCOA1", x = "Microbial traits")+
  theme_bw()+
  guides(fill=guide_legend(title=NULL,reverse=TRUE))+
  theme(
    panel.background = element_rect(fill = "gray99" ),
    legend.title = element_text(colour="black", size=13, face="bold"),
    legend.text = element_text(colour="black", size=12),
    axis.text=element_text(size=13,face="bold",colour="black"),
    axis.title=element_text(size=16,face="bold"),
    legend.position="top")+
  coord_flip()
p1

##MCOA2
rm(list=ls())
dat<-read.csv("data/2microbial traits contribution_amplicon.csv",head=T,row.names=1)

dat_plot <- reshape2::melt(dat, id = 'SynVar2')

p <- ggplot(dat_plot, aes(value, SynVar2)) +
  geom_point() +
  facet_wrap(~variable, ncol = 3, scale = 'free') +
  geom_smooth(method = 'lm')
p

env <- c('bac_shannon', 'KO_shannon', 'vir_shannon_mgm', 'fun_shannon', 'HGT_per_assembled_genes', 
         'growth_rate',"Ags","rrn","Resfam","GC")
R2_adj <- c()
p_value <- c()

for (i in env) {
  fit_stat <- summary(lm(dat[['SynVar2']]~dat[[i]]))  #lm
  R2_adj <- c(R2_adj, fit_stat$adj.r.squared)  #R2
  p_value <- c(p_value, fit_stat$coefficients[2,4])  # p value
}

env_stat <- data.frame(row.names = env, R2_adj, p_value)
env_stat  

fit1 <- lm(SynVar2~ bac_shannon+KO_shannon+vir_shannon_mgm+fun_shannon+HGT_per_assembled_genes+Ags
           +rrn+Resfam,data = dat)
summary(fit1)  



library(relaimpo)
crf <- calc.relimp(fit1, rela = T, 
                   type = c('lmg', 'last', 'first', 'betasq', 'pratt', 'genizi', 'car'))
crf
plot(crf)

#lmg
con<-read.csv("data/microbial traits to MCOA2_amplicon.csv",head=T,row.names=1)

con1<-within(con,{
  type<-factor(type,levels=c("rrn","GC","S.Fungi","S.Bacteria", "Average.Genome.Size","HGT.events","S.Phage","S.ARG"))
})

p1 <- ggplot(con1,aes(type,lmg))+
  geom_col(aes(fill=group))+
  scale_fill_manual(values=c("steelblue1","darksalmon","#00c16e","#7552cc"),
                    breaks=c("n", "p"),
                    labels=c("Trait negatively associated with MCOA2", "Trait positively associated with MCOA2"))+#自定义颜色
  labs(y = "Contribution to MCOA2", x = "Microbial traits")+
  
  #scale_y_continuous( limits=c(-0.4, 0.2))+
  theme_bw()+
  guides(fill=guide_legend(title=NULL,reverse=TRUE))+
  theme(
    panel.background = element_rect(fill = "gray99" ),
    legend.title = element_text(colour="black", size=13, face="bold"),
    legend.text = element_text(colour="black", size=12),
    axis.text=element_text(size=13,face="bold",colour="black"),
    axis.title=element_text(size=16,face="bold"),
    legend.position="top")+
  coord_flip()
p1


####Fig.S2a
otu <- read.csv("data/16S.amplicon.otutable.csv")
bac.nmds<-metaMDS(vegdist(otu[1:96,4:8542],method = 'bray'),trymax=100)
bac.nmds

#permanova
permanova<-adonis(otu[1:96,4:8542] ~ otu$treat)###perm anova### 
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
#write.csv(df_points,"bacterial_taxa.NMDS.csv")

####读入新的旋转后的坐标
df_points2 <- read.csv("bacterial_taxa.NMDS.csv",header = T,row.names = 1)
df <-df_points2[,c("MDS1", "MDS2", "Type")]
find_hull <- function(df_points2) df_points2[chull(df_points2$MDS1, df_points2$MDS2), ]

hulls <- ddply(df, "Type", find_hull )###border



#plot
p1<-ggplot(data=df_points2, aes(x= MDS1 , y= MDS2))+ 
  geom_point(aes(colour=samples),size=6,alpha=0.7)+ 
  geom_polygon(data=hulls, alpha=.1,aes(fill = Type))+
  labs(x = "NMDS1", y = "NMDS2")+
  scale_colour_manual(values=c("#C7E2FF","lightblue1","#95E9FF","#5EB5FF","#339CFF","blue",
                                        "#FFE1FF","#FFC8FF","#FFAFFF","#FF63FF","#FF6371","red"))+
  scale_fill_manual(values=c("blue","red"))+
  #geom_segment(data=mcoa_env3,aes(x=0,y=0,xend=NMDS1*0.5,yend=NMDS2*0.5),arrow=arrow(length=unit(0.2,'cm')),size=0.8,color='black')+
  #geom_text(data=mcoa_env3,aes(NMDS1*0.6,NMDS2*0.6,label=group),color='black',size=3)+
  ggtitle(expression('Stress= 0.086'~~~~'Treatment:'~italic(R)^2~'= 0.355***'))+
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
ggsave("Fig.S2a.pdf", p1,height=5,width =6,limitsize = FALSE )


####Fig.S2b
#4.2.3
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


data<-read.table("data/16S.ANAU50.phylum.txt",header=TRUE,sep="\t")
combined =data##用于计算显著性

AN50_UP <- data[,c(1,2,4)]
ratio <- data.frame(cbind(AN50_UP,AN50_UP$AN50_up/AN50_UP$total))
dat <- ratio[,c(1,4)]
colnames(dat) <-c("Annotation","AN50_UP_ratio")

merAN<- dat


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
  "Verrucomicrobiota",
  "Planctomycetota",
  "Armatimonadota",
  "Gemmatimonadota",
  "Acidobacteriota",
  "Proteobacteria",
  "Bacteroidota",
  "Chloroflexota",
  "Actinobacteriota"
  
))


samples_new <- sapply(as.character(dat$variable), function(x) strsplit(x, "[_]"))

Timepoint <- {}
for(i in 1:18){
  T <- samples_new[[i]][1]
  Timepoint <- append(Timepoint,T)
}

SampleType <- {}
for(i in 1:18){
  T <- samples_new[[i]][2]
  SampleType <- append(SampleType,T)
}

dat$SampleType = SampleType
dat$Timepoint = Timepoint


combined<-read.table("data/16S.ANAU50.phylum.combined.txt",header=TRUE,sep="\t")

combined$AN50_UP_pvalue <- phyper(combined$AN50_up-1,combined$total,sum(combined$total)-combined$total,sum(combined$AN50_up),
                                  lower.tail = FALSE, log.p = FALSE)


combined$AU50_UP_pvalue <- phyper(combined$AU50_up-1,combined$total,sum(combined$total)-combined$total,sum(combined$AU50_up),
                                  lower.tail = FALSE, log.p = FALSE)


combined[which(combined$AN50_UP_pvalue < 0.05),'sig.AN'] <- 'Significant'
combined[which(combined$AU50_UP_pvalue < 0.05),'sig.AU'] <- 'Significant'

dat[which(dat$Annotation <= 1 ),'sig'] <- 'Non-significant(1)'
dat[which(dat$value > 1),'sig'] <- 'Significant(1)'
dat1 <- subset(dat,Timepoint=="AN50")
dat2 <- subset(dat,Timepoint=="AU50")



ggplot(dat, aes(x=value, y=Annotation,colour=Timepoint,shape=sig)) +
  scale_colour_manual(name="",values = c("royalblue1","#F24F50","#9B4F0F","#6FB98F"))+
  scale_shape_manual(name="",values=c(1,1,19,19))+facet_wrap(~SampleType,scales="free_x",ncol=2)+
  scale_shape_manual(values  = c(21,16), limits = c('Non-significant(1)', 'Significant(1)')) +
  geom_point(size =6,stroke = 1)+theme_bw()+
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


#Fig.S2cd
#R4.0.4
library(rfPermute)
library(randomForest)
library(rfUtilities)

fun_env<-read.csv("data/random_forest.csv",head=T,row.names=1)

###NMDS1
fun.rf <- randomForest(fun_env$NMDS1_species~IN+pH+plantrichness+NH3g,
                       data=fun_env,num.cores = 4, nrep = 1000)
rf.significance(fun.rf,fun_env[,c(7:10)])
fun.rf##the details of model

fun.rfP1 <- rfPermute( fun_env$NMDS1_species~IN+pH+plantrichness+NH3g,
                       data=fun_env,num.cores = 4, nrep = 1000)
importance(fun.rfP1)
plotImportance(fun.rfP1,imp.type = NULL,plot = TRUE,alpha = 0.05,ranks = TRUE,ylab = "Increase in MSE (%)",scale = TRUE,
               sig.only = FALSE)
###NMDS2
fun.rf2 <- randomForest(fun_env$NMDS2_species~IN+pH+plantrichness+NH3g,
                        data=fun_env,num.cores = 4, nrep = 1000)
rf.significance(fun.rf2,fun_env[,c(7:10)])
fun.rf2

fun.rfP12 <- rfPermute( fun_env$NMDS2_species~IN+pH+plantrichness+NH3g,
                        data=fun_env,num.cores = 4, nrep = 1000)
importance(fun.rfP12)
plotImportance(fun.rfP12,imp.type = NULL,plot = TRUE,alpha = 0.05,ranks = TRUE,ylab = "Increase in MSE (%)",scale = TRUE,
               sig.only = FALSE)


#Fig.S2e & Fig.S5
#R3.6.2
library(readxl)
library(psych)
library(tidyr)
library(ggplot2)
library(vegan)
library(colorRamps)
library(ggrepel)
library(grid)

ko.sel0 <-read.csv("data/OTU-NMDS.csv",header=T,row.names=1)
ko.sel.ID0 <-  ko.sel0[, c(1:3)]

ko.sel <-ko.sel0[, c(-1:-3)]
ko.sel2 <- log(ko.sel,10)
ko.sel3 <- scale(ko.sel2)
ko.sel4 <- sqrt(ko.sel)
env.amp1 <-read.csv("data/random_forest.csv",header=T,row.names=1)
colnames(ko.sel) == env.amp1$name


env.tmp1<-env.amp1[,c("NMDS1_species","NMDS2_species")]
#Fig.S5
#env.tmp1<-env.amp1[,c("plantrichness","NH3g","pH","IN")]


env.tmp2 <- data.frame(t(ko.sel))
spman.d12 = corr.test(env.tmp1, env.tmp2,use="pairwise",method="pearson",adjust="fdr",alpha=.05,ci=FALSE)


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
                   levels = c("NMDS2_species","NMDS1_species"))

cor.out$X<- factor(cor.out$X, 
                   labels = c("NMDS2", "NMDS1"))
#Fig.S5
#cor.out$X<- factor(cor.out$X, 
#levels = c("plantrichness","NH3g","pH","IN"))

#cor.out$X<- factor(cor.out$X, 
#labels = c("Plant richness","NH3 stress","pH","IN"))




merged_df <- merge(cor.out, ko.sel.ID0, by = "Y")
cor.out<-merged_df


library(splitstackshape)

#cor.out<-cSplit(cor.out, "Y", "_")

cor.out$category <- factor(cor.out$Phylum, levels = c("Acidobacteriota","Actinobacteriota","Armatimonadota","Bacteroidota","Chloroflexota", "Gemmatimonadota",
                                                      "Planctomycetota","Proteobacteria","Verrucomicrobiota"))



cor.out$category <- factor(cor.out$Phylum,labels = c("Acidobacteriota","Actinobacteriota","ARM","Bacteroidota","CHL", "GEM",
                                                     "Planctomycetota","Proteobacteria","VER"))




p<- ggplot(cor.out, aes(Y,X)) +
  facet_grid(.~ category,  scales = "free", space = "free" )+
  geom_tile(aes(fill = r), size=0)+
  #geom_hline(yintercept = c(8.5,9.5), color = "yellow")+
  #scale_fill_gradient(guide = "legend", high='green', low='blue',name="rho")+
  scale_fill_gradient2(guide = "legend", high='red3',mid="white", low='dodgerblue3',name="rho")+
  #scale_fill_gradientn(colours = c('dodgerblue3',"white","red3"),name="rho")+
  theme(strip.text = element_text(size = 12,face="bold", color = c("black")),
        strip.background = element_rect(fill = "red"),
        panel.spacing = unit(0.2, "lines"),
        axis.title= element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y=element_text(colour="black", size=12, face="bold"))
p


g <- ggplot_gtable(ggplot_build(p))
strips <- which(grepl('strip-', g$layout$name))
fills <- c("lightskyblue","lightgoldenrod2",'sienna1',"violet","mediumpurple1","seagreen1","lightpink2", "paleturquoise2", "wheat2")
for (i in seq_along(strips)) {
  k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  l <- which(grepl('titleGrob', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- fills[i]
}
plot(g)

ggsave("Fig.S2d.pdf", g,height=2,width =31,limitsize = FALSE )


####Fig.S3b
library(EnhancedVolcano)
library(patchwork)
gene5<-read.csv("data/LOGFC.16S.amplicon.OTU.ANAU50.csv",row.names=1,head=T)
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
                      subtitle = "16S OTU: Total = 9183 variables",
                      caption = bquote("AN50 up = 98 variables; AU50 up = 92 variables"),
                      pCutoff = 0.05,
                      FCcutoff = 0.5,
                      pointSize = 2.0,
                      labSize = 6,
                      col=c('gray70', 'gray70', 'gray70', 'red3'),
                      colAlpha = 0.5,
                      ylim = c(0, 20),
                      xlim = c(-5,8),
                      legendPosition = 'right',
                      selectLab = row.names(gene5)[which(names(keyvals5) %in% c("Up","Down","None"))],colCustom = keyvals5)

p5
ggsave("Fig.S3b.pdf", p5,height=6,width =7,limitsize = FALSE )


####Fig.S4
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

dat<-read.csv("data/relative abundance.16S.phylum.amplicon.csv",head=T,row.names=1)

p2 <-ggplot()+
  geom_jitter(data=dat, aes(y=Abundance, x= log(level3+1)),color=dat$samples2,size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=Abundance, x= log(level3+1),color=type ),method="loess", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  
  facet_wrap(~group,ncol=5,scales = "free")+
  theme_bw()+
  xlab(bquote('Ln(N addition rate ('*g~N~m ^-2~yr^-1*'))')) +
  ylab(bquote('Relative abundance of dominant phyla'))+
  scale_color_manual(values = c("blue3","red3", "gold", "#00ff00"))+ 
  scale_y_continuous(labels = scales::percent)+
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18),
        strip.text = element_text(face = "bold", colour = "black", size = 16))+
  labs(color="Nitrogen type")

p2
ggsave("Fig.S4.pdf", p2, width = 28, height = 11)



####Fig.S6
#Fig.S6a
library(vegan)
otu <- read.csv("data/16S.shotgun.otutable.csv",header=T,row.names = 1)

bac.nmds<-metaMDS(vegdist(otu[1:96,3:12156],method = 'bray'),trymax=1000)
bac.nmds

#####permanova
permanova<-adonis(otu[1:96,3:12156] ~ otu$treat)###perm anova### 
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
#write.csv(df_points,"shotgun.bacterial_taxa.NMDS.csv")

####读入新的旋转后的坐标
df_points2 <- read.csv("shotgun.bacterial_taxa.NMDS.csv",header = T,row.names = 1)
df <-df_points2[,c("MDS1", "MDS2", "Type")]
find_hull <- function(df_points2) df_points2[chull(df_points2$MDS1, df_points2$MDS2), ]

hulls <- ddply(df, "Type", find_hull )###border


#plot
p6<-ggplot(data=df_points2, aes(x= MDS1 , y= MDS2,group=Type))+ 
  geom_point(aes(colour=samples),size=6,alpha=0.7)+ 
  geom_polygon(data=hulls, alpha=.1,aes(fill = Type))+
  labs(x = "NMDS1", y = "NMDS2")+
  scale_colour_manual(values=c("#C7E2FF","lightblue1","#95E9FF","#5EB5FF","#339CFF","blue",
                               "#FFE1FF","#FFC8FF","#FFAFFF","#FF63FF","#FF6371","red"))+
  scale_fill_manual(values=c("blue","red"))+
  ggtitle(expression('Stress= 0.059'~~~~'Treatment:'~italic(R)^2~'= 0.398***'))+
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

p6
ggsave("Fig.S6a.pdf", p6,height=5,width =6,limitsize = FALSE )



#Fig.S6b
#4.2.3
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


data<-read.table("data/16S.shotgun.AN50AU50.phylum.txt",header=TRUE,sep="\t")
combined =data##用于计算显著性

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
  "Verrucomicrobiota",
  "Planctomycetota",
  "Armatimonadota",
  
  "Acidobacteriota",
  "Proteobacteria",
  
  "Gemmatimonadota",
  "Bacteroidota",
  "Chloroflexota",
  "Actinobacteriota"
  
))


samples_new <- sapply(as.character(dat$variable), function(x) strsplit(x, "[_]"))

Timepoint <- {}
for(i in 1:18){
  T <- samples_new[[i]][1]
  Timepoint <- append(Timepoint,T)
}

SampleType <- {}
for(i in 1:18){
  T <- samples_new[[i]][2]
  SampleType <- append(SampleType,T)
}

dat$SampleType = SampleType
dat$Timepoint = Timepoint


combined<-read.table("data/16S.shotgun.AN50AU50.phylum.combined.txt",header=TRUE,sep="\t")

combined$AN50_UP_pvalue <- phyper(combined$AN50_up-1,combined$total,sum(combined$total)-combined$total,sum(combined$AN50_up),
                                  lower.tail = FALSE, log.p = FALSE)


combined$AU50_UP_pvalue <- phyper(combined$AU50_up-1,combined$total,sum(combined$total)-combined$total,sum(combined$AU50_up),
                                  lower.tail = FALSE, log.p = FALSE)


combined[which(combined$AN50_UP_pvalue < 0.05),'sig.AN'] <- 'Significant'
combined[which(combined$AU50_UP_pvalue < 0.05),'sig.AU'] <- 'Significant'

dat[which(dat$value <= 1 ),'sig'] <- 'Non-significant'
dat[which(dat$value > 1),'sig'] <- 'Significant'
dat1 <- subset(dat,Timepoint=="AN50")
dat2 <- subset(dat,Timepoint=="AU50")


ggplot(dat, aes(x=value, y=Annotation,colour=Timepoint,shape=sig)) +
  scale_colour_manual(name="",values = c("royalblue1","#F24F50","#9B4F0F","#6FB98F"))+
  scale_shape_manual(name="",values=c(1,1,19,19))+facet_wrap(~SampleType,scales="free_x",ncol=2)+
  scale_shape_manual(values  = c(21,21), limits = c('Non-significant', 'Significant')) +
  geom_point(size =6,stroke = 1)+theme_bw()+
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


#Fig.S6cd
#R4.0.4
library(rfPermute)
library(randomForest)
library(rfUtilities)
fun_env<-read.csv("data/random_forest.csv",head=T,row.names=1)

#NMDS1
fun.rf <- randomForest(fun_env$NMDS1_species_shotgun~IN+pH+plantrichness+NH3g,
                       data=fun_env,num.cores = 4, nrep = 1000)
rf.significance(fun.rf,fun_env[,c(7:10)])
fun.rf##the details of model

fun.rfP1 <- rfPermute( fun_env$NMDS1_species_shotgun~IN+pH+plantrichness+NH3g,
                       data=fun_env,num.cores = 4, nrep = 1000)
importance(fun.rfP1)
plotImportance(fun.rfP1,imp.type = NULL,plot = TRUE,alpha = 0.05,ranks = TRUE,ylab = "Increase in MSE (%)",scale = TRUE,
               sig.only = FALSE)


#NMDS2
fun.rf2 <- randomForest(fun_env$NMDS2_species_shotgun~IN+pH+plantrichness+NH3g,
                        data=fun_env,num.cores = 4, nrep = 1000)##ntree=500 mtry = 6 default
rf.significance(fun.rf2,fun_env[,c(7:10)])##change the combined table with different columns
fun.rf2##the details of model

fun.rfP12 <- rfPermute( fun_env$NMDS2_species_shotgun~IN+pH+plantrichness+NH3g,
                        data=fun_env,num.cores = 4, nrep = 1000)
importance(fun.rfP12)
plotImportance(fun.rfP12,imp.type = NULL,plot = TRUE,alpha = 0.05,ranks = TRUE,ylab = "Increase in MSE (%)",scale = TRUE,
               sig.only = FALSE)



#Fig.S6e & Fig.S8
library(readxl)
library(psych)
library(tidyr)
library(ggplot2)
library(vegan)
library(colorRamps)
library(ggrepel)
library(grid)

ko.sel0 <-read.csv("data/16S.shotgun.otutable_NMDS.csv",header=T,row.names=1)
ko.sel.ID0 <-  ko.sel0[, c(1:3)]
ko.sel <-ko.sel0[, c(-1:-3)]
ko.sel2 <- log(ko.sel,10)
ko.sel3 <- scale(ko.sel)
ko.sel4 <- scale(ko.sel3)
ko.sel5 <- sqrt(ko.sel)

env.amp1 <-read.csv("data/random_forest.csv",header=T,row.names=1)
colnames(ko.sel) == env.amp1$name

env.tmp1<-env.amp1[,c("NMDS1_species_shotgun","NMDS2_species_shotgun")]
#Fig.S8
#env.tmp1<-env.amp1[,c("plantrichness","NH3g","pH","IN")]

env.tmp2 <- data.frame(t(ko.sel))
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
                   levels = c("NMDS2_species_shotgun","NMDS1_species_shotgun"))

cor.out$X<- factor(cor.out$X, 
                   labels = c("NMDS2", "NMDS1"))

#Fig.S8
#cor.out$X<- factor(cor.out$X, 
# levels = c("plantrichness","NH3g","pH","IN"))

#cor.out$X<- factor(cor.out$X, 
#labels = c("Plant richness","NH3 stress","pH","IN"))

merged_df <- merge(cor.out, ko.sel.ID0, by = "Y")
cor.out<-merged_df

library(splitstackshape)
cor.out$category <- factor(cor.out$Phylum, levels = c("Acidobacteriota","Actinobacteriota","Armatimonadota","Bacteroidota","Chloroflexota", "Gemmatimonadota", "Planctomycetota","Proteobacteria","Verrucomicrobiota")) 
  
                                                      

cor.out$category <- factor(cor.out$Phylum,labels = c("Acidobacteriota","Actinobacteriota","A","Bacteroidota","CHL", "GEM",
                                                     "PLA","Proteobacteria","VER"))

p8<- ggplot(cor.out, aes(Y,X)) +
  facet_grid(.~ category,  scales = "free", space = "free" )+
  geom_tile(aes(fill = r), size=0)+
  scale_fill_gradient2(guide = "legend", high='red3',mid="white", low='dodgerblue3',name="rho")+
  theme(strip.text = element_text(size = 12,face="bold", color = c("black")),
        strip.background = element_rect(fill = "red"),
        panel.spacing = unit(0.2, "lines"),
        axis.title= element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y=element_text(colour="black", size=12, face="bold"))
p8


g <- ggplot_gtable(ggplot_build(p8))
strips <- which(grepl('strip-', g$layout$name))
fills <- c("lightskyblue","lightgoldenrod2",'sienna1',"violet","mediumpurple1","seagreen1","lightpink2", "paleturquoise2", "wheat2")
for (i in seq_along(strips)) {
  k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  l <- which(grepl('titleGrob', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- fills[i]
}
plot(g)

ggsave("Fig.S6e.pdf", g,height=2,width =31,limitsize = FALSE )


####Fig.S7
#Fig.S7b
library(EnhancedVolcano)
library(patchwork)

gene5<-read.csv("data/LOGFC.16SrRNA.AN50AU50.csv",row.names=1,head=T)
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
                      subtitle = "16S Shotgun: Total = 12154 variables",
                      caption = bquote("AN50 up = 941 variables; AU50 up = 1012 variables"),
                      pCutoff = 0.05,
                      FCcutoff = 0.5,
                      pointSize = 2.0,
                      labSize = 6,
                      col=c('gray70', 'gray70', 'gray70', 'red3'),
                      colAlpha = 0.5,
                      ylim = c(0, 20),
                      xlim = c(-7,8),
                      legendPosition = 'right',
                      selectLab = row.names(gene5)[which(names(keyvals5) %in% c("Up","Down","None"))],colCustom = keyvals5)
p5
ggsave("Fig.S7b.pdf", p5,height=6,width =7,limitsize = FALSE )


#Fig.S7c
dat<-read.csv("data/relative abundance.16SrRNA.phylum.shotgun.csv",head=T,row.names=1)

p6 <-ggplot()+
  geom_jitter(data=dat, aes(y=Abundance, x= log(level3+1)),color=dat$samples2,size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=Abundance, x= log(level3+1),color=type ),method="loess", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  facet_wrap(~group,ncol=5,scales = "free")+
  theme_bw()+
  xlab(bquote('Ln(N addition rate ('*g~N~m ^-2~yr^-1*'))')) +
  ylab(bquote('Relative abundance of dominant phyla'))+
  scale_color_manual(values = c("blue3","red3", "gold", "#00ff00"))+ 
  scale_y_continuous(labels = scales::percent)+
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18),
        strip.text = element_text(face = "bold", colour = "black", size = 16))+
  labs(color="Nitrogen type")

p6
ggsave("Fig.S7c.pdf", p6, width = 28, height = 11)




####Fig.S25-S30
#R3.6.2
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

dat<-read.csv("data/bacteria taxa_contribution_amplicon.csv",head=T,row.names=1)
dat.AN <- dat[1:48,]
dat.AU <- dat[49:96,]
#Fig.S25
#MCOA1
p1 <-ggplot()+
  geom_jitter(data=dat, aes(y=Acidobacteria, x= SynVar1,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=Acidobacteria, x= SynVar1),color="black",method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA1')) +
  ylab(bquote('Acidobacteriota (Reads)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.420***'))

p1
summary(lm(dat$Acidobacteria~dat$SynVar1))


p2 <-ggplot()+
  geom_jitter(data=dat, aes(y=Actinobacteria, x= SynVar1,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=Actinobacteria, x= SynVar1),color="black",method="lm",linewidth=2,linetype=1)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA1')) +
  ylab(bquote('Actinobacteria (Reads)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.052**'))

p2
summary(lm(dat$Actinobacteria~dat$SynVar1))


p3 <-ggplot()+
  geom_jitter(data=dat, aes(y=Bacteroidetes, x= SynVar1,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=Bacteroidetes, x= SynVar1),color="black",method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA1')) +
  ylab(bquote('Bacteroidota (Reads)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.287***'))

p3
summary(lm(dat$Bacteroidetes~dat$SynVar1))


p4 <-ggplot()+
  geom_jitter(data=dat, aes(y=Chloroflexi, x= SynVar1,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=Chloroflexi, x= SynVar1),color="black",method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA1')) +
  ylab(bquote('Chloroflexota (Reads)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.434***'))

p4
summary(lm(dat$Chloroflexi~dat$SynVar1))


p5 <-ggplot()+
  geom_jitter(data=dat, aes(y=Armatimonadetes, x= SynVar1,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=Armatimonadetes, x= SynVar1),color="black",method="lm", span=0.9,linewidth=2,linetype=1)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA1')) +
  ylab(bquote('Armatimonadota (Reads)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.069**'))

p5
summary(lm(dat$Armatimonadetes~dat$SynVar1))

p6 <-ggplot()+
  geom_jitter(data=dat, aes(y=Gemmatimonadetes, x= SynVar1,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=Gemmatimonadetes, x= SynVar1),color="black",method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA1')) +
  ylab(bquote('Gemmatimonadota (Reads)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.181***'))

p6
summary(lm(dat$Gemmatimonadetes~dat$SynVar1))


p7 <-ggplot()+
  geom_jitter(data=dat, aes(y=Planctomycetes, x= SynVar1,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=Planctomycetes, x= SynVar1),color="black",method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA1')) +
  ylab(bquote('Planctomycetota (Reads)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.129***'))

p7
summary(lm(dat$Planctomycetes~dat$SynVar1))



p8 <-ggplot()+
  geom_jitter(data=dat, aes(y=Verrucomicrobia, x= SynVar1,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=Verrucomicrobia, x= SynVar1),color="black",method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA1')) +
  ylab(bquote('Verrucomicrobiota (Reads)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.216***'))

p8
summary(lm(dat$Verrucomicrobia~dat$SynVar1))


p <- 
  p1 + p2 + p3 + p4+p5+p6+p7+p8+
  plot_layout(guides = "collect",ncol=4,nrow=2) +
  theme(legend.position='right')
p

ggsave("Fig.S25.pdf", p, height=12,width =25,limitsize = FALSE )


####Fig.S26
#MCOA2
dat<-read.csv("data/bacteria taxa_contribution_amplicon.csv",head=T,row.names=1)
dat.AN <- dat[1:48,]
dat.AU <- dat[49:96,]

p1 <-ggplot()+
  geom_jitter(data=dat, aes(y=Acidobacteria, x= SynVar2,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=Acidobacteria, x= SynVar2),color="black",method="lm", span=0.9,linewidth=2,linetype=1)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA2')) +
  ylab(bquote('Acidobacteriota (Reads)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.183***'))

p1
summary(lm(dat$Acidobacteria~dat$SynVar2))


p2 <-ggplot()+
  geom_jitter(data=dat, aes(y=Actinobacteria, x= SynVar2,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=Actinobacteria, x= SynVar2),color="black",method="lm",linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA2')) +
  ylab(bquote('Actinobacteria (Reads)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.129***'))

p2
summary(lm(dat$Actinobacteria~dat$SynVar2))


p3 <-ggplot()+
  geom_jitter(data=dat, aes(y=Bacteroidetes, x= SynVar2,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=Bacteroidetes, x= SynVar2),color="black",method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA2')) +
  ylab(bquote('Bacteroidota (Reads)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.426***'))

p3
summary(lm(dat$Bacteroidetes~dat$SynVar2))



p4 <-ggplot()+
  geom_jitter(data=dat, aes(y=Armatimonadetes, x= SynVar2,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=Armatimonadetes, x= SynVar2),color="black",method="lm", span=0.9,linewidth=2,linetype=1)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA2')) +
  ylab(bquote('Armatimonadota (Reads)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.174***'))

p4
summary(lm(dat$Armatimonadetes~dat$SynVar2))




p5 <-ggplot()+
  geom_jitter(data=dat, aes(y=Gemmatimonadetes, x= SynVar2,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=Gemmatimonadetes, x= SynVar2),color="black",method="lm", span=0.9,linewidth=2,linetype=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA2')) +
  ylab(bquote('Gemmatimonadota (Reads)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.030'^ns))

p5
summary(lm(dat$Gemmatimonadetes~dat$SynVar2))


p6 <-ggplot()+
  geom_jitter(data=dat, aes(y=Planctomycetes, x= SynVar2,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=Planctomycetes, x= SynVar2),color="black",method="lm", span=0.9,linewidth=2,se=F,linetype=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA2')) +
  ylab(bquote('Planctomycetota (Reads)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= -0.002'^ns))

p6
summary(lm(dat$Planctomycetes~dat$SynVar2))


p7 <-ggplot()+
  geom_jitter(data=dat, aes(y=Proteobacteria, x= SynVar2,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=Proteobacteria, x= SynVar2),color="black",method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA2')) +
  ylab(bquote('Proteobacteria (Reads)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.429***'))

p7
summary(lm(dat$Proteobacteria~dat$SynVar2))

p8<-ggplot()+
  geom_jitter(data=dat, aes(y=Verrucomicrobia, x= SynVar2,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=Verrucomicrobia, x= SynVar2),color="black",method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA2')) +
  ylab(bquote('Verrucomicrobia (Reads)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.429***'))

p8
summary(lm(dat$Verrucomicrobia~dat$SynVar2))


p <- 
  p1 + p2 + p3 + p4 + p5+p6+p7+
  plot_layout(guides = "collect",ncol=4,nrow=2) +
  theme(legend.position='right')
p

ggsave("Fig.S26.pdf", p,height=12,width =25,limitsize = FALSE )



##Fig.S27
#MCOA1
dat<-read.csv("data/1functional potential contribution_amplicon.csv",head=T,row.names=1)
dat.AN <- dat[1:48,]
dat.AU <- dat[49:96,]
fit1 <- lm(SynVar1~YAS5+YAS6+YAS7+YAS9+YAS11+YAS12+YAS14+YAS16+YAS17+YAS20+YAS22+YAS23+YAS24+YAS25+YAS26+YAS27,data = dat)
p1 <-ggplot()+
  geom_jitter(data=dat, aes(y=YAS5, x= SynVar1,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=YAS5, x= SynVar1),color="black",method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA1')) +
  ylab(bquote('Carbohydrate metabolism (Reads)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.339***'))

p1
summary(lm(dat$YAS5~dat$SynVar1))


p2 <-ggplot()+
  geom_jitter(data=dat, aes(y=YAS6, x= SynVar1,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=YAS6, x= SynVar1),color="black",method="lm",linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA1')) +
  ylab(bquote('Carbon fixation pathways in prokaryotes (Reads)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.044*'))

p2
summary(lm(dat$YAS6~dat$SynVar1))


p3 <-ggplot()+
  geom_jitter(data=dat, aes(y=YAS7, x= SynVar1,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=YAS7, x= SynVar1),color="black",method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA1')) +
  ylab(bquote('Cell growth and death (Reads)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.573***'))

p3
summary(lm(dat$YAS7~dat$SynVar1))


p4 <-ggplot()+
  geom_jitter(data=dat, aes(y=YAS9, x= SynVar1,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=YAS9, x= SynVar1),color="black",method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA1')) +
  ylab(bquote('Energy metabolism (Reads)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.138***'))

p4
summary(lm(dat$YAS9~dat$SynVar1))


p5 <-ggplot()+
  geom_jitter(data=dat, aes(y=YAS11, x= SynVar1,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=YAS11, x= SynVar1),color="black",method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA1')) +
  ylab(bquote('Fatty acid biosynthesis (Reads)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.136***'))

p5
summary(lm(dat$YAS11~dat$SynVar1))


p6 <-ggplot()+
  geom_jitter(data=dat, aes(y=YAS12, x= SynVar1,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=YAS12, x= SynVar1),color="black",method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA1')) +
  ylab(bquote('Folding, sorting and degradation (Reads)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.842***'))

p6
summary(lm(dat$YAS12~dat$SynVar1))


p7 <-ggplot()+
  geom_jitter(data=dat, aes(y=YAS14, x= SynVar1,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=YAS14, x= SynVar1),color="black",method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA1')) +
  ylab(bquote('Lipid metabolism (Reads)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.744***'))

p7
summary(lm(dat$YAS14~dat$SynVar1))


p8 <-ggplot()+
  geom_jitter(data=dat, aes(y=YAS16, x= SynVar1,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=YAS16, x= SynVar1),color="black",method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA1')) +
  ylab(bquote('Membrane transport (Reads)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.506***'))

p8
summary(lm(dat$YAS16~dat$SynVar1))


p9 <-ggplot()+
  geom_jitter(data=dat, aes(y=YAS17, x= SynVar1,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=YAS17, x= SynVar1),color="black",method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA1')) +
  ylab(bquote('Metabolism of cofactors and vitamins (Reads)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.750***'))

p9
summary(lm(dat$YAS17~dat$SynVar1))


p11 <-ggplot()+
  geom_jitter(data=dat, aes(y=YAS20, x= SynVar1,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=YAS20, x= SynVar1),color="black",method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA1')) +
  ylab(bquote('Nitrogen metabolism (reads)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.586***'))

p11
summary(lm(dat$YAS20~dat$SynVar1))


p12 <-ggplot()+
  geom_jitter(data=dat, aes(y=YAS22, x= SynVar1,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=YAS22, x= SynVar1),color="black",method="lm",linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA1')) +
  ylab(bquote('Quorum sensing (Reads)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.126***'))

p12
summary(lm(dat$YAS22~dat$SynVar1))


p13 <-ggplot()+
  geom_jitter(data=dat, aes(y=YAS23, x= SynVar1,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=YAS23, x= SynVar1),color="black",method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA1')) +
  ylab(bquote('Replication and repair (Reads)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.871***'))

p13
summary(lm(dat$YAS23~dat$SynVar1))


p14 <-ggplot()+
  geom_jitter(data=dat, aes(y=YAS24, x= SynVar1,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=YAS24, x= SynVar1),color="black",method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA1')) +
  ylab(bquote('Signal transduction (Reads)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.480***'))

p14
summary(lm(dat$YAS24~dat$SynVar1))


p15 <-ggplot()+
  geom_jitter(data=dat, aes(y=YAS25, x= SynVar1,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=YAS25, x= SynVar1),color="black",method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA1')) +
  ylab(bquote('Sulfur metabolism (Reads)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.755***'))

p15
summary(lm(dat$YAS25~dat$SynVar1))


p16 <-ggplot()+
  geom_jitter(data=dat, aes(y=YAS26, x= SynVar1,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=YAS26, x= SynVar1),color="black",method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA1')) +
  ylab(bquote('Transcription (Reads)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.409***'))

p16
summary(lm(dat$YAS26~dat$SynVar1))


p17 <-ggplot()+
  geom_jitter(data=dat, aes(y=YAS27, x= SynVar1,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=YAS27, x= SynVar1),color="black",method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA1')) +
  ylab(bquote('Translation (Reads)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.850***'))

p17
summary(lm(dat$YAS27~dat$SynVar1))



p <- p1 + p2+ p3 + p4 +
  p5+p7+p6+p8+p9+
  p11 + p12 + p13 + p14 + p15+p16+p17+
  plot_layout(guides = "collect",ncol=4,nrow=4) +
  theme(legend.position='right')

p

ggsave("Fig.S27.pdf", p,height=30,width =30,limitsize = FALSE )

#Fig.S28
####MCOA2
dat<-read.csv("data/2functional potential contribution_amplicon.csv",head=T,row.names=1)
dat.AN <- dat[1:48,]
dat.AU <- dat[49:96,]
fit1 <- lm(SynVar2~YAS2+YAS3+YAS4+YAS6+YAS8+YAS9+YAS10+YAS12+YAS13+YAS15+YAS16+YAS17+YAS20+YAS22+YAS24,data = dat)


p1 <-ggplot()+
  geom_jitter(data=dat, aes(y=YAS2, x= SynVar2,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=YAS2, x= SynVar2),color="black",method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA2')) +
  ylab(bquote('Bacterial chemotaxis (Reads)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.335***'))

p1
summary(lm(dat$YAS2~dat$SynVar2))


p2 <-ggplot()+
  geom_jitter(data=dat, aes(y=YAS3, x= SynVar2,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=YAS3, x= SynVar2),color="black",method="lm",linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA2')) +
  ylab(bquote('Biofilm formation (Reads)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.112***'))

p2
summary(lm(dat$YAS3~dat$SynVar2))


p3 <-ggplot()+
  geom_jitter(data=dat, aes(y=YAS4, x= SynVar2,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=YAS4, x= SynVar2),color="black",method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA2')) +
  ylab(bquote('Biosynthesis of other secondary metabolites (Reads)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.453***'))

p3
summary(lm(dat$YAS4~dat$SynVar2))


p4 <-ggplot()+
  geom_jitter(data=dat, aes(y=YAS6, x= SynVar2,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=YAS6, x= SynVar2),color="black",method="lm", span=0.9,linetype=1,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA2')) +
  ylab(bquote('Carbon fixation pathways in prokaryotes (Reads)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.290***'))

p4
summary(lm(dat$YAS6~dat$SynVar2))


p5 <-ggplot()+
  geom_jitter(data=dat, aes(y=YAS8, x= SynVar2,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=YAS8, x= SynVar2),color="black",method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA2')) +
  ylab(bquote('Cell motility (Reads)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.487***'))

p5
summary(lm(dat$YAS8~dat$SynVar2))


p6 <-ggplot()+
  geom_jitter(data=dat, aes(y=YAS9, x= SynVar2,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=YAS9, x= SynVar2),color="black",method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA2')) +
  ylab(bquote('Energy metabolism (Reads)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.269***'))

p6
summary(lm(dat$YAS9~dat$SynVar2))


p7 <-ggplot()+
  geom_jitter(data=dat, aes(y=YAS10, x= SynVar2,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=YAS10, x= SynVar2),color="black",method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA2')) +
  ylab(bquote('Exopolysaccharide biosynthesis (Reads)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.495***'))

p7
summary(lm(dat$YAS10~dat$SynVar2))


p8 <-ggplot()+
  geom_jitter(data=dat, aes(y=YAS12, x= SynVar2,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=YAS12, x= SynVar2),color="black",method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA2')) +
  ylab(bquote('Folding, sorting and degradation (Reads)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.044*'))

p8
summary(lm(dat$YAS12~dat$SynVar2))


p9 <-ggplot()+
  geom_jitter(data=dat, aes(y=YAS13, x= SynVar2,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=YAS13, x= SynVar2),color="black",method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA2')) +
  ylab(bquote('Glycan biosynthesis and metabolism (Reads)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.117***'))

p9
summary(lm(dat$YAS13~dat$SynVar2))


p10 <-ggplot()+
  geom_jitter(data=dat, aes(y=YAS15, x= SynVar2,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=YAS15, x= SynVar2),color="black",method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA2')) +
  ylab(bquote('Lipopolysaccharide biosynthesis (Reads)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.378***'))

p10
summary(lm(dat$YAS15~dat$SynVar2))


p11 <-ggplot()+
  geom_jitter(data=dat, aes(y=YAS16, x= SynVar2,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=YAS16, x= SynVar2),color="black",method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA2')) +
  ylab(bquote('Membrane transport (Reads)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.123***'))

p11
summary(lm(dat$YAS16~dat$SynVar2))


p12 <-ggplot()+
  geom_jitter(data=dat, aes(y=YAS17, x= SynVar2,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=YAS17, x= SynVar2),color="black",method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA2')) +
  ylab(bquote('Metabolism of cofactors and vitamins (reads)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.031*'))

p12
summary(lm(dat$YAS17~dat$SynVar2))


p13 <-ggplot()+
  geom_jitter(data=dat, aes(y=YAS20, x= SynVar2,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=YAS20, x= SynVar2),color="black",method="lm",linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA2')) +
  ylab(bquote('Nitrogen metabolism (Reads)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.046*'))

p13
summary(lm(dat$YAS20~dat$SynVar2))


p14 <-ggplot()+
  geom_jitter(data=dat, aes(y=YAS22, x= SynVar2,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=YAS22, x= SynVar2),color="black",method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA2')) +
  ylab(bquote('Quorum sensing (Reads)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.556***'))

p14
summary(lm(dat$YAS22~dat$SynVar2))


p15 <-ggplot()+
  geom_jitter(data=dat, aes(y=YAS24, x= SynVar2,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=YAS24, x= SynVar2),color="black",method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA2')) +
  ylab(bquote('Signal transduction (Reads)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.061**'))

p15
summary(lm(dat$YAS24~dat$SynVar2))



p <- 
  p1 + p2 + p3 + p4 +p7+p6+p5+p8+p9+p10+p11 + p12 + p13 + p14 +p15+
  plot_layout(guides = "collect",ncol=4,nrow=4) +
  theme(legend.position='right')
p

ggsave("Fig.S28.pdf", p,height=30,width =30,limitsize = FALSE )




#Fig.S29
#MCOA1
dat<-read.csv("data/1microbial traits contribution_amplicon.csv",head=T,row.names=1)
dat.AN <- dat[1:48,]
dat.AU <- dat[49:96,]

p1 <-ggplot()+
  geom_jitter(data=dat, aes(y=bac_shannon, x= SynVar1,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=bac_shannon, x= SynVar1),color="black",method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA1')) +
  ylab(bquote('S.Bacteria'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.774***'))

p1
summary(lm(dat$bac_shannon~dat$SynVar1))


p2 <-ggplot()+
  geom_jitter(data=dat, aes(y=KO_shannon, x= SynVar1,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=KO_shannon, x= SynVar1),color="black",method="lm",linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA1')) +
  ylab(bquote('S.KO'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.733***'))

p2
summary(lm(dat$KO_shannon~dat$SynVar1))


p3 <-ggplot()+
  geom_jitter(data=dat, aes(y=vir_shannon_mgm, x= SynVar1,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=vir_shannon_mgm, x= SynVar1),color="black",method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA1')) +
  ylab(bquote('S.Phage'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.376***'))

p3
summary(lm(dat$vir_shannon_mgm~dat$SynVar1))


p4 <-ggplot()+
  geom_jitter(data=dat, aes(y=fun_shannon, x= SynVar1,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=fun_shannon, x= SynVar1),color="black",method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA1')) +
  ylab(bquote('S.Fungi'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.313***'))

p4
summary(lm(dat$fun_shannon~dat$SynVar1))


p5 <-ggplot()+
  geom_jitter(data=dat, aes(y=HGT_per_assembled_genes, x= SynVar1,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=HGT_per_assembled_genes, x= SynVar1),color="black",method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA1')) +
  ylab(bquote('HGT Events \n (HGT Events per 1K Assembled Genes)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.263***'))

p5
summary(lm(dat$HGT_per_assembled_genes~dat$SynVar1))


p6 <-ggplot()+
  geom_jitter(data=dat, aes(y=-1*growth_rate, x= SynVar1,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=-1*growth_rate, x= SynVar1),color="black",method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA1')) +
  ylab(bquote('Growth Rate \n (-1 * Minimal Doubling Time (hours))'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.646***'))

p6
summary(lm(dat$growth_rate~dat$SynVar1))


p7 <-ggplot()+
  geom_jitter(data=dat, aes(y=Ags, x= SynVar1,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=Ags, x= SynVar1),color="black",method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA1')) +
  ylab(bquote('Average Genome Size (Mbp)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.201***'))

p7
summary(lm(dat$Ags~dat$SynVar1))


p8 <-ggplot()+
  geom_jitter(data=dat, aes(y=rrn, x= SynVar1,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=rrn, x= SynVar1),color="black",method="lm", span=0.9,linewidth=2,linetype=2,se=F)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA1')) +
  ylab(bquote('Average 16S Copy Number'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= -0.001'^ns))

p8
summary(lm(dat$rrn~dat$SynVar1))


p9 <-ggplot()+
  geom_jitter(data=dat, aes(y=GC, x= SynVar1,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=GC, x= SynVar1),color="black",method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA1')) +
  ylab(bquote('GC Content (%)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.544***'))

p9
summary(lm(dat$GC~dat$SynVar1))


p10 <-ggplot()+
  geom_jitter(data=dat, aes(y=Resfam, x= SynVar1,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=Resfam, x= SynVar1),color="black",method="lm", span=0.9,linewidth=2,se=F,linetype=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA1')) +
  ylab(bquote('S.ARG'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.014'^ns))

p10
summary(lm(dat$Resfam~dat$SynVar1))



p <- 
  p1 + p2 + p3 + p4 +
  p5+p6+p7+p8+p9+p10+
  plot_layout(guides = "collect",ncol=5,nrow=2) +
  theme(legend.position='right')
p

ggsave("Fig.S29.pdf", p, height=14,width =33,limitsize = FALSE )

#Fig.S30
#MCOA2
dat<-read.csv("data/2microbial traits contribution_amplicon.csv",head=T,row.names=1)
dat.AN <- dat[1:48,]
dat.AU <- dat[49:96,]

p1 <-ggplot()+
  geom_jitter(data=dat, aes(y=bac_shannon, x= SynVar2,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=bac_shannon, x= SynVar2),color="black",method="lm", span=0.9,linewidth=2,se=F,linetype=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA2')) +
  ylab(bquote('S.Bacteria'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.013'^ns))

p1
summary(lm(dat$bac_shannon~dat$SynVar2))


p2 <-ggplot()+
  geom_jitter(data=dat, aes(y=vir_shannon_mgm, x= SynVar2,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=vir_shannon_mgm, x= SynVar2),color="black",method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA2')) +
  ylab(bquote('S.Phage'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.494***'))

p2
summary(lm(dat$vir_shannon_mgm~dat$SynVar2))


p3 <-ggplot()+
  geom_jitter(data=dat, aes(y=fun_shannon, x= SynVar2,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=fun_shannon, x= SynVar2),color="black",method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA2')) +
  ylab(bquote('S.Fungi'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.278***'))

p3
summary(lm(dat$fun_shannon~dat$SynVar2))


p4 <-ggplot()+
  geom_jitter(data=dat, aes(y=HGT_per_assembled_genes, x= SynVar2,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=HGT_per_assembled_genes, x= SynVar2),color="black",method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA2')) +
  ylab(bquote('HGT Events \n (HGT Events per 1K Assembled Genes)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.481***'))

p4
summary(lm(dat$HGT_per_assembled_genes~dat$SynVar2))


p5 <-ggplot()+
  geom_jitter(data=dat, aes(y=Ags, x= SynVar2,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=Ags, x= SynVar2),color="black",method="lm", span=0.9,linewidth=2,se=F,linetype=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA2')) +
  ylab(bquote('Average Genome Size (Mbp)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.029'^ns))

p5
summary(lm(dat$Ags~dat$SynVar2))


p6 <-ggplot()+
  geom_jitter(data=dat, aes(y=rrn, x= SynVar2,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=rrn, x= SynVar2),color="black",method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA2')) +
  ylab(bquote('Average 16S Copy Number'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.784***'))

p6
summary(lm(dat$rrn~dat$SynVar2))


p7<-ggplot()+
  geom_jitter(data=dat, aes(y=GC, x= SynVar2,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=GC, x= SynVar2),color="black",method="lm", span=0.9,linewidth=2,linetype=1)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA2')) +
  ylab(bquote('GC Content (%)'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.284***'))

p7
summary(lm(dat$GC~dat$SynVar2))


p8 <-ggplot()+
  geom_jitter(data=dat, aes(y=Resfam, x= SynVar2,color=type),size=10,alpha=0.7)+
  geom_smooth( data=dat,aes(y=Resfam, x= SynVar2),color="black",method="lm", span=0.9,linewidth=2)+ 
  geom_point( )+ 
  theme_bw()+
  xlab(bquote('MCOA2')) +
  ylab(bquote('S.ARG'))+
  scale_color_manual(values = c("blue","red", "gold", "#00ff00"))+ 
  theme(axis.text=element_text(size=50,face="bold"),
        axis.title=element_text(size=24,face="bold"),
        axis.text.x =element_text(color="black",size = 22),
        axis.text.y =element_text(color="black",size = 22),
        legend.key = element_rect(fill = 'transparent'),
        legend.text = element_text(face = "bold", colour = "black", size = 16),
        plot.title=element_text(face = "bold", colour = "black", size = 22),
        legend.title = element_text(face = "bold", colour = "black", size = 18))+
  labs(color="Nitrogen type")+
  ggtitle(expression(~italic(R) ^2~'= 0.599***'))

p8
summary(lm(dat$Resfam~dat$SynVar2))


p <- 
  p1 + p2 + p3 + p4 +
  p5+p6+p7+p8+
  plot_layout(guides = "collect",ncol=4,nrow=2) +
  theme(legend.position='right')
p

ggsave("Fig.S30.pdf", p,height=14,width =28,limitsize = FALSE )

####Fig.S23
#Fig.S23a
#R3.6.2
library(omicade4)
functional.traits<-read.csv("data/functional traits_MCOA.csv",head=T,row.names=1)
functional<-read.csv("data/functional potential_MCOA.csv",head=T,row.names=1)
bacteria<-read.csv("data/bacterial taxa 16SrRNA_MCOA.csv",head=T,row.names=1)
#log10
functional.traits2 <- log(functional.traits,10)
functional2<-log(functional,10)
bacteria2<-log(bacteria,10)
#z-score
functional.traits3 <- scale(functional.traits2)
functional3 <-scale(functional2)
bacteria3 <- scale(bacteria2)
#list
lista <- list(functional.traits3,functional3,bacteria3)
names(lista)<-c("functional.traits","functional","bacteria")

head(lista$functional.traits[1:6]) #Agilent platform.
head(lista$functional[1:6])        #H-GU133 platform
head(lista$bacteria[1:6])    #H-GU133 plus 2.0 platform


#MCOA
mcoin <- mcia(lista, cia.nf = 3, cia.scan = FALSE, nsc = TRUE)
mcoin

summary(mcoin)

#mian result
names(mcoin$mcoa)

mcoin$mcoa$pseudoeig

mcoin$mcoa$RV

mcoin$mcoa

sample_type <- colnames(lista$functional.traits)
sample_type <- sapply(strsplit(sample_type, split="\\."), function(x) x[1])
sample_type

library(vegan)
library(ggplot2)
df_points <- as.data.frame(mcoin$mcoa$SynVar)
df_points$samples <- sample_type
df_points$Type <-c("AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN",
                   "AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN",
                   "AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU",
                   "AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU")
library(plyr)
df <-df_points[,c("SynVar1", "SynVar2", "Type")]
find_hull <- function(df_points) df_points[chull(df_points$SynVar1, df_points$SynVar2), ]
hulls <- ddply(df, "Type", find_hull )#border

#准备构建旋转后的坐标
#write.csv(df_points,"MCOA_amplicon.NMDS_shotgun.csv")


####读入新的旋转后的坐标
df_points2 <- read.csv("MCOA_amplicon.NMDS_shotgun.csv",header = T,row.names = 1)
df <-df_points2[,c("SynVar1", "SynVar2", "Type")]
find_hull <- function(df_points2) df_points2[chull(df_points2$SynVar1, df_points2$SynVar2), ]

hulls <- ddply(df, "Type", find_hull )###border


#读取环境数据
env <- read.csv("data/only_randomforest.csv",head=T,row.names=1)

mcoa_env<-envfit(df_points2[,c(6,7)]~.,data=env,perm=999,choices=c(1,2),display='sites')


mcoa_env_adj<-mcoa_env
mcoa_env_adj$vectors$pvals<-p.adjust(mcoa_env_adj$vectors$pvals,method='bonferroni')
mcoa_env_adj

mcoa_env2<-data.frame(cbind(mcoa_env_adj$vectors$arrows,mcoa_env_adj$vectors$r,mcoa_env_adj$vectors$pvals))
names(mcoa_env2)<-c('SynVar1','SynVar2','r2','p.adj')
#write.csv(mcoa_env2,'mcoa_env_16SrRNA_shotgun.csv',quote=FALSE)

mcoa_env3<- read.csv("mcoa_env_16SrRNA_shotgun.csv",head=T,row.names=1)


####绘图
p1<-ggplot(data=df_points2, aes(x= SynVar1 , y= SynVar2))+ 
  
  geom_point(aes(colour=samples),size=6,alpha=0.7)+ 
  geom_polygon(data=hulls, alpha=.1,aes(fill = Type))+
  labs(x = "MCOA1", y = "MCOA2")+
  scale_colour_manual(values=c("#C7E2FF","lightblue1","#95E9FF","#5EB5FF","#339CFF","blue",
                               "#FFE1FF","#FFC8FF","#FFAFFF","#FF63FF","#FF6371","red"))+
  scale_fill_manual(values=c("blue","red"))+
  geom_segment(data=mcoa_env3,aes(x=0,y=0,xend=SynVar1*1.5,yend=SynVar2*1.5),arrow=arrow(length=unit(0.2,'cm')),size=0.8,color='black')+
  geom_text(data=mcoa_env3,aes(SynVar1*2.5,SynVar2*2.5,label=group),color='black',size=3)+
  # scale_x_continuous(limits = c(-3, 3.7))+
  theme_bw()+
  guides(color=guide_legend(title= "Treatment"))+
  theme(legend.position = 'left',legend.direction = 'vertical',
        legend.title = element_text(colour="black", size=13, face="bold"),
        legend.text = element_text(colour="black", size=12),
        axis.text=element_text(size=15,face="bold",colour="black"),
        axis.title=element_text(size=17,face="bold"))+
  geom_vline(xintercept = 0, color = 'gray34', linetype = 1, linewidth = 1.5) +
  geom_hline(yintercept = 0, color = 'gray34', linetype = 1, linewidth = 1.5) 


p1

ggsave("Fig.S23a.pdf", p1,height=5,width =6,limitsize = FALSE )



#Fig.S23b
eigen_value2<- read.csv("data/variation and dimension_MCOA2.csv",head=T,row.names=1)
p2 <- ggplot(data = eigen_value2, aes(x =dimension, y = variation100))+
  geom_bar(stat = "identity", 
           width = 0.8, colour = "black", linewidth = 0.8,
           fill = "thistle3", alpha = 1) +
  ylim(0, 80) + 
  scale_x_continuous(breaks = seq(0, 20, by = 1),minor_breaks=NULL)+
  labs(y = "% of variation", x = "MCOA dimension")+
  theme_bw()+
  theme(
    axis.title = element_text(size = 17, face = "plain", color = "black"), 
    axis.text = element_text(size = 12, face = "plain", color = "black") )
p2
ggsave("Fig.S23b.pdf", p2,height=5,width =6,limitsize = FALSE )



#Fig.S23c
con2 <- read.csv("data/contribution of each dataset_MCOA_16srRNA2.csv",head=T,row.names=1)
#MCOA1
p1 <- ggplot(con2,aes(Dataset,mcoa1))+
  geom_col(aes(fill=Dataset))+
  scale_fill_manual(values=c("#7552cc","royalblue1","#00c16e","#7552cc"),
                    breaks=c("Functional potential","Microbial traits","Bacterial taxa"),
                    labels=c("Functional potential","Microbial traits",  "Bacterial taxa"))+#自定义颜色
  labs(y = "Contribution to MCOA1", x = "Dateset Category")+
  theme_bw()+
  guides(fill=guide_legend(title=NULL,reverse=TRUE))+
  theme(
    panel.background = element_rect(fill = "gray99" ),
    legend.title = element_text(colour="black", size=13, face="bold"),
    legend.text = element_text(colour="black", size=12),
    axis.text=element_text(size=13,face="bold",colour="black"),
    axis.title=element_text(size=16,face="bold"),
    legend.position="top")+
  coord_flip()
p1


#MCOA2
p2 <- ggplot(con2,aes(Dataset,mcoa2))+
  geom_col(aes(fill=Dataset))+
  scale_fill_manual(values=c("#7552cc","royalblue1","#00c16e","#7552cc"),
                    breaks=c("Bacterial taxa","Functional potential","Microbial traits"),
                    labels=c("Bacterial taxa","Functional potential",  "Microbial traits"))+#自定义颜色
  labs(y = "Contribution to MCOA2", x = "Dateset Category")+
  theme_bw()+
  guides(fill=guide_legend(title=NULL,reverse=TRUE))+
  theme(
    panel.background = element_rect(fill = "gray99" ),
    legend.title = element_text(colour="black", size=13, face="bold"),
    legend.text = element_text(colour="black", size=12),
    axis.text=element_text(size=13,face="bold",colour="black"),
    axis.title=element_text(size=16,face="bold"),
    legend.position="top")+
  coord_flip()
p2



####Fig. S24
#Fig.S24a
library(vegan)
library(rmarkdown)

otu <- read.csv("data/Fungi_NMDS.csv")

bac.nmds<-metaMDS(vegdist(otu[1:96,4:3982],method = 'bray'),trymax=100)
bac.nmds

#####permanova分析
permanova<-adonis(otu[1:96,4:3982] ~ otu$treat)###perm anova### 
permanova


x=bac.nmds$points[,1]
y=bac.nmds$points[,2]

max(x)
min(x)
max(y)
min(y)

library(vegan)
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
hulls <- ddply(df, "Type", find_hull )###寻找边界


#准备构建旋转后的坐标
#write.csv(df_points,"FUngi_amplicon.NMDS.csv")

####读入新的旋转后的坐标
df_points2 <- read.csv("data/Fungi_amplicon.NMDS.csv",header = T,row.names = 1)
df <-df_points2[,c("MDS1", "MDS2", "Type")]
find_hull <- function(df_points2) df_points2[chull(df_points2$MDS1, df_points2$MDS2), ]

hulls <- ddply(df, "Type", find_hull )###border

#######画样品二维空间分布

p1<-ggplot(data=df_points2, aes(x= MDS1 , y= MDS2,group=Type))+ 
  geom_point(aes(colour=samples),size=6,alpha=0.7)+ 
  geom_polygon(data=hulls, alpha=.1,aes(fill = Type))+
  labs(x = "NMDS1", y = "NMDS2")+
  #scale_shape_manual(values=c(0, 6, 3, 1))+
  scale_colour_manual(values=c("#C7E2FF","lightblue1","#95E9FF","#5EB5FF","#339CFF","blue",
                               "#FFE1FF","#FFC8FF","#FFAFFF","#FF63FF","#FF6371","red"))+
  scale_fill_manual(values=c("blue","red"))+
  
  ggtitle(expression('Stress= 0.114'~~~~'Treatment:'~italic(R)^2~'= 0.354***'))+
  
  theme_bw()+
  guides(color=guide_legend(title= "Treatment"))+
  #shape=guide_legend(title= "Compartment"))+
  theme(legend.position = 'left',legend.direction = 'vertical',
        legend.title = element_text(colour="black", size=13, face="bold"),
        legend.text = element_text(colour="black", size=12),
        axis.text=element_text(size=15,face="bold",colour="black"),
        axis.title=element_text(size=17,face="bold"),
        plot.title=element_text(face = "bold", colour = "black", size = 15))+
  geom_vline(xintercept = 0, color = 'gray34', linetype = 1, linewidth = 1.5) +
  geom_hline(yintercept = 0, color = 'gray34', linetype = 1, linewidth = 1.5) 
p1
ggsave("Fig.S24a.pdf", p1,height=5,width =6,limitsize = FALSE )




#Fig.S24b
#4.0.4
library(rfPermute)
library(randomForest)
library(rfUtilities)

fun_env<-read.csv("data/random_forest_fungi.csv",head=T,row.names=1)


#NMDS1
fun.rf <- randomForest(fun_env$NMDS1~IN+pH+plantrichness+NH3g,
                       data=fun_env,num.cores = 4, nrep = 1000)##ntree=500 mtry = 6 default
rf.significance(fun.rf,fun_env[,c(7:10)])##change the combined table with different columns,计算显著性 the same result P=0.001,R2=0.399
fun.rf##the details of model

fun.rfP1 <- rfPermute( fun_env$NMDS1~IN+pH+plantrichness+NH3g,
                       data=fun_env,num.cores = 4, nrep = 1000)
importance(fun.rfP1)
plotImportance(fun.rfP1,imp.type = NULL,plot = TRUE,alpha = 0.05,ranks = TRUE,ylab = "Increase in MSE (%)",scale = TRUE,
               sig.only = FALSE)



###NMDS2
fun.rf2 <- randomForest(fun_env$NMDS2~IN+pH+plantrichness+NH3g,
                        data=fun_env,num.cores = 4, nrep = 1000)##ntree=500 mtry = 6 default
rf.significance(fun.rf2,fun_env[,c(7:10)])##change the combined table with different columns,计算显著性 the same result P=0.001,R2=0.399
fun.rf2##the details of model

fun.rfP12 <- rfPermute( fun_env$NMDS2~IN+pH+plantrichness+NH3g,
                        data=fun_env,num.cores = 4, nrep = 1000)
importance(fun.rfP12)
plotImportance(fun.rfP12,imp.type = NULL,plot = TRUE,alpha = 0.05,ranks = TRUE,ylab = "Increase in MSE (%)",scale = TRUE,
               sig.only = FALSE)


#Fig.S24c
library(omicade4)
functional.traits<-read.csv("data/functional traits33.csv",head=T,row.names=1)
functional<-read.csv("data/function2.csv",head=T,row.names=1)
bacteria<-read.csv("data/bacteria2.csv",head=T,row.names=1)
fungi<-read.csv("data/fungi2.csv",head=T,row.names=1)
#log10
functional.traits2 <- log(functional.traits,10)
functional2<-log(functional,10)
bacteria2<-log(bacteria,10)
fungi2<-log(fungi,10)
#z-score标准化：观测值减去平均值，然后再除以标准差，得到均值为0， 标准差为1的数据，且数据符合正太分布
functional.traits3 <- scale(functional.traits2)
functional3 <-scale(functional2)
bacteria3 <- scale(bacteria2)
fungi3 <- scale(fungi2)
#构建列表
lista <- list(functional.traits3,functional3,bacteria3,fungi3)
names(lista)<-c("functional.traits","functional","bacteria","fungi")

#行为样品名称，列为基因，功能性状或者OTU
head(lista$functional.traits[1:6]) #Agilent platform.
head(lista$functional[1:6])        #H-GU133 platform
head(lista$bacteria[1:6])    #H-GU133 plus 2.0 platform
head(lista$fungi[1:6]) 

#执行 MCoIA，详情 ?mcia
mcoin <- mcia(lista, cia.nf = 3, cia.scan = FALSE, nsc = TRUE)
mcoin

summary(mcoin)



#提取主要结果，例如
names(mcoin$mcoa)

#MCoIA 的伪特征值，代表了各 MCoIA 轴承载的协惯量
mcoin$mcoa$pseudoeig

#RV 值，代表了各数据集间结构的一致性
mcoin$mcoa$RV

mcoin$mcoa
#作图观测，以观测 MCoIA 前两轴的特征为例
#图中对象的颜色，根据其所属的样品类型着色
sample_type <- colnames(lista$functional.traits)
sample_type <- sapply(strsplit(sample_type, split="\\."), function(x) x[1])
sample_type

library(vegan)
library(ggplot2)
df_points <- as.data.frame(mcoin$mcoa$SynVar)
df_points$samples <- sample_type
df_points$Type <-c("AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN",
                   "AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN","AN",
                   "AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU",
                   "AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU","AU")
library(plyr)
df <-df_points[,c("SynVar1", "SynVar2", "Type")]
find_hull <- function(df_points) df_points[chull(df_points$SynVar1, df_points$SynVar2), ]
hulls <- ddply(df, "Type", find_hull )###寻找边界

#准备构建旋转后的坐标
#write.csv(df_points,"MCOA_amplicon add fungi.NMDS.csv")


####读入新的旋转后的坐标
df_points2 <- read.csv("data/MCOA_amplicon add fungi.NMDS.csv",header = T,row.names = 1)
df <-df_points2[,c("SynVar1", "SynVar2", "Type")]
find_hull <- function(df_points2) df_points2[chull(df_points2$SynVar1, df_points2$SynVar2), ]

hulls <- ddply(df, "Type", find_hull )###border



#读取环境数据
env <- read.csv("data/only_randomforest.csv",head=T,row.names=1)

mcoa_env<-envfit(df_points2[,c(6,7)]~.,data=env,perm=999,choices=c(1,2),display='sites')


mcoa_env_adj<-mcoa_env
mcoa_env_adj$vectors$pvals<-p.adjust(mcoa_env_adj$vectors$pvals,method='bonferroni')
mcoa_env_adj

mcoa_env2<-data.frame(cbind(mcoa_env_adj$vectors$arrows,mcoa_env_adj$vectors$r,mcoa_env_adj$vectors$pvals))
names(mcoa_env2)<-c('SynVar1','SynVar2','r2','p.adj')
#write.csv(mcoa_env2,'mcoa_env_ITS add fungi.csv',quote=FALSE)


mcoa_env3<- read.csv("data/mcoa_env_ITS add fungi.csv",head=T,row.names=1)


####绘图
p1<-ggplot(data=df_points2, aes(x= SynVar1 , y= SynVar2))+ 
  
  geom_point(aes(colour=samples),size=6,alpha=0.7)+ 
  geom_polygon(data=hulls, alpha=.1,aes(fill = Type))+
  labs(x = "MCOA1", y = "MCOA2")+
  scale_colour_manual(values=c("#C7E2FF","lightblue1","#95E9FF","#5EB5FF","#339CFF","blue",
                               "#FFE1FF","#FFC8FF","#FFAFFF","#FF63FF","#FF6371","red"))+
  scale_fill_manual(values=c("blue","red"))+
  geom_segment(data=mcoa_env3,aes(x=0,y=0,xend=SynVar1*1.5,yend=SynVar2*1.5),arrow=arrow(length=unit(0.2,'cm')),size=0.8,color='black')+
  geom_text(data=mcoa_env3,aes(SynVar1*2.5,SynVar2*2.5,label=group),color='black',size=3)+
  #scale_x_continuous(limits = c(-3, 3.7))+
  theme_bw()+
  guides(color=guide_legend(title= "Treatment"))+
  theme(legend.position = 'left',legend.direction = 'vertical',
        legend.title = element_text(colour="black", size=13, face="bold"),
        legend.text = element_text(colour="black", size=12),
        axis.text=element_text(size=15,face="bold",colour="black"),
        axis.title=element_text(size=17,face="bold"))+
  geom_vline(xintercept = 0, color = 'gray34', linetype = 1, linewidth = 1.5) +
  geom_hline(yintercept = 0, color = 'gray34', linetype = 1, linewidth = 1.5) 


p1

ggsave("Fig.S24c.pdf", p1,height=5,width =5.6,limitsize = FALSE )


#Fig.S24e
#####画轴的解释量
eigen_value <- as.data.frame(mcoin$mcoa$pseudoeig)
eigen_value$number <-c(1:20)
names(eigen_value)[1:2] <- c('eigenvalues', 'dimension')

#write.csv(eigen_value,"variation and dimension.csv")

eigen_value2<- read.csv("data/variation and dimension_fungi.csv",head=T,row.names=1)
p2 <- ggplot(data = eigen_value2, aes(x =dimension, y = variation100))+
  geom_bar(stat = "identity", 
           width = 0.8, colour = "black", linewidth = 0.8,
           fill = "thistle3", alpha = 1) +
  ylim(0, 80) + # 设置y轴范围
  scale_x_continuous(breaks = seq(0, 20, by = 1),minor_breaks=NULL)+
  labs(y = "% of variation", x = "MCOA dimension")+
  theme_bw()+
  theme(
    axis.title = element_text(size = 17, face = "plain", color = "black"), # 设置标题的字体及大小
    axis.text = element_text(size = 12, face = "plain", color = "black") # 设置坐标轴的字体及大小
  )
p2
ggsave("Fig.S24e.pdf", p2,height=5,width =6,limitsize = FALSE )


#Fig.S24f
con2 <- read.csv("data/contribution_MCOA_fungi.csv",head=T,row.names=1)
#1轴
p1 <- ggplot(con2,aes(Dataset,mcoa1))+
  geom_col(aes(fill=Dataset))+
  scale_fill_manual(values=c("#7552cc","royalblue1","#00c16e","#7552cc"),
                    breaks=c("aFungal.taxa","Bacterial.taxa","Functional.potential","Microbial.traits"),
                    labels=c("Fungal taxa","Bacterial taxa","Functional potential","Microbial traits"))+#自定义颜色
  labs(y = "Contribution to MCOA1", x = "Dateset Category")+
  theme_bw()+
  guides(fill=guide_legend(title=NULL,reverse=TRUE))+
  #scale_fill_discrete(breaks=c("n", "p"),
  #labels=c("Trait negatively associated with MCOA1", "Trait positively associated with MCOA1"))+
  #shape=guide_legend(title= "Compartment"))+
  theme(
    panel.background = element_rect(fill = "gray99" ),
    legend.title = element_text(colour="black", size=13, face="bold"),
    legend.text = element_text(colour="black", size=12),
    axis.text=element_text(size=13,face="bold",colour="black"),
    axis.title=element_text(size=16,face="bold"),
    legend.position="top")+
  
  coord_flip()
p1


#2轴
p2 <- ggplot(con2,aes(Dataset,mcoa2))+
  geom_col(aes(fill=Dataset))+
  scale_fill_manual(values=c("#7552cc","royalblue1","#00c16e","#7552cc"),
                    breaks=c("aFungal.taxa","Bacterial.taxa","Functional.potential","Microbial.traits"),
                    labels=c("Fungal taxa","Bacterial taxa","Functional potential","Microbial traits"))+#自定义颜色
  labs(y = "Contribution to MCOA2", x = "Dateset Category")+
  theme_bw()+
  guides(fill=guide_legend(title=NULL,reverse=TRUE))+
  #scale_fill_discrete(breaks=c("n", "p"),
  #labels=c("Trait negatively associated with MCOA1", "Trait positively associated with MCOA1"))+
  #shape=guide_legend(title= "Compartment"))+
  theme(
    panel.background = element_rect(fill = "gray99" ),
    legend.title = element_text(colour="black", size=13, face="bold"),
    legend.text = element_text(colour="black", size=12),
    axis.text=element_text(size=13,face="bold",colour="black"),
    axis.title=element_text(size=16,face="bold"),
    legend.position="top")+
  
  coord_flip()
p2

