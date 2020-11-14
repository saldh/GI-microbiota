setwd("D:/Data/16S_data/GI-microbiota/result/")
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(plyr)
library(ggrepel)
library(dplyr)
library(agricolae)
### alpha diversity --- qiime2 results
Shannon <- read.table('alpha/shannon.tsv',header = T,sep = '')
Evenness <- read.table('alpha/evenness.tsv',header = T,sep = '')
Richness <- read.table('alpha/observed_otus.tsv',header = T,sep = '')
Faith_pd <- read.table('alpha/faith_pd.tsv',header = T,sep = '')
Chao1 <- read.table("alpha/chao1.tsv",header = T,sep = '')
Chao1 <- Chao1[rownames(Chao1) %in% rownames(Shannon),]
Simpson <- read.table("alpha/simpson.tsv",header = T,sep = '')
Simpson <- Simpson[rownames(Simpson) %in% rownames(Shannon),]
alpha <- cbind(Shannon,Evenness,Richness,Faith_pd,Chao1,Simpson)
colnames(alpha) <- c("Shannon","Evenness","Richness","Faith_pd","Chao1","Simpson")
# 合并metadata
alpha$SampleID <- rownames(alpha)
metadata <- read.table("metadata.txt",header = T,sep = "\t")
submeta <- metadata[metadata$SampleID %in% alpha$SampleID,]
alpha <- merge(alpha, submeta,by = "SampleID")
write.csv(alpha, 'alpha/alpha.csv', quote = FALSE)

### 首先观察幽门螺杆菌对胃肠道菌群的影响
alpha <- read.csv("alpha/alpha.csv",row.names = 1,header = T)
# hp infection
alpha$SampleType <- factor(alpha$SampleType ,levels = c("Gastric biopsies","Gastric juice","Stool samples"))
my_comparisons = list(c('HpN','HpP'))
# Shannon
(p <- ggplot(alpha, aes(x= HpInfection , y=Shannon))+
    geom_boxplot(aes(fill = HpInfection ),show.legend = F,outlier.colour = NA)+
    stat_boxplot(geom = 'errorbar',width=0.3, colour='black')+
    geom_jitter(alpha = 0.5,width = 0.1,size =1) + theme_bw() +
    facet_wrap( ~ SampleType,nrow = 1,scales = 'free_y')+
    scale_fill_manual(values = c("#4DBBD5B2","#E64B35B2"))+
    theme(axis.text.y = element_text(face = 'bold',colour = 'black')) +
    theme(axis.text = element_text(colour = 'black'))+
    theme(strip.text = element_text(face = 'bold'), 
          panel.grid=element_blank(),
          panel.background = element_rect(fill='white', colour='black'),
          strip.background = element_rect(fill = 'gray90',colour = 'black')) + 
    geom_signif(comparisons = my_comparisons, step_increase = 0.2, 
                map_signif_level = T, test = wilcox.test)  +
    labs(y="Shannon Index",x='Hp Infection'))
ggsave(paste0("alpha/boxplot_Shannon_hp.pdf"), p, width=180, height=150, units="mm")
ggsave(paste0("alpha/boxplot_Shannon_hp.jpg"), p, width=180, height=150, units="mm")

# 胃肠道菌群与疾病的关联   
# 木村竹本分类
# 胃镜诊断结果
alpha$GD <- factor(alpha$GD ,levels = c("SG","C1","C2","C3","O1","O2"))
alpha$SampleType <- factor(alpha$SampleType ,levels = c("Gastric biopsies","Gastric juice","Stool samples"))
my_comparisons = list(c('SG',"C1"),c("SG","C2"),c("SG","C3"),c("SG","O1"),c("SG","O2"))
(p <- ggplot(alpha, aes(x= GD , y= Shannon))+
  geom_boxplot(aes(fill = GD ),show.legend = F,outlier.colour = NA)+
  stat_boxplot(geom = 'errorbar',width=0.3, colour='black')+
  geom_jitter(alpha = 0.5,width = 0.1,size =1) + theme_bw() +
  facet_wrap( ~ SampleType,nrow = 1,scales = 'free_y')+
  scale_fill_brewer(palette = 'Set2')+
  theme(axis.text.y = element_text(face = 'bold',colour = 'black')) +
  theme(axis.text = element_text(colour = 'black'))+
  theme(strip.text = element_text(face = 'bold'), 
        panel.grid=element_blank(),
        panel.background = element_rect(fill='white', colour='black'),
        strip.background = element_rect(fill = 'gray90',colour = 'black')) + 
  geom_signif(comparisons = my_comparisons, step_increase = 0.2, map_signif_level = T, test = wilcox.test)  +
  stat_compare_means(method = "anova")+ # Add global p-value
  labs(y="Shannon Index",x='Endoscopic diagnosis results'))

ggsave(paste0("alpha/boxplot_Shannon_GD.pdf"), p, width=180, height=150, units="mm")
ggsave(paste0("alpha/boxplot_Shannon_GD.jpg"), p, width=180, height=150, units="mm")

# 去除幽门螺杆菌患者
sub_alpha <- subset(alpha, HpInfection != "HpP")
sub_alpha$GD <- factor(sub_alpha$GD ,levels = c("SG","C1","C2","C3","O1","O2"))
my_comparisons = list(c('SG',"C1"),c("SG","C2"),c("SG","C3"),c("SG","O1"),c("SG","O2"))
# Faith_pd
(p <- ggplot(sub_alpha, aes(x= GD , y= Faith_pd))+
  geom_boxplot(aes(fill = GD ),show.legend = F,outlier.colour = NA)+
  stat_boxplot(geom = 'errorbar',width=0.3, colour='black')+
  geom_jitter(alpha = 0.5,width = 0.1,size =1) + theme_bw() +
  facet_wrap( ~ SampleType,nrow = 1,scales = 'free_y')+
  scale_fill_brewer(palette = 'Set2')+
  theme(axis.text.y = element_text(face = 'bold',colour = 'black')) +
  theme(axis.text = element_text(colour = 'black'))+
  theme(strip.text = element_text(face = 'bold'), 
        panel.grid=element_blank(),
        panel.background = element_rect(fill='white', colour='black'),
        strip.background = element_rect(fill = 'gray90',colour = 'black')) + 
    stat_compare_means(method = "anova")+ # Add global p-value
      labs(y="Faith_pd Index",x='Endoscopic diagnosis results'))
ggsave(paste0("alpha/boxplot_Faith_pd_EA_-HP.pdf"), p, width=160, height=150, units="mm")
ggsave(paste0("alpha/boxplot_Faith_pd_EA_-HP.jpg"), p, width=160, height=150, units="mm")

### 按病理诊断结果分类
## 去除HP影响因素
sub_alpha <- subset(alpha, HpInfection != "HpP")
my_comparisons = list(c('SG',"AG"),c("SG","IM"),c("SG","IN"))
sub_alpha$PD <- factor(sub_alpha$PD,levels = c("SG","AG","IM","IN"))
sub_alpha$SampleType <- factor(sub_alpha$SampleType ,levels = c("Gastric biopsies","Gastric juice","Stool samples"))
## 分组计算 Faith _ pd正态性，方差齐性
sub_alpha_mucosa <- subset(sub_alpha, SampleType == "Gastric biopsies")
shapiro.test(sub_alpha_mucosa$Chao1)
bartlett.test(Chao1 ~PD,sub_alpha_mucosa) # p <0.05 ，k-square < chi-square
qchisq(0.95,3)
fligner.test(Chao1 ~ PD,sub_alpha_mucosa)
# Chao1
# 方差分析 -- fecal
summary(aov(Chao1 ~ PD ,subset(sub_alpha, SampleType == "Stool samples")))
# 方差分析 -- juice
summary(aov(Chao1 ~ PD ,subset(sub_alpha, SampleType == "Gastric juice")))
# 方差分析 -- mucosa
summary(aov(Chao1 ~ PD ,subset(sub_alpha, SampleType == "Gastric biopsies")))
tukey.reuslt <- as.data.frame(TukeyHSD(anova.result)$PD) # 多重检验及p值校正
tukey.reuslt$Qvalue <- p.adjust(tukey.reuslt$`p adj` , method = "BH"  )
tukey.reuslt # SG-IN : P 0.02466519
# ggplot2
(p <- ggplot(sub_alpha, aes(x= PD , y= Chao1))+
  geom_boxplot(aes(fill = PD ),show.legend = F,outlier.colour = NA)+
  stat_boxplot(geom = 'errorbar',width=0.3, colour='black')+
  geom_jitter(alpha = 0.5,width = 0.1,size =1) + theme_bw() +
  facet_wrap( ~ SampleType,nrow = 1,scales = 'free_y')+
  scale_fill_manual(values = c("#0073C2FF","#EFC000FF","#E18727FF","#A73030FF"))+
    theme(axis.text.y = element_text(face = 'bold',colour = 'black')) +
    theme(axis.text = element_text(colour = 'black'))+
    theme(strip.text = element_text(face = 'bold'), 
          panel.grid=element_blank(),
          panel.background = element_rect(fill='white', colour='black'),
          strip.background = element_rect(fill = 'gray90',colour = 'black')) +
    #geom_signif(comparisons = my_comparisons, step_increase = 0.2, map_signif_level = T, test = wilcox.test)  +
    stat_compare_means(method = "anova",hjust=0,label.y = 600)+
    labs(y="Chao1 Index",x='Pathological Classification'))
ggsave(paste0("alpha/boxplot_Chao1_PD_-HP.pdf"), p, width=150, height=100, units="mm")
ggsave(paste0("alpha/boxplot_Chao1_PD_-HP.jpg"), p, width=150, height=100, units="mm")

## 分组计算 Faith _ pd正态性，方差齐性
sub_alpha_mucosa <- subset(sub_alpha, SampleType == "Gastric biopsies")
shapiro.test(sub_alpha_mucosa$Chao1)
bartlett.test(Faith_pd ~PD,sub_alpha_mucosa) # p <0.05 ，k-square < chi-square
qchisq(0.95,3)
fligner.test(Faith_pd ~ PD,sub_alpha_mucosa)
# 方差分析 -- mucosa
anova.result <- aov(Faith_pd ~ PD ,sub_alpha_mucosa)
summary(anova.result)
tukey.reuslt <- as.data.frame(TukeyHSD(anova.result)$PD) # 多重检验及p值校正
tukey.reuslt$Qvalue <- p.adjust(tukey.reuslt$`p adj` , method = "BH"  )
tukey.reuslt # SG-IN P 0.04297743
# 方差分析 -- fecal
summary(aov(Faith_pd ~ PD ,subset(sub_alpha, SampleType == "Stool samples")))
# 方差分析 -- mucosa
summary(aov(Faith_pd ~ PD ,subset(sub_alpha, SampleType == "Gastric juice")))
# 根据显著组设置比较组
my_comparisons = list(c("SG","IN"))
# Faith_pd
(p <- ggplot(sub_alpha, aes(x= PD , y= Faith_pd))+
    geom_boxplot(aes(fill = PD ),show.legend = F,outlier.colour = NA)+
  stat_boxplot(geom = "errorbar",width=0.3, colour='black')+
  geom_jitter(alpha = 0.5,width = 0.1,size =1) + theme_bw() +
    facet_wrap( ~ SampleType,nrow = 1,scales = 'free_y')+
    scale_fill_manual(values = c("#0073C2FF","#EFC000FF","#E18727FF","#A73030FF"))+
  theme(legend.position = 'none')+
  theme(axis.text.y = element_text(face = 'bold',colour = 'black')) +
  theme(axis.text = element_text(colour = 'black'))+
  theme(strip.text = element_text(face = 'bold'), 
        panel.grid=element_blank(),
        panel.background = element_rect(fill='white', colour='black'),
        strip.background = element_rect(fill = 'gray90',colour = 'black')) +
    stat_compare_means(method = "anova",label.y = 35,hjust = 0)+
        labs(y="Faith PD Index",x='Pathological Classification'))
ggsave(paste0("alpha/boxplot_Faith_pd_PD_-HP.pdf"), p, width=150, height=100, units="mm")
ggsave(paste0("alpha/boxplot_Faith_pd_PD_-HP.jpg"), p, width=150, height=100, units="mm")

### AG,IM,IN合并为胃癌前病变组PLGC，SG为普通炎症组
alpha$SampleType <- factor(alpha$SampleType ,levels = c("Gastric biopsies","Gastric juice","Stool samples"))
alpha$PathogenicType <- factor(alpha$PathogenicType,levels = c("SG","PLGC"))
my_comparisons = list(c('SG',"PLGC"))
# Faith_pd
(p <- ggplot(alpha, aes(x= PathogenicType , y= Faith_pd))+
    geom_boxplot(aes(fill = PathogenicType ),show.legend = F,outlier.colour = NA)+
    stat_boxplot(geom = 'errorbar',width=0.3, colour='black')+
    geom_jitter(alpha = 0.5,width = 0.1,size =1) + 
    facet_wrap( ~ SampleType,nrow = 1,scales = 'free_y')+
    scale_fill_manual(values = c("#0073C2FF","#E77728"))+
    theme(legend.position = 'none')+
    theme(axis.text.y = element_text(face = 'bold',colour = 'black')) +
    theme(axis.text = element_text(colour = 'black'))+
    theme(strip.text = element_text(face = 'bold'),  
          panel.background = element_rect(fill='white', colour='black'),
          strip.background = element_rect(fill = 'gray90',colour = 'black')) +
    geom_signif(comparisons = my_comparisons, step_increase = 0.2,
                map_signif_level = T, test = wilcox.test)  +
    labs(y="Faith_pd Index",x='PLGC Classification'))
ggsave(paste0("alpha/boxplot_Faith_pd_PLGC.pdf"), p, width=110, height=150, units="mm")
ggsave(paste0("alpha/boxplot_Faith_pd_PLGC.jpg"), p, width=110, height=150, units="mm")
# Faith_pd
sub_alpha <- subset(alpha, HpInfection != "HpP" )
my_comparisons = list(c('SG',"PLGC"))
sub_alpha$PathogenicType <- factor(sub_alpha$PathogenicType,levels = c("SG","PLGC"))
sub_alpha$SampleType <- factor(sub_alpha$SampleType ,levels = c("Gastric biopsies","Gastric juice","Stool samples"))
# -hp
(p <- ggplot(sub_alpha, aes(x= PathogenicType , y= Faith_pd))+
  geom_boxplot(aes(fill = PathogenicType ),show.legend = F,outlier.colour = NA)+
  stat_boxplot(geom = 'errorbar',width=0.3, colour='black')+
  geom_jitter(alpha = 0.5,width = 0.1,size =1) + 
  facet_wrap( ~ SampleType,nrow = 1,scales = 'free_y')+
  scale_fill_manual(values = c("#0073C2FF","#E77728"))+
  theme(legend.position = 'none')+
  theme(axis.text.y = element_text(face = 'bold',colour = 'black')) +
  theme(axis.text = element_text(colour = 'black'))+
  theme(strip.text = element_text(face = 'bold'),  
        panel.background = element_rect(fill='white', colour='black'),
        strip.background = element_rect(fill = 'gray90',colour = 'black')) +
  geom_signif(comparisons = my_comparisons, step_increase = 0.2,
              map_signif_level = T, test = wilcox.test)  +
    labs(y="Faith_pd Index",x='PLGC Classification'))
ggsave(paste0("alpha/boxplot_Faith_pd_PLGC_-HP.pdf"), p, width=110, height=150, units="mm")
ggsave(paste0("alpha/boxplot_Faith_pd_PLGC_-HP.jpg"), p, width=110, height=150, units="mm")

#Chao1
# + hp
(p <- ggplot(alpha, aes(x= PathogenicType , y= Chao1))+
  geom_boxplot(aes(fill = PathogenicType ),show.legend = F,outlier.colour = NA)+
  stat_boxplot(geom = 'errorbar',width=0.3, colour='black')+
  geom_jitter(alpha = 0.5,width = 0.1,size =1) + theme_bw() +
  facet_wrap( ~ SampleType,nrow = 1,scales = 'free_y')+
  scale_fill_manual(values = c("#0073C2FF","#E77728"))+
    theme(legend.position = 'none')+
    theme(axis.text.y = element_text(face = 'bold',colour = 'black')) +
    theme(axis.text = element_text(colour = 'black'))+
    theme(strip.text = element_text(face = 'bold'),  
          panel.background = element_rect(fill='white', colour='black'),
          strip.background = element_rect(fill = 'gray90',colour = 'black')) +
    geom_signif(comparisons = my_comparisons, step_increase = 0.2,
                map_signif_level = T, test = wilcox.test)  +
    labs(y="Chao1 Index",x='PLGC Classification'))
ggsave(paste0("alpha/boxplot_Chao1_PLGC.pdf"), p, width=120, height=100, units="mm")
ggsave(paste0("alpha/boxplot_Chao1_PLGC.jpg"), p, width=120, height=100, units="mm")
#- hp
(p <- ggplot(sub_alpha, aes(x= PathogenicType , y= Chao1))+
  geom_boxplot(aes(fill = PathogenicType ),show.legend = F,outlier.colour = NA)+
  stat_boxplot(geom = 'errorbar',width=0.3, colour='black')+
  geom_jitter(alpha = 0.5,width = 0.1,size =1) + theme_bw() +
  facet_wrap( ~ SampleType,nrow = 1,scales = 'free_y')+
  scale_fill_manual(values = c("#0073C2FF","#E77728"))+
    theme(legend.position = 'none')+
    theme(axis.text.y = element_text(face = 'bold',colour = 'black')) +
    theme(axis.text = element_text(colour = 'black'))+
    theme(strip.text = element_text(face = 'bold'),  
          panel.background = element_rect(fill='white', colour='black'),
          strip.background = element_rect(fill = 'gray90',colour = 'black')) +
    geom_signif(comparisons = my_comparisons, step_increase = 0.2,
                map_signif_level = T, test = wilcox.test)  +
    labs(y="Chao1 Index",x='PLGC Classification'))
ggsave(paste0("alpha/boxplot_Chao1_PLGC_-HP.pdf"), p, width=120, height=100, units="mm")
ggsave(paste0("alpha/boxplot_Chao1_PLGC_-HP.jpg"), p, width=120, height=100, units="mm")

## OLGA
sub_alpha <- subset(alpha, HpInfection != "HpP")
sub_alpha$SampleType <- factor(sub_alpha$SampleType ,levels = c("Gastric biopsies","Gastric juice","Stool samples"))
sub_alpha$OLGA <- factor(sub_alpha$OLGA,levels = c("0","1","2","3"))
# anova
sub_alpha_mucosa$OLGA <- factor(sub_alpha_mucosa$OLGA,levels = c("0","1","2","3"))
anova.result <- aov(Faith_pd ~ OLGA ,sub_alpha_mucosa)
summary(anova.result)
tukey.reuslt <- as.data.frame(TukeyHSD(anova.result)$OLGA) # 多重检验及p值校正
my_comparisons = list(c('0',"1"),c('0',"2"),c('1',"2"))
# Faith_pd
(p <- ggplot(sub_alpha, aes(x= OLGA , y= Faith_pd))+
    geom_boxplot(aes(fill = OLGA ),show.legend = F,outlier.colour = NA)+
    stat_boxplot(geom = 'errorbar',width=0.3, colour='black')+
    geom_jitter(alpha = 0.5,width = 0.1,size =1) + theme_bw() +
    facet_wrap( ~ SampleType,nrow = 1,scales = 'free_y')+
    scale_fill_manual(values = c("#0073C2FF","#EFC000FF","#E18727FF","#A73030FF"))+
    theme(legend.position = 'none')+
    theme(axis.text.y = element_text(face = 'bold',colour = 'black')) +
    theme(axis.text = element_text(colour = 'black'))+
    theme(strip.text = element_text(face = 'bold'),  
          panel.background = element_rect(fill='white', colour='black'),
          strip.background = element_rect(fill = 'gray90',colour = 'black')) +
    stat_compare_means(method = "anova",label.y = 35,hjust = 0)+
    labs(y="Faith_pd Index",x='OLGA scores'))
ggsave(paste0("alpha/boxplot_Faith_pd_OLGA_-HP.pdf"), p, width=150, height=100, units="mm")
ggsave(paste0("alpha/boxplot_Faith_pd_OLGA_-HP.jpg"), p, width=150, height=100, units="mm")

## CHAO1
# anova
sub_alpha_mucosa$OLGA <- factor(sub_alpha_mucosa$OLGA,levels = c("0","1","2","3"))
anova.result <- aov(Chao1 ~ OLGA ,sub_alpha_mucosa)
summary(anova.result)
tukey.reuslt <- as.data.frame(TukeyHSD(anova.result)$OLGA) # 多重检验及p值校正
# ggplot
(p <- ggplot(sub_alpha, aes(x= OLGA , y= Chao1))+
    geom_boxplot(aes(fill = OLGA ),show.legend = F,outlier.colour = NA)+
    stat_boxplot(geom = 'errorbar',width=0.3, colour='black')+
    geom_jitter(alpha = 0.5,width = 0.1,size =1) + theme_bw() +
    facet_wrap( ~ SampleType,nrow = 1,scales = 'free_y')+
    scale_fill_manual(values = c("#0073C2FF","#EFC000FF","#E18727FF","#A73030FF"))+
    theme(legend.position = 'none')+
    theme(axis.text.y = element_text(face = 'bold',colour = 'black')) +
    theme(axis.text = element_text(colour = 'black'))+
    theme(strip.text = element_text(face = 'bold'),  
          panel.background = element_rect(fill='white', colour='black'),
          strip.background = element_rect(fill = 'gray90',colour = 'black')) +
    stat_compare_means(method = "anova",label.y = 600,hjust = 0)+
    labs(y="Chao1 Index",x='OLGA'))
ggsave(paste0("alpha/boxplot_Chao1_OLGA_-HP.pdf"), p, width=120, height=120, units="mm")
ggsave(paste0("alpha/boxplot_Chao1_OLGA_-HP.jpg"), p, width=120, height=120, units="mm")

## OLGIM
sub_alpha <- subset(alpha, HpInfection != "HpP")
sub_alpha$SampleType <- factor(sub_alpha$SampleType ,levels = c("Gastric biopsies","Gastric juice","Stool samples"))
sub_alpha$OLGIM <- factor(sub_alpha$OLGIM,levels = c("0","1","2","3","4"))
# anova MUCOSA
sub_alpha_mucosa$OLGIM <- factor(sub_alpha_mucosa$OLGIM,levels = c("0","1","2","3"))
anova.result <- aov(Faith_pd ~ OLGIM ,sub_alpha_mucosa)
summary(anova.result)
tukey.reuslt <- as.data.frame(TukeyHSD(anova.result)$OLGIM) # 多重检验及p值校正
# anova
sub_alpha_juice <- subset(sub_alpha,SampleType == "Gastric juice")
sub_alpha_juice$OLGIM <- factor(sub_alpha_juice$OLGIM,levels = c("0","1","2","3"))
anova.result <- aov(Faith_pd ~ OLGIM ,sub_alpha_juice)
summary(anova.result)
tukey.reuslt <- as.data.frame(TukeyHSD(anova.result)$OLGIM) # 多重检验及p值校正
my_comparisons = list(c('0',"2"))
# Faith_pd
(p <- ggplot(sub_alpha, aes(x= OLGIM , y= Faith_pd))+
    geom_boxplot(aes(fill = OLGIM ),show.legend = F,outlier.colour = NA)+
    stat_boxplot(geom = 'errorbar',width=0.3, colour='black')+
    geom_jitter(alpha = 0.5,width = 0.1,size =1) + theme_bw() +
    facet_wrap( ~ SampleType,nrow = 1,scales = 'free_y')+
    scale_fill_manual(values = c("#0073C2FF","#EFC000FF","#E18727FF","#A73030FF"))+
    theme(legend.position = 'none')+
    theme(axis.text.y = element_text(face = 'bold',colour = 'black')) +
    theme(axis.text = element_text(colour = 'black'))+
    theme(strip.text = element_text(face = 'bold'),  
          panel.background = element_rect(fill='white', colour='black'),
          strip.background = element_rect(fill = 'gray90',colour = 'black')) +
    stat_compare_means(method = "anova",label.y = 35,hjust = 0)+
    labs(y="Faith_pd Index",x='OLGIM scores'))
ggsave(paste0("alpha/boxplot_Faith_pd_OLGIM_-HP.pdf"), p, width=150, height=100, units="mm")
ggsave(paste0("alpha/boxplot_Faith_pd_OLGIM_-HP.jpg"), p, width=150, height=100, units="mm")
## CHAO1
# anova
sub_alpha_mucosa$OLGIM <- factor(sub_alpha_mucosa$OLGIM,levels = c("0","1","2","3"))
anova.result <- aov(Chao1 ~ OLGIM ,sub_alpha_mucosa)
summary(anova.result)
tukey.reuslt <- as.data.frame(TukeyHSD(anova.result)$OLGIM) # 多重检验及p值校正
# GGPLOT
(p <- ggplot(sub_alpha, aes(x= OLGIM , y= Chao1))+
    geom_boxplot(aes(fill = OLGIM ),show.legend = F,outlier.colour = NA)+
    stat_boxplot(geom = 'errorbar',width=0.3, colour='black')+
    geom_jitter(alpha = 0.5,width = 0.1,size =1) + theme_bw() +
    facet_wrap( ~ SampleType,nrow = 1,scales = 'free_y')+
    scale_fill_manual(values = c("#0073C2FF","#EFC000FF","#E18727FF","#A73030FF"))+
    theme(legend.position = 'none')+
    theme(axis.text.y = element_text(face = 'bold',colour = 'black')) +
    theme(axis.text = element_text(colour = 'black'))+
    theme(strip.text = element_text(face = 'bold'),  
          panel.background = element_rect(fill='white', colour='black'),
          strip.background = element_rect(fill = 'gray90',colour = 'black')) +
    stat_compare_means(method = "anova",label.y = 600,hjust = 0)+
    labs(y="Chao1 Index",x='OLGIM'))
ggsave(paste0("alpha/boxplot_Chao1_OLGIM_-HP.pdf"), p, width=120, height=120, units="mm")
ggsave(paste0("alpha/boxplot_Chao1_OLGIM_-HP.jpg"), p, width=120, height=120, units="mm")


