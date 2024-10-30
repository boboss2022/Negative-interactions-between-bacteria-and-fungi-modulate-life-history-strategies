library(rfPermute)
library(ggplot2)
#Read data
df = read.csv("sem1.csv",header=T)
names(df)
df <- df[, -c(1, 8,18,28:36)]

set.seed(123)
#Running random forest
rf_results1<-rfPermute(CopyNumber~., data =df, importance = TRUE, ntree = 1000)

#Look at the main results of the random forest
rf_results1$rf
#The interpretation rate of the predictor was extracted
predictor_var1<- data.frame(importance(rf_results1, scale = TRUE), check.names = FALSE)
#write.csv(predictor_var1,"随机结果-1.csv")
#Extract the significance of predictor variables
predictor_sig1<-as.data.frame((rf_results1$pval)[,,2])
colnames(predictor_sig1)<-c("sig","other")
#Combine significance factor and interpretation rate tables
df_pre1<-cbind(predictor_var1,predictor_sig1)
df_pre1$sig1 <- ifelse(df_pre1$sig >= 0.01 & df_pre1$sig <= 0.05, "*", 
                       ifelse(df_pre1$sig >= 0.001 & df_pre1$sig < 0.01, "**", 
                              ifelse(df_pre1$sig < 0.001, "***", " ")))
write.csv(df_pre1,"Random result2.csv")

#k1 <- df_pre1$IncNodePurity[df_pre1$sig1=="*"]<-"#3679D1"

df_pre1 <- read.csv("Random result2.csv",row.names = 1)

df_pre1$row <- factor(rownames(df_pre1), levels = rev(rownames(df_pre1)))
#df_pre1$row <- factor(df_pre1$row, levels = rev(levels(df_pre1$row)))
df_pre1 <- df_pre1[nrow(df_pre1):1, ]

custom_labels <- rev(c("MAP" = "MAP", "TN" = "TN",
                   "TC" = "TC", "bf_com" = "Bacteria-fungus competition",
                   #"nw1" = "Niche width(Bacterial)", "nw2" = "Niche width(fungus)",
                   "bac_com" = "Bacterial competition","fun_com" = "Fungus competition",
                   "pH" = "pH", "Sand" = "Sand","Shannon1" = "Shannon(Bacterial)", "TP" = "TP",
                   "WC" = "WC","silt" = "Silt",
                   "bac" = "Network complexity(Bacterial)","fun" = "Network complexity(fungus)",
                   "bf" = "Network complexity(interaction)",
                   "Richness1" = "Richness(Bacterial)","Clay" = "Clay",
                   "pcoa1" = "Pcoa(Bacterial)", #"type" = "Habitat type",
                   "PD1" = "PD(Bacterial)","Richness2" = "Richness(fungus)","pd2" = "PD(fungus)",
                   "pcoa2" = "Pcoa(fungus)","Shannon2" = "Shannon(fungus)"
                   ))
#Draw a bar chart
names(df_pre1)
p1 <- ggplot(data=df_pre1,aes(x=df_pre1$row,y=df_pre1$X.IncMSE,fill=type))+
  geom_bar(stat = 'identity',width=0.4)+
  
  theme_classic()+ labs(x='',y='Increase in MSE(%)')+
  
  scale_y_continuous(expand = c(0,0),limit = c(0, 30))+
  
  geom_text(aes(label = df_pre1$sig1,y= df_pre1$X.IncMSE+ 0.8, x =df_pre1$row),vjust = 0.5,size = 6)+
  
  theme(axis.text = element_text(color = "black",size = 11))+ 
  annotate('text', label = sprintf('italic(R^2) == %s', "0.752"), x = 15, y = 22, size = 5, parse = TRUE)+
    guides(fill = guide_legend(title = NULL))+
  theme(legend.position = c(0.8, 0.7)) +
  #theme(legend.direction = "horizontal") +
  #scale_fill_discrete(labels = c("Abiotic factors","Biological factors"))
  theme(axis.text.x = element_text( hjust = 1)) +
  theme(axis.text.y=element_text(size=12))+#调整y轴刻度
  theme(axis.text.x=element_text(size=12))+#调整x轴刻度
  theme(axis.title.x=element_text(size=12))+#调整x轴标签
  theme(axis.title.y=element_text(size=12))+#调整x轴标签
  scale_fill_discrete(labels = c("Abiotic factors","Biological factors"))+
  theme(legend.text = element_text(size = 12))+
  scale_x_discrete(labels = custom_labels)+ coord_flip()
p1

#Construction of structural equation model
library(tidyverse)
library(nlme)
library(glmmTMB)
library(piecewiseSEM)
#remove.packages('glmmTMB')#xfun，htmltools，stringi，purrr，
#tzdb，minqa，tidyr，lme4，mvtnorm，readr，igraph，quantreg
#install.packages('TMB', type = 'source')

df = read.csv("sem1.csv",header=T)
#.libPaths("C:/Program Files/R/R-4.3.0/library")
#devtools::install_github("jslefche/piecewiseSEM@devel")

scale2 <- function(x) as.numeric(scale(x))
#install.packages("purrr")


df <- df |> 
  mutate_if(is.numeric, scale2) |>
  as.data.frame()
#soil
soil <- lm(CopyNumber ~ TC+TN+TP+pH+WC+Sand, df)

summary(soil)
coefs(soil, standardize = 'scale')
#Extract coefficients to generate factor scores to predict target variables
S_TC <-  summary(soil)$coefficients[2, 1]
S_TN <-  summary(soil)$coefficients[3, 1]
S_TP <-  summary(soil)$coefficients[4, 1]
S_pH <-  summary(soil)$coefficients[5, 1]
S_WC <-  summary(soil)$coefficients[6, 1]
S_Sand <-  summary(soil)$coefficients[7, 1]

soil <- S_TC*df$TC + S_TN*df$TN+ S_TP*df$TP+ S_pH*df$pH+ S_WC*df$WC+ S_Sand*df$Sand
df$soil <- soil
summary(lm(CopyNumber ~ soil, df))
#Use the function in PecewiseSEM to get the normalized coefficient: Coefs
coefs(lm(CopyNumber ~ soil, df))
names(df)

#Diversity index
div <- lm(CopyNumber ~ Shannon1+pcoa1+PD1, df)

summary(div)
coefs(div, standardize = 'scale')

d_s1 <-  summary(div)$coefficients[2, 1]
d_s2 <-  summary(div)$coefficients[3, 1]
d_s3 <-  summary(div)$coefficients[4, 1]

div <- d_s1*df$Shannon1 + d_s2*df$pcoa1 +d_s3*df$PD1
df$div <- div
summary(lm(CopyNumber ~ div, df))
#Use the function in PecewiseSEM to get the normalized coefficient: Coefs
coefs(lm(CopyNumber ~ div, df))
names(df)

site.cv0 <- psem(
  lme(soil ~ MAP,
      random = ~ 1 | type/site ,method="REML",
      #random = ~ 1 | Reference/ID,method="REML",
      na.action = na.omit,
      data = df),
  lme(div ~ soil + MAP,
      random = ~ 1 | type/site ,method="REML",
      #random = ~ 1 | Reference/ID,method="REML",
      na.action = na.omit,
      data = df),
  lme(bf_com ~ soil + div,
      random = ~ 1 | type/site ,method="REML",
      #random = ~ 1 | Reference/ID,method="REML",
      na.action = na.omit,
      data = df),
  lme(bf ~ soil + div,
      random = ~ 1 | type/site ,
      #random = ~ 1 | Reference/ID,method="REML",
      na.action = na.omit,method="REML",
      data = df),
  lme(CopyNumber ~ soil + bf_com+bf +div + MAP,
      random = ~ 1 | type/site ,method="REML",
      #random = ~ 1 | Reference/ID,method="REML",
      na.action = na.omit,
      data = df),
  bf %~~% bf_com
)
# 
summary(site.cv0, .progressBar = FALSE)

#SEM

plot(site.cv0)

