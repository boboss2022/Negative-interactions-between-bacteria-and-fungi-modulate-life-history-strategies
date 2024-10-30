##rrnDBcorrectOTU 包
#devtools::install_github("Listen-Lii/rrnDBcorrectOTU".force = TRUE)

##使用 rrnDBcorrectOTU 包，读取上述数据表格后校正 OTU 表的 16S 丰度
library(rrnDBcorrectOTU)

otu = read.csv("otu.csv", row.names = 1)
#RDP classifer 
classifer = read.csv(file="tax.csv",header=T, row.names = 1)
classifer <- read.table("taxonomy.txt", header = TRUE, row.names = 1)
#rrnDB
rrnDB = read.csv(file="rrnDB-5.8_pantaxa_stats_RDP.tsv",header=T,sep="\t")
rrnDB = read.table(file="rrnDB-5.8_pantaxa_stats_RDP.txt",header=T,sep="\t")

result = rco(otu,classifer,rrnDB)
write.table(result$whole.res, 'whole.res.txt', sep = '\t', col.names = NA, quote = FALSE)
write.table(result$correct.table, 'correct.table.txt', sep = '\t', col.names = NA, quote = FALSE)

A <-read.table(file="whole.res.txt",header=T,sep="\t")
averages <- A[,c(1,4)]
averages[averages < 0] <- 0

otu = A[,-2:-4]

rownames(otu) <- otu[, 1]
otu <- otu[, -1]


row = as.numeric(length(row.names(otu)))

col = as.numeric(length(colnames(otu)))

col_sum = rep(colSums(otu), row)

col_sum = matrix(col_sum, nrow = col, ncol = row)

ASV_relative=otu/t(col_sum)
#Check whether the result of each column is 1
colSums(ASV_relative)
# CK1 CK2 CK3 LT1 LT2 LT3 MT1 MT2 MT3 HT1 HT2 HT3 
# 1   1   1   1   1   1   1   1   1   1   1   1 

write.csv(ASV_relative,"ASV_relative.csv")

ASV <- read.csv("ASV_relative.csv", header=T)
library(dplyr)
weighted_CopyNumber <- ASV %>%
  inner_join(averages, by = "ID") %>%
  mutate(across(starts_with("Bacteria_"), ~ . * CopyNumber)) %>%
  rowwise() %>%
  mutate(weighted_CopyNumber = sum(c_across(starts_with("Bacteria_")))) %>%
  ungroup()

a <- colSums(ASV[,-1], na.rm = TRUE)
head(a)

col_names_a <- colnames(ASV[,-1])

b <- colSums(weighted_CopyNumber[,2:109], na.rm = TRUE)
col_names_b <- colnames(weighted_CopyNumber[,2:109])
c <- b / a


result_df <- data.frame(sample_id = col_names_a, weighted_CopyNumber = b)

write.table(result_df, file = "weighted_CopyNumber.txt", sep = "\t", row.names = FALSE)

library(ggplot2)
setwd("E:/要完成的文章/第二篇文章/data")
df = read.csv('function.csv',header = TRUE)
df <- df[,-1]

result1=c()#Create an empty vector to store the extracted results
for (i  in  c(2:5)){#lmer(MAP~df[,i]+(1|Type/site),data=data)
  fit=lm(MAP~df[,i],data=df)
  result1=rbind(result1,c(colnames(df)[i],coef(summary(fit))[2,c(1,2,4)]))
}
result1


R2_adj <- c()
p_value <- c()

for (i in c(2:5)) {
  fit_stat <- summary(lm(MAP~df[,i],data=df))  #Unitary linear regression
  R2_adj <- c(R2_adj, fit_stat$adj.r.squared)  #The corrected R2 was extracted
  p_value <- c(p_value, fit_stat$coefficients[2,4])  #The significance P-value was extracted
}

env_stat <- data.frame(row.names = names(df[, 2:5]), R2_adj, p_value)
env_stat  

df1 <- df[, c(1,5)]
dat_plot <- reshape2::melt(df1, id = 'MAP')

library(ggplot2)
library(ggpmisc)
library(reshape2)

# Customize the content of the faceted label
custom_labels <- c("weighted_genome_size" = "Genome size", "weighted_gc_percentage" = "GC content (%)",
                   "weighted_protein_count" = "Protein counts", "weighted_CopyNumber" = "Average rRNA operon copy number"
)

p <- ggplot(dat_plot, aes(MAP, value,color=MAP)) +
  #geom_point(shape = 21, size = 4, colour = "black", stroke = 1) +
  #scale_fill_gradient(low = "#BAC8F5", high = "#5E7CF2")+
  #geom_smooth(method = lm, se = TRUE, color = "black", fill = "grey")+
  geom_jitter(position=position_jitter(0.17), size=3,alpha=0.5)+
  
  theme_bw(base_line_size = 1.05,base_rect_size = 1.05)+
  
  geom_smooth(method=lm,level=0.95,color="#3D3E83",size=1.4,se=T)+#
  xlab("MAP") + ylab("")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))+
  theme_bw()+
  #theme(legend.position = "none")+#
  theme(strip.text = element_text(size = 12,colour = "black"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+ 
  facet_wrap(~variable, ncol = 2, scale = 'free',
             labeller = labeller(variable = custom_labels)) +
  scale_color_gradientn(values=seq(0,1,0.2),colors = c("gray98","#6191BF","#3D3E83"))+
  stat_poly_eq(method = 'lm',label.x.npc = 0.08, label.y.npc = 0.1,
               aes(label=paste(after_stat(rr.label),
                               after_stat(p.value.label),sep = "*\", \"*")),
               size=4)
p

library(ggplot2)
library(ggsignif)
df = read.csv('function.csv',header = TRUE)
df <- df[,-1]
df1 <- df[, c(1,5,8)]
dat_plot1 <- reshape2::melt(df1, id = c('type',"MAP"))
#Define comparison group
my.comparisons <- list(c("Desert", "Grassland"), c("Desert", "Mountain"), c("Grassland", "Mountain")) 
unique(dat_plot1$type)
unique(dat_plot1$variable)

p1 <- ggplot(dat_plot1, aes(type, value))+
  stat_boxplot(geom = "errorbar", width=0.1,size=0.8)+
  geom_boxplot(aes(fill=type), outlier.color = NA)+
  ggsignif::geom_signif(comparisons = my.comparisons, test = "wilcox.test", map_signif_level = FALSE,
                        step_increase = 0.15, size = 0.2, textsize = 4
                        , vjust = 1.5) + 
  geom_jitter(aes(fill=type), shape=21,size=2.5,alpha=0.6,width=0.3)+
  scale_fill_manual(values = c("#f85a40","#00c16e","#037ef3"))+
  facet_wrap(~variable, ncol = 2, scale = 'free',
             labeller = labeller(variable = custom_labels)) +
  theme_bw()+
  #theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  
  #theme(strip.background = element_rect(fill=c("#FAE3AD")))+
  xlab("") + ylab("")+
  theme(strip.text = element_text(size = 12,colour = "black"))+
  theme(legend.position = "none",
        axis.text = element_text(size=10),
        #strip.text = element_text(margin = margin(b = 5)),
        axis.text.x = element_text(vjust = 1, hjust = 1))
p1


