setwd("E:/data")
env = read.csv("env.csv", row.names = 1)

env$a <- (env$BG + env$CBH)/(env$BG + env$CBH + env$AKP)
env$b <- (env$BG + env$CBH)/(env$BG + env$CBH + env$NAG)

env$vector_length <- sqrt((env$a^2+env$b^2)) # Computed vector length

# Computed vector Angle
env$vector_angle <- atan2(env$a, env$b) * 180 / pi - 45# Computed vector Angle
# Add the result as a new column to the env data box

env$vector_length1 <- sqrt((log(env$BG) / log(env$NAG + env$LAP))^2 + (log(env$BG) / log(env$AKP))^2)
env$vector_angle1 <- atan2(log(env$BG) / log(env$AKP), log(env$BG) / log(env$NAG + env$LAP)) * 180 / pi- 45

write.csv(env,"env1.csv")

env = read.csv("env1.csv", row.names = 1)

env$MCN <- (env$BG)/(env$NAG)
env$MNP <- (env$NAG)/(env$AKP)

env$a1 <- (env$BG)/(env$LAP + env$NAG)
env$b1 <- (env$LAP + env$NAG)/(env$AKP)
env <- env[!(rownames(env) %in% c('Bacteria_F281_raw', 'Bacteria_F282_raw',
                                  'Bacteria_F283_raw')), ]
library(ggpmisc)
library(ggplot2)
unique(env$type)

p1 <- ggplot(env, aes(MAP, MCN)) +
  #geom_point(shape = 21, size = 4, colour = "black", stroke = 1) +
  #scale_fill_gradient(low = "#BAC8F5", high = "#5E7CF2")+
  #scale_color_manual(values = c("#f85a40","#00c16e","#037ef3"))+
  geom_jitter(aes(fill=type), shape=21,size=6,alpha=0.6)+
  scale_fill_manual(values = c("#f85a40","#00c16e","#037ef3"))+
  
  theme_bw(base_line_size = 1.05,base_rect_size = 1.05)+
  xlab("MAP") + ylab("BG:NAG")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))+
  theme_bw()+
  theme(legend.position = "none")+
  theme(strip.text = element_text(size = 12,colour = "black"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+ 
  geom_smooth(method=lm,level=0.95,color="black",size=1.4,se=T)+#
  #facet_wrap(~"variable", scale = 'free')+
  #labeller = labeller(variable = custom_labels)) +
    stat_poly_eq(method = 'lm',label.x.npc = 0.2, label.y.npc = 0.7,
               aes(label=paste(after_stat(rr.label),
                               after_stat(p.value.label),sep = "*\", \"*")),
               size=4)
#geom_vline(xintercept = 1, linetype = "dashed", color = "black",size = 1)

p1

p2 <- ggplot(env, aes(MAP, MNP)) +
  #geom_point(shape = 21, size = 4, colour = "black", stroke = 1) +
  #scale_fill_gradient(low = "#BAC8F5", high = "#5E7CF2")+
  #scale_color_manual(values = c("#f85a40","#00c16e","#037ef3"))+
  geom_jitter(aes(fill=type), shape=21,size=6,alpha=0.6)+
  scale_fill_manual(values = c("#f85a40","#00c16e","#037ef3"))+
  
  theme_bw(base_line_size = 1.05,base_rect_size = 1.05)+
  xlab("MAP") + ylab("NAG:AP")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))+
  theme_bw()+
  #theme(legend.position = "none")+#
  theme(strip.text = element_text(size = 12,colour = "black"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        legend.position = c(0.3, 0.85),
        legend.direction = "horizontal",
        legend.title = element_blank())+ 
  geom_smooth(method=lm,level=0.95,color="black",size=1.4,se=T)+#
  #facet_wrap(~variable, ncol = 2, scale = 'free',
  #labeller = labeller(variable = custom_labels)) +
  stat_poly_eq(method = 'lm',label.x.npc = 0.1, label.y.npc = 0.7,
               aes(label=paste(after_stat(rr.label),
                               after_stat(p.value.label),sep = "*\", \"*")),
               size=4)
#geom_vline(xintercept = 1, linetype = "dashed", color = "black",size = 1)

p2

library(gridExtra)

# 组合图形
combined <- grid.arrange(p1, p2, ncol = 2)
