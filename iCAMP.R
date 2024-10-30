
# install.packages('NST')
library(NST)
library(iCAMP)
dir.create("./data/iCAMP.data/")
save.wd="./result_and_plot/Data.mining//iCAMP.ITs/"
if(!dir.exists(save.wd)){fs::dir_create(save.wd)}



# nworker is the number of threads and memory.g is the running memory, which is set according to the amount of analysis data, computing resources, and the desired run time.
prefix="Test"
rand.time=10
nworker=4
memory.G=30

library(iCAMP)
library(ape)
library(picante)

# Prepare data: otu table, evolution tree, grouping file, environment factor file #-------
library(tidyverse)
library(phyloseq)
library(ggClusterNet)

ps01 = base::readRDS("E:/amplicon_v.3/data/Rawdata//ps_16s001.rds")%>%
  filter_taxa(function(x) sum(x ) > 0 , TRUE)

ps01 = base::readRDS("E:/amplicon_v.3/data/Rawdata_fa//ITs.rds")%>%
  filter_taxa(function(x) sum(x ) > 0 , TRUE)
comm= ps01 %>% vegan_otu()
comm[1:5,1:5]
tree=ps01 %>% 
  filter_OTU_ps(500) %>% phy_tree() %>% multi2di() 
treat=as.tibble(sample_data(ps01))
row.names(treat) = treat$SampleID
env=read.csv("E:/amplicon_v.3/data/d1/env.csv",row.names = 1)
env

#  This part of the code is less important than phylsoeq
sampid.check = match.name(rn.list=list(comm=comm,treat=treat,env=env))
treat=sampid.check$treat
comm=sampid.check$comm
comm=comm[,colSums(comm)>0,drop=FALSE]
env=sampid.check$env
spid.check = match.name(cn.list=list(comm=comm),tree.list=list(tree=tree))
comm=spid.check$comm
tree=spid.check$tree


# Phylogenetic distance matrix between phylogenetic species
if(!file.exists("E:/result_and_plot/Data.mining//iCAMP.ITs/pd.desc")) {
  pd.big=iCAMP::pdist.big(tree = tree, wd=save.wd, nworker = nworker, memory.G = memory.G)
}else{
  pd.big=list()
  pd.big$tip.label=read.csv(paste0(save.wd,"/pd.taxon.name.csv"),row.names = 1,stringsAsFactors = FALSE)[,1]
  pd.big$pd.wd=save.wd
  pd.big$pd.file="pd.desc"
  pd.big$pd.name.file="pd.taxon.name.csv"}




# today <- Sys.Date()
# mydate <- strftime(today,format='%Y-%m-%d')

#BiocManager::install("random")
library("random")
tem0 = as.vector(randomStrings(n = 1, len = 5))
tem = paste0("./result_and_plot/Data.mining//iCAMP.16s/",tem0)


# The ecological niche difference between species was calculated according to the environmental variables of species
niche.dif=iCAMP::dniche(env = env,
                        comm = comm,
                        method = "niche.value",
                        nworker = nworker,out.dist=FALSE,
                        bigmemo=TRUE,
                        nd.wd= tem)

# Evaluation of within-bin phylogenetic signals
ds = 0.2
bin.size.limit = 5
phylobin=taxa.binphy.big(tree = tree, pd.desc = pd.big$pd.file,pd.spname = pd.big$tip.label,
                         pd.wd = pd.big$pd.wd, ds = ds, bin.size.limit = bin.size.limit,
                         nworker = nworker)


sp.bin=phylobin$sp.bin[,3,drop=FALSE]
sp.ra=colMeans(comm/rowSums(comm))
abcut=3
commc=comm[,colSums(comm)>=abcut,drop=FALSE]
dim(commc)
spname.use=colnames(commc)

binps=iCAMP::ps.bin(sp.bin = sp.bin,sp.ra = sp.ra,spname.use = spname.use,
                    pd.desc = pd.big$pd.file, pd.spname = pd.big$tip.label, pd.wd = pd.big$pd.wd,
                    nd.list = niche.dif$nd,nd.spname = niche.dif$names,ndbig.wd = niche.dif$nd.wd,
                    cor.method = "pearson",r.cut = 0.1, p.cut = 0.05, min.spn = 5)


if(file.exists(paste0("./result_and_plot/Data.mining//iCAMP.16s//",prefix,".PhyloSignalSummary.csv"))){
  appendy=TRUE;col.namesy=FALSE
}else{
  appendy=FALSE;col.namesy=TRUE}


write.table(data.frame(ds=ds,n.min=bin.size.limit,binps$Index),file = paste0("./result_and_plot/Data.mining//iCAMP.16s//",prefix,".PhyloSignalSummary.csv"),
            append = appendy, quote=FALSE, sep=",", row.names = FALSE,col.names = col.namesy)

if(file.exists(paste0("./result_and_plot/Data.mining//iCAMP.16s//",prefix,".PhyloSignalDetail.csv"))){appendy2=TRUE;col.namesy2=FALSE
}else{
  appendy2=FALSE;col.namesy2=TRUE}


write.table(data.frame(ds=ds,n.min=bin.size.limit,binID=rownames(binps$detail),binps$detail),
            file = paste0(prefix,".PhyloSignalDetail.csv"),
            append = appendy2, quote = FALSE, sep = ",", row.names = FALSE, col.names = col.namesy2)



bin.size.limit = 12
sig.index="Confidence"
# Analysis and quantification of relative importance of different processes based on phylogenetic zero model #--------
#Analysis of community construction based on Ning et al. (2020) method
#ses.cut =1.96, βNRI=1.96 as the threshold of homogeneous and heterogeneous selection; Rc.cut =0.95, RC=0.95 as the diffusion and drift threshold
# bin.sies.limit is used to specify the number of taxon in bin. This dataset is too small, so 5 is used. In practice, a reasonable bin can be selected according to phylogenetic signal testing or some Settings are tried
#rand specifies the number of randomizations to build a zero distribution, which nworker uses for multithreading to increase computation speed
icres=iCAMP::icamp.big(comm=comm, 
                       pd.desc = pd.big$pd.file,
                       pd.spname=pd.big$tip.label,
                       pd.wd = pd.big$pd.wd, 
                       rand = rand.time, 
                       tree=tree,
                       prefix = prefix,
                       ds = 0.2, pd.cut = NA, sp.check = TRUE,
                       phylo.rand.scale = "within.bin", taxa.rand.scale = "across.all",
                       phylo.metric = "bMPD", sig.index=sig.index, bin.size.limit = bin.size.limit, 
                       nworker = nworker, memory.G = memory.G, rtree.save = FALSE, detail.save = TRUE, 
                       qp.save = FALSE, detail.null = FALSE, ignore.zero = TRUE, output.wd = save.wd, 
                       correct.special = TRUE, unit.sum = rowSums(comm), special.method = "depend",
                       ses.cut = 1.96, rc.cut = 0.95, conf.cut=0.975, omit.option = "no",meta.ab = NULL)

head(icres$CbMPDiCBraya)


write.csv(icres$CbMPDiCBraya, './result_and_plot/Data.mining//iCAMP.16s/icamp.out.csv', row.names = FALSE)




map = sample_data(ps01)
# head(map)

treat.use = data.frame( as.numeric(as.factor(map$Group))) %>%as.matrix()
row.names(treat.use) = row.names(map)

pnstout=NST::pNST(comm=comm, 
                  pd.desc=pd.big$pd.file, 
                  pd.wd=pd.big$pd.wd, 
                  pd.spname=pd.big$tip.label,
                  group=treat.use,
                  abundance.weighted=TRUE,
                  rand=rand.time, phylo.shuffle=TRUE, 
                  nworker=nworker,
                  output.rand = TRUE, SES=FALSE, RC=FALSE)

write.csv(pnstout$index.grp,file = paste0("./result_and_plot/Data.mining//iCAMP.16s//",prefix,".pNST.summary.",colnames(treat)[1],".csv"))
write.csv(pnstout$index.pair.grp,file = paste0("./result_and_plot/Data.mining//iCAMP.16s//",prefix,".pNST.pairwise.",colnames(treat)[1],".csv"))

pnst.bt=NST::nst.boot(nst.result=pnstout, group=treat.use,
                      rand=rand.time, nworker=nworker)
write.csv(pnst.bt$NST.summary,file = paste0("./result_and_plot/Data.mining//iCAMP.16s//",prefix,".pNST.bootstr.",colnames(treat)[1],".csv"))
write.csv(pnst.bt$NST.compare,file = paste0("./result_and_plot/Data.mining//iCAMP.16s//",prefix,".pNST.compare.",colnames(treat)[1],".csv"))



i=2
treat.use = treat[,i,drop=FALSE] %>% as.data.frame()
# treat.use = data.frame( as.numeric(as.factor(map$Group))) %>%as.matrix()
row.names(treat.use) = treat$SampleID


icamp.result=icres$CbMPDiCBraya
icboot=iCAMP::icamp.boot(icamp.result = icamp.result,
                         treat = treat.use,
                         rand.time = rand.time,
                         compare = TRUE,
                         silent = FALSE,
                         between.group = TRUE,
                         ST.estimation = TRUE
                         )
save(icboot,file=paste0(prefix,".iCAMP.Boot.",colnames(treat)[i],".rda"))
write.csv(icboot$summary,file = paste0("./result_and_plot/Data.mining//iCAMP.16s/",prefix,".iCAMP.BootSummary.",colnames(treat)[i],".csv"),row.names = FALSE)
write.csv(icboot$compare,file = paste0("./result_and_plot/Data.mining//iCAMP.16s/",prefix,".iCAMP.Compare.",colnames(treat)[i],".csv"),row.names = FALSE)







# The qpen() function implements Stegen's community construction analysis 
library(iCAMP)

#Phylogenetic distance is obtained from evolutionary tree
pd <- cophenetic(tree) 


set.seed(123)
qpen.out <- qpen(comm = comm, pd = pd, sig.bNTI = 2, sig.rc = 0.95, rand.time = 100, nworker = 4)

qpen.out$ratio  #Heterogeneous.Selection, Homogeneous.Selection, Dispersal.Limitation, Homogenizing.Dispersal and Un dominated) sample pair dominated
head(qpen.out$result)  # bMNTD, BC, bNTI, RC values of each sample pair and the most important ecological processes


write.csv(qpen.out, './result_and_plot/Data.mining//iCAMP.16s/qpen.out.csv', row.names = FALSE)
write.csv(qpen.out, './result_and_plot/Data.mining//iCAMP.ITs/qpen.out.csv', row.names = FALSE)

library(tidyverse)
icamp=read.csv("E:/result_and_plot/Data.mining//iCAMP.16s/qpen.out.csv")
icamp

its=read.csv("E:/result_and_plot/Data.mining//iCAMP.ITs/qpen.out.csv")
its

env=read.csv("E:/data/env1.csv")
group1=env[,1:2]

group1$sample1 = group1$sample_id
head(group1$type)
#Treat
icamp %>%
  left_join(group1, by = "sample1") %>%
  left_join(group1, by = c("sample2" = "sample1")) %>%
  filter(type.x %in% "Grassland",type.y %in% "Grassland") ->icamp_Grassland

icamp %>%
  left_join(group1, by = "sample1") %>%
  left_join(group1, by = c("sample2" = "sample1")) %>%
  filter(type.x %in% "Desert",type.y %in% "Desert") ->icamp_Desert

icamp %>%
  left_join(group1, by = "sample1") %>%
  left_join(group1, by = c("sample2" = "sample1")) %>%
  filter(type.x %in% "Mountain",type.y %in% "Mountain") ->icamp_Mountain

d1 <- bind_rows(
  icamp_Grassland %>% mutate(ecosystem = "Grassland"),
  icamp_Desert %>% mutate(ecosystem = "Desert"),
  icamp_Mountain %>% mutate(ecosystem = "Mountain")
)

names(d1)

its %>%
  left_join(group1, by = "sample1") %>%
  left_join(group1, by = c("sample2" = "sample1")) %>%
  filter(type.x %in% "Grassland",type.y %in% "Grassland") ->its_Grassland

its %>%
  left_join(group1, by = "sample1") %>%
  left_join(group1, by = c("sample2" = "sample1")) %>%
  filter(type.x %in% "Desert",type.y %in% "Desert") ->its_Desert

its %>%
  left_join(group1, by = "sample1") %>%
  left_join(group1, by = c("sample2" = "sample1")) %>%
  filter(type.x %in% "Mountain",type.y %in% "Mountain") ->its_Mountain

d2 <- bind_rows(
  its_Grassland %>% mutate(ecosystem = "Grassland"),
  its_Desert %>% mutate(ecosystem = "Desert"),
  its_Mountain %>% mutate(ecosystem = "Mountain")
)


d3 <- bind_rows(
  d2 %>% mutate(ty = "Fungi"),
  d1 %>% mutate(ty = "Bacteria"))

names(d3)
v1 = c("#f85a40","#00c16e","#037ef3")
# Mapping of bacterial β-NTI distribution
p1 <- ggplot(d1, aes(ecosystem, bNTI, fill = ecosystem)) + 
  geom_violin(trim = FALSE, position = position_dodge(0.9), show.legend = FALSE) + 
  geom_boxplot(width = 0.2, position = position_dodge(0.9), show.legend = FALSE) + 
  geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "red") + # 添加两条虚线
  scale_fill_manual(values = v1)+
  labs(x = "", 
       y = "β-NTI")+
  facet_grid( ~"Bacteria", drop=TRUE,scale="free", space="free_x")+
  theme_bw(base_size = 15)

p2 <- ggplot(d2, aes(ecosystem, bNTI, fill = ecosystem)) + 
  geom_violin(trim = FALSE, position = position_dodge(0.9), show.legend = FALSE) + 
  geom_boxplot(width = 0.2, position = position_dodge(0.9), show.legend = FALSE) + 
  geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "red") + # 添加两条虚线
  scale_fill_manual(values = v1)+
  labs(x = "", 
       y = "β-NTI")+
  facet_grid( ~"Fungi", drop=TRUE,scale="free", space="free_x")+
  theme_bw(base_size = 15)

p3 <- ggplot(d3, aes(ecosystem, bNTI, fill = ty)) + 
  geom_violin(trim = FALSE, position = position_dodge(0.9), show.legend = FALSE) + 
  geom_boxplot(width = 0.2, position = position_dodge(0.9), show.legend = FALSE) + 
  geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "red") + # 添加两条虚线
  scale_fill_manual(values = v1)+
  labs(x = "", 
       y = "β-NTI")+
  #facet_grid( ~"Fungi", drop=TRUE,scale="free", space="free_x")+
  theme_bw(base_size = 15)
library(ggplot2)
library(ggpubr)

p <-ggplot(data = d3, aes(ecosystem, bNTI, fill = ty))+
  geom_violin(position = position_dodge(width = 1),scale = "width",lwd=0.8,alpha=0.85)+
  geom_boxplot(position = position_dodge(1),width=.4,lwd=0.8,alpha=1)+
  geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "red") + # 添加两条虚线
  #geom_point(shape=21,size=2,color = "black",alpha=0.6) +
             #position = position_jitterdodge(jitter.width = 0.8,
                                             #jitter.height = 0,
                                             #dodge.width = 1)) +
  stat_summary(fun.y = "median",geom = "point",shape = 16, 
               size = 2, color = "white",position =position_dodge(1)) +
  scale_fill_manual(values = c("#C7C7C7", "#F5BB5D"))+
  scale_color_manual(values = c("#C7C7C7", "#F5BB5D"))+
  #theme_classic()+
  theme_bw(base_size = 15)+
  theme(legend.position = c(0.9,0.7))+
  theme(panel.grid = element_blank())+
  theme(legend.title = element_blank())+#删除所有图例
  #theme(legend.direction = 'horizontal')+
  stat_compare_means(method = "wilcox.test", 
                     method.args = list(p.adjust.method = "bonferroni", stars = c("***", "**", "*", ".", "NS")),
                     aes(label = paste0("p : ", after_stat(p.signif))),
                     size = 5)+
  facet_wrap(.~ecosystem,scales = "free_x",nrow = 1)+
   labs(x="",y="β-NTI")+
  scale_y_continuous(expand = expansion(mult = c(0,0.1)))
p

ggplot(data = d3, aes(ecosystem, bNTI, fill = ty))+
  geom_violin(position = position_dodge(width = 1),scale = "width",lwd=0.8,alpha=0.85)+
  geom_boxplot(position = position_dodge(1),width=.4,lwd=0.8,alpha=1)+
  #, outlier.shape = NA,outlier.shape = 21)+
  geom_point(aes(fill=ty),shape=21,size=2,color = "black",alpha=0.6,
             position = position_jitterdodge(jitter.width = 0.8,
                                             jitter.height = 0,
                                             dodge.width = 1)) +
  stat_summary(fun.y = "median",geom = "point",shape = 16, 
               size = 2, color = "white",position =position_dodge(1)) +
  scale_fill_manual(values = c("#849871","#ddc5b2"))+
  scale_color_manual(values = c("#849871","#ddc5b2"))+
  theme_bw()+
  stat_compare_means(#comparisons = list( c("NR","R")),
    method = "wilcox.test", size=5,
    aes( label = paste0("p : ", after_stat(p.signif)),)
  )+
  facet_wrap(.~ecosystem,scales = "free",nrow = 2)+
  theme(strip.text.x = element_text(size=9,color="black"))+labs(x="",y="")+
  scale_y_continuous(expand = expansion(mult = c(0,0.1)))

icamp_Grassland %>%
  count(process) %>%
  mutate(proportion = n/528 *100) ->icamp_Grassland1

icamp_Desert %>%
  count(process) %>%
  mutate(proportion = n/741 *100) ->icamp_Desert1

icamp_Mountain %>%
  count(process) %>%
  mutate(proportion = n/528 *100) ->icamp_Mountain1

its_Grassland %>%
  count(process) %>%
  mutate(proportion = n/528 *100) ->its_Grassland1

its_Desert %>%
  count(process) %>%
  mutate(proportion = n/741 *100) ->its_Desert1

its_Mountain %>%
  count(process) %>%
  mutate(proportion = n/528 *100) ->its_Mountain1
library(dplyr)

# Combined data
df1 <- bind_rows(
  icamp_Grassland1 %>% mutate(ecosystem = "Grassland"),
  icamp_Desert1 %>% mutate(ecosystem = "Desert"),
  icamp_Mountain1 %>% mutate(ecosystem = "Mountain")
)

df2 <- bind_rows(
  its_Grassland1 %>% mutate(ecosystem = "Grassland"),
  its_Desert1 %>% mutate(ecosystem = "Desert"),
  its_Mountain1 %>% mutate(ecosystem = "Mountain")
)

names(df1)

library(ggalluvial)
v2 <- c("#D7DBE3", "#88B1DB", "#F5DF1E", "#87D5D6", "#A86634")

v2 <-c("#5994CF", "#3BB0EB", "#F2DBA4", "#A86634")
v3 <-c("#5994CF", "#3BB0EB", "#F2DBA4", "#F08B27","#A86634")
p1 <- ggplot(df1, aes(x =ecosystem, y = proportion, 
                fill =process,stratum =process, alluvium =process)) +
  geom_stratum(width = 0.8,alpha=1,color="white") +geom_flow(alpha = 0.01) + 
  scale_fill_manual(values=v2)+
  theme(axis.text=element_text(colour='black',size=9))+
  labs(x = '', y = 'Proportion',fill=" ")+theme_bw()+
  #theme(legend.position = "none")+
  #theme(axis.text=element_text(colour='black',size=11))+
  #theme(strip.text = element_text(size = 12))+
  facet_grid(~"Bacteria", drop=TRUE,scale="free",space="free_x")+
  theme_bw(base_size = 15)+
  theme(legend.position = "none")
p1
#"Fungi", 
p2 <- ggplot(df2, aes(x =ecosystem, y = proportion, 
                fill =process,stratum =process, alluvium =process)) +
  geom_stratum(width = 0.8,alpha=1,color="white") +geom_flow(alpha = 0.01) + 
  scale_fill_manual(values=v3)+
  theme(axis.text=element_text(colour='black',size=9))+
  labs(x = '', y = '',fill=" ")+theme_bw()+
  #theme(axis.text=element_text(colour='black',size=11))+
  #theme(legend.position = c(0.95,0.6))+
  #theme(legend.position = 'bottom')+
  #theme(legend.direction = 'horizontal')+
  #theme(strip.text = element_text(size = 12))+
  facet_grid(~"Fungi", drop=TRUE,scale="free",space="free_x")+
  theme_bw(base_size = 15)
p2


library(ggpubr)
com <- ggarrange(
  p, 
  ggarrange(p1, p2, nrow = 1, ncol = 2, labels = c("B"),
            widths = c(0.85, 1.55),
            heights = c(1, 1)), 
  nrow = 2,
  labels = c("A", ""),
  widths = c(2, 2),
  heights = c(1, 1)
)
com
library(export)
graph2ppt(com,file = "cor_mantel.ppt", width = 13.5, height = 5.6 )
ggsave("cor_mantel.pdf",com, device = "pdf",family = "Times",width = 13.5,height = 5.6)


#
icamp_Grassland %>%
  mutate(total_Heterogeneous = sum(Heterogeneous.Selection)/n(),
         total_Homogeneous = sum(Homogeneous.Selection)/n(),
         total_Dispersal = sum(Dispersal.Limitation)/n(),
         total_Homogenizing = sum(Homogenizing.Dispersal)/n(),
         total_Drift = sum(Drift.and.Others)/n()) -> icamp_Grassland

icamp_Desert %>%
  mutate(total_Heterogeneous = sum(Heterogeneous.Selection)/n(),
         total_Homogeneous = sum(Homogeneous.Selection)/n(),
         total_Dispersal = sum(Dispersal.Limitation)/n(),
         total_Homogenizing = sum(Homogenizing.Dispersal)/n(),
         total_Drift = sum(Drift.and.Others)/n()) -> icamp_Desert

icamp_Mountain %>%
  mutate(total_Heterogeneous = sum(Heterogeneous.Selection)/n(),
         total_Homogeneous = sum(Homogeneous.Selection)/n(),
         total_Dispersal = sum(Dispersal.Limitation)/n(),
         total_Homogenizing = sum(Homogenizing.Dispersal)/n(),
         total_Drift = sum(Drift.and.Others)/n()) -> icamp_Mountain
#输出
write.csv(icamp_Grassland, './result_and_plot/Data.mining//iCAMP.16s/qpen.out.csv', row.names = FALSE)


