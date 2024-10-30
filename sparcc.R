#SpiecEasi ：https://github.com/zdk123/SpiecEasi
#install.packages('devtools')
#devtools::install_github('zdk123/SpiecEasi')

library(SpiecEasi)

setwd("E:/要完成的文章/第二篇文章/Life-history")

path <- "./20240326/16s_sparcc/"
path <- "./20240326/ITs_sparcc/"
path <- "./20240326/BF_sparcc/"
dir.create(path, recursive = TRUE)

amgut1.filt = read.csv("E:/要完成的文章/第二篇文章/Life-history/16s_1.csv",  row.names=1)
amgut1.filt = read.csv("E:/要完成的文章/第二篇文章/Life-history/16s_2.csv",  row.names=1)
amgut1.filt = read.csv("E:/要完成的文章/第二篇文章/Life-history/16s_3.csv",  row.names=1)

amgut1.filt = read.csv("E:/要完成的文章/第二篇文章/Life-history/ITs_1.csv",  row.names=1)
amgut1.filt = read.csv("E:/要完成的文章/第二篇文章/Life-history/ITs_2.csv",  row.names=1)
amgut1.filt = read.csv("E:/要完成的文章/第二篇文章/Life-history/ITs_3.csv",  row.names=1)

amgut1.filt = read.csv("E:/要完成的文章/第二篇文章/Life-history/BF_g1.csv",  row.names=1)

amgut1.filt <- data.frame(t(amgut1.filt))
amgut1.filt[1:6,1:6]

library(doParallel)

# Set the number of cores to use
cores <- 6  


registerDoParallel(cores)
#Perform sparcc analysis
set.seed(123)
amgut1.filt.sparcc <- sparcc(amgut1.filt, iter = 20, inner_iter = 10,th = 0.1)
sparcc0 <- amgut1.filt.sparcc$Cor  #Sparse correlation matrix



#Sparse correlation matrix of output observations
colnames(sparcc0) <- colnames(amgut1.filt)
rownames(sparcc0) <- colnames(amgut1.filt)

#tablename <- paste(path,'sparcc0.txt', sep = '\t', col.names = NA, quote = FALSE)
#write.csv(data,tablename)
#filename <- file.path(path, "sparcc0.txt")


write.table(sparcc0, './20240326/16s_sparcc/sparcc1.txt', sep = '\t', col.names = NA, quote = FALSE)
write.table(sparcc0, './20240326/ITs_sparcc/sparcc1.txt', sep = '\t', col.names = NA, quote = FALSE)
write.table(sparcc0, './20240326/ITs_sparcc/sparcc2.txt', sep = '\t', col.names = NA, quote = FALSE)
write.table(sparcc0, './20240326/ITs_sparcc/sparcc3.txt', sep = '\t', col.names = NA, quote = FALSE)

#The random matrix was obtained by 100 bootstrap sampling
set.seed(123)
n = 100
#bacteria
for (i in 1:n) {
  amgut1.filt.boot <- sample(amgut1.filt, replace = TRUE)  #bootstrap
  amgut1.filt.sparcc_boot <- sparcc(amgut1.filt.boot, iter = 20, inner_iter = 10, th = 0.1)  #sparcc 参数设置和上文保持一致
  sparcc_boot <- amgut1.filt.sparcc_boot$Cor
  colnames(sparcc_boot) <- colnames(amgut1.filt.boot)
  rownames(sparcc_boot) <- colnames(amgut1.filt.boot)
  write.table(sparcc_boot, paste('./20240326/16s_sparcc/sparcc_boot_1', i, '.txt', sep = ''), sep = '\t', col.names = NA, quote = FALSE)  #输出随机值的稀疏相关性矩阵
}
#fungus
for (i in 1:n) {
  amgut1.filt.boot <- sample(amgut1.filt, replace = TRUE)  #bootstrap
  amgut1.filt.sparcc_boot <- sparcc(amgut1.filt.boot, iter = 20, inner_iter = 10, th = 0.1)  #sparcc 参数设置和上文保持一致
  sparcc_boot <- amgut1.filt.sparcc_boot$Cor
  colnames(sparcc_boot) <- colnames(amgut1.filt.boot)
  rownames(sparcc_boot) <- colnames(amgut1.filt.boot)
  write.table(sparcc_boot, paste('./20240326/ITs_sparcc/sparcc_boot_3', i, '.txt', sep = ''), sep = '\t', col.names = NA, quote = FALSE)  #输出随机值的稀疏相关性矩阵
}
#Based on the sparse correlation matrix of the above observed values and the results of bootstrap 100 times, the pseudo-P value of sparse correlation is obtained
p <- sparcc0
p[p!=0] <- 0

for (i in 1:n) {
  p_boot <- read.delim(paste('./20240325/16s_sparcc/sparcc_boot', i, '.txt', sep = ''), sep = '\t', row.names = 1)
  p[abs(p_boot)>=abs(sparcc0)] <- p[abs(p_boot)>=abs(sparcc0)] + 1
}

for (i in 1:n) {
  p_boot <- read.delim(paste('./20240326/ITs_sparcc/sparcc_boot_3', i, '.txt', sep = ''), sep = '\t', row.names = 1)
  p[abs(p_boot)>=abs(sparcc0)] <- p[abs(p_boot)>=abs(sparcc0)] + 1
} 

p <- p / n
write.table(p, './20240325/16s_sparcc/pvals.two_sided.txt', sep = '\t', col.names = NA, quote = FALSE)
write.table(p, './20240326/ITs_sparcc/pvals.two_sided1.txt', sep = '\t', col.names = NA, quote = FALSE)
write.table(p, './20240326/ITs_sparcc/pvals.two_sided2.txt', sep = '\t', col.names = NA, quote = FALSE)
write.table(p, './20240326/ITs_sparcc/pvals.two_sided3.txt', sep = '\t', col.names = NA, quote = FALSE)

#Correlation matrix of observations
cor_sparcc <- read.delim('./20240325/16s_sparcc/sparcc0.txt', row.names = 1, sep = '\t', check.names = FALSE)
cor_sparcc <- read.delim('./20240326/ITs_sparcc/sparcc1.txt', row.names = 1, sep = '\t', check.names = FALSE)
cor_sparcc <- read.delim('./20240326/ITs_sparcc/sparcc3.txt', row.names = 1, sep = '\t', check.names = FALSE)

#Pseudo-p-valued matrix
pvals <- read.delim('./20240325/16s_sparcc/pvals.two_sided.txt', row.names = 1, sep = '\t', check.names = FALSE)
pvals <- read.delim('./20240326/ITs_sparcc/pvals.two_sided1.txt', row.names = 1, sep = '\t', check.names = FALSE)
pvals <- read.delim('./20240326/ITs_sparcc/pvals.two_sided3.txt', row.names = 1, sep = '\t', check.names = FALSE)

#Retain the value of | correlation |>0.3 and p<0.05 (all correlations that do not meet the conditions are assigned 0)
cor_sparcc[abs(cor_sparcc)<=0.3 | pvals>=0.05] <- 0

#Converts the value in the diagonal of the correlation matrix (which represents the autocorrelation) to 0
diag(cor_sparcc) <- 0



#Output a network file of the adjacency matrix type
write.table(cor_sparcc, './20240325/16s_sparcc/neetwork.adj.txt', col.names = NA, sep = '\t', quote = FALSE)
write.table(cor_sparcc, './20240326/ITs_sparcc/neetwork1.adj.txt', col.names = NA, sep = '\t', quote = FALSE)
write.table(cor_sparcc, './20240326/ITs_sparcc/neetwork3.adj.txt', col.names = NA, sep = '\t', quote = FALSE)


library(igraph)


neetwork_adj <- read.delim('./20240325/16s_sparcc/neetwork.adj.txt', row.names = 1, sep = '\t', check.names = FALSE)
neetwork_adj <- read.delim('./20240326/ITs_sparcc/neetwork1.adj.txt', row.names = 1, sep = '\t', check.names = FALSE)
neetwork_adj <- read.delim('./20240326/ITs_sparcc/neetwork3.adj.txt', row.names = 1, sep = '\t', check.names = FALSE)

head(neetwork_adj)[1:6]   


g <- graph_from_adjacency_matrix(as.matrix(neetwork_adj), mode = 'undirected', weighted = TRUE, diag = FALSE)
g    #igraph 


E(g)$sparcc <- E(g)$weight
E(g)$weight <- abs(E(g)$weight)


edge <- data.frame(as_edgelist(g))

edge_list <- data.frame(
  source = edge[[1]],
  target = edge[[2]],
  weight = E(g)$weight,
  sparcc = E(g)$sparcc
)
head(edge_list)

write.table(edge_list, './20240325/16s_sparcc/network.edge_list1.txt', sep = '\t', row.names = FALSE, quote = FALSE)

write.table(edge_list, './20240326/ITs_sparcc/network.edge_list3.txt', sep = '\t', row.names = FALSE, quote = FALSE)

#A list of node attributes, corresponding to an edge list, records node attributes
node_list <- data.frame(
  id = V(g)$name,
  label = V(g)$name,
  degree = degree(g)
)
head(node_list)

write.table(node_list, './20240325/16s_sparcc/network.node_list1.txt', sep = '\t', row.names = FALSE, quote = FALSE)

write.table(node_list, './20240326/ITs_sparcc/network.node_list3.txt', sep = '\t', row.names = FALSE, quote = FALSE)

d1 <- read.delim('./20240326/16s_sparcc/network.edge_list1.txt',  sep = '\t', check.names = FALSE)
d1 <- read.delim('./20240326/16s_sparcc/network.edge_list2.txt',  sep = '\t', check.names = FALSE)
d1 <- read.delim('./20240326/16s_sparcc/network.edge_list3.txt',  sep = '\t', check.names = FALSE)

d2 <- read.delim('./20240326/ITs_sparcc/network.edge_list1.txt',  sep = '\t', check.names = FALSE)
d2 <- read.delim('./20240326/ITs_sparcc/network.edge_list2.txt',  sep = '\t', check.names = FALSE)
d2 <- read.delim('./20240326/ITs_sparcc/network.edge_list3.txt',  sep = '\t', check.names = FALSE)

d3 <- read.delim('./20240326/BF_sparcc/network.edge_list1.txt',  sep = '\t', check.names = FALSE)
d3 <- read.delim('./20240326/BF_sparcc/network.edge_list2.txt',  sep = '\t', check.names = FALSE)
d3 <- read.delim('./20240326/BF_sparcc/network.edge_list3.txt',  sep = '\t', check.names = FALSE)

names(d2)
d2$cor <- ifelse(d2$sparcc >= 0, "+", "-")

write.csv(d1, "./20240326/16s_sparcc/network.edge_list1.csv", row.names = F)
write.csv(d1, "./20240326/16s_sparcc/network.edge_list2.csv", row.names = F)
write.csv(d1, "./20240326/16s_sparcc/network.edge_list3.csv", row.names = F)

write.csv(d2, "./20240326/ITs_sparcc/network.edge_list1.csv", row.names = F)
write.csv(d2, "./20240326/ITs_sparcc/network.edge_list2.csv", row.names = F)
write.csv(d2, "./20240326/ITs_sparcc/network.edge_list3.csv", row.names = F)

write.csv(d3, "./20240326/BF_sparcc/network.edge_list1.csv", row.names = F)
write.csv(d3, "./20240326/BF_sparcc/network.edge_list2.csv", row.names = F)
write.csv(d3, "./20240326/BF_sparcc/network.edge_list3.csv", row.names = F)