otu = read.csv("E:/要完成的文章/第二篇文章/Life-history/16s_1.csv",  row.names=1)
otu = read.csv("E:/要完成的文章/第二篇文章/Life-history/ITs_1.csv",  row.names=1)

otu = read.csv("E:/要完成的文章/第二篇文章/Life-history/16s_2.csv",  row.names=1)
otu = read.csv("E:/要完成的文章/第二篇文章/Life-history/ITs_2.csv",  row.names=1)

otu = read.csv("E:/要完成的文章/第二篇文章/Life-history/16s_3.csv",  row.names=1)
otu = read.csv("E:/要完成的文章/第二篇文章/Life-history/ITs_3.csv",  row.names=1)

otu = read.csv("E:/要完成的文章/第二篇文章/Life-history/BF_g1.csv",  row.names=1)
otu = read.csv("E:/要完成的文章/第二篇文章/Life-history/BF_g2.csv",  row.names=1)
otu = read.csv("E:/要完成的文章/第二篇文章/Life-history/BF_g3.csv",  row.names=1)

library(igraph)

#输入数据，邻接矩阵

#注：最好把其中的孤立点删除（也就是不存在和其它任何节点具有显著相关性的节点）。我这个示例中没有剔除它们，某些分析中可能会产生干扰。

neetwork_adj <- read.delim('./20240326/16s_sparcc/neetwork1.adj.txt', row.names = 1, sep = '\t', check.names = FALSE)
neetwork_adj <- read.delim('./20240326/16s_sparcc/neetwork2.adj.txt', row.names = 1, sep = '\t', check.names = FALSE)
neetwork_adj <- read.delim('./20240326/16s_sparcc/neetwork3.adj.txt', row.names = 1, sep = '\t', check.names = FALSE)

neetwork_adj <- read.delim('./20240326/ITs_sparcc/neetwork1.adj.txt', row.names = 1, sep = '\t', check.names = FALSE)
neetwork_adj <- read.delim('./20240326/ITs_sparcc/neetwork2.adj.txt', row.names = 1, sep = '\t', check.names = FALSE)
neetwork_adj <- read.delim('./20240326/ITs_sparcc/neetwork3.adj.txt', row.names = 1, sep = '\t', check.names = FALSE)

neetwork_adj <- read.delim('./20240326/BF_sparcc/neetwork1.adj.txt', row.names = 1, sep = '\t', check.names = FALSE)
neetwork_adj <- read.delim('./20240326/BF_sparcc/neetwork2.adj.txt', row.names = 1, sep = '\t', check.names = FALSE)
neetwork_adj <- read.delim('./20240326/BF_sparcc/neetwork3.adj.txt', row.names = 1, sep = '\t', check.names = FALSE)

head(neetwork_adj)[1:6]    #邻接矩阵类型的网络文件

#neetwork_adj <- neetwork_adj[rowSums(neetwork_adj) != 0, colSums(neetwork_adj) != 0]

#邻接矩阵 -> igraph 的邻接列表，获得含权的无向网络
g <- graph_from_adjacency_matrix(as.matrix(neetwork_adj), mode = 'undirected', weighted = TRUE, diag = FALSE)
g    #igraph 的邻接列表


#并将所有子网络存储到一个列表（list）里面
sub_graph <- list()
for (i in names(otu)) {
  sample_i <- otu[i]
  select_node <- rownames(sample_i)[which(sample_i != 0)]
  sub_graph[[i]] <- subgraph(g, select_node)
}
sub_graph

# 初始化向量
sample_name <- c()
nodes_num <- c()
degree <- c()
average_path_length <- c()
betweenness_centralization <- c()
linkage_density <- c()
connectance <- c()
average_neighborhood <- c()

for(i in 1:length(sub_graph)){
  sample_name <- c(sample_name, names(sub_graph[i]))
  nodes <- V(sub_graph[[i]])
  
  # 节点数量
  nodes_num <- c(nodes_num, length(nodes))
  
  # 平均度
  degree <- c(degree, mean(degree(sub_graph[[i]])))
  
  # 平均路径长度
  # 检查权重是否为非负
  if(any(E(sub_graph[[i]])$weight < 0)){
    # 如果存在负权重，使用绝对值作为权重
    abs_weights <- abs(E(sub_graph[[i]])$weight)
    average_path_length <- c(average_path_length, average.path.length(sub_graph[[i]], weights = abs_weights, directed = FALSE))
  } else {
    average_path_length <- c(average_path_length, average.path.length(sub_graph[[i]], directed = FALSE))
  }
  
  # 介数中心性
  betweenness_centralization <- c(betweenness_centralization, centralization.betweenness(sub_graph[[i]])$centralization)
  
  # 连接密度 (linkage density)
  e <- sum(sapply(nodes, function(n) sum(degree(sub_graph[[i]], n))/2)) / (length(nodes)*(length(nodes)-1)/2)
  linkage_density <- c(linkage_density, e)
  
  # 连接率 (connectance)
  f <- mean(degree(sub_graph[[i]])) / (length(nodes)-1)
  connectance <- c(connectance, f)
  
  # 平均邻居度 (average neighborhood)
  # 检查每个节点是否有邻居，邻居度是否大于零
  g_values <- sapply(nodes, function(n) {
    neighbors_degrees <- degree(sub_graph[[i]], neighbors(sub_graph[[i]], n))
    if(length(neighbors_degrees) > 0){
      mean(neighbors_degrees)
    } else {
      0
    }
  })
  
  average_neighborhood <- c(average_neighborhood, mean(g_values))
}


# 创建数据框
sub_graph_stat <- data.frame(
  nodes_num,
  degree,
  average_path_length,
  betweenness_centralization,
  linkage_density,
  connectance,
  average_neighborhood
)

# 设置行名
rownames(sub_graph_stat) <- sample_name

# 显示前几行数据
head(sub_graph_stat)

write.csv(sub_graph_stat, './20240326/16s_sparcc/Network_complexity1.csv', quote = FALSE)
write.csv(sub_graph_stat, './20240326/16s_sparcc/Network_complexity2.csv', quote = FALSE)
write.csv(sub_graph_stat, './20240326/16s_sparcc/Network_complexity3.csv', quote = FALSE)

write.csv(sub_graph_stat, './20240326/ITs_sparcc/Network_complexity1.csv', quote = FALSE)
write.csv(sub_graph_stat, './20240326/ITs_sparcc/Network_complexity2.csv', quote = FALSE)
write.csv(sub_graph_stat, './20240326/ITs_sparcc/Network_complexity3.csv', quote = FALSE)

write.csv(sub_graph_stat, './20240326/BF_sparcc/Network_complexity1.csv', quote = FALSE)
write.csv(sub_graph_stat, './20240326/BF_sparcc/Network_complexity2.csv', quote = FALSE)
write.csv(sub_graph_stat, './20240326/BF_sparcc/Network_complexity3.csv', quote = FALSE)


##计算各个子网络的拓扑指数，随便选了 4 个指数为例
sample_name <- c()
nodes_num <- c()
degree <- c()
Positive_edge <- c()
Negative_edge <- c()
#average_path_length <- c()
betweenness_centralization <- c()
library(tidyverse)
for(i in 1:length(sub_graph)){
  sample_name <- c(sample_name, names(sub_graph[i]))
  nodes_num <- c(nodes_num, length(V(sub_graph[[i]])))  #节点数量
  degree <- c(degree, mean(degree(sub_graph[[i]])))  #平均度
  Positive_edge <- c(Positive_edge, nrow(get.edge.attribute(sub_graph[[i]])[[1]] %>% 
                                           as_tibble() %>% 
                                           dplyr::filter(value > 0)))
  Negative_edge <- c(Negative_edge, nrow(get.edge.attribute(sub_graph[[i]])[[1]] %>% 
                                           as_tibble() %>% 
                                           dplyr::filter(value < 0)))
  #average_path_length <- c(average_path_length, average.path.length(sub_graph[[i]], directed = FALSE))  #平均路径长度
  betweenness_centralization <- c(betweenness_centralization, centralization.betweenness(sub_graph[[i]])$centralization)  #介数中心性
}

sub_graph_stat <- data.frame(nodes_num, degree, Positive_edge, Negative_edge,  betweenness_centralization)
rownames(sub_graph_stat) <- sample_name
head(sub_graph_stat)  #与各样本有关的“子网络”的拓扑指数，本示例以节点数量、平均度、平均路径长度、介数中心性为例

#输出统计表
write.csv(sub_graph_stat, './20240326/16s_sparcc/16s1_subnet.csv', quote = FALSE)
write.csv(sub_graph_stat, './20240326/16s_sparcc/16s2_subnet.csv', quote = FALSE)
write.csv(sub_graph_stat, './20240326/16s_sparcc/16s3_subnet.csv', quote = FALSE)

write.csv(sub_graph_stat, './20240326/ITs_sparcc/ITs1_subnet.csv', quote = FALSE)
write.csv(sub_graph_stat, './20240326/ITs_sparcc/ITs2_subnet.csv', quote = FALSE)
write.csv(sub_graph_stat, './20240326/ITs_sparcc/ITs3_subnet.csv', quote = FALSE)

write.csv(sub_graph_stat, './20240326/BF_sparcc/BFs1_subnet.csv', quote = FALSE)
write.csv(sub_graph_stat, './20240326/BF_sparcc/BFs2_subnet.csv', quote = FALSE)
write.csv(sub_graph_stat, './20240326/BF_sparcc/BFs3_subnet.csv', quote = FALSE)
