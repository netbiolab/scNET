#!/usr/bin/env Rscript
#usage Rscript.R -bs [net1 BS] -sv [net2 SV] -mc [net SC] 

#this code takes three network (bin selected from BS, SV, SC) and check overlap of the two network


# Seungbyn Baek mod : 2021.05.11
# Rscript.R [BS_fit] [SV_fit] [outdir] [filename]


suppressPackageStartupMessages(library(argparse))
## input files
parser <- ArgumentParser()
parser$add_argument("-b", "--BSnet", help="cut net from BS")
parser$add_argument("-s", "--SVnet", help="cut net from SV")
parser$add_argument("-c", "--SCnet", help="cut net from SC")

parser$add_argument("-d", "--outdir", help="saving directory")
parser$add_argument("-o", "--output", help="filename, current directory will be added to prefix")

args <- parser$parse_args()

library(igraph)
library(ggraph)
library(graphlayouts)
library(VennDiagram)
library(extrafont)
library(gridExtra)
library(Vennerable)


path.to.BS <- args$BSnet
path.to.SV <- args$SVnet
path.to.SC <- args$SCnet
out_dir <- args$outdir
output.file <- args$output

setwd(out_dir)


#remove negative link if it exists
edgelist.bs <- read.table(path.to.BS, sep = '\t', header=F)
edgelist.bs <- edgelist.bs[edgelist.bs$V3 > 0,]

edgelist.saver <- read.table(path.to.SV, sep='\t', header=F)
edgelist.saver <- edgelist.saver[edgelist.saver$V3 > 0,]

edgelist.sc <- read.table(path.to.SC, sep='\t', header=F)
edgelist.sc <- edgelist.sc[edgelist.sc$V3 > 0,]


#sort the networks for edge overlap
sort_edgelist <- function(edgelist){
  tmp <- t(apply(edgelist[,1:2],1,sort))
  x <- cbind(tmp, edgelist[,3])
  sorted.net <- as.data.frame(x)
  sorted.net[,3] <- as.numeric(sorted.net[,3])
  
  return(sorted.net)
}

edgelist.bs <- sort_edgelist(edgelist.bs)
edgelist.saver <- sort_edgelist(edgelist.saver) 
edgelist.sc <- sort_edgelist(edgelist.sc)


#make into igraph format
sc_net_saver <- graph.data.frame(edgelist.saver, directed = FALSE)
sc_net_bs <- graph.data.frame(edgelist.bs, directed=F)
sc_net_sc <- graph.data.frame(edgelist.sc, directed=F)



####plot 
# pdf(paste0('./',output.file,'_overall_net_topo_SV.pdf'))
# ggraph::autograph(sc_net_saver)
# dev.off()

# pdf(paste0('./',output.file,'_overall_net_topo_BS.pdf'))
# ggraph::autograph(sc_net_bs)
# dev.off()

# define a custom color palette
got_palette <- c("#1A5878", "#C44237", "#AD8941", "#E99093", "#50594B")

cat('drawing networks...\n')

setwd(out_dir)
pdf(paste0(output.file,'network_topo.pdf'), width = 10, height=10)
LLS <- E(sc_net_saver)$V3
ggraph(sc_net_saver,layout = "stress")+
  geom_edge_link(aes(edge_width = LLS),edge_colour = "grey66")+
  geom_node_point(color = "#1A5878")+
  scale_edge_width(range = c(0.2,3))+
  scale_size(range = c(1,6))+
  #theme_graph(base_family = "FreeSans",title_size = 18)+
  theme(legend.position = "bottom")+
  ggtitle('Saver Network from Benchmark')

LLS <- E(sc_net_bs)$V3
ggraph(sc_net_bs, layout = "stress")+
  geom_edge_link0(aes(edge_width = LLS),edge_colour = "grey66")+
  geom_node_point(color = "#C44237")+
  scale_edge_width(range = c(0.2,3))+
  scale_size(range = c(1,6))+
  #theme_graph(base_family = "FreeSans", title_size = 18)+
  theme(legend.position = "bottom")+
  ggtitle('BS Network from Benchmark')

LLS <- E(sc_net_sc)$V3
ggraph(sc_net_sc, layout = "stress")+
  geom_edge_link0(aes(edge_width = LLS),edge_colour = "grey66")+
  geom_node_point(color = "#0CB702")+
  scale_edge_width(range = c(0.2,3))+
  scale_size(range = c(1,6))+
  #theme_graph(base_family = "FreeSans", title_size = 18)+
  theme(legend.position = "bottom")+
  ggtitle('SC Network from Benchmark')
dev.off()



nodes.bs <- V(sc_net_bs)$name
nodes.saver <- V(sc_net_saver)$name
nodes.sc <- V(sc_net_sc)$name


#see overlap
set1 <- nodes.saver
set2 <- nodes.bs
set3 <- nodes.sc

nodes.list <- list(set1,set2,set3)
names(nodes.list) <- c('SV', 'BS', 'SC')
Vnodes <- Venn(nodes.list)


#get overlap of edges
set4 <- as_ids(E(sc_net_saver))
set5 <- as_ids(E(sc_net_bs))
set6 <- as_ids(E(sc_net_sc))

edges.list <- list(set4,set5,set6)
names(edges.list) <- c('SV','BS','SC')
Vedges <- Venn(edges.list)


cat('Drawing networks....\n')
pdf(paste0(output.file,'_VennDiagram_NetworkOverlap.pdf'))
plot(Vnodes, doWeights = TRUE)
plot(Vedges, doWeights = TRUE)
###############################################
#inter <- E(intersection(intersection(sc_net_bs,sc_net_saver),sc_net_sc)) #use intersection not intersect!!!


#what does integrated graph look like?
net_all <- union(union(sc_net_bs, sc_net_saver), sc_net_sc)
#net_inter <- intersection(sc_net_bs, sc_net_saver, sc_net_sc)

V(sc_net_saver)$platform = "SV"
V(sc_net_bs)$platform = "BS"
V(sc_net_sc)$platform = "SC"

net_all <- graph.union(sc_net_bs, sc_net_saver, sc_net_sc, byname=TRUE)

V(net_all)$platform <- ifelse(V(net_all)$name %in% V(sc_net_bs)$name & V(net_all)$name %in% V(sc_net_saver)$name & V(net_all)$name %in% V(sc_net_sc)$name, 'union',
                           ifelse(V(net_all)$name %in% V(sc_net_bs)$name, 'BS',
                                  ifelse(V(net_all)$name %in% V(sc_net_saver)$name, 'SV', 
                                         ifelse(V(net_all)$name %in% V(sc_net_sc)$name, 'SC', NA))))



ggraph(net_all,layout = "stress")+
  geom_edge_link0(edge_colour = "grey66")+
  geom_node_point(shape=21,aes(fill = V(net_all)$platform))+
  scale_edge_width(range = c(0.2,3))+
  scale_size(range = c(1,6))+
  #theme_graph(base_family = "FreeSans", title_size = 18)+
  theme(legend.position = "bottom") +
  ggtitle('Combined Network') +
  guides(fill=guide_legend(title="Network Methods"))


dev.off()



#top 50 degree / btw nodes overlap
degree.bs <- igraph::strength(sc_net_bs, weights = E(sc_net_bs)$V3)
between.bs <- igraph::betweenness(sc_net_bs) # no weight

deg.bs.top100 <- degree.bs[order(degree.bs, decreasing = T)][1:100]
length(grep('^RP', names(deg.bs.top100)))
cat(paste('ribosome related:',length(grep('^RP', names(deg.bs.top100))), 'genes out of 100\n'))
write.table(names(deg.bs.top100), file=paste0(output.file,'_deg_BS_top100'), row.names = F, col.names = F, quote = F)

btw.bs.top100 <- between.bs[order(between.bs, decreasing = T)][1:100]
cat(paste('ribosome related:',length(grep('^RP', names(btw.bs.top100))), 'genes out of 100\n'))
write.table(names(deg.bs.top100), file=paste0(output.file,'_btw_BS_top100'), row.names = F, col.names = F, quote = F)


degree.sv <- igraph::strength(sc_net_saver, weights = E(sc_net_saver)$V3)
between.sv <- igraph::betweenness(sc_net_saver)
deg.sv.top100 <- degree.sv[order(degree.sv, decreasing = T)][1:100]
write.table(names(deg.bs.top100), file=paste0(output.file,'_deg_SV_top100'), row.names = F, col.names = F, quote = F)
btw.sv.top100 <- between.sv[order(between.sv, decreasing = T)][1:100]
write.table(names(deg.bs.top100), file=paste0(output.file,'_btw_SV_top100'), row.names = F, col.names = F, quote = F)


degree.sc <- igraph::strength(sc_net_sc, weights = E(sc_net_sc)$V3)
between.sc <- igraph::betweenness(sc_net_sc)
deg.sc.top100 <- degree.sc[order(degree.sc, decreasing = T)][1:100]
write.table(names(deg.sc.top100), file=paste0(output.file,'_deg_SC_top100'), row.names = F, col.names = F, quote = F)
btw.sc.top100 <- between.sc[order(between.sc, decreasing = T)][1:100]
write.table(names(deg.sc.top100), file=paste0(output.file,'_btw_SC_top100'), row.names = F, col.names = F, quote = F)


#see overlap
set7 <- deg.sv.top100
set8 <- deg.bs.top100
set9 <- deg.sc.top100


deg.list <- list(names(set7),names(set8),names(set9))
names(deg.list) <- c('SV','BS','SC')
Vdeg <- Venn(deg.list)

set10 <- btw.sv.top100
set11 <- btw.bs.top100
set12 <- btw.sc.top100


btw.list <- list(names(set10), names(set11), names(set12))
names(btw.list) <- c('SV','BS','SC')
Vbtw <- Venn(btw.list)

###############################################
pdf(paste0(output.file,'_VennDiagram_deg_btw_top100.pdf'))
plot(Vdeg, doWeights = TRUE)
plot(Vbtw, doWeights = TRUE)
dev.off()




