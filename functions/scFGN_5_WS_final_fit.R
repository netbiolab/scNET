#!/usr/bin/env Rscript

#this code draws Weighted SUM network, run regression.py, fit to curve and produce final integrated network
#make a new column for the integrated column with LLS score
#if edge exist in both network, take the largest LLS of the two as integrated link
#check first if the LLS difference is big

suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser(description = 'this code takes WSnet, fit LLS and make the final integrated network')
parser$add_argument("-i", "--input", help="WSnet with BS SV SC")
parser$add_argument("-o", "--output", help="filename, current directory will be added to prefix")
parser$add_argument("-d", "--outdir", help="saved directory")
parser$add_argument("-p", "--package", help="provide path to where scNet package is")
parser$add_argument("-gs", "--goldstandard", default='example/input/gold_standard_symbol_HNv3',help = "gold standard to evaluate LLS")

args <- parser$parse_args()


library(drc)
library(igraph)


out_dir <- args$outdir
output.file <- args$output
input.path <- args$input
path.to.package <- args$package
path.to.gs <- args$goldstandard

setwd(out_dir)
ws <- read.table(input.path, sep='\t')

ws.net <- ws[,c(1,2,ncol(ws))]
colnames(ws.net) <- c('gene1','gene2',"WS")

write.table(ws.net,paste0(output.file,'_WS_net_final'), quote=F, row.names = F, col.names = F, sep='\t')

#run regredssion.py on WS_net

system(paste0('python3 ',path.to.package,'/functions/LLS.py ', out_dir,output.file, '_WS_net_final ',path.to.package,path.to.gs))

lls <- read.table(paste0(output.file,'_WS_net_final.binlls'), header=T)

#remove negative and fit
if( tail(lls$BinLLS,1)<0){
    lls <- lls[-nrow(lls),]
}

model <- drm(BinLLS ~ MeanBinStatistics, fct = L.4(), data = lls) #parameter 4 log logistic model

lls.fit <- predict(model, data.frame(ws.net$WS))
cat(paste('max fitted LLS value: ', print(max(lls.fit))))

ws.net$LLS_final <- lls.fit

ws.net$WS <- NULL
ws.net <- ws.net[order(ws.net$LLS_final, decreasing=T), ]

ws.net <- ws.net[ws.net[,3] > 1, ]

write.table(ws.net, paste0(output.file,'_Final_integrated_network.tsv'), quote = F, row.names = F, col.names = T, sep = '\t')

#get node and edges of final network
igraph_net <- graph.data.frame(ws.net, directed=F)
edge.num <- length(E(igraph_net))
node.num <- length(V(igraph_net))

cat(paste('The final network has:\n', node.num, 'nodes and', edge.num, 'edges\n'))
