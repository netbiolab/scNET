#!/usr/bin/env Rscript

#this code draws Weighted SUM network, run regression.py, fit to curve and produce final integrated network
#make a new column for the integrated column with LLS score
#if edge exist in both network, take the largest LLS of the two as integrated link
#check first if the LLS difference is big

# Seungbyn Baek mod : 2021.05.11 
# script.R [WS_net] [out_dir] [filename]

# Junha Cha mod : 2022.06.30
# script.R -i [WS_net] -d [out_dir] -o [filename] 

suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser(description = 'this code takes WSnet, fit LLS and make the final integrated network')
parser$add_argument("-i", "--input", help="WSnet with BS SV SC") #no dir prefix
parser$add_argument("-o", "--output", help="filename, current directory will be added to prefix") #no dir prefix
parser$add_argument("-d", "--outdir", help="saved directory")
parser$add_argument("-p", "--package", help="provide path to where scNet package is")

args <- parser$parse_args()


library(drc)
library(igraph)


# test if there is at least one argument: if not, return an error
# if (length(args)==0) {
#   stop("one argument must be supplied (./WS_net)", call.=FALSE)
# }

out_dir <- args$outdir
output.file <- args$output
input.path <- args$input
path.to.package <- args$package

setwd(out_dir)
ws <- read.table(input.path, sep='\t')

#colnames(ws) <- c('gene1', 'gene2', 'LLS1', 'LLS2','WS')

ws.net <- ws[,c(1,2,ncol(ws))]
colnames(ws.net) <- c('gene1','gene2',"WS")
#if file exists for WS_net
# if (!file.exists('./WS_net')){
#   cat("Error: WS_net file already exsits in this directory!! make sure you make a separate directory for each network")
#   quit()
# }

write.table(ws.net,paste0(output.file,'_WS_net_final'), quote=F, row.names = F, col.names = F, sep='\t')

#run regredssion.py on WS_net
system(paste('python3', paste0(path.to.package,'functions/LLS.py'), paste0(out_dir,'/',output.file, '_WS_net_final')))

lls <- read.table(paste0(output.file,'_WS_net_final.binlls'), header=T)

#remove negative and fit
if( tail(lls$BinLLS,1)<0){
    lls <- lls[-nrow(lls),]
}

model <- drm(BinLLS ~ MeanBinStatistics, fct = L.4(), data = lls) #parameter 4 log logistic model
pdf(paste0(output.file,'_WS_regression_curve_fit.pdf'))
plot(model, log="", main = "4 param Logistic function", pch = 19)
dev.off()

lls.fit <- predict(model, data.frame(ws.net$WS))
cat(paste('max fitted LLS value: ', print(max(lls.fit))))

ws.net$LLS_final <- lls.fit

ws.net$WS <- NULL
ws.net <- ws.net[order(ws.net$LLS_final, decreasing=T), ]

#take only 1 or higher LLS value
ws.net <- ws.net[ws.net[,3] > 0.3, ]

write.table(ws.net, paste0(output.file,'_Final_integrated_scFGN.tsv'), quote = F, row.names = F, col.names = T, sep = '\t')

#get node and edges of final network
igraph_net <- graph.data.frame(ws.net, directed=F)
edge.num <- length(E(igraph_net))
node.num <- length(V(igraph_net))

cat(paste('The final network has:\n', node.num, 'nodes and', edge.num, 'edges\n'))
