#!/usr/bin/env Rscript

#ex) ./scFGN_1_SC.R -d /home/sej4926/Bottom_up_network/Benchmarking/SuperCell/OVC/Macrophage/ -o OVC_Macrophage_10 -i /home/sej4926/scNET_OVC/tumor_only/Macrophage/OVC_Macrophage_countmatrix.tsv 
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser(description = 'this code takes gene by cell countmatrix, filter by cell type, and runs supercell and PCC correlation matrix calculation')
parser$add_argument("-s", "--sort", default='sort', help="run absort to get clearer top bins [default %(default)s]")
parser$add_argument("-c", "--cutoff", type="integer", default=0.05, help="cutoff threshold for zero enriched gene")
parser$add_argument("-i", "--input", help="PATH to exprs matrix or seurat rds")
#parser$add_argument("-r", "--reuse", default='F', help="type path of calcuated bs object. This will be used instead [default %(default)s]")
parser$add_argument("-o", "--output", help="filename, current directory will be added to prefix") #no dir prefix
parser$add_argument("-d", "--outdir", help="path to saved directory")
parser$add_argument("-g", "--gamma", default=20, help="supercell gamma parameter")


args <- parser$parse_args()

data.path <- args$input
cutoff <- args$cutoff
sort.type <- args$sort
data.name <- args$output #preferably celltype name
folder.path <- args$outdir
gamma <- args$gamma
#path.to.reuse <- args$reuse

ReadData <- function(expr.file) {
  cat("Reading Data\n")
  datatype <- tail(unlist(strsplit(expr.file, "\\.")), n=1)
  if (datatype == 'csv'){
    seperate = ','
    data <- read.table(file = expr.file, header =T, sep = seperate)
  }
  else if (datatype == 'tsv'){
    seperate = '\t'
    data <- read.table(file = expr.file, header =T, sep = seperate)
  }
  else if(datatype == 'rds'){
    seurat <- readRDS(expr.file)
    data <- as.matrix(seurat@assays$RNA@counts)
    rm(seurat)
  }
  else{
    print('unknown file type: check input file')
    quit()
  }
  #read.table and make appropriate adjustment for bigscale input
  return(data)
}

library(data.table)
library(SuperCell)
library(ggplot2)
library(ggpubr)
#library(anndata)

args = commandArgs(trailingOnly=TRUE)

GE <- ReadData(data.path)

SC <- SCimplify(GE,  # gene expression matrix 
                k.knn = 10, # number of nearest neighbors to build kNN network
                gamma = gamma, # graining level
                n.var.genes = 1000 # number of the top variable genes to use for dimentionality reduction 
)

SC.GE <- supercell_GE(GE, SC$membership)
SC.GE <- as.data.frame(as.matrix(SC.GE))
SC.GE <- SC.GE[apply(SC.GE, 1, function(x) sum(x > 0) / length(x)) > cutoff, ]

setwd(folder.path)

write.table(SC.GE,file='zero_SuperCell_expression_matrix.tsv',sep='\t',row.names=T,col.names=T,quote=F)
SC.GE <- t(SC.GE)

hnv3.genes <- readRDS('/home3/junhacha/HumanNetv3/Benchmark/HNv3_coverage_genes.rds')
metacells.f <- SC.GE[,colnames(SC.GE) %in% hnv3.genes]

corr.mat <- cor(metacells.f, method='pearson')
print(paste0( dim(corr.mat)[1],"*",dim(corr.mat)[1]," correlation matrix is made!"))
corr.mat[!lower.tri(corr.mat)] <- NA

corr.net <- reshape2::melt(corr.mat, na.rm = T)

#take only positive values
corr.net <- corr.net[corr.net[,3] > quantile(corr.net[,3], 0.90), ]


output1 <- paste0(data.name, '_SC_PCCnet')
cat(paste('network name:', output1,'\n'))

#write network
write.table(corr.net, file = output1, quote = F, row.names = F, sep = '\t', col.names = F)

sort.file.name <- paste0(output1, '_possorted')
system(paste('sort -nrk 3,3', output1, '>', sort.file.name))


#run regression analysis for net1 and net2
cat("Running LLS.py ...\n")
system(paste('python3 /home3/junhacha/bin/LLS.py', sort.file.name))


#remove unsorted network
cat("removing unsorted network\n")
system(paste('rm', output1))


#read LLS benchmark output and save regression plots
LLS.prefiltered <- read.table(file = paste0(sort.file.name, '.binlls'), sep='\t',header = T)


#save plots
cat('draw to pdf...done!\n')
pdf(paste0(data.name, '_SC_benchmark.pdf'), width = 14)
par(mfrow = c(1,2))
plot(LLS.prefiltered$GeneCoverage / 18802 *100, LLS.prefiltered$cumLLS, main = 'Benchmark', xlab = 'coverage', ylab = 'cumLLS', type = 'l')
plot(LLS.prefiltered$MeanBinStatistics,LLS.prefiltered$BinLLS, main = 'Regression 1000bin', xlab = 'avg PCC', ylab = "binLLS", pch=19)
dev.off()