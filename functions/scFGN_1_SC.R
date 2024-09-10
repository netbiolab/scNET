#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser(description = 'this code takes gene by cell countmatrix, filter by cell type, and runs supercell and PCC correlation matrix calculation')
parser$add_argument("-s", "--sort", default='sort', help="run absort to get clearer top bins [default %(default)s]")
parser$add_argument("-c", "--cutoff", type="integer", default=0.05, help="cutoff threshold for zero enriched gene")
parser$add_argument("-i", "--input", help="PATH to exprs matrix or seurat rds")
parser$add_argument("-o", "--output", help="filename, current directory will be added to prefix") #no dir prefix
parser$add_argument("-d", "--outdir", help="path to saved directory")
parser$add_argument("-gc", "--genes", help="path to ccds or wanted genes, saved as a character vector RDS file")
parser$add_argument("-g", "--gamma", default=20, help="supercell gamma parameter")
parser$add_argument("-p", "--package", help = "provide path to where scNet package is")
parser$add_argument("-gs", "--goldstandard", default='example/input/gold_standard_symbol_Hnv3.rds',help = "gold standard to evaluate LLS")


args <- parser$parse_args()

data.path <- args$input
cutoff <- args$cutoff
sort.type <- args$sort
data.name <- args$output #preferably celltype name
folder.path <- args$outdir
path.to.genes <- args$genes
gamma <- args$gamma
path.to.package <- args$package
path.to.gs <- args$goldstandard

ReadCoverage <- function() {
  cat("Reading Coverage\n")
  datatype <- tail(unlist(strsplit(path.to.genes, "\\.")), n=1)
  if (datatype == 'csv'){
    seperate = ','
    coverage <- read.table(file = paste0(path.to.package,path.to.genes), header =T, sep = seperate)
  }
  else if (datatype == 'tsv'){
    seperate = '\t'
    coverage <- read.table(file =paste0(path.to.package,path.to.genes), header =T, sep = seperate)
  }
  else if(datatype == 'rds'){
    coverage <- readRDS(file=paste0(path.to.package,path.to.genes))
  }
  else{
    coverage <- read.table(file = paste0(path.to.package,path.to.genes), colClasses = 'character')[,1]

  }
  #read coverage and map data
  #check if the input for file exists was properly made
  if (!file.exists(paste0(path.to.package,path.to.genes))){
    cat("Error: input gene file not found. Check the directory or whether it is an RDS file")
    quit()
  }

  return(coverage)
}

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

  return(data)
}

library(data.table)
library(SuperCell)
library(ggplot2)
library(ggpubr)


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

write.table(SC.GE,file='SuperCell_expression_matrix.tsv',sep='\t',row.names=T,col.names=T,quote=F)
SC.GE <- t(SC.GE)

coverage <- ReadCoverage()

metacells.f <- SC.GE[,colnames(SC.GE) %in% coverage]

corr.mat <- cor(metacells.f, method='pearson')
print(paste0( dim(corr.mat)[1],"*",dim(corr.mat)[1]," correlation matrix is made!"))
corr.mat[!lower.tri(corr.mat)] <- NA

corr.net <- reshape2::melt(corr.mat, na.rm = T)

corr.net <- corr.net[corr.net[,3] > quantile(corr.net[,3], 0.90), ]


output1 <- paste0(data.name, '_SC_PCCnet')
cat(paste('network name:', output1,'\n'))

#write network
write.table(corr.net, file = output1, quote = F, row.names = F, sep = '\t', col.names = F)

sort.file.name <- paste0(output1, '_possorted')
system(paste('sort -nrk 3,3', output1, '>', sort.file.name))


#run regression analysis for net1 and net2
cat("Running LLS.py ...\n")
system(paste('python3 ',path.to.package,'/functions/LLS.py', sort.file.name, ' ', path.to.package,path.to.gs))

#remove unsorted network
cat("removing unsorted network\n")
system(paste('rm', output1))


#read LLS benchmark output and save regression plots
LLS.prefiltered <- read.table(file = paste0(sort.file.name, '.binlls'), sep='\t',header = T)
