#!/usr/bin/env Rscript
#this code takes an expression data or seurat as an input, runs SAVER multicore -> sort network -> run LLS.py
#rds file must be a seurat file]
#./Rscript.R --help for detail
#Run Saver -> coverage filter -> cut correlation -> LLS.py -> draw benchmark graph



#this code takes saver data, filter by cell type, and runs parallel PCC correlation matrix calculation
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser(description = 'this code takes seurat data, filter by cell type, and runs saver and PCC correlation matrix calculation')
parser$add_argument("-s", "--sort", default='sort', help="run absort to get clearer top bins [default %(default)s]")
parser$add_argument("-c", "--cutoff", type = 'double', default=99, help="99 percent zero genes are removed [default %(default)s]")
parser$add_argument("-n", "--nCore", type = 'integer', default=12, help="number of threads for saver computation [default %(default)s]")
parser$add_argument("-i", "--input", help="FULL PATH to exprs matrix or seurat rds")
parser$add_argument("-r", "--reuse", default='F', help="type path to calculated saver rds object. This will be used instead [default %(default)s]")
parser$add_argument("-o", "--output", help="filename, current directory will be added to prefix") #no dir prefix
parser$add_argument("-d", "--outdir", help="saved directory")
parser$add_argument("-gc", "--genes", help="path to ccds or wanted genes, saved as a character vector RDS file")
parser$add_argument("-p", "--package", help = "provide path to where scNet package is")
parser$add_argument("-gs", "--goldstandard", default='example/input/gold_standard_symbol_HNv3',help = "gold standard to evaluate LLS")


args <- parser$parse_args()

library(Seurat)
library(tictoc)
library(SAVER)
library(reshape2)



path.to.exprs <- args$input
sort.type <- args$sort
ncores <- as.numeric(args$nCore)
output.file <- args$output
cutoff <- 1 - (as.numeric(args$cutoff) / 100)
out_dir <- args$outdir
path.to.package <- args$package
path.to.genes <- args$genes
path.to.reuse <- args$reuse
path.to.gs <- args$goldstandard

setwd(out_dir)

if (args$reuse != 'F'){
  reuse.sv <- TRUE
}else if (args$reuse == 'F'){
  reuse.sv <- FALSE
}

#test arguments
args = commandArgs(trailingOnly=TRUE)
if (sort.type != 'absort' & sort.type != 'sort'){
  cat('Wrong sorting method')
  quit()
} 

if (reuse.sv != F & reuse.sv != T){
  cat('wrong input of saver object usage')
  quit()
}

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
    data <- as.matrix(read.table(file = expr.file, header =T, sep = seperate))

  }
  else if (datatype == 'tsv'){
    seperate = '\t'
    data <- as.matrix(read.table(file = expr.file, header =T, sep = seperate))
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
  #read.table and make appropriate adjustment for saver input
  return(data)
}


#read data and filter by coding genes

coverage <- ReadCoverage()
data <- ReadData(path.to.exprs)
data.prefiltered <- data[rownames(data) %in% coverage, ]

coverage.filtered.genes <- nrow(data.prefiltered)
cat(paste('coverage genes filtering left:', nrow(data.prefiltered), '\n'))

#additionally filter genes that are >99% zeros (nonzero values of 1 percent or higher)
cat(paste('additional filtering to retain over', cutoff * 100, 'percent non-zero values...\n'))
data.prefiltered <- data.prefiltered[apply(data.prefiltered, 1, function(x) sum(x > 0) / length(x)) > cutoff, ]
data.prefiltered <- data.prefiltered[,! colSums(data.prefiltered)== 0 ]
#print how many genes were filtered
genes.left <- nrow(data.prefiltered)

#if total filtered genes are less then 20 percent of entire coverage genes..re-cut! 
if((genes.left / coverage.filtered.genes) < 0.2){
  cat(paste('only', genes.left, 'genes are left..'))
  cat('too many genes were filtered!! increase cutoff parameter')
  quit()
}

cat(paste('removed additional',(coverage.filtered.genes - genes.left), 'genes...\n'))
cat(paste('total genes left:'), genes.left, 'genes\n')

#run saver
if (reuse.sv == F){
  cat('running saver...\n')
  tic('saver running time...')
  saver <- saver(data.prefiltered, ncores = ncores)
  estimate <- saver$estimate
  toc()
  saveRDS(saver, file = paste0(output.file, '_saver.rds'))
 
} else if (reuse.sv == T){
  #read from the saved objects
  cat("Calling saver object...\n")
  saver <- readRDS(file = path.to.reuse)
}


#save memory this takes a lot of space!
rm(data)

gc()

cat('get corr Matrix...\n')
#get corr matrix adjusting for estimtae uncertainty froms saver calculation
corm <- cor.genes(saver)

#convert NaN to NA
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
corm[is.nan(corm)] <- NA

#reshape into edgelist
corm[!lower.tri(corm)] <- NA
net <- melt(corm, na.rm = T)

#cut saver net (for faster calcluation to top 5 percent based on sorting technique)
if (sort.type == 'absort'){
  net <- net[abs(corm[,3]) > quantile(net(corm[,3]), 0.95), ]
}
if (sort.type == 'sort'){
  net <- net[net[,3] > quantile(net[,3], 0.95), ]
}

#set name for network
output1 <- paste0(output.file, '_SV_PCCnet')
cat(paste('saver network name:', output1,'\n'))


#write network
write.table(net, file = output1, quote = F, row.names = F, sep = '\t', col.names = F)

#sort network based on sorting technique
if (sort.type == 'absort'){
  sort.file.name <- paste0(output1, '_absorted')
  system(paste0(path.to.package,'functions/absort ', output1, ncores))
  
}else if(sort.type == 'sort'){
  sort.file.name <- paste0(output1, '_possorted')
  system(paste('sort -nrk 3,3', output1, '>', sort.file.name))
}


#check if the input for regression.py was properly made
if (!file.exists(sort.file.name)){
  cat("Error: Problem with sorting. Read SAVER.rds from file and try again")
  quit()
}


#run regression analysis for net1 and net2
cat("Running LLS.py ...\n")
system(paste0('python3 ',path.to.package,'functions/LLS.py ', sort.file.name, ' ', path.to.package,path.to.gs))


#remove unsorted network
cat("removing unsorted network\n")
system(paste('rm', output1))


#check if the output for regression was properly made
if(!file.exists(paste0(sort.file.name, '.binlls'))){
  cat("Error: Problem with LLS analysis prefiltered. Read SAVER.rds from current directory and try again\n")
  quit()
}
