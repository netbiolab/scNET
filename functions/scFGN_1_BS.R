#!/usr/bin/env Rscript
#this code takes an expression data as an input, runs bigSCale correlation calculaion, and produce a unsorted PCC graph
#rds file must be a seurat file
#./Rscript.R [-i FULL path to exprs UMI data] [-s sort] [-o 20K_monocytes]
#./Rscript.R --help for detail
#RunbigSCale -> coverage filter -> cut correlation -> LLS.py -> draw benchmark graph
#Junha Cha
#07-12-2019

suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser(description = 'this code takes seurat data, filter by cell type, and runs bigSCale2 and PCC correlation matrix calculation')
parser$add_argument("-s", "--sort", default='sort', help="run absort to get clearer top bins [default %(default)s]")
parser$add_argument("-n", "--nCore", type="integer", default=4, help="number of threads for sorting [default %(default)s]")
parser$add_argument("-c", "--cutoff", type="integer", default=0.95, help="cutoff threshold for compute.network [default %(default)s]")
parser$add_argument("-i", "--input", help="PATH to exprs matrix or seurat rds")
parser$add_argument("-r", "--reuse", default='F', help="type path of calcuated bs object. This will be used instead [default %(default)s]")
parser$add_argument("-o", "--output", help="filename, current directory will be added to prefix") #no dir prefix
parser$add_argument("-d", "--outdir", help="path to saved directory")
parser$add_argument("-gc", "--genes", help="path to ccds or wanted genes, saved as a character vector RDS file")
parser$add_argument("-p", "--package", help = "provide path to where scNet package is")

args <- parser$parse_args()

library(bigSCale)
library(reshape2)
library(float)
library(Seurat)



path.to.exprs <- args$input
cutoff <- args$cutoff
sort.type <- args$sort
ncores <- args$nCore
output.file <- args$output #preferably celltype name
out_dir <- args$outdir
path.to.package <- args$package
path.to.genes <- args$genes
path.to.reuse <- args$reuse

setwd(path.to.package)

if (args$reuse != 'F'){
  reuse.bs <- TRUE
}else if (args$reuse == 'F'){
  reuse.bs <- FALSE
}


#test arguments
args = commandArgs(trailingOnly=TRUE)
if (sort.type != 'absort' & sort.type != 'sort'){
  cat('Wrong sorting method')
  quit()
} 

#if (reuse.bs != F & reuse.bs != T){
#  cat('wrong input of bigscale object usage')
#  quit()
#}


#pubdata.path <- '/home2/bsb0613/Research/RawData/Network_DB'


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
  #read.table and make appropriate adjustment for bigscale input
  return(data)
}


#get coverage genes, and filter the input scRNA seq data
coverage <- ReadCoverage()
data <- ReadData(path.to.exprs)
data.prefiltered <- data[rownames(data) %in% coverage, ] #filter the data to the input genes


#run Bigscale and saveRDS bigscale object
if (reuse.bs == F){
  cat("Running BigScale...\n")
  
  cat(paste('start...\n'))
  bigscale.prefiltered <- compute.network(expr.data = data.prefiltered, gene.names = rownames(data.prefiltered), clustering = 'recursive')
  saveRDS(bigscale.prefiltered, file = paste0(output.file, '_bs.rds'))
} else if (reuse.bs == T){
  #read from the saved objects
  cat("Calling bigscale objects...\n")
  bigscale.prefiltered <- readRDS(file = path.to.reuse)
}



#process corr
zscore <-dbl(bigscale.prefiltered$tot.scores)
colnames(zscore) <- colnames(bigscale.prefiltered$correlations)
#write genespace to know what gene we are making network with (they are filtered)
write.table(colnames(zscore), file = paste0(outfir,"/",output.file,'_Genespace_BS.txt'), quote=F, row.names=F, col.names=F)

cat('calculating PCC...\n')
corr.prefiltered <- cor(zscore, method = 'pearson')


#make PCC unsorted prefiltered using for each 
cat('Melting correlation matrix\n')
corr.prefiltered[!lower.tri(corr.prefiltered)] <- NA
net <- melt(corr.prefiltered, na.rm = T)

if (sort.type == 'absort'){
  net <- net[abs(net$value) > quantile(abs(net$value), cutoff), ] 
}
if (sort.type == 'sort'){
  net <- net[net$value > quantile(net$value, cutoff), ]
  #net <- net[net$value > 0,]
}


output1 <- paste0(out_dir,'/', output.file, '_BS_PCCnet')
cat(paste('BS network name:', output1,'\n'))

#writeoutput in inprogress folder
write.table(net, file = output1, quote = F, row.names = F, sep = '\t', col.names = F)


#system absort commaned for net
cat(paste('sorting with...', sort.type, '\n'))
cat("Sorting network...\n")

if (sort.type == 'absort'){
  sort.file.name <- paste0(output1, '_absorted')
  system(paste(path.to.package,'/functions/absort', output1, ncores))
  
}else if(sort.type == 'sort'){
  sort.file.name <- paste0(output1, '_possorted')
  system(paste('sort -nrk 3,3', output1, '>', sort.file.name))
}



#check if the input for regression.py was properly made
if (!file.exists(sort.file.name)){
  cat("Error: Problem with sorting. Read bigscale.rds from file and try again")
  quit()
}


#run regression analysis for net1 and net2
cat("Running LLS.py ...\n")
system(paste('python3', paste0(path.to.package,'/functions/LLS.py'), sort.file.name))
# system(paste('python3.6 ~/bin/condLLS.py', sort.file.name, paste0('Genespace_BS_', output.file,'.txt')))

#remove unsorted network
cat("removing unsorted network\n")
system(paste('rm', output1))



#check if the output for regression was properly made
if(!file.exists(paste0(sort.file.name, '.binlls'))){
  cat("Error: Problem with LLS.py analysis prefiltered. Read bigscale.rds from inprogress file and try again\n")
  quit()
}

#read LLS benchmark output and save regression plots
LLS.prefiltered <- read.table(file = paste0(sort.file.name, '.binlls'), sep='\t',header = T)


#save plots
cat('draw to pdf...done!\n')
pdf(paste0(output.file, '_BS_benchmark.pdf'), width = 14)
par(mfrow = c(1,2))
plot(LLS.prefiltered$GeneCoverage / 18802 *100, LLS.prefiltered$cumLLS, main = 'Benchmark', xlab = 'coverage', ylab = 'cumLLS', type = 'l')
plot(LLS.prefiltered$MeanBinStatistics,LLS.prefiltered$BinLLS, main = paste0('BS Regression 1000bin (',sort.type,')'), xlab = 'avg PCC', ylab = "LLS", pch=19)
dev.off()