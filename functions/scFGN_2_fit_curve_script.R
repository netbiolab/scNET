#!/usr/bin/env Rscript
#this code takes BS and SAVER and MC network and .binLLS file and fits sigmoidal or logistic regression model
#make a new column for the integrated column with LLS score
#Rscript script.R [possorted.binlls] [possorted]
#Junha Cha 2021.04

# Seungbyn Baek mod 2021.05.11
# script.R -bl [BS_possorted.binlls] -b [BS_possorted] -sl [SV_possorted.binlls] -s [SV_possorted] [MC_ossorted.binlls] [MC_possorted] [BS_bin] [SV_bin] [MC_bin] [out_dir] [outputfile]

#modified Junha Cha 2022.06.30

suppressPackageStartupMessages(library(argparse))
## input files
parser <- ArgumentParser()
parser$add_argument("-bl", "--BSlls", help="binlls file from BS result")
parser$add_argument("-b", "--BSpossorted", help="Sorted file from BS result")
#parser$add_argument("-bb", "--BSbin", type="integer", help="Bin count from BS result")

parser$add_argument("-scl", "--SClls", help="binlls file from SC result")
parser$add_argument("-sc", "--SCpossorted", help="Sorted file from SC result")
#parser$add_argument("-scb", "--SCbin", type="integer",help="Bin count from SC result")

parser$add_argument("-sl", "--SVlls", help="binlls file from SV result")
parser$add_argument("-s", "--SVpossorted", help="Sorted file from SV result")
#parser$add_argument("-sb", "--SVbin", type="integer",help="Bin count from SV result")


parser$add_argument("-d", "--outdir", help="saving directory")
parser$add_argument("-o", "--output", help="filename, current directory will be added to prefix")

args <- parser$parse_args()


#set parameters
output.file <- args$output
out_dir <- args$outdir




#Library loading
library(drc)

bs_lls_dir <- args$BSlls 
bs_net_dir <- args$BSpossorted
#bs_bin <- args$BSbin

sv_lls_dir <- args$SVlls
sv_net_dir <- args$SVpossorted
#sv_bin <- args$SVbin

mc_lls_dir <- args$SClls
mc_net_dir <- args$SCpossorted
#mc_bin <- args$SCbin

out_dir <- args$outdir
celltype.name <- args$output
setwd(out_dir)

bs_bin <- as.numeric(as.matrix(read.table('BS_bin_tmp',sep='\t',header=F)))
sv_bin <- as.numeric(as.matrix(read.table('SV_bin_tmp',sep='\t',header=F)))
mc_bin <- as.numeric(as.matrix(read.table('SC_bin_tmp',sep='\t',header=F)))


#BS files
bs_lls <- read.table(bs_lls_dir, sep='\t', header=TRUE)
bs_net <- read.table(bs_net_dir, sep='\t')
#SV files
sv_lls <- read.table(sv_lls_dir, sep='\t', header=TRUE)
sv_net <- read.table(sv_net_dir, sep='\t')

#MC files
mc_lls <- read.table(mc_lls_dir, sep='\t', header=TRUE)
mc_net <- read.table(mc_net_dir, sep='\t')


#take sigma plot and derived sigmoidal parameter-4 equation
#BS_cut files
bs_lls.cut <- bs_lls[1:bs_bin,]
bs_net.cut <- bs_net[1:(1000*bs_bin),]
#SV_cut files
sv_lls.cut <- sv_lls[1:sv_bin,]
sv_net.cut <- sv_net[1:(1000*sv_bin),]
#MC_cut files
mc_lls.cut <- mc_lls[1:mc_bin,]
mc_net.cut <- mc_net[1:(1000*mc_bin),]


#BS_plot (before,after)
pdf(file=paste0(output.file,"_BS_benchmark_compare.pdf"), width=14)
par(mfrow = c(1,2))
plot(bs_lls$MeanBinStatistics, bs_lls$BinLLS,main = paste0('BS Regression 1000bin (before)'), xlab = 'avg PCC', ylab = "LLS", pch=19)
plot(bs_lls.cut$MeanBinStatistics, bs_lls.cut$BinLLS,main = paste0('BS Regression 1000bin (',bs_bin,' bins)'), xlab = 'avg PCC', ylab = "LLS", pch=19)
dev.off()

#SV_plot (before,after)
pdf(file=paste0(output.file,"_SV_benchmark_compare.pdf"), width=14)
par(mfrow = c(1,2))
plot(sv_lls$MeanBinStatistics, sv_lls$BinLLS,main = paste0('SAVER Regression 1000bin (before)'), xlab = 'avg PCC', ylab = "LLS", pch=19)
plot(sv_lls.cut$MeanBinStatistics, sv_lls.cut$BinLLS,main = paste0('SAVER Regression 1000bin (',sv_bin,' bins)'), xlab = 'avg PCC', ylab = "LLS", pch=19)
dev.off()


#SC_plot (before,after)
pdf(file=paste0(output.file,"_MC_benchmark_compare.pdf"), width=14)
par(mfrow = c(1,2))
plot(mc_lls$MeanBinStatistics, mc_lls$BinLLS, main = paste0('MetaCell Regression 1000bin (before)'), xlab = 'avg PCC', ylab = "LLS", pch=19)
plot(mc_lls.cut$MeanBinStatistics, mc_lls.cut$BinLLS, main = paste0('MetaCell Regression 1000bin (',mc_bin,' bins)'), xlab = 'avg PCC', ylab = "LLS", pch=19)
dev.off()



fit_LLS <- function(lls.cut, net.cut, output.file, method){
  ## BS_data
  data <- lls.cut[,c('MeanBinStatistics', 'BinLLS')]
  print(data)
  
  #option 1 log logistic
  model <- drm(BinLLS ~ MeanBinStatistics, fct = L.4(), data = data) #parameter 4 log logistic model
  print(model)
  
  #check
  pdf(file=paste0(output.file,"_",method,"_model.pdf"))
  plot(model, log="", main = "BS 4 param Logistic function",pch=19)
  dev.off()
  
  #take LLS values from constructed model
  lls.fit <- predict(model, data.frame(net.cut[,3]))
  max(lls.fit)
  
  
  #make LLS converted column
  net.cut$lls.fit <- lls.fit
  colnames(net.cut) <- c('gene1', 'gene2', 'PCC', 'LLS_fit')
  
  #remove PCC column
  net.cut$PCC <- NULL
  
  #filter to 1 or higher LLS
  net.cut <- net.cut[net.cut[,3] > 1, ]
  
  return(net.cut)
}


bs_net.cut <- fit_LLS(bs_lls.cut,bs_net.cut, output.file = output.file, method = 'BS')
write.table(bs_net.cut, sep='\t',quote = F, row.names = F, col.names = F, file = paste0(output.file,'_BS_fit'))

sv_net.cut <- fit_LLS(sv_lls.cut, sv_net.cut, output.file = output.file, method = 'SV')
write.table(sv_net.cut, sep='\t',quote = F, row.names = F, col.names = F, file = paste0(output.file,'_SV_fit'))

mc_net.cut <- fit_LLS(mc_lls.cut,mc_net.cut, output.file = output.file, method = 'SC')
write.table(mc_net.cut, sep='\t',quote = F, row.names = F, col.names = F, file = paste0(output.file,'_SC_fit'))


