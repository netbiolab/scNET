#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(patchwork)


suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser(description = 'this code takes WSnet, fit LLS and make the final integrated network')
parser$add_argument("-i", "--input", help="WSnet with BS SV SC")
parser$add_argument("-o", "--output", help="filename, will be added as prefix to all output files") 
parser$add_argument("-d", "--outdir", help="saved directory")
parser$add_argument("-bs", "--bigscale2", help="binlls file for bigscale2")
parser$add_argument("-sv", "--saver", help="binlls file for saver")
parser$add_argument("-sc", "--supercell", help="binlls file for supercell")


#set argument
args <- parser$parse_args()


out_dir <- args$outdir
output.file <- args$output
input.path <- args$input
BS_binlls <- args$bigscale2
SV_binlls <- args$saver
SC_binlls <- args$supercell

setwd(out_dir)



tryCatch({
bm_res <- fread(file = BS_binlls)
# Get optimal PCC threshold
prev_lls <- bm_res[1, BinLLS]
pcc_thres <- NA
for (row in seq_len(nrow(bm_res))) {
  cur_lls <- bm_res[row, BinLLS]
  if (row > 1) {
      if (abs(prev_lls - cur_lls) >= 2 || cur_lls <= -0.5) {
        last_pos_row <- row - 1
        while (bm_res[last_pos_row, BinLLS] < 0) {
          last_pos_row <- last_pos_row - 1
        }
        pcc_thres <- bm_res[row - 1, MeanBinStatistics]
        bin <- last_pos_row
        cat(paste0("Detected valid PCC threshold: ", pcc_thres, "\n"))
        cat(paste0("Detected valid Bin: ", bin, "\n"))
        break
    }
  }
  if (row > 2) {
    test_df <- bm_res[c(row - 2, row - 1, row)]
    test_df <- test_df[BinLLS < 0]
    if (nrow(test_df) > 1) {
      last_pos_row <- row - 2
      while(bm_res[last_pos_row, BinLLS] < 0) {
        last_pos_row <- last_pos_row - 1
      }
      pcc_thres <- bm_res[last_pos_row, MeanBinStatistics]
      bin <- last_pos_row
      cat(paste0("Detected valid PCC threshold: ", pcc_thres, "\n"))
      cat(paste0("Detected valid Bin: ", bin, "\n"))
      break
    }
  }
  prev_lls <- cur_lls
  }

BS_bin_suggest <- last_pos_row

p1 <- ggplot(data = bm_res, aes(x = MeanBinStatistics, y = BinLLS)) +
geom_point(show.legend = F) + theme_bw() +
labs(title = paste0('Regression_auto_threshold (BigScale Bin = ',BS_bin_suggest,')'), x = "avg PCC", y = "LLS")
if (!is.na(pcc_thres)) {
pcc_lab <- as.character(round(pcc_thres, 6))
y_pos <- min(bm_res[, BinLLS]) + 1
p1 <- p1 + geom_vline(xintercept = pcc_thres, color = "red", linetype = "dashed")
p1 <- p1 + annotate(geom = "vline", xintercept = pcc_thres, linetype = "dashed") +
    annotate(geom = "text", label = pcc_lab, x = pcc_thres, y = y_pos, color = "red", angle = 90, vjust = 1.2)
}
}, error=function(e) {
    p1 <- ggplot() + theme_void() +labs(title = "error for BigScale")
})

tryCatch({
bm_res <- fread(file = SV_binlls)
# Get optimal PCC threshold
prev_lls <- bm_res[1, BinLLS]
pcc_thres <- NA
for (row in seq_len(nrow(bm_res))) {
  cur_lls <- bm_res[row, BinLLS]
  if (row > 1) {
      if (abs(prev_lls - cur_lls) >= 2 || cur_lls <= -0.5) {
        last_pos_row <- row - 1
        while (bm_res[last_pos_row, BinLLS] < 0) {
          last_pos_row <- last_pos_row - 1
        }
        pcc_thres <- bm_res[row - 1, MeanBinStatistics]
        bin <- last_pos_row
        cat(paste0("Detected valid PCC threshold: ", pcc_thres, "\n"))
        cat(paste0("Detected valid Bin: ", bin, "\n"))
        break
    }
  }
  if (row > 2) {
    test_df <- bm_res[c(row - 2, row - 1, row)]
    test_df <- test_df[BinLLS < 0]
    if (nrow(test_df) > 1) {
      last_pos_row <- row - 2
      while(bm_res[last_pos_row, BinLLS] < 0) {
        last_pos_row <- last_pos_row - 1
      }
      pcc_thres <- bm_res[last_pos_row, MeanBinStatistics]
      bin <- last_pos_row
      cat(paste0("Detected valid PCC threshold: ", pcc_thres, "\n"))
      cat(paste0("Detected valid Bin: ", bin, "\n"))
      break
    }
  }
  prev_lls <- cur_lls
  }


SV_bin_suggest <- last_pos_row

p2 <- ggplot(data = bm_res, aes(x = MeanBinStatistics, y = BinLLS)) +
geom_point(show.legend = F) + theme_bw() +
labs(title = paste0('Regression_auto_threshold (SAVER Bin = ',SV_bin_suggest,')'), x = "avg PCC", y = "LLS")
if (!is.na(pcc_thres)) {
pcc_lab <- as.character(round(pcc_thres, 6))
y_pos <- min(bm_res[, BinLLS]) + 1
p2 <- p2 + geom_vline(xintercept = pcc_thres, color = "red", linetype = "dashed")
p2 <- p2 + annotate(geom = "vline", xintercept = pcc_thres, linetype = "dashed") +
    annotate(geom = "text", label = pcc_lab, x = pcc_thres, y = y_pos, color = "red", angle = 90, vjust = 1.2)
}
}, error=function(e) {
    p2 <- ggplot() + theme_void() +labs(title = "error for SAVER")
})




tryCatch({
  bm_res <- fread(file = SC_binlls)
  # Get optimal PCC threshold
  prev_lls <- bm_res[1, BinLLS]
  pcc_thres <- NA
  for (row in seq_len(nrow(bm_res))) {
    cur_lls <- bm_res[row, BinLLS]
    if (row > 1) {
      if (abs(prev_lls - cur_lls) >= 2 || cur_lls <= -0.5) {
        last_pos_row <- row - 1
        while (bm_res[last_pos_row, BinLLS] < 0) {
          last_pos_row <- last_pos_row - 1
        }
        pcc_thres <- bm_res[row - 1, MeanBinStatistics]
        bin <- last_pos_row
        cat(paste0("Detected valid PCC threshold: ", pcc_thres, "\n"))
        cat(paste0("Detected valid Bin: ", bin, "\n"))
        break
      }
    }
    if (row > 2) {
      test_df <- bm_res[c(row - 2, row - 1, row)]
      test_df <- test_df[BinLLS < 0]
      if (nrow(test_df) > 1) {
        last_pos_row <- row - 2
        while(bm_res[last_pos_row, BinLLS] < 0) {
          last_pos_row <- last_pos_row - 1
        }
        pcc_thres <- bm_res[last_pos_row, MeanBinStatistics]
        bin <- last_pos_row
        cat(paste0("Detected valid PCC threshold: ", pcc_thres, "\n"))
        cat(paste0("Detected valid Bin: ", bin, "\n"))
        break
      }
    }
    prev_lls <- cur_lls
  }
  
  
  SC_bin_suggest <- last_pos_row
  
  p3 <- ggplot(data = bm_res, aes(x = MeanBinStatistics, y = BinLLS)) +
    geom_point(show.legend = F) + theme_bw() +
    labs(title = paste0('Regression_auto_threshold (SuperCell Bin = ',SC_bin_suggest,')'), x = "avg PCC", y = "LLS")
  if (!is.na(pcc_thres)) {
    pcc_lab <- as.character(round(pcc_thres, 6))
    y_pos <- min(bm_res[, BinLLS]) + 1
    p3 <- p3 + geom_vline(xintercept = pcc_thres, color = "red", linetype = "dashed")
    p3 <- p3 + annotate(geom = "vline", xintercept = pcc_thres, linetype = "dashed") +
      annotate(geom = "text", label = pcc_lab, x = pcc_thres, y = y_pos, color = "red", angle = 90, vjust = 1.2)
  }
}, error=function(e) {
  p3 <- ggplot() + theme_void() +labs(title = "error for SuperCell")
})





out_p <- p1 | p2 | p3
ggsave(paste0(output.file,'_Bin_automation.pdf'), plot = out_p, width = 21, height = 7)



#write bin suggest file
write.table(BS_bin_suggest, file="BS_bin_tmp",sep='\t',row.names=F,col.names=F)
write.table(SV_bin_suggest, file="SV_bin_tmp",sep='\t',row.names=F,col.names=F)
write.table(SC_bin_suggest, file="SC_bin_tmp",sep='\t',row.names=F,col.names=F)



  