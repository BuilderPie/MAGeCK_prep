#!/usr/bin/env Rscript
# --------------
# Date:  2021-02-28 10:23:10
# Author:Dian Li
# Email: lidian@zju.edu.cn
# --------------
# About project: extract RRA / MLE results to LFC file.
# rm(list = ls())
# cat('\014')
# ======================================================= #
load_package = function(pkgs){
  if(!is.element('BiocManager', installed.packages()[,1])){
    install.packages('BiocManager')
  }
  for(el in pkgs){
    if (!is.element(el, installed.packages()[,1]))BiocManager::install(el)
    suppressWarnings(suppressMessages(invisible(require(el, character.only=TRUE))))
  }
}
# load_package(c('optparse'))
# load_package(c('ggplot2'))
# load_package(c('ggrepel'))
source('utils.R')
# load_package(c('MAGeCKFlute'))
# ================================================================= #
# ================================================================= #
# option_list = list(
#   make_option(c("-d", "--dir"), type="character", default=NULL, 
#               help="directory name that include 'contrast_table.txt' and 'median_lfc' results. Usually follow the pattern Model_Condition_Category_median_lfc.txt", metavar="character"),
#   make_option(c("--max"), type="character", default=3, 
#               help="max value for normalization", metavar="character"))
#   
# opt_parser = OptionParser(option_list=option_list);
# opt = parse_args(opt_parser);
# 
# if (is.null(opt$dir)){
#   print_help(opt_parser)
#   stop("At least one argument must be supplied (-d / --dir)", call.=FALSE)
# }
# ======================================================= #
# ======================================================= #
load_gdata_median = function(gdata_median, max_limit){
  gdata = read.table(gdata_median, header = TRUE, sep = "\t", na.strings = "Empty", stringsAsFactors = FALSE)
  max_limit = as.numeric(max_limit)
  
  # if(grepl("High|Low", condition, ignore.case = TRUE) & grepl('Sorting', category, ignore.case = TRUE)){
  #   gdata$LFC = gdata$LFC * max_limit / max(abs(gdata$LFC), na.rm = T)
  # }else if (dim(gdata)[1] < 5000){
  #   gdata$LFC = gdata$LFC * max_limit / max(abs(gdata$LFC), na.rm = T)
  # }else{
    quant_99 = quantile(gdata$LFC, probs = seq(0,1,0.01))[100]
    # quant_99 = max(abs(gdata$LFC), na.rm = T)
    if (quant_99 < 1e-3){
      quant_99 = quantile(gdata$LFC[gdata$LFC!=0], probs = seq(0,1,0.01))[100]
    }
    gdata$LFC = gdata$LFC * max_limit / quant_99
  # }
  return(gdata)
}

# ======================================================= #
# ======================================================= #
step5_normalize_lfc = function(folder, max_limit){
  # ===== read gdata and match with contrast_table.txt
  print(paste0("Start: normalize median lfc for ", folder))
  
  gdata_median = list.files(folder, pattern = "median_lfc.txt", recursive = T, ignore.case = T, full.names = T)

  if (length(gdata_median)>0){
    gdata_merge = Reduce(function(x, y) merge(x, y, by="gene"),  lapply(gdata_median, FUN = load_gdata_median, max_limit = max_limit))
    colnames(gdata_merge)[-1] = gsub(pattern = "_median_lfc.txt", "", basename(gdata_median))
    write.table(x = gdata_merge, file = file.path(folder, "results", paste0(basename(folder), '_normalized_lfc.txt')),
                sep = '\t', quote = FALSE, row.names = F, col.names = T)

  }

  print(paste0("Finish: normalize median lfc for ", folder))
}
# ======================================================= #
# step5_normalize_lfc(folder = opt$dir, max_limit = opt$max)

