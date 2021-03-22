#!/usr/bin/env Rscript
# --------------
# Date:  2021-02-28 10:23:10
# Author:Dian Li
# Email: lidian@zju.edu.cn
# --------------
# About project: extract RRA / MLE results to LFC file.
rm(list = ls())
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
load_package(c('optparse'))
# load_package(c('ggplot2'))
# load_package(c('ggrepel'))
source('utils.R')
# load_package(c('MAGeCKFlute'))
# ================================================================= #
# ================================================================= #
option_list = list(
  make_option(c("-d", "--dir"), type="character", default=NULL, 
              help="directory name that include 'contrast_table.txt' and 'median_lfc' results. Usually follow the pattern Model_Condition_Category_median_lfc.txt", metavar="character"),
  make_option(c("--max"), type="character", default=3, 
              help="max value for normalization", metavar="character"))
  
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$dir)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (-d / --dir)", call.=FALSE)
}
# ======================================================= #
# ======================================================= #
load_gdata = function(gene_sum, folder, max_limit){
  gdata = read.table(file.path(folder, gene_sum), header = TRUE, sep = "\t", na.strings = "Empty", stringsAsFactors = FALSE)
  
  gene_sum_name = gsub('_median_lfc.txt', '', basename(gene_sum))
  condition = strsplit(gene_sum_name, "_")[[1]][2]
  category = strsplit(gene_sum_name, "_")[[1]][3]
  max_limit = as.numeric(max_limit)
  
  if(grepl("High|Low", condition, ignore.case = TRUE) & grepl('Sorting', category, ignore.case = TRUE)){
    gdata$LFC = gdata$LFC * max_limit / max(abs(gdata$LFC), na.rm = T)
  }else if (dim(gdata)[1] < 5000){
    gdata$LFC = gdata$LFC * max_limit / max(abs(gdata$LFC), na.rm = T)
  }else{
    quant_99 = quantile(gdata$LFC, probs = seq(0,1,0.01))[100]
    # quant_99 = max(abs(gdata$LFC), na.rm = T)
    if (quant_99 < 1e-3){
      quant_99 = quantile(gdata$LFC[gdata$LFC!=0], probs = seq(0,1,0.01))[100]
    }
    gdata$LFC = gdata$LFC * max_limit / quant_99
  }
  return(gdata)
}
# ======================================================= #
# ======================================================= #
merge_contrast = function(folder, gene_sum_name){
  contrast_tmp = file.path(folder, paste0(gene_sum_name, "_tmp_contrast.txt"))
  if (any(file.exists(contrast_tmp))){
    contrast_merge = do.call(rbind, lapply(contrast_tmp, FUN = function(x) read.table(x, sep = "\t", header = T)))
    write.table(x = contrast_merge, file = file.path(folder, paste0('contrast_table_merge.txt')),
                sep = '\t', quote = FALSE, row.names = F, col.names = T)
    file.remove(contrast_tmp)
  }
}
# ======================================================= #
# ======================================================= #
step5_normalize_lfc = function(folder, max_limit){
  # ===== read gdata and match with contrast_table.txt
  print(paste0("Start: normalize median lfc for ", folder))
  
  gene_sum = list.files(folder, pattern = "median_lfc.txt", recursive = T)
  
  if (length(gene_sum)>0){
    gene_sum_name = gsub('_median_lfc.txt', '', basename(gene_sum))
    # names(gene_sum) = gene_sum_name
    gdata_merge = Reduce(function(x, y) merge(x, y, by="gene"),  lapply(gene_sum, FUN = load_gdata, folder = folder, max_limit = max_limit))
    colnames(gdata_merge)[-1] = gene_sum_name
    write.table(x = gdata_merge, file = file.path(folder, paste0(basename(folder), '_normalized_lfc.txt')), 
                sep = '\t', quote = FALSE, row.names = F, col.names = T)
    # === merge contrast tables and delete tmp files
    merge_contrast(folder, gene_sum_name)
  }

  print(paste0("Finish: merge lfc for ", folder))
}
# ======================================================= #
step5_normalize_lfc(folder = opt$dir, max_limit = opt$max)

