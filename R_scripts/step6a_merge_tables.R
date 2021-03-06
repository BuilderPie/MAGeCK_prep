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
load_package(c('ggplot2'))
load_package(c('ggrepel'))
load_package(c('ComplexHeatmap'))
load_package(c('circlize'))
source('utils.R')
# load_package(c('MAGeCKFlute'))
# ================================================================= #
# ================================================================= #
# option_list = list(
#   make_option(c("-d", "--dir"), type="character", default=NULL, 
#               help="directory name that include 'contrast_table.txt' and 'median_lfc' results. Usually follow the pattern Model_Condition_Category_median_lfc.txt", metavar="character"),
#   make_option(c("-o", "--output"), type="character", default=3, 
#               help="output folder for merged lfc table and contrast table", metavar="character"))
# 
# opt_parser = OptionParser(option_list=option_list);
# opt = parse_args(opt_parser);
# 
# if (is.null(opt$dir)){
#   print_help(opt_parser)
#   stop("At least one argument must be supplied (-d / --dir)", call.=FALSE)
# }
# 
# if (is.null(opt$output)){
#   print_help(opt_parser)
#   stop("At least one argument must be supplied (-o / --output)", call.=FALSE)
# }


# ======================================================= #
# ======================================================= #
step6a_merge_tables = function(folder, output_dir = NULL, output_file = NULL){
  # ===== read gdata and match with contrast_table.txt
  print(paste0("Start: merge normalized median lfc for ", folder))
  
  norm_lfc = list.files(folder, pattern = "normalized_lfc.txt", recursive = T, ignore.case = T, full.names = T)
  median_lfc = list.files(folder, pattern = "median_lfc_output.txt", recursive = T, ignore.case = T, full.names = T)
  # ===== normalized lfc files for each study
  gdata_merge= c()
  for(file in norm_lfc){
    gdata = read.table(file, header = TRUE, sep = "\t", na.strings = "Empty", stringsAsFactors = FALSE, check.names = F)
    colnames(gdata)[2:dim(gdata)[2]] = paste0(dir_check_merge(dirname(file))[1], "_",colnames(gdata)[2:dim(gdata)[2]])
    if (length(gdata_merge)==0) {
      gdata_merge = gdata
      # print(head(gdata_merge))
    }
    else{
      gdata_merge = merge(gdata_merge, gdata, by = "gene", all=T)
      # print(head(gdata))
      # print(head(gdata_merge))
    }
  }
  # ===== un-normalized lfc files for each study for test purposes
  gdata_median_merge= c() 
  for(file in median_lfc){
    gdata = read.table(file, header = TRUE, sep = "\t", na.strings = "Empty", stringsAsFactors = FALSE, check.names = F)
    colnames(gdata)[2:dim(gdata)[2]] = paste0(dir_check_merge(dirname(file))[1], "_",colnames(gdata)[2:dim(gdata)[2]])
    if (length(gdata_median_merge)==0) {
      gdata_median_merge = gdata
    }
    else{
      gdata_median_merge = merge(gdata_median_merge, gdata, by = "gene", all=T)
    }
  }
  # ===== combine new data to existing files
  if (!is.null(output_file)){
    if(file.exists(output_file)){
      gdata_merge_pre = read.table(output_file, header = TRUE, na.strings = "Empty", stringsAsFactors = FALSE, check.names = F,  quote = "", comment.char = "")
      ind_col = unlist(lapply(gsub("\\+", "\\\\+", colnames(gdata_merge)), FUN = function(x) grep(x, colnames(gdata_merge_pre))))
      gdata_merge = merge(gdata_merge, gdata_merge_pre[, -ind_col[-1], drop = F], by = "gene", all=T)
    }else{
      if(!dir.exists(dirname(output_file))) dir.create(dirname(output_file), recursive = T)
    }
    if (grepl(".txt$", output_file)) output_file = gsub(".txt", "", output_file)
    if (grepl(".RDS$", output_file)) output_file = gsub(".RDS", "", output_file)
    write.table(x = gdata_merge, file = paste0(output_file, ".txt"), sep = '\t', quote = FALSE, row.names = F, col.names = T)
    saveRDS(gdata_merge, file = paste0(output_file, ".RDS"))
  }
  
  if (!is.null(output_dir)){
    if(!dir.exists(output_dir)) dir.create(output_dir, recursive = T)
    write.table(x = gdata_merge, file = file.path(output_dir, paste0('all_lfc_normalized.txt')), sep = '\t', quote = FALSE, row.names = F, col.names = T)
    saveRDS(gdata_merge, file = file.path(output_dir, paste0('all_lfc_normalized.RDS')))
    
    write.table(x = gdata_median_merge, file = file.path(output_dir, paste0('all_lfc_median.txt')), sep = '\t', quote = FALSE, row.names = F, col.names = T)
    saveRDS(gdata_median_merge, file = file.path(output_dir, paste0('all_lfc_median.RDS')))
  }
  # print(head(gdata_merge))
  # if (length(norm_lfc)>0){
  #   # names(gene_sum) = gene_sum_name
  #   gdata_merge = Reduce(function(x, y) merge(x, y, by="gene", all=T),  lapply(norm_lfc, FUN = function(file){
  #     gdata = read.table(file, header = TRUE, sep = "\t", na.strings = "Empty", stringsAsFactors = FALSE, check.names = F)
  #     colnames(gdata)[2:dim(gdata)[2]] = paste0(dir_check_merge(dirname(file))[1], "_",colnames(gdata)[2:dim(gdata)[2]])
  #     # print(colnames(gdata))
  #     print(file)
  #   }))
  #   # print(head(gdata_merge))
  #   colnames(gdata_merge)[1] = 'gene'
  #   print(head(gdata_merge))
  #   # write.table(x = gdata_merge, file = file.path(output_dir, paste0('all_normalized_lfc.txt')),
  #   #             sep = '\t', quote = FALSE, row.names = F, col.names = T)
  #   
  # }
  
  print(paste0("Finish: merge normalized median lfc for ", folder))
}
# ======================================================= #
# step6_merge_tables(folder = opt$dir, output_dir = opt$output, signature, plot_group)

