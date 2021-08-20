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
# load_package(c('ggrepel'))
load_package(c('ComplexHeatmap'))
load_package(c('ggfortify'))
load_package(c('circlize'))
source('utils.R')
# load_package(c('MAGeCKFlute'))
# ======================================================= #
# option_list = list(
#   make_option(c("-c", "--contrast"), type="character", default=NULL, 
#               help="path to contrast table file. Usually follow the pattern of contrast_table.txt", metavar="character"),
#   make_option(c("-r", "--rawcount"), type="character", default=NULL, 
#               help="path to rawcount file. Usually follow the pattern of rawcount.txt", metavar="character"))
#   
# opt_parser = OptionParser(option_list=option_list);
# opt = parse_args(opt_parser);
# 
# if (is.null(opt$rawcount)){
#   print_help(opt_parser)
#   stop("At least one argument must be supplied (-f / --file)", call.=FALSE)
# }
# ======================================================= #
# ======================================================= #
#' Sub function, which read contrast table and rawcount data
#' @param contrast input from command line; path to contrast table.
#' @param rawcount input from command line; path to rawcount/rawcount.txt.
read_rawcount = function(contrast, rawcount){
  contrast = read.table(contrast, sep = "\t", header = TRUE, na.strings = "Empty", stringsAsFactors = FALSE, check.names = F,  quote = "", comment.char = "")
  rawcount = read.table(rawcount, sep = "\t", header = TRUE, na.strings = "Empty", stringsAsFactors = FALSE, check.names = F,  quote = "", comment.char = "")
  # print(head(rawcount))
  df_out = list(rawcount = rawcount, 
                treatment = unlist(strsplit(as.vector(contrast$Sample), split = '\\,|\\;')),
                control = unlist(strsplit(as.vector(contrast$Control), split = '\\,|\\;')))
  return(df_out)
}
# ======================================================= #
# ======================================================= #
#' Sub function, which perform PCA and correlation analysis for rawcount data
#' @param contrast input from command line; path to contrast table.
step3a_qc_pre_pca = function(contrast){
  main_dir = dir_check(dirname(contrast))[1]
  contrast_dir = dirname(contrast)
  print(paste0("Start: pre QC for ", main_dir))
  if (file.exists(contrast) & length(dir(file.path(contrast_dir, "rawcount")))>0 & length(dir(file.path(contrast_dir, "lib")))>0 ){
    contrast = read.table(contrast, sep = "\t", header = T, check.names = F)
    unique_count_file = unique(contrast[, c("Count_File")])
    unique_count_file_ind = lapply(1:length(unique_count_file), function(x) {which(contrast$Count_File == unique_count_file[x])})

    lapply(unique_count_file_ind, FUN = function(count_ind , contrast, contrast_dir){
      treatment = unlist(strsplit(as.vector(contrast$Sample[count_ind]), split = '\\,|\\;'))
      control = unlist(strsplit(as.vector(contrast$Control[count_ind]), split = '\\,|\\;'))
      
      rawcount =  file.path(contrast_dir, paste0("rawcount/", contrast$Count_File[count_ind[1]]))
      rawcount = read.table(rawcount, sep = "\t", header = TRUE, na.strings = "Empty", stringsAsFactors = FALSE, check.names = F,  quote = "", comment.char = "#", fill=TRUE)
      
      df = as.data.frame(t(rawcount[, c(treatment, control)]))
      df = df[ , colSums(is.na(df)) == 0]
      
      pca_res <- prcomp(df)
      
      df$condition = c(rep("Treatment", length(treatment)), rep("Control", length(control)))
      
      fileName = paste(dir_check(contrast_dir)[1], contrast$Count_File[count_ind[1]], sep = "\n")
      fileSave_cor = file.path(contrast_dir, 'qc', paste0(gsub(".txt", "", contrast$Count_File[count_ind[1]]), "_corrPlot.png"))
      fileSave_pca = file.path(contrast_dir, 'qc', paste0(gsub(".txt", "", contrast$Count_File[count_ind[1]]), "_pcaPlot.png"))
      
      p1 = autoplot(pca_res, data = df, colour = 'condition', main = paste0(fileName, "  rawCount_PCA"))
      ggsave(filename = fileSave_pca, plot = p1, width = 5, height = 4, dpi = 300)

      df = as.data.frame((rawcount[, union(treatment, control)]))
      cor_res = cor(df, use = "complete.obs")
      png(filename = fileSave_cor, width = (6 + 6*(dim(df)[2])/20), height = (4 + 4*(dim(df)[2])/10), units = 'in', res = 300)
      draw(Heatmap(cor_res, col = colorRamp2(c(min(cor_res, na.rm = T), 1), c("white", "red")), column_title = paste0(fileName), name = "Pearson\nCorr"))
      dev.off()
      
    }, contrast = contrast, contrast_dir = contrast_dir)
    
    # df_list = read_rawcount(contrast, rawcount)
    # 
    # 
    # df = as.data.frame(t(df_list[["rawcount"]][, c(df_list[["treatment"]], df_list[["control"]])]))
    # df = df[ , colSums(is.na(df)) == 0]
    # pca_res <- prcomp(df)
    # 
    # df$condition = c(rep("Treatment", length(df_list[["treatment"]])), rep("Control", length(df_list[["control"]])))
    # 
    # check_name = strsplit(basename(dirname(contrast)), "_")[[1]]
    # if (length(check_name) == 4){
    #   if (grepl("^\\d", check_name[1]) & grepl("^\\d", check_name[4])){
    #     fileName = basename(dirname(contrast))
    #     fileSave_cor = file.path(dirname(contrast), 'qc', paste0("rawCount_corrPlot.png"))
    #     fileSave_pca = file.path(dirname(contrast), 'qc', paste0("rawCount_pcaPlot.png"))
    #   }
    # } else{
    #   check_name = strsplit(basename(dirname(dirname(contrast))), "_")[[1]]
    #   if (length(check_name) == 4){
    #     if (grepl("^\\d", check_name[1]) & grepl("^\\d", check_name[4])){
    #       fileName = paste0(basename(dirname(dirname(contrast))), '\n', basename(dirname(contrast)))
    #       fileSave_cor = file.path(dirname(contrast), 'qc', paste0("rawCount_corrPlot.png"))
    #       fileSave_pca = file.path(dirname(contrast), 'qc', paste0("rawCount_pcaPlot.png"))
    #     }
    #   } else{
    #     fileName = paste0(basename(dirname(contrast)))
    #     fileSave_cor = file.path(dirname(contrast), 'qc', paste0("rawCount_corrPlot.png"))
    #     fileSave_pca = file.path(dirname(contrast), 'qc', paste0("rawCount_pcaPlot.png"))
    #   }
    # }
    # 
    # p1 = autoplot(pca_res, data = df, colour = 'condition', main = paste0(fileName, "\nrawCount_PCA"))
    # ggsave(filename = fileSave_pca, plot = p1, width = 5, height = 4, dpi = 300)
    # 
    # cor_res = cor(t(df[, -which(colnames(df)=="condition" | colnames(df)=="batch")]))
    # png(filename = fileSave_cor, width = (6 + 6*length(df$condition)/40), height = (4 + 4*length(df$condition)/30), units = 'in', res = 300)
    # draw(Heatmap(cor_res, col = colorRamp2(c(min(cor_res, na.rm = T), 1), c("white", "red")), column_title = paste0(fileName), name = "Pearson\nCorr"))
    # dev.off()
    
    
    
    # fileName =  paste0(gsub(".gene_summary.txt", "", basename(gene_sum)), "_", tolower(run))
    # p1 = autoplot(pca_res, data = df, colour = 'condition', size = 5, main = paste0(fileName, "\nrawCount_PCA"))
    # ggsave(filename = file.path((gsub(paste0("/", tolower(run), "/.*"), "", gene_sum)), 'qc', paste0(fileName, "_rawCount_pcaPlot.png")), plot = p1, width = 5, height = 4, dpi = 300)
    # 
    # cor_res = cor(t(df[, -which(colnames(df)=="condition" | colnames(df)=="batch")]))
    # png(filename = , width = 6, height = 4, units = 'in', res = 300)
    # draw(Heatmap(cor_res, col = colorRamp2(c(min(cor_res, na.rm = T), 1), c("white", "red")), column_title = paste0(fileName), name = "Pearson\nCorr"))
    # dev.off()
  }
  print(paste0("Finish: pre QC for ", main_dir))
  
}

# ======================================================= #
# ======================================================= #
#' Main function, which passes command line args to step3b_qc_post_visual
#' @param contrast input from command line; path to contrast table.
# step3a_qc_pre_pca(contrast = opt$contrast, rawcount = opt$rawcount)

