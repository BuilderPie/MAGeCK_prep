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
load_package(c('tidyverse'))
# load_package(c('ComplexHeatmap'))
# load_package(c('ggfortify'))
# load_package(c('circlize'))
source('utils.R')
# load_package(c('MAGeCKFlute'))
# ======================================================= #
# option_list = list(
#   make_option(c("-c", "--contrast"), type="character", default=NULL, 
#               help="path to contrast table file. Usually follow the pattern of contrast_table.txt", metavar="character"),
#   make_option(c("-g", "--geneSum"), type="character", default=NULL, 
#               help="path to gene summary file. Usually follow the pattern of gene_summary.txt", metavar="character"))
# 
# opt_parser = OptionParser(option_list=option_list);
# opt = parse_args(opt_parser);
# 
# if (is.null(opt$geneSum)){
#   print_help(opt_parser)
#   stop("At least one argument must be supplied (-g / --geneSum)", call.=FALSE)
# }

# ======================================================= #
# ======================================================= #
#' Sub function, which perform visual quality check for post rra run results
#' @param contrast input from command line; path to contrast table.
#' @param geneSum input from command line; path to gene_summary.txt from MAGeCK rra run.
step3b_qc_post_visual = function(contrast, geneSum){
  # ===== read gdata and match with contrast_table.txt
  print(paste0("Start: post QC visual for ", geneSum))
  gdata = read_geneSum(contrast, geneSum)

  # ===== prepare rank data 
  rankdata = gdata[["gdata"]]$Score
  names(rankdata) = gdata[["gdata"]]$id
  # ===== 
  fileName =  paste0(gsub(".gene_summary.txt", "", basename(geneSum)), "_", tolower(gdata[["run"]]))
  fileTitle = paste0(gsub(".gene_summary.txt", "", basename(geneSum)), "_", tolower(gdata[["run"]]))  
  
  # ===== volcano plot
  VolcanoView(gdata[["gdata"]], x = "Score", y = "FDR", Label = "id", main = fileTitle,
              top = 10, topnames = gdata[["topnames"]],
              filename = file.path(dirname(dirname(geneSum)), 'qc', paste0(fileName, "_volcanoPlot.png")),
              width = 3, height = 1.8)
  # ===== rank plot
  RankView(rankdata = rankdata, main = fileTitle,
           top = 10, bottom = 10, genelist = gdata[["topnames"]],
           filename = file.path(dirname(dirname(geneSum)), 'qc', paste0(fileName, "_rankPlot.png")),
           width = 3, height = 1.8)
  # ===== scatter plot
  set.seed(200)
  gdata[["gdata"]]$RandomIndex = sample(dim(gdata[["gdata"]])[1])
  
  p1 = ScatterView(data = gdata[["gdata"]][gdata[["gdata"]]$Score<=0, ], x = "RandomIndex", y = "Score",
                   label = 'id', y_cut = -CutoffCalling(gdata[["gdata"]]$Score, 2),
                   groups = "bottom", group_col = 'blue',
                   toplabels = gdata[["topnames"]], ylab = "Log2FC",
                   main = fileTitle, display_cut = TRUE)
  ggsave(plot = p1, filename = file.path(dirname(dirname(geneSum)), 'qc', paste0(fileName, "_scatterBottom.png")),
         width = 3, height = 1.8, units = 'in', dpi = 300)
  
  p2 = ScatterView(data = gdata[["gdata"]][gdata[["gdata"]]$Score>=0, ], x = "RandomIndex", y = "Score",
                   label = 'id', y_cut = CutoffCalling(gdata[["gdata"]]$Score, 2),
                   groups = "top", group_col = 'red',
                   toplabels = gdata[["topnames"]], ylab = "Log2FC",
                   main = fileTitle, display_cut = TRUE)
  ggsave(plot = p2, filename = file.path(dirname(dirname(geneSum)), 'qc', paste0(fileName, "_scatterTop.png")),
         width = 3, height = 1.8, units = 'in', dpi = 300)
  # ===== histogram
  gdata[["gdata"]]$sign = rep("0", dim(gdata[["gdata"]])[1])
  gdata[["gdata"]]$sign[gdata[["gdata"]]$Score<0] = "neg"
  gdata[["gdata"]]$sign[gdata[["gdata"]]$Score>0] = "pos"
  gdata[["gdata"]]$sign = as.factor(gdata[["gdata"]]$sign)
  p3 = ggplot(gdata[["gdata"]], aes(x=Score, color = sign)) +
    geom_histogram(fill="white",binwidth=0.2, alpha=0.5, position="identity") +
    labs(title = fileTitle, x = "Log2FC", y = "Count") +
    theme(text = element_text(size=6))
  ggsave(plot = p3, filename = file.path(dirname(dirname(geneSum)), 'qc', paste0(fileName, "_histogram.png")),
         width = 3, height = 1.8, units = 'in', dpi = 300)
  # ===== histogram for rawcount.count_normalized.txt
  rawcount_normalized = dir(path = file.path(dirname(contrast), "count"), pattern = "_normalized.txt", full.names = T)
  if (length(rawcount_normalized)>=1){
    ind_rawcount = grep(gsub(pattern = ".txt", "", gdata[["annorow"]]$Count_File),rawcount_normalized)
    rawcount_normalized = rawcount_normalized[ind_rawcount]
    
    rawcount_normalized = read.table(file=rawcount_normalized, header = T, sep = "\t", check.names = F, stringsAsFactors = F, quote="")
    rawcount_normalized = rawcount_normalized[, -grep(pattern = "Gene|sgRNA", x = colnames(rawcount_normalized), ignore.case = T)]
    p4 = ggplot(pivot_longer(data = rawcount_normalized, cols = colnames(rawcount_normalized), values_to = "Normalized_LFC"), aes(x = Normalized_LFC)) +
      geom_histogram(fill = "white", colour = "black", bins = 60) +
      facet_wrap(name ~ ., scales="free", ncol = 2)
    ggsave(plot = p4, filename = file.path(dirname(dirname(geneSum)), 'qc', paste0(fileName, "_histogram_rawcount_normalized.png")),
           units = 'in', dpi = 300)
  }
  # =====
  print(paste0("Finish: post QC visual for ", geneSum))
}
# ======================================================= #
# ======================================================= #
#' Main function, which passes command line args to step3b_qc_post_visual
#' @param contrast input from command line; path to contrast table.
#' @param geneSum input from command line; path to gene_summary.txt from MAGeCK rra run.
# step3b_qc_post_visual(contrast = opt$contrast, geneSum = opt$geneSum)
