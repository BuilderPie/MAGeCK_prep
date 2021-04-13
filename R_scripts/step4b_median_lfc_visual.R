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
source('utils.R')
# load_package(c('MAGeCKFlute'))
# ================================================================= #
# ================================================================= #
# option_list = list(
#   make_option(c("-d", "--dir"), type="character", default=NULL, 
#               help="directory name that include 'contrast_table.txt' and 'rra/mle' results. Usually follow the pattern Model_Condition_Category_r#.gene_summary.txt", metavar="character"))
# 
# opt_parser = OptionParser(option_list=option_list);
# opt = parse_args(opt_parser);
# 
# if (is.null(opt$dir)){
#   print_help(opt_parser)
#   stop("At least one argument must be supplied (-d / --dir)", call.=FALSE)
# }

# ================================================================= #
# ================================================================= #

# ================================================================= #
# ================================================================= #
#' Sub function, which run median lfc step
#' @param folder input from command line; path to the folder which contains studies.
#' @param search_type optional; default search pattern is rra.

step4b_median_lfc_visual = function(folder, search_type="rra"){
  print(paste0("Start: visualize median lfc for ", folder))
  
  gdata_median = list.files(path = folder, pattern = "_median_lfc.txt", recursive = T, full.names = T)
  gdata_median = gdata_median[grep(paste0("_", search_type, "_"), gdata_median, ignore.case = T)] # to filter rra or mle based on search_type
  if (length(gdata_median) > 0){
    gdata_median_df = as.data.frame(do.call(rbind, strsplit(basename(gdata_median), "_")))
    
    contrast = file.path(folder,'contrast_table.txt')
    contrast = read.table(contrast, sep = "\t", header = T, check.names = F)
    
    lapply(1:length(gdata_median), FUN = function(x){
      ind = which(gdata_median_df[x,1] == contrast$Model & gdata_median_df[x,2] == contrast$Condition & gdata_median_df[x,3] == contrast$Category)
      
      df = read.table(gdata_median[x], header = TRUE, sep = "\t", na.strings = "Empty", stringsAsFactors = FALSE)
      topnames = unlist(strsplit(as.character(contrast[ind, "PosControls"]), split = "\\,|\\;"))
      if(any(grepl("mouse", contrast[ind, "Organism"], ignore.case = TRUE))){
        topnames = unname(TransGeneID(topnames, "symbol", "symbol", fromOrg = "mmu", toOrg = "hsa"))
      }
      if (length(topnames) < 3){
        topnames_merge = c(topnames, "JAK1", "JAK2", "STAT1", "IFNGR2", "TAP1")
      }
      # ===== visualize
      rankdata = df$LFC
      names(rankdata) = df$gene
      fileName = paste0(apply(gdata_median_df[x,1:4] , 1 , paste , collapse = "_" ), "_median_lfc_rankPlot.png")
      fileTitle = paste0(apply(gdata_median_df[x,1:4] , 1 , paste , collapse = "_" ), "_median_lfc_rank")
      RankView(rankdata = rankdata, main = fileTitle, top = 0, bottom = 0, genelist = topnames,
               filename = file.path(folder, 'qc',fileName), width = 3, height = 1.8)
    })
  }
  
  
  print(paste0("Finish: visualize median lfc for ", folder))
}

# ======================================================= #
# step4_median_lfc(folder = opt$dir, search_type = 'rra')

