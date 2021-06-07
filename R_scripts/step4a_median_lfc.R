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
#' Sub function, which run median lfc step
#' @param folder input from command line; path to the folder which contains studies.
#' @param search_type optional; default search pattern is rra.

step4a_median_lfc = function(folder, search_type="rra"){
  print(paste0("Start: median lfc for ", folder))
  # gdata = read_geneSum(contrast, geneSum)
  
  contrast = file.path(folder,'contrast_table.txt')
  
  geneSum = list.files(path = folder, pattern = "gene_summary.txt", recursive = T)
  geneSum = geneSum[grep(paste0("^", search_type), geneSum, ignore.case = T)] # to filter rra or mle based on search_type
  geneSum = file.path(folder, geneSum)
  
  qc_quantile_chebyshev_raw = qc_quantile_chebyshev = read.csv(file.path(folder, 'qc', 'quantile_chebyshev.csv'), row.names = 1)
  qc_names_raw_df = as.data.frame(do.call(rbind, strsplit(rownames(qc_quantile_chebyshev), "_")))
  qc_names_raw_df = cbind(qc_names_raw_df, qc_quantile_chebyshev_raw$autoQC)
  colnames(qc_names_raw_df) = c("Model", "Cell_Type", "Condition", "Category", "Replicate", "run", "autoQC")
  
  qc_quantile_chebyshev = qc_quantile_chebyshev[qc_quantile_chebyshev$autoQC==1, ]

  if (dim(qc_quantile_chebyshev)[1]>0){
    qc_names_df = as.data.frame(do.call(rbind, strsplit(rownames(qc_quantile_chebyshev), "_")))
    qc_names_unique = apply(unique(qc_names_df[, c(1:4)]),1,paste,collapse="_")
    # qc_names_unique = gsub("\\+", "\\\\+", qc_names_unique)

    lapply(qc_names_unique, FUN = function(x){
      ind = grep(gsub("\\+", "\\\\+", x), geneSum) # search for + sign, as it needs esccape sign
      if (length(ind) > 0){
        gdata_merge = lapply(geneSum[ind], FUN = read_geneSum, contrast = contrast)
        # print(gdata_merge[[1]][["gdata"]][1:20,1:3])
        gdata_merge = do.call(rbind, lapply(gdata_merge, function(x) return(x[["gdata"]])))[, 1:2] # delete FDR column
        colnames(gdata_merge) = c("gene", "LFC")
        
        gdata_merge = gdata_merge[!(is.na(gdata_merge$gene) | gdata_merge$gene == ""), ]
        # ===== median for gdata_merge
        gdata_merge[, "gene"] = as.factor(gdata_merge[, "gene"])
        gdata_median = by(gdata_merge[, "LFC"], gdata_merge[, "gene"], median)
        gdata_median = data.frame(gene=names(gdata_median), LFC=as.numeric(unname(gdata_median)))
        # ===== save median LFC 
        gdata_median = gdata_median[order(gdata_median[,'gene']), ]
        gdata_median = gdata_median[!(duplicated(gdata_median$gene)|is.na(gdata_median$gene) | gdata_median$gene == ""), ]
        
        fileName = paste0(paste(x, collapse = "_"), "_", search_type, "_median_lfc.txt")
        write.table(x = gdata_median, file = file.path(folder, 'results',fileName),
                    sep = '\t', quote = FALSE, row.names = F, col.names = T)
      }
    })
  }
  
  
  # ===== extend contrast table and save as tmp file for next step merge

  dir_info =  t(as.vector(dir_check(folder)))
  contrast = read.table(contrast, sep = "\t", header = T, check.names = F)
  contrast$Replicate = rep("r1", dim(contrast)[1])
  unique_anno = unique(contrast[, c("Model", "Condition", "Category")])
  lapply(1:dim(unique_anno)[1], function(x) {
    ind = which(contrast$Model == unique_anno$Model[x] & contrast$Condition == unique_anno$Condition[x] & contrast$Category == unique_anno$Category[x])
    contrast$Replicate[ind] <<- paste0("r", c(1:length(ind)))
  })
  
  # print(qc_names_raw_df)
  contrast = (merge(contrast, qc_names_raw_df, by = c("Model", "Condition", "Category", "Replicate"), all.x = TRUE))
  contrast = contrast[, c(dim(contrast)[2], 1:(dim(contrast)[2]-1))]
  contrast = cbind(dir_info, contrast,
              t(unlist(strsplit(dir_info[1], split = "_"))), Sys.Date())
  colnames(contrast)[1:2] = c('DirName', 'SubDir')
  colnames(contrast)[(dim(contrast)[2]-4):dim(contrast)[2]] = c('PMID', 'Last_Author', 'Journal', "Year", "Data_QC")
  write.table(x = contrast, file = file.path(folder, "results", "contrast_extend_tmp.txt"),
              sep = '\t', quote = FALSE, row.names = F, col.names = T)
  
  print(paste0("Finish: median lfc for ", folder))
}
# ======================================================= #
#' Main function, which passes command line args to step0_call_functions
#' @param folder input from command line; path to the folder which contains studies.
#' @param search_type optional; default search pattern is rra.
# step4a_median_lfc(folder = opt$dir, search_type = 'rra')

