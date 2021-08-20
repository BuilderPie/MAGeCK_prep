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
# load_package(c('MAGeCKFlute'))
# ======================================================= #
# option_list = list(
#   make_option(c("-d", "--dir"), type="character", default=NULL, 
#               help="directory name that contains data and contrast table. Usually follow the pattern PMID_LastAuthor_Journal_Year", metavar="character"))
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
#' Prepare folders.
#' @param folder path to the study folder which contains one contrast_table.txt file.
prepare_folder = function(folder){
      print(paste0("Working directory is: ", folder))
      dir.create(file.path(folder, 'count'), showWarnings = F)
      dir.create(file.path(folder, 'mle'), showWarnings = F)
      dir.create(file.path(folder, 'rra'), showWarnings = F)
      dir.create(file.path(folder, 'qc'), showWarnings = F)
      dir.create(file.path(folder, 'logs'), showWarnings = F)
      dir.create(file.path(folder, 'results'), showWarnings = F)
      dir.create(file.path(folder, 'rawcount'), showWarnings = F)
      dir.create(file.path(folder, 'lib'), showWarnings = F)
      # Sys.chmod(list.dirs(folder, recursive = T), mode = as.character(mode), use_umask = FALSE)
      # Sys.chmod(list.files(folder, recursive = T, full.names = T), mode = as.character(mode), use_umask = FALSE)
}
# ======================================================= #
# ======================================================= #
#' Prepare bash commend to run count module.
#' @param unique_ind row indices of contrast table for all replicates of one study design.
#' @param contrast contrast table
cmd_count_prep = function(unique_count_file_ind, contrast){
  # file_name = sapply(X = 1:length(unique_ind), FUN = function(i) paste0(contrast[unique_ind[i], "Model"], "_", contrast[unique_ind[i], "Condition"], "_", contrast[unique_ind[i], "Category"], "_r", i))
  # cmd_count = paste0("mageck count -k rawcount/", contrast$Count_File[unique_ind]," -l lib/library.csv --norm-method ", contrast$Norm_Method[unique_ind], " -n count/", file_name)
  cmd_count = paste0("mageck count -k rawcount/", contrast$Count_File[unique_count_file_ind[1]]," -l lib/library.csv --norm-method ", contrast$Norm_Method[1], " -n count/", gsub(".txt", "", contrast$Count_File[unique_count_file_ind[1]]))
  
  return(cmd_count)
}
# ======================================================= #
# ======================================================= #
#' Prepare bash commend to run MAGeCK rra module.
#' @param unique_ind row indices of contrast table for all replicates of one study design.
#' @param contrast contrast table
cmd_rra_prep = function(unique_ind, contrast){
  file_name = sapply(X = 1:length(unique_ind), FUN = function(i) paste0(contrast[unique_ind[i], "Model"], "_", contrast[unique_ind[i], "Cell_Type"], "_", contrast[unique_ind[i], "Condition"], "_", contrast[unique_ind[i], "Category"], "_r", i))
  count_name = gsub(".txt", "", contrast$Count_File[unique_ind])
  # cmd_rra = paste0("mageck test -k count/", file_name, ".count_normalized.txt -c ", contrast$Control[unique_ind], " -t ", contrast$Sample[unique_ind], " --norm-method ", contrast$Norm_Method[unique_ind], " -n rra/", file_name)
  cmd_rra = paste0("mageck test -k count/", count_name, ".count_normalized.txt -c ", contrast$Control[unique_ind], " -t ", contrast$Sample[unique_ind], " --norm-method ", contrast$Norm_Method[unique_ind], " -n rra/", file_name)
  
  return(cmd_rra)
}
# ======================================================= #
# ======================================================= #
#' Prepare bash commend to run MAGeCK mle module.
#' @param unique_ind row indices of contrast table for all replicates of one study design.
#' @param contrast contrast table
cmd_mle_prep = function(unique_ind, contrast){
  file_name = sapply(X = 1:length(unique_ind), FUN = function(i) paste0(contrast[unique_ind[i], "Model"], "_", contrast[unique_ind[i], "Cell_Type"], "_", contrast[unique_ind[i], "Condition"], "_", contrast[unique_ind[i], "Category"], "_r", i))
  count_name = gsub(".txt", "", contrast$Count_File[unique_ind])
  # cmd_mle = sapply(X = 1:length(unique_ind), FUN = function(i) paste0("mageck mle -k count/", file_name[i], ".count_normalized.txt -d designmatrix/designmatrix_", i, ".txt", " --norm-method ", contrast$Norm_Method[unique_ind[i]], " -n mle/", file_name[i]))
  cmd_mle = paste0("mageck mle -k count/", count_name, ".count_normalized.txt -d designmatrix/designmatrix_", ".txt", " --norm-method ", contrast$Norm_Method[unique_ind], " -n mle/", file_name)
  
  return(cmd_mle)
}
# ======================================================= #
# ======================================================= #
#' Prepare bash files to run MAGeCK rra and mle modules.
#' @param folder path to the study folder which contains one contrast_table.txt file.
prepare_command = function(folder){
  contrast = read.table(file.path(folder, "contrast_table.txt"), sep = "\t", header = TRUE, na.strings = "Empty", stringsAsFactors = FALSE, check.names = F,  quote = "", comment.char = "")
  unique_count_file = unique(contrast[, c("Count_File")])
  unique_count_file_ind = lapply(1:length(unique_count_file), function(x) {which(contrast$Count_File == unique_count_file[x])})
  
  # cmd_count = unlist(lapply(unique_ind, FUN = cmd_count_prep, contrast = contrast))
  cmd_count = unlist(lapply(unique_count_file_ind, FUN = cmd_count_prep, contrast = contrast))
  file.create(file.path(folder, "run_mageck_count.sh"))
  writeLines(cmd_count, con = file.path(folder, "run_mageck_count.sh"))
  
  unique_anno = unique(contrast[, c("Model", "Condition", "Category")])
  unique_ind = lapply(1:dim(unique_anno)[1], function(x) {which(contrast$Model == unique_anno$Model[x] & contrast$Condition == unique_anno$Condition[x] & contrast$Category == unique_anno$Category[x])})
  
  cmd_rra = unlist(lapply(unique_ind, FUN = cmd_rra_prep, contrast = contrast))
  file.create(file.path(folder, "run_mageck_rra.sh"))
  writeLines(cmd_rra, con = file.path(folder, "run_mageck_rra.sh"))
  
  cmd_mle = unlist(lapply(unique_ind, FUN = cmd_mle_prep, contrast = contrast))
  file.create(file.path(folder, "run_mageck_mle.sh"))
  writeLines(cmd_mle, con = file.path(folder, "run_mageck_mle.sh"))
}

# ======================================================= #
#' Prepare folders and bash files to run MAGeCK rra and mle modules.
#' @param folder path to the study folder which contains one contrast_table.txt file.
step1_prepare_files = function(folder){
  print(paste0("Start: step 1 prepare folder and command for ", folder))
  prepare_folder(folder)
  prepare_command(folder)
  print(paste0("Finish: step 1 prepare folder and command for ", folder))
}
# ======================================================= #
#' Main function, which passes command line args to step1_prepare_files
#' @param folder input from command line; path to the study folder which contains one contrast_table.txt file.
# step1_prepare_files(folder = opt$dir, mode = 774)
# step1_prepare_files(folder = opt$dir)
