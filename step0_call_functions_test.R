#!/usr/bin/env Rscript
# --------------
# Date:  2021-04-26 19:00:10
# Author:Dian Li
# Email: lidian@zju.edu.cn
# --------------
# About project: to organize and call each step function
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
load_package(c('optparse'))
# load_package(c('MAGeCKFlute'))
# source("step3b_qc_post_visual.")
# ======================================================= #
option_list = list(
  make_option(c("-d", "--dir"), type="character", default=NULL, 
              help="directory name that contains all the first level study folder", metavar="character"),
  make_option(c("-p", "--rpath"), type="character", default=NULL, 
              help="directory name that contains all the R function steps", metavar="character"),
  make_option(c("-o", "--output_dir"), type="character", default=NULL, 
              help="directory name that save the merged lfc table and final qc plots", metavar="character"),
  make_option(c("-f", "--output_file"), type="character", default=NULL, 
              help="path to the existing or to be created merged lfc table; optional; if defined then output files will be saved to its folder", metavar="character"))

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$dir)){
  print_help(opt_parser)
  stop("The following argument must be supplied (-d / --dir)", call.=FALSE)
}

if (is.null(opt$rpath)){
  print_help(opt_parser)
  stop("The following argument must be supplied (-p / --rPath)", call.=FALSE)
}

# ======================================================= #
#' Detect folders and call each step function
#' @param folder input from command line; path to the folder which contains all studies.
#' @param study input from command line; name of specific studies inside the folder.
step0_call_functions = function(folder, rPath, output_dir=NULL, output_file=NULL){
  cur_wd = getwd()
  setwd(rPath)
  lapply(list.files(pattern = "[.]R$", recursive = TRUE), source)
  setwd(cur_wd)
  
  contrast = list.files(folder, pattern = "contrast_table.txt", recursive = T, full.names = T, ignore.case = T)
  
  if (length(contrast) ==0) stop("no active studies has been found")
  # for (i in contrast){
    # step1_prepare_files(dirname(i))
  # }
  
  
  # ============================== #
  # ======= test for parallel
  # https://stackoverflow.com/questions/38318139/run-a-for-loop-in-parallel-in-r
  # library(foreach)
  # library(doParallel)
  
  #setup parallel backend to use many processors
  # cores=detectCores()
  # cl <- makeCluster(cores[1]-1) #not to overload your computer
  # registerDoParallel(cl)
  # 
  # foreach(i=1:length(contrast)) %dopar% {
  #   step1_prepare_files(dirname(contrast[i]))
  #   # tempMatrix = functionThatDoesSomething() #calling a function
  #   #do other things if you want
  #   
  #   # tempMatrix #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
  # }
  # ====================================== #
  # for (i in contrast){
    # step2_run_mageck_vispr(dirname(i), run = "rra")
    # step3a_qc_pre_pca(contrast = i)
  # }
  
  # foreach(i=1:length(contrast)) %dopar% {
    # step2_run_mageck_vispr(dirname(contrast[i]), run = "rra")
    # step3a_qc_pre_pca(contrast = contrast[i])
  # }
  
  
  # geneSumList = lapply(contrast, FUN = function(x) {
  #   list.files(dirname(x), pattern = "gene_summary.txt", recursive = T, full.names = T, ignore.case = T)
  # })

  # for (i in 1:length(geneSumList)){
  #   if (length(geneSumList[[i]]) > 0){
  #     geneSum = geneSumList[[i]]
      # lapply(geneSum, FUN = step3b_qc_post_visual, contrast = contrast[i])
      # lapply(geneSum, FUN = step3c_qc_post_table, contrast = contrast[i])
  #   }
  # }
  
  # for (i in 1:length(contrast)){
  #   step4a_median_lfc(folder = dirname(contrast[i]))
  #   step4b_median_lfc_visual(folder = dirname(contrast[i]))
  #   step5_normalize_lfc(folder = dirname(contrast[i]), max_limit = 2)
  # }
  

  
  # step6a_merge_tables(folder = folder, output_dir = output_dir, output_file = output_file)

  if (!is.null(output_file)){
    step6b_merge_tables_visual_heatmap(folder = folder, output_dir = dirname(output_file))
    step6c_merge_tables_visual_rankplot(folder = folder, output_dir = dirname(output_file))
  }
  if (!is.null(output_dir)){
    step6b_merge_tables_visual_heatmap(folder = folder, output_dir = output_dir)
    step6c_merge_tables_visual_rankplot(folder = folder, output_dir = output_dir)
  }

}
# ======================================================= #
#' Main function, which passes command line args to step0_call_functions
#' @param folder input from command line; path to the folder which contains studies.
step0_call_functions(folder = opt$dir, rPath = opt$rpath, output_dir = opt$output_dir, output_file = opt$output_file)
