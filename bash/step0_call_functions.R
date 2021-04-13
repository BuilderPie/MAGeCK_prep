#!/usr/bin/env Rscript
# --------------
# Date:  2021-04-26 19:00:10
# Author:Dian Li
# Email: lidian@zju.edu.cn
# --------------
# About project: to organize and call each step function
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
# load_package(c('MAGeCKFlute'))
source("step3b_qc_post_visual.")
# ======================================================= #
option_list = list(
  make_option(c("-d", "--dir"), type="character", default=NULL, 
              help="directory name that contains all the first level study folder", metavar="character"),
  make_option(c("-s", "--study"), type="character", default=NULL, 
              help="directory name that contains specific studies", metavar="character"))

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$dir)){
  print_help(opt_parser)
  stop("At least one of the following argument must be supplied (-d / --dir)", call.=FALSE)
}


# ======================================================= #
#' Detect folders and call each step function
#' @param folder input from command line; path to the folder which contains all studies.
#' @param study input from command line; name of specific studies inside the folder.
step0_call_functions = function(folder, study){
  if (!is.null(study)){
    contrast = list.files(file.path(folder, study), pattern = "contrast_table.txt", recursive = T, full.names = T, ignore.case = T)
  } else{
    contrast = list.files(folder, pattern = "contrast_table.txt", recursive = T, full.names = T, ignore.case = T)
  }
  print(contrast)
  geneSum = lapply(contrast, FUN = function(x) {
    list.files(dirname(x), pattern = "gene_summary.txt", recursive = T, full.names = T, ignore.case = T)
    })
  
  # geneSum = list.files(dirname(contrast), pattern = "gene_summary.txt", recursive = T, full.names = T, ignore.case = T)
  print(geneSum)
}
# ======================================================= #
#' Main function, which passes command line args to step0_call_functions
#' @param folder input from command line; path to the folder which contains all studies.
#' @param study input from command line; name of specific studies inside the folder.
step0_call_functions(folder = opt$dir, study = opt$study)
