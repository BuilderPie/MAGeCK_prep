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
# load_package(c('MAGeCKFlute'))
# ======================================================= #
option_list = list(
  make_option(c("-d", "--dir"), type="character", default=NULL, 
              help="directory name that contains data and contrast table. Usually follow the pattern PMID_LastAuthor_Journal_Year", metavar="character"),
  make_option(c("--run"), type="character", default="rra", 
              help="which mode to run MAGeCK: could be either rra or mle", metavar="character"))
  
  
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$dir)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (-d / --dir)", call.=FALSE)
}
# ======================================================= #
step2_run_mageck_vispr = function(folder, run, dryrun){
  cur_wd = getwd()
  setwd(folder)
  
  system(paste0("bash run_mageck_count.sh"))
  if (any(grepl("rra", run, ignore.case = T))) system(paste0("bash run_mageck_rra.sh"))
  if (any(grepl("mle", run, ignore.case = T))) system(paste0("bash run_mageck_mle.sh"))

  setwd(cur_wd)
}
# ======================================================= #
step2_run_mageck_vispr(folder = opt$dir, run = opt$run, dryrun = FALSE)

