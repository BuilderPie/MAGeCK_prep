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
              help="directory name that contains data and contrast table. Usually follow the pattern PMID_LastAuthor_Journal_Year", metavar="character"))

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$dir)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (-d / --dir)", call.=FALSE)
}
# ======================================================= #
prepare_folder = function(folder, mode){
      print(paste0("Working directory is: ", folder))
      dir.create(file.path(folder, 'count'), showWarnings = F)
      dir.create(file.path(folder, 'mle'), showWarnings = F)
      dir.create(file.path(folder, 'rra'), showWarnings = F)
      dir.create(file.path(folder, 'qc'), showWarnings = F)
      dir.create(file.path(folder, 'logs'), showWarnings = F)
      dir.create(file.path(folder, 'results'), showWarnings = F)
      Sys.chmod(list.dirs(folder, recursive = T), mode = as.character(mode), use_umask = FALSE)
      Sys.chmod(list.files(folder, recursive = T, full.names = T), mode = as.character(mode), use_umask = FALSE)
}

prepare_command = function(folder, mode){
  tbl = read.table(file.path(folder, "contrast_table.txt"), sep = "\t", header = T)
  print(tbl)
  unique_anno = unique(tbl[, c("Model", "Condition", "Category")])
  if (dim(unique_anno)[1] == 1){
    file_name = sapply(X = 1:dim(tbl)[1], FUN = function(i) paste0(tbl[i, "Model"], "_", tbl[i, "Condition"], "_", tbl[i, "Category"], "_r", i))
    
    cmd_count = paste0("mageck count -k rawcount/rawcount.txt -l lib/library.csv --norm-method ", tbl$Norm_method[1], " -n count/", file_name[1])
    file.create(file.path(folder, "run_mageck_count.sh"))
    writeLines(cmd_count, con = file.path(folder, "run_mageck_count.sh"))
    
    cmd_rra = sapply(X = 1:dim(tbl)[1], FUN = function(i) paste0("mageck test -k count/", file_name[1], ".count_normalized.txt -c ", tbl$Control[i], " -t ", tbl$Sample[i], " --norm-method ", tbl$Norm_method[i], " -n rra/", file_name[i]))
    file.create(file.path(folder, "run_mageck_rra.sh"))
    writeLines(cmd_rra, con = file.path(folder, "run_mageck_rra.sh"))
    
    cmd_mle = sapply(X = 1:dim(tbl)[1], FUN = function(i) paste0("mageck mle -k count/", file_name[1], ".count_normalized.txt -d designmatrix/designmatrix_", i, ".txt", " --norm-method ", tbl$Norm_method[i], " -n mle/", file_name[i]))
    file.create(file.path(folder, "run_mageck_mle.sh"))
    writeLines(cmd_mle, con = file.path(folder, "run_mageck_mle.sh"))
    
    Sys.chmod(list.files(folder, recursive = T, full.names = T), mode = as.character(mode), use_umask = FALSE)
    }
}


step1_prepare_files = function(folder, mode){
  prepare_folder(folder, mode)
  prepare_command(folder, mode)
}
# ======================================================= #
step1_prepare_files(folder = opt$dir, mode = 774)

