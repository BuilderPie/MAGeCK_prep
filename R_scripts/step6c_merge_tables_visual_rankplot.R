# #!/usr/bin/env Rscript
# # --------------
# # Date:  2021-02-28 10:23:10
# # Author:Dian Li
# # Email: lidian@zju.edu.cn
# # --------------
# # About project: extract RRA / MLE results to LFC file.
# # rm(list = ls())
# # cat('\014')
# # ======================================================= #
# load_package = function(pkgs){
#   if(!is.element('BiocManager', installed.packages()[,1])){
#     install.packages('BiocManager')
#   }
#   for(el in pkgs){
#     if (!is.element(el, installed.packages()[,1]))BiocManager::install(el)
#     suppressWarnings(suppressMessages(invisible(require(el, character.only=TRUE))))
#   }
# }
# # load_package(c('optparse'))
# load_package(c('ggplot2'))
# load_package(c('ggrepel'))
# load_package(c('ComplexHeatmap'))
# load_package(c('circlize'))
source('utils.R')
# # load_package(c('MAGeCKFlute'))
# # ================================================================= #
# # ================================================================= #
# # option_list = list(
# #   make_option(c("-d", "--dir"), type="character", default=NULL, 
# #               help="directory name that include 'contrast_table.txt' and 'median_lfc' results. Usually follow the pattern Model_Condition_Category_median_lfc.txt", metavar="character"),
# #   make_option(c("-o", "--output"), type="character", default=3, 
# #               help="output folder for merged lfc table and contrast table", metavar="character"))
# # 
# # opt_parser = OptionParser(option_list=option_list);
# # opt = parse_args(opt_parser);
# # 
# # if (is.null(opt$dir)){
# #   print_help(opt_parser)
# #   stop("At least one argument must be supplied (-d / --dir)", call.=FALSE)
# # }
# # 
# # if (is.null(opt$output)){
# #   print_help(opt_parser)
# #   stop("At least one argument must be supplied (-o / --output)", call.=FALSE)
# # }
# 
# # ======================================================= #
# # ======================================================= #
norm_lfc_plot_rankplot = function(plotGroup, output_dir, gdata_merge, posControl){
  each_rankplot = function(ind, output_dir, gdata_merge, posControl, ann_category){
    rankdata = as.numeric(gdata_merge[, ind])
    names(rankdata) = gdata_merge$gene
    rankdata = rankdata[!is.na(rankdata)]
    
    fileName = colnames(gdata_merge)[ind]
    ind_underscore = gregexpr(pattern ='_',fileName)
    fileTitle = sub(paste0("^(.{",ind_underscore[[1]][4]-1,"})."), "\\1\n", fileName)
    RankView(rankdata = rankdata, main = fileTitle,
             top = 0, bottom = 0, genelist = posControl,
             filename = file.path(output_dir, 'qc_rankplot', paste0(ann_category[ind], "_", fileName, ".png")),
             width = 3, height = 1.8)
  }
  
  ann_category = sapply(strsplit(colnames(gdata_merge), '_'), FUN = function(x) x[length(x)-1])

  ind = unlist(lapply(plotGroup, function(x) grep(x, ann_category)))
  if (length(ind > 0)){
    lapply(ind, each_rankplot, output_dir = output_dir, gdata_merge = gdata_merge, posControl = posControl, ann_category = ann_category)
  }

}
# ======================================================= #
# ======================================================= #
step6c_merge_tables_visual_rankplot = function(folder, output_dir){
  print(paste0("Start: normalize lfc QC rankplot for ", folder))
  # if(!dir.exists(file.path(output_dir, 'qc_heatmap'))) dir.create(file.path(output_dir, 'qc_heatmap'), recursive = T)
  if(!dir.exists(file.path(output_dir, 'qc_rankplot'))) dir.create(file.path(output_dir, 'qc_rankplot'), recursive = T)
  # if(!dir.exists(file.path(output_dir, 'qc_tables'))) dir.create(file.path(output_dir, 'qc_tables'), recursive = T)
  
  gdata_merge = read.table(file.path(output_dir, "all_lfc_normalized.txt"), header = TRUE, na.strings = "Empty", stringsAsFactors = FALSE, check.names = F,  quote = "", comment.char = "")
  # === map positive control genes to merged normalized table
  lapply(plotGroup, norm_lfc_plot_rankplot, output_dir=output_dir, gdata_merge=gdata_merge, posControl=posControl)
  
  print(paste0("Finish: normalize lfc QC rankplot for ", folder))
  
}

# # ======================================================= #
# # step6_merge_tables(folder = opt$dir, output_dir = opt$output, signature, plot_group)
# 
