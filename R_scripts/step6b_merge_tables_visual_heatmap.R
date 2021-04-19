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
load_package(c('ComplexHeatmap'))
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
norm_lfc_plot_heatmap = function(plot_group, output_dir, LFC_heatmap){
  annotation_row = data.frame(Category = rep(0, dim(LFC_heatmap)[2]),
                              stringsAsFactors = FALSE)
  annotation_row$Category = sapply(strsplit(colnames(LFC_heatmap), '_'), FUN = function(x) x[length(x)-1])
  # col_tmp = colnames(LFC_heatmap)
  colnames(LFC_heatmap) = sapply(strsplit(colnames(LFC_heatmap), '_'), FUN = function(x){
    if (length(x) == 8){
      # paste(x[2], x[3], x[4], x[5], x[6],sep = '_')
      paste(x[2], x[3], x[4], x[5], x[6],sep = '_')
    }
  })
  rownames(annotation_row) = colnames(LFC_heatmap)
  LFC_heatmap[is.na(LFC_heatmap)] <- as.double("NA")
  ind = unlist(lapply(plot_group, function(x) grep(x, annotation_row$Category)))

  if (length(ind > 0)){
    df = t(LFC_heatmap[, ind, drop = FALSE])
    png(filename = file.path(output_dir, 'qc_heatmap', paste0("Heatmap_", paste(unlist(plot_group), collapse = "_"), ".png")),
        width = 12, height = ifelse(length(ind) < 12, 1.8+0.3*length(ind), 0.5*length(ind)), units = 'in', res = 300)
    col_map = colorRamp2(seq(-4, 4, length = 3), c("blue", "#EEEEEE", "red"), space = "RGB")
    draw(Heatmap(df, na_col = "grey", col = col_map, name = "LFC", rect_gp = gpar(col = "white", lwd = 2),
                 cluster_rows = FALSE, cluster_columns = F, row_names_max_width = max_text_width(rownames(df), gp = gpar(fontsize = 12)),
                 show_column_names = T, show_row_names = T, column_title = paste(unlist(plot_group), collapse = "_"),
                 cell_fun = function(j, i, x, y, w, h, fill) {
                   if (!is.na(df[i, j])){
                     if (abs(df[i, j]) >= 2) {
                       grid.text("*", x, y)
                     }
                   }
                 }
                 ))
    dev.off()
  }

}

# ======================================================= #
# ======================================================= #
step6b_merge_tables_visual_heatmap = function(folder, output_dir){
  print(paste0("Start: normalized lfc QC heatmap for ", folder))
  if(!dir.exists(file.path(output_dir, 'qc_heatmap'))) dir.create(file.path(output_dir, 'qc_heatmap'), recursive = T)
  # if(!dir.exists(file.path(output_dir, 'qc_rankplot'))) dir.create(file.path(output_dir, 'qc_rankplot'), recursive = T)
  # if(!dir.exists(file.path(output_dir, 'qc_tables'))) dir.create(file.path(output_dir, 'qc_tables'), recursive = T)
  
  gdata_merge = read.table(file.path(output_dir, "all_lfc_normalized.txt"), header = TRUE, na.strings = "Empty", stringsAsFactors = FALSE, check.names = F,  quote = "", comment.char = "")
  # === map positive control genes to merged normalized table
  ind_gene = grep(paste0('^', (paste(posControl, collapse = '$|^')), '$'), gdata_merge$gene)
  LFC_heatmap = gdata_merge[ind_gene, -1, drop = FALSE]
  rownames(LFC_heatmap) = gdata_merge[ind_gene, "gene"]
  LFC_heatmap <- LFC_heatmap[rowMeans(is.na(LFC_heatmap)) <= 0.9,, drop = FALSE]
  LFC_heatmap[] <- sapply(LFC_heatmap, as.numeric)
  
  LFC_means <- rowMeans(LFC_heatmap,na.rm = T)
  LFC_heatmap = LFC_heatmap[order(LFC_means, decreasing = T), ]

  filter_colSum = colSums(abs(LFC_heatmap), na.rm = TRUE)>0.2
  LFC_heatmap = LFC_heatmap[, filter_colSum, drop = FALSE]
  
  lapply(plotGroup, FUN = norm_lfc_plot_heatmap, output_dir=output_dir, LFC_heatmap=LFC_heatmap)
  
  print(paste0("Finish: normalized lfc QC heatmap for ", folder))
  
}
# # ======================================================= #
# # step6_merge_tables(folder = opt$dir, output_dir = opt$output, signature, plot_group)
# 
