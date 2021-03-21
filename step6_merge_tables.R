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
load_package(c('ggplot2'))
load_package(c('ggrepel'))
load_package(c('ComplexHeatmap'))
load_package(c('circlize'))
source('utils.R')
# load_package(c('MAGeCKFlute'))
# ================================================================= #
# ================================================================= #
option_list = list(
  make_option(c("-d", "--dir"), type="character", default=NULL, 
              help="directory name that include 'contrast_table.txt' and 'median_lfc' results. Usually follow the pattern Model_Condition_Category_median_lfc.txt", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=3, 
              help="output folder for merged lfc table and contrast table", metavar="character"))
  
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$dir)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (-d / --dir)", call.=FALSE)
}

if (is.null(opt$output)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (-o / --output)", call.=FALSE)
}

# ================================================================= #
# ================================================================= #
signature = c("HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAPBP", "HLA-E", "HLA-F",
              "STAT1", "IRF1", "CD274", "IFNGR1", "IFNGR2", "CARM1",
              "CASP8", "TRADD", "FADD", "TRAPPC12", "MPV17L2",
              "CD40", "GATA3", "XBP1", "IRF4", "IL4", "IL13", "FOXP3", "LAG3",
              "ADAR", "PTPN2", "DCAF15", "CD3D", "LCP2", "CD5", "CBLB", "UBE2N", "MYC")

plot_group = list(essentiality = c('Essentiality', 'SyntheticLethal'),
                  sorting = c('Sorting'),
                  coculture = c('Coculture-NK', 'Coculture-T', 'Coculture-Tumor'),
                  invivo = c('Invivo'))
# ======================================================= #
# ======================================================= #
load_gdata = function(gene_sum, folder, max_limit){
  gdata = read.table(file.path(folder, gene_sum), header = TRUE, sep = "\t", na.strings = "Empty", stringsAsFactors = FALSE)
  gene_sum_name = gsub('normalized_lfc.txt', '', basename(gene_sum))
  colnames(gdata) = paste0(gene_sum_name, colnames(gdata))
  return(gdata)
}
# ======================================================= #
# ======================================================= #
convert_name = function(x){
  if (length(x) == 7){
    paste(x[2], x[3], x[4], x[5], x[6],sep = '_')
  }
}
# ======================================================= #
# ======================================================= #
plot_heatmap = function(plot_group, output_dir, LFC_heatmap){
  annotation_row = data.frame(Category = rep(0, dim(LFC_heatmap)[2]),
                              stringsAsFactors = FALSE)
  annotation_row$Category = sapply(strsplit(colnames(LFC_heatmap), '_'), FUN = function(x) x[length(x)])
  # col_tmp = colnames(LFC_heatmap)
  colnames(LFC_heatmap) = sapply(strsplit(colnames(LFC_heatmap), '_'), FUN = convert_name)
  rownames(annotation_row) = colnames(LFC_heatmap)
  LFC_heatmap[is.na(LFC_heatmap)] <- as.double("NA")
  
  ind = unlist(lapply(plot_group, function(x) grep(x, annotation_row$Category)))
  if (length(ind > 0)){
    png(filename = file.path(output_dir, 'qc_heatmap', paste0("Heatmap_", paste(unlist(plot_group), collapse = "_"), ".png")), 
        width = 24, height = ifelse(length(ind) < 12, 1.5+0.3*length(ind), 0.5*length(ind)), units = 'in', res = 300)
    draw(Heatmap(t(LFC_heatmap[, ind]), col = colorRamp2(c(-4, 4), c("blue", "red")), column_title = '', name = "LFC",
                 show_column_names = T, show_row_names = T))
    dev.off() 
  }

}
# ======================================================= #
# ======================================================= #
qc_visual = function(gdata, output_dir, signature, plot_group){
  ind_gene = grep(paste0('^', (paste(signature, collapse = '$|^')), '$'), gdata$gene)
  LFC_heatmap = gdata[ind_gene, -1]
  rownames(LFC_heatmap) = gdata[ind_gene, "gene"]
  LFC_heatmap <- LFC_heatmap[rowMeans(is.na(LFC_heatmap)) <= 0.5,]
  
  LFC_means <- rowMeans(LFC_heatmap,na.rm = T)
  LFC_heatmap = LFC_heatmap[order(LFC_means, decreasing = T), ]
  
  filter_colSum = colSums(abs(LFC_heatmap), na.rm = TRUE)>0.2
  LFC_heatmap = LFC_heatmap[, filter_colSum]
  
  lapply(plot_group, plot_heatmap, output_dir=output_dir, LFC_heatmap=LFC_heatmap)
  # for (i in 1:length(plot_group)){
  #   # idx_group = which(annotation_row$Category %in% plot_group[[i]])
  #   idx_group = vector()
  #   for (j in 1:length(plot_group[[i]])){
  #     idx_group = c(idx_group, grep(), which(annotation_row$Category %in% plot_group[[i]][j]))
  #   }
  #   colPal = rev(colorRampPalette(c("#c12603", "white", "#0073B6"), space = "Lab")(199))
    # pheatmap::pheatmap(sig_LFC_scaled, limit = c(-6,6),color=colPal, breaks=breaks, border_color=NA,
    # pheatmap::pheatmap(t(sig_LFC_scaled[, idx_group]), limit = c(-4,4), color=colPal,
    #                    breaks = seq(-4, 4, length.out = 200),border_color=NA,
    #                    cluster_row = FALSE, cluster_col = FALSE,
    #                    annotation_row = annotation_row[idx_group, ],
    #                    filename = file.path(path[['LFC_scaled_qc_path']], paste0("QC_Heatmap_scaled_", names(plot_group)[i], ".png")),
    #                    width = 24, height = ifelse(length(idx_group) < 12, 0.3*12, 0.5*length(idx_group)),
    #                    fontsize_row = 22, fontsize_col = 18, fontsize = 18,
    #                    na_col = "grey")
    # filename = file.path(path[['LFC_scaled_qc_path']], paste0("QC_Heatmap_scaled_", names(plot_group)[i], ".svg"))
    # svglite(filename,
    #         width = 24,
    #         height = ifelse(length(idx_group) < 12, 0.3*12, 0.5*length(idx_group)),
    #         system_fonts = list(sans = "Arial"))
    # par(family = "sans")
    # pheatmap::pheatmap(t(sig_LFC_scaled[, idx_group]), limit = c(-4,4), color=colPal,
    #                    breaks = seq(-4, 4, length.out = 200),border_color=NA,
    #                    cluster_row = FALSE, cluster_col = FALSE,
    #                    annotation_row = annotation_row[idx_group, ],
    #                    # filename = file.path(path[['LFC_scaled_qc_path']], paste0("QC_Heatmap_scaled_", names(plot_group)[i], ".png")),
    #                    # width = 24, height = ifelse(length(idx_group) < 12, 0.3*12, 0.5*length(idx_group)),
    #                    fontsize_row = 22, fontsize_col = 18, fontsize = 18,
    #                    na_col = "grey")
    # dev.off()
  # }
}
# ======================================================= #
# ======================================================= #
step6_merge_tables = function(folder, output_dir, signature, plot_group){
  # ===== read gdata and match with contrast_table.txt
  print(paste0("Start: normalize median lfc for ", folder))
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = T)
  if(!dir.exists(file.path(output_dir, 'qc_heatmap'))) dir.create(file.path(output_dir, 'qc_heatmap'), recursive = T)
  if(!dir.exists(file.path(output_dir, 'qc_rankplot'))) dir.create(file.path(output_dir, 'qc_rankplot'), recursive = T)
  if(!dir.exists(file.path(output_dir, 'qc_tables'))) dir.create(file.path(output_dir, 'qc_tables'), recursive = T)
  
  gene_sum = list.files(folder, pattern = "normalized_lfc.txt", recursive = T)
  
  if (length(gene_sum)>0){
    gene_sum_name = gsub('_normalized_lfc.txt', '', basename(gene_sum))
    print(gene_sum_name)
    print(gene_sum)
    # names(gene_sum) = gene_sum_name
    gdata_merge = Reduce(function(x, y) merge(x, y, by="gene"),  lapply(gene_sum, FUN = load_gdata, folder = folder, max_limit = max_limit))
    # print(head(gdata_merge))
    colnames(gdata_merge)[1] = 'gene'
    write.table(x = gdata_merge, file = file.path(output_dir, paste0('all_normalized_lfc.txt')),
                sep = '\t', quote = FALSE, row.names = F, col.names = T)
    
    qc_visual(gdata_merge, output_dir, signature, plot_group)
  }

  print(paste0("Finish: merge lfc for ", folder))
}
# ======================================================= #
step6_merge_tables(folder = opt$dir, output_dir = opt$output, signature, plot_group)

