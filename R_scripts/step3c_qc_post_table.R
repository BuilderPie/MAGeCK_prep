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
# load_package(c('ggplot2'))
# load_package(c('ggrepel'))
# load_package(c('ComplexHeatmap'))
# load_package(c('ggfortify'))
# load_package(c('circlize'))
source('utils.R')
# load_package(c('MAGeCKFlute'))
# ======================================================= #
# option_list = list(
#   make_option(c("-c", "--contrast"), type="character", default=NULL, 
#               help="path to contrast table file. Usually follow the pattern of contrast_table.txt", metavar="character"),
#   make_option(c("-g", "--geneSum"), type="character", default=NULL, 
#               help="path to gene summary file. Usually follow the pattern of gene_summary.txt", metavar="character"))
# 
# opt_parser = OptionParser(option_list=option_list);
# opt = parse_args(opt_parser);
# 
# if (is.null(opt$geneSum)){
#   print_help(opt_parser)
#   stop("At least one argument must be supplied (-g / --geneSum)", call.=FALSE)
# }

# ======================================================= #
# ======================================================= #
#' Sub function, which perform auto quality check for post rra run results, and generate evaluation tables
#' @param contrast input from command line; path to contrast table.
#' @param geneSum input from command line; path to gene_summary.txt from MAGeCK rra run.
step3c_qc_post_table = function(contrast, geneSum, output="quantile_chebyshev.csv"){
  # ===== read gdata and match with contrast_table.txt
  print(paste0("Start: post QC table for ", geneSum))
  gdata = read_geneSum(contrast, geneSum)
  # =====
  fileName =  paste0(gsub(".gene_summary.txt", "", basename(geneSum)), "_", tolower(gdata[["run"]]))
  
  # ===== generate random index for shapiro.test
  set.seed(100)
  rand_neg = sample(which(gdata[["gdata"]]$Score < 0))
  if (length(rand_neg)>5000) {
    rand_neg = rand_neg[1:5000]
  } else{
    rand_neg = rand_neg[1:length(rand_neg)]
  }
  
  set.seed(100)
  rand_pos = sample(which(gdata[["gdata"]]$Score > 0))
  if (length(rand_pos)>5000) {
    rand_pos = rand_pos[1:5000]
  } else{
    rand_pos = rand_pos[1:length(rand_pos)]
  }
  # ===== shapiro.test for normal distribution
  if(grepl("High", gdata[["annorow"]]$Condition, ignore.case = TRUE) &
     grepl('Sorting', gdata[["annorow"]]$Category, ignore.case = TRUE)){
    shapiro_W_neg = 0
    shapiro_p_neg = 1
    shapiro_W_pos = shapiro.test(gdata[["gdata"]]$Score[rand_pos])$statistic
    shapiro_p_pos = shapiro.test(gdata[["gdata"]]$Score[rand_pos])$p.value
  } else if (grepl("Low", gdata[["annorow"]]$Condition, ignore.case = TRUE) &
             grepl('Sorting', gdata[["annorow"]]$Category, ignore.case = TRUE)){
    shapiro_W_pos = 0
    shapiro_p_pos = 1
    shapiro_W_neg = shapiro.test(gdata[["gdata"]]$Score[rand_neg])$statistic
    shapiro_p_neg = shapiro.test(gdata[["gdata"]]$Score[rand_neg])$p.value
  } else{
    shapiro_W_pos = shapiro.test(gdata[["gdata"]]$Score[rand_pos])$statistic
    shapiro_p_pos = shapiro.test(gdata[["gdata"]]$Score[rand_pos])$p.value
    shapiro_W_neg = shapiro.test(gdata[["gdata"]]$Score[rand_neg])$statistic
    shapiro_p_neg = shapiro.test(gdata[["gdata"]]$Score[rand_neg])$p.value
  }
  # ===== 
  gdata_neg = gdata[["gdata"]][gdata[["gdata"]]$Score <= 0, ]
  score_sd_neg = sd(gdata_neg$Score)
  gdata_pos = gdata[["gdata"]][gdata[["gdata"]]$Score > 0, ]
  score_sd_pos = sd(gdata_pos$Score)
  # ===== calculate the median rank / position of positive control genes in the gdata. The smaller the value, the higher the rank in both direction
  topnames_order_neg = grep(paste0('^', paste(gdata[["topnames"]],collapse="$|^"), '$'), gdata_neg$id[order(gdata_neg$Score, decreasing = F)], ignore.case = T)
  topnames_order_pos = grep(paste0('^', paste(gdata[["topnames"]],collapse="$|^"), '$'), gdata_pos$id[order(gdata_pos$Score, decreasing = T)], ignore.case = T)
  topnames_order = median(c(topnames_order_neg, topnames_order_pos))
  # ===== assemble the QC evalutation output in the form of data.frame
  tmp2 = data.frame(
    dim(gdata[["gdata"]])[1],
    # ifelse(grepl(annorow$PassQC, "Y", ignore.case = T), 1, 0),
    shapiro_W_neg, shapiro_p_neg,shapiro_W_pos, shapiro_p_pos,
    sum(gdata[["gdata"]]$Score<0), sum(gdata[["gdata"]]$Score<0) / dim(gdata[["gdata"]])[1],
    sum(gdata[["gdata"]]$Score==0), sum(gdata[["gdata"]]$Score==0) / dim(gdata[["gdata"]])[1],
    sum(gdata[["gdata"]]$Score>0), sum(gdata[["gdata"]]$Score>0) / dim(gdata[["gdata"]])[1],
    quantile(gdata[["gdata"]]$Score, probs = seq(0,1,0.05))[20],
    quantile(gdata[["gdata"]]$Score, probs = seq(0,1,0.01))[100],
    
    sum(gdata_neg$Score > (-1 * score_sd_neg)) / dim(gdata[["gdata"]])[1],
    sum(gdata_neg$Score > (-2 * score_sd_neg)) / dim(gdata[["gdata"]])[1],
    # sum(gdata_neg$Score > (-3 * score_sd_neg)) / dim(gdata)[1],
    # sum(gdata_neg$Score > (-4 * score_sd_neg)) / dim(gdata)[1],
    sum(gdata_pos$Score < (1 * score_sd_pos)) / dim(gdata[["gdata"]])[1],
    sum(gdata_pos$Score < (2 * score_sd_pos)) / dim(gdata[["gdata"]])[1],
    # sum(gdata_pos$Score < (3 * score_sd_pos)) / dim(gdata)[1],
    # sum(gdata_pos$Score < (4 * score_sd_pos)) / dim(gdata)[1],
    # sd(score_sorted[1:100]) / mean(score_sorted[1:100]),
    # sd(score_sorted[(length(score_sorted)-99):length(score_sorted)]) / mean(score_sorted[(length(score_sorted)-99):length(score_sorted)]),
    # median(topnames_order[topnames_order <= (dim(gdata)[1] / 2)]) / dim(gdata)[1],
    # median(topnames_order[topnames_order > (dim(gdata)[1] / 2)]) / dim(gdata)[1]
    topnames_order / dim(gdata[["gdata"]])[1]
  )
  # ===== name the QC data frame
  # rownames(quantile_output_chebyshev)[dim(quantile_output_chebyshev)[1]] = fileName
  rownames(tmp2) = fileName
  
  colnames(tmp2) = c('numGenes',
                     'shapiro_W_neg', 'shapiro_p_neg', 'shapiro_W_pos', 'shapiro_p_pos',
                     '#_neg', 'ratio_neg', "#_0", 'ratio_0', "#_pos", 'ratio_pos',
                     "95%", "99%",
                     # '1_SD_neg', '2_SD_neg', '3_SD_neg', '4_SD_neg',
                     # '1_SD_pos', '2_SD_pos', '3_SD_pos', '4_SD_pos',
                     '1_SD_neg', '2_SD_neg',
                     '1_SD_pos', '2_SD_pos',
                     # 'CV_neg100', 'CV_pos100',
                     # 'controlGene_median_rank_bottom', 'controlGene_median_rank_top')
                     'controlGene_median_rank')
  # ===== calculate autoQC result
  if(grepl("High", gdata[["annorow"]]$Condition, ignore.case = TRUE) &
     grepl('Sorting', gdata[["annorow"]]$Category, ignore.case = TRUE)){
    tmp2_quantile_1 = quantile(gdata[["gdata"]]$Score, probs = seq(0,1,0.01))[2]

    autoQC = (tmp2[, 'controlGene_median_rank'] < 0.35) &
      (tmp2_quantile_1 < 0) &
      (tmp2[, 'ratio_pos'] < 0.75) &
      (tmp2[, 'numGenes'] > 100)
  } else if (grepl("Low", gdata[["annorow"]]$Condition, ignore.case = TRUE) &
             grepl('Sorting', gdata[["annorow"]]$Category, ignore.case = TRUE)){
    autoQC = (tmp2[, 'controlGene_median_rank'] < 0.35) &
      (tmp2[, '99%'] > 0) &
      (tmp2[, 'ratio_neg'] < 0.75) &
      (tmp2[, 'numGenes'] > 100)
  } else {
    autoQC = (tmp2[, 'controlGene_median_rank'] < 0.35) &
      (tmp2[, '99%'] > 0) &
      (tmp2[, 'ratio_pos'] < 0.75) &
      (tmp2[, 'ratio_neg'] < 0.75) &
      (tmp2[, 'numGenes'] > 100)
  }
  
  
  autoQC = ifelse(unname(autoQC), 1, 0)
  tmp2 = cbind(autoQC, tmp2)
  # ===== append QC table to the existing one, or create new one if not existed.
  # print(tmp2)
  # if (dim(quantile_output_chebyshev)[1] > 0) {
  if (length(list.files(path = file.path(dirname(dirname(geneSum)), 'qc'), output)) > 0){
    output_file = read.csv(file = file.path(dirname(dirname(geneSum)), 'qc', output), stringsAsFactors = FALSE, na.strings = c("","NA"), check.names = F, row.names = 1, header = T)
    if (!is.vector(output_file)) {
      rowname_tmp2 = rownames(tmp2)
      rowname_tmp2 = paste0("^", gsub("\\+", "\\\\+", rowname_tmp2), "$")
      rowname_check2 = grep(paste(rowname_tmp2, collapse = '|'), rownames(output_file))
      
      if (length(rowname_check2)>0){
        output_file[rowname_check2, ] = tmp2
      } else{
        output_file = rbind(output_file, tmp2)
      }
      
      output_file = output_file[order(row.names(output_file)), ]  # order the output based on row name
    }
  } else{
    output_file = tmp2
  }
  # ===== write output file to local qc folder
  write.csv(output_file, file = file.path(dirname(dirname(geneSum)), 'qc', output), row.names = T)
  # ===== 
  print(paste0("Finish: post QC table for ", geneSum))
}
# ======================================================= #
#' Main function, which perform auto quality check for post rra run results, and generate evaluation tables
#' @param contrast input from command line; path to contrast table.
#' @param geneSum input from command line; path to gene_summary.txt from MAGeCK rra run.
# step3c_qc_post_table(contrast = opt$contrast, geneSum = opt$geneSum, output="quantile_chebyshev.csv")

