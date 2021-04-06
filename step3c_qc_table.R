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
load_package(c('ggfortify'))
load_package(c('circlize'))
source('utils.R')
# load_package(c('MAGeCKFlute'))
# ======================================================= #
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="file name that include 'gene_summary.txt'. Usually follow the pattern Model_Condition_Category_r#.gene_summary.txt", metavar="character"))

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (-f / --file)", call.=FALSE)
}

# ======================================================= #
# ======================================================= #
qc_quantile_chebyshev = function(gdata, annorow, topnames, gene_sum, run, output){
  fileName =  paste0(gsub(".gene_summary.txt", "", basename(gene_sum)), "_", tolower(run))
  # =============== #
  # f = paste(annorow$Dir_Name, annorow$SubDir, annorow$Condition, sep = '_')
  # =============== #
  set.seed(100)
  rand_neg = sample(which(gdata$Score < 0))
  if (length(rand_neg)>5000) {
    rand_neg = rand_neg[1:5000]
  } else{
    rand_neg = rand_neg[1:length(rand_neg)]
  }
  
  set.seed(100)
  rand_pos = sample(which(gdata$Score > 0))
  if (length(rand_pos)>5000) {
    rand_pos = rand_pos[1:5000]
  } else{
    rand_pos = rand_pos[1:length(rand_pos)]
  }
  
  if(grepl("High", annorow$Condition, ignore.case = TRUE) &
     grepl('Sorting', annorow$Category, ignore.case = TRUE)){
    shapiro_W_neg = 0
    shapiro_p_neg = 1
    shapiro_W_pos = shapiro.test(gdata$Score[rand_pos])$statistic
    shapiro_p_pos = shapiro.test(gdata$Score[rand_pos])$p.value
  } else if (grepl("Low", annorow$Condition, ignore.case = TRUE) &
             grepl('Sorting', annorow$Category, ignore.case = TRUE)){
    shapiro_W_pos = 0
    shapiro_p_pos = 1
    shapiro_W_neg = shapiro.test(gdata$Score[rand_neg])$statistic
    shapiro_p_neg = shapiro.test(gdata$Score[rand_neg])$p.value
  } else{
    shapiro_W_pos = shapiro.test(gdata$Score[rand_pos])$statistic
    shapiro_p_pos = shapiro.test(gdata$Score[rand_pos])$p.value
    shapiro_W_neg = shapiro.test(gdata$Score[rand_neg])$statistic
    shapiro_p_neg = shapiro.test(gdata$Score[rand_neg])$p.value
  }
  
  gdata_neg = gdata[gdata$Score <= 0, ]
  score_sd_neg = sd(gdata_neg$Score)
  gdata_pos = gdata[gdata$Score > 0, ]
  score_sd_pos = sd(gdata_pos$Score)
  
  # topnames_order_neg = order(gdata_neg$Score, decreasing = F)[grep(paste0('^', paste(topnames,collapse="$|^"), '$'), gdata_neg$id, ignore.case = T)]
  topnames_order_neg = grep(paste0('^', paste(topnames,collapse="$|^"), '$'), gdata_neg$id[order(gdata_neg$Score, decreasing = F)], ignore.case = T)
  # topnames_order_pos = order(gdata_pos$Score, decreasing = T)[grep(paste0('^', paste(topnames,collapse="$|^"), '$'), gdata_pos$id, ignore.case = T)]
  topnames_order_pos = grep(paste0('^', paste(topnames,collapse="$|^"), '$'), gdata_pos$id[order(gdata_pos$Score, decreasing = T)], ignore.case = T)
  topnames_order = median(c(topnames_order_neg, topnames_order_pos))
  
  # quantile_output_chebyshev = rbind(quantile_output_chebyshev, c(
  tmp2 = data.frame(
    dim(gdata)[1],
    # ifelse(grepl(annorow$PassQC, "Y", ignore.case = T), 1, 0),
    shapiro_W_neg, shapiro_p_neg,shapiro_W_pos, shapiro_p_pos,
    sum(gdata$Score<0), sum(gdata$Score<0) / dim(gdata)[1],
    sum(gdata$Score==0), sum(gdata$Score==0) / dim(gdata)[1],
    sum(gdata$Score>0), sum(gdata$Score>0) / dim(gdata)[1],
    quantile(gdata$Score, probs = seq(0,1,0.05))[20],
    quantile(gdata$Score, probs = seq(0,1,0.01))[100],
    
    sum(gdata_neg$Score > (-1 * score_sd_neg)) / dim(gdata)[1],
    sum(gdata_neg$Score > (-2 * score_sd_neg)) / dim(gdata)[1],
    # sum(gdata_neg$Score > (-3 * score_sd_neg)) / dim(gdata)[1],
    # sum(gdata_neg$Score > (-4 * score_sd_neg)) / dim(gdata)[1],
    sum(gdata_pos$Score < (1 * score_sd_pos)) / dim(gdata)[1],
    sum(gdata_pos$Score < (2 * score_sd_pos)) / dim(gdata)[1],
    # sum(gdata_pos$Score < (3 * score_sd_pos)) / dim(gdata)[1],
    # sum(gdata_pos$Score < (4 * score_sd_pos)) / dim(gdata)[1],
    # sd(score_sorted[1:100]) / mean(score_sorted[1:100]),
    # sd(score_sorted[(length(score_sorted)-99):length(score_sorted)]) / mean(score_sorted[(length(score_sorted)-99):length(score_sorted)]),
    # median(topnames_order[topnames_order <= (dim(gdata)[1] / 2)]) / dim(gdata)[1],
    # median(topnames_order[topnames_order > (dim(gdata)[1] / 2)]) / dim(gdata)[1]
    topnames_order / dim(gdata)[1]
  )
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
  # ======================= #
  # calculate autoQC result
  autoQC = (tmp2[, 'controlGene_median_rank'] < 0.35) &
                        (tmp2[, '99%'] > 0) &
                        (tmp2[, 'ratio_pos'] < 0.65) &
                        (tmp2[, 'ratio_neg'] < 0.65) &
                        (tmp2[, 'numGenes'] > 100)

  autoQC = ifelse(autoQC, 1, 0)
  tmp2 = cbind(autoQC, tmp2)
  # print(tmp2)
  # if (dim(quantile_output_chebyshev)[1] > 0) {
  if (length(list.files(path = file.path((gsub(paste0("/", tolower(run), "/.*"), "", gene_sum)), 'qc'), output)) > 0){
    output_file = read.csv(file = file.path((gsub(paste0("/", tolower(run), "/.*"), "", gene_sum)), 'qc', output), stringsAsFactors = FALSE, na.strings = c("","NA"), check.names = F, row.names = 1)
    # print(class(output_file))
    if (!is.vector(output_file)) {
      rowname_check2 = grep(paste(rownames(tmp2), collapse = '|'), rownames(output_file))
      # print(rowname_check2)
      # print(rowname_check2)
      if (length(rowname_check2)>0){
        output_file[rowname_check2, ] = tmp2
      } else{
        output_file = rbind(output_file, tmp2)
      }
    }
  } else{
    output_file = tmp2
  }
  # ========================================= #
  write.csv(output_file, file = file.path((gsub(paste0("/", tolower(run), "/.*"), "", gene_sum)), 'qc', output), row.names = T)
}

# ======================================================= #
# ======================================================= #
step3_qc = function(gene_sum){
  # ===== read gdata and match with contrast_table.txt
  print(paste0("Start: QC for ", gene_sum))
  run = strsplit(dirname(gene_sum), split = "/")
  run = run[[1]][length(run[[1]])]
  tbl = read.table(file.path((gsub(paste0("/", tolower(run), "/.*"), "", gene_sum)),'contrast_table.txt'), sep = "\t", header = T)
  ann = strsplit(gsub(".gene_summary.txt", "", basename(gene_sum)), '_')
  annorow = tbl[which(tbl$Model == ann[[1]][1] & tbl$Condition == ann[[1]][2] & tbl$Category == ann[[1]][3])[1], ]
  
  gdata = read.table(gene_sum, header = TRUE, sep = "\t", na.strings = "Empty", stringsAsFactors = FALSE)
  gdata$id = gsub(";.*", "", gdata$id)
  gdata = ReadRRA(gdata, score = 'lfc')
  
  # ===== organize data for Sorting and high / low condition
  if(grepl("High", annorow$Condition, ignore.case = TRUE) & grepl('Sorting', annorow$Category, ignore.case = TRUE)){
    set.seed(200)
    gdata$Score[gdata$Score<0] = runif(sum(gdata$Score<0),-0.05,0)
  }
  if(grepl("Low", annorow$Condition, ignore.case = TRUE) & grepl('Sorting', annorow$Category, ignore.case = TRUE)){
    set.seed(200)
    gdata$Score[gdata$Score<0] = runif(sum(gdata$Score<0),-0.05,0)
    gdata$Score = -gdata$Score
  }
  # ===== extract positive control genes and convert to hsa symbols
  if (!is.na(annorow$PosControls)){
    topnames = unlist(strsplit(as.character(annorow$PosControls), split = "\\,|\\;"))
    
    if(any(grepl("mouse", annorow$Organism, ignore.case = TRUE))){
      map = TransGeneID(topnames, "symbol", "symbol", fromOrg = "mmu", toOrg = "hsa")
      topnames = unname(map)
      # ===== in case there is Entrez id (0.2 as cutoff)
      if (sum(!is.na(as.numeric(gdata$id))) > 0.2*dim(gdata)[1]){
        num_ind = !is.na(as.numeric(gdata$id))
        gdata$id[num_ind] = unname(TransGeneID(as.numeric(gdata$id[num_ind]), fromType = "Entrez", toType = "Symbol", fromOrg = "mmu", toOrg = "hsa"))
      } else{
        gdata$id = TransGeneID(gdata$id, "symbol", "symbol", fromOrg = "mmu", toOrg = "hsa")
      }
    }
    if (length(topnames) < 3){
      topnames = c(topnames, "JAK1", "JAK2", "STAT1", "IFNGR2", "TAP1")
    }
  } else{
    topnames = c("JAK1", "JAK2", "STAT1", "IFNGR2", "TAP1")
  }
  # =====
  qc_quantile_chebyshev(gdata, annorow, topnames, gene_sum, run, output="quantile_chebyshev.csv")
  # ===== 
  # ===== 
  # ===== 
  print(paste0("Finish: QC for ", gene_sum))
}
# ======================================================= #
step3_qc(gene_sum = opt$file)

