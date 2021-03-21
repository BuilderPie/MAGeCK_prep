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
source('utils.R')
# load_package(c('MAGeCKFlute'))
# ================================================================= #
# ================================================================= #
option_list = list(
  make_option(c("-d", "--dir"), type="character", default=NULL, 
              help="directory name that include 'contrast_table.txt' and 'rra/mle' results. Usually follow the pattern Model_Condition_Category_r#.gene_summary.txt", metavar="character"))

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$dir)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (-d / --dir)", call.=FALSE)
}

# ================================================================= #
# ================================================================= #

# ======================================================= #
# ======================================================= #
load_gdata = function(gene_sum, folder){
  gdata = read.table(file.path(folder, gene_sum), header = TRUE, sep = "\t", na.strings = "Empty", stringsAsFactors = FALSE)
  gdata$id = gsub(";.*", "", gdata$id)
  gdata = ReadRRA(gdata, score = 'lfc')
  gdata = gdata[, c("id", "Score")]
  return(gdata)
}

# ======================================================= #
# ======================================================= #
step4_median_lfc = function(folder){
  # ===== read gdata and match with contrast_table.txt
  print(paste0("Start: merge lfc for ", folder))
  
  tbl = read.table(file.path(folder,'contrast_table.txt'), sep = "\t", header = T)
  qc_quantile_chebyshev = read.csv(file.path(folder, 'qc', 'quantile_chebyshev.csv'), row.names = 1)
  
  gene_sum = list.files(path = folder, pattern = "gene_summary.txt", recursive = T)
  all_run =  dirname(gene_sum)
  study = paste0(gsub(".gene_summary.txt", "", basename(gene_sum)), "_", all_run)
  gene_sum = (gene_sum[which(qc_quantile_chebyshev[study, "autoQC"] == 1)])
  
  if (length(gene_sum)>0){
    gdata_merge = do.call(rbind, (lapply(gene_sum, FUN = load_gdata, folder = folder)))
    colnames(gdata_merge) = c("gene", "LFC")
    topnames = unique(unlist(strsplit(as.character(tbl$PosControls[1]), split = "\\,|\\;")))
    
    if(any(grepl("mouse", tbl$Organism[1], ignore.case = TRUE))){
      # === in case there is Entrez id (0.2 as cutoff)
      if (sum(!is.na(as.numeric(gdata_merge[,'gene']))) > 0.2*dim(gdata_merge)[1]){
        num_ind = !is.na(as.numeric(gdata_merge[,'gene']))
        map = TransGeneID(as.numeric(gdata_merge[num_ind,'gene']), fromType = "Entrez", toType = "Symbol", fromOrg = "mmu", toOrg = "hsa")
      } else{
        map = TransGeneID(gdata_merge[,'gene'], "symbol", "symbol", fromOrg = "mmu", toOrg = "hsa")
      }
      # topnames = map$hsa
      mmu_gene = names(map)
      hsa_gene = unname(map)
      map = data.frame(mmu_gene = mmu_gene, hsa_gene = hsa_gene)
      map = map[!(is.na(map$hsa_gene) | (duplicated(map[,c("mmu_gene", "hsa_gene")]))),]
      subcount = merge(gdata_merge, map, by.x = "gene", by.y = "mmu_gene", all = TRUE)
      subcount = subcount[, c(3,2)]
      colnames(subcount)[1] = 'gene'
      # gdata_merge = subcount
      subcount = subcount[!(is.na(subcount$gene) | subcount$gene == ""), ]
      gdata_merge = subcount
      
      topnames = unname(TransGeneID(topnames, "symbol", "symbol", fromOrg = "mmu", toOrg = "hsa"))
    }
    
    if (length(topnames) < 3){
      topnames_merge = c(topnames, "JAK1", "JAK2", "STAT1", "IFNGR2", "TAP1")}
    # === median for gdata_merge
    gdata_merge[, "gene"] = as.factor(gdata_merge[, "gene"])
    gdata_median = by(gdata_merge[, "LFC"], gdata_merge[, "gene"], median)
    gdata_median = data.frame(gene=names(gdata_median), LFC=as.numeric(unname(gdata_median)))
    # === convert values for sorting + high/low conditon
    if(grepl("High", tbl$Condition, ignore.case = TRUE) & grepl('Sorting', tbl$Category, ignore.case = TRUE)){
      set.seed(200)
      gdata_median$LFC[gdata_median$LFC<0] = runif(sum(gdata$Score<0),-0.05,0)
    }
    if(grepl("Low", tbl$Condition, ignore.case = TRUE) & grepl('Sorting', tbl$Category, ignore.case = TRUE)){
      set.seed(200)
      gdata_median$LFC[gdata_median$LFC<0] = runif(sum(gdata$Score<0),-0.05,0)
      gdata_median$LFC = -gdata_median$LFC
    }
    # =========================================== #
    # save raw LFC with QC = Y
    gdata_median = gdata_median[order(gdata_median[,'gene']), ]
    gdata_median = gdata_median[!(duplicated(gdata_median$gene)|is.na(gdata_median$gene) | gdata_median$gene == ""), ]
    
    fileName = paste0(basename(folder), "_median_lfc.txt")
    write.table(x = gdata_median, file = file.path(folder, 'results',fileName),
                sep = '\t', quote = FALSE, row.names = F, col.names = T)
  }
  
  # =========================================== #
  # === visualize
  rankdata = gdata_median$LFC
  names(rankdata) = gdata_median$gene
  fileName = paste0(basename(folder), "_median_lfc_rankPlot.png")
  fileTitle = paste0(basename(folder), "_median_lfc")
  RankView(rankdata = rankdata, main = fileTitle, top = 0, bottom = 0, genelist = topnames,
           filename = file.path(folder, 'results',fileName), width = 3, height = 1.8)
  
  print(paste0("Finish: merge lfc for ", folder))
}
# ======================================================= #
step4_median_lfc(folder = opt$dir)

