# --------------
# Date:  2020-09-29 22:00:00
# Author:Dian Li
# Email: lidian@zju.edu.cn
# --------------
# utils functions for rra data post process

# ========================================================================= #
# ========================================================================= #
ReadRRA = function (gene_summary, score = c("lfc", "rra")[1]) 
{
  if (is.null(dim(gene_summary))) {
    gene_summary = read.table(file = gene_summary, sep = "\t", 
                              header = TRUE, quote = "", comment.char = "", check.names = FALSE, 
                              stringsAsFactors = FALSE)
  }
  if (all(c("id", "Score", "FDR") %in% colnames(gene_summary))) {
    dd = as.data.frame(gene_summary[, c("id", "Score", "FDR")], 
                       stringsAsFactors = FALSE)
    dd$id = as.character(dd$id)
    return(dd)
  }
  gene_summary = gene_summary[, c(1, 3, 9, 8, 14, 5, 11)]
  colnames(gene_summary) = c("id", "negscore", "poscore", "neglfc", 
                             "poslfc", "negfdr", "posfdr")
  dd = gene_summary
  if ("lfc" %in% tolower(score)) {
    dd$LFC = dd$poslfc
    dd$FDR = dd$posfdr
    dd$LFC[abs(dd$neglfc) > dd$poslfc] = dd$neglfc[abs(dd$neglfc) > 
                                                     dd$poslfc]
    dd$FDR[abs(dd$neglfc) > dd$poslfc] = dd$negfdr[abs(dd$neglfc) > 
                                                     dd$poslfc]
    dd = dd[, c("id", "LFC", "FDR")]
  }
  else if ("rra" %in% tolower(score)) {
    idx_neg = dd$negscore < dd$poscore
    dd$LFC = apply(-log10(dd[, 2:3]), 1, max)
    dd$LFC[idx_neg] = -dd$LFC[idx_neg]
    dd$FDR = dd$posfdr
    dd$FDR[idx_neg] = dd$negfdr[idx_neg]
    dd = dd[, c("id", "LFC", "FDR")]
  }
  colnames(dd) = c("id", "Score", "FDR")
  dd$id = as.character(dd$id)
  return(dd)
}

# ========================================================================= #
# ========================================================================= #
# customized RankView to replace parameters of original RankView in MAGeCKFlute
RankView = function (rankdata, genelist = NULL, top = 10, bottom = 10, cutoff = NULL, 
                     main = NULL, filename = NULL, width = 5, height = 4, ...) 
{
  requireNamespace("ggrepel", quietly = TRUE) || stop("need ggrepel package")
  if (length(cutoff) == 0) 
    cutoff = CutoffCalling(rankdata, 1)
  if (length(cutoff) == 1) 
    cutoff = sort(c(-cutoff, cutoff))
  data = data.frame(Gene = names(rankdata), diff = rankdata, stringsAsFactors = FALSE)
  data$Rank = rank(data$diff)
  data$group = "no"
  data$group[data$diff > cutoff[2]] = "up"
  data$group[data$diff < cutoff[1]] = "down"
  idx = (data$Rank <= bottom) | (data$Rank > (max(data$Rank) - top)) | (data$Gene %in% genelist)
  mycolour = c(no = "gray80", up = "#e41a1c", down = "#377eb8")
  p = ggplot(data)
  p = p + geom_jitter(aes_string(x = "diff", y = "Rank", color = "group"), size = 0.1)
  if (!all(cutoff == 0)) 
    p = p + geom_vline(xintercept = cutoff, linetype = "dotted")
  if (sum(idx) > 0) 
    p = p + geom_label_repel(aes_string(x = "diff", y = "Rank", fill = "group", label = "Gene"), 
                             data = data[idx,], fontface = "bold", color = "white", size = 1.5, 
                             box.padding = unit(0.15, "lines"), segment.color = "grey50",
                             label.padding = unit(0.1, "lines"),
                             point.padding = unit(0.2, "lines"), segment.size = 0.2)
  p = p + scale_color_manual(values = mycolour)
  p = p + scale_fill_manual(values = mycolour)
  p = p + theme(panel.background = element_rect(fill = "white", colour = "black"))
  p = p + theme(text = element_text(colour = "black", size = 6, family = "Helvetica"), 
                plot.title = element_text(hjust = 0.5, size = 6), axis.text = element_text(colour = "gray10"))
  p = p + theme(axis.line = element_line(size = 0.5, colour = "black"), 
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                panel.border = element_blank(), panel.background = element_blank())
  # p = p + labs(x = "Treatment-Control beta score", y = "Rank", 
  #              title = main)
  p = p + labs(x = "LFC", y = "Rank", title = main)
  p = p + theme(legend.position = "none")
  if (!is.null(filename)) {
    ggsave(plot = p, filename = filename, units = "in", width = width, height = height, ...)
  }
  return(p)
}


# ========================================================================= #
# ========================================================================= #
# customized VolcanoView to replace parameters of original VolcanoView in MAGeCKFlute
VolcanoView = function (df, x = "logFC", y = "adj.P.Val", Label = NA, top = 5, 
          topnames = NULL, x_cutoff = log2(1.5), y_cutoff = 0.05, mycolour = c("gray80", 
                                                                               "#e41a1c", "#377eb8"), alpha = 0.6, force = 0.1, main = NULL, 
          xlab = "Log2 Fold Change", ylab = "-Log10(Adjust.P)", filename = NULL, 
          width = 4, height = 2.5, ...) 
{
  requireNamespace("ggrepel", quietly = TRUE) || stop("need ggrepel package")
  gg = df[, c(x, y)]
  gg$group = "no"
  gg$group[gg[, x] > x_cutoff & gg[, y] < y_cutoff] = "up"
  gg$group[gg[, x] < -x_cutoff & gg[, y] < y_cutoff] = "down"
  gg[, y] = -log10(gg[, y])
  if (!(top == 0 & is.null(topnames))) {
    gg$Label = rownames(gg)
    if (!is.na(Label)) 
      gg$Label = df[, Label]
    gg = gg[order(gg[, y], abs(gg[, x]), decreasing = TRUE), 
            ]
    idx1 = idx2 = c()
    if (top > 0) {
      idx1 = which(gg$group == "up")[1:min(top, sum(gg$group == "up"))]
      idx2 = which(gg$group == "down")[1:min(top, sum(gg$group == "down"))]
    }
    idx = unique(c(idx1, idx2, which(gg$Label %in% topnames)))
    gg$Label = as.character(gg$Label)
    gg$Label[setdiff(1:nrow(gg), idx)] = ""
  }
  gg$color = gg$group
  gg$color[gg$Label != ""] = "black"
  mycolour = c(mycolour, "black")
  names(mycolour) = c("no", "up", "down", "black")
  p = ggplot(gg, aes_string(x = x, y = y, label = "Label"))
  p = p + geom_point(aes_string(fill = "group"), shape = 21, alpha = alpha, show.legend = FALSE)
  p = p + geom_point(aes_string(colour = "color"), shape = 21, alpha = alpha, show.legend = FALSE)
  p = p + scale_color_manual(values = mycolour)
  p = p + scale_fill_manual(values = mycolour)
  p = p + theme(text = element_text(colour = "black", size = 6, family = "Helvetica"), 
                plot.title = element_text(hjust = 0.5, size = 6), axis.text = element_text(colour = "gray10"))
  p = p + theme(axis.line = element_line(size = 0.5, colour = "black"), 
                panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                panel.border = element_blank(), panel.background = element_blank())
  p = p + geom_hline(yintercept = -log10(y_cutoff), linetype = "dotted")
  p = p + geom_vline(xintercept = c(-x_cutoff, x_cutoff), linetype = "dotted")
  p = p + labs(x = xlab, y = ylab, title = main)
  if (!(top == 0 & is.null(topnames))) 
    p = p + ggrepel::geom_text_repel(force = force, fontface = "bold", 
                                     size = 1.5, segment.color = "grey50", segment.size = 0.2)
  p = p + theme(legend.position = "none")
  if (!is.null(filename)) {
    ggsave(plot = p, filename = filename, width = width, 
           height = height, units = "in", ...)
  }
  return(p)
}
# ========================================================================= #
# ========================================================================= #
TransGeneID = function (genes, fromType = "Symbol", toType = "Entrez", organism = "hsa", 
          fromOrg = organism, toOrg = organism, ensemblHost = "www.ensembl.org", 
          unique = TRUE, update = FALSE) 
{
  genes = as.character(genes)
  fromType = tolower(fromType)
  toType = tolower(toType)
  if (length(genes) < 1) 
    return(c())
  keggcode = rep(c("hsa", "mmu", "rno", "bta", "cfa", "ptr", 
                   "ssc"), 2)
  names(keggcode) = c(tolower(c("Human", "Mouse", "Rat", "Bovine", 
                                "Canine", "Chimp", "Pig")), c("hsa", "mmu", "rno", "bta", 
                                                              "cfa", "ptr", "ssc"))
  if (!tolower(organism) %in% names(keggcode)) 
    stop("Organism error ...")
  if (!tolower(fromOrg) %in% names(keggcode)) 
    stop("fromOrg error ...")
  if (!tolower(toOrg) %in% names(keggcode)) 
    stop("toOrg error ...")
  organism = keggcode[tolower(organism)]
  fromOrg = keggcode[tolower(fromOrg)]
  toOrg = keggcode[tolower(toOrg)]
  datasets = paste0(c("hsapiens", "mmusculus", "btaurus", "cfamiliaris", 
                      "ptroglodytes", "rnorvegicus", "sscrofa"), "_gene_ensembl")
  if (fromOrg == toOrg) {
    if (all(c(fromType, toType) %in% c("entrez", "symbol", "ensembl"))) {
      ann <- getGeneAnn(organism, update = update)$Gene
      if ("symbol" %in% c(fromType, toType) & any(!genes %in% ann$symbol)) {
        ann = rbind(ann, ann)
        idx = (nrow(ann)/2 + 1):nrow(ann)
        ann$symbol[idx] = ann$synonyms[idx]
      }
      ann = ann[, c(fromType, toType)]
    }
    else if (all(c(fromType, toType) %in% c("uniprot", "ensembl", "refseq", "symbol"))) {
      ann <- getGeneAnn(organism, update = update)$Protein
      ann = ann[, c(fromType, toType)]
    }
    else {
      ds = datasets[grepl(organism, datasets)]
      ensembl <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = ds, host = ensemblHost)
      attrs = biomaRt::listAttributes(ensembl)$name
      if (sum(attrs == fromType) == 0) {
        idx1 = grepl(tolower(fromType), attrs)
        idx = idx1
        if (sum(idx1) > 2) 
          idx = idx1 & grepl("_id", attrs)
        fromType = ifelse(sum(idx) > 0, attrs[idx][1], 
                          attrs[idx1][1])
        if (fromType == "hgnc_symbol" & fromOrg == "mmu") 
          fromType = "mgi_symbol"
      }
      if (sum(attrs == toType) == 0) {
        idx1 = grepl(tolower(toType), attrs)
        idx = idx1
        if (sum(idx1) > 2) 
          idx = idx1 & grepl("_id", attrs)
        toType = ifelse(sum(idx) > 0, attrs[idx][1], attrs[idx1])
        if (toType == "hgnc_symbol" & toOrg == "mmu") 
          toType = "mgi_symbol"
      }
      ann = biomaRt::getBM(attributes = c(fromType, toType), 
                           mart = ensembl, filters = fromType, values = genes)
    }
    idx = ann[, toType] == "" | is.na(ann[, toType])
    ann = ann[!idx, ]
    idx = ann[, fromType] == "" | is.na(ann[, fromType])
    ann = ann[!idx, ]
    tmp = ann
    tmp[, fromType] = gsub("\\..*|-.*", "", tmp[, fromType])
    ann = rbind.data.frame(ann, tmp)
    ann = ann[ann[, fromType] %in% genes, ]
  }
  else {
    if (all(c(fromType, toType) %in% c("symbol", "entrez"))) {
      ann = getOrtAnn(fromOrg, toOrg, update)
      ann = ann[, c(paste0(fromOrg, "_", fromType), paste0(toOrg, 
                                                           "_", toType))]
      colnames(ann) = c(fromOrg, toOrg)
    }
    else {
      from = biomaRt::useMart("ensembl", dataset = datasets[grepl(fromOrg, 
                                                                  datasets)])
      to = biomaRt::useMart("ensembl", dataset = datasets[grepl(toOrg, 
                                                                datasets)])
      attrs_1 = biomaRt::listAttributes(from)$name
      attrs_2 = biomaRt::listAttributes(to)$name
      if (sum(attrs_1 == fromType) == 0) {
        idx1 = grepl(tolower(fromType), attrs_1)
        idx = idx1
        if (sum(idx1) > 2) 
          idx = idx1 & grepl("_id", attrs_1)
        fromType = ifelse(sum(idx) > 0, attrs_1[idx][1], 
                          attrs_1[idx1][1])
        if (fromType == "hgnc_symbol" & fromOrg == "mmu") 
          fromType = "mgi_symbol"
      }
      if (sum(attrs_2 == toType) == 0) {
        idx1 = grepl(tolower(toType), attrs_2)
        idx = idx1
        if (sum(idx1) > 2) 
          idx = idx1 & grepl("_id", attrs_2)
        toType = ifelse(sum(idx) > 0, attrs_2[idx][1], 
                        attrs_2[idx1])
        if (toType == "hgnc_symbol" & toOrg == "mmu") 
          toType = "mgi_symbol"
      }
      ann = biomaRt::getLDS(attributes = fromType, mart = from, 
                            filters = fromType, values = genes, attributesL = toType, 
                            martL = to)
      colnames(ann) = c(fromOrg, toOrg)
    }
    idx = ann[, toOrg] == "" | is.na(ann[, toOrg])
    ann = ann[!idx, ]
    ann = ann[ann[, fromOrg] %in% genes, ]
  }
  ann = ann[!duplicated(paste0(ann[, 1], ann[, 2])), ]
  if (unique) {
    ann = ann[!duplicated(ann[, 1]), ]
    rownames(ann) = ann[, 1]
    tmp = ann[genes, 2]
    names(tmp) = genes
    ann = tmp
  }
  return(ann)
}
# ========================================================================= #
# ========================================================================= #
getOrtAnn = function (fromOrg = "mmu", toOrg = "hsa", update = FALSE) 
{
  rdsann = file.path(system.file("extdata", package = "MAGeCKFlute"), 
                     paste0("HOM_GeneID_Annotation_", fromOrg, "_", toOrg, 
                            ".rds"))
  if (file.exists(rdsann) & !update) 
    return(readRDS(rdsann))
  keggcode = c("hsa", "mmu", "rno", "bta", "cfa", "ptr", "ssc")
  names(keggcode) = c("human", "mouse", "rat", "bovine", "canine", 
                      "chimp", "pig")
  locfname <- file.path(system.file("extdata", package = "MAGeCKFlute"), 
                        "HOM_MouseHumanSequence.rpt.gz")
  if ((!file.exists(locfname)) | update) {
    refname <- "http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt"
    download.file(refname, locfname, quiet = TRUE)
  }
  mgi_ann <- read.table(locfname, sep = "\t", header = TRUE, 
                        stringsAsFactors = FALSE)
  mgi_ann = mgi_ann[, c(1, 2, 4, 5)]
  colnames(mgi_ann) = c("homoloid", "org", "symbol", "entrez")
  mgi_ann$org = gsub(", laboratory", "", mgi_ann$org)
  mgi_ann$org = keggcode[mgi_ann$org]
  locfname <- file.path(system.file("extdata", package = "MAGeCKFlute"), 
                        "homologene.data.gz")
  if ((!file.exists(locfname)) | update) {
    refname <- "ftp://ftp.ncbi.nih.gov/pub/HomoloGene/current/homologene.data"
    download.file(refname, locfname, quiet = TRUE)
  }
  ncbi_ann <- read.table(locfname, sep = "\t", stringsAsFactors = FALSE, 
                         quote = "")
  ncbi_ann = ncbi_ann[, c(1, 2, 4, 3)]
  colnames(ncbi_ann) = c("homoloid", "org", "symbol", "entrez")
  names(keggcode) = c(9606, 10090, 10116, 9913, 9615, 9598, 
                      9823)
  ncbi_ann$org = keggcode[as.character(ncbi_ann$org)]
  ncbi_ann = ncbi_ann[!is.na(ncbi_ann$org), ]
  ann = rbind.data.frame(mgi_ann, ncbi_ann)
  genes = unique(ann$entrez)
  idx1 = ann$entrez %in% genes
  idx2 = ann$org == toOrg
  idx3 = ann$homoloid %in% ann$homoloid[idx1]
  tmp1 = ann[idx1, c("homoloid", "symbol", "entrez")]
  tmp2 = ann[(idx2 & idx3), c("homoloid", "symbol", "entrez")]
  colnames(tmp1)[2:3] = paste0(fromOrg, c("_symbol", "_entrez"))
  colnames(tmp2)[2:3] = paste0(toOrg, c("_symbol", "_entrez"))
  ann = merge(tmp1, tmp2, by = "homoloid")[, -1]
  datasets = paste0(c("hsapiens", "mmusculus", "btaurus", "cfamiliaris", 
                      "ptroglodytes", "rnorvegicus", "sscrofa"), "_gene_ensembl")
  from = biomaRt::useMart("ensembl", dataset = datasets[grepl(fromOrg, 
                                                              datasets)])
  to = biomaRt::useMart("ensembl", dataset = datasets[grepl(toOrg, 
                                                            datasets)])
  from_symbol <- ifelse(fromOrg == "mmu", "mgi_symbol", "hgnc_symbol")
  to_symbol <- ifelse(toOrg == "mmu", "mgi_symbol", "hgnc_symbol")
  ensembl_ann = biomaRt::getLDS(attributes = c(from_symbol, 
                                               "entrezgene_id"), mart = from, attributesL = c(to_symbol, 
                                                                                              "entrezgene_id"), martL = to)
  colnames(ensembl_ann) = c(paste0(fromOrg, c("_symbol", "_entrez")), 
                            paste0(toOrg, c("_symbol", "_entrez")))
  ann = rbind.data.frame(ann, ensembl_ann)
  idx = duplicated(paste(ann[, 1], ann[, 2], ann[, 3], ann[, 
                                                           4], sep = "_"))
  ann = ann[!idx, ]
  saveRDS(ann, rdsann)
  return(ann)
}
# ========================================================================= #
# ========================================================================= #
CutoffCalling = function (d, scale = 1) 
{
  param = 1
  if (is.logical(scale) & scale) {
    param = round(length(d)/20000, digits = 1)
  }
  else if (is.numeric(scale)) {
    param = scale
  }
  Control_mean = 0
  sorted_beta = sort(abs(d))
  temp = quantile(sorted_beta, 0.68)
  temp_2 = qnorm(0.84)
  cutoff = round(temp/temp_2, digits = 3)
  names(cutoff) = NULL
  cutoff = cutoff * param
  return(cutoff)
}
# ========================================================================= #
# ========================================================================= #
ScatterView = function (data, x = "x", y = "y", label = 0, model = c("none", "ninesquare", "volcano", "rank")[1], x_cut = NULL, y_cut = NULL, 
          slope = 1, intercept = NULL, auto_cut = FALSE, auto_cut_x = auto_cut, 
          auto_cut_y = auto_cut, auto_cut_diag = auto_cut, groups = NULL, 
          group_col = NULL, groupnames = NULL, label.top = TRUE, top = 0, 
          toplabels = NULL, display_cut = FALSE, color = NULL, shape = 16, 
          size = 1, alpha = 0.6, main = NULL, xlab = x, ylab = y, legend.position = "none", 
          ...) 
{
  requireNamespace("ggplot2", quietly = TRUE) || stop("need ggplot package")
  requireNamespace("ggrepel", quietly = TRUE) || stop("need ggrepel package")
  requireNamespace("ggpubr", quietly = TRUE) || stop("need ggpubr package")
  data = as.data.frame(data, stringsAsFactors = FALSE)
  data = data[!(is.na(data[, x]) | is.na(data[, y])), ]
  if (label == 0) 
    data$Label = rownames(data)
  else data$Label = as.character(data[, label])
  if (!is.null(groupnames)) 
    legend.position = "right"
  model = tolower(model)
  if (model == "ninesquare") {
    if (length(x_cut) == 0) 
      x_cut = c(-CutoffCalling(data[, x], 1.5), CutoffCalling(data[, x], 1.5))
    if (length(y_cut) == 0) 
      y_cut = c(-CutoffCalling(data[, y], 1.5), CutoffCalling(data[, y], 1.5))
    if (length(intercept) == 0) 
      intercept = c(-CutoffCalling(data[, y] - data[, x], 1.5), CutoffCalling(data[, y] - data[, x], 1.5))
  }
  if (model == "volcano") {
    if (length(x_cut) == 0) 
      x_cut = c(-CutoffCalling(data[, x], 1.5), CutoffCalling(data[, x], 1.5))
    if (length(y_cut) == 0) 
      y_cut = -log10(0.05)
  }
  if (model == "rank") {
    if (length(x_cut) == 0) 
      x_cut = c(-CutoffCalling(data[, x], 1.5), CutoffCalling(data[, x], 1.5))
  }
  if (model == "none") {
    if (auto_cut_x) 
      x_cut = c(-CutoffCalling(data[, x], 1.5), CutoffCalling(data[, x], 1.5))
    if (auto_cut_y) 
      y_cut = c(-CutoffCalling(data[, y], 1.5), CutoffCalling(data[, y], 1.5))
    if (auto_cut_diag) 
      intercept = c(-CutoffCalling(data[, y] - data[, x], 1.5), CutoffCalling(data[, y] - data[, x], 1.5))
  }
  avail_groups = c("topleft", "topright", "bottomleft", "bottomright", 
                   "midleft", "topcenter", "midright", "bottomcenter", "midcenter", 
                   "top", "mid", "bottom", "left", "center", "right", "none")
  mycolour = c("#1f78b4", "#fb8072", "#33a02c", "#ff7f00", 
               "#bc80bd", "#66c2a5", "#6a3d9a", "#fdb462", "#ffed6f", 
               "#e78ac3", "#fdb462", "#8da0cb", "#66c2a5", "#fccde5", 
               "#fc8d62", "#d9d9d9")
  names(mycolour) = avail_groups
  if (model == "ninesquare") 
    groups = c("midleft", "topcenter", "midright", "bottomcenter")
  if (model == "volcano") 
    groups = c("topleft", "topright")
  if (model == "rank") 
    groups = c("left", "right")
  groups = intersect(groups, avail_groups)
  if (length(x_cut) > 0) {
    idx1 = data[, x] < min(x_cut)
    idx2 = data[, x] > max(x_cut)
  }
  else {
    idx1 = NA
    idx2 = NA
  }
  if (length(y_cut) > 0) {
    idx3 = data[, y] < min(y_cut)
    idx4 = data[, y] > max(y_cut)
  }
  else {
    idx3 = NA
    idx4 = NA
  }
  if (length(intercept) > 0) {
    idx5 = data[, y] < slope * data[, x] + min(intercept)
    idx6 = data[, y] > slope * data[, x] + max(intercept)
  }
  else {
    idx5 = NA
    idx6 = NA
  }
  data$group = "none"
  for (gr in groups) {
    if (gr == "topleft") 
      idx = cbind(idx1, idx4, idx6)
    if (gr == "topcenter") 
      idx = cbind(!idx1, !idx2, idx4, idx6)
    if (gr == "topright") 
      idx = cbind(idx2, idx4, idx6)
    if (gr == "midleft") 
      idx = cbind(idx1, idx6, !idx3, !idx4)
    if (gr == "midcenter") 
      idx = cbind(!idx1, !idx2, !idx3, !idx4, !idx5, !idx6)
    if (gr == "midright") 
      idx = cbind(idx2, !idx3, !idx4, idx5)
    if (gr == "bottomleft") 
      idx = cbind(idx1, idx3, idx5)
    if (gr == "bottomcenter") 
      idx = cbind(!idx1, !idx2, idx3, idx5)
    if (gr == "bottomright") 
      idx = cbind(idx2, idx3, idx5)
    if (gr == "top") {
      if (length(y_cut) > 0 & length(intercept) > 0) 
        idx = idx4 & idx6
      else if (length(y_cut) > 0) 
        idx = idx4
      else idx = idx6
    }
    if (gr == "mid") 
      idx = (!idx3) & (!idx4)
    if (gr == "bottom") {
      if (length(y_cut) > 0 & length(intercept) > 0) 
        idx = idx3 & idx5
      else if (length(y_cut) > 0) 
        idx = idx3
      else idx = idx5
    }
    if (gr == "left") {
      if (length(x_cut) > 0 & length(intercept) > 0) 
        if (slope > 0) 
          idx = idx1 & idx6
        else idx = idx1 & idx5
        else if (length(x_cut) > 0) 
          idx = idx1
        else if (slope > 0) 
          idx = idx6
        else idx = idx5
    }
    if (gr == "center") 
      idx = (!idx1) & (!idx2)
    if (gr == "right") {
      if (length(x_cut) > 0 & length(intercept) > 0) 
        if (slope > 0) 
          idx = idx2 & idx5
        else idx = idx2 & idx6
        else if (length(x_cut) > 0) 
          idx = idx2
        else if (slope > 0) 
          idx = idx5
        else idx = idx6
    }
    if (is.null(ncol(idx))) {
      if (sum(!is.na(idx)) > 0) 
        data$group[idx] = gr
      else warning("No cutpoint for group:", gr)
    }
    else {
      idx = idx[, !is.na(idx[1, ])]
      if (is.null(ncol(idx))) 
        warning("No cutpoint for group:", gr)
      else if (ncol(idx) < 4 & gr == "midcenter") 
        warning("No cutpoint for group:", gr)
      else data$group[rowSums(idx) == ncol(idx)] = gr
    }
  }
  data$group = factor(data$group, levels = unique(c(groups, "none")))
  if (length(groupnames) != length(groups)) 
    groupnames = groups
  if (length(groups) > 0) 
    names(groupnames) = groups
  if (length(group_col) == length(groups)) 
    mycolour[groups] = group_col
  if (length(groups) == 0) 
    mycolour["none"] = "#FF6F61"
  data$rank = top + 1
  for (g in groups) {
    idx1 = data$group == g
    x_symb = 0
    y_symb = 0
    if (g == "topleft") {
      x_symb = 1
      y_symb = -1
    }
    if (g == "topcenter") {
      x_symb = 0
      y_symb = -1
    }
    if (g == "topright") {
      x_symb = -1
      y_symb = -1
    }
    if (g == "midleft") {
      x_symb = 1
      y_symb = 0
    }
    if (g == "midright") {
      x_symb = -1
      y_symb = 0
    }
    if (g == "bottomleft") {
      x_symb = 1
      y_symb = 1
    }
    if (g == "bottomcenter") {
      x_symb = 0
      y_symb = 1
    }
    if (g == "bottomright") {
      x_symb = -1
      y_symb = 1
    }
    if (g == "top") {
      x_symb = 0
      y_symb = -1
    }
    if (g == "bottom") {
      x_symb = 0
      y_symb = 1
    }
    if (g == "left") {
      x_symb = 1
      y_symb = 0
    }
    if (g == "right") {
      x_symb = -1
      y_symb = 0
    }
    tmp = data[, c(x, y)]
    tmp[, x] = (tmp[, x] - min(tmp[, x]))/(max(tmp[, x]) - min(tmp[, x]))
    tmp[, y] = (tmp[, y] - min(tmp[, y]))/(max(tmp[, y]) - min(tmp[, y]))
    data$rank[idx1] = rank((x_symb * tmp[, x] + y_symb * tmp[, y])[idx1])
  }
  data$rank[data$rank == 0] = Inf
  if (mode(toplabels) == "list") {
    data$Label[data$rank > top & !(data$Label %in% unlist(toplabels))] = ""
    data$group = data$Label
    if (length(toplabels) > 0) {
      tmp = stack(toplabels)
      tmp = tmp[!duplicated(tmp[, 1]), ]
      rownames(tmp) = tmp[, 1]
      data$group[data$group %in% tmp[, 1]] = as.character(tmp[data$group[data$group %in% tmp[, 1]], 2])
      data$group[!(data$group %in% tmp[, 2]) & data$group != ""] = "Top hits"
    }
  }
  else {
    data$Label[data$rank > top & !(data$Label %in% toplabels)] = ""
  }
  if (is.null(color)) 
    color = "group"
  gg = data
  gg = gg[order(gg[, color]), ]
  p = ggplot(gg, aes_string(x, y, label = "Label", color = color))
  if (all(c(shape, size) %in% colnames(gg))) 
    p = p + geom_point(aes_string(shape = shape, size = size), alpha = alpha)
  else if (shape %in% colnames(gg)) 
    p = p + geom_point(aes_string(shape = shape), size = size, alpha = alpha)
  else if (size %in% colnames(gg)) 
    p = p + geom_point(aes_string(size = size), shape = shape, alpha = alpha)
  else p = p + geom_point(size = size, shape = shape, alpha = alpha)
  if (color == "group") {
    if (mode(toplabels) != "list") 
      p = p + scale_color_manual(values = mycolour, labels = groupnames)
    else p = p + scale_color_manual(values = c("#d9d9d9", 
                                               "#fb8072", "#80b1d3", "#fdb462", "#bc80bd", "#b3de69", 
                                               "#bebada", "#8dd3c7", "#ffffb3", "#fccde5", "#ccebc5", 
                                               "#ffed6f"))
  }
  else if (color %in% colnames(gg)) {
    if (mode(gg[, color]) == "numeric") 
      p = p + scale_color_gradient2(low = "#377eb8", high = "#e41a1c", 
                                    midpoint = 0)
    else if (!"try-error" %in% class(try(col2rgb(gg[1, color]), silent = TRUE))) {
      mycolour = unique(gg[, color])
      names(mycolour) = mycolour
      p = p + scale_color_manual(values = mycolour)
    }
  }
  if (label.top) 
    p = p + ggrepel::geom_text_repel(size = 2, segment.size = 0.2)
  if (display_cut) {
    if (length(x_cut) > 0) 
      p = p + geom_vline(xintercept = x_cut, linetype = "dotted")
    if (length(y_cut) > 0) 
      p = p + geom_hline(yintercept = y_cut, linetype = "dotted")
    if (length(intercept) > 0) 
      p = p + geom_abline(slope = slope, intercept = intercept, linetype = "dotted")
  }
  p = p + labs(x = xlab, y = ylab, title = main, color = NULL)
  p = p + ggpubr::theme_pubr() + theme(plot.title = element_text(hjust = 0.5))
  p = p + theme(legend.position = legend.position,
                text = element_text(colour = "black", size = 6, family = "Helvetica"), 
                plot.title = element_text(hjust = 0.5, size = 6), axis.text = element_text(colour = "gray10")) 
  return(p)
}