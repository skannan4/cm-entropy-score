#Transcriptomic entropy enables quantification of cardiomyocyte maturation at single cell level
#Chulan Kwon Laboratory
#Primary author: Suraj Kannan
#February 18,2020
#Document: Figures

#All necessary import files will be in this folder
setwd("~/Documents/Research/Reference")
library(ggplot2)
library(reshape2)
library(Matrix)
library(DropletUtils)
library(grid)
library(stringr)
library(monocle)

#So initially, the figures here were made with my computer screen in mind. As it turns out, with svg formatting, it is better to start small and then go bigger, rather than the other way around. Thus, the size of everything needed to be changed. The fig_factor below allows for this - if set to 1, the figures will be appropriate for a computer screen/ppt. If set to 2.5, they will be appropriate for the manuscript.
fig_factor = 1

#####Fig 01
#This figure computes the raw Shannon entropy on our reference data. This uses the oldest version of our entropy function, which has no gene filtering, no subselection of top genes, and also includes none of our QC (filtering out poor quality cells, filtering out non-CMs etc). We define and run this function here (in all its original glory as written by Michael Farid) but remove it subsequently since we never use it again.

michael_entropy = function(expressionMatrix){
  #feel free to check out Michael's spotify at mfd50 for some fire playlists
  p=sweep(expressionMatrix,MARGIN=2,FUN="/",STATS=colSums(expressionMatrix));
  logp=log(p)
  logp[logp==-Inf]=0
  plogp=p*logp
  result=-1*colSums(plogp)
  result
}

fig1_ref_df = data.frame(row.names = colnames(kannan_ref_data), timepoint = combined_datasets[combined_datasets$data == "kannan_ref_data", ]$timepoint, entropy = michael_entropy(kannan_ref_data))

ggplot(fig1_ref_df, aes(x = timepoint, y = entropy, fill = timepoint)) + geom_violin(scale = "width", position = position_dodge(width = 1), lwd =  1/fig_factor) + geom_boxplot(position = position_dodge(width = 1), lwd = 1/fig_factor, outlier.size = 1.5/fig_factor) + xlab("Timepoint") + ylab(expression('Shannon Entropy'~italic(S))) + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.position = "none") 

rm(fig1_ref_df, michael_entropy)

#####Fig 02
#This figure compares entropy score for various mappings of both our reference data and the 10x chromium data. The purpose of this figure is to highlight broadly that there is consistency cross-mapping; however, for genomic mappings methods (we use zUMIs here as the example) or for transcriptomic mapping methods that are given pseudogenes (we used kallisto as an example here), we show that correction of mitochondrial pseudogenes is sufficient to correct differences in entropy. Currently, the mito correct function is hard-coded into master_entropy (oops); thus, here we temporarily define a master_entropy function that doesn't have the mito correct. However, note that this version still includes all the other aspects of the master_entropy function (top 1000 genes and gene filtering).

master_entropy_nomito = function(expressionMatrix, n = 1000, species = "mouse"){
  #This portion splits up large datasets to go easier on memory usage.
  if(ncol(expressionMatrix) > 5000){
    lst = split(colnames(expressionMatrix), (seq(ncol(expressionMatrix))-1) %/% 5000)
    result = unlist(lapply(lst, function(x){master_entropy(expressionMatrix[, x], n = n, species = species)}))
  }
  else{
    if(species == "mouse"){
      #expressionMatrix = mito_correct(expressionMatrix)
      #begone mito correct!
    }
    if(species == "human"){
      G_list = G_list_human
    }
    expressionMatrix = expressionMatrix[rownames(expressionMatrix) %in% G_list[G_list$use_for_entropy == TRUE, ]$symbol & !startsWith(rownames(expressionMatrix), "Rps") & !startsWith(rownames(expressionMatrix), "Rpl") & !startsWith(rownames(expressionMatrix), "RPS") & !startsWith(rownames(expressionMatrix), "RPL"), ]
    expressionMatrix = apply(expressionMatrix, 2, function(x) tail(sort(x), n))
    p=sweep(expressionMatrix,MARGIN=2,FUN="/",STATS=colSums(expressionMatrix));
    logp=log(p)
    logp[logp==-Inf]=0
    plogp=p*logp
    result=-1*colSums(plogp)
  }
  return(result)
}

kannan_ref_kalcr = readMM("Fig02_data/ref_cellranger.mtx")
rownames(kannan_ref_kalcr) = read.table("Fig02_data/ref_cellranger.barcodes.txt", as.is = TRUE)$V1
colnames(kannan_ref_kalcr) = read.table("Fig02_data/ref_cellranger.genes.txt", as.is = TRUE)$V1
kannan_ref_kalcr = kannan_ref_kalcr[colnames(kannan_ref_data), ]
colnames(kannan_ref_kalcr) = substr(colnames(kannan_ref_kalcr), 1, 18)
kannan_ref_kalcr = rename_genes(t(kannan_ref_kalcr))

kannan_ref_full = readMM("Fig02_data/ref_full.mtx")
rownames(kannan_ref_full) = read.table("Fig02_data/ref_full.barcodes.txt", as.is = TRUE)$V1
colnames(kannan_ref_full) = read.table("Fig02_data/ref_full.genes.txt", as.is = TRUE)$V1
kannan_ref_full = kannan_ref_full[colnames(kannan_ref_data), ]
colnames(kannan_ref_full) = substr(colnames(kannan_ref_full), 1, 18)
kannan_ref_full = rename_genes(t(kannan_ref_full))

kannan_ref_entropy_comp = data.frame(row.names = colnames(kannan_ref_data), zUMIs_uncorrected = master_entropy_nomito(kannan_ref_data), zUMIs_corrected = master_entropy(kannan_ref_data), kb.cellranger_uncorrected = master_entropy_nomito(kannan_ref_kalcr), kb.cellranger_corrected = master_entropy(kannan_ref_kalcr), kb.full_uncorrected = master_entropy_nomito(kannan_ref_full), kb.full_corrected = master_entropy(kannan_ref_full), cellnames = colnames(kannan_ref_data), timepoint = combined_datasets[combined_datasets$data == "kannan_ref_data", ]$timepoint)
kannan_ref_entropy_comp = melt(kannan_ref_entropy_comp, id.vars = c("cellnames", "timepoint"))
kannan_ref_entropy_comp$method = factor(sapply(strsplit(as.character(kannan_ref_entropy_comp$variable), "_"), "[[", 1), levels = c("zUMIs", "kb.full", "kb.cellranger"))
kannan_ref_entropy_comp$correction = factor(sapply(strsplit(as.character(kannan_ref_entropy_comp$variable), "_"), "[[", 2), levels = c("uncorrected", "corrected"))
kannan_ref_mito_comp = data.frame(row.names = colnames(kannan_ref_data), zUMIs_uncorrected = mito(kannan_ref_data, correct = FALSE), zUMIs_corrected = mito(kannan_ref_data), kb.cellranger_uncorrected = mito(kannan_ref_kalcr, correct = FALSE), kb.cellranger_corrected = mito(kannan_ref_kalcr), kb.full_uncorrected = mito(kannan_ref_full, correct = FALSE), kb.full_corrected = mito(kannan_ref_full), cellnames = colnames(kannan_ref_data), timepoint = combined_datasets[combined_datasets$data == "kannan_ref_data", ]$timepoint)
kannan_ref_mito_comp = melt(kannan_ref_mito_comp, id.vars = c("cellnames", "timepoint"))
kannan_ref_mito_comp$method = factor(sapply(strsplit(as.character(kannan_ref_mito_comp$variable), "_"), "[[", 1), levels = c("zUMIs", "kb.full", "kb.cellranger"))
kannan_ref_mito_comp$correction = factor(sapply(strsplit(as.character(kannan_ref_mito_comp$variable), "_"), "[[", 2), levels = c("uncorrected", "corrected"))

###Fig 02a
ggplot(kannan_ref_entropy_comp, aes(x = method, y = value, fill = correction)) + geom_violin(scale = "width", position = position_dodge(width = 1), lwd =  1/fig_factor) + geom_boxplot(position = position_dodge(width = 1), lwd = 1/fig_factor, outlier.size = 1.5/fig_factor) + xlab("Mapping Method") + ylab("Entropy Score") + labs(fill = "Correction") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor), strip.text = element_text(size = 18/fig_factor)) + facet_wrap(~timepoint, scales = "free_x")

###Fig 02b
ggplot(kannan_ref_mito_comp, aes(x = method, y = value, fill = correction)) + geom_violin(scale = "width", position = position_dodge(width = 1), lwd =  1/fig_factor) + geom_boxplot(position = position_dodge(width = 1), lwd = 1/fig_factor, outlier.size = 1.5/fig_factor) + xlab("Mapping Method") + ylab("Mitochondrial Percentage") + labs(fill = "Correction") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor), strip.text = element_text(size = 18/fig_factor)) + facet_wrap(~timepoint, scales = "free_x")

chromium_kalcr = readMM("Fig02_data/chrom_cellranger.mtx")
rownames(chromium_kalcr) = read.table("Fig02_data/chrom_cellranger.barcodes.txt", as.is = TRUE)$V1
colnames(chromium_kalcr) = read.table("Fig02_data/chrom_cellranger.genes.txt", as.is = TRUE)$V1
chromium_kalcr = chromium_kalcr[colnames(chromium_data), ]
colnames(chromium_kalcr) = substr(colnames(chromium_kalcr), 1, 18)
chromium_kalcr = rename_genes(t(chromium_kalcr))

chromium_full = readMM("Fig02_data/chrom_full.mtx")
rownames(chromium_full) = read.table("Fig02_data/chrom_full.barcodes.txt", as.is = TRUE)$V1
colnames(chromium_full) = read.table("Fig02_data/chrom_full.genes.txt", as.is = TRUE)$V1
chromium_full = chromium_full[colnames(chromium_data), ]
colnames(chromium_full) = substr(colnames(chromium_full), 1, 18)
chromium_full = rename_genes(t(chromium_full))

chromium_zumis = readRDS("Fig02_data/chrom_dgecounts.rds")
chromium_zumis = rename_genes(chromium_zumis$umicount$inex$all)
chromium_zumis = chromium_zumis[, colnames(chromium_data)]

chromium_entropy_comp = data.frame(row.names = colnames(chromium_data), zUMIs_uncorrected = master_entropy_nomito(chromium_zumis), zUMIs_corrected = master_entropy(chromium_zumis), kb.cellranger_uncorrected = master_entropy_nomito(chromium_kalcr), kb.cellranger_corrected = master_entropy(chromium_kalcr), kb.full_uncorrected = master_entropy_nomito(chromium_full), kb.full_corrected = master_entropy(chromium_full), cellranger_uncorrected = master_entropy_nomito(chromium_data), cellranger_corrected = master_entropy(chromium_data), cellnames = colnames(chromium_data))
chromium_entropy_comp = melt(chromium_entropy_comp, id.vars = c("cellnames"))
chromium_entropy_comp$method = factor(sapply(strsplit(as.character(chromium_entropy_comp$variable), "_"), "[[", 1), levels = c("zUMIs", "kb.full", "kb.cellranger", "cellranger"))
chromium_entropy_comp$correction = factor(sapply(strsplit(as.character(chromium_entropy_comp$variable), "_"), "[[", 2), levels = c("uncorrected", "corrected"))
chromium_mito_comp = data.frame(row.names = colnames(chromium_data), zUMIs_uncorrected = mito(chromium_zumis, correct = FALSE), zUMIs_corrected = mito(chromium_zumis), kb.cellranger_uncorrected = mito(chromium_kalcr, correct = FALSE), kb.cellranger_corrected = mito(chromium_kalcr), kb.full_uncorrected = mito(chromium_full, correct = FALSE), kb.full_corrected = mito(chromium_full), cellranger_uncorrected = mito(chromium_data, correct = FALSE), cellranger_corrected = mito(chromium_data), cellnames = colnames(chromium_data))
chromium_mito_comp = melt(chromium_mito_comp, id.vars = c("cellnames"))
chromium_mito_comp$method = factor(sapply(strsplit(as.character(chromium_mito_comp$variable), "_"), "[[", 1), levels = c("zUMIs", "kb.full", "kb.cellranger", "cellranger"))
chromium_mito_comp$correction = factor(sapply(strsplit(as.character(chromium_mito_comp$variable), "_"), "[[", 2), levels = c("uncorrected", "corrected"))

###Fig 02c
ggplot(chromium_entropy_comp, aes(x = method, y = value, fill = correction)) + geom_violin(scale = "width", position = position_dodge(width = 1), lwd =  1/fig_factor) + geom_boxplot(position = position_dodge(width = 1), lwd = 1/fig_factor, outlier.size = 1.5/fig_factor) + xlab("Mapping Method") + ylab("Entropy Score") + labs(fill = "Correction") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor))

###Fig 02d
ggplot(chromium_mito_comp, aes(x = method, y = value, fill = correction)) + geom_violin(scale = "width", position = position_dodge(width = 1), lwd =  1/fig_factor) + geom_boxplot(position = position_dodge(width = 1), lwd = 1/fig_factor, outlier.size = 1.5/fig_factor) + xlab("Mapping Method") + ylab("Mitochondrial Percentage") + labs(fill = "Correction") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor))

rm(chromium_entropy_comp, chromium_full, chromium_kalcr, chromium_mito_comp, chromium_zumis, kannan_ref_entropy_comp, kannan_ref_full, kannan_ref_kalcr, kannan_ref_mito_comp, master_entropy_nomito)

#####Fig 03
#This figure computes the entropy score for 4 datasets - two counts based, and two UMI (Dueck, Jia e9.5, Hill e10.5 [first 100 good cells], and Duan 10x [first 100 good cells]) at various subsamplings. We selected these studies to span a range of starting depths and methods. For each, we make two plots - the entropy score at each subsampling as well as the accuracy (computed as 1 - normalized error, where normalized error is (subsampled entropy - full entropy)/full entropy). We also compute the minimum depth that achieves 98% accuracy (selected as a random cutoff, since it corresponds to ~0.1 difference in entropy). Note that these figures are generated using our default entropy score, which uses 1000 genes as cutoff. This explains the increasing entropy as the number of genes in the subsampled datasets approaches 1000; this would shift if we selected a different gene cutoff, and entirely disappear if we had no cutoff (for example, with the raw Shannon entropy). Since the purpose of this figure was to look at the stability of our entropy score used elsewhere in the manuscript, we kept our defaults. Because of the characteristics of this curve, we selected the point of 98% accuracy a depth where the median genes is greater than 1000; we also only included cells whose number of expressed genes > 1000.

prop_calc = function(data, reads){
  prop = reads/colSums(data)
  prop[prop > 1] = 1
  return(prop)
}

#This code is general purpose; to test any dataset, just replace the data name in the first two lines, and the maximum depth in the third line.

good_cells = combined_datasets[combined_datasets$data == "dueck_data" & combined_datasets$good_cell == TRUE, ]$cellname
downsampling_data = dueck_data[, good_cells]
downsamplings = 10^seq(3, 7.4, 0.05)
downsampled_umis = data.frame(row.names = colnames(downsampling_data), cells = colnames(downsampling_data))
for(i in downsamplings){
  downsampled_umis[, paste("entropy", trunc(i, 0), sep = "_")] = master_entropy(downsampleMatrix(mito_correct(downsampling_data), prop = prop_calc(downsampling_data, i)))
  downsampled_umis[, paste("accuracy", trunc(i, 0), sep = "_")] = 1-abs((master_entropy(downsampling_data) - master_entropy(downsampleMatrix(mito_correct(downsampling_data), prop = prop_calc(downsampling_data, i))))/(master_entropy(downsampling_data)))
  downsampled_umis[, paste("genes", trunc(i, 0), sep = "_")] = colSums(downsampleMatrix(mito_correct(downsampling_data), prop = prop_calc(downsampling_data, i)) > 0)
}
downsampled_umis = melt(downsampled_umis, id.vars = "cells" )
downsampled_umis$grouping = sapply(strsplit(as.character(downsampled_umis$variable), "_"), "[[", 1)
downsampled_umis$downsampling = as.numeric(sapply(strsplit(as.character(downsampled_umis$variable), "_"), "[[", 2))
genes = vector()
for(i in rownames(downsampled_umis)){
  cell = downsampled_umis[i, ]$cells
  samp = downsampled_umis[i, ]$downsampling
  genes[i] = downsampled_umis[downsampled_umis$cells == cell & downsampled_umis$downsampling == samp & downsampled_umis$grouping == "genes", ]$value
}
downsampled_umis$genes = genes

#To make our plots, we want to include only cells with genes > 1000 (since our gene cutoff is 1000). However, if we plot only cells with genes > 1000, at the lowest subsamplings only cells with naturally higher entropy will dominate. Thus, we plot using the lowest depth such that the median genes expressed in all cells is >1000.
gene_count = tapply(downsampled_umis[downsampled_umis$grouping == "genes", ]$value, downsampled_umis[downsampled_umis$grouping == "genes", ]$downsampling, median)
gene_min = as.numeric(rownames(gene_count)[gene_count > 1000][1])

###Fig 03a
ggplot(downsampled_umis[downsampled_umis$grouping == "entropy" & downsampled_umis$genes > 1000 & downsampled_umis$downsampling > gene_min, ], aes(x = log10(downsampling), y = value, group = downsampling)) + geom_boxplot(lwd = 1/fig_factor, outlier.size = 1.5/fig_factor) + xlab(expression(log[10](Downsampling))) + ylab("Entropy Score") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor))

#We now wish to plot the accuracy, and identify the minimum depth such that >98% accuracy is achieved. Here, we will compute this using the 1000 gene cutoff as well.

accuracy_count = tapply(downsampled_umis[downsampled_umis$grouping == "accuracy" & downsampled_umis$genes > 1000 & downsampled_umis$downsampling > gene_min, ]$value, downsampled_umis[downsampled_umis$grouping == "accuracy" & downsampled_umis$genes > 1000 & downsampled_umis$downsampling > gene_min, ]$downsampling, median)
accuracy_min = as.numeric(rownames(accuracy_count)[accuracy_count > 0.98][1])

grob = grobTree(textGrob(paste("98% Accuracy at Depth >=", accuracy_min), x=0.35,  y=0.05, hjust=0,
                          gp=gpar(col="red", fontsize=24/fig_factor)))

###Fig 03b
ggplot(downsampled_umis[downsampled_umis$grouping == "accuracy" & downsampled_umis$genes > 1000 & downsampled_umis$downsampling > gene_min, ], aes(x = log10(downsampling), y = value, group = downsampling)) + geom_boxplot(lwd = 1/fig_factor, outlier.size = 1.5/fig_factor) + xlab(expression(log[10](Downsampling))) + ylab("Accuracy") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor)) + geom_hline(yintercept = 0.98, color = "red") + geom_vline(xintercept = log10(accuracy_min)) + annotation_custom(grob)

rm(downsampled_umis, downsampling_data, grob, accuracy_count, accuracy_min, cell, downsamplings, gene_count, gene_min, good_cells, i, samp, prop_calc)

#####Fig 04
#This figure computes entropy score for mouse embryonic stem cells from six different library preparation methods, taken from Ziegenhain et al. (Molecular Cell 2017). We selected this particular benchmarking paper because a) it includes one well-defined cell type; b) it includes a range of library preparation techniques (droplet, plate, chip/machine; counts, UMI); and c) the samples are all sequenced to high depth, thus eliminating differences in depth as a potential variable. The purpose of this figure is to demonstrate that for a defined celltype, our entropy score is relatively constant across multiple sequencing methods. Please note: since we didn't use cardiomyocytes here, we didn't compute an SCN score. Thus, we selected good cells using the other metrics (normally compiled into "good_cell").

ziegenhain_data = read.table("Fig04_data/ziegenhain_data.txt", as.is = TRUE, row.names = 1, header = TRUE)
ziegenhain_pheno = read.table("Fig04_data/ziegenhain_pheno.txt", as.is = TRUE, row.names = 1, header = TRUE)
ziegenhain_data = rename_genes(ziegenhain_data)
#The Ziegenhain data has *very* high expression of a gene called Gm42418. As discussed by Valentine Svensson on his blog, this is likely a mismapped rRNA similar to 18S. Given that it is likely a contaminant, we remove it here.
ziegenhain_data = ziegenhain_data[rownames(ziegenhain_data) != "Gm42418", ]
ziegenhain_temp = data_qc("ziegenhain_data", study = "Ziegenhain", timepoint_list = ziegenhain_pheno$batch, scn_calc = FALSE, species = "mouse", other_meta = ziegenhain_pheno$method)

###Fig 04a
ggplot(ziegenhain_temp[ziegenhain_temp$top5_norm < 1.3 & ziegenhain_temp$depth_norm > -0.5 & ziegenhain_temp$genes > 1000, ], aes(x = other_meta, y = entropy, fill = other_meta)) + geom_violin(scale = "width", position = position_dodge(width = 1), lwd =  1/fig_factor) + geom_boxplot(position = position_dodge(width = 1), lwd = 1/fig_factor, outlier.size = 1.5/fig_factor) + xlab("Sequencing Method") + ylab("Entropy Score") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.position = "none") 

###Fig 04b
#We also plot the mitochondrial percentage - CELseq2 appears to have a lower entropy than all the other datasets, but also a much higher mitochondrial percentage, suggesting that this data may be less optimal in quality.
ggplot(ziegenhain_temp[ziegenhain_temp$top5_norm < 1.3 & ziegenhain_temp$depth_norm > -0.5 & ziegenhain_temp$genes > 1000, ], aes(x = other_meta, y = mito, fill = other_meta)) + geom_violin(scale = "width", position = position_dodge(width = 1), lwd =  1/fig_factor) + geom_boxplot(position = position_dodge(width = 1), lwd = 1/fig_factor, outlier.size = 1.5/fig_factor) + xlab("Sequencing Method") + ylab("Mitochondrial Percentage") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.position = "none") 

rm(ziegenhain_data, ziegenhain_pheno, ziegenhain_temp)

#####Fig 05
#This figure plots the celltype classification scores for all of the cells in all datasets. To classify, we used SingleCellNet against the Tabula Muris reference (see code for the SCN calculation method). Here, we plot the score for the CM classification across timepoint, and color the cells by whether their highest classification was for CM or for another celltype. Plots are made for mouse and human in vivo as well as human directed differentiation.

###Fig 05a
ggplot(combined_datasets[combined_datasets$full_label == "mouse in vivo", ], aes(x = timepoint, y = cm_score, color = max_celltype == "cardiac muscle cell", group = study)) + geom_jitter(size = 1/fig_factor) + xlab("Timepoint") + ylab("SingleCellNet CM Score") + scale_color_discrete(name = "Celltype \nClassification", labels = c("Non-CM", "CM")) + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), axis.text.x = element_text(angle = 45), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 18/fig_factor), legend.key.size = unit(1.2/fig_factor, "lines"))

###Fig 05b
ggplot(combined_datasets[combined_datasets$full_label == "human in vivo", ], aes(x = timepoint, y = cm_score, color = max_celltype == "cardiac muscle cell", group = study)) + geom_jitter(size = 1/fig_factor) + xlab("Timepoint") + ylab("SingleCellNet CM Score") + scale_color_discrete(name = "Celltype \nClassification", labels = c("Non-CM", "CM")) + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), axis.text.x = element_text(angle = 45), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 18/fig_factor), legend.key.size = unit(1.2/fig_factor, "lines"))

###Fig 05c
ggplot(combined_datasets[combined_datasets$full_label == "human directed differentiation", ], aes(x = timepoint, y = cm_score, color = max_celltype == "cardiac muscle cell", group = study)) + geom_jitter(size = 1/fig_factor) + xlab("Timepoint") + ylab("SingleCellNet CM Score") + scale_color_discrete(name = "Celltype \nClassification", labels = c("Non-CM", "CM")) + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), axis.text.x = element_text(angle = 45), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 18/fig_factor), legend.key.size = unit(1.2/fig_factor, "lines"))

#####Fig 06
#This figure plots the percentage of counts going to ribosomal protein coding genes (e.g. Rps, Rpl) in the various mouse in vivo datasets. We plot both by timepoint and by sequencing protocol; we observe that the abundance of ribosomal protein coding genes appears higher in certain protocols over others. Please note - we focus on the ribosomal protein coding genes here because they are the ones that make it through the gene filtering; however, some datasets have other ribosomal contaminants (for example, rRNA [Rn45, mt-Rnr1, mt-Rnr2] or mismapped rRNAs) that we are not plotting here.

###Fig 06a
ggplot(combined_datasets[combined_datasets$full_label == "mouse in vivo", ], aes(x = timepoint, y = ribo, fill = study)) + geom_violin(scale = "width", position = position_dodge(width = 1), lwd =  1/fig_factor) + geom_boxplot(position = position_dodge(width = 1), lwd = 1/fig_factor, outlier.size = 1.5/fig_factor) + xlab("Timepoint") + ylab("Ribosomal Protein Coding Gene %") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), axis.text.x = element_text(angle = 45), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor), legend.key.size = unit(1.2/fig_factor, "lines")) + scale_fill_discrete(name = "Study") + guides(fill = guide_legend(ncol = 1))

###Fig 06b
ggplot(combined_datasets[combined_datasets$full_label == "mouse in vivo", ], aes(x = sequencing, y = ribo, fill = study)) + geom_violin(scale = "width", position = position_dodge(width = 1), lwd =  1/fig_factor) + geom_boxplot(position = position_dodge(width = 1), lwd = 1/fig_factor, outlier.size = 1.5/fig_factor) + xlab("Sequencing Protocol") + ylab("Ribosomal Protein Coding Gene %") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor), legend.key.size = unit(1.2/fig_factor, "lines")) + scale_fill_discrete(name = "Study") + guides(fill = guide_legend(ncol = 1))

#####Fig 07
#This figure highlights our QC approach - namely, removing cells based on low depth or high percentage of genes skewed towards the top 5 genes. Both are metrics of poor quality cells (and affect the entropy score). We compute our metrics per timepoint per dataset, since both can affect the appropriate benchmarking. Fig 07a-b show examples of how depth/top5% can affect entropy (we use the Churko data because it conveniently demonstrates both phenomena). Fig 07c-d show our metrics plotted for all the datasets as well as the cutoffs (note that we use very conservative cutoffs).

###Fig 07a
ggplot(combined_datasets[combined_datasets$data == "churko_data", ], aes(x = log10(depth), y = entropy, color = timepoint)) + geom_point(size = 1/fig_factor) + xlab(expression(log[10](Depth))) + ylab("Entropy Score") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor)) + scale_color_discrete(name = "Timepoint")

###Fig 07b
ggplot(combined_datasets[combined_datasets$data == "churko_data", ], aes(x = top5, y = entropy, color = timepoint)) + geom_point(size = 1/fig_factor) + xlab("Fraction Counts to Top 5 Genes") + ylab("Entropy Score") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor)) + scale_color_discrete(name = "Timepoint")

###Fig 07c
ggplot(combined_datasets[combined_datasets$max_celltype == "cardiac muscle cell" & combined_datasets$genes > 1000, ], aes(x = timepoint, y = depth_norm, fill = study)) + geom_boxplot(outlier.size = 1.5/fig_factor) + xlab("Timepoint") + ylab("Normalized Depth QC Metric") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), axis.text.x = element_text(size = 14/fig_factor,angle = 45), legend.title = element_text(size = 14/fig_factor), legend.text = element_text(size = 8/fig_factor), legend.key.size = unit(1.2/fig_factor, "lines")) + scale_fill_discrete(name = "Study") + guides(fill = guide_legend(ncol = 1)) + geom_hline(yintercept = -0.5, color = "red", size = 1.5/fig_factor)

###Fig 07d
ggplot(combined_datasets[combined_datasets$max_celltype == "cardiac muscle cell" & combined_datasets$genes > 1000, ], aes(x = timepoint, y = top5_norm, fill = study)) + geom_boxplot(outlier.size = 1.5/fig_factor) + xlab("Timepoint") + ylab("Normalized Top5% QC Metric") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), axis.text.x = element_text(size = 14/fig_factor,angle = 45), legend.title = element_text(size = 14/fig_factor), legend.text = element_text(size = 8/fig_factor), legend.key.size = unit(1.2/fig_factor, "lines")) + scale_fill_discrete(name = "Study") + guides(fill = guide_legend(ncol = 1)) + geom_hline(yintercept = 1.3, color = "red", size = 1.5/fig_factor)

#####Fig 08
#This figure plots the percentage of counts going to mitochondrial transcripts in the various mouse in vivo datasets. We plot per timepoint, including only cells that passed the quality filters (top5, depth, and CM classification filters). The purpose of this figure is two-fold - to show that the mitochondrial read percentage increases over CM development, and to highlight certain poor quality datasets. In the case of the latter, cells pass the previous QC filters because the entire dataset is affected (usually due to difficulty in isolating postnatal CMs). We identify these datasets manually, based on empirical observation of high mitochondrial percentage. The datasets on this criteria are outlined in red in Fig 08a; in Fig 08b they are removed. In Fig 08c, we show mitochondral read percentages only for datasets included in the final analysis (i.e. eliminating low depth studies).

###Fig 08a
ggplot(combined_datasets[combined_datasets$full_label == "mouse in vivo" & combined_datasets$good_cell == TRUE, ], aes(x = timepoint, y = mito, fill = study, color = (reason != "high mitochondrial percentage" | is.na(reason)))) + geom_violin(scale = "width", position = position_dodge(width = 1), lwd =  1/fig_factor) + geom_boxplot(position = position_dodge(width = 1), lwd = 1/fig_factor, outlier.size = 1.5/fig_factor) + xlab("Timepoint") + ylab("Mito %") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), axis.text.x = element_text(angle = 45), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor),  legend.key.size = unit(1.2/fig_factor, "lines")) + scale_fill_discrete(name = "Study") + guides(fill = guide_legend(ncol = 1)) + scale_color_manual(guide = FALSE, values = c("red", "black"))

###Fig 08b
ggplot(combined_datasets[combined_datasets$full_label == "mouse in vivo" & combined_datasets$good_cell == TRUE & (combined_datasets$reason != "high mitochondrial percentage" | is.na(combined_datasets$reason)), ], aes(x = timepoint, y = mito, fill = study)) + geom_violin(scale = "width", position = position_dodge(width = 1), lwd =  1/fig_factor) + geom_boxplot(position = position_dodge(width = 1), lwd = 1/fig_factor, outlier.size = 1.5/fig_factor) + xlab("Timepoint") + ylab("Mito %") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), axis.text.x = element_text(angle = 45), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor), legend.key.size = unit(1.2/fig_factor, "lines")) + scale_fill_discrete(name = "Study") + guides(fill = guide_legend(ncol = 1))

###Fig 08c
ggplot(combined_datasets[combined_datasets$full_label == "mouse in vivo" & combined_datasets$good_cell == TRUE & combined_datasets$include_dataset == TRUE & combined_datasets$genes > 1000, ], aes(x = timepoint, y = mito, fill = study)) + geom_violin(scale = "width", position = position_dodge(width = 1), lwd =  1/fig_factor) + geom_boxplot(position = position_dodge(width = 1), lwd = 1/fig_factor, outlier.size = 1.5/fig_factor) + xlab("Timepoint") + ylab("Mito %") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), axis.text.x = element_text(angle = 45), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor), legend.key.size = unit(1.2/fig_factor, "lines")) + scale_fill_discrete(name = "Study") + guides(fill = guide_legend(ncol = 1))

#####Fig 09
#The purpose of this figure is to compare the entropy scores of our reference data to three pseudotime methods - the well-known Monocle 2 as well Slingshot and SCORPIUS, the latter two which have scored notably well for linear trajectories (see Saelens and Cannoodt, Nature Biotech 2019). For each method, we compute and compare pseudotimes; we then use Monocle's differential gene approach, simply substituting in the corresponding pseudotime values. Note - we are aware that Monocle 2 is deprecated; however, for simple trajectories such as the one in this dataset, DDRTree is more than sufficient (and franky a little easier to use) - hence our choice of Monocle 2 over 3.

library(slingshot)
library(SingleCellExperiment)
library(mclust)
library(RColorBrewer)
library(SCORPIUS)
library(VennDiagram)

#used for everything
good_cells = combined_datasets[combined_datasets$data == "kannan_ref_data" & combined_datasets$good_cell == TRUE & combined_datasets$genes > 1000, ]$cellname
trajectory_data = mito_correct(kannan_ref_data[, good_cells])
trajectory_pheno = combined_datasets[combined_datasets$data == "kannan_ref_data" & combined_datasets$good_cell == TRUE & combined_datasets$genes > 1000, ]
rownames(trajectory_pheno) = trajectory_pheno$cellname

#Monocle 2 - We used the dpFeature approach to be as unbiased as possible
trajectory_monocle = newCellDataSet(trajectory_data, phenoData = new("AnnotatedDataFrame", data = trajectory_pheno), featureData = new("AnnotatedDataFrame", data = data.frame(row.names = rownames(trajectory_data), gene1 = rownames(trajectory_data), gene2 = rownames(trajectory_data))), lowerDetectionLimit = 0.5, expressionFamily = negbinomial.size())
trajectory_monocle = estimateSizeFactors(trajectory_monocle)
trajectory_monocle = estimateDispersions(trajectory_monocle)
trajectory_monocle = detectGenes(trajectory_monocle, min_expr = 0.1)
fData(trajectory_monocle)$use_for_ordering = fData(trajectory_monocle)$num_cells_expressed > 0.05 * ncol(trajectory_monocle)
trajectory_monocle = reduceDimension(trajectory_monocle,
                            max_components = 2,
                            norm_method = 'log',
                            num_dim = 3,
                            reduction_method = 'tSNE',
                            verbose = T)
trajectory_monocle <- clusterCells(trajectory_monocle, verbose = F)
trajectory_monocle <- detectGenes(trajectory_monocle, min_expr = 0.1)
#selected only genes expressed in 10% of the cells; we'll reuse this for all the methods
expressed_genes = row.names(subset(fData(trajectory_monocle), num_cells_expressed >= 67))
trajectory_monocle = trajectory_monocle[expressed_genes, ]
clustering_DEG_genes = differentialGeneTest(trajectory_monocle, fullModelFormulaStr = '~Cluster', cores = 1)
ordering_genes = row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]
trajectory_monocle = setOrderingFilter(trajectory_monocle, ordering_genes = ordering_genes)
trajectory_monocle = reduceDimension(trajectory_monocle, method = 'DDRTree')
trajectory_monocle = orderCells(trajectory_monocle, reverse = TRUE)

###Fig 09a
plot_cell_trajectory(trajectory_monocle, color_by = "timepoint", cell_size = 1.5/fig_factor) + xlab("Component 1") + ylab("Component 2") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor)) + scale_color_discrete(name = "Timepoint")

#Slingshot - using PCA
trajectory_slingshot = SingleCellExperiment(assays = List(counts = trajectory_data[expressed_genes, ]), colData = trajectory_pheno)
FQnorm = function(counts){
  rk = apply(counts,2,rank,ties.method='min')
  counts.sort = apply(counts,2,sort)
  refdist = apply(counts.sort,1,median)
  norm = apply(rk,2,function(r){ refdist[r] })
  rownames(norm) = rownames(counts)
  return(norm)
}
assays(trajectory_slingshot)$norm = FQnorm(assays(trajectory_slingshot)$counts)
pca = prcomp(t(log1p(assays(trajectory_slingshot)$norm)), scale. = FALSE)
rd1 = pca$x[,1:2]
reducedDims(trajectory_slingshot) = SimpleList(PCA = rd1)
cl1 = Mclust(rd1)$classification
colData(trajectory_slingshot)$GMM = cl1
trajectory_slingshot <- slingshot(trajectory_slingshot, clusterLabels = 'GMM', reducedDim = 'PCA')
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(12)
names(colors) = levels(trajectory_slingshot$timepoint)[levels(trajectory_slingshot$timepoint) %in% trajectory_slingshot$timepoint]
colors2 = colors[as.character(trajectory_slingshot$timepoint)]

###Fig 09b
#Would I be happier if this weren't in base R? Of course, but I'm not going through the headache of figuring it out.
plot(reducedDims(trajectory_slingshot)$PCA, col = colors2, pch=16, asp = 1)
legend("bottomleft",legend = names(colors), fill = colors)
lines(SlingshotDataSet(trajectory_slingshot), lwd=2/fig_factor, col = c("black"))

#SCORPIUS
space <- reduce_dimensionality(t(trajectory_data[expressed_genes, ]), "spearman", ndim = 3)
trajectory_scorpius <- infer_trajectory(space)

###Fig 09c
#Please see above point about Base R
draw_trajectory_plot(
  space, 
  progression_group = trajectory_pheno$timepoint,
  progression_group_palette = colors2,
  path = trajectory_scorpius$path,
  contour = FALSE
)

#Plotting and computing the correlations
trajectory_pheno$monocle = trajectory_monocle$Pseudotime
trajectory_pheno$slingshot = trajectory_slingshot$slingPseudotime_1
trajectory_pheno$scorpius = trajectory_scorpius$time

###Used for Fig 09d
cor(trajectory_pheno$entropy, trajectory_pheno$monocle)
cor(trajectory_pheno$entropy, trajectory_pheno$slingshot)
cor(trajectory_pheno$entropy, trajectory_pheno$scorpius)

###Fig 09e
ggplot(trajectory_pheno, aes(x = monocle, y = entropy, color = timepoint)) + geom_point(size = 1/fig_factor) + xlab("Monocle Pseudotime") + ylab("Entropy Score") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor)) + scale_color_discrete(name = "Timepoint")

###Fig 09f
ggplot(trajectory_pheno, aes(x = slingshot, y = entropy, color = timepoint)) + geom_point(size = 1/fig_factor) + xlab("Slingshot Pseudotime") + ylab("Entropy Score") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor)) + scale_color_discrete(name = "Timepoint")

###Fig 09g
ggplot(trajectory_pheno, aes(x = scorpius, y = entropy, color = timepoint)) + geom_point(size = 1/fig_factor) + xlab("Scorpius Pseudotime") + ylab("Entropy Score") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 18/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor)) + scale_color_discrete(name = "Timepoint")

#We compute the differentially expressed genes based on each of the pseudotimes outputted by each method (as well as for entropy). Most of these linear pseudotime approaches use a similar approach to compute differential genes (some type of GAM fitting); to simplify, we use Monocle 2's built-in differential gene calculator, just replacing pseudotime with each value.

pData(trajectory_monocle)$monocle = trajectory_monocle$Pseudotime
pData(trajectory_monocle)$slingshot = trajectory_slingshot$slingPseudotime_1
pData(trajectory_monocle)$scorpius = trajectory_scorpius$time
 
diff_monocle = differentialGeneTest(trajectory_monocle, fullModelFormulaStr = "~sm.ns(Pseudotime, df=3)")
diff_slingshot = differentialGeneTest(trajectory_monocle, fullModelFormulaStr = "~sm.ns(slingshot, df=3)")
diff_scorpius = differentialGeneTest(trajectory_monocle, fullModelFormulaStr = "~sm.ns(scorpius, df=3)")
diff_entropy = differentialGeneTest(trajectory_monocle, fullModelFormulaStr = "~sm.ns(entropy, df=3)")

diff_monocle_genes = rownames(diff_monocle)[diff_monocle$qval < 0.05]
diff_slingshot_genes = rownames(diff_slingshot)[diff_slingshot$qval < 0.05]
diff_scorpius_genes = rownames(diff_scorpius)[diff_scorpius$qval < 0.05]
diff_entropy_genes = rownames(diff_entropy)[diff_entropy$qval < 0.05]

###Fig 09h (This one will save directly)
venn.diagram(x = list(diff_entropy_genes, diff_monocle_genes, diff_slingshot_genes, diff_scorpius_genes), category.names = c("Entropy \nScore", "Monocle 2", "Slingshot", "SCORPIUS"), filename = "FirstDraftFigures/Ref09h.tiff", output = TRUE, imagetype = "tiff", col = "black",
             lty = "dotted",
             lwd = 4/fig_factor,
             fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),
             alpha = 0.50,
             label.col = c("orange", "white", "darkorchid4", "white", "white", "red",
                           "white", "white", "darkblue", "white",
                           "white", "white", "white", "darkgreen", "white"),
             cex = 1.5,
             fontfamily = "sans",
             fontface = "bold",
             cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),
             cat.cex = 1.5,
             cat.fontfamily = "sans", width = 3000)

rm(clustering_DEG_genes, diff_entropy, diff_monocle, diff_scorpius, diff_slingshot, pca, rd1, space, cl1, colors, colors2, diff_entropy_genes, diff_monocle_genes, diff_scorpius_genes, diff_slingshot_genes, good_cells, ordering_genes, FQnorm, trajectory_data, trajectory_pheno, trajectory_scorpius, trajectory_slingshot, expressed_genes)
#We temporarily hold onto the trajectory_monocle object so we can use in another figure

#####Fig 10
#This is the $$$ figure. If I graduate with a PhD, it will almost entirely be on the merits of this figure. Months of my life were lost to make this one line of code. Here, we plot the entropy score for all of our in vivo datasets across timepoint. For supplement, we also plot them while coloured by other parameters (isolation, sequencing, mapping).

supp.labs = c("Mouse", "Human")
names(supp.labs) = c("mouse", "human")

###Fig 10a
ggplot(combined_datasets[combined_datasets$sample_type == "in vivo" & combined_datasets$good_cell == TRUE & combined_datasets$genes > 1000 & combined_datasets$include_dataset == TRUE, ], aes(x = timepoint, y = entropy, fill = study)) + geom_violin(scale = "width", position = position_dodge(width = 1), lwd =  1/fig_factor) + geom_boxplot(position = position_dodge(width = 1), lwd = 1/fig_factor, outlier.size = 1.5/fig_factor) + xlab("Timepoint") + ylab("Entropy Score") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 16/fig_factor), axis.text.x = element_text(angle = 45), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor), legend.key.size = unit(1.2/fig_factor, "lines"), strip.text = element_text(size = 16/fig_factor)) + scale_fill_discrete(name = "Study") + guides(fill = guide_legend(ncol = 1)) + facet_grid(~species, scales = "free_x", space = "free_x", labeller = labeller(species = supp.labs))

###Fig 10b
ggplot(combined_datasets[combined_datasets$sample_type == "in vivo" & combined_datasets$good_cell == TRUE & combined_datasets$genes > 1000 & combined_datasets$include_dataset == TRUE, ], aes(x = timepoint, y = entropy, dodge = study, fill = isolation)) + geom_violin(scale = "width", position = position_dodge(width = 1), lwd =  1/fig_factor) + geom_boxplot(position = position_dodge(width = 1), lwd = 1/fig_factor, outlier.size = 1.5/fig_factor) + xlab("Timepoint") + ylab("Entropy Score") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 16/fig_factor), axis.text.x = element_text(angle = 45), legend.key.size = unit(1.2/fig_factor, "lines"), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor), strip.text = element_text(size = 16/fig_factor)) + scale_fill_discrete(name = "Isolation \nMethod") + guides(fill = guide_legend(ncol = 1)) + facet_grid(~species, scales = "free_x", space = "free_x", labeller = labeller(species = supp.labs))

###Fig 10c
ggplot(combined_datasets[combined_datasets$sample_type == "in vivo" & combined_datasets$good_cell == TRUE & combined_datasets$genes > 1000 & combined_datasets$include_dataset == TRUE, ], aes(x = timepoint, y = entropy, dodge = study, fill = sequencing)) + geom_violin(scale = "width", position = position_dodge(width = 1), lwd =  1/fig_factor) + geom_boxplot(position = position_dodge(width = 1), lwd = 1/fig_factor, outlier.size = 1.5/fig_factor) + xlab("Timepoint") + ylab("Entropy Score") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 16/fig_factor), axis.text.x = element_text(angle = 45), legend.key.size = unit(1.2/fig_factor, "lines"), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor), strip.text = element_text(size = 16/fig_factor)) + scale_fill_discrete(name = "Sequencing \nProtocol") + guides(fill = guide_legend(ncol = 1)) + facet_grid(~species, scales = "free_x", space = "free_x", labeller = labeller(species = supp.labs))

###Fig 10d
ggplot(combined_datasets[combined_datasets$sample_type == "in vivo" & combined_datasets$good_cell == TRUE & combined_datasets$genes > 1000 & combined_datasets$include_dataset == TRUE, ], aes(x = timepoint, y = entropy, dodge = study, fill = mapping)) + geom_violin(scale = "width", position = position_dodge(width = 1), lwd =  1/fig_factor) + geom_boxplot(position = position_dodge(width = 1), lwd = 1/fig_factor, outlier.size = 1.5/fig_factor) + xlab("Timepoint") + ylab("Entropy Score") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 16/fig_factor), axis.text.x = element_text(angle = 45), legend.key.size = unit(1.2/fig_factor, "lines"), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor), strip.text = element_text(size = 16/fig_factor)) + scale_fill_discrete(name = "Mapping \nMethod") + guides(fill = guide_legend(ncol = 1)) + facet_grid(~species, scales = "free_x", space = "free_x", labeller = labeller(species = supp.labs))

rm(supp.labs)

#####Fig 11
#The purpose of this figure is to show gene expression patterns for well known genes in cardiac maturation across entropy. We selected four classes of genes - sarcomeric (Myh6, Tnni3, Myom2, Ttn), cell cycle (Ccnd1, Ccnb1, Cdk1, Cdk4), calcium handling (Casq2, Ryr2, Pln, S100a1), and metabolic/energetic (mt-Co1, mt-Co3, Atp5a1, Atp2a2, Ckmt2, Cox5b, Cox6a2, Slc2a1). We selected the particular genes based on past literature. We plotted the genes using Monocle 2's built-in features for simplicity, using Monocle's normalization for relative expression. We simply replace "Pseudotime" in Monocle 2 with -1 * Entropy (done so that genes can still be interpretted left-to-right over development). While the main figure uses our reference, we also include the Chromium, Goodyer PF, and Duan 10x datasets in the supplement.

#This code is general purpose - can replace with any dataset you want.
dataset = "kannan_ref_data"
good_cells = combined_datasets[combined_datasets$data == dataset & combined_datasets$good_cell == TRUE & combined_datasets$genes > 1000, ]$cellname
plot_data = mito_correct(get(dataset)[, good_cells])
plot_pheno = combined_datasets[combined_datasets$data == dataset & combined_datasets$good_cell == TRUE & combined_datasets$genes > 1000, ]
rownames(plot_pheno) = plot_pheno$cellname
plot_monocle = newCellDataSet(plot_data, phenoData = new("AnnotatedDataFrame", data = plot_pheno), featureData = new("AnnotatedDataFrame", data = data.frame(row.names = rownames(plot_data), gene1 = rownames(plot_data), gene2 = rownames(plot_data))), lowerDetectionLimit = 0.5, expressionFamily = negbinomial.size())
plot_monocle = estimateSizeFactors(plot_monocle)
plot_monocle = estimateDispersions(plot_monocle)
plot_monocle$Pseudotime = -1 * plot_monocle$entropy

###Fig 11a
plot_genes_in_pseudotime(plot_monocle[c("Myh6", "Tnni3", "Myom2", "Ttn", "Ccnd1", "Ccnb1", "Cdk1", "Cdk4", "Casq2", "Ryr2", "Pln", "S100a1", "mt-Co1", "mt-Co3", "Atp5a1", "Atp2a2", "Ckmt2", "Cox5b", "Cox6a2", "Slc2a1"), ], nrow = 5, ncol = 4, color = "timepoint") + xlab("-1 * Entropy Score") + theme_bw() + scale_color_discrete(name = "Timepoint")
#This one doesn't scale with fig_factor - too annoying for me to figure out

#Needed to make some cuts for the version in the paper, and reorganize a bit
plot_genes_in_pseudotime(plot_monocle[c("Myh6", "Tnni3", "Myom2", "Ryr2", "Ccnb1", "Cdk1", "Cdk4", "S100a1", "mt-Co1", "mt-Co3", "Casq2", "Atp2a2", "Ckmt2", "Cox6a2", "Slc2a1", "Pln"), ], nrow = 4, ncol = 4, color = "timepoint") + xlab("-1 * Entropy Score") + theme_bw() + scale_color_discrete(name = "Timepoint") + theme(strip.text.x = element_text(size = 18))

rm(dataset, good_cells, plot_data, plot_pheno, plot_monocle, trajectory_monocle)

#####Fig 12
#Another $$$ figure. This figure plots the entropy score for the human in vitro directed differentiation datasets; we plot along with the human in vivo datasets to enable direct comparison.

supp.labs = c("In Vivo", "Directed Differentiation")
names(supp.labs) = c("in vivo", "directed differentiation")

###Fig 12a
ggplot(combined_datasets[combined_datasets$species == "human" & combined_datasets$good_cell == TRUE & combined_datasets$genes > 1000 & combined_datasets$include_dataset == TRUE & !combined_datasets$timepoint %in% c("D0", "D2", "D5"), ], aes(x = timepoint, y = entropy, fill = study)) + geom_violin(scale = "width", position = position_dodge(width = 1), lwd =  1/fig_factor) + geom_boxplot(position = position_dodge(width = 1), lwd = 1/fig_factor, outlier.size = 1.5/fig_factor) + xlab("Timepoint") + ylab("Entropy Score") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 16/fig_factor), axis.text.x = element_text(angle = 45), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor), strip.text = element_text(size = 16/fig_factor), legend.key.size = unit(1.2/fig_factor, "lines")) + scale_fill_discrete(name = "Study") + guides(fill = guide_legend(ncol = 1)) + facet_grid(~sample_type, scales = "free_x", space = "free_x", labeller = labeller(sample_type = supp.labs)) # I don't know why D0-5 cells make it past SCN filter, but we just kick them out

rm(supp.labs)

#####Fig 13
#Here, we plot the entropy score for the direct reprogramming dataset; we plot along with the mouse in vivo datasets to enable direct comparison. Because CM-like cells do not emerge until D3, we plot only cells D3 and after. However, we perform a second analysis here as well. In the original manuscript (Stone and Gifford, Cell Stem Cell 2019), the authors classified cells into different groups based on their point in a reprogramming trajectory (Fig 1F). We used their labeling approach to classify cells as "Ccnb1 Alternative" (clusters 4, 10 in their analysis), "Mmp3 Alternative" (clusters 1, 3, 7 in their analysis), and "Tnni3 Main" (clusters 5, 0, 2 in their analysis). We then plot entropy scores for all of the cells (including the non-CMs) based on their classification into those three groups.

supp.labs = c("In Vivo", "Reprogramming")
names(supp.labs) = c("in vivo", "direct reprogramming")

###Fig 13a
ggplot(combined_datasets[combined_datasets$species == "mouse" & combined_datasets$good_cell == TRUE & combined_datasets$genes > 1000 & combined_datasets$include_dataset == TRUE & !combined_datasets$timepoint %in% c("DM1", "D1", "D2"), ], aes(x = timepoint, y = entropy, fill = study)) + geom_violin(scale = "width", position = position_dodge(width = 1), lwd =  1/fig_factor) + geom_boxplot(position = position_dodge(width = 1), lwd = 1/fig_factor, outlier.size = 1.5/fig_factor) + xlab("Timepoint") + ylab("Entropy Score") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 16/fig_factor), axis.text.x = element_text(angle = 45), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor), strip.text = element_text(size = 16/fig_factor)) + scale_fill_discrete(name = "Study") + guides(fill = guide_legend(ncol = 1)) + facet_grid(~sample_type, scales = "free_x", space = "free_x", labeller = labeller(sample_type = supp.labs))

#I initally misinterpretted the authors' claims when I was labeling the metadata. I thought that the Ccnb1 population represented the start of the trajectory; the authors actually claim it is a third branch from the trajectory (which starts in the middle of all three trajectories). To plot appropriately, we rename the groups here.
reprogramming = combined_datasets[combined_datasets$full_label == "mouse direct reprogramming" & combined_datasets$top5_norm < 1.3 & combined_datasets$depth_norm > -0.5 & combined_datasets$include_dataset == TRUE, ]
reprogramming$label = ""
reprogramming$other_meta = sapply(strsplit(reprogramming$other_meta, " "), "[", 1)
reprogramming[reprogramming$other_meta == "start", ]$label = "Ccnb1 Alternative"
reprogramming[reprogramming$other_meta == "refractory", ]$label = "Mmp3 Alternative"
reprogramming[reprogramming$other_meta == "main", ]$label = "Tnni3 Main"
reprogramming$label = factor(reprogramming$label, levels = c("Ccnb1 Alternative", "Mmp3 Alternative", "Tnni3 Main"))

###Fig 13b
ggplot(reprogramming, aes(x = timepoint, y = entropy, fill = label)) + geom_violin(scale = "width", position = position_dodge(width = 0.9), lwd =  1/fig_factor) + geom_boxplot(position = position_dodge(width = 0.9), lwd = 1/fig_factor, outlier.size = 1.5/fig_factor) + xlab("Timepoint") + ylab("Entropy Score") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 16/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor), strip.text = element_text(size = 16/fig_factor), legend.key.size = unit(1.2/fig_factor, "lines")) + scale_fill_discrete(name = "Group") + guides(fill = guide_legend(ncol = 1))

rm(reprogramming, supp.labs)

#####Fig 14
#This figure shows the distribution of gene abundances per timepoint; we compute, for each timepoint in the reference data, the CPM of each gene and then make density plots for the gene abundances. This code is heavily adapted from Michael Farid's code.

kannan_ref_pheno = combined_datasets[combined_datasets$data == "kannan_ref_data", ]
kannan_plot_data = sweep(kannan_ref_data, 2, colSums(kannan_ref_data), "/")
#average over all cells in a timepoint and wrangle the data into a numeric data.frame
aggByTimepoint=t(aggregate(t(as.matrix(kannan_plot_data)), by=list(kannan_ref_pheno$timepoint), mean))
aggByTimepoint=aggByTimepoint[2:nrow(aggByTimepoint),]
aggByTimepoint=matrix(as.numeric(aggByTimepoint),nrow(aggByTimepoint),ncol(aggByTimepoint)) 
aggByTimepoint=data.frame(aggByTimepoint)
colnames(aggByTimepoint)=c("e14", "e18","p0","p4","p8","p11","p14", "p18", "p22", "p28", "p35", "p56")
row.names(aggByTimepoint)=row.names(kannan_ref_data)
maybe=melt(aggByTimepoint)

#Please note - this figure has a heavily adjusted bandwidth as well as a hard tail cutoff on the x axis (the maximum expressed gene has CPM ~72000). This needs to be taken into context when interpretting this figure. We do these changes to make the gene trends visualizable. However, there are obviously genes expressed above 5000 CPM; the y-axis is also kind of uninterpettable.

ggplot(maybe, aes(value * 10^6, colour = variable)) +
  geom_density(adjust=1000, size=1/fig_factor)+scale_x_continuous(limits = c(0,5000))+
  xlab("CPM")+ylab("Adjusted Density")+labs(color="Timepoint") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 16/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor), legend.key.size = unit(1.2/fig_factor, "lines"))

rm(aggByTimepoint, maybe, kannan_plot_data, kannan_ref_pheno)

#####Fig15
#This figure shows the ratio of entropy computed on 3' counts vs UMIs (i.e. the data pre- and post-UMI collapsing) for four of our datasets (Kannan Reference, Kannan Live/Fixed, Kannan Chamber, and Murphy), where we had both sets of data as output from zUMIs.

umi_ratio1 = combined_datasets[combined_datasets$data == "kannan_ref_data" & combined_datasets$good_cell == TRUE & combined_datasets$genes > 1000, ]
umi_ratio2 = combined_datasets[combined_datasets$data == "kannan_fixed_data" & combined_datasets$good_cell == TRUE & combined_datasets$genes > 1000, ]
umi_ratio3 = combined_datasets[combined_datasets$data == "kannan_chamber_data" & combined_datasets$good_cell == TRUE & combined_datasets$genes > 1000, ]
umi_ratio4 = combined_datasets[combined_datasets$data == "murphy_data" & combined_datasets$good_cell == TRUE & combined_datasets$genes > 1000, ]

umi_ratio1$ratio = (umi_ratio1$entropy)/(combined_datasets[combined_datasets$data == "kannan_ref_reads_data" & combined_datasets$cellname %in% umi_ratio1$cellname, ]$entropy)
umi_ratio2$ratio = (umi_ratio2$entropy)/(combined_datasets[combined_datasets$data == "kannan_fixed_reads_data" & combined_datasets$cellname %in% umi_ratio2$cellname, ]$entropy)
umi_ratio3$ratio = (umi_ratio3$entropy)/(combined_datasets[combined_datasets$data == "kannan_chamber_reads_data" & combined_datasets$cellname %in% umi_ratio3$cellname, ]$entropy)
umi_ratio4$ratio = (umi_ratio4$entropy)/(combined_datasets[combined_datasets$data == "murphy_reads_data" & combined_datasets$cellname %in% umi_ratio4$cellname, ]$entropy)

umi_ratio = rbind(umi_ratio1, umi_ratio2, umi_ratio3, umi_ratio4)

#Please note - I plotted the inverse of the ratio calculated above. There was no particular reason to choose either direction - I went with this one (3' Counts:UMIs) because it was easier to write in the main body of the text that the reverse. Results should be comparable either way, obviously.

ggplot(umi_ratio, aes(x = timepoint, y = 1/ratio, fill = study)) + geom_violin(scale = "width", position = position_dodge(width = 1), lwd =  1/fig_factor) + geom_boxplot(position = position_dodge(width = 1), lwd = 1/fig_factor, outlier.size = 1.5/fig_factor) + xlab("Timepoint") + ylab("3' Counts:UMI Entropy Score Ratio") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 16/fig_factor), axis.text.x = element_text(angle = 45), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor), legend.key.size = unit(1.2/fig_factor, "lines"), strip.text = element_text(size = 16/fig_factor)) + scale_fill_discrete(name = "Study") + guides(fill = guide_legend(ncol = 1))
median(umi_ratio$ratio)

rm(umi_ratio, umi_ratio1, umi_ratio2, umi_ratio3, umi_ratio4)

#####Fig 16
#This figure shows integration of five perinatal datasets using functionality in Seurat V3. The code is adapted from the following Seurat vigneete: https://satijalab.org/seurat/v3.1/integration.html. The purpose of this figure was to demonstrate that, even following integration, CMs aren't always well-aggregated by timepoint, suggesting persistent batch effects. Despite this, entropy score is still able to capture the appropriate maturation status. We used these five datasets: 10x Chromium, Duan 10x, Kannan Reference, Murphy, and Tabula Muris Senis (because they all had overlap around the perinatal timepoints). We used the SCTransform method, but set the Kannan reference as the reference dataset. The rationale was to mimic a sort of alternative approach to predicting maturation status, whereby test datasets are projected onto a reference (such as our maturation reference) and then assigned some metric.

makeSeurat = function(data){
  good_cells = combined_datasets[combined_datasets$data == data & combined_datasets$good_cell == TRUE, ]$cellname
  tab = mito_correct(get(data)[, good_cells])
  meta.data = combined_datasets[combined_datasets$data == data & combined_datasets$good_cell == TRUE, ]
  colnames(tab) = paste(good_cells, data, sep = "_")
  rownames(meta.data) = paste(good_cells, data, sep = "_")
  return(CreateSeuratObject(tab, meta.data = meta.data))
}

kannan_ref_seurat = makeSeurat("kannan_ref_data")
murphy_seurat = makeSeurat("murphy_data")
duan_seurat = makeSeurat("duan10_data")
chromium_seurat = makeSeurat("chromium_data")
tabula_seurat = makeSeurat("tabula_data")

#To run this code, you may need quite a bit of RAM (see the changed settings to future).
options(future.globals.maxSize = 1500 * 1024^2)
heart.list = list(kannan_ref_seurat, murphy_seurat, duan_seurat, chromium_seurat, tabula_seurat)
for (i in 1:length(heart.list)) {
  heart.list[[i]] <- SCTransform(heart.list[[i]], verbose = FALSE)
}
heart.features <- SelectIntegrationFeatures(object.list = heart.list, nfeatures = 3000)
heart.list <- PrepSCTIntegration(object.list = heart.list, anchor.features = heart.features, 
                                    verbose = FALSE)
heart.anchors <- FindIntegrationAnchors(object.list = heart.list, normalization.method = "SCT", 
                                           anchor.features = heart.features, verbose = FALSE, k.filter = 100, dims = 1:10, reference = 1)
heart.integrated <- IntegrateData(anchorset = heart.anchors, normalization.method = "SCT", 
                                     verbose = FALSE, dims = 1:10)
DefaultAssay(heart.integrated) <- "integrated"
heart.integrated <- RunPCA(heart.integrated, npcs = 10, verbose = FALSE)
heart.integrated <- RunUMAP(heart.integrated, reduction = "pca", dims = 1:10)
heart.integrated$timepoint = factor(heart.integrated$timepoint, levels = lev3)

###Fig 16a
DimPlot(heart.integrated, reduction = "umap", group.by = "timepoint", pt.size = 2/fig_factor) + xlab("UMAP 1")+ylab("UMAP 2")+labs(color="Timepoint") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 16/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor), legend.key.size = unit(1.2/fig_factor, "lines"))

###Fig 16b
DimPlot(heart.integrated, reduction = "umap", group.by = "study", pt.size = 2/fig_factor) + xlab("UMAP 1")+ylab("UMAP 2")+labs(color="Study") + theme_linedraw() + theme(axis.title = element_text(size = 24/fig_factor), axis.text = element_text(size = 16/fig_factor), legend.title = element_text(size = 24/fig_factor), legend.text = element_text(size = 12/fig_factor), legend.key.size = unit(1.2/fig_factor, "lines"))

###Fig16c
#Note - I was too lazy to format this figure (would have to handle each ggplot individually and it wasn't worth it)
FeaturePlot(heart.integrated, features = c("Tnni3", "mt-Co1", "Ckmt2", "Cdk1"))

rm(heart.list, heart.features, heart.anchors, heart.integrated, kannan_ref_seurat, duan10_data, chromium_seurat, murphy_seurat, tabula_seurat, makeSeurat)
