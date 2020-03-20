#Transcriptomic entropy enables quantification of cardiomyocyte maturation at single cell level
#Chulan Kwon Laboratory
#Primary author: Suraj Kannan
#February 13, 2020
#Document: Functions

#The following document contains the codebase for functions used in our transcriptomic entropy project. It is intended to be paired with the corresponding R workspace, which contains several objects necessary for the successful running of these scripts. The functions are automatically loaded as part of the workspace as well.

#Functions to include: data_qc
#get_gene_entropy, get_gene_entropy_list, plot_gene_entropy, plot_gene_entropy_list, plot_gene_entropy_list_groups

library(ggplot2)
library(pheatmap)
library(singleCellNet)
library(Matrix)

####GENERAL DATA HANDLING FUNCTIONS####

#mito_correct 
#(requires mito_map)
#This function is used to correct datasets that have been mapped to the genome or used a full reference including pseudogenes. We observed that many mitochondrial genes can be mismapped by these methods to pseudogenes, which may in turn affect the calculated entropy (given the high abundance of mitochondrial genes in CMs). We developed a correction method that adds pseudogenes to the corresponding mitochondrial gene. Note that currently, this function only works for mouse datasets. This is because we did not observe significant pseudogene contamination in the human datasets not handled with CellRanger; thus, we felt a human map at this time was unnecessary.

mito_correct = function(data){
  data = as.data.frame(as.matrix(data))
  rownames(data) = str_replace(rownames(data), "[.]", "-") #This stupid correction came about when I found out that sometimes, converting to dataframe replaces hyphens with periods for no reason. y u do dis r
  for (mito in names(mito_map)){
    mito_bin = unlist(mito_map[mito])
    if(any(mito_bin %in% rownames(data))){
      mito_bin = mito_bin[mito_bin %in% rownames(data)]
      data[mito, ] = colSums(data[mito_bin, ])
      data = data[rownames(data) == mito | !rownames(data) %in% mito_bin, ]
    }
  }
  return(Matrix(as.matrix(data, sparse = TRUE)))
}

#mito 
#(requires G_list, G_list_human)
#This function calculates the percentage of mitochondrial reads in the dataset. By default, we only include genes that would be used in the entropy calculation - protein coding, antisense, and lncRNAs. This will notably exclude the mitochondrial tRNAs and, more significantly, the mitochondrial rRNAs. If the "all" flag is set to TRUE, this function will compute the full percentage of mitochondrial reads. By default, this function will also perform a mitochondrial correct; this can be turned off with the "correct" flag.

mito = function(dataset, species = "mouse", all = FALSE, correct = TRUE){
  if(ncol(dataset) > 5000){
    lst = split(colnames(dataset), (seq(ncol(dataset))-1) %/% 5000)
    result = unlist(lapply(lst, function(x){mito(dataset[, x], species = species, all = all, correct = correct)}))
  }
  else{
    mitostring = "mt-"
    if(species == "mouse" & correct == TRUE){
      dataset = mito_correct(dataset)
    }
    if(species == "human"){
      mitostring = "MT-"
      G_list = G_list_human
    }
    if (all == FALSE){
      dataset = dataset[rownames(dataset) %in% G_list[G_list$use_for_entropy == TRUE, ]$symbol, ]
    }
    result = (colSums(dataset[startsWith(rownames(dataset), mitostring), ])/colSums(dataset))
  }
  return(result)
}

#ribo 
#(requires G_list, G_list_human)
#This function calculates the percentage of ribosomal protein-coding reads in the dataset. By default, we only include genes that would be used in the entropy calculation - protein coding, antisense, and lncRNAs. If the "all" flag is set to TRUE, this function will compute the full percentage of ribosomal reads.

ribo = function(dataset, species = "mouse", all = FALSE){
  if(ncol(dataset) > 5000){
    lst = split(colnames(dataset), (seq(ncol(dataset))-1) %/% 5000)
    result = unlist(lapply(lst, function(x){ribo(dataset[, x], species = species, all = all)}))
  }
  else{
    ribostring1 = "Rps"
    ribostring2 = "Rpl"
    if(species == "human"){
      ribostring1 = "RPS"
      ribostring2 = "RPL"
      G_list = G_list_human
    }
    if (all == FALSE){
      dataset = dataset[rownames(dataset) %in% G_list[G_list$use_for_entropy == TRUE, ]$symbol, ]
    }
    result = (colSums(dataset[startsWith(rownames(dataset), ribostring1) | startsWith(rownames(dataset), ribostring2), ])/colSums(dataset))
  }
  return(result)
}

#rename_genes 
#(requires G_list, G_list_human)
#This function renames genes from ENSEMBL ids to gene symbol (mgi for mouse, hgcn for human). Our approach is fairly blunt in that we filter out ENSEMBL ids with no known names, as well as duplicated gene symbols. However, this usually suffices to capture the vast majority of interesting data. Please note that this function assumes no version numbers appended to the id.

rename_genes = function(data, species = "mouse") {
  if(species == "human"){
    G_list = G_list_human
  }
  data = as.data.frame(as.matrix(data))
  data$ensembl_gene_id = rownames(data)
  data = merge(data, G_list, by = "ensembl_gene_id", all.x = TRUE)
  data$ensembl_gene_id = NULL
  data$gene_biotype = NULL
  data$use_for_entropy = NULL
  data = data[!duplicated(data$symbol) & !is.na(data$symbol) & data$symbol != "", ]
  rownames(data) = data$symbol
  data$symbol = NULL
  data = Matrix(as.matrix(data), sparse = TRUE)
  return(data)
}

#scn and scn_human 
#(requires expT, expTrain, oTab, cgenesA, cgenesA_human, grps, xpairs, and xpairs_human; also needs singleCellNet loaded)
#These functions compute cell type classifications of the input dataset using SingleCellNet (Tan et al., Cell Systems 2019) against a Tabula Muris reference. The function takes as input data (the data table), pheno (a phenotype table containing, at minimum, a column for timepoint), study (a name for the study), and combined_datasets (named based on a legacy version of this function - it is basically any dataframe that will house the output of this function). For the sake of simplicity, a separate mouse and human version of this function are provided.

scn = function(data, pheno, study, combined_datasets){
  cgenesB = cgenesA[cgenesA %in% rownames(data)]
  xpairs_list = as.data.frame(matrix(unlist(strsplit(xpairs, "_")), nrow=length(xpairs), byrow=T))
  xpairsB = xpairs[xpairs_list$V1 %in% cgenesB & xpairs_list$V2 %in% cgenesB]
  system.time(pdTrain<-query_transform(expT[cgenesB, ], xpairsB))
  system.time(rf_tspAll<-sc_makeClassifier(pdTrain[xpairsB,], genes=xpairsB, groups=grps, nRand=10, ntrees=500))
  print("Done making classifier")
  system.time(kidTransAll<-query_transform(data[cgenesB,], xpairsB))
  print("A")
  system.time(crParkall<-rf_classPredict(rf_tspAll, kidTransAll, numRand=10))
  print("B")
  sgrp<-as.vector(pheno$timepoint)
  names(sgrp)<-rownames(pheno)
  sc_hmClass(crParkall, sgrp, max=5000, isBig=TRUE, cCol=F, font=8)
  print("C")
  crParkall = crParkall[, names(sgrp)]
  combined_datasets[combined_datasets$study == study, ]$cm_score = crParkall["cardiac muscle cell", ]
  combined_datasets[combined_datasets$study == study, ]$max_score = apply(crParkall, 2, max)
  combined_datasets[combined_datasets$study == study, ]$max_celltype = rownames(crParkall)[apply(crParkall,2,which.max)]
  return(combined_datasets)
}

scn_human = function(data, pheno, study, combined_datasets){
  cgenesB = cgenesA_human[cgenesA_human %in% rownames(data)]
  xpairs_list = as.data.frame(matrix(unlist(strsplit(xpairs_human, "_")), nrow=length(xpairs_human), byrow=T))
  xpairsB = xpairs_human[xpairs_list$V1 %in% cgenesB & xpairs_list$V2 %in% cgenesB]
  system.time(pdTrain<-query_transform(expTrain[cgenesB, ], xpairsB))
  system.time(rf_tspAll<-sc_makeClassifier(pdTrain[xpairsB,], genes=xpairsB, groups=grps, nRand=10, ntrees=500))
  print("Done making classifier")
  system.time(kidTransAll<-query_transform(data[cgenesB,], xpairsB))
  print("A")
  system.time(crParkall<-rf_classPredict(rf_tspAll, kidTransAll, numRand=10))
  print("B")
  sgrp<-as.vector(pheno$timepoint)
  names(sgrp)<-rownames(pheno)
  sc_hmClass(crParkall, sgrp, max=5000, isBig=TRUE, cCol=F, font=8)
  print("C")
  crParkall = crParkall[, names(sgrp)]
  combined_datasets[combined_datasets$study == study, ]$cm_score = crParkall["cardiac muscle cell", ]
  combined_datasets[combined_datasets$study == study, ]$max_score = apply(crParkall, 2, max)
  combined_datasets[combined_datasets$study == study, ]$max_celltype = rownames(crParkall)[apply(crParkall,2,which.max)]
  return(combined_datasets)
}

#top_n_genes
#This function returns the percentage of counts that go the top n genes, as specified by the user. As of now, this function computes based on *all* counts, not just those in the categories used for entropy calculation. This function is used for QC filtering.

top_n_genes = function(expressionMatrix, n){
  if(ncol(expressionMatrix) > 10000){
    lst = split(colnames(expressionMatrix), (seq(ncol(expressionMatrix))-1) %/% 10000)
    return(unlist(lapply(lst, function(x){top_n_genes(expressionMatrix[, x], n)})))
  }
  else{
    p=sweep(expressionMatrix,MARGIN=2,FUN="/",STATS=colSums(expressionMatrix));
    return(apply(p, 2, function(x) sum(tail(sort(x), n))))
  }
}

#top_n_table
#This function returns a matrix of the names of the top n genes, as specified by the user, expressed in each cell. This function has been immensely invaluable for troubleshooting and exploring datasets quickly.

top_n_table = function(data, n) {
  readout = matrix(nrow = n, ncol = ncol(data))
  for (i in 1:ncol(data)) {
    readout[, i] = tail(names(sort(structure(data[, i], names = rownames(data), class = "numeric"))), n)
  }
  colnames(readout) = colnames(data)
  return(readout)
}

#I wrote a top_n_table2 that returns only genes used in the entropy calculation
top_n_table2 = function(data, n, species = "mouse") {
  G_list = G_list
  if(species == "human"){
    G_list = G_list_human
  }
  data = data[rownames(data) %in% G_list[G_list$use_for_entropy == TRUE, ]$symbol & !startsWith(rownames(data), "RPS") & !startsWith(rownames(data), "RPL") & !startsWith(rownames(data), "Rps") & !startsWith(rownames(data), "Rpl"), ]
  readout = matrix(nrow = n, ncol = ncol(data))
  for (i in 1:ncol(data)) {
    readout[, i] = tail(names(sort(structure(data[, i], names = rownames(data), class = "numeric"))), n)
  }
  colnames(readout) = colnames(data)
  return(readout)
}

#master_entropy
#(requires G_list and G_list_human)
#This function is the primary entropy-calculating function used in this project (unless otherwise noted). The name reflects that this function went through many changes and iterations; the "master" reflecting that this was supposed to be finalized version. (Ironically, the master_entropy function has itself gone through many changes, but given how habituated we have become with this name, it continues to be the working name of this function). This function automatically performs mitochondrial correction using mito_correct and subselects protein coding, antisense, and lincRNA genes. Additionally, genes coding for ribosomal proteins are removed. Entropy is calculated up to n genes, as specified by the user; the default is 1000.

master_entropy = function(expressionMatrix, n = 1000, species = "mouse"){
  #This portion splits up large datasets to go easier on memory usage.
  if(ncol(expressionMatrix) > 5000){
    lst = split(colnames(expressionMatrix), (seq(ncol(expressionMatrix))-1) %/% 5000)
    result = unlist(lapply(lst, function(x){master_entropy(expressionMatrix[, x], n = n, species = species)}))
  }
  else{
    if(species == "mouse"){
      expressionMatrix = mito_correct(expressionMatrix) 
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

#data_qc
#(requires just about every previous function to run properly)
#This function wraps together many of the other functions provided above to compute a consolidated QC profile for the dataset. The function takes the following main inputs: dataset (provided here as a character, rather than as the dataset itself!), study (the name of the study), timepoint_list (a vector containing all of the timepoints), scn_calc (a logical on whether the scn classification function should be run; could be skipped, for example, if the authors want to save some time but it is advised not to skip), species, sample type (e.g. in vivo, directed differentiation), isolation, sequencing, mapping, datatype, doi (all characteristics of the study), other_meta (a vector containing any other information that the authors want to include), append_to_alldata (in case the study has not been yet documented in the alldata object). A number of outputs are provided, including important QC metrics and whether the cell is usable, mitochodrial and ribosomal percentage readouts, and the entropy.

data_qc = function(dataset, study, timepoint_list = NA, scn_calc = TRUE, species = "mouse", sample_type = NA, append_to_alldata = FALSE, append_only = FALSE, isolation = NA, sequencing = NA, mapping = NA, datatype = NA, doi = NA, other_meta = NA){
  data = get(dataset)
  qc = as.data.frame(matrix(nrow = ncol(data), ncol = 27))
  colnames(qc) = c("timepoint", "study", "depth", "genes", "entropy", "top5", "top5_norm", "depth_norm", "mito", "mito_all", "ribo", "cm_score", "max_score", "max_celltype", "good_cell", "data", "cellname", "other_meta", "include_dataset", "reason", "isolation", "sequencing", "mapping", "datatype", "species", "sample_type", "full_label")
  rownames(qc) = colnames(data)
  qc$study = study
  qc$timepoint = timepoint_list
  qc$depth = colSums(data)
  qc$genes = colSums(data > 0)
  if(append_to_alldata == TRUE){
    alldata[study, ] <<- c(isolation, sequencing, mapping, datatype, as.numeric(median(qc$depth)), as.numeric(median(qc$genes)), as.numeric(nrow(qc)), paste(unique(qc$timepoint), collapse = ", "), doi, dataset, species, sample_type)
  }
  if(append_only == FALSE){
    qc$entropy = master_entropy(data, species = species)
    
    qc$top5 = top_n_genes(data, 5)
    qc$top5_norm = 0
    if(!is.na(unique(qc$timepoint))){
      for (j in unique(qc$timepoint)){
        qc[qc$timepoint == j, ]$top5_norm = qc[qc$timepoint == j, ]$top5/median(qc[qc$timepoint == j, ]$top5)
      }
    }
    
    qc$depth_norm = 0
    if (!is.na(unique(qc$timepoint))){
      for (j in unique(qc$timepoint)){
        qc[qc$timepoint == j, ]$depth_norm = log(qc[qc$timepoint == j, ]$depth/median(qc[qc$timepoint == j, ]$depth))
      }
    }
    
    qc$mito = mito(data, species = species)
    qc$mito_all = mito(data, species = species, all = TRUE)
    qc$ribo = ribo(data, species = species)
    
    if(scn_calc == TRUE){
      if(species == "human"){
        scn_data = scn_human(data, qc, study, qc)
      }
      else {
        scn_data = scn(data, qc, study, qc)
      }
      qc$cm_score = scn_data$cm_score
      qc$max_score = scn_data$max_score
      qc$max_celltype = scn_data$max_celltype
    }
    
    #Please note that the gene number isn't hardcoded into the computation of "good cell" - this is because other uses may want to compute other entropy scores with different gene numbers but still relevantly compute a good cell
    qc$good_cell = ((qc$top5_norm < 1.3) & (qc$depth_norm > -0.5) & (qc$max_celltype == "cardiac muscle cell"))
    qc$data = dataset
    qc$cellname = colnames(data)
    qc$include_dataset = TRUE
    qc$isolation = isolation
    qc$sequencing = sequencing
    qc$mapping = mapping
    qc$datatype = datatype
    qc$species = species
    qc$sample_type = sample_type
    qc$full_label = paste(qc$species, qc$sample_type)
    qc$other_meta = other_meta
    return(qc)
  }
}
