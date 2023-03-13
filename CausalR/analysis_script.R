#! /usr/bin/env Rscript

set.seed(123)

library("readxl")
library("readr")
library("XML")
library("igraph")
library("foreach")
library("doParallel")
library("ranger")
library("palmerpenguins")
library("tidyverse")
library("kableExtra")
library("CARNIVAL")
library("OmnipathR")
library("CausalR")

dir.create("sif")
dir.create("output")

load(file = "../CARNIVAL/ppi_enio.RData")
nodes <- unique(c(ppi$source, ppi$target))
ppi$sign[which(ppi$sign=="1")] <- "Activates"
ppi$sign[which(ppi$sign=="-1")] <- "Inhibits"
write.table(x = ppi, file = "sif/net.sif", quote = FALSE, sep = "\t", 
            row.names = FALSE, col.names = FALSE)

ccg <- CreateCCG(filename = "sif/net.sif")

load(file = "../Perturb_Seq_Benchmarking/gene_list.RData")
cases <- names(gene_list)[which(grepl(pattern = "-", x = names(gene_list), fixed = TRUE))]
cases <- sapply(strsplit(x = cases, split = " - ", fixed = TRUE), "[", 1)
ff <- paste0("ttop_", tolower(cases), "_k562.RData")

## K562

dir.create("output/k562")

registerDoParallel(cores=15)

foreach(ii = 1:length(ff)) %dopar% {
  
  print(paste0("Step ---- ", ii, "/", length(ff)))
  
  load(file = paste0("../Gene_Expression/output/k562/", ff[ii]))
  ttop <- ttop[which(ttop$ID%in%nodes), ]
  
  df <- matrix(data = , nrow = nrow(ttop), ncol = 2)
  df[, 1] <- ttop$ID
  df[, 2] <- 0
  
  idx1 <- intersect(x = which(ttop$FDR<=0.05), y = which(ttop$logFC>0))
  if(length(idx1)>0){
    df[idx1, 2] <- 1
  }
  
  idx2 <- intersect(x = which(ttop$FDR<=0.05), y = which(ttop$logFC<0))
  if(length(idx2)>0){
    df[idx2, 2] <- -1
  }
  
  scanResults <- runSCANR(ccg, df, numberOfDeltaToScan=5,
                          topNumGenes=100,correctPredictionsThreshold=1,
                          writeResultFiles = FALSE, writeNetworkFiles = "none",quiet=FALSE)
  
  save(scanResults, file = paste0("output/k562/scan_results_", 
                                  sapply(strsplit(x = ff[ii], split = "_", fixed = TRUE), "[", 2), 
                                  ".RData"))
  
}

resCausalR <- list()
for(ii in 1:length(cases)){
  
  load(file = paste0("output/k562/scan_results_", tolower(cases[ii]), ".RData"))
  nodes <- scanResults$uniqueNodes$names
  nn <- c()
  for(jj in 1:length(nodes)){
    nn <- c(nn, substr(x = nodes[jj], start = 1, stop = nchar(nodes[jj])-1))
  }
  
  resCausalR[[length(resCausalR)+1]] <- nn
}
names(resCausalR) <- cases

save(resCausalR, file = "output/k562/resCausalR.RData")
