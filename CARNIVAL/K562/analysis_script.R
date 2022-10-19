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

dir.create("output")

# interactions <- import_omnipath_interactions(resources = "BioGRID")
# interactions <- interactions[which(interactions$is_directed==1), ]
# interactions <- interactions[which((interactions$is_stimulation+interactions$is_inhibition)==1), ]
# ppi <- matrix(data = "1", nrow = nrow(interactions), ncol = 3)
# ppi[, 1] <- interactions$source_genesymbol
# ppi[which(interactions$is_inhibition==1), 2] <- "-1"
# ppi[, 3] <- interactions$target_genesymbol
# colnames(ppi) <- c("source", "sign", "target")
# ppi <- as.data.frame(ppi)
load(file = "../ppi.RData")

ff1 <- list.files("../../Gene_Expression/output/k562/")[
  which(grepl(pattern = "ttop_", x = list.files("../../Gene_Expression/output/k562/"), fixed = TRUE))]
ff2 <- list.files("../../Gene_Expression/output/k562/")[
  which(grepl(pattern = "tf_act_", x = list.files("../../Gene_Expression/output/k562/"), fixed = TRUE))]

cases <- intersect(x = sapply(strsplit(x = ff1, split = "_", fixed = TRUE), "[", 2), 
                   y = sapply(strsplit(x = ff2, split = "_", fixed = TRUE), "[", 3))

for(ii in 1:length(cases)){
  
  load(file = paste0("../../Gene_Expression/output/k562/ttop_", cases[ii], "_k562.RData"))
  bn <- ppi[intersect(x = which(ppi$source%in%ttop$ID), 
                      y = which(ppi$target%in%ttop$ID)), ]
  
  load(file = paste0("../../Gene_Expression/output/k562/tf_act_", cases[ii], "_k562.RData"))
  topTF <- input.scores[which(input.scores$pval<=0.05), ]
  measObj <- matrix(data = , nrow = 1, ncol = nrow(topTF))
  measObj[1, ] <- topTF$nes
  colnames(measObj) <- topTF$id
  measObj <- as.data.frame(measObj)
  
  res <- runCARNIVAL(inputObj = NULL,
                     measObj = measObj,
                     netObj = bn,
                     weightObj = NULL,
                     solverPath = "/beegfs/homes/egjerga/cplex",
                     solver = "cplex",
                     timelimit = 10800,
                     limitPop = 100,
                     poolCap = 100,
                     mipGAP = 0.1,
                     poolrelGAP = 0.1,
                     poolIntensity = 4, 
                     threads = 20)
  
  save(res, file = paste0("output/res_", cases[ii], ".RData"))
  
}

ff <- list.files("output/")
for(ii in 1:length(ff)){

  load(file = paste0("output/", ff[ii]))
  case <- gsub(pattern = "res_", replacement = "", fixed = TRUE,
               x = gsub(pattern = ".RData", replacement = "", x = ff[ii], fixed = TRUE))

  write.table(x = res$weightedSIF, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE,
              file = paste0("../encore_network_counts/k562/network_", case, ".txt"))

  write.table(x = res$nodesAttributes, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE,
              file = paste0("../encore_network_counts/k562/attributes_", case, ".txt"))

}
