#! /usr/bin/env Rscript

set.seed(1234)

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
library("aggregation")
library("LINDA")
library("fgsea")

dir.create("output")

ff1 <- list.files("../../Gene_Expression/output/hepg2/")[
  which(grepl(pattern = "ttop_", x = list.files("../../Gene_Expression/output/hepg2/"), fixed = TRUE))]
ff2 <- list.files("../../Gene_Expression/output/hepg2/")[
  which(grepl(pattern = "tf_act_", x = list.files("../../Gene_Expression/output/hepg2/"), fixed = TRUE))]
ff3 <- list.files("../../Transcript_Expression/output/hepg2/")[
  which(grepl(pattern = "es_data_", x = list.files("../../Transcript_Expression/output/hepg2/"), fixed = TRUE))]

cases <- intersect(x = intersect(x = sapply(strsplit(x = ff1, split = "_", fixed = TRUE), "[", 2), 
                                 y = sapply(strsplit(x = ff2, split = "_", fixed = TRUE), "[", 3)), 
                   y = sapply(strsplit(x = ff3, split = "_", fixed = TRUE), "[", 3))

load(file = system.file("extdata", "digger_human_transcripts.RData", package = "LINDA"))

##Set the LINDA optimization parameters
lambda1 <- 100
lambda2 <- 1
input.node <- NULL
mipgap = 0
relgap = 0
populate = 500
nSolutions = 100
intensity = 1
timelimit = 7200
process_log = FALSE
replace = 1
solverPath <- "/beegfs/homes/egjerga/cplex"
pValThresh <- 0.05
threads <- 8

registerDoParallel(cores=5)

ctrlList <- list()
kdList <- list()
foreach(ii = 1:length(cases)) %dopar% {
  
  load(file = paste0("../../Gene_Expression/output/hepg2/ttop_", cases[ii], "_hepg2.RData"))
  load(file = paste0("../../Gene_Expression/output/hepg2/tf_act_", cases[ii], "_hepg2.RData"))
  load(file = paste0("../../Transcript_Expression/output/hepg2/es_data_", cases[ii], "_hepg2.RData"))
  
  background.network <- bg[intersect(x = which(bg$gene_source%in%ttop$ID), 
                                     y = which(bg$gene_target%in%ttop$ID)), ]
  
  input.scores$nes <- 0
  input.scores$nes[which(input.scores$pval<=0.05)] <- 1
  input.scores <- input.scores[, c(1, 2)]
  top <- length(which(input.scores$nes==1))
  
  # Ctrl Case
  res_ctrl <- runLINDA(input.scores = input.scores, 
                       as.input = NULL, 
                       background.network = background.network, 
                       solverPath = solverPath, 
                       pValThresh = pValThresh, 
                       top = top, 
                       lambda1 = lambda1, 
                       lambda2 = lambda2, 
                       mipgap = mipgap, 
                       relgap = relgap, 
                       populate = populate, 
                       nSolutions = nSolutions, 
                       timelimit = timelimit, 
                       replace = replace, 
                       threads = threads, 
                       condition = ii, 
                       intensity = intensity, 
                       solver = "cplex", 
                       save_res = TRUE)
  
  ctrlList[[length(ctrlList)+1]] <- res_ctrl
  
  
  # KD Case
  res_kd <- runLINDA(input.scores = input.scores, 
                     as.input = es_data, 
                     background.network = background.network, 
                     solverPath = solverPath, 
                     pValThresh = pValThresh, 
                     top = top, 
                     lambda1 = lambda1, 
                     lambda2 = lambda2, 
                     mipgap = mipgap, 
                     relgap = relgap, 
                     populate = populate, 
                     nSolutions = nSolutions, 
                     timelimit = timelimit, 
                     replace = replace, 
                     threads = threads, 
                     condition = ii+length(cases), 
                     intensity = intensity, 
                     solver = "cplex", 
                     save_res = TRUE)
  
  kdList[[length(kdList)+1]] <- res_kd
  
}

names(ctrlList) <- cases
names(kdList) <- cases

save(ctrlList, file = "output/res_hepg2_ctrl.RData")
save(kdList, file = "output/res_hepg2_kd.RData")


