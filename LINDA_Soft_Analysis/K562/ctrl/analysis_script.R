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

ff1 <- list.files("../../../Gene_Expression/output/k562/")[
  which(grepl(pattern = "ttop_", x = list.files("../../../Gene_Expression/output/k562/"), fixed = TRUE))]
ff2 <- list.files("../../../Gene_Expression/output/k562/")[
  which(grepl(pattern = "tf_act_", x = list.files("../../../Gene_Expression/output/k562/"), fixed = TRUE))]
ff3 <- list.files("../../../Transcript_Expression/output/k562/")[
  which(grepl(pattern = "es_data_", x = list.files("../../../Transcript_Expression/output/k562/"), fixed = TRUE))]

cases <- c("hnrnpa2b1", "hnrnpc", "magoh", "hnrnpu", "ncbp2", "pabpn1", 
           "papola", "pcbp1", "polr2g", "ppil4", "prpf6", "prpf8", "puf60", 
           "sf1", "srsf1", "smndc1")

load(file = system.file("extdata", "digger_human_transcripts.RData", package = "LINDA"))

##Set the LINDA optimization parameters
lambda1 <- 100
lambda2 <- 1
input.node <- NULL
mipgap = 0.02
relgap = 0.02
populate = 500
nSolutions = 100
intensity = 1
timelimit = 5400
process_log = FALSE
replace = 1
solverPath <- "/beegfs/homes/egjerga/cplex"
pValThresh <- NULL
threads <- 8

registerDoParallel(cores=5)

ctrlList <- list()
kdList <- list()
foreach(ii = 1:length(cases)) %dopar% {
  
  if(!file.exists(paste0("res_", ii, ".RData"))){
    
    load(file = paste0("../../../Gene_Expression/output/k562/ttop_", cases[ii], "_k562.RData"))
    load(file = paste0("../../../Gene_Expression/output/k562/tf_act_", cases[ii], "_k562.RData"))
    load(file = paste0("../../../Transcript_Expression/output/k562/es_data_", cases[ii], "_k562.RData"))
    
    background.network <- bg[intersect(x = which(bg$gene_source%in%ttop$ID), 
                                       y = which(bg$gene_target%in%ttop$ID)), ]
    
    background.network <- background.network[complete.cases(background.network), ]
    
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
    
  }
  
}

names(ctrlList) <- cases

save(ctrlList, file = "output/res_k562_ctrl.RData")

nn <- c()
for(ii in 1:length(cases)){

  ff <- paste0("res_", ii, ".RData")
  if(file.exists(ff)){

    load(file = ff)
    ctrlList[[length(ctrlList)+1]] <- res
    nn <- c(nn, cases[ii])

  }

}
names(ctrlList) <- nn
res_k562_ctrl <- ctrlList

save(res_k562_ctrl, file = "output/res_k562_ctrl.RData")
