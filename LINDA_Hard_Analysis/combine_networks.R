library("readr")
library("XML")
library("igraph")

dir.create("encore_network_counts/")

source("prepare_cytoscape_visualization.R")

load(file = system.file("extdata", "digger_human_transcripts.RData", package = "LINDA"))

## HepG2
dir.create("encore_network_counts/hepg2_networks")

load(file = "HepG2/ctrl/output/res_hepg2_ctrl.RData")
load(file = "HepG2/kd/output/res_hepg2_kd.RData")
all(names(res_hepg2_ctrl)==names(res_hepg2_kd))

cases <- intersect(x = names(res_hepg2_ctrl), y = names(res_hepg2_kd))

for(ii in 1:length(cases)){
  
  print(paste0("Step ---- ", ii, "/", length(cases)))
  
  load(file = paste0("../Gene_Expression/output/hepg2/ttop_", cases[ii], "_hepg2.RData"))
  load(file = paste0("../Gene_Expression/output/hepg2/tf_act_", cases[ii], "_hepg2.RData"))
  load(file = paste0("../Transcript_Expression/output/hepg2/es_data_", cases[ii], "_hepg2.RData"))
  
  background.network <- bg[intersect(x = which(bg$gene_source%in%ttop$ID), 
                                     y = which(bg$gene_target%in%ttop$ID)), ]
  background.network <- background.network[complete.cases(background.network), ]
  
  tf_act <- input.scores
  tf_act$nes <- 0
  tf_act$nes[which(tf_act$pval<=0.05)] <- 1
  
  res_ctrl <- res_hepg2_ctrl[[ii]]
  res_kd <- res_hepg2_kd[[ii]]
  
  net <- prepare_cytoscape_visualization(netA = res_ctrl$combined_interactions, 
                                         netB = res_kd$combined_interactions, 
                                         spliceA = NULL, 
                                         spliceB = es_data, 
                                         pValThresh = 0.05, 
                                         background.network = background.network, 
                                         tf = tf_act$id[which(tf_act$nes==1)], 
                                         sources = "Perturbation")
  
  write.table(x = net$network, file = paste0("encore_network_counts/hepg2_networks/network_", cases[ii], ".txt"), 
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  write.table(x = net$attributes, file = paste0("encore_network_counts/hepg2_networks/attributes_", cases[ii], ".txt"), 
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
}

## K562
dir.create("encore_network_counts/k562_networks")

load(file = "K562/ctrl/output/res_k562_ctrl.RData")
load(file = "K562/kd/output/res_k562_kd.RData")
all(names(res_k562_ctrl)==names(res_k562_kd))

cases <- intersect(x = names(res_k562_ctrl), y = names(res_k562_kd))

for(ii in 1:length(cases)){
  
  print(paste0("Step ---- ", ii, "/", length(cases)))
  
  load(file = paste0("../Gene_Expression/output/k562/ttop_", cases[ii], "_k562.RData"))
  load(file = paste0("../Gene_Expression/output/k562/tf_act_", cases[ii], "_k562.RData"))
  load(file = paste0("../Transcript_Expression/output/k562/es_data_", cases[ii], "_k562.RData"))
  
  background.network <- bg[intersect(x = which(bg$gene_source%in%ttop$ID), 
                                     y = which(bg$gene_target%in%ttop$ID)), ]
  background.network <- background.network[complete.cases(background.network), ]
  
  tf_act <- input.scores
  tf_act$nes <- 0
  tf_act$nes[which(tf_act$pval<=0.05)] <- 1
  
  res_ctrl <- res_k562_ctrl[[ii]]
  res_kd <- res_k562_kd[[ii]]
  
  net <- prepare_cytoscape_visualization(netA = res_ctrl$combined_interactions, 
                                         netB = res_kd$combined_interactions, 
                                         spliceA = NULL, 
                                         spliceB = es_data, 
                                         pValThresh = 0.05, 
                                         background.network = background.network, 
                                         tf = tf_act$id[which(tf_act$nes==1)], 
                                         sources = "Perturbation")
  
  write.table(x = net$network, file = paste0("encore_network_counts/k562_networks/network_", cases[ii], ".txt"), 
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  write.table(x = net$attributes, file = paste0("encore_network_counts/k562_networks/attributes_", cases[ii], ".txt"), 
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
}
