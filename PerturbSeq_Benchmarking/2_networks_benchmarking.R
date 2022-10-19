set.seed(1234)

library(readr)
library(fgsea)
library(LINDA)
library(ggplot2)
library(ggpubr)
library(ggrepel)

dir.create("output")

# load(file = "gene_list.RData")
load(file = "../../PerturbSeq/Phenotypes/gene_list.RData")

load(file = system.file("extdata", "digger_human_transcripts.RData", package = "LINDA"))
universe1 <- unique(c(bg$gene_source, bg$gene_target))

load(file = "../CARNIVAL/ppi.RData")
universe2 <- unique(c(ppi$source, ppi$target))

load(file = "../LINDA_Soft_Analysis/K562/ctrl/output/res_k562_ctrl.RData")
linda_su <- res_k562_ctrl
rm(res_k562_ctrl)

load(file = "../LINDA_Soft_Analysis/K562/kd/output/res_k562_kd.RData")
linda_sa <- res_k562_kd
rm(res_k562_kd)

setCases <- names(gene_list)[which(grepl(pattern = " - ", x = names(gene_list), fixed = TRUE))]
setCases <- unique(sapply(strsplit(x = setCases, split = " - ", fixed = TRUE), "[", 1))

ff <- list.files(path = "../CARNIVAL/K562/output/")

scores_su <- c()
scores_sa <- c()
scores_carn <- c()
cc <- c()

for(ii in 1:length(setCases)){
  
  idx1 <- which(names(linda_su)==tolower(setCases[ii]))
  
  idx2 <- which(names(linda_sa)==tolower(setCases[ii]))
  
  idx3 <- which(ff==paste0("res_", tolower(setCases[ii]), ".RData"))
  
  if(length(idx1)>0 && length(idx2)>0 && length(idx3)>0){
    
    net_su <- linda_su[[idx1]]$combined_interactions
    su <- unique(c(net_su[, 1], net_su[, 3]))
    
    net_sa <- linda_sa[[idx2]]$combined_interactions
    sa <- unique(c(net_sa[, 1], net_sa[, 3]))
    
    load(file = paste0("../CARNIVAL/K562/output/", ff[idx3]))
    net_carn <- res$weightedSIF
    carn <- unique(c(net_carn[, 1], net_carn[, 3]))
    
    fgseaResSU <- fora(pathways = gene_list, genes = su, universe = universe1, minSize = 1, maxSize = Inf)
    fgseaResSA <- fora(pathways = gene_list, genes = sa, universe = universe1, minSize = 1, maxSize = Inf)
    fgseaResCARN <- fora(pathways = gene_list, genes = carn, universe = universe2, minSize = 1, maxSize = Inf)
    
    px <- names(gene_list)[which(grepl(pattern = paste0(setCases[ii], " - "), x = names(gene_list), fixed = TRUE))]
    
    for(kk in 1:length(px)){
      
      ind_su <- which(fgseaResSU$pathway==px[kk])
      if(length(ind_su)==0){
        scores_su <- c(scores_su, 1)
      } else {
        scores_su <- c(scores_su, (-1)*log10(fgseaResSU$pval[which(fgseaResSU$pathway==px[kk])]))
      }
      
      ind_sa <- which(fgseaResSA$pathway==px[kk])
      if(length(ind_sa)==0){
        scores_sa <- c(scores_sa, 1)
      } else {
        scores_sa <- c(scores_sa, (-1)*log10(fgseaResSA$pval[which(fgseaResSA$pathway==px[kk])]))
      }
      
      ind_carn <- which(fgseaResCARN$pathway==px[kk])
      if(length(ind_carn)==0){
        scores_carn <- c(scores_carn, 1)
      } else {
        scores_carn <- c(scores_carn, (-1)*log10(fgseaResCARN$pval[which(fgseaResCARN$pathway==px[kk])]))
      }
      
      # cc <- c(cc, paste0(toupper(cases[jj]), " = ", px[kk]))
      cc <- c(cc, px[kk])
      
    }
  
  }
  
}

# SA vs SU
df <- matrix(data = , nrow = length(scores_sa)*2, ncol = 3)
df[, 1] <- c(scores_sa, scores_su)
df[, 2] <- c(rep("splice_aware", length(scores_sa)), rep("splice_unaware", length(scores_su)))
df[, 3] <- rep(cc, 2)
colnames(df) <- c("enrichment_score", "network", "case")
df <- as.data.frame(df)
df$enrichment_score <- as.numeric(df$enrichment_score)
df1 <- df

pdf(file = paste0("output/splice_aware_vs_splice_unaware_thresh_", thresh[ii], ".pdf"), width = 8, height = 8)
p <- ggpaired(df, x = "network", y = "enrichment_score",
              color = "network", line.color = "gray", line.size = 0.8,
              palette = "jco", label = "case", 
              repel = TRUE, label.rectangle = TRUE, font.label = 5)+
  stat_compare_means(paired = TRUE, method = "t.test")
plot(p)
dev.off()

# SA vs CARNIVAL
df <- matrix(data = , nrow = length(scores_sa)*2, ncol = 3)
df[, 1] <- c(scores_sa, scores_carn)
df[, 2] <- c(rep("splice_aware", length(scores_sa)), rep("carnival", length(scores_carn)))
df[, 3] <- rep(cc, 2)
colnames(df) <- c("enrichment_score", "network", "case")
df <- as.data.frame(df)
df$enrichment_score <- as.numeric(df$enrichment_score)
df2 <- df

# df <- df[-which(df$case=="HNRNPU - Spliceosome"), ]
# df <- df[-which(df$case=="PUF60 - Spliceosome"), ]

pdf(file = paste0("output/splice_aware_vs_carnival_thresh_", thresh[ii], ".pdf"), width = 8, height = 8)
p <- ggpaired(df, x = "network", y = "enrichment_score",
              color = "network", line.color = "gray", line.size = 0.4,
              palette = "jco", label = "case", 
              repel = TRUE, label.rectangle = TRUE)+
  stat_compare_means(paired = TRUE, method = "t.test")
plot(p)
dev.off()

tt_carn <- t.test(x = scores_sa, y = scores_carn, paired = TRUE)
tt_su <- t.test(x = scores_sa, y = scores_su, paired = TRUE)

for(ii in 1:length(thresh)){
  
  load(file = paste0("../Soft_Analysis/saRes_", thresh[ii], ".RData"))
  load(file = paste0("../Soft_Analysis/suRes_", thresh[ii], ".RData"))
  
  scores_su <- c()
  scores_sa <- c()
  scores_carn <- c()
  cc <- c()
  
  pertCases <- toupper(names(saRes))
  cases <- intersect(x = setCases, y = pertCases)
  
  for(jj in 1:length(cases)){
    
    idx1 <- which(names(suRes)==tolower(cases[jj]))
    idx2 <- which(names(saRes)==tolower(cases[jj]))
    
    if(length(idx1)>0 && length(idx2)>0){
      
      su <- suRes[[idx1]]
      su <- unique(c(su[, 1], su[, 3]))
      
      sa <- saRes[[idx2]]
      sa <- unique(c(sa[, 1], sa[, 3]))
      
      load(file = paste0("../CARNIVAL/K562/output/res_", tolower(cases[jj]), ".RData"))
      carn <- unique(c(res$weightedSIF[, 1], res$weightedSIF[, 3]))
      
      px <- names(gene_list)[which(grepl(pattern = paste0(cases[jj], " - "), x = names(gene_list), fixed = TRUE))]
      
      fgseaResSU <- fora(pathways = gene_list, genes = su, universe = universe1, minSize = 1, maxSize = Inf)
      fgseaResSA <- fora(pathways = gene_list, genes = sa, universe = universe1, minSize = 1, maxSize = Inf)
      fgseaResCARN <- fora(pathways = gene_list, genes = carn, universe = universe2, minSize = 1, maxSize = Inf)
      
      for(kk in 1:length(px)){
        
        ind_su <- which(fgseaResSU$pathway==px[kk])
        if(length(ind_su)==0){
          scores_su <- c(scores_su, 1)
        } else {
          scores_su <- c(scores_su, (-1)*log10(fgseaResSU$pval[which(fgseaResSU$pathway==px[kk])]))
        }
        
        ind_sa <- which(fgseaResSA$pathway==px[kk])
        if(length(ind_sa)==0){
          scores_sa <- c(scores_sa, 1)
        } else {
          scores_sa <- c(scores_sa, (-1)*log10(fgseaResSA$pval[which(fgseaResSA$pathway==px[kk])]))
        }
        
        ind_carn <- which(fgseaResCARN$pathway==px[kk])
        if(length(ind_carn)==0){
          scores_carn <- c(scores_carn, 1)
        } else {
          scores_carn <- c(scores_carn, (-1)*log10(fgseaResCARN$pval[which(fgseaResCARN$pathway==px[kk])]))
        }
        
        # cc <- c(cc, paste0(toupper(cases[jj]), " = ", px[kk]))
        cc <- c(cc, px[kk])
        
      }
      
    }
    
    scores_su_all <- c(scores_su_all, scores_su)
    scores_sa_all <- c(scores_sa_all, scores_sa)
    scores_carn_all <- c(scores_carn_all, scores_carn)
    cases_all <- c(cases_all, cc)
    
  }
  
  # SA vs SU
  df <- matrix(data = , nrow = length(scores_sa)*2, ncol = 3)
  df[, 1] <- c(scores_sa, scores_su)
  df[, 2] <- c(rep("splice_aware", length(scores_sa)), rep("splice_unaware", length(scores_su)))
  df[, 3] <- rep(cc, 2)
  colnames(df) <- c("enrichment_score", "network", "case")
  df <- as.data.frame(df)
  df$enrichment_score <- as.numeric(df$enrichment_score)
  df1 <- df
  
  # df <- df[-which(df$case=="HNRNPU - Spliceosome"), ]
  # df <- df[-which(df$case=="PUF60 - Spliceosome"), ]
  
  pdf(file = paste0("output/splice_aware_vs_splice_unaware_thresh_", thresh[ii], ".pdf"), width = 8, height = 8)
  p <- ggpaired(df, x = "network", y = "enrichment_score",
                color = "network", line.color = "gray", line.size = 0.8,
                palette = "jco", label = "case", 
                repel = TRUE, label.rectangle = TRUE, font.label = 5)+
    stat_compare_means(paired = TRUE, method = "t.test")
  plot(p)
  dev.off()
  
  # SA vs CARNIVAL
  df <- matrix(data = , nrow = length(scores_sa)*2, ncol = 3)
  df[, 1] <- c(scores_sa, scores_carn)
  df[, 2] <- c(rep("splice_aware", length(scores_sa)), rep("carnival", length(scores_carn)))
  df[, 3] <- rep(cc, 2)
  colnames(df) <- c("enrichment_score", "network", "case")
  df <- as.data.frame(df)
  df$enrichment_score <- as.numeric(df$enrichment_score)
  df2 <- df
  
  # df <- df[-which(df$case=="HNRNPU - Spliceosome"), ]
  # df <- df[-which(df$case=="PUF60 - Spliceosome"), ]
  
  pdf(file = paste0("output/splice_aware_vs_carnival_thresh_", thresh[ii], ".pdf"), width = 8, height = 8)
  p <- ggpaired(df, x = "network", y = "enrichment_score",
                color = "network", line.color = "gray", line.size = 0.4,
                palette = "jco", label = "case", 
                repel = TRUE, label.rectangle = TRUE)+
    stat_compare_means(paired = TRUE, method = "t.test")
  plot(p)
  dev.off()
  
  tt_carn <- t.test(x = scores_sa, y = scores_carn, paired = TRUE)
  tt_su <- t.test(x = scores_sa, y = scores_su, paired = TRUE)
  
  mm[cnt, 1] <- tt_carn$statistic
  mm[cnt, 2] <- tt_carn$p.value
  mm[cnt, 3] <- tt_su$statistic
  mm[cnt, 4] <- tt_su$p.value
  mm[cnt, 5] <- thresh[ii]
  
  cnt <- cnt + 1
  
}

## Check significance progress
colnames(mm) <- c("t_carn", "p_carn", "t_su", "p_su", "thresh")
mm <- as.data.frame(mm)

df <- matrix(data = , nrow = nrow(mm)*2, ncol = 3)
df[, 1] <- c(mm$p_carn, mm$p_su)
df[, 2] <- c(10, 20, 50, 100)
df[, 3] <- c(rep("carnival", nrow(mm)), rep("splice_unaware", nrow(mm)))
colnames(df) <- c("significance", "threshold", "case")
df <- as.data.frame(df)
df$significance <- as.numeric(df$significance)

pdf(file = "output/significance_progress.pdf", width = 8, height = 6)
p <- ggplot(df, aes(x = threshold, y = significance, fill = case)) +
  geom_line(aes(group=case)) +
  geom_point(size = 4, shape = 21)
plot(p)
dev.off()


## Compare enrichment scores for all the cases
# SA vs SU
df <- matrix(data = , nrow = length(scores_sa_all)*2, ncol = 3)
df[, 1] <- c(scores_sa_all, scores_su_all)
df[, 2] <- c(rep("splice_aware", length(scores_sa_all)), rep("splice_unaware", length(scores_su_all)))
df[, 3] <- rep(cases_all, 2)
colnames(df) <- c("enrichment_score", "network", "case")
df <- as.data.frame(df)
df$enrichment_score <- as.numeric(df$enrichment_score)

pdf(file = paste0("output/splice_aware_vs_splice_unaware_all.pdf"), width = 8, height = 8)
p <- ggpaired(df, x = "network", y = "enrichment_score",
              color = "network", line.color = "gray", line.size = 0.4,
              palette = "jco")+
  stat_compare_means(paired = TRUE, method = "t.test")
plot(p)
dev.off()

# SA vs CARNIVAL
df <- matrix(data = , nrow = length(scores_sa_all)*2, ncol = 3)
df[, 1] <- c(scores_sa_all, scores_carn_all)
df[, 2] <- c(rep("splice_aware", length(scores_sa_all)), rep("carnival", length(scores_carn_all)))
df[, 3] <- rep(cases_all, 2)
colnames(df) <- c("enrichment_score", "network", "case")
df <- as.data.frame(df)
df$enrichment_score <- as.numeric(df$enrichment_score)

pdf(file = paste0("output/splice_aware_vs_carnival_all.pdf"), width = 8, height = 8)
p <- ggpaired(df, x = "network", y = "enrichment_score",
              color = "network", line.color = "gray", line.size = 0.4,
              palette = "jco")+
  stat_compare_means(paired = TRUE, method = "t.test")
plot(p)
dev.off()


##Do a Volcano-Plot
# Splice-Aware vs CARNIVAL
df <- mm[, c(1, 2, 5)]
colnames(df) <- c("t", "pval", "case")
df$t <- as.numeric(df$t)
df$pval <- as.numeric(df$pval)
marker <- rep("stable", nrow(df))
marker[which(df$pval<=0.2)] <- "significant"

pdf(file = "output/volcano_splice_aware_vs_carnival.pdf", width = 10, height = 10)
p <- ggplot(data = df,aes(x=t,y=-log10(pval))) +
  geom_point(aes(col=marker),alpha=0.5,size=3) +
  scale_color_manual(values=c("forestgreen","black")) + 
  ggtitle(paste0("Volcano Plot PP-data: Effect Size vs -log10(pval) : n=",nrow(df)))+
  geom_label_repel(size = 5,aes(label = case),
                   # box.padding   = 0.005,
                   # point.padding = 0.01,
                   segment.color = 'grey80') +
  theme(axis.text.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x =  element_text(color = "grey20", size = 20, face = "plain") ,
        axis.title.y =  element_text(color = "grey20", size = 20, face = "plain"), 
        legend.text =  element_text(color = "grey20", size = 10, face = "plain"),
        legend.title =  element_text(color = "grey20", size = 12, face = "plain"),
        plot.title =  element_text(color = "grey20", size = 15, face = "plain"))
plot(p)
dev.off()

# Splice-Aware vs Splice-Unaware
df <- mm[, c(3, 4, 5)]
colnames(df) <- c("t", "pval", "case")
df$t <- as.numeric(df$t)
df$pval <- as.numeric(df$pval)
marker <- rep("stable", nrow(df))
marker[which(df$pval<=0.2)] <- "significant"

pdf(file = "output/volcano_splice_aware_vs_splice_unaware.pdf", width = 10, height = 10)
p <- ggplot(data = df,aes(x=t,y=-log10(pval))) +
  geom_point(aes(col=marker),alpha=0.5,size=3) +
  scale_color_manual(values=c("forestgreen","black")) + 
  ggtitle(paste0("Volcano Plot PP-data: Effect Size vs -log10(pval) : n=",nrow(df)))+
  geom_label_repel(size = 5,aes(label = case),
                   # box.padding   = 0.005,
                   # point.padding = 0.01,
                   segment.color = 'grey80') +
  theme(axis.text.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x =  element_text(color = "grey20", size = 20, face = "plain") ,
        axis.title.y =  element_text(color = "grey20", size = 20, face = "plain"), 
        legend.text =  element_text(color = "grey20", size = 10, face = "plain"),
        legend.title =  element_text(color = "grey20", size = 12, face = "plain"),
        plot.title =  element_text(color = "grey20", size = 15, face = "plain"))
plot(p)
dev.off()

## Plot together
df <- unique(rbind(df1, df2))
map <- matrix(data = , nrow = length(unique(df$case)), ncol = 2)
map[, 1] <- unique(df$case)
map[, 2] <- paste0("[", 1:nrow(map), "]")

for(ii in 1:nrow(df)){
  
  idx <- which(map[, 1]==df$case[ii])
  df$case[ii] <- map[idx, 2]
  
}
df$network <- factor(x = df$network, levels = c("splice_aware", "carnival", "splice_unaware"))

# pdf(file = "output/boxplot_combined.pdf", width = 10, height = 6)
# p <- ggpaired(df, x = "network", y = "enrichment_score",
#               color = "network", line.color = "gray", line.size = 0,
#               palette = "jco", label = "case",
#               repel = TRUE, label.rectangle = FALSE)+
#   stat_compare_means(paired = TRUE, method = "t.test")
# plot(p)
# dev.off()

pdf(file = "output/boxplot_combined.pdf", width = 5, height = 3)
p <- ggpaired(df, x = "network", y = "enrichment_score",
              color = "network", line.color = "gray", line.size = 0,
              palette = c("darkblue", "darkred", "grey"),
              repel = TRUE, label.rectangle = FALSE)+
  stat_compare_means(paired = TRUE, method = "t.test")
plot(p)
dev.off()
