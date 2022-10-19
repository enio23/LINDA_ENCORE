#! /usr/bin/env Rscript

set.seed(1234)

library("readr")
library("org.Hs.eg.db")
library("biomaRt")
library("stringr")
library("dplyr")
library("edgeR")
library("dorothea")
library("foreach")
library("doParallel")
library("piano")
library("BiRewire")

dir.create("output")

source("estimate_significance.R")

mart=useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host="apr2019.archive.ensembl.org")
ensg2symbol=getBM(attributes=c("ensembl_gene_id","external_gene_name"), mart=mart)

splicing_factors <- read.table("reactome_mrna_splicing_factors.txt", quote="\"", comment.char="")

data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B","C"))

## HepG2 samples
dir.create("output/hepg2")

encode_samples <- read.delim(file = "encode_samples.txt", header = FALSE)

encode_samples <- encode_samples[intersect(x = which(encode_samples$V4=="hepg2"), y = which(encode_samples$V1%in%splicing_factors$V1)), 1:3]

koGenes <- unique(encode_samples$V1)

registerDoParallel(cores=50)

foreach(ii = 1:length(koGenes)) %dopar% {

  currGene <- koGenes[ii]

  idxKO <- intersect(x = which(encode_samples$V1==currGene), y = which(encode_samples$V2=="ko"))
  idxCtrl <- intersect(x = which(encode_samples$V1==currGene), y = which(encode_samples$V2=="ctrl"))

  download.file(url = encode_samples$V3[idxKO], destfile = paste0(getwd(), "/temp_ko_", ii, ".txt"))
  download.file(url = encode_samples$V3[idxCtrl], destfile = paste0(getwd(), "/temp_ctrl_", ii, ".txt"))

  temp_ko <- read.table(file = paste0("temp_ko_", ii, ".txt"), sep = " ")
  temp_ctrl <- read.table(file = paste0("temp_ctrl_", ii, ".txt"), sep = " ")

  # Download KO
  experiments <- temp_ko$V1[which(grepl(pattern = "gene_quantifications_GRCh38.tsv.gz", x = temp_ko$V1, fixed = TRUE))]
  experiments <- gsub(pattern = "\t", replacement = " ", x = experiments)
  experiments <- unlist(x = strsplit(x = experiments, split = " ", fixed = TRUE))
  experiments <- experiments[which(grepl(pattern = ".tsv", x = experiments, fixed = TRUE))]
  experiments <- experiments[which(grepl(pattern = "GRCh38", x = experiments, fixed = TRUE))]
  cnt <- 1
  for(jj in 1:length(experiments)){
    if(grepl(pattern = "gene_quantifications", x = experiments[jj])){
      download.file(url = experiments[jj], destfile = paste0(getwd(), "/temp", ii, ".tsv"))
      temp <- read_tsv(file = paste0(getwd(), "/temp", ii, ".tsv"))
      if("gene_id"%in%colnames(temp)){
        write_tsv(x = temp, file = paste0("ko", cnt, "_", ii, ".tsv"))
        cnt <- cnt + 1
      }
      file.remove(paste0(getwd(), "/temp", ii, ".tsv"))
    }
  }

  # Download Ctrl
  experiments <- temp_ctrl$V1[which(grepl(pattern = "gene_quantifications_GRCh38.tsv.gz", x = temp_ctrl$V1, fixed = TRUE))]
  experiments <- gsub(pattern = "\t", replacement = " ", x = experiments)
  experiments <- unlist(x = strsplit(x = experiments, split = " ", fixed = TRUE))
  experiments <- experiments[which(grepl(pattern = ".tsv", x = experiments, fixed = TRUE))]
  experiments <- experiments[which(grepl(pattern = "GRCh38", x = experiments, fixed = TRUE))]
  cnt <- 1
  for(jj in 1:length(experiments)){
    if(grepl(pattern = "gene_quantifications", x = experiments[jj])){
      download.file(url = experiments[jj], destfile = paste0(getwd(), "/temp", ii, ".tsv"))
      temp <- read_tsv(file = paste0(getwd(), "/temp", ii, ".tsv"))
      if("gene_id"%in%colnames(temp)){
        write_tsv(x = temp, file = paste0("ctrl", cnt, "_", ii, ".tsv"))
        cnt <- cnt + 1
      }
      file.remove(paste0(getwd(), "/temp", ii, ".tsv"))
    }
  }

  # Build Counts data matrix
  ko1 <- read.table(file = paste0("ko1_", ii, ".tsv"), header = TRUE)
  ko2 <- read.table(file = paste0("ko2_", ii, ".tsv"), header = TRUE)
  ctrl1 <- read.table(file = paste0("ctrl1_", ii, ".tsv"), header = TRUE)
  ctrl2 <- read.table(file = paste0("ctrl2_", ii, ".tsv"), header = TRUE)

  commonGenes <- intersect(x = intersect(x = ko1$gene_id, y = ko2$gene_id),
                           y = intersect(x = ctrl1$gene_id, y = ctrl2$gene_id))

  if(length(commonGenes)>0){

    ko1 <- ko1[which(ko1$gene_id%in%commonGenes), ]
    ko2 <- ko2[which(ko2$gene_id%in%commonGenes), ]
    ctrl1 <- ctrl1[which(ctrl1$gene_id%in%commonGenes), ]
    ctrl2 <- ctrl2[which(ctrl2$gene_id%in%commonGenes), ]

    ko1 <- ko1[order(ko1$gene_id), ]
    ko2 <- ko2[order(ko2$gene_id), ]
    ctrl1 <- ctrl1[order(ctrl1$gene_id), ]
    ctrl2 <- ctrl2[order(ctrl2$gene_id), ]

    idx2rem <- which(duplicated(ko1$gene_id)); if(length(idx2rem)>0){ko1 <- ko1[-idx2rem, ]}
    idx2rem <- which(duplicated(ko2$gene_id)); if(length(idx2rem)>0){ko2 <- ko2[-idx2rem, ]}
    idx2rem <- which(duplicated(ctrl1$gene_id)); if(length(idx2rem)>0){ctrl1 <- ctrl1[-idx2rem, ]}
    idx2rem <- which(duplicated(ctrl2$gene_id)); if(length(idx2rem)>0){ctrl2 <- ctrl2[-idx2rem, ]}

    df <- matrix(data = , nrow = nrow(ko1), ncol = 4)
    df[, 1] <- ko1$expected_count
    df[, 2] <- ko2$expected_count
    df[, 3] <- ctrl1$expected_count
    df[, 4] <- ctrl2$expected_count
    df <- as.data.frame(df)
    rownames(df) <- sapply(strsplit(x = ko1$gene_id, split = ".", fixed = TRUE), '[', 1)
    colnames(df) <- c(paste0("ko_", 1:2), paste0("ctrl_", 1:2))

    uGenes <- unique(ensg2symbol$external_gene_name)
    dfMapped <- matrix(data = , nrow = length(uGenes), ncol = ncol(df))
    for(ll in 1:length(uGenes)){
      ensg <- ensg2symbol$ensembl_gene_id[which(ensg2symbol$external_gene_name==uGenes[ll])]
      for(mm in 1:ncol(df)){
        dfMapped[ll, mm] <- median(x = df[which(rownames(df)%in%ensg), mm], na.rm = TRUE)
      }
    }
    rownames(dfMapped) <- uGenes
    colnames(dfMapped) <- colnames(df)

    dfMapped <- dfMapped[complete.cases(dfMapped), ]

    # Do differential gene expression analysis
    conditions<-factor(c("ko", "ko", "ctrl", "ctrl"))
    design <- model.matrix(~ conditions)

    y <- DGEList(counts=dfMapped, group=conditions)
    keep <- filterByExpr(y, group=conditions)
    y <- y[keep, , keep.lib.sizes=FALSE]
    y <- calcNormFactors(y)
    y <- estimateGLMCommonDisp(y, design, verbose=TRUE)
    y <- estimateGLMTagwiseDisp(y, design)

    fit <- glmQLFit(y,design,robust=TRUE)

    res=glmQLFTest(fit, coef=2)
    ttop=as.data.frame(topTags(res,n=nrow(df)))
    ttop$ID <- rownames(ttop)

    save(ttop, file = paste0("output/hepg2/ttop_", tolower(currGene), "_hepg2.RData"))

    # TF Activities
    ss <- ttop$logFC
    names(ss) <- ttop$ID

    input.scores <- estimate_significance(expr = as.matrix(ss), regulons = regulons, nperm = 1000)

    save(input.scores, file = paste0("output/hepg2/tf_act_", tolower(koGenes)[ii], "_hepg2.RData"))

    file.remove(paste0("ko1_", ii, ".tsv"))
    file.remove(paste0("ko2_", ii, ".tsv"))
    file.remove(paste0("ctrl1_", ii, ".tsv"))
    file.remove(paste0("ctrl2_", ii, ".tsv"))
    file.remove(paste0("temp_ko_", ii, ".txt"))
    file.remove(paste0("temp_ctrl_", ii, ".txt"))

  }

}

ff <- list.files()
ff <- ff[c(which(grepl(pattern = ".tsv", x = ff, fixed = TRUE)), which(grepl(pattern = "temp_", x = ff, fixed = TRUE)))]
file.remove(ff)



## K562 samples
dir.create("output/k562")

encode_samples <- read.delim(file = "encode_samples.txt", header = FALSE)

encode_samples <- encode_samples[intersect(x = which(encode_samples$V4=="k562"), y = which(encode_samples$V1%in%splicing_factors$V1)), 1:3]

koGenes <- unique(encode_samples$V1)

registerDoParallel(cores=50)

foreach(ii = 1:length(koGenes)) %dopar% {
  
  currGene <- koGenes[ii]
  
  idxKO <- intersect(x = which(encode_samples$V1==currGene), y = which(encode_samples$V2=="ko"))
  idxCtrl <- intersect(x = which(encode_samples$V1==currGene), y = which(encode_samples$V2=="ctrl"))
  
  download.file(url = encode_samples$V3[idxKO], destfile = paste0(getwd(), "/temp_ko_", ii, ".txt"))
  download.file(url = encode_samples$V3[idxCtrl], destfile = paste0(getwd(), "/temp_ctrl_", ii, ".txt"))
  
  temp_ko <- read.table(file = paste0("temp_ko_", ii, ".txt"), sep = " ")
  temp_ctrl <- read.table(file = paste0("temp_ctrl_", ii, ".txt"), sep = " ")
  
  # Download KO
  experiments <- temp_ko$V1[which(grepl(pattern = "gene_quantifications_GRCh38.tsv.gz", x = temp_ko$V1, fixed = TRUE))]
  experiments <- gsub(pattern = "\t", replacement = " ", x = experiments)
  experiments <- unlist(x = strsplit(x = experiments, split = " ", fixed = TRUE))
  experiments <- experiments[which(grepl(pattern = ".tsv", x = experiments, fixed = TRUE))]
  experiments <- experiments[which(grepl(pattern = "GRCh38", x = experiments, fixed = TRUE))]
  cnt <- 1
  for(jj in 1:length(experiments)){
    if(grepl(pattern = "gene_quantifications", x = experiments[jj])){
      download.file(url = experiments[jj], destfile = paste0(getwd(), "/temp", ii, ".tsv"))
      temp <- read_tsv(file = paste0(getwd(), "/temp", ii, ".tsv"))
      if("gene_id"%in%colnames(temp)){
        write_tsv(x = temp, file = paste0("ko", cnt, "_", ii, ".tsv"))
        cnt <- cnt + 1
      }
      file.remove(paste0(getwd(), "/temp", ii, ".tsv"))
    }
  }
  
  # Download Ctrl
  experiments <- temp_ctrl$V1[which(grepl(pattern = "gene_quantifications_GRCh38.tsv.gz", x = temp_ctrl$V1, fixed = TRUE))]
  experiments <- gsub(pattern = "\t", replacement = " ", x = experiments)
  experiments <- unlist(x = strsplit(x = experiments, split = " ", fixed = TRUE))
  experiments <- experiments[which(grepl(pattern = ".tsv", x = experiments, fixed = TRUE))]
  experiments <- experiments[which(grepl(pattern = "GRCh38", x = experiments, fixed = TRUE))]
  cnt <- 1
  for(jj in 1:length(experiments)){
    if(grepl(pattern = "gene_quantifications", x = experiments[jj])){
      download.file(url = experiments[jj], destfile = paste0(getwd(), "/temp", ii, ".tsv"))
      temp <- read_tsv(file = paste0(getwd(), "/temp", ii, ".tsv"))
      if("gene_id"%in%colnames(temp)){
        write_tsv(x = temp, file = paste0("ctrl", cnt, "_", ii, ".tsv"))
        cnt <- cnt + 1
      }
      file.remove(paste0(getwd(), "/temp", ii, ".tsv"))
    }
  }
  
  # Build Counts data matrix
  ko1 <- read.table(file = paste0("ko1_", ii, ".tsv"), header = TRUE)
  ko2 <- read.table(file = paste0("ko2_", ii, ".tsv"), header = TRUE)
  ctrl1 <- read.table(file = paste0("ctrl1_", ii, ".tsv"), header = TRUE)
  ctrl2 <- read.table(file = paste0("ctrl2_", ii, ".tsv"), header = TRUE)
  
  commonGenes <- intersect(x = intersect(x = ko1$gene_id, y = ko2$gene_id), 
                           y = intersect(x = ctrl1$gene_id, y = ctrl2$gene_id))
  
  if(length(commonGenes)>0){
    
    ko1 <- ko1[which(ko1$gene_id%in%commonGenes), ]
    ko2 <- ko2[which(ko2$gene_id%in%commonGenes), ]
    ctrl1 <- ctrl1[which(ctrl1$gene_id%in%commonGenes), ]
    ctrl2 <- ctrl2[which(ctrl2$gene_id%in%commonGenes), ]
    
    ko1 <- ko1[order(ko1$gene_id), ]
    ko2 <- ko2[order(ko2$gene_id), ]
    ctrl1 <- ctrl1[order(ctrl1$gene_id), ]
    ctrl2 <- ctrl2[order(ctrl2$gene_id), ]
    
    idx2rem <- which(duplicated(ko1$gene_id)); if(length(idx2rem)>0){ko1 <- ko1[-idx2rem, ]}
    idx2rem <- which(duplicated(ko2$gene_id)); if(length(idx2rem)>0){ko2 <- ko2[-idx2rem, ]}
    idx2rem <- which(duplicated(ctrl1$gene_id)); if(length(idx2rem)>0){ctrl1 <- ctrl1[-idx2rem, ]}
    idx2rem <- which(duplicated(ctrl2$gene_id)); if(length(idx2rem)>0){ctrl2 <- ctrl2[-idx2rem, ]}
    
    df <- matrix(data = , nrow = nrow(ko1), ncol = 4)
    df[, 1] <- ko1$expected_count
    df[, 2] <- ko2$expected_count
    df[, 3] <- ctrl1$expected_count
    df[, 4] <- ctrl2$expected_count
    df <- as.data.frame(df)
    rownames(df) <- sapply(strsplit(x = ko1$gene_id, split = ".", fixed = TRUE), '[', 1)
    colnames(df) <- c(paste0("ko_", 1:2), paste0("ctrl_", 1:2))
    
    uGenes <- unique(ensg2symbol$external_gene_name)
    dfMapped <- matrix(data = , nrow = length(uGenes), ncol = ncol(df))
    for(ll in 1:length(uGenes)){
      ensg <- ensg2symbol$ensembl_gene_id[which(ensg2symbol$external_gene_name==uGenes[ll])]
      for(mm in 1:ncol(df)){
        dfMapped[ll, mm] <- median(x = df[which(rownames(df)%in%ensg), mm], na.rm = TRUE)
      }
    }
    rownames(dfMapped) <- uGenes
    colnames(dfMapped) <- colnames(df)
    
    dfMapped <- dfMapped[complete.cases(dfMapped), ]
    
    # Do differential gene expression analysis
    conditions<-factor(c("ko", "ko", "ctrl", "ctrl"))
    design <- model.matrix(~ conditions)
    
    y <- DGEList(counts=dfMapped, group=conditions)
    keep <- filterByExpr(y, group=conditions)
    y <- y[keep, , keep.lib.sizes=FALSE]
    y <- calcNormFactors(y)
    y <- estimateGLMCommonDisp(y, design, verbose=TRUE)
    y <- estimateGLMTagwiseDisp(y, design)
    
    fit <- glmQLFit(y,design,robust=TRUE)
    
    res=glmQLFTest(fit, coef=2)
    ttop=as.data.frame(topTags(res,n=nrow(df)))
    ttop$ID <- rownames(ttop)
    
    save(ttop, file = paste0("output/k562/ttop_", tolower(currGene), "_k562.RData"))
    
    # TF Activities
    ss <- ttop$logFC
    names(ss) <- ttop$ID
    
    input.scores <- estimate_significance(expr = as.matrix(ss), regulons = regulons, nperm = 1000)
    
    save(input.scores, file = paste0("output/k562/tf_act_", tolower(koGenes)[ii], "_k562.RData"))
    
    file.remove(paste0("ko1_", ii, ".tsv"))
    file.remove(paste0("ko2_", ii, ".tsv"))
    file.remove(paste0("ctrl1_", ii, ".tsv"))
    file.remove(paste0("ctrl2_", ii, ".tsv"))
    file.remove(paste0("temp_ko_", ii, ".txt"))
    file.remove(paste0("temp_ctrl_", ii, ".txt"))
    
  }
  
}

ff <- list.files()
ff <- ff[c(which(grepl(pattern = ".tsv", x = ff, fixed = TRUE)), which(grepl(pattern = "temp_", x = ff, fixed = TRUE)))]
file.remove(ff)
