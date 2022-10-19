library("readr")
library("GSA")
library("httr")
library("clusterProfiler")
library("biomaRt")

gmt <- GSA.read.gmt("src/human_ontology.gmt")
gene_list <- gmt$genesets
names(gene_list) <- gmt$geneset.names

load(file = "../LINDA_Soft_Analysis/K562/kd/output/res_k562_kd.RData")

cases <- toupper(names(res_k562_kd))

perturbation_clusters <- read.delim("src/perturbation_clusters.txt")

cases_list <- list()
cases_list_names <- c()
for(ii in 1:length(cases)){
  
  phenotypes <- c()
  cnt <- 0
  for(jj in 1:nrow(perturbation_clusters)){
    
    nearby_genes <- unique(c(unlist(strsplit(x = perturbation_clusters$members[jj], 
                                             split = ",", fixed = TRUE)),
                             unlist(strsplit(x = perturbation_clusters$nearby_genes[jj], 
                                             split = ",", fixed = TRUE))))
    
    if(cases[ii]%in%nearby_genes){
      
      phenotypes <- unique(c(phenotypes, perturbation_clusters$best_description[jj]))
      cnt <- cnt + 1
      
    }
    
  }
  
  if(cnt > 0){
    
    cases_list_names <- c(cases_list_names, cases[ii])
    cases_list[[length(cases_list)+1]] <- phenotypes
    
  }
  
}
names(cases_list) <- cases_list_names


## Now manually annotate the phenotypes
phenotype_list <- list()
phenotype_list_names <- c()

mart=useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host="apr2019.archive.ensembl.org")
map=getBM(attributes=c("ensembl_gene_id","external_gene_name", "uniprot_gn"), mart=mart)

# HNRNPA2B1
df <- read.delim("src/histone_acetylation.txt", header=FALSE)
my_protein_ids <- sapply(strsplit(x = df$V1, split = ":", fixed = TRUE), "[", 2)
my_gene_ids <- map$external_gene_name[which(map$uniprot_gn%in%my_protein_ids)]
phenotype_list[[length(phenotype_list)+1]] <- my_gene_ids
phenotype_list_names <- c(phenotype_list_names, "HNRNPA2B1 - Histone Acetylation")

# HNRNPC
df <- read.delim("src/nua4_histone_acetyltransverase_complex.txt", header=FALSE)
my_protein_ids <- sapply(strsplit(x = df$V1, split = ":", fixed = TRUE), "[", 2)
my_gene_ids <- map$external_gene_name[which(map$uniprot_gn%in%my_protein_ids)]
phenotype_list[[length(phenotype_list)+1]] <- my_gene_ids
phenotype_list_names <- c(phenotype_list_names, "HNRNPC - NuA4 histone acetyltransferase complex")

# HNRNPU
df <- read.delim("src/spliceosome.txt", header=FALSE)
my_protein_ids <- sapply(strsplit(x = df$V1, split = ":", fixed = TRUE), "[", 2)
my_gene_ids <- map$external_gene_name[which(map$uniprot_gn%in%my_protein_ids)]
phenotype_list[[length(phenotype_list)+1]] <- my_gene_ids
phenotype_list_names <- c(phenotype_list_names, "HNRNPU - Spliceosome")

# MAGOH
df <- read.delim("src/nonsense_mediated_decay.txt", header=FALSE)
my_protein_ids <- sapply(strsplit(x = df$V1, split = ":", fixed = TRUE), "[", 2)
my_gene_ids <- map$external_gene_name[which(map$uniprot_gn%in%my_protein_ids)]
phenotype_list[[length(phenotype_list)+1]] <- my_gene_ids
phenotype_list_names <- c(phenotype_list_names, "MAGOH - Non-sense Mediated Decay")

# NCBP2 (GO:0006370)
df <- read.delim("src/mrna_capping.txt", header=FALSE)
my_protein_ids <- sapply(strsplit(x = df$V1, split = ":", fixed = TRUE), "[", 2)
my_gene_ids <- map$external_gene_name[which(map$uniprot_gn%in%my_protein_ids)]
phenotype_list[[length(phenotype_list)+1]] <- my_gene_ids
phenotype_list_names <- c(phenotype_list_names, "NCBP2 - mRNA Capping")

# PABPN1 (GO:0000178)
df <- read.delim("src/exosome.txt", header=FALSE)
my_protein_ids <- sapply(strsplit(x = df$V1, split = ":", fixed = TRUE), "[", 2)
my_gene_ids <- map$external_gene_name[which(map$uniprot_gn%in%my_protein_ids)]
phenotype_list[[length(phenotype_list)+1]] <- my_gene_ids
phenotype_list_names <- c(phenotype_list_names, "PABPN1 - Exosome and mRNA turnover")

# PAPOLA (GO:0016592)
df <- read.delim("src/mediator_complex.txt", header=FALSE)
my_protein_ids <- sapply(strsplit(x = df$V1, split = ":", fixed = TRUE), "[", 2)
my_gene_ids <- map$external_gene_name[which(map$uniprot_gn%in%my_protein_ids)]
phenotype_list[[length(phenotype_list)+1]] <- my_gene_ids
phenotype_list_names <- c(phenotype_list_names, "PAPOLA - Mediator Complex")

# PCBP1 (GO:0008180)
df <- read.delim("src/cop9_signalosome.txt", header=FALSE)
my_protein_ids <- sapply(strsplit(x = df$V1, split = ":", fixed = TRUE), "[", 2)
my_gene_ids <- map$external_gene_name[which(map$uniprot_gn%in%my_protein_ids)]
phenotype_list[[length(phenotype_list)+1]] <- my_gene_ids
phenotype_list_names <- c(phenotype_list_names, "PCBP1 - COP9 Signalosome")

# POLR2G (GO:0006289)
df <- read.delim("src/nucleotide_excision_repair.txt", header=FALSE)
my_protein_ids <- sapply(strsplit(x = df$V1, split = ":", fixed = TRUE), "[", 2)
my_gene_ids <- map$external_gene_name[which(map$uniprot_gn%in%my_protein_ids)]
phenotype_list[[length(phenotype_list)+1]] <- my_gene_ids
phenotype_list_names <- c(phenotype_list_names, "POLR2G - Nucleotide Excision Repair")

# PPIL4 (GO:0000398)
df <- read.delim("src/spliceosome.txt", header=FALSE)
my_protein_ids <- sapply(strsplit(x = df$V1, split = ":", fixed = TRUE), "[", 2)
my_gene_ids <- map$external_gene_name[which(map$uniprot_gn%in%my_protein_ids)]
phenotype_list[[length(phenotype_list)+1]] <- my_gene_ids
phenotype_list_names <- c(phenotype_list_names, "PPIL4 - Spliceosome")

# PRPF6 (GO:0000398)
df <- read.delim("src/spliceosome.txt", header=FALSE)
my_protein_ids <- sapply(strsplit(x = df$V1, split = ":", fixed = TRUE), "[", 2)
my_gene_ids <- map$external_gene_name[which(map$uniprot_gn%in%my_protein_ids)]
phenotype_list[[length(phenotype_list)+1]] <- my_gene_ids
phenotype_list_names <- c(phenotype_list_names, "PRPF6 - Spliceosome")

# PRPF8 (GO:0000398)
df <- read.delim("src/spliceosome.txt", header=FALSE)
my_protein_ids <- sapply(strsplit(x = df$V1, split = ":", fixed = TRUE), "[", 2)
my_gene_ids <- map$external_gene_name[which(map$uniprot_gn%in%my_protein_ids)]
phenotype_list[[length(phenotype_list)+1]] <- my_gene_ids
phenotype_list_names <- c(phenotype_list_names, "PRPF8 - Spliceosome")

# PUF60 (GO:0000398)
df <- read.delim("src/spliceosome.txt", header=FALSE)
my_protein_ids <- sapply(strsplit(x = df$V1, split = ":", fixed = TRUE), "[", 2)
my_gene_ids <- map$external_gene_name[which(map$uniprot_gn%in%my_protein_ids)]
phenotype_list[[length(phenotype_list)+1]] <- my_gene_ids
phenotype_list_names <- c(phenotype_list_names, "PUF60 - Spliceosome")

# SF1
df <- read.delim("src/nua4_histone_acetyltransverase_complex.txt", header=FALSE)
my_protein_ids <- sapply(strsplit(x = df$V1, split = ":", fixed = TRUE), "[", 2)
my_gene_ids <- map$external_gene_name[which(map$uniprot_gn%in%my_protein_ids)]
phenotype_list[[length(phenotype_list)+1]] <- my_gene_ids
phenotype_list_names <- c(phenotype_list_names, "SF1 - NuA4 histone acetyltransferase complex")

# SMNDC1 (GO:0000398)
df <- read.delim("src/spliceosome.txt", header=FALSE)
my_protein_ids <- sapply(strsplit(x = df$V1, split = ":", fixed = TRUE), "[", 2)
my_gene_ids <- map$external_gene_name[which(map$uniprot_gn%in%my_protein_ids)]
phenotype_list[[length(phenotype_list)+1]] <- my_gene_ids
phenotype_list_names <- c(phenotype_list_names, "SMNDC1 - Spliceosome")

# SRSF1 (GO:0016592)
df <- read.delim("src/mediator_complex.txt", header=FALSE)
my_protein_ids <- sapply(strsplit(x = df$V1, split = ":", fixed = TRUE), "[", 2)
my_gene_ids <- map$external_gene_name[which(map$uniprot_gn%in%my_protein_ids)]
phenotype_list[[length(phenotype_list)+1]] <- my_gene_ids
phenotype_list_names <- c(phenotype_list_names, "SRSF1 - Mediator Complex")


## Save everything
names(phenotype_list) <- phenotype_list_names
save(phenotype_list, file = "phenotype_list.RData")

temp <- gene_list
nn1 <- names(temp)
nn2 <- names(phenotype_list)

gene_list <- list()
for(ii in 1:length(temp)){
  gene_list[[length(gene_list)+1]] <- temp[[ii]]
}

for(ii in 1:length(phenotype_list)){
  gene_list[[length(gene_list)+1]] <- phenotype_list[[ii]]
}

names(gene_list) <- c(nn1, nn2)
save(gene_list, file = "gene_list.RData")

