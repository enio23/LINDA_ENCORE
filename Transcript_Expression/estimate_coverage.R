library(LINDA)
library(foreach)
library(doParallel)

load(file = system.file("extdata", "digger_human_transcripts.RData", package = "LINDA"))

domains <- unique(c(paste0(bg$gene_source, "_", bg$pfam_source),
                    paste0(bg$gene_target, "_", bg$pfam_target)))

map_table <- matrix(data = , nrow = length(domains), ncol = 2)
map_table[, 1] <- domains
for(ii in 1:nrow(map_table)){
  
  # print(paste0("Step ---- ", ii, "/", nrow(map_table)))
  
  transcripts <- c()
  
  idx <- which(paste0(bg$gene_source, "_", bg$pfam_source)==domains[ii])
  if(length(idx)>0){
    
    transcripts <- c(transcripts, unique(unlist(strsplit(x = bg$id_source[idx], split = "_", fixed = TRUE))))
    
  }
  
  idx <- which(paste0(bg$gene_target, "_", bg$pfam_target)==domains[ii])
  if(length(idx)>0){
    
    transcripts <- c(transcripts, unique(unlist(strsplit(x = bg$id_target[idx], split = "_", fixed = TRUE))))
    
  }
  
  map_table[ii, 2] <- paste(unique(transcripts), collapse = "_")
  
}
colnames(map_table) <- c("domains", "transcripts")
map_table <- as.data.frame(map_table)
map_table <- map_table[-which(map_table$transcripts=="NA"), ]

# hepg2
registerDoParallel(cores=20)

ff <- list.files(path = "output/hepg2/")
ff <- ff[which(grepl(pattern = "ttop_", x = ff, fixed = TRUE))]
perc1 <- rep(0, length(ff))
foreach(ii = 1:length(ff)) %dopar% {
  
  print(paste0("Step ---- ", ii, "/", length(ff)))
  
  load(file = paste0("output/hepg2/", ff[ii]))
  transcripts <- ttop$ID
  
  if(file.exists(paste0("../Gene_Expression/output/hepg2/", ff[ii]))){
    
    load(file = paste0("../Gene_Expression/output/hepg2/", ff[ii]))
    genes <- sapply(strsplit(x = map_table$domains, split = "_", fixed = TRUE), "[", 1)
    tmp <- map_table[which(genes%in%ttop$ID), ]
    
    cases <- rep(0, nrow(tmp))
    for(jj in 1:length(transcripts)){
      
      idx <- which(grepl(pattern = transcripts[jj], x = tmp$transcripts, fixed = TRUE))
      if(length(idx)>0){
        
        cases[idx] <- 1
        
      }
      
    }
    
    perc1[ii] <- sum(cases)/length(cases)
    
  }
  
}

# perc1 <- c(0.7613316, 0.7417979, 0.7222643, 0.7293486, 0.7363643, 0.7328565,
#            0.7403535, 0.750877, 0.7412477, 0.7522526, 0.7334755, 0.7227457,
#            0.7366394, 0.7363643, 0.7429672, 0.7253594, 0.7424857, 0.7409038,
#            0.7393906, 0.7387028, 0.7432423, 0.7221267, 0.7415916, 0.7456496,
#            0.7378087, 0.7349199, 0.7308618, 0.7348511, 0.7382901, 0.7541096,
#            0.753697, 0.7304491, 0.7558291, 0.7395282, 0.743036, 0.7411789,
#            0.7477818, 0.7318247, 0.7250155, 0.7406287, 0.750533, 0.7517023,
#            0.7435862, 0.735195, 0.7354701, 0.7294862, 0.7464062, 0.7303804)

perc1 <- c(0.9597445, 0.9612926, 0.9553719, 0.9548141, 0.9577045, 0.9594361,
           0.9598343, 0.9618259, 0.9577427, 0.9568212, 0.9573254, 0.9544371,
           0.9562827, 0.959649, 0.9537351, 0.9578136, 0.9608954, 0.962031,
           0.9601371, 0.9600181, 0.9506353, 0.9637119, 0.9597711, 0.9605751,
           0.9560231, 0.9596098, 0.9567103, 0.9569873, 0.9596004, 0.9568054,
           0.9566166, 0.9637084, 0.9552386, 0.9595371, 0.9556888, 0.960389,
           0.9586063, 0.9524765, 0.9566036, 0.959936, 0.9560838, 0.9601147,
           0.9612171, 0.9600145, 0.9608632, 0.9594896, 0.9528857)

# k562
registerDoParallel(cores=20)

ff <- list.files(path = "output/k562/")
ff <- ff[which(grepl(pattern = "ttop_", x = ff, fixed = TRUE))]
perc2 <- rep(0, length(ff))
foreach(ii = 1:length(ff)) %dopar% {
  
  print(paste0("Step ---- ", ii, "/", length(ff)))
  
  load(file = paste0("output/k562/", ff[ii]))
  transcripts <- ttop$ID
  
  if(file.exists(paste0("../Gene_Expression/output/k562/", ff[ii]))){
    
    load(file = paste0("../Gene_Expression/output/k562/", ff[ii]))
    genes <- sapply(strsplit(x = map_table$domains, split = "_", fixed = TRUE), "[", 1)
    tmp <- map_table[which(genes%in%ttop$ID), ]
    
    cases <- rep(0, nrow(tmp))
    for(jj in 1:length(transcripts)){
      
      idx <- which(grepl(pattern = transcripts[jj], x = tmp$transcripts, fixed = TRUE))
      if(length(idx)>0){
        
        cases[idx] <- 1
        
      }
      
    }
    
    perc2[ii] <- sum(cases)/length(cases)
    
  }
  
}

# perc2 <- c(0.7241213, 0.6758374, 0.676594, 0.682234, 0.6964716, 0.6889057,
#            0.6945457, 0.6923447, 0.6700598, 0.6955086, 0.6915881, 0.7023867,
#            0.7064447, 0.7025242, 0.7122223, 0.6846413, 0.685054, 0.7127725,
#            0.6818901, 0.7085769, 0.7124286, 0.6929638, 0.693514, 0.6949584,
#            0.6658642, 0.6858106, 0.6857418, 0.6713667, 0.7004608, 0.6942706,
#            0.6915194, 0.6899374, 0.7358828, 0.7135291, 0.6934452, 0.6987413,
#            0.7207511, 0.6889057, 0.6697847, 0.6931013, 0.7005296, 0.694133,
#            0.6884242, 0.6849852, 0.6764564, 0.6582296, 0.6642135)
perc2 <- c(0.9370644, 0.9466615, 0.9473787, 0.9473278, 0.9554729, 0.9448204,
           0.9510942, 0.9538712, 0.9427873, 0.9527753, 0.9493369, 0.950188,
           0.9542809, 0.9517755, 0.9462028, 0.9485181, 0.94437, 0.9563514,
           0.9507704, 0.9500513, 0.9546428, 0.9496238, 0.9508525, 0.9501282,
           0.9394296, 0.9450507, 0.9463036, 0.9433483, 0.9496843, 0.9522448,
           0.9388852, 0.9516144, 0.9592244, 0.9600599, 0.9472486, 0.9522592,
           0.9480662, 0.9348451, 0.9347426, 0.9537241, 0.9511182, 0.9505186,
           0.9488653, 0.9488117, 0.9472865, 0.9282006, 0.9317872)
