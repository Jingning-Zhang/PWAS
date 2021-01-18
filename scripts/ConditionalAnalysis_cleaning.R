suppressMessages(library("readr"))

##############################
# 1. Analyze by each tissue

tissue_list <- readLines("../GTex_V7_tissue_list.txt")


for (tissue in tissue_list){
    
    load(paste0("../Results/ConditionalAnalysis/RDat/",tissue,".RDat"))
    PWAS_hit=dat.sentinel.pwas$ID

    res <- list()
    for (i in 1:length(PWAS_hit)){
        
        TWAS_hit <- twas.hit[i]
        TWAS_p <- twas.p[i]
        Dist_of_hits <- dist[i]
        Corr_of_hits <- corr[i]
        PcT_p <- PcT.p[i]
        TcP_p <- TcP.p[i]
        
        res[[PWAS_hit[i]]] <- list(TWAS_hit=TWAS_hit, TWAS_p=TWAS_p,
                                   Dist_of_hits=Dist_of_hits, Corr_of_hits=Corr_of_hits,
                                   PcT_p=PcT_p, TcP_p=TcP_p)
    }

    p.twas <- 0.05/readRDS("../GTex_V7_n_gene.rds")[tissue]
    
    TWAS_p <- numeric()
    TWAS_hit <- character()
    Dist_of_hits <- integer()
    Corr_of_hits <- numeric()
    PcT_p <- numeric()
    TcP_p <- numeric()
    sig <- logical()
    for (i in 1:length(res)){
        m <- (res[[i]]$TWAS_p < p.twas)
        m[is.na(m)] <- FALSE
        sig[i] <- m
        
        TWAS_p[i] <- res[[i]]$TWAS_p
        TWAS_hit[i] <- res[[i]]$TWAS_hit
        Dist_of_hits[i] <- res[[i]]$Dist_of_hits
        Corr_of_hits[i] <- res[[i]]$Corr_of_hits
        PcT_p[i] <- res[[i]]$PcT_p
        TcP_p[i] <- res[[i]]$TcP_p
    }
    
    a <- data.frame(PWAS_hit=dat.sentinel.pwas$ID,
                    PWAS_p=signif(dat.sentinel.pwas$PWAS.P, 3),
                    TWAS_hit=ifelse(sig,paste0(TWAS_hit,"*"),TWAS_hit),
                    TWAS_p=ifelse(sig,paste0(signif(TWAS_p, 3),"*"),TWAS_p),
                    Dist_of_hits=Dist_of_hits,
                    Corr_of_hits=signif(Corr_of_hits,3),
                    PcT_p=signif(PcT_p,3),
                    TcP_p=signif(TcP_p,3))
    
    write_tsv(a, paste0("../Results/ConditionalAnalysis/Table/",tissue,".txt"))
    
}

##############################
# 2. All-tissue analysis

tissue_list <- readLines("../GTex_V7_tissue_list.txt")

# get all the regional sentinel PWAS genes
load(paste0("../Results/ConditionalAnalysis/RDat/",tissue_list[1],".RDat"))
PWAS_hit=dat.sentinel.pwas$ID

# multiple testing correction
p.twas <- 0.05/sum(readRDS("../GTex_V7_n_gene.rds"))

# for each regional sentinel PWAS gene, pick out the most significant nearby TWAS gene across all tissues
res <- list()
for (i in 1:length(PWAS_hit)){
    
    TWAS_hit=character()
    TWAS_p=numeric()
    Dist_of_hits=numeric()
    Corr_of_hits=numeric()
    PcT_p=numeric()
    TcP_p=numeric()
    
    for (j in 1:length(tissue_list)){
        tissue <- tissue_list[j]
        
        load(paste0("../Results/ConditionalAnalysis/RDat/",tissue,".RDat"))
        TWAS_hit[j] <- twas.hit[i]
        TWAS_p[j] <- twas.p[i]
        Dist_of_hits[j] <- dist[i]
        Corr_of_hits[j] <- corr[i]
        PcT_p[j] <- PcT.p[i]
        TcP_p[j] <- TcP.p[i]
    }
    
    res[[PWAS_hit[i]]] <- list(TWAS_hit=TWAS_hit, TWAS_p=TWAS_p,
                               Dist_of_hits=Dist_of_hits, Corr_of_hits=Corr_of_hits,
                               PcT_p=PcT_p, TcP_p=TcP_p)
    
}

min_TWAS_p <- numeric()
min_tissue <- character()
min_TWAS_hit <- character()
min_Dist_of_hits <- integer()
min_Corr_of_hits <- numeric()
min_PcT_p <- numeric()
min_TcP_p <- numeric()
sig_tiss <- character()
N_tiss <- integer()
for (i in 1:length(res)){
    m <- (res[[i]]$TWAS_p < p.twas)
    m[is.na(m)] <- FALSE
    if(sum(m)>0){
        sig_tiss[i] <- paste(tissue_list[m],collapse=",")
        N_tiss[i] <- length(tissue_list[m])
        min_TWAS_p[i] <- min(res[[i]]$TWAS_p[m])
        min_tissue[i] <- tissue_list[m][which.min(res[[i]]$TWAS_p[m])]
        min_TWAS_hit[i] <- res[[i]]$TWAS_hit[m][which.min(res[[i]]$TWAS_p[m])]
        min_Dist_of_hits[i] <- res[[i]]$Dist_of_hits[m][which.min(res[[i]]$TWAS_p[m])]
        min_Corr_of_hits[i] <- res[[i]]$Corr_of_hits[m][which.min(res[[i]]$TWAS_p[m])]
        min_PcT_p[i] <- res[[i]]$PcT_p[m][which.min(res[[i]]$TWAS_p[m])]
        min_TcP_p[i] <- res[[i]]$TcP_p[m][which.min(res[[i]]$TWAS_p[m])]
    }else{
        sig_tiss[i] <- NA
        N_tiss[i] <- 0
        min_TWAS_p[i] <- min(res[[i]]$TWAS_p, na.rm = T)
        min_tissue[i] <- tissue_list[which.min(res[[i]]$TWAS_p)]
        min_TWAS_hit[i] <- res[[i]]$TWAS_hit[which.min(res[[i]]$TWAS_p)]
        min_Dist_of_hits[i] <- res[[i]]$Dist_of_hits[which.min(res[[i]]$TWAS_p)]
        min_Corr_of_hits[i] <- res[[i]]$Corr_of_hits[which.min(res[[i]]$TWAS_p)]
        min_PcT_p[i] <- res[[i]]$PcT_p[which.min(res[[i]]$TWAS_p)]
        min_TcP_p[i] <- res[[i]]$TcP_p[which.min(res[[i]]$TWAS_p)]
    }
}

a <- data.frame(min_TWAS_Tissue=ifelse(N_tiss!=0,paste0(min_tissue,"*"),min_tissue),
                min_TWAS_TWAS_hit=ifelse(N_tiss!=0,paste0(min_TWAS_hit,"*"),min_TWAS_hit),
                min_TWAS_p=ifelse(N_tiss!=0,paste0(signif(min_TWAS_p,3),"*"),min_TWAS_p),
                Dist_of_hits=min_Dist_of_hits,
                Corr_of_hits=signif(min_Corr_of_hits,3),
                PcT_p=signif(min_PcT_p,3),
                TcP_p=signif(min_TcP_p,3),
                N_significant_tissues_in_TWAS=N_tiss,
                all_significant_tissues_in_TWAS=sig_tiss)

write_tsv(a, paste0("../Results/ConditionalAnalysis/Table/all-tissue.txt"))





