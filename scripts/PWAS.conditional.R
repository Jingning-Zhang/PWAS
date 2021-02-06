suppressMessages(library("optparse"))
suppressMessages(library("readr"))
suppressMessages(library("stringr"))
suppressMessages(library("dplyr"))
suppressMessages(library("plink2R"))

option_list = list(
    make_option("--PWAS", action="store", default=NA, type='character',
                help="Path to PWAS output table. [required]"),
    make_option("--tissue_list", action="store", default=NA, type='character',
                help="Path to tissue list for which the conditional analysis is performed [required]"),
    make_option("--TWAS", action="store", default=NA, type='character',
                help="Path to the folder containing TWAS output tables. [required]"),
    make_option("--out", action="store", default=NA, type='character',
                help="Path to output files [required]"),
    make_option("--imputed_P", action="store", default=NA, type='character',
                help="Path to imputed plasma protein table is stored [required]"),
    make_option("--imputed_T", action="store", default=NA, type='character',
                help="Path to imputed gene expressions table is stored [required]")
    
)

opt = parse_args(OptionParser(option_list=option_list))

tissue_list <- readLines(opt$tissue_list)

dat.pwas <- suppressMessages(read_tsv(opt$PWAS))
p.pwas <- 0.05/nrow(dat.pwas)
dat.pwas <- dat.pwas[!(is.na(dat.pwas$PWAS.P)),]

## get sentinel PWAS genes
dat.sentinel <- tibble()
dat <- dat.pwas[dat.pwas$PWAS.P < p.pwas,]
for(i in 1:22){
    tmp0 <- dat[dat$CHR == i,]
    if(nrow(tmp0)==0){
        next()
    }
    tmp <- data.frame(file=tmp0$FILE,
                      pval=tmp0$PWAS.P,
                      tss=tmp0$P0,
                      num=1:nrow(tmp0))
    res <- integer()
    while (nrow(tmp)>0) {
        flag <- which(tmp$pval==min(tmp$pval))[1]
        res <- c(res, tmp$num[flag])
        tmp <- tmp[which((tmp$tss > tmp$tss[flag] + 1000000) | tmp$tss < tmp$tss[flag] - 1000000),]
    }
    tmp0 <- tmp0[res,]
    dat.sentinel <- rbind(dat.sentinel, tmp0)
}
dat.sentinel <- dat.sentinel[!is.na(dat.sentinel$PWAS.P),]
dat.sentinel.pwas <- dat.sentinel

cat(paste0("There are ", nrow(dat.sentinel.pwas), " significant PWAS loci.\n"))

pred_prot <- suppressMessages(read_tsv(opt$imputed_P)) # load imputed cis-regulated protein levels for reference individuals

cat(paste0("Starting to perform conditional analysis --\n"))

for (tissue in tissue_list){
    
    pred_ge <- suppressMessages(read_tsv(paste0(opt$imputed_T,"/", tissue,".txt"))) # load imputed cis-regulated gene expression levels for reference individuals
    
    dat.twas <- suppressMessages(read_tsv(paste0(opt$TWAS,"/",tissue,".out"))) # load TWAS result table
    
    dat.twas <- dat.twas[!(is.na(dat.twas$P0)),]
    dat.twas <- dat.twas[!(is.na(dat.twas$TWAS.P)),]
    
    
    ## perform conditional anlaysis for each sentinel PWAS gene and its nearby TWAS genes
    PcT.z <- numeric()
    PcT.p <- numeric()
    TcP.z <- numeric()
    TcP.p <- numeric()
    twas.p <- numeric()
    twas.hit <- character()
    dist <- integer()
    corr <- numeric()
    for (i in 1:nrow(dat.sentinel.pwas)) {
        chr <- dat.sentinel.pwas$CHR[i]
        dat.twas.tmp <- dat.twas[dat.twas$CHR == chr, ]
        
        ## extract out nearby TWAS genes
        
        tmp <- dat.twas.tmp[(dat.twas.tmp$P0 < dat.sentinel.pwas$P0[i] + 500000) & (dat.twas.tmp$P0 > dat.sentinel.pwas$P0[i] - 500000),]
        if(nrow(tmp)==0){
            PcT.z[i] <- NA
            TcP.z[i] <- NA
            PcT.p[i] <- NA
            TcP.p[i] <- NA
            twas.p[i] <- NA
            twas.hit[i] <- NA
            dist[i] <- NA
            corr[i] <- NA
        }else{
            dat.twas.tmp <- tmp[which.min(tmp$TWAS.P),]
            dist[i] <- dat.twas.tmp$P0 - dat.sentinel.pwas$P0[i]
            
            twas.p[i] <- dat.twas.tmp$TWAS.P
            twas.hit[i] <- dat.twas.tmp$ID
            
            pred_mat = matrix(nrow=498,ncol=2)
            tmp <- strsplit(dat.sentinel.pwas$FILE[i], "/")[[1]]; tmp <- tmp[length(tmp)]; tmp <- substr(tmp, start=1, stop=nchar(tmp)-9)
            pred_mat[,1] = pred_prot[[tmp]]
            tmp <- strsplit(dat.twas.tmp$FILE[1], "/")[[1]]; tmp <- tmp[length(tmp)]; tmp <- substr(tmp, start=nchar(dat.twas.tmp$PANEL[1])+2, stop=nchar(tmp)-9)
            pred_mat[,2] = pred_ge[[tmp]]
            
            pred_mat = scale( pred_mat )
            
            # compute cis-regulated genetic correlation
            corr_mat = cor(pred_mat)
            corr[i] = corr_mat[1,2]
            
            b.se = sqrt(1 - corr_mat[1,2]^2)
            
            # estimate conditional effect size (P conditional on T, P|T)
            b = dat.sentinel.pwas$PWAS.Z[i] - corr[i] * dat.twas.tmp$TWAS.Z
            PcT.z[i] = b / b.se
            PcT.p[i] = 2*(pnorm( abs( PcT.z[i] ) , lower.tail=F))
            
            # estimate conditional effect size (T conditional on P, T|P)
            b = dat.twas.tmp$TWAS.Z - corr[i] * dat.sentinel.pwas$PWAS.Z[i]
            TcP.z[i] = b / b.se
            TcP.p[i] = 2*(pnorm( abs( TcP.z[i] ) , lower.tail=F))
            
        }
        
    }
    
    save(dat.sentinel.pwas,
         PcT.z, PcT.p,
         TcP.z, TcP.p,
         twas.p, twas.hit,
         dist, corr, 
         file = paste0(opt$out, "/",tissue,".RDat"))

    cat(paste0(tissue, " is completed.\n"))
    
}


