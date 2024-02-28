library(dplyr)

args = commandArgs(trailingOnly=TRUE)
#FAM, COVARS, COVARS_COLS, PCS, PC_DIMS, WITHCOVARS, !{params.collection_name}
#args <- c("antiTNF_QCed_VCF_noATCG.fam", "antiTNF.covars", "antiTNF.covar_cols", "pcs.txt", 0, TRUE, "antiTNF")

fam <- read.table(args[1])
colnames(fam) <- c("FID","IID","PAT","MAT","SEX","PHENO")
covars <- read.table(args[2], header=T)
covar_cols <- args[3] %>% strsplit(.,split=",") %>% unlist()

number_pcs <- args[5]
withcovars <- args[6] %>% as.logical()
collection_name <- args[7]

if(number_pcs != 0){
  pcs <- read.table(args[4], header=T)
  pc_cols <- paste0("PC", seq(1:number_pcs))
  
  fam_pcs <- merge(fam, pcs, by=c("FID", "IID")) %>%
              select(FID, IID, all_of(pc_cols))
  
} else{
  pc_cols <- ""
  fam_pcs <- fam %>%
    select(FID, IID)
}

if(withcovars){
  fam_pcs_covars <- merge(fam_pcs, covars[,c("FID","IID",covar_cols)], by=c("FID", "IID"))
}else{
  fam_pcs_covars <- fam_pcs
}


#covar_cols; a text line with all used column names comma separated
write(colnames(fam_pcs_covars)[c(-1,-2)], 
      sep = ",", 
      file = paste0(collection_name, ".covar_cols"), 
      ncolumns=length(colnames(fam_pcs_covars)[c(-1,-2)])
      )
#covars; a tab separated table with FID, IID and all used covariate columns, in the same order as the input fam file.
write.table(fam_pcs_covars, 
            sep="\t", 
            file = paste0(collection_name, ".covars"), quote = F, row.names = F)
