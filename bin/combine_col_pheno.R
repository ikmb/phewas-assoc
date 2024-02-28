library(dplyr)

args = commandArgs(trailingOnly=TRUE)

# Inputs
file1 <- args[1]
file2 <- args[2]
traittype <- args[3]
output_file <- args[4]



# Read the phenofile and covariate file into dataframes
df1 <- read.table(file1, header = T)
df2 <- read.table(file2, header = T, sep="\t")

if(traittype == "binary"){
  print("binary trait")
  df1[,3] <- df1[,3]-1
}else if (traittype == "quantitative"){
  print("quantitave trait")
}else{
  print( "Invalid trait type")
  stop()
}

merged_df <- inner_join(df1, df2, by = c("FID", "IID"))

write.table(merged_df, file = output_file, row.names = FALSE, quote=F, sep=" ")