#utils
# Gene symbol search function
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
gene = select(org.Hs.eg.db, keys=keys(org.Hs.eg.db), columns=c("SYMBOL","GENETYPE"))
pcg = gene$SYMBOL[which(gene$GENETYPE=="protein-coding")]

ENTREZtoSYMBOL <- function(x){
  vec=mapIds(org.Hs.eg.db,keys=x,column="SYMBOL",keytype="ENTREZID",multiVals="first")
  repl=which(is.na(vec))
  for(y in repl){vec[y]=x[y]}
  return(vec)
}

# Gene symbol search function
ENSEMBLIDtoSYMBOL <- function(x){
  vec=mapIds(org.Hs.eg.db,keys=x,column="SYMBOL",keytype="ENSEMBL",multiVals="first")
  repl=which(is.na(vec))
  for(y in repl){vec[y]=x[y]}
  return(vec)
}

ENSEMBLPROTtoSYMBOL <- function(x){
  vec=mapIds(org.Hs.eg.db,keys=x,column="SYMBOL",keytype="ENSEMBLPROT",multiVals="first")
  repl=which(is.na(vec))
  for(y in repl){vec[y]=x[y]}
  return(vec)
}


#########
##fpkmToTpm
fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
##countToTPM
countToTpm <- function(counts, effLen)
{
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}

##TCGA annotation GTF file
library("rtracklayer")
gtf_data = import('gencode.v22.annotation.gtf') #gtf的路径
#这里使用import导入gtf文件， 生成一个GRangs对象
gtf_data = as.data.frame(gtf_data)
gtf_data = gtf_data[which(gtf_data$type == "gene" & gtf_data$gene_type=="protein_coding"),]



