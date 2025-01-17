dataset_landscape = function(data,name){
  ##only save the samples with both expression data and clinical information
  expMarker = data.frame(Sample=intersect(colnames(data$TPM),data$Clinical$Sample))
  ##PD-L1
  expression = data$TPM["CD274",]
  expMarker = merge(expMarker,data.frame(Sample=names(expression),PD_L1=expression),by="Sample",all.x = T)
  
  ##PD-1
  expression = data$TPM["PDCD1",]
  expMarker = merge(expMarker,data.frame(Sample=names(expression),PD_1=expression),by="Sample",all.x = T)
  
  ##PD-L2
  expression = data$TPM["PDCD1LG2",]
  expMarker = merge(expMarker,data.frame(Sample=names(expression),PD_L2=expression),by="Sample",all.x = T)
  
  ##CX3CL1
  expression = data$TPM["CX3CL1",]
  expMarker = merge(expMarker,data.frame(Sample=names(expression),CX3CL1=expression),by="Sample",all.x = T)
  
  ##CTLA4 
  expression = data$TPM["CTLA4",]
  expMarker = merge(expMarker,data.frame(Sample=names(expression),CTLA4=expression),by="Sample",all.x = T)
  
  ###cytolytic activity (CYT) scores
  CYT_gene = c("GZMA","PRF1")
  CYT_gene = intersect(CYT_gene,rownames(data$TPM))
  expression = data$TPM[CYT_gene,]
  mean_expression = apply(expression, 2, mean)
  expMarker = merge(expMarker,data.frame(Sample=names(mean_expression),CYT_score=as.numeric(mean_expression)),by="Sample",all.x = T)
  
  ##HLA_DRA 
  expression = data$TPM["HLA-DRA",]
  expMarker = merge(expMarker,data.frame(Sample=names(expression),HLA_DRA=as.numeric(expression)),by="Sample",all.x = T)
  
  ##INF-gamma signature
  IFN_gamma = c("IDO1","CXCL10","CXCL9","HLA-DRA","STAT1","IFNG")
  expression = data$TPM[intersect(IFN_gamma,rownames(data$TPM)),]
  expression = as.data.frame(t(expression))
  expression["IFN_gamma"] = apply(expression, 1, function(x){return(mean(as.numeric(x)))})
  expression["Sample"] = rownames(expression)
  expMarker = merge(expMarker,expression[c("Sample","IFN_gamma")],by="Sample",all.x = T)
  
  ####Expanded_immune_gene_signature
  Expanded_immune_gene_signature = c("CD3D","IDO1","CIITA","CD3E","CCL5","GZMK","CD2","HLA-DRA","CXCL13","IL2RG","NKG7","HLA-E","CXCR6","LAG3","TAGAP","CXCL10","STAT1","GZMB")
  expression = data$TPM[intersect(Expanded_immune_gene_signature,rownames(data$TPM)),]
  expression = as.data.frame(t(expression))
  expression["Expanded_immune_gene_signature"] = apply(expression, 1, function(x){return(mean(as.numeric(x)))})
  expression["Sample"] = rownames(expression)
  expMarker = merge(expMarker,expression[c("Sample","Expanded_immune_gene_signature")],by="Sample",all.x = T)
  
  ####T_cell_inflamed_GEP_score
  T_cell_inflamed_GEP = c("TIGIT","PDCD1LG2","CD27","CD8A","LAG3","CD274","CXCR6","CMKLR1","NKG7","CCL5","PSMB10","IDO1","CXCL9","HLA-DQA1","CD276","STAT1","HLA-DRB1","HLA-E")
  expression = data$TPM[intersect(T_cell_inflamed_GEP,rownames(data$TPM)),]
  expression = as.data.frame(t(expression))
  expression["T_cell_inflamed_GEP_score"] = apply(expression, 1, function(x){return(mean(as.numeric(x)))})
  expression["Sample"] = rownames(expression)
  expMarker = merge(expMarker,expression[c("Sample","T_cell_inflamed_GEP_score")],by="Sample",all.x = T)
  
  ##Immunophenoscore (IPS)
  source("marker/Immunophenogram/IPS.R")
  DF = calculateIPS(data$TPM)
  expMarker = merge(expMarker,data.frame(Sample=DF$SAMPLE,Immunophenoscore=DF$IPS),by="Sample",all.x = T)
  
  ##immuno-predictive score (IMPRES)
  Gene1 = c("PDCD1","CD27","CTLA4","CD40","CD86","CD28","CD80","CD274","CD86","CD40","CD86","CD40","CD28","CD40","TNFRSF14")
  Gene2 = c("TNFSF4","PDCD1","TNFSF4","CD28","TNFSF4","CD86","TNFSF9","VSIR","HAVCR2","PDCD1","CD200","CD80","CD276","CD274","CD86")
  IMPRES_pairs = data.frame(Gene1,Gene2)
  FS = apply(IMPRES_pairs, 1, function(x){
    if(x[1] %in% rownames(data$TPM) & x[2] %in% rownames(data$TPM)){
      expression1 = data$TPM[x[1],]
      expression2 = data$TPM[x[2],]
      return(as.numeric(expression1<expression2))
    }else{
      return(rep(0,ncol(data$TPM)))
    }
  })
  rownames(FS) = colnames(data$TPM)
  IMPRES_score = apply(FS,1,sum)
  expMarker = merge(expMarker,data.frame(Sample=names(IMPRES_score),IMPRES_score=as.numeric(IMPRES_score)),by="Sample",all.x = T)

  ##CRMA
  MGAEA_genes = c("MAGEA3", "CSAG3", "CSAG2","MAGEA2", "MAGEA2B", "CSAG1", "MAGEA12", "MAGEA6")
  intersect_genes = intersect(MGAEA_genes,rownames(data$TPM))
  expression = data$TPM[intersect_genes,]
  expression = as.data.frame(t(expression))
  expression["CRMA_score"] = apply(expression, 1, function(x){return(mean(as.numeric(x)))})
  expression["Sample"] = rownames(expression)
  expMarker = merge(expMarker,expression[c("Sample","CRMA_score")],by="Sample",all.x = T)
  
  ####The immune resistance program
  discretize<-function(v,n.cat){
    q1<-quantile(v,seq(from = (1/n.cat),to = 1,by = (1/n.cat)))
    u<-matrix(nrow = length(v))
    for(i in 2:n.cat){
      u[(v>=q1[i-1])&(v<q1[i])]<-i
    }
    return(u)
  }
  get.semi.random.OE <- function(r,genes.dist.q,b.sign,num.rounds = 1000,full.flag = F){
    # Previous name: get.random.sig.scores
    sign.q<-as.matrix(table(genes.dist.q[b.sign]))
    q<-rownames(sign.q)
    idx.all<-c()
    B<-matrix(data = F,nrow = length(genes.dist.q),ncol = num.rounds)
    Q<-matrix(data = 0,nrow = length(genes.dist.q),ncol = num.rounds)
    for (i in 1:nrow(sign.q)){
      num.genes<-sign.q[i]
      if(num.genes>0){
        idx<-which(is.element(genes.dist.q,q[i]))
        for (j in 1:num.rounds){
          idxj<-sample(idx,num.genes) 
          Q[i,j]<-sum(B[idxj,j]==T)
          B[idxj,j]<-T
        }  
      }
    }
    rand.scores<-apply(B,2,function(x) colMeans(r$zscores[x,]))
    if(full.flag){return(rand.scores)}
    rand.scores<-rowMeans(rand.scores)
    return(rand.scores)
  }
  get.OE.bulk <- function(r,gene.sign = NULL,num.rounds = 1000,full.flag = F){
    set.seed(1234)
    r$genes.mean<-rowMeans(r$tpm)
    r$zscores<-sweep(r$tpm,1,r$genes.mean,FUN = '-')
    r$genes.dist<-r$genes.mean
    r$genes.dist.q<-discretize(r$genes.dist,n.cat = 50)
    r$sig.scores<-matrix(data = 0,nrow = ncol(r$tpm),ncol = length(gene.sign))
    sig.names<-names(gene.sign)
    colnames(r$sig.scores)<-sig.names
    r$sig.scores.raw<-r$sig.scores
    rand.flag<-is.null(r$rand.scores)|!all(is.element(names(gene.sign),colnames(r$rand.scores)))
    if(rand.flag){
      print("Computing also random scores.")
      r$rand.scores<-r$sig.scores
    }
    for (i in sig.names){
      b.sign<-is.element(r$genes,gene.sign[[i]])
      if(sum(b.sign)<2){next()}
      if(rand.flag){
        rand.scores<-get.semi.random.OE(r,r$genes.dist.q,b.sign,num.rounds = num.rounds)
      }else{
        rand.scores<-r$rand.scores[,i]
      }
      raw.scores<-colMeans(r$zscores[b.sign,])
      final.scores<-raw.scores-rand.scores
      r$sig.scores[,i]<-final.scores
      r$sig.scores.raw[,i]<-raw.scores
      r$rand.scores[,i]<-rand.scores
    }
    if(full.flag){return(r)}
    sig.scores<-r$sig.scores
    return(sig.scores)
  }
  data$tpm = data$TPM
  data$genes = rownames(data$TPM)
  load("marker/resistance.program.RData")
  resistance.scores = get.OE.bulk(data, gene.sign = res.sig)
  ims = resistance.scores[,"resu.up"] - resistance.scores[,"resu.down"]
  #overall expression
  OE = data.frame(Sample=colnames(data$TPM),The_immune_resistance_program=scale(as.numeric(ims)))
  expMarker = merge(expMarker,OE,by="Sample",all.x = T)
  
  ##### EMT_Stroma_core_signature
  EMT_Stroma_core_signature = c("FLNA","EMP3","CALD1","FN1","FOXC2","LOX","FBN1","TNC")
  expression = data$TPM[intersect(EMT_Stroma_core_signature,rownames(data$TPM)),]
  expression = as.data.frame(t(expression))
  expression["EMT_Stroma_core_signature"] = apply(expression, 1, function(x){return(mean(as.numeric(x)))})
  expression["Sample"] = rownames(expression)
  expMarker = merge(expMarker,expression[c("Sample","EMT_Stroma_core_signature")],by="Sample",all.x = T)
  
  #####pan-fibroblast TGFβ response signature (F-TBRS)
  F_TBRS_genes = c("ACTA2", "ACTG2", "ADAM12", "ADAM19", "CNN1", "COL4A1", "CCN2", "CTPS1",
                   "RFLNB", "FSTL3", "HSPB1", "IGFBP3", "PXDC1", "SEMA7A", "SH3PXD2A", "TAGLN", 
                   "TGFBI", "TNS1", "TPM1")
  m = data$TPM
  m <- t(scale( t( m ),
                center=TRUE, 
                scale=TRUE)
  )
  m2 = m[intersect(F_TBRS_genes,rownames(data$TPM)),]
  ##' Calculate score across genes and samples
  gsScore <- function(gm, summarizationFunction="PC") {
    if (summarizationFunction == "PC") {
      pc <- prcomp(t(gm),
                   retx=TRUE)
      gss <- pc$x[,1] * sign(cor(pc$x[,1], colMeans(gm)))
    } else {
      gss <- colMeans(gm)
    }
    return(gss)
  }
  F_TBRS_score = gsScore(m2) 
  expression = data.frame(Sample=names(F_TBRS_score),F_TBRS=as.numeric(F_TBRS_score))
  expMarker = merge(expMarker,expression,by="Sample",all.x = T)
  
  ##TMEscore
  library(readxl)
  TME_signature <- read_delim("marker/TMEScore.txt", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)
  TME_gene_immune = TME_signature$Symbol[which(TME_signature$`TME-signature-group`=="TME-gene-A")]
  TME_gene_stroma = TME_signature$Symbol[which(TME_signature$`TME-signature-group`=="TME-gene-B")]
  m = data$TPM
  m <- t(scale( t( m ),
                center=TRUE, 
                scale=TRUE)
  )
  m1 = m[intersect(TME_gene_immune,rownames(m)),]
  m2 = m[intersect(TME_gene_stroma,rownames(m)),]
  ##filter gene with 0 expression
  m1 = m1[which(!is.nan(rowSums(m1))),]
  m2 = m2[which(!is.nan(rowSums(m2))),]
  pc1 = prcomp(t(m1),retx=TRUE)
  # pc1$x[,1]
  pc2 = prcomp(t(m2),retx=TRUE)
  # pc2$x[,1]
  TMEscore = pc1$x[,1] - pc2$x[,1]
  expression = data.frame(Sample=names(TMEscore),TMEscore=as.numeric(TMEscore))
  expMarker = merge(expMarker,expression,by="Sample",all.x = T)
  
  ##risk score
  IRGs = c("LEPR","PRLHR","NR2F2","PRL","NRP1","ANGPTL5","IGF1","TNFRSF10B","TNFRSF10A","PLAU","IFI30")
  coeff = c(0.32196,-0.64921,-0.32677,0.23573,0.39005,0.38166,-0.03522,0.02975,0.39830,0.14607,-0.68625)
  ##only compute the gene expression in the dataset
  expression = data$TPM[IRGs[which(is.element(IRGs,rownames(data$TPM)))],]
  riskScore = apply(expression, 2, function(x){return(sum(x*coeff[which(is.element(IRGs,rownames(data$TPM)))]))})
  expMarker = merge(expMarker,data.frame(Sample=names(riskScore),RiskScore=riskScore),by="Sample",all.x = T)
  
  ##Tertiary lymphoid structures (TLS) 
  TLS_gene = c("CD79B","CD1D","CCR6","LAT","SKAP1","CETP","EIF1AY","RBP5","PTGDS")
  TLS_gene = intersect(TLS_gene,rownames(data$TPM))
  expression = data$TPM[TLS_gene,]
  mean_expression = apply(expression, 2, mean)
  expMarker = merge(expMarker,data.frame(Sample=names(mean_expression),TLS_score=as.numeric(mean_expression)),by="Sample",all.x = T)
  
  ##CXCL9
  expression = data$TPM["CXCL9",]
  expMarker = merge(expMarker,data.frame(Sample=names(expression),CXCL9=as.numeric(expression)),by="Sample",all.x = T)
  
  ##MPS 
  library(readr)
  MPS_gene_list <- read_delim("marker/MPS_gene_list.txt", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)
  expression = data$TPM[intersect(MPS_gene_list$`Gene Symbol`,rownames(data$TPM)),]
  MPS_postive_gene = MPS_gene_list$`Gene Symbol`[which(MPS_gene_list$`Sign in the signature`==1)]
  MPS_negative_gene = MPS_gene_list$`Gene Symbol`[which(MPS_gene_list$`Sign in the signature`==-1)]
  MPS_score = apply(expression, 2, function(x){return(sum(x[which(is.element(names(x),MPS_postive_gene))]) - sum(x[which(is.element(names(x),MPS_negative_gene))]))})
  expMarker = merge(expMarker,data.frame(Sample=names(MPS_score),MPS_score=as.numeric(MPS_score)),by="Sample",all.x = T)
  
  ##Renal 101 Immuno signature
  Renal_101_Immuno_signature = c("CD3G","CD3E","CD8B","THEMIS","TRAT1","GRAP2","CD247",
                                 "CD2","CD96","PRF1","CD6","IL7R","ITK","GPR18","EOMES",
                                 "SIT1","NLRC3","CD244","KLRD1","SH2D1A","CCL5","XCL2",
                                 "CST7","GFI1","KCNA3","PSTPIP1")
  Renal_101_Immuno_signature = intersect(Renal_101_Immuno_signature,rownames(data$TPM))
  expression = data$TPM[Renal_101_Immuno_signature,]
  mean_expression = apply(expression, 2, mean)
  expMarker = merge(expMarker,data.frame(Sample=names(mean_expression),Renal_101_Immuno_signature=as.numeric(mean_expression)),by="Sample",all.x = T)
  
  ##HRH1
  expression = data$TPM["HRH1",]
  expMarker = merge(expMarker,data.frame(Sample=names(expression),HRH1=as.numeric(expression)),by="Sample",all.x = T)
  
  ##TIDE
  library(readr)
  file_path = paste0("TIDE_output/",name,".csv")
  res <- read_csv(file_path)
  expMarker = merge(expMarker,data.frame(Sample=res$Patient,TIDE=res$TIDE),by="Sample",all.x = T)
  
  ##TIS_score
  library(readr)
  IIS_TIS_signature <- read_delim("marker/IIS_TIS_signature.txt", 
                                  delim = "\t", escape_double = FALSE, 
                                  trim_ws = TRUE)
  all = unique(IIS_TIS_signature$`Cell type`)
  IIS_TIS_geneset = list()
  for (i in all) {
    b = IIS_TIS_signature$Symbol[which(IIS_TIS_signature$`Cell type`==i)]
    IIS_TIS_geneset[[i]] = b
  }
  
  library("GSVA")
  gsva.es <- gsva(data$TPM, IIS_TIS_geneset, method="ssgsea", verbose=T)
  # TIS_score
  TIS_sinature = c("CD8 T cells","T helper cells","Tcm cells","Tem cells","Th1 cells","Th2 cells","Th17 cells","Treg cells")
  TIS_score = apply(gsva.es,2,function(x){return(sum(x[TIS_sinature]))})
  IIS_score = colSums(gsva.es)
  expMarker = merge(expMarker,data.frame(Sample=names(TIS_score),IIS_score=IIS_score,TIS_score=TIS_score),by="Sample",all.x = T)
  
  ## APM_score
  APM_signature = c("HLA-A","HLA-B","HLA-C","B2M","TAP1","TAP2","TAPBP")
  APM_score = gsva(data$TPM, list(APM_signature), method="ssgsea", verbose=T)
  expMarker = merge(expMarker,data.frame(Sample=colnames(APM_score),APM_score=as.numeric(APM_score)),by="Sample",all.x = T)
  
  ####innate anti-PD-1 resistance (IPRES)
  library("qusage")
  IPRES_signatures = read.gmt("marker/IPRES_signatures.gmt")
  gsva.es <- gsva(data$TPM, IPRES_signatures, method="ssgsea", verbose=T)
  gsva.es = scale(t(gsva.es))
  IPRES_score = apply(gsva.es,1,mean)
  expMarker = merge(expMarker, data.frame(Sample=names(IPRES_score),IPRES_score=IPRES_score),by="Sample",all.x = T)
  
  
  ##58 C-ECM ssGSEA
  library(readr)
  C_ECM_genes <- read_delim("marker/C_ECM.txt", delim = ";", 
                            escape_double = FALSE, col_names = FALSE, 
                            trim_ws = TRUE)
  geneSet = list()
  geneSet[["C-ECM"]] = C_ECM_genes$X1
  C_ECM_score <- gsva(data$TPM, geneSet, method="ssgsea", verbose=T)
  expMarker = merge(expMarker,data.frame(Sample=colnames(C_ECM_score),C_ECM_score=as.numeric(C_ECM_score[1,])),by="Sample",all.x = T)
  
  ##immune microenvironment score (IMS)
  IMS_signature <- read_delim("marker/IMS_signature.txt", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)
  IMS_signature_geneset = list()
  all = unique(IMS_signature$`Immune cell type`)
  for (i in all) {
    b = IMS_signature$Gene[which(IMS_signature$`Immune cell type` == i)]
    IMS_signature_geneset[[i]] = b
  }
  IMS_signature_meta <- read_delim("marker/IMS_signature_meta.txt", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)
  IMS_signature_HR_more = IMS_signature_meta$`Immune Cell`[which(IMS_signature_meta$HR>1)]
  IMS_signature_HR_less = IMS_signature_meta$`Immune Cell`[which(IMS_signature_meta$HR<1)]
  
  gsva.es <- gsva(data$TPM, IMS_signature_geneset, method="ssgsea", verbose=T)
  
  IMS_score <- apply(gsva.es, 2, function(x){
    NES1 = sum(x[which(is.element(names(x),IMS_signature_HR_more))])
    NES2 = sum(x[which(is.element(names(x),IMS_signature_HR_less))])
    return(NES2-NES1)
  })
  expMarker = merge(expMarker,data.frame(Sample=names(IMS_score),IMS_score=IMS_score),by="Sample",all.x = T)
  
  ##Super pathway signatures
  load('marker/Pathway_Singatures.Rdata')
  #Prepare signaure to ssGSEA
  prepare_sig <- function(sig){
    leadingEdge.list <- sig$leadingEdge
    names(leadingEdge.list) <- sig$pathway
    return(leadingEdge.list)
  }
  PASS.PRE.Sigs <- prepare_sig(Pathway.Sigs$PASS_PRE)
  ssgsea <- gsva(expr = as.matrix(data$TPM), gset.idx.list = PASS.PRE.Sigs, method='ssgsea')
  load('marker/PASS_PRE_Coefficient.Rdata')
  PASS_PRE =  apply(ssgsea[coeft.df$Pathway[-1],],2,function(x){return(sum(x*as.numeric(coeft.df$Weight[-1])))})
  expMarker = merge(expMarker,data.frame(Sample=names(PASS_PRE),PASS_PRE=PASS_PRE),by="Sample",all.x = T)
  
  PASS.ON.Sigs <- prepare_sig(Pathway.Sigs$PASS_ON)
  ssgsea <- gsva(expr = as.matrix(data$TPM), gset.idx.list = PASS.ON.Sigs, method='ssgsea')
  load('marker/PASS_ON_Coefficient.Rdata')
  PASS_ON =  apply(ssgsea[coeft.df$Pathway[-1],],2,function(x){return(sum(x*as.numeric(coeft.df$Weight[-1])))})
  expMarker = merge(expMarker,data.frame(Sample=names(PASS_ON),PASS_ON=PASS_ON),by="Sample",all.x = T)
  
  ###MIAS
  library(GSVA)
  library(readr)
  MIAS_signatures <- list()
  a = read_csv("marker/MIAS.txt")
  MIAS_signatures[["MIAS"]] <- a$x
  MIAS_score = gsva(data$TPM, MIAS_signatures, method="ssgsea", verbose=T)[1,]
  expMarker = merge(expMarker,data.frame(Sample=names(MIAS_score),MIAS_score=as.numeric(MIAS_score)),by="Sample",all.x = T)
  
  # CD8 T cells 
  library(xCell)
  xCell = xCellAnalysis(data$TPM)
  expMarker = merge(expMarker,data.frame(Sample=names(xCell["CD8+ T-cells",]),CD8T_xCell=as.numeric(xCell["CD8+ T-cells",])),by="Sample",all.x = T)
  
  library(MCPcounter)
  MCPcounter_score = MCPcounter.estimate(data$TPM,featuresType="HUGO_symbols")
  expMarker = merge(expMarker,data.frame(Sample=names(MCPcounter_score["CD8 T cells",]),CD8T_MCPcounter=as.numeric(MCPcounter_score["CD8 T cells",])),by="Sample",all.x = T)
  
  file_path = paste0("CIBERSORTx_output/CIBERSORTx_",name,"_Results.txt")
  library(readr)
  CIBERSORTx_result <- read_delim(file_path,
                                  delim = "\t",
                                  escape_double = FALSE, 
                                  trim_ws = TRUE)
  expMarker = merge(expMarker,data.frame(Sample=CIBERSORTx_result$Mixture,CD8T_CIBERSORTx=as.numeric(CIBERSORTx_result$`T cells CD8`)),by="Sample",all.x = T)
  
  ##immunoscore construct based on CIBERSORTx_output
  file_path = paste0("CIBERSORTx_output/CIBERSORTx_",name,"_Results.txt")
  library(readr)
  CIBERSORTx_result <- read_delim(file_path,
                                  delim = "\t",
                                  escape_double = FALSE, 
                                  trim_ws = TRUE)
  
  immunoscore = 1.13*CIBERSORTx_result$`B cells naive` + 1.36*CIBERSORTx_result$`B cells memory` + 5.92*CIBERSORTx_result$Eosinophils + 
    9.70*CIBERSORTx_result$`T cells follicular helper` + 15.34*CIBERSORTx_result$`T cells regulatory (Tregs)` -
    1.14*CIBERSORTx_result$`Macrophages M0` - 2.31*CIBERSORTx_result$`Plasma cells` - 4.52*CIBERSORTx_result$`T cells gamma delta`
  expMarker = merge(expMarker,data.frame(Sample=CIBERSORTx_result$Mixture,Immunoscore_CIBERSORTx=immunoscore),by="Sample",all.x = T)
  
  ##EcoTyper: CE9
  file_path = paste0("EcoTyper_output/",name,"_ecotyper_output/Carcinoma_Ecotypes/Ecotype_Assignment.txt")
  library(readr)
  Ecotype_Assignment <- read_delim(file_path, 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)
  expMarker = merge(expMarker,data.frame(Sample=Ecotype_Assignment$ID,Ecotype=Ecotype_Assignment$`Carcinoma Ecotype`),by="Sample",all.x = T)
  
  
  file_path = paste0("MFP_output/",name,".txt")
  library(readr)
  MFP_result <- read_delim(file_path,
                           delim = "\t", escape_double = FALSE, 
                           col_names = FALSE, trim_ws = TRUE)
  
  expMarker = merge(expMarker,data.frame(Sample=MFP_result$X1,MFP=MFP_result$X2),by="Sample",all.x = T)

  
  ##averager2ssGSEA
  ##IFN_gamma_marker_ssGSEA
  library(GSVA)
  IFN_gamma = c("ID01","CXCL10","CXCL9","HLA-DRA","STAT1","INFG")
  gene_set = list()
  gene_set[["IFN_gamma"]] = IFN_gamma
  gsva.es <- gsva(data$TPM, gene_set, method="ssgsea", verbose=T)
  expMarker = merge(expMarker,data.frame(Sample=colnames(gsva.es),IFN_gamma_ssGSEA=as.numeric(gsva.es)),by="Sample")
  
  ##Expanded_immune_gene_signature_ssGSEA
  Expanded_immune_gene_signature = c("CD3D","IDO1","CIITA","CD3E","CCL5","GZMK","CD2","HLA-DRA","CXCL13","IL2RG","NKG7","HLA-E","CXCR6","LAG3","TAGAP","CXCL10","STAT1","GZMB")
  gene_set = list()
  gene_set[["Expanded_immune_gene_signature"]] = Expanded_immune_gene_signature
  gsva.es <- gsva(data$TPM, gene_set, method="ssgsea", verbose=T)
  expMarker = merge(expMarker,data.frame(Sample=colnames(gsva.es),Expanded_immune_gene_ssGSEA=as.numeric(gsva.es)),by="Sample")
  
  ##T_cell_inflamed_GEP_ssGSEA
  T_cell_inflamed_GEP = c("TIGIT","PDCD1LG2","CD27","CD8A","LAG3","CD274","CXCR6","CMKLR1","NKG7","CCL5","PSMB10","IDO1","CXCL9","HLA-DQA1","CD276","STAT1","HLA-DRB1","HLA-E")
  gene_set = list()
  gene_set[["T_cell_inflamed_GEP"]] = T_cell_inflamed_GEP
  gsva.es <- gsva(data$TPM, gene_set, method="ssgsea", verbose=T)
  expMarker = merge(expMarker,data.frame(Sample=colnames(gsva.es),T_cell_inflamed_GEP_ssGSEA=as.numeric(gsva.es)),by="Sample")
  
  
  ##CRMA_ssGSEA
  MGAEA_genes = c("MAGEA3", "CSAG3", "CSAG2","MAGEA2", "MAGEA2B", "CSAG1", "MAGEA12", "MAGEA6")
  gene_set = list()
  gene_set[["CRMA"]] = MGAEA_genes
  gsva.es <- gsva(data$TPM, gene_set, method="ssgsea", verbose=T)
  expMarker = merge(expMarker,data.frame(Sample=colnames(gsva.es),CRMA_ssGSEA=as.numeric(gsva.es)),by="Sample")
  
  ##EMT_Stroma_core_ssGSEA
  EMT_Stroma_core_signature = c("FLNA","EMP3","CALD1","FN1","FOXC2","LOX","FBN1","TNC")
  gene_set = list()
  gene_set[["EMT_Stroma_core_signature"]] = EMT_Stroma_core_signature
  gsva.es <- gsva(data$TPM, gene_set, method="ssgsea", verbose=T)
  expMarker = merge(expMarker,data.frame(Sample=colnames(gsva.es),EMT_Stroma_core_ssGSEA=as.numeric(gsva.es)),by="Sample")
  
  ##F_TBRS_ssGSEA
  F_TBRS_genes = c("ACTA2", "ACTG2", "ADAM12", "ADAM19", "CNN1", "COL4A1", "CTGF", "CTPS1",
                   "FAM101B", "FSTL3", "HSPB1", "IGFBP3", "PXDC1", "SEMA7A", "SH3PXD2A", "TAGLN", 
                   "TGFBI", "TNS1", "TPM1")
  
  gene_set = list()
  gene_set[["F_TBRS_genes"]] = F_TBRS_genes
  gsva.es <- gsva(data$TPM, gene_set, method="ssgsea", verbose=T)
  expMarker = merge(expMarker,data.frame(Sample=colnames(gsva.es),F_TBRS_ssGSEA=as.numeric(gsva.es)),by="Sample")
  
  ##Risk_ssGSEA
  IRGs = c("LEPR","PRLHR","NR2F2","PRL","NRP1","ANGPTL5","IGF1","TNFRSF10B","TNFRSF10A","PLAU","IFI30")
  gene_set = list()
  gene_set[["IRGs"]] = IRGs
  gsva.es <- gsva(data$TPM, gene_set, method="ssgsea", verbose=T)
  expMarker = merge(expMarker,data.frame(Sample=colnames(gsva.es),RiskScore_ssGSEA=as.numeric(gsva.es)),by="Sample")
  
  ##TLS_score_ssGSEA
  TLS_gene = c("CD79B","CD1D","CCR6","LAT","SKAP1","CETP","EIF1AY","RBP5","PTGDS")
  gene_set = list()
  gene_set[["TLS_gene"]] = TLS_gene
  gsva.es <- gsva(data$TPM, gene_set, method="ssgsea", verbose=T)
  expMarker = merge(expMarker,data.frame(Sample=colnames(gsva.es),TLS_score_ssGSEA=as.numeric(gsva.es)),by="Sample")
  
  ##Renal_101_Immuno_ssGSEA
  Renal_101_Immuno_signature = c("CD3G","CD3E","CD8B","THEMIS","TRAT1","GRAP2","CD247",
                                 "CD2","CD96","PRF1","CD6","IL7R","ITK","GPR18","EOMES",
                                 "SIT1","NLRC3","CD244","KLRD1","SH2D1A","CCL5","XCL2",
                                 "CST7","GFI1","KCNA3","PSTPIP1")
  gene_set = list()
  gene_set[["Renal_101_Immuno_signature"]] = Renal_101_Immuno_signature
  gsva.es <- gsva(data$TPM, gene_set, method="ssgsea", verbose=T)
  expMarker = merge(expMarker,data.frame(Sample=colnames(gsva.es),Renal_101_Immuno_ssGSEA=as.numeric(gsva.es)),by="Sample")
  
  
  return(expMarker)
}

load("./ICB_data/Braun et al/Braun_data.Rdata")
data1=Braun_data
name1="Braun_2020"
Landscape = dataset_landscape(data1,name1)
Braun_data[["Landscape"]] = Landscape
save(Braun_data,file = 'ICB_data/Braun et al/Braun_2020.Rdata')

load("./ICB_data/Du et al/MGH_PRE_data.RData")
data2=MGH_PRE_data
name2="MGH_PRE_2021"
Landscape = dataset_landscape(data2,name2)
MGH_PRE_data[["Landscape"]] = Landscape
save(MGH_PRE_data,file = 'ICB_data/Du et al/MGH_PRE_2021.Rdata')

load("./ICB_data/Du et al/MGH_ON_data.RData")
data3=MGH_ON_data
name3="MGH_ON_2021"
Landscape = dataset_landscape(data3,name3)
MGH_ON_data[["Landscape"]] = Landscape
save(MGH_ON_data,file = 'ICB_data/Du et al/MGH_ON_2021.Rdata')

load("./ICB_data/Gide et al/Gide_data.RData")
data4=Gide_data
name4="Gide_2019"
Landscape = dataset_landscape(data4,name4)
Gide_data[["Landscape"]] = Landscape
save(Gide_data,file = 'ICB_data/Gide et al/Gide_2019.Rdata')

load("./ICB_data/Hugo et al/Hugo_data.RData")
data5=Hugo_data
name5="Hugo_2016"
Landscape = dataset_landscape(data5,name5)
Hugo_data[["Landscape"]] = Landscape
save(Hugo_data,file = 'ICB_data/Hugo et al/Hugo_2016.Rdata')

load("./ICB_data/Jung et al/Jung_data.Rdata")
data6=Jung_data
name6="Jung_2019"
Landscape = dataset_landscape(data6,name6)
Jung_data[["Landscape"]] = Landscape
save(Jung_data,file = 'ICB_data/Jung et al/Jung_2019.Rdata')

load("./ICB_data/Kim et al/Kim_data.RData")
data7=Kim_data
name7="Kim_2018"
Landscape = dataset_landscape(data7,name7)
Kim_data[["Landscape"]] = Landscape
save(Kim_data,file = 'ICB_data/Kim et al/Kim_2018.Rdata')

load("./ICB_data/Lee et al/Lee_data.RData")
data8=Lee_data
name8="Lee_2020"
Landscape = dataset_landscape(data8,name8)
Lee_data[["Landscape"]] = Landscape
save(Lee_data,file = 'ICB_data/Lee et al/Lee_2020.Rdata')

load("./ICB_data/Liu et al/Liu_data.Rdata")
data9=Liu_data
name9="Liu_2019"
Landscape = dataset_landscape(data9,name9)
Liu_data[["Landscape"]] = Landscape
save(Liu_data,file = 'ICB_data/Liu et al/Liu_2019.Rdata')

load("./ICB_data/Mariathasan et al/Mariathasan_data.Rdata")
data10=Mariathasan_data
name10="Mariathasan_2018"
Landscape = dataset_landscape(data10,name10)
Mariathasan_data[["Landscape"]] = Landscape
save(Mariathasan_data,file = 'ICB_data/Mariathasan et al/Mariathasan_2018.Rdata')

load("./ICB_data/Miao et al/Miao_data.Rdata")
data11=Miao_data
name11="Miao_2018"
Landscape = dataset_landscape(data11,name11)
Miao_data[["Landscape"]] = Landscape
save(Miao_data,file = 'ICB_data/Miao et al/Miao_2018.Rdata')

load("./ICB_data/Motzer et al/Motzer_data.Rdata")
data12=Motzer_data
name12="Motzer_2020"
Landscape = dataset_landscape(data12,name12)
Motzer_data[["Landscape"]] = Landscape
save(Motzer_data,file = 'ICB_data/Motzer et al/Motzer_2020.Rdata')

load("./ICB_data/Nathanson et al/Nathanson_data.Rdata")
data13=Nathanson_data
name13="Nathanson_2017"
Landscape = dataset_landscape(data13,name13)
Nathanson_data[["Landscape"]] = Landscape
save(Nathanson_data,file = 'ICB_data/Nathanson et al/Nathanson_2017.Rdata')

load("./ICB_data/Riaz et al/Riaz_data.RData")
data14=Riaz_data
name14="Riaz_2017"
Landscape = dataset_landscape(data14,name14)
Riaz_data[["Landscape"]] = Landscape
save(Riaz_data,file = 'ICB_data/Riaz et al/Riaz_2017.Rdata')

load("./ICB_data/Snyder et al/Snyder_data.Rdata")
data15=Snyder_data
name15="Snyder_2017"
Landscape = dataset_landscape(data15,name15)
Snyder_data[["Landscape"]] = Landscape
save(Snyder_data,file = 'ICB_data/Snyder et al/Snyder_2017.Rdata')

load("./ICB_data/VanAllen et al/VanAllen_data.RData")
data16=VanAllen_data
name16="VanAllen_2015"
Landscape = dataset_landscape(data16,name16)
VanAllen_data[["Landscape"]] = Landscape
save(VanAllen_data,file = 'ICB_data/VanAllen et al/VanAllen_2015.Rdata')


###sub-datasets
load("./ICB_data/Gide et al/Gide_PRE_data.Rdata")
data17=Gide_PRE_data
name17="Gide_PRE_2019"
Landscape = dataset_landscape(data17,name17)
Gide_PRE_data[["Landscape"]] = Landscape
save(Gide_PRE_data,file = 'ICB_data/Gide et al/Gide_PRE_2019.Rdata')

load("./ICB_data/Gide et al/Gide_ON_data.Rdata")
data18=Gide_ON_data
name18="Gide_ON_2019"
Landscape = dataset_landscape(data18,name18)
Gide_ON_data[["Landscape"]] = Landscape
save(Gide_ON_data,file = 'ICB_data/Gide et al/Gide_ON_2019.Rdata')

load("./ICB_data/Riaz et al/Riaz_PRE_data.Rdata")
data19=Riaz_PRE_data
name19="Riaz_PRE_2017"
Landscape = dataset_landscape(data19,name19)
Riaz_PRE_data[["Landscape"]] = Landscape
save(Riaz_PRE_data,file = 'ICB_data/Riaz et al/Riaz_PRE_2017.Rdata')

load("./ICB_data/Riaz et al/Riaz_ON_data.Rdata")
data20=Riaz_ON_data
name20="Riaz_ON_2017"
Landscape = dataset_landscape(data20,name20)
Riaz_ON_data[["Landscape"]] = Landscape
save(Riaz_ON_data,file = 'ICB_data/Riaz et al/Riaz_ON_2017.Rdata')

load("./ICB_data/Lee et al/Lee_PRE_data.Rdata")
data21=Lee_PRE_data
name21="Lee_PRE_2020"
Landscape = dataset_landscape(data21,name21)
Lee_PRE_data[["Landscape"]] = Landscape
save(Lee_PRE_data,file = 'ICB_data/Lee et al/Lee_PRE_2020.Rdata')

load("./ICB_data/Lee et al/Lee_ON_data.Rdata")
data22=Lee_ON_data
name22="Lee_ON_2020"
Landscape = dataset_landscape(data22,name22)
Lee_ON_data[["Landscape"]] = Landscape
save(Lee_ON_data,file = 'ICB_data/Lee et al/Lee_ON_2020.Rdata')


load("./ICB_data/Gide et al/Gide_MONO_data.Rdata")
data23=Gide_MONO_data
name23="Gide_MONO_2019"
Landscape = dataset_landscape(data23,name23)
Gide_MONO_data[["Landscape"]] = Landscape
save(Gide_MONO_data,file = 'ICB_data/Gide et al/Gide_MONO_2019.Rdata')

load("./ICB_data/Gide et al/Gide_COMBINE_data.Rdata")
data24=Gide_COMBINE_data
name24="Gide_COMBINE_2019"
Landscape = dataset_landscape(data24,name24)
Gide_COMBINE_data[["Landscape"]] = Landscape
save(Gide_COMBINE_data,file = 'ICB_data/Gide et al/Gide_COMBINE_2019.Rdata')

load("./ICB_data/Riaz et al/Riaz_NAIVE_data.Rdata")
data25=Riaz_NAIVE_data
name25="Riaz_NAIVE_2017"
Landscape = dataset_landscape(data25,name25)
Riaz_NAIVE_data[["Landscape"]] = Landscape
save(Riaz_NAIVE_data,file = 'ICB_data/Riaz et al/Riaz_NAIVE_2017.Rdata')

load("./ICB_data/Riaz et al/Riaz_EXPOSURE_data.Rdata")
data26=Riaz_EXPOSURE_data
name26="Riaz_EXPOSURE_2017"
Landscape = dataset_landscape(data26,name26)
Riaz_EXPOSURE_data[["Landscape"]] = Landscape
save(Riaz_EXPOSURE_data,file = 'ICB_data/Riaz et al/Riaz_EXPOSURE_2017.Rdata')

load("./ICB_data/Liu et al/Liu_NAIVE_data.Rdata")
data27=Liu_NAIVE_data
name27="Liu_NAIVE_2019"
Landscape = dataset_landscape(data27,name27)
Liu_NAIVE_data[["Landscape"]] = Landscape
save(Liu_NAIVE_data,file = 'ICB_data/Liu et al/Liu_NAIVE_2019.Rdata')

load("./ICB_data/Liu et al/Liu_EXPOSURE_data.Rdata")
data28=Liu_EXPOSURE_data
name28="Liu_EXPOSURE_2019"
Landscape = dataset_landscape(data28,name28)
Liu_EXPOSURE_data[["Landscape"]] = Landscape
save(Liu_EXPOSURE_data,file = 'ICB_data/Liu et al/Liu_EXPOSURE_2019.Rdata')


