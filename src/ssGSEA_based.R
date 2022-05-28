##ssGSEA based method
##


TIS_score = function(data,name){
  library(readr)
  IIS_TIS_signature <- read_delim("marker/IIS_TIS_signature.txt", 
                                  delim = "\t", escape_double = FALSE, 
                                  trim_ws = TRUE)
  all = unique(IIS_TIS_signature$`Cell type`)
  IIS_TIS_geneset = list()
  for (i in all) {
    b = IIS_TIS_signature$Symbol[which(IIS_TIS_signature$`Cell type`==i)]
    #assign(i,b)
    #if(length(b)>1){
    IIS_TIS_geneset[[i]] = b
    #}
  }
  
  library("GSVA")
  gsva.es <- gsva(data$TPM, IIS_TIS_geneset, method="ssgsea", verbose=T)
  # TIS_score
  TIS_sinature = c("CD8 T cells","T helper cells","Tcm cells","Tem cells","Th1 cells","Th2 cells","Th17 cells","Treg cells")
  TIS_score = apply(gsva.es,2,function(x){return(sum(x[TIS_sinature]))})
  IIS_score = colSums(gsva.es)
  
  expMarker = merge(data$Samples,data.frame(Sample=names(TIS_score),IIS_score=IIS_score,TIS_score=TIS_score),by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  
  library(ggpubr)
  g = ggboxplot(expMarker,x="Resp_NoResp",y="TIS_score",color = "Resp_NoResp",
                  palette = "npg",add = "jitter") + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "less"))
  ggsave(paste0("Results/Results_wilcox_test/",name,"/TIS_score.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(TIS_score ~ Resp_NoResp, data = expMarker,alternative = "less")
  auc = roc(Resp_NoResp ~ TIS_score, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("Results/Results_AUC/",name,"/TIS_score.png"),width = 6,height = 6, bg="white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
}

IIS_score = function(data,name){
  library(readr)
  IIS_TIS_signature <- read_delim("marker/IIS_TIS_signature.txt", 
                                  delim = "\t", escape_double = FALSE, 
                                  trim_ws = TRUE)
  all = unique(IIS_TIS_signature$`Cell type`)
  IIS_TIS_geneset = list()
  for (i in all) {
    b = IIS_TIS_signature$Symbol[which(IIS_TIS_signature$`Cell type`==i)]
    #assign(i,b)
    #if(length(b)>1){
    IIS_TIS_geneset[[i]] = b
    #}
  }
  
  library("GSVA")
  gsva.es <- gsva(data$TPM, IIS_TIS_geneset, method="ssgsea", verbose=T)
  # TIS_score
  TIS_sinature = c("CD8 T cells","T helper cells","Tcm cells","Tem cells","Th1 cells","Th2 cells","Th17 cells","Treg cells")
  TIS_score = apply(gsva.es,2,function(x){return(sum(x[TIS_sinature]))})
  IIS_score = colSums(gsva.es)
  
  expMarker = merge(data$Samples,data.frame(Sample=names(TIS_score),IIS_score=IIS_score,TIS_score=TIS_score),by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  
  g = ggboxplot(expMarker,x="Resp_NoResp",y="IIS_score",color = "Resp_NoResp",
                  palette = "npg",add = "jitter") + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "less"))
  ggsave(paste0("Results/Results_wilcox_test/",name,"/IIS_score.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(IIS_score ~ Resp_NoResp, data = expMarker,alternative = "less")
  auc = roc(Resp_NoResp ~ IIS_score, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("Results/Results_AUC/",name,"/IIS_score.png"),width = 6,height = 6, bg="white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
  
}



# APM_score
APM_score = function(data,name){
  APM_signature = c("HLA-A","HLA-B","HLA-C","B2M","TAP1","TAP2","TAPBP")
  APM_score = gsva(data$TPM, list(APM_signature), method="ssgsea", verbose=T)
  expMarker = merge(data$Samples,data.frame(Sample=colnames(APM_score),APM_score=as.numeric(APM_score)),by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  
  g = ggboxplot(expMarker,x="Resp_NoResp",y="APM_score",color = "Resp_NoResp",
                  palette = "npg",add = "jitter") + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "less"))
  ggsave(paste0("Results/Results_wilcox_test/",name,"/APM_score.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(APM_score ~ Resp_NoResp, data = expMarker,alternative = "less")
  auc = roc(Resp_NoResp ~ APM_score, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("Results/Results_AUC/",name,"/APM_score.png"),width = 6,height = 6, bg="white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
}




##################
##innate anti-PD-1 resistance (IPRES) signature
# We row-normalized
# the GSVA scores of each gene set in the IPRES signature across the samples
# from the four cohorts
# The IPRES (enrichment) score was defined as the average
# Z score across all gene sets in the IPRES signature

IPRES_score = function(data,name){
  library("qusage")
  IPRES_signatures = read.gmt("marker/IPRES_signatures.gmt")
  gsva.es <- gsva(data$TPM, IPRES_signatures, method="ssgsea", verbose=T)
  gsva.es = scale(t(gsva.es))
  IPRES_score = apply(gsva.es,1,mean)
  expMarker = merge(data$Samples,data.frame(Sample=names(IPRES_score),IPRES_score=IPRES_score),by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  library(ggpubr)
  g = ggboxplot(expMarker,x="Resp_NoResp",y="IPRES_score",color = "Resp_NoResp",
                  palette = "npg",add = "jitter") +  stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "greater"))
  ggsave(paste0("Results/Results_wilcox_test/",name,"/IPRES_score.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(IPRES_score ~ Resp_NoResp, data = expMarker,alternative = "greater")
  auc = roc(Resp_NoResp ~ IPRES_score, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("Results/Results_AUC/",name,"/IPRES_score.png"),width = 6,height = 6, bg="white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
}





#####
##58 C-ECM ssGSEA
C_ECM_score = function(data,name){
  library(readr)
  C_ECM_genes <- read_delim("marker/C_ECM.txt", delim = ";", 
                            escape_double = FALSE, col_names = FALSE, 
                            trim_ws = TRUE)
  geneSet = list()
  geneSet[["C-ECM"]] = C_ECM_genes$X1
  C_ECM_score <- gsva(data$TPM, geneSet, method="ssgsea", verbose=T)
  expMarker = merge(data$Samples,data.frame(Sample=colnames(C_ECM_score),C_ECM_score=as.numeric(C_ECM_score[1,])),by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  library(ggpubr)
  g = ggboxplot(expMarker,x="Resp_NoResp",y="C_ECM_score",color = "Resp_NoResp",
                  palette = "npg",add = "jitter") +  stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "greater"))
  ggsave(paste0("Results/Results_wilcox_test/",name,"/C_ECM_score.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(C_ECM_score ~ Resp_NoResp, data = expMarker,alternative = "greater")
  auc = roc(Resp_NoResp ~ C_ECM_score, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("Results/Results_AUC/",name,"/C_ECM_score.png"),width = 6,height = 6, bg="white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
}


#############
##immune microenvironment score (IMS) from gastric cancer
IMS_score = function(data,name){
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
  expMarker = merge(data$Samples,data.frame(Sample=names(IMS_score),IMS_score=IMS_score),by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  library(ggpubr)
  g = ggboxplot(expMarker,x="Resp_NoResp",y="IMS_score",color = "Resp_NoResp",
                  palette = "npg",add = "jitter") + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "less"))
  ggsave(paste0("Results/Results_wilcox_test/",name,"/IMS_score.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(IMS_score ~ Resp_NoResp, data = expMarker,alternative = "less")
  auc = roc(Resp_NoResp ~ IMS_score, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("Results/Results_AUC/",name,"/IMS_score.png"),width = 6,height = 6, bg="white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
}



#############################
##Super pathway signatures
PASS_PRE = function(data,name){
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
  expMarker = merge(data$Samples,data.frame(Sample=names(PASS_PRE),PASS_PRE=PASS_PRE),by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  library(ggpubr)
  g = ggboxplot(expMarker,x="Resp_NoResp",y="PASS_PRE",color = "Resp_NoResp",
                  palette = "npg",add = "jitter") + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "less"))
  ggsave(paste0("Results/Results_wilcox_test/",name,"/PASS_PRE.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(PASS_PRE ~ Resp_NoResp, data = expMarker,alternative = "less")
  auc = roc(Resp_NoResp ~ PASS_PRE, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("Results/Results_AUC/",name,"/PASS_PRE.png"),width = 6,height = 6, bg="white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
}

PASS_ON = function(data,name){
  load('marker/Pathway_Singatures.Rdata')
  #Prepare signaure to ssGSEA
  prepare_sig <- function(sig){
    leadingEdge.list <- sig$leadingEdge
    names(leadingEdge.list) <- sig$pathway
    return(leadingEdge.list)
  }
  PASS.ON.Sigs <- prepare_sig(Pathway.Sigs$PASS_ON)
  ssgsea <- gsva(expr = as.matrix(data$TPM), gset.idx.list = PASS.ON.Sigs, method='ssgsea')
  load('marker/PASS_ON_Coefficient.Rdata')
  PASS_ON =  apply(ssgsea[coeft.df$Pathway[-1],],2,function(x){return(sum(x*as.numeric(coeft.df$Weight[-1])))})
  expMarker = merge(data$Samples,data.frame(Sample=names(PASS_ON),PASS_ON=PASS_ON),by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  library(ggpubr)
  g = ggboxplot(expMarker,x="Resp_NoResp",y="PASS_ON",color = "Resp_NoResp",
                  palette = "npg",add = "jitter") +  stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "less"))
  ggsave(paste0("Results/Results_wilcox_test/",name,"/PASS_ON.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(PASS_ON ~ Resp_NoResp, data = expMarker,alternative = "less")
  auc = roc(Resp_NoResp ~ PASS_ON, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("Results/Results_AUC/",name,"/PASS_ON.png"),width = 6,height = 6, bg="white")
  result = c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
}


###MIAS
MIAS_score = function(data,name){
  library(GSVA)
  library(readr)
  MIAS_signatures <- list()
  a = read_csv("marker/MIAS.txt")
  MIAS_signatures[["MIAS"]] <- a$x
  MIAS_score = gsva(data$TPM, MIAS_signatures, method="ssgsea", verbose=T)[1,]
  expMarker = merge(data$Samples,data.frame(Sample=names(MIAS_score),MIAS_score=as.numeric(MIAS_score)),by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  library(ggpubr)
  g = ggboxplot(expMarker,x="Resp_NoResp",y="MIAS_score",color = "Resp_NoResp",
                  palette = "npg",add = "jitter") +  stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "less"))
  ggsave(paste0("Results/Results_wilcox_test/",name,"/MIAS_score.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(MIAS_score ~ Resp_NoResp, data = expMarker,alternative = "less")
  auc = roc(Resp_NoResp ~ MIAS_score, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("Results/Results_AUC/",name,"/MIAS_score.png"),width = 6,height = 6, bg="white")
  result = c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
}



