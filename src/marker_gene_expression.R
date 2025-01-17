#marker gene expression can predict ICB response

library(gridExtra)
library(grid)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(pROC)

# scale_tpm = scale(data$TPM)
# write.table(scale_tpm,file = paste0("ICB_data/",name,"_scale_tpm.txt"),quote = F,sep = "\t")


my_comparisons = list(c("No_Response","Response"))

##PD-L1
PD_L1_marker = function(data,name){
  expression = data$TPM["CD274",]
  expMarker = merge(data$Samples,data.frame(Sample=names(expression),PDL1=expression),by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  g = ggboxplot(expMarker,x="Resp_NoResp",y="PDL1",color = "Resp_NoResp",palette = "npg",add = "jitter") + 
    stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "less"))
  ggsave(paste0("Results/Results_wilcox_test/",name,"/PD_L1.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(PDL1 ~ Resp_NoResp, data = expMarker,alternative = "less")
  auc = roc(Resp_NoResp ~ PDL1, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("Results/Results_AUC/",name,"/PD_L1.png"),width = 6,height = 6,bg = "white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
}


##PD-1
PD_1_marker = function(data,name){
  expression = data$TPM["PDCD1",]
  #names(expression)
  expMarker = merge(data$Samples,data.frame(Sample=names(expression),PD1=expression),by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  g = ggboxplot(expMarker,x="Resp_NoResp",y="PD1",color = "Resp_NoResp",palette = "npg",add = "jitter") +
    stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "less"))
  ggsave(paste0("Results/Results_wilcox_test/",name,"/PD_1.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(PD1 ~ Resp_NoResp, data = expMarker,alternative = "less")
  auc = roc(Resp_NoResp ~ PD1, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("Results/Results_AUC/",name,"/PD_1.png"),width = 6,height = 6, bg="white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
}



##PD-L2
PD_L2_marker = function(data,name){
  expression = data$TPM["PDCD1LG2",]
  #names(expression)
  expMarker = merge(data$Samples,data.frame(Sample=names(expression),PDL2=expression),by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  g = ggboxplot(expMarker,x="Resp_NoResp",y="PDL2",color = "Resp_NoResp",palette = "npg",add = "jitter") + 
    stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "less"))
  ggsave(paste0("Results/Results_wilcox_test/",name,"/PD_L2.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(PDL2 ~ Resp_NoResp, data = expMarker,alternative = "less")
  auc = roc(Resp_NoResp ~ PDL2, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("Results/Results_AUC/",name,"/PD_L2.png"),width = 6,height = 6, bg="white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
}


##CX3CL1
CX3CL1_marker = function(data,name){
  expression = data$TPM["CX3CL1",]
  #names(expression)
  expMarker = merge(data$Samples,data.frame(Sample=names(expression),CX3CL1=expression),by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  g = ggboxplot(expMarker,x="Resp_NoResp",y="CX3CL1",color = "Resp_NoResp",palette = "npg",add = "jitter") + 
    stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "less"))
  ggsave(paste0("Results/Results_wilcox_test/",name,"/CX3CL1.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(CX3CL1 ~ Resp_NoResp, data = expMarker,alternative = "less")
  auc = roc(Resp_NoResp ~ CX3CL1, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("Results/Results_AUC/",name,"/CX3CL1.png"),width = 6,height = 6, bg="white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
}


##CTLA4 
CTLA4_marker = function(data,name){
  expression = data$TPM["CTLA4",]
  #names(expression)
  expMarker = merge(data$Samples,data.frame(Sample=names(expression),CTLA4=expression),by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  g = ggboxplot(expMarker,x="Resp_NoResp",y="CTLA4",color = "Resp_NoResp",palette = "npg",add = "jitter") + 
    stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "less"))
  ggsave(paste0("Results/Results_wilcox_test/",name,"/CTLA4.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(CTLA4 ~ Resp_NoResp, data = expMarker,alternative = "less")
  auc = roc(Resp_NoResp ~ CTLA4, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("Results/Results_AUC/",name,"/CTLA4.png"),width = 6,height = 6, bg="white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
}


############
###cytolytic activity (CYT) scores
CYT_marker = function(data,name){
  CYT_gene = c("GZMA","PRF1")
  CYT_gene = intersect(CYT_gene,rownames(data$TPM))
  expression = data$TPM[CYT_gene,]
  mean_expression = apply(expression, 2, mean)
  expMarker = merge(data$Samples,data.frame(Sample=names(mean_expression),CYT_score=as.numeric(mean_expression)),by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  g = ggboxplot(expMarker,x="Resp_NoResp",y='CYT_score',color = "Resp_NoResp",
                palette = "npg",add = "jitter") + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "less"))
  ggsave(paste0("Results/Results_wilcox_test/",name,"/CYT_score.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(CYT_score ~ Resp_NoResp, data = expMarker,alternative = "less")
  auc = roc(Resp_NoResp ~ CYT_score, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("Results/Results_AUC/",name,"/CYT_score.png"),width = 6,height = 6, bg="white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
}

##############
##HLA_DRA 
HLA_DRA_marker = function(data,name){
  expression = data$TPM["HLA-DRA",]
  expMarker = merge(data$Samples,data.frame(Sample=names(expression),HLA_DRA=as.numeric(expression)),by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  g = ggboxplot(expMarker,x="Resp_NoResp",y='HLA_DRA',color = "Resp_NoResp",
                palette = "npg",add = "jitter") + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "less"))
  ggsave(paste0("Results/Results_wilcox_test/",name,"/HLA_DRA.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(HLA_DRA ~ Resp_NoResp, data = expMarker,alternative = "less")
  auc = roc(Resp_NoResp ~ HLA_DRA, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("Results/Results_AUC/",name,"/HLA_DRA.png"),width = 6,height = 6, bg="white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
}



##INF-gamma signature
IFN_gamma_marker = function(data,name){
  IFN_gamma = c("IDO1","CXCL10","CXCL9","HLA-DRA","STAT1","IFNG")
  expression = data$TPM[intersect(IFN_gamma,rownames(data$TPM)),]
  #names(expression)
  expression = as.data.frame(t(expression))
  expression["IFN_gamma"] = apply(expression, 1, function(x){return(mean(as.numeric(x)))})
  expression["Sample"] = rownames(expression)
  expMarker = merge(data$Samples,expression[c("Sample","IFN_gamma")],by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  g = ggboxplot(expMarker,x="Resp_NoResp",y="IFN_gamma",color = "Resp_NoResp",
                palette = "npg",add = "jitter") + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "less"))
  ggsave(paste0("Results/Results_wilcox_test/",name,"/IFN_gamma.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(IFN_gamma ~ Resp_NoResp, data = expMarker,alternative = "less")
  auc = roc(Resp_NoResp ~ IFN_gamma, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("Results/Results_AUC/",name,"/IFN_gamma.png"),width = 6,height = 6, bg="white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
}



####Expanded_immune_gene_signature
Expanded_immune_gene_signature_marker = function(data,name){
  Expanded_immune_gene_signature = c("CD3D","IDO1","CIITA","CD3E","CCL5","GZMK","CD2","HLA-DRA","CXCL13","IL2RG","NKG7","HLA-E","CXCR6","LAG3","TAGAP","CXCL10","STAT1","GZMB")
  expression = data$TPM[intersect(Expanded_immune_gene_signature,rownames(data$TPM)),]
  #names(expression)
  expression = as.data.frame(t(expression))
  expression["Expanded_immune_gene_signature"] = apply(expression, 1, function(x){return(mean(as.numeric(x)))})
  expression["Sample"] = rownames(expression)
  expMarker = merge(data$Samples,expression[c("Sample","Expanded_immune_gene_signature")],by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  g = ggboxplot(expMarker,x="Resp_NoResp",y="Expanded_immune_gene_signature",color = "Resp_NoResp",
                palette = "npg",add = "jitter") + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "less"))
  ggsave(paste0("Results/Results_wilcox_test/",name,"/Expanded_immune_gene_signature.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(Expanded_immune_gene_signature ~ Resp_NoResp, data = expMarker,alternative = "less")
  auc = roc(Resp_NoResp ~ Expanded_immune_gene_signature, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("Results/Results_AUC/",name,"/Expanded_immune_gene_signature.png"),width = 6,height = 6, bg="white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
}

####T_cell_inflamed_GEP_score
T_cell_inflamed_GEP_score = function(data,name){
  T_cell_inflamed_GEP = c("TIGIT","PDCD1LG2","CD27","CD8A","LAG3","CD274","CXCR6","CMKLR1","NKG7","CCL5","PSMB10","IDO1","CXCL9","HLA-DQA1","CD276","STAT1","HLA-DRB1","HLA-E")
  expression = data$TPM[intersect(T_cell_inflamed_GEP,rownames(data$TPM)),]
  #names(expression)
  expression = as.data.frame(t(expression))
  expression["T_cell_inflamed_GEP_score"] = apply(expression, 1, function(x){return(mean(as.numeric(x)))})
  expression["Sample"] = rownames(expression)
  expMarker = merge(data$Samples,expression[c("Sample","T_cell_inflamed_GEP_score")],by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  g = ggboxplot(expMarker,x="Resp_NoResp",y="T_cell_inflamed_GEP_score",color = "Resp_NoResp",
                palette = "npg",add = "jitter") + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "less"))
  ggsave(paste0("Results/Results_wilcox_test/",name,"/T_cell_inflamed_GEP_score.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(T_cell_inflamed_GEP_score ~ Resp_NoResp, data = expMarker,alternative = "less")
  auc = roc(Resp_NoResp ~ T_cell_inflamed_GEP_score, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("Results/Results_AUC/",name,"/T_cell_inflamed_GEP_score.png"),width = 6,height = 6, bg="white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
  
}


##Immunophenoscore (IPS)
Immunophenoscore = function(data,name){
  source("marker/Immunophenogram/IPS.R")
  DF = calculateIPS(data$TPM)
  expMarker = merge(data$Samples,data.frame(Sample=DF$SAMPLE,IPS=DF$IPS),by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  g = ggboxplot(expMarker,x="Resp_NoResp",y="IPS",color = "Resp_NoResp",
                palette = "npg",add = "jitter") + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "less"))
  ggsave(paste0("Results/Results_wilcox_test/",name,"/Immunophenoscore.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(IPS ~ Resp_NoResp, data = expMarker,alternative = "less")
  auc = roc(Resp_NoResp ~ IPS, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("Results/Results_AUC/",name,"/Immunophenoscore.png"),width = 6,height = 6, bg="white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
}




##immuno-predictive score (IMPRES) coding in MATLAB
#15 pairs of gene,0-15, High scores predict good response
IMPRES = function(data,name){
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
  expMarker = merge(data$Samples,data.frame(Sample=names(IMPRES_score),IMPRES_score=as.numeric(IMPRES_score)),by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  g = ggboxplot(expMarker,x="Resp_NoResp",y="IMPRES_score",color = "Resp_NoResp",
                palette = "npg",add = "jitter") + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "less"))
  ggsave(paste0("Results/Results_wilcox_test/",name,"/IMPRES_score.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(IMPRES_score ~ Resp_NoResp, data = expMarker,alternative = "less")
  auc = roc(Resp_NoResp ~ IMPRES_score, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("Results/Results_AUC/",name,"/IMPRES_score.png"),width = 6,height = 6, bg="white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
}


##################
##CRMA : (MAGEA3, CSAG3, CSAG2,MAGEA2, MAGEA2B, CSAG1, MAGEA12, MAGEA6)
##CTLA-4 Blockade
CRMA_score = function(data,name){
  MGAEA_genes = c("MAGEA3", "CSAG3", "CSAG2","MAGEA2", "MAGEA2B", "CSAG1", "MAGEA12", "MAGEA6")
  intersect_genes = intersect(MGAEA_genes,rownames(data$TPM))
  expression = data$TPM[intersect_genes,]
  expression = as.data.frame(t(expression))
  expression["CRMA_score"] = apply(expression, 1, function(x){return(mean(as.numeric(x)))})
  expression["Sample"] = rownames(expression)
  expMarker = merge(data$Samples,expression[c("Sample","CRMA_score")],by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  g = ggboxplot(expMarker,x="Resp_NoResp",y='CRMA_score',color = "Resp_NoResp",
                palette = "npg",add = "jitter") + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "greater"))
  ggsave(paste0("Results/Results_wilcox_test/",name,"/CRMA_score.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(CRMA_score ~ Resp_NoResp, data = expMarker,alternative = "greater")
  auc = roc(Resp_NoResp ~ CRMA_score, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("Results/Results_AUC/",name,"/CRMA_score.png"),width = 6,height = 6, bg="white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
}


#######The immune resistance program
Immune_resistance_program_score = function(data,name){
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
  expMarker = merge(data$Samples,OE,by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  g = ggboxplot(expMarker,x="Resp_NoResp",y='The_immune_resistance_program',color = "Resp_NoResp",
                palette = "npg",add = "jitter") + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "greater"))
  ggsave(paste0("Results/Results_wilcox_test/",name,"/The_immune_resistance_program.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(The_immune_resistance_program ~ Resp_NoResp, data = expMarker,alternative = "greater")
  auc = roc(Resp_NoResp ~ The_immune_resistance_program, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("Results/Results_AUC/",name,"/The_immune_resistance_program.png"),width = 6,height = 6, bg="white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
}



############ EMT/Stroma_core signature: 8 genes were included in the EdgeSeq expression panel 
##Patients with high CD8 infiltration and low EMT/Stromal core gene expression had the highest response rates and longest PFS and OS
EMT_Stroma_core_marker = function(data,name){
  EMT_Stroma_core_signature = c("FLNA","EMP3","CALD1","FN1","FOXC2","LOX","FBN1","TNC")
  expression = data$TPM[intersect(EMT_Stroma_core_signature,rownames(data$TPM)),]
  expression = as.data.frame(t(expression))
  expression["EMT_Stroma_core_signature"] = apply(expression, 1, function(x){return(mean(as.numeric(x)))})
  expression["Sample"] = rownames(expression)
  expMarker = merge(data$Samples,expression[c("Sample","EMT_Stroma_core_signature")],by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  g = ggboxplot(expMarker,x="Resp_NoResp",y="EMT_Stroma_core_signature",color = "Resp_NoResp",
                palette = "npg",add = "jitter") + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "greater"))
  ggsave(paste0("Results/Results_wilcox_test/",name,"/EMT_Stroma_core_signature.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(EMT_Stroma_core_signature ~ Resp_NoResp, data = expMarker,alternative = "greater")
  auc = roc(Resp_NoResp ~ EMT_Stroma_core_signature, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("Results/Results_AUC/",name,"/EMT_Stroma_core_signature.png"),width = 6,height = 6, bg="white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
}


##########pan-fibroblast TGFβ response signature (F-TBRS), based on PCA
F_TBRS_score = function(data,name){
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
  expMarker = merge(data$Samples,expression,by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  g = ggboxplot(expMarker,x="Resp_NoResp",y="F_TBRS",color = "Resp_NoResp",
                palette = "npg",add = "jitter") + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "greater"))
  ggsave(paste0("Results/Results_wilcox_test/",name,"/F_TBRS.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(F_TBRS ~ Resp_NoResp, data = expMarker,alternative = "greater")
  auc = roc(Resp_NoResp ~ F_TBRS, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("Results/Results_AUC/",name,"/F_TBRS.png"),width = 6,height = 6, bg="white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
}




##########
##TMEscore
TMEscore = function(data,name){
  library(readr)
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
  expMarker = merge(data$Samples,expression,by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  g = ggboxplot(expMarker,x="Resp_NoResp",y="TMEscore",color = "Resp_NoResp",
                palette = "npg",add = "jitter") + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "less"))
  ggsave(paste0("Results/Results_wilcox_test/",name,"/TMEscore.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(TMEscore ~ Resp_NoResp, data = expMarker,alternative = "less")
  auc = roc(Resp_NoResp ~ TMEscore, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("Results/Results_AUC/",name,"/TMEscore.png"),width = 6,height = 6, bg="white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
}





#############
##risk score, 11 IRGs (immune-related genes)
Risk_score = function(data,name){
  IRGs = c("LEPR","PRLHR","NR2F2","PRL","NRP1","ANGPTL5","IGF1","TNFRSF10B","TNFRSF10A","PLAU","IFI30")
  coeff = c(0.32196,-0.64921,-0.32677,0.23573,0.39005,0.38166,-0.03522,0.02975,0.39830,0.14607,-0.68625)
  ##only compute the gene expression in the dataset
  expression = data$TPM[IRGs[which(is.element(IRGs,rownames(data$TPM)))],]
  riskScore = apply(expression, 2, function(x){return(sum(x*coeff[which(is.element(IRGs,rownames(data$TPM)))]))})
  expMarker = merge(data$Samples,data.frame(Sample=names(riskScore),RiskScore=riskScore),by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  g = ggboxplot(expMarker,x="Resp_NoResp",y="RiskScore",color = "Resp_NoResp",
                palette = "npg",add = "jitter") + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "greater"))
  ggsave(paste0("Results/Results_wilcox_test/",name,"/RiskScore.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(RiskScore ~ Resp_NoResp, data = expMarker,alternative = "greater")
  auc = roc(Resp_NoResp ~ RiskScore, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("Results/Results_AUC/",name,"/RiskScore.png"),width = 6,height = 6, bg="white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
}



############
#Tertiary lymphoid structures (TLS) gene signature score
TLS_score = function(data,name){
  TLS_gene = c("CD79B","CD1D","CCR6","LAT","SKAP1","CETP","EIF1AY","RBP5","PTGDS")
  TLS_gene = intersect(TLS_gene,rownames(data$TPM))
  expression = data$TPM[TLS_gene,]
  mean_expression = apply(expression, 2, mean)
  expMarker = merge(data$Samples,data.frame(Sample=names(mean_expression),TLS_score=as.numeric(mean_expression)),by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  g = ggboxplot(expMarker,x="Resp_NoResp",y='TLS_score',color = "Resp_NoResp",
                palette = "npg",add = "jitter") + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "less"))
  ggsave(paste0("Results/Results_wilcox_test/",name,"/TLS_score.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(TLS_score ~ Resp_NoResp, data = expMarker,alternative = "less")
  auc = roc(Resp_NoResp ~ TLS_score, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("Results/Results_AUC/",name,"/TLS_score.png"),width = 6,height = 6, bg="white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
}



##########
##CXCL9
CXCL9_marker = function(data,name){
  expression = data$TPM["CXCL9",]
  expMarker = merge(data$Samples,data.frame(Sample=names(expression),CXCL9=as.numeric(expression)),by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  g = ggboxplot(expMarker,x="Resp_NoResp",y='CXCL9',color = "Resp_NoResp",
                palette = "npg",add = "jitter") + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "less"))
  ggsave(paste0("Results/Results_wilcox_test/",name,"/CXCL9.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(CXCL9 ~ Resp_NoResp, data = expMarker,alternative = "less")
  auc = roc(Resp_NoResp ~ CXCL9, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("Results/Results_AUC/",name,"/CXCL9.png"),width = 6,height = 6, bg="white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
}




##MPS melanocytic plasticity signature, 45 genes
MPS_score = function(data,name){
  library(readr)
  MPS_gene_list <- read_delim("marker/MPS_gene_list.txt", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)
  expression = data$TPM[intersect(MPS_gene_list$`Gene Symbol`,rownames(data$TPM)),]
  MPS_postive_gene = MPS_gene_list$`Gene Symbol`[which(MPS_gene_list$`Sign in the signature`==1)]
  MPS_negative_gene = MPS_gene_list$`Gene Symbol`[which(MPS_gene_list$`Sign in the signature`==-1)]
  MPS_score = apply(expression, 2, function(x){return(sum(x[which(is.element(names(x),MPS_postive_gene))]) - sum(x[which(is.element(names(x),MPS_negative_gene))]))})
  expMarker = merge(data$Samples,data.frame(Sample=names(MPS_score),MPS_score=as.numeric(MPS_score)),by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  g = ggboxplot(expMarker,x="Resp_NoResp",y='MPS_score',color = "Resp_NoResp",
                palette = "npg",add = "jitter") + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "greater"))
  ggsave(paste0("Results/Results_wilcox_test/",name,"/MPS_score.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(MPS_score ~ Resp_NoResp, data = expMarker,alternative = "greater")
  auc = roc(Resp_NoResp ~ MPS_score, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("Results/Results_AUC/",name,"/MPS_score.png"),width = 6,height = 6, bg="white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
}



##########
#Renal 101 Immuno signature
Renal_101_Immuno_score = function(data,name){
  Renal_101_Immuno_signature = c("CD3G","CD3E","CD8B","THEMIS","TRAT1","GRAP2","CD247",
                                 "CD2","CD96","PRF1","CD6","IL7R","ITK","GPR18","EOMES",
                                 "SIT1","NLRC3","CD244","KLRD1","SH2D1A","CCL5","XCL2",
                                 "CST7","GFI1","KCNA3","PSTPIP1")
  Renal_101_Immuno_signature = intersect(Renal_101_Immuno_signature,rownames(data$TPM))
  expression = data$TPM[Renal_101_Immuno_signature,]
  mean_expression = apply(expression, 2, mean)
  expMarker = merge(data$Samples,data.frame(Sample=names(mean_expression),Renal_101_Immuno_signature=as.numeric(mean_expression)),by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  g = ggboxplot(expMarker,x="Resp_NoResp",y='Renal_101_Immuno_signature',color = "Resp_NoResp",
                palette = "npg",add = "jitter") + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "less"))
  ggsave(paste0("Results/Results_wilcox_test/",name,"/Renal_101_Immuno_signature.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(Renal_101_Immuno_signature ~ Resp_NoResp, data = expMarker,alternative = "less")
  auc = roc(Resp_NoResp ~ Renal_101_Immuno_signature, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("Results/Results_AUC/",name,"/Renal_101_Immuno_signature.png"),width = 6,height = 6, bg="white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
}




##########
##HRH1
HRH1_marker = function(data,name){
  expression = data$TPM["HRH1",]
  expMarker = merge(data$Samples,data.frame(Sample=names(expression),HRH1=as.numeric(expression)),by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  g = ggboxplot(expMarker,x="Resp_NoResp",y='HRH1',color = "Resp_NoResp",
                palette = "npg",add = "jitter") + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "greater"))
  ggsave(paste0("Results/Results_wilcox_test/",name,"/HRH1.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(HRH1 ~ Resp_NoResp, data = expMarker,alternative = "greater")
  auc = roc(Resp_NoResp ~ HRH1, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("Results/Results_AUC/",name,"/HRH1.png"),width = 6,height = 6, bg="white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
}



# ## generate TIDE input files
# tide_input = function(data,name){
#   data = log2(data$TPM+1)
#   data = scale(data)
#   Gene = rownames(data)
#   data = cbind(Gene,data)
#   write.table(data,file = paste0("./TIDE_input/",name,".txt"),quote = F,sep = "\t",row.names = F,col.names = T)
# }
# tide_input(data1,name1)
# tide_input(data2,name2)
# tide_input(data3,name3)
# tide_input(data4,name4)
# tide_input(data5,name5)
# tide_input(data6,name6)
#Gene = rownames(data7$TPM)
#write.table(cbind(Gene,data7$TPM),file = paste0("./TIDE_input/",name7,".txt"),quote = F,sep = "\t",row.names = F,col.names = T)
# tide_input(data8,name8)
# tide_input(data9,name9)
# tide_input(data10,name10)
# tide_input(data11,name11)
# tide_input(data12,name12)
# tide_input(data13,name13)
# tide_input(data14,name14)
# tide_input(data15,name15)
# tide_input(data16,name16)
# tide_input(Gide_PRE_data,"Gide_PRE_data")
# tide_input(Gide_ON_data,"Gide_ON_data")
# tide_input(Lee_PRE_data,"Lee_PRE_data")
# tide_input(Lee_ON_data,"Lee_ON_data")
# tide_input(Nathanson_PRE_data,"Nathanson_PRE_data")
# tide_input(Nathanson_ON_data,"Nathanson_ON_data")
# tide_input(Riaz_PRE_data,"Riaz_PRE_data")
# tide_input(Riaz_ON_data,"Riaz_ON_data")
# tide_input(Gide_MONO_data,"Gide_MONO_data")
# tide_input(Gide_COMBINE_data,"Gide_COMBINE_data")
# tide_input(Riaz_NAIVE_data,"Riaz_NAIVE_data")
# tide_input(Riaz_EXPOSURE_data,"Riaz_EXPOSURE_data")
# tide_input(Liu_NAIVE_data,"Liu_NAIVE_data")
# tide_input(Liu_EXPOSURE_data,"Liu_EXPOSURE_data")


TIDE = function(data,name){
  library(readr)
  file_path = paste0("TIDE_output/",name,".csv")
  res <- read_csv(file_path)
  expMarker = merge(data$Samples,data.frame(Sample=res$Patient,TIDE=res$TIDE),by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  g = ggboxplot(expMarker,x="Resp_NoResp",y='TIDE',color = "Resp_NoResp",
                palette = "npg",add = "jitter") + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "greater"))
  ggsave(paste0("Results/Results_wilcox_test/",name,"/TIDE.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(TIDE ~ Resp_NoResp, data = expMarker,alternative = "greater")
  auc = roc(Resp_NoResp ~ TIDE, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("Results/Results_AUC/",name,"/TIDE.png"),width = 6,height = 6, bg="white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
}

