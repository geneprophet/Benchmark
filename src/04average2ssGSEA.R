##select signatures which contains >5 genes to compute ssGSEA instead of average expression
library("GSVA")
library("ggpubr")
library(pROC)
my_comparisons = list(c("No_Response","Response"))

IFN_gamma_marker_ssGSEA = function(data,name){
  IFN_gamma = c("IDO1","CXCL10","CXCL9","HLA-DRA","STAT1","IFNG")
  gene_set = list()
  gene_set[["IFN_gamma"]] = IFN_gamma
  gsva.es <- gsva(data$TPM, gene_set, method="ssgsea", verbose=T)
  expMarker = merge(data$Samples,data.frame(Sample=colnames(gsva.es),IFN_gamma_ssGSEA=as.numeric(gsva.es)),by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  g = ggboxplot(expMarker,x="Resp_NoResp",y="IFN_gamma_ssGSEA",color = "Resp_NoResp",
                  palette = "npg",add = "jitter") + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "less"))
  ggsave(paste0("Results/Results_wilcox_test/",name,"/IFN_gamma_ssGSEA.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(IFN_gamma_ssGSEA ~ Resp_NoResp, data = expMarker,alternative = "less")
  auc = roc(Resp_NoResp ~ IFN_gamma_ssGSEA, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("Results/Results_AUC/",name,"/IFN_gamma_ssGSEA.png"),width = 6,height = 6, bg="white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
}

####Expanded_immune_gene_signature
Expanded_immune_gene_signature_ssGSEA= function(data,name){
  Expanded_immune_gene_signature = c("CD3D","IDO1","CIITA","CD3E","CCL5","GZMK","CD2","HLA-DRA","CXCL13","IL2RG","NKG7","HLA-E","CXCR6","LAG3","TAGAP","CXCL10","STAT1","GZMB")
  gene_set = list()
  gene_set[["Expanded_immune_gene_signature"]] = Expanded_immune_gene_signature
  gsva.es <- gsva(data$TPM, gene_set, method="ssgsea", verbose=T)
  expMarker = merge(data$Samples,data.frame(Sample=colnames(gsva.es),Expanded_immune_gene_ssGSEA=as.numeric(gsva.es)),by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  g = ggboxplot(expMarker,x="Resp_NoResp",y="Expanded_immune_gene_ssGSEA",color = "Resp_NoResp",
                palette = "npg",add = "jitter") + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "less"))
  ggsave(paste0("Results/Results_wilcox_test/",name,"/Expanded_immune_gene_ssGSEA.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(Expanded_immune_gene_ssGSEA ~ Resp_NoResp, data = expMarker,alternative = "less")
  auc = roc(Resp_NoResp ~ Expanded_immune_gene_ssGSEA, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("Results/Results_AUC/",name,"/Expanded_immune_gene_ssGSEA.png"),width = 6,height = 6, bg="white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
}

####T_cell_inflamed_GEP_score
T_cell_inflamed_GEP_ssGSEA = function(data,name){
  T_cell_inflamed_GEP = c("TIGIT","PDCD1LG2","CD27","CD8A","LAG3","CD274","CXCR6","CMKLR1","NKG7","CCL5","PSMB10","IDO1","CXCL9","HLA-DQA1","CD276","STAT1","HLA-DRB1","HLA-E")
  gene_set = list()
  gene_set[["T_cell_inflamed_GEP"]] = T_cell_inflamed_GEP
  gsva.es <- gsva(data$TPM, gene_set, method="ssgsea", verbose=T)
  expMarker = merge(data$Samples,data.frame(Sample=colnames(gsva.es),T_cell_inflamed_GEP_ssGSEA=as.numeric(gsva.es)),by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  g = ggboxplot(expMarker,x="Resp_NoResp",y="T_cell_inflamed_GEP_ssGSEA",color = "Resp_NoResp",
                palette = "npg",add = "jitter") + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "less"))
  ggsave(paste0("Results/Results_wilcox_test/",name,"/T_cell_inflamed_GEP_ssGSEA.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(T_cell_inflamed_GEP_ssGSEA ~ Resp_NoResp, data = expMarker,alternative = "less")
  auc = roc(Resp_NoResp ~ T_cell_inflamed_GEP_ssGSEA, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("Results/Results_AUC/",name,"/T_cell_inflamed_GEP_ssGSEA.png"),width = 6,height = 6, bg="white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
  
}

CRMA_ssGSEA = function(data,name){
  MGAEA_genes = c("MAGEA3", "CSAG3", "CSAG2","MAGEA2", "MAGEA2B", "CSAG1", "MAGEA12", "MAGEA6")
  gene_set = list()
  gene_set[["CRMA"]] = MGAEA_genes
  gsva.es <- gsva(data$TPM, gene_set, method="ssgsea", verbose=T)
  expMarker = merge(data$Samples,data.frame(Sample=colnames(gsva.es),CRMA_ssGSEA=as.numeric(gsva.es)),by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  g = ggboxplot(expMarker,x="Resp_NoResp",y="CRMA_ssGSEA",color = "Resp_NoResp",
                palette = "npg",add = "jitter") + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "greater"))
  ggsave(paste0("Results/Results_wilcox_test/",name,"/CRMA_ssGSEA.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(CRMA_ssGSEA ~ Resp_NoResp, data = expMarker,alternative = "greater")
  auc = roc(Resp_NoResp ~ CRMA_ssGSEA, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("Results/Results_AUC/",name,"/CRMA_ssGSEA.png"),width = 6,height = 6, bg="white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
}

EMT_Stroma_core_ssGSEA = function(data,name){
  EMT_Stroma_core_signature = c("FLNA","EMP3","CALD1","FN1","FOXC2","LOX","FBN1","TNC")
  gene_set = list()
  gene_set[["EMT_Stroma_core_signature"]] = EMT_Stroma_core_signature
  gsva.es <- gsva(data$TPM, gene_set, method="ssgsea", verbose=T)
  expMarker = merge(data$Samples,data.frame(Sample=colnames(gsva.es),EMT_Stroma_core_ssGSEA=as.numeric(gsva.es)),by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  g = ggboxplot(expMarker,x="Resp_NoResp",y="EMT_Stroma_core_ssGSEA",color = "Resp_NoResp",
                palette = "npg",add = "jitter") + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "greater"))
  ggsave(paste0("Results/Results_wilcox_test/",name,"/EMT_Stroma_core_ssGSEA.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(EMT_Stroma_core_ssGSEA ~ Resp_NoResp, data = expMarker,alternative = "greater")
  auc = roc(Resp_NoResp ~ EMT_Stroma_core_ssGSEA, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("Results/Results_AUC/",name,"/EMT_Stroma_core_ssGSEA.png"),width = 6,height = 6, bg="white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
}

F_TBRS_ssGSEA = function(data,name){
  F_TBRS_genes = c("ACTA2", "ACTG2", "ADAM12", "ADAM19", "CNN1", "COL4A1", "CCN2", "CTPS1",
                   "RFLNB", "FSTL3", "HSPB1", "IGFBP3", "PXDC1", "SEMA7A", "SH3PXD2A", "TAGLN", 
                   "TGFBI", "TNS1", "TPM1")
  
  gene_set = list()
  gene_set[["F_TBRS_genes"]] = F_TBRS_genes
  gsva.es <- gsva(data$TPM, gene_set, method="ssgsea", verbose=T)
  expMarker = merge(data$Samples,data.frame(Sample=colnames(gsva.es),F_TBRS_ssGSEA=as.numeric(gsva.es)),by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  g = ggboxplot(expMarker,x="Resp_NoResp",y="F_TBRS_ssGSEA",color = "Resp_NoResp",
                palette = "npg",add = "jitter") + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "greater"))
  ggsave(paste0("Results/Results_wilcox_test/",name,"/F_TBRS_ssGSEA.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(F_TBRS_ssGSEA ~ Resp_NoResp, data = expMarker,alternative = "greater")
  auc = roc(Resp_NoResp ~ F_TBRS_ssGSEA, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("Results/Results_AUC/",name,"/F_TBRS_ssGSEA.png"),width = 6,height = 6, bg="white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
  
}

Risk_ssGSEA = function(data,name){
  IRGs = c("LEPR","PRLHR","NR2F2","PRL","NRP1","ANGPTL5","IGF1","TNFRSF10B","TNFRSF10A","PLAU","IFI30")
  gene_set = list()
  gene_set[["IRGs"]] = IRGs
  gsva.es <- gsva(data$TPM, gene_set, method="ssgsea", verbose=T)
  expMarker = merge(data$Samples,data.frame(Sample=colnames(gsva.es),RiskScore_ssGSEA=as.numeric(gsva.es)),by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  g = ggboxplot(expMarker,x="Resp_NoResp",y="RiskScore_ssGSEA",color = "Resp_NoResp",
                palette = "npg",add = "jitter") + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "greater"))
  ggsave(paste0("Results/Results_wilcox_test/",name,"/RiskScore_ssGSEA.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(RiskScore_ssGSEA ~ Resp_NoResp, data = expMarker,alternative = "greater")
  auc = roc(Resp_NoResp ~ RiskScore_ssGSEA, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("Results/Results_AUC/",name,"/RiskScore_ssGSEA.png"),width = 6,height = 6, bg="white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
}

TLS_score_ssGSEA = function(data,name){
  TLS_gene = c("CD79B","CD1D","CCR6","LAT","SKAP1","CETP","EIF1AY","RBP5","PTGDS")
  gene_set = list()
  gene_set[["TLS_gene"]] = TLS_gene
  gsva.es <- gsva(data$TPM, gene_set, method="ssgsea", verbose=T)
  expMarker = merge(data$Samples,data.frame(Sample=colnames(gsva.es),TLS_score_ssGSEA=as.numeric(gsva.es)),by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  g = ggboxplot(expMarker,x="Resp_NoResp",y="TLS_score_ssGSEA",color = "Resp_NoResp",
                palette = "npg",add = "jitter") + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "less"))
  ggsave(paste0("Results/Results_wilcox_test/",name,"/TLS_score_ssGSEA.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(TLS_score_ssGSEA ~ Resp_NoResp, data = expMarker,alternative = "less")
  auc = roc(Resp_NoResp ~ TLS_score_ssGSEA, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("Results/Results_AUC/",name,"/TLS_score_ssGSEA.png"),width = 6,height = 6, bg="white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
}

Renal_101_Immuno_ssGSEA = function(data,name){
  Renal_101_Immuno_signature = c("CD3G","CD3E","CD8B","THEMIS","TRAT1","GRAP2","CD247",
                                 "CD2","CD96","PRF1","CD6","IL7R","ITK","GPR18","EOMES",
                                 "SIT1","NLRC3","CD244","KLRD1","SH2D1A","CCL5","XCL2",
                                 "CST7","GFI1","KCNA3","PSTPIP1")
  gene_set = list()
  gene_set[["Renal_101_Immuno_signature"]] = Renal_101_Immuno_signature
  gsva.es <- gsva(data$TPM, gene_set, method="ssgsea", verbose=T)
  expMarker = merge(data$Samples,data.frame(Sample=colnames(gsva.es),Renal_101_Immuno_ssGSEA=as.numeric(gsva.es)),by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  g = ggboxplot(expMarker,x="Resp_NoResp",y="Renal_101_Immuno_ssGSEA",color = "Resp_NoResp",
                palette = "npg",add = "jitter") + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "less"))
  ggsave(paste0("Results/Results_wilcox_test/",name,"/Renal_101_Immuno_ssGSEA.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(Renal_101_Immuno_ssGSEA ~ Resp_NoResp, data = expMarker,alternative = "less")
  auc = roc(Resp_NoResp ~ Renal_101_Immuno_ssGSEA, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("Results/Results_AUC/",name,"/Renal_101_Immuno_ssGSEA.png"),width = 6,height = 6, bg="white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
}


##compare average and ssGSEA, which is better?
source("src/marker_gene_expression.R")

benchmark_comparision = function(data,name){
  print(paste0("Begin benchmark of dataset: ",name))
  dir.create(paste0("Results/Results_wilcox_test/",name))
  dir.create(paste0("Results/Results_AUC/",name))
  benchmark_result = rbind(IFN_gamma_marker_ssGSEA(data,name),
                       IFN_gamma_marker(data,name),
                       Expanded_immune_gene_signature_ssGSEA(data,name),
                       Expanded_immune_gene_signature_marker(data,name),
                       T_cell_inflamed_GEP_ssGSEA(data,name),
                       T_cell_inflamed_GEP_score(data,name),
                       CRMA_ssGSEA(data,name),
                       CRMA_score(data,name),
                       EMT_Stroma_core_ssGSEA(data,name),
                       EMT_Stroma_core_marker(data,name),
                       F_TBRS_ssGSEA(data,name),
                       F_TBRS_score(data,name),
                       Risk_ssGSEA(data,name),
                       Risk_score(data,name),
                       TLS_score_ssGSEA(data,name),
                       TLS_score(data,name),
                       Renal_101_Immuno_ssGSEA(data,name),
                       Renal_101_Immuno_score(data,name))
  
  rownames(benchmark_result) = c("IFN_gamma_ssGSEA","IFN_gamma","Expanded_immune_gene_ssGSEA","Expanded_immune_gene_signature",
                                 "T_cell_inflamed_GEP_ssGSEA","T_cell_inflamed_GEP_score","CRMA_ssGSEA","CRMA_score",
                                 "EMT_Stroma_core_ssGSEA","EMT_Stroma_core_signature",
                                 "F_TBRS_ssGSEA","F_TBRS","RiskScore_ssGSEA","RiskScore","TLS_score_ssGSEA","TLS_score",
                                 "Renal_101_Immuno_ssGSEA","Renal_101_Immuno_signature")
  
  Dataset=rep(name,18)
  result = as.data.frame(cbind(benchmark_result,Dataset))
  result$Biomarker = rownames(benchmark_result)
  return(result)
}

load("./ICB_data/Braun et al/Braun_data.Rdata")
data1=Braun_data; name1="Braun_2020"
load("./ICB_data/Du et al/MGH_PRE_data.RData")
data2=MGH_PRE_data; name2="MGH_PRE_2021"
load("./ICB_data/Du et al/MGH_ON_data.RData")
data3=MGH_ON_data; name3="MGH_ON_2021"
load("./ICB_data/Gide et al/Gide_data.RData")
data4=Gide_data; name4="Gide_2019"
load("./ICB_data/Hugo et al/Hugo_data.RData")
data5=Hugo_data; name5="Hugo_2016"
load("./ICB_data/Jung et al/Jung_data.Rdata")
data6=Jung_data; name6="Jung_2019"
load("./ICB_data/Kim et al/Kim_data.RData")
data7=Kim_data; name7="Kim_2018"
load("./ICB_data/Lee et al/Lee_data.RData")
data8=Lee_data; name8="Lee_2020"
load("./ICB_data/Liu et al/Liu_data.Rdata")
data9=Liu_data; name9="Liu_2019"
load("./ICB_data/Mariathasan et al/Mariathasan_data.Rdata")
data10=Mariathasan_data; name10="Mariathasan_2018"
load("./ICB_data/Miao et al/Miao_data.Rdata")
data11=Miao_data; name11="Miao_2018"
load("./ICB_data/Motzer et al/Motzer_data.Rdata")
data12=Motzer_data; name12="Motzer_2020"
load("./ICB_data/Nathanson et al/Nathanson_data.Rdata")
data13=Nathanson_data; name13="Nathanson_2017"
load("./ICB_data/Riaz et al/Riaz_data.RData")
data14=Riaz_data; name14="Riaz_2017"
load("./ICB_data/Snyder et al/Snyder_data.Rdata")
data15=Snyder_data; name15="Snyder_2017"
load("./ICB_data/VanAllen et al/VanAllen_data.RData")
data16=VanAllen_data; name16="VanAllen_2015"

load("ICB_data/Gide et al/Gide_PRE_data.Rdata")
load("ICB_data/Gide et al/Gide_ON_data.Rdata")
data17=Gide_PRE_data; name17="Gide_PRE_2019"
data18=Gide_ON_data; name18="Gide_ON_2019"
load("ICB_data/Lee et al/Lee_PRE_data.Rdata")
load("ICB_data/Lee et al/Lee_ON_data.Rdata")
data19=Lee_PRE_data; name19="Lee_PRE_2020"
data20=Lee_ON_data; name20="Lee_ON_2020"
load("ICB_data/Riaz et al/Riaz_PRE_data.RData")
load("ICB_data/Riaz et al/Riaz_ON_data.Rdata")
data21=Riaz_PRE_data; name21="Riaz_PRE_2017"
data22=Riaz_ON_data; name22="Riaz_ON_2017"

load("ICB_data/Gide et al/Gide_MONO_data.Rdata")
load("ICB_data/Gide et al/Gide_COMBINE_data.Rdata")
data23 = Gide_MONO_data; name23 = "Gide_MONO_2019"
data24 = Gide_COMBINE_data; name24 = "Gide_COMBINE_2019"
load("ICB_data/Riaz et al/Riaz_NAIVE_data.Rdata")
load("ICB_data/Riaz et al/Riaz_EXPOSURE_data.Rdata")
data25 = Riaz_NAIVE_data; name25 = "Riaz_NAIVE_2017"
data26 = Riaz_EXPOSURE_data; name26 = "Riaz_EXPOSURE_2017"
load("ICB_data/Liu et al/Liu_NAIVE_data.Rdata")
load("ICB_data/Liu et al/Liu_EXPOSURE_data.Rdata")
data27 = Liu_NAIVE_data; name27 = "Liu_NAIVE_2019"
data28 = Liu_EXPOSURE_data; name28 = "Liu_EXPOSURE_2019"




res = rbind(benchmark_comparision(data1,name1),
            benchmark_comparision(data2,name2),
            benchmark_comparision(data3,name3),
            benchmark_comparision(data4,name4),
            benchmark_comparision(data5,name5),
            benchmark_comparision(data6,name6),
            benchmark_comparision(data7,name7),
            benchmark_comparision(data8,name8),
            benchmark_comparision(data9,name9),
            benchmark_comparision(data10,name10),
            benchmark_comparision(data11,name11),
            benchmark_comparision(data12,name12),
            benchmark_comparision(data13,name13),
            benchmark_comparision(data14,name14),
            benchmark_comparision(data15,name15),
            benchmark_comparision(data16,name16),
            benchmark_comparision(data17,name17),
            benchmark_comparision(data18,name18),
            benchmark_comparision(data19,name19),
            benchmark_comparision(data20,name20),
            benchmark_comparision(data21,name21),
            benchmark_comparision(data22,name22),
            benchmark_comparision(data23,name23),
            benchmark_comparision(data24,name24),
            benchmark_comparision(data25,name25),
            benchmark_comparision(data26,name26),
            benchmark_comparision(data27,name27),
            benchmark_comparision(data28,name28))

res$Dataset[which(res$Dataset=="Braun_2020")] = "ccRCC_124_Braun"
res$Dataset[which(res$Dataset=="MGH_PRE_2021")] = "Melanoma_19_MGH_PRE"
res$Dataset[which(res$Dataset=="MGH_ON_2021")] = "Melanoma_31_MGH_ON"
res$Dataset[which(res$Dataset=="Gide_2019")] = "Melanoma_90_Gide"
res$Dataset[which(res$Dataset=="Hugo_2016")] = "Melanoma_26_Hugo"
res$Dataset[which(res$Dataset=="Jung_2019")] = "NSCLC_27_Jung"
res$Dataset[which(res$Dataset=="Kim_2018")] = "Gastric-Cancer_45_Kim"
res$Dataset[which(res$Dataset=="Lee_2020")] = "Melanoma_79_Lee"
res$Dataset[which(res$Dataset=="Liu_2019")] = "Melanoma_121_Liu"
res$Dataset[which(res$Dataset=="Mariathasan_2018")] = "Urothelial-Cancer_348_Mariathasan"
res$Dataset[which(res$Dataset=="Miao_2018")] = "ccRCC_33_Miao"
res$Dataset[which(res$Dataset=="Motzer_2020")] = "NSCLC_354_Motzer"
res$Dataset[which(res$Dataset=="Nathanson_2017")] = "Melanoma_24_Nathanson"
res$Dataset[which(res$Dataset=="Riaz_2017")] = "Melanoma_103_Riaz"
res$Dataset[which(res$Dataset=="Snyder_2017")] = "Urothelial-Cancer_26_Snyder"
res$Dataset[which(res$Dataset=="VanAllen_2015")] = "Melanoma_42_VanAllen"

res$Dataset[which(res$Dataset=="Gide_PRE_2019")] = "Gide_PRE_Melanoma_72"
res$Dataset[which(res$Dataset=="Gide_ON_2019")] = "Gide_ON_Melanoma_18"
res$Dataset[which(res$Dataset=="Riaz_PRE_2017")] = "Riaz_PRE_Melanoma_49"
res$Dataset[which(res$Dataset=="Riaz_ON_2017")] = "Riaz_ON_Melanoma_54"
res$Dataset[which(res$Dataset=="Lee_PRE_2020")] = "Lee_PRE_Melanoma_44"
res$Dataset[which(res$Dataset=="Lee_ON_2020")] = "Lee_ON_Melanoma_35"
res$Dataset[which(res$Dataset=="Gide_MONO_2019")] = "Gide_MONO_Melanoma_50"
res$Dataset[which(res$Dataset=="Gide_COMBINE_2019")] = "Gide_COMBINE_Melanoma_40"
res$Dataset[which(res$Dataset=="Riaz_NAIVE_2017")] = "Riaz_NAIVE_Melanoma_45"
res$Dataset[which(res$Dataset=="Riaz_EXPOSURE_2017")] = "Riaz_EXPOSURE_Melanoma_58"
res$Dataset[which(res$Dataset=="Liu_NAIVE_2019")] = "Liu_NAIVE_Melanoma_74"
res$Dataset[which(res$Dataset=="Liu_EXPOSURE_2019")] = "Liu_EXPOSURE_Melanoma_47"


res$Dataset = as.factor(res$Dataset)
res$Biomarker = as.factor(res$Biomarker)
res$p_value = as.numeric(res$p_value)
res$AUC = as.numeric(res$AUC)

source("./src/visualization.R")
p = dotplot(res)
ggsave(p,filename = "./figures/average2ssgsea_benchmark2.pdf",width = 16,height = 8)

save(list = ls(),file = "Results/average2ssgsea_results2.Rdata")



