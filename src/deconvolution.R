##deconvolution
# ##############
# CD8 T cells 
# ##xCell: 64 cell type signatures
# # devtools::install_github('dviraran/xCell')
CD8T_xCell = function(data,name){
  print("xCell")
  library(xCell)
  xCell = xCellAnalysis(data$TPM)
  expMarker = merge(data$Samples,data.frame(Sample=names(xCell["CD8+ T-cells",]),CD8T_xCell=as.numeric(xCell["CD8+ T-cells",])),by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  g = ggboxplot(expMarker,x="Resp_NoResp",y="CD8T_xCell",color = "Resp_NoResp",palette = "npg",add = "jitter") + 
    stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "less"))
  ggsave(paste0("Results/Results_wilcox_test/",name,"/CD8T_xCell.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(CD8T_xCell ~ Resp_NoResp, data = expMarker,alternative = "less")
  auc = roc(Resp_NoResp ~ CD8T_xCell, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("Results/Results_AUC/",name,"/CD8T_xCell.png"),width = 6,height = 6, bg="white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
}


CD8T_MCPcounter = function(data,name){
  print("MCPcounter")
  library(MCPcounter)
  MCPcounter_score = MCPcounter.estimate(data$TPM,featuresType="HUGO_symbols")
  expMarker = merge(data$Samples,data.frame(Sample=names(MCPcounter_score["CD8 T cells",]),CD8T_MCPcounter=as.numeric(MCPcounter_score["CD8 T cells",])),by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  g = ggboxplot(expMarker,x="Resp_NoResp",y="CD8T_MCPcounter",color = "Resp_NoResp",palette = "npg",add = "jitter") + 
    stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "less"))
  ggsave(paste0("Results/Results_wilcox_test/",name,"/CD8T_MCPcounter.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(CD8T_MCPcounter ~ Resp_NoResp, data = expMarker,alternative = "less")
  auc = roc(Resp_NoResp ~ CD8T_MCPcounter, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("Results/Results_AUC/",name,"/CD8T_MCPcounter.png"),width = 6,height = 6, bg="white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
}



#generate CIBERSORTx input file
# CIBERSORTx_input = function(data,name){
#   data = as.data.frame(data$TPM)
#   Gene = rownames(data)
#   data = cbind(Gene,data)
#   write.table(data,file = paste0("./CIBERSORTx_input/",name,".txt"),quote = F,sep = "\t",row.names = F,col.names = T)
# }
# CIBERSORTx_input(data1,name1)
# CIBERSORTx_input(data2,name2)
# CIBERSORTx_input(data3,name3)
# CIBERSORTx_input(data4,name4)
# CIBERSORTx_input(data5,name5)
# CIBERSORTx_input(data6,name6)
# CIBERSORTx_input(data7,name7)
# CIBERSORTx_input(data8,name8)
# CIBERSORTx_input(data9,name9)
# CIBERSORTx_input(data10,name10)
# CIBERSORTx_input(data11,name11)
# CIBERSORTx_input(data12,name12)
# CIBERSORTx_input(data13,name13)
# CIBERSORTx_input(data14,name14)
# CIBERSORTx_input(data15,name15)
# CIBERSORTx_input(data16,name16)
# CIBERSORTx_input(Gide_PRE_data,"Gide_PRE_data")
# CIBERSORTx_input(Gide_ON_data,"Gide_ON_data")
# CIBERSORTx_input(Lee_PRE_data,"Lee_PRE_data")
# CIBERSORTx_input(Lee_ON_data,"Lee_ON_data")
# CIBERSORTx_input(Nathanson_PRE_data,"Nathanson_PRE_data")
# CIBERSORTx_input(Nathanson_ON_data,"Nathanson_ON_data")
# CIBERSORTx_input(Riaz_PRE_data,"Riaz_PRE_data")
# CIBERSORTx_input(Riaz_ON_data,"Riaz_ON_data")
# CIBERSORTx_input(Gide_MONO_data,"Gide_MONO_data")
# CIBERSORTx_input(Gide_COMBINE_data,"Gide_COMBINE_data")
# CIBERSORTx_input(Riaz_NAIVE_data,"Riaz_NAIVE_data")
# CIBERSORTx_input(Riaz_EXPOSURE_data,"Riaz_EXPOSURE_data")
# CIBERSORTx_input(Liu_NAIVE_data,"Liu_NAIVE_data")
# CIBERSORTx_input(Liu_EXPOSURE_data,"Liu_EXPOSURE_data")


##read the CIBERSORTx_output files

CD8T_CIBERSORTx = function(data,name){
  file_path = paste0("CIBERSORTx_output/CIBERSORTx_",name,"_Results.txt")
  library(readr)
  CIBERSORTx_result <- read_delim(file_path,
                                  delim = "\t",
                                  escape_double = FALSE, 
                                  trim_ws = TRUE)
  expMarker = merge(data$Samples,data.frame(Sample=CIBERSORTx_result$Mixture,CD8T_CIBERSORTx=as.numeric(CIBERSORTx_result$`T cells CD8`)),by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  g = ggboxplot(expMarker,x="Resp_NoResp",y="CD8T_CIBERSORTx",color = "Resp_NoResp",palette = "npg",add = "jitter") + 
    stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "less"))
  ggsave(paste0("Results/Results_wilcox_test/",name,"/CD8T_CIBERSORTx.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(CD8T_CIBERSORTx ~ Resp_NoResp, data = expMarker,alternative = "less")
  auc = roc(Resp_NoResp ~ CD8T_CIBERSORTx, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("Results/Results_AUC/",name,"/CD8T_CIBERSORTx.png"),width = 6,height = 6, bg="white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
}



##############
##immunoscore construct based on CIBERSORTx_output
##immunoscore = (1.13 × fraction level of naive B cells) + (1.36 × fraction level of memory B cells) + (5.92 × fraction level of eosinophils) + (9.70 × fraction level of follicular helper  T cells) + (15.34 × fraction level of Tregs) - (1.14 × fraction level of M0 macrophages) - (2.31 × fraction level of plasma cells) - (4.52 × fraction level of γδT  cell).
Immunoscore_CIBERSORTx = function(data,name){
  file_path = paste0("CIBERSORTx_output/CIBERSORTx_",name,"_Results.txt")
  library(readr)
  CIBERSORTx_result <- read_delim(file_path,
                                  delim = "\t",
                                  escape_double = FALSE, 
                                  trim_ws = TRUE)
  
  immunoscore = 1.13*CIBERSORTx_result$`B cells naive` + 1.36*CIBERSORTx_result$`B cells memory` + 5.92*CIBERSORTx_result$Eosinophils + 
                  9.70*CIBERSORTx_result$`T cells follicular helper` + 15.34*CIBERSORTx_result$`T cells regulatory (Tregs)` -
                  1.14*CIBERSORTx_result$`Macrophages M0` - 2.31*CIBERSORTx_result$`Plasma cells` - 4.52*CIBERSORTx_result$`T cells gamma delta`
  expMarker = merge(data$Samples,data.frame(Sample=CIBERSORTx_result$Mixture,Immunoscore_CIBERSORTx=immunoscore),by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  g = ggboxplot(expMarker,x="Resp_NoResp",y="Immunoscore_CIBERSORTx",color = "Resp_NoResp",palette = "npg",add = "jitter") + 
    stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "less"))
  ggsave(paste0("Results/Results_wilcox_test/",name,"/Immunoscore_CIBERSORTx.png"),width = 6,height = 6, bg="white")
  res = wilcox.test(Immunoscore_CIBERSORTx ~ Resp_NoResp, data = expMarker,alternative = "less")
  auc = roc(Resp_NoResp ~ Immunoscore_CIBERSORTx, data = expMarker,auc = T)
  ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
  ggsave(paste0("Results/Results_AUC/",name,"/Immunoscore_CIBERSORTx.png"),width = 6,height = 6, bg="white")
  result=c(res$p.value,as.numeric(auc$auc))
  names(result) = c("p_value","AUC")
  return(result)
}



#############
##EcoTyper: CE9
EcoTyper = function(data,name){
  file_path = paste0("EcoTyper_output/",name,"_ecotyper_output/Carcinoma_Ecotypes/Ecotype_Assignment.txt")
  library(readr)
  Ecotype_Assignment <- read_delim(file_path, 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)
  expMarker = merge(data$Samples,data.frame(Sample=Ecotype_Assignment$ID,Ecotype=Ecotype_Assignment$`Carcinoma Ecotype`),by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  expMarker = expMarker[which(!is.na(expMarker$Ecotype)),]
  
  a = length(expMarker$Sample[which(expMarker$Resp_NoResp=="Response" & expMarker$Ecotype=="CE9")])
  b = length(expMarker$Sample[which(expMarker$Resp_NoResp=="No_Response" & expMarker$Ecotype=="CE9")])
  c = length(expMarker$Sample[which(expMarker$Resp_NoResp=="Response" & expMarker$Ecotype!="CE9")])
  d = length(expMarker$Sample[which(expMarker$Resp_NoResp=="No_Response" & expMarker$Ecotype!="CE9")])
  matrix(c(a,b,c,d),nrow = 2)
  res = fisher.test(matrix(c(a,b,c,d),nrow = 2),alternative = "greater")
  res$p.value
  ## no AUC for EcoTyper, then return AUC as 0.5
  result=c(res$p.value,0.5)
  names(result) = c("p_value","AUC")
  return(result)
}


# ##MFP
# ##generate MFP input
# MFP_input = function(data,name){
#   ## if the expression matrix exist the minus value,
#   if(length(which(data$TPM<0))>0){
#     tpm = t(data$TPM)
#     tpm = as.data.frame(tpm)
#     genes = rownames(tpm)
#     tpm = cbind(genes,tpm)
#     write.table(tpm,file = paste0("./MFP_input/",name,".txt"),quote = F,sep = "\t",row.names = F,col.names = T)
# 
#   }else{
#     tpm = t(log2(data$TPM+1))
#     tpm = as.data.frame(tpm)
#     genes = rownames(tpm)
#     tpm = cbind(genes,tpm)
#     write.table(tpm,file = paste0("./MFP_input/",name,".txt"),quote = F,sep = "\t",row.names = F,col.names = T)
# 
#   }
# }
# MFP_input(data1,name1)
# MFP_input(data2,name2)
# MFP_input(data3,name3)
# MFP_input(data4,name4)
# MFP_input(data5,name5)
# MFP_input(data6,name6)
# MFP_input(data7,name7)
# MFP_input(data8,name8)
# MFP_input(data9,name9)
# MFP_input(data10,name10)
# MFP_input(data11,name11)
# MFP_input(data12,name12)
# MFP_input(data13,name13)
# MFP_input(data14,name14)
# MFP_input(data15,name15)
# MFP_input(data16,name16)
# MFP_input(Gide_PRE_data,"Gide_PRE_data")
# MFP_input(Gide_ON_data,"Gide_ON_data")
# MFP_input(Lee_PRE_data,"Lee_PRE_data")
# MFP_input(Lee_ON_data,"Lee_ON_data")
# MFP_input(Nathanson_PRE_data,"Nathanson_PRE_data")
# MFP_input(Nathanson_ON_data,"Nathanson_ON_data")
# MFP_input(Riaz_PRE_data,"Riaz_PRE_data")
# MFP_input(Riaz_ON_data,"Riaz_ON_data")
# MFP_input(Gide_MONO_data,"Gide_MONO_data")
# MFP_input(Gide_COMBINE_data,"Gide_COMBINE_data")
# MFP_input(Riaz_NAIVE_data,"Riaz_NAIVE_data")
# MFP_input(Riaz_EXPOSURE_data,"Riaz_EXPOSURE_data")
# MFP_input(Liu_NAIVE_data,"Liu_NAIVE_data")
# MFP_input(Liu_EXPOSURE_data,"Liu_EXPOSURE_data")

MFP = function(data,name){
  file_path = paste0("MFP_output/",name,".txt")
  library(readr)
  MFP_result <- read_delim(file_path,
                    delim = "\t", escape_double = FALSE, 
                    col_names = FALSE, trim_ws = TRUE)
  
  expMarker = merge(data$Samples,data.frame(Sample=MFP_result$X1,MFP=MFP_result$X2),by="Sample")
  expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
  expMarker = expMarker[which(!is.na(expMarker$MFP)),]
  
  a = length(expMarker$Sample[which(expMarker$Resp_NoResp=="Response" & expMarker$MFP=="IE")])
  b = length(expMarker$Sample[which(expMarker$Resp_NoResp=="No_Response" & expMarker$MFP=="IE")])
  c = length(expMarker$Sample[which(expMarker$Resp_NoResp=="Response" & expMarker$MFP!="IE")])
  d = length(expMarker$Sample[which(expMarker$Resp_NoResp=="No_Response" & expMarker$MFP!="IE")])
  matrix(c(a,b,c,d),nrow = 2)
  res = fisher.test(matrix(c(a,b,c,d),nrow = 2),alternative = "greater")
  res$p.value
  ## no AUC for MFP, then return AUC as 0.5
  result=c(res$p.value,0.5)
  names(result) = c("p_value","AUC")
  return(result)
  
}


{
# ###########
# ####SIC  MCP-Counter
# # install.packages(c("devtools","curl")) ##Installs devtools and the MCPcounter dependancy 'curl'
# # library(devtools)
# # install_github("ebecht/MCPcounter",ref="master", subdir="Source")
# ??MCPcounter.estimate
# library(MCPcounter)
# ## expressionMatrix of SRAC, focus on UPS,DDLPS,LMS
# library(readr)
# TCGA_SARC_GDC_phenotype <- read_delim("D:/TCGA _GDC/TCGA-SARC.GDC_phenotype.tsv", 
#                                       delim = "\t", escape_double = FALSE, 
#                                       trim_ws = TRUE)
# TCGA_SARC = TCGA_SARC_GDC_phenotype %>% filter(primary_diagnosis.diagnoses %in% c("Leiomyosarcoma, NOS","Dedifferentiated liposarcoma","Undifferentiated sarcoma"))
# expressionMatrix = expressionMatrix[,intersect(TCGA_SARC$submitter_id.samples,colnames(expressionMatrix))]
# expressionMatrix[1:3,1:3]
# expressionMatrix = expressionMatrix[,intersect(TCGA_SARC$submitter_id.samples,colnames(expressionMatrix))]
# MCPcounter_score = MCPcounter.estimate(log2(expressionMatrix+1),featuresType="HUGO_symbols")
# MCPcounter_score[1:3,1:3]
# ##gene signatures for functional orientation
# immunosuppression = c("CXCL12", "TGFB1", "TGFB3", "LGALS1")
# T_cell_activation = c("CXCL9", "CXCL10", "CXCL16", "IFNG", "IL15")
# T_cell_survival = c("CD70","CD27")
# regulatory_T_cells = c("FOXP3","TNFRSF18")
# major_histocompatibility_complex_class_I=c("HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G" , "B2M")
# meyeloid_cells_chemotactism = c("CCL2")
# TLS = c("CXCL13")
# ##immune checkpoints related genes
# ICs = c("PDCD1","CD274","PDCD1LG2","CTLA4","HAVCR2","LAG3")
# ##compute the total score matrix and cluster
# immunosuppression_score = colMeans(expressionMatrix[immunosuppression,])
# T_cell_activation_score = colMeans(expressionMatrix[T_cell_activation,])
# T_cell_survival_score = colMeans(expressionMatrix[T_cell_survival,])
# regulatory_T_cells_score = colMeans(expressionMatrix[regulatory_T_cells,])
# major_histocompatibility_complex_class_I_score = colMeans(expressionMatrix[major_histocompatibility_complex_class_I,])
# meyeloid_cells_chemotactism_score = expressionMatrix[meyeloid_cells_chemotactism,]
# TLS_score = expressionMatrix[TLS,]
# ICs_score = expressionMatrix[ICs,]
# scoreMatrix = rbind(MCPcounter_score[-10,],immunosuppression_score,T_cell_activation_score,
#                     T_cell_survival_score,regulatory_T_cells_score,major_histocompatibility_complex_class_I_score,
#                     TLS_score,ICs_score)
# 
# library('pheatmap')
# pheatmap(MCPcounter_score[-10,],cluster_rows =F,show_colnames = F,
#          cutree_cols = 5,scale = "row",breaks = seq(-4,4,by=0.1),legend_breaks = seq(-4,4,by=1)) -> res
# pheatmap(scoreMatrix,cluster_rows =F,show_colnames = F,
#          cutree_cols = 5,scale = "row",breaks = seq(-4,4,by=0.1),legend_breaks = seq(-4,4,by=1)) -> res
# 
# library(gplots)
# heatmap.2(MCPcounter_score[-10,],dendrogram ="column",scale = "row")
# 
# 
# 
# 
# 
# 
# ###########
# ##prognostic  32-gene signature for gastric cancer
# signatures32 = c("FHL2","PML","BRCA1","WT1","AREG","TP63","ESR1","BEST1","ACTA2",
#                  "HIPK2","IGSF9","ASCC2","JUN","PPP2R5A","SMAD3","CREBBP","EP300",
#                  "DDX5","TP53","HSF1","TGS1","PAWR","FAM96A","WTAP","PCNA",
#                  "GNL3","WRN","SMARCA4","NCOA6","RPA1","MSH6","PARP1")
# intersec_gene = intersect(signatures32,rownames(data$TPM))
# library(NMF)
# # res = nmf(data$TPM[intersec_gene,],2:7,nrun=10)
# # plot(res)
# res <- nmf(data$TPM[intersec_gene,],4,nrun=1000)
# coefmap(res)
# xx = consensusmap(res,hclustfun="average")
# predict(res)
}

