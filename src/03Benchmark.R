##Benchnark_results

benchmark_result = function(data,name){
  print(paste0("Begin benchmark of dataset: ",name))
  dir.create(paste0("Results/Results_wilcox_test/",name))
  dir.create(paste0("Results/Results_AUC/",name))
  source("./src/marker_gene_expression.R")
  marker_benchmark_results = rbind(PD_L1_marker(data,name),
                                   PD_1_marker(data,name),
                                   PD_L2_marker(data,name),
                                   CX3CL1_marker(data,name),
                                   CTLA4_marker(data,name),
                                   CYT_marker(data,name),
                                   HLA_DRA_marker(data,name),
                                   IFN_gamma_marker(data,name),
                                   Expanded_immune_gene_signature_marker(data,name),
                                   T_cell_inflamed_GEP_score(data,name),
                                   Immunophenoscore(data,name),
                                   IMPRES(data,name),
                                   CRMA_score(data,name),
                                   Immune_resistance_program_score(data,name),
                                   EMT_Stroma_core_marker(data,name),
                                   F_TBRS_score(data,name),
                                   TMEscore(data,name),
                                   Risk_score(data,name),
                                   TLS_score(data,name),
                                   CXCL9_marker(data,name),
                                   MPS_score(data,name),
                                   Renal_101_Immuno_score(data,name),
                                   HRH1_marker(data,name),
                                   TIDE(data,name))
  rownames(marker_benchmark_results) = c("PD_L1","PD_1","PD_L2","CX3CL1","CTLA4","CYT_score","HLA_DRA",
                                         "IFN_gamma","Expanded_immune_gene_signature","T_cell_inflamed_GEP_score",
                                         "Immunophenoscore","IMPRES_score","CRMA_score",
                                         "The_immune_resistance_program","EMT_Stroma_core_signature","F_TBRS","TMEscore",
                                         "RiskScore","TLS_score","CXCL9","MPS_score","Renal_101_Immuno_signature","HRH1","TIDE")
  
  source("src/ssGSEA_based.R")
  ssGSEA_based_benchmark = rbind(TIS_score(data,name),
                                 IIS_score(data,name),
                                 APM_score(data,name),
                                 IPRES_score(data,name),
                                 C_ECM_score(data,name),
                                 IMS_score(data,name),
                                 PASS_PRE(data,name),
                                 PASS_ON(data,name),
                                 MIAS_score(data,name))
  rownames(ssGSEA_based_benchmark) = c("TIS_score","IIS_score","APM_score","IPRES_score","C_ECM_score",
                                       "IMS_score","PASS_PRE","PASS_ON","MIAS_score")
  source("src/deconvolution.R")
  deconvolution_benchmark_results = rbind(CD8T_CIBERSORTx(data,name),
                                          CD8T_MCPcounter(data,name),
                                          CD8T_xCell(data,name),
                                          Immunoscore_CIBERSORTx(data,name),
                                          MFP(data,name),
                                          EcoTyper(data,name))
  rownames(deconvolution_benchmark_results)  = c("CD8T_CIBERSORTx","CD8T_MCPcounter","CD8T_xCell",
                                                 "Immunoscore_CIBERSORTx","MFP","Ecotype")
  
  
  Dataset=rep(name,39)
  result = as.data.frame(cbind(rbind(marker_benchmark_results,ssGSEA_based_benchmark,deconvolution_benchmark_results),Dataset))
  result$Biomarker = c(rownames(marker_benchmark_results),rownames(ssGSEA_based_benchmark),rownames(deconvolution_benchmark_results))
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

res = rbind(benchmark_result(data1,name1),
            benchmark_result(data2,name2),
            benchmark_result(data3,name3),
            benchmark_result(data4,name4),
            benchmark_result(data5,name5),
            benchmark_result(data6,name6),
            benchmark_result(data7,name7),
            benchmark_result(data8,name8),
            benchmark_result(data9,name9),
            benchmark_result(data10,name10),
            benchmark_result(data11,name11),
            benchmark_result(data12,name12),
            benchmark_result(data13,name13),
            benchmark_result(data14,name14),
            benchmark_result(data15,name15),
            benchmark_result(data16,name16))
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
res$Dataset = as.factor(res$Dataset)
res$Biomarker = as.factor(res$Biomarker)
res$p_value = as.numeric(res$p_value)
res$AUC = as.numeric(res$AUC)


source("./src/visualization.R")
p = dotplot(res)
ggsave(p,filename = "./figures/benchmark.pdf",width = 16,height = 12)

save(list = ls(),file = "Results/benchmark_results.Rdata")


