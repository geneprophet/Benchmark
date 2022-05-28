##construct the mysql table and save to the server
library(RMySQL)
# mydb = dbConnect(MySQL(), user='kkk', password='kkk@123456', 
#                  dbname='kanghe_db', host='192.168.164.82',port=33060)
mydb = dbConnect(MySQL(), user='root', password='kang123456', 
                 dbname='kanghe_db', host='127.0.0.1',port=3306)
dbListTables(mydb)

library(readxl)
dataset <- read_excel("statistical.xlsx", 
                      sheet = "Dataset")
dbWriteTable(mydb,"dataset",dataset, overwrite = F, append= T,
             row.names = F)


biomarker <- read_excel("statistical.xlsx", 
                        sheet = "Biomarker")
biomarker$category[which(biomarker$category==1)] = "simple operation"
biomarker$category[which(biomarker$category==2)] = "ssGSEA based"
biomarker$category[which(biomarker$category==3)] = "deconvolution based"
load("marker/all_biomarker_genes.Rdata")
biomarker2gene = c()
for (i in unique(biomarker_genes$biomarker)) {
  a = biomarker_genes$gene[which(biomarker_genes$biomarker==i)]
  biomarker2gene = rbind(biomarker2gene,c(i,paste0(a,collapse = ', ')))
}
tmp = as.data.frame(biomarker2gene)
colnames(tmp) = c("name","genes")
all_biomarkers = merge.data.frame(biomarker,tmp,by = "name",all.x = T,sort = F)

dbWriteTable(mydb,"biomarker",all_biomarkers, overwrite = F, append= T,
             row.names = F)




kk_extract = function(data,name){
  cnames = c("Sample","Resp_NoResp","OS","OS_CNSR","PFS","PFS_CNSR","Pre_On","Treatment","Sex","Age")
  Clinical = data$Clinical[,intersect(colnames(data$Clinical),cnames)]
  if(is.element("Age",colnames(Clinical))){
    Clinical$Age = as.numeric(Clinical$Age)
  }
  if(is.element("OS",colnames(Clinical))){
    Clinical$OS = as.numeric(Clinical$OS)
  }
  if(is.element("OS_CNSR",colnames(Clinical))){
    Clinical$OS_CNSR = as.numeric(Clinical$OS_CNSR)
  }
  if(is.element("PFS",colnames(Clinical))){
    Clinical$PFS = as.numeric(Clinical$PFS)
  }
  if(is.element("PFS_CNSR",colnames(Clinical))){
    Clinical$PFS_CNSR = as.numeric(Clinical$PFS_CNSR)
  }
  Clinical$Dataset = name
  return(Clinical)
}
{
load("./ICB_data/Braun et al/Braun_2020.Rdata")
data=Braun_data; name="Braun_2020"
clinical1 = kk_extract(data,name)
load("./ICB_data/Du et al/MGH_PRE_2021.RData")
data=MGH_PRE_data; name="MGH_PRE_2021"
clinical2 = kk_extract(data,name)
load("./ICB_data/Du et al/MGH_ON_2021.RData")
data=MGH_ON_data; name="MGH_ON_2021"
clinical3 = kk_extract(data,name)
load("./ICB_data/Gide et al/Gide_2019.RData")
data=Gide_data; name="Gide_2019"
clinical4 = kk_extract(data,name)
load("./ICB_data/Hugo et al/Hugo_2016.RData")
data=Hugo_data; name="Hugo_2016"
clinical5 = kk_extract(data,name)
load("./ICB_data/Jung et al/Jung_2019.Rdata")
data=Jung_data; name="Jung_2019"
clinical6 = kk_extract(data,name)
load("./ICB_data/Kim et al/Kim_2018.RData")
data=Kim_data; name="Kim_2018"
clinical7 = kk_extract(data,name)
load("./ICB_data/Lee et al/Lee_2020.RData")
data=Lee_data; name="Lee_2020"
clinical8 = kk_extract(data,name)
load("./ICB_data/Liu et al/Liu_2019.Rdata")
data=Liu_data; name="Liu_2019"
clinical9 = kk_extract(data,name)
load("./ICB_data/Mariathasan et al/Mariathasan_2018.Rdata")
data=Mariathasan_data; name="Mariathasan_2018"
clinical10 = kk_extract(data,name)
load("./ICB_data/Miao et al/Miao_2018.Rdata")
data=Miao_data; name="Miao_2018"
clinical11 = kk_extract(data,name)
load("./ICB_data/Motzer et al/Motzer_2020.Rdata")
data=Motzer_data; name="Motzer_2020"
clinical12 = kk_extract(data,name)
load("./ICB_data/Nathanson et al/Nathanson_2017.Rdata")
data=Nathanson_data; name="Nathanson_2017"
clinical13 = kk_extract(data,name)
load("./ICB_data/Riaz et al/Riaz_2017.RData")
data=Riaz_data; name="Riaz_2017"
clinical14 = kk_extract(data,name)
load("./ICB_data/Snyder et al/Snyder_2017.Rdata")
data=Snyder_data; name="Snyder_2017"
clinical15 = kk_extract(data,name)
load("./ICB_data/VanAllen et al/VanAllen_2015.RData")
data=VanAllen_data; name="VanAllen_2015"
clinical16 = kk_extract(data,name)

load("ICB_data/Gide et al/Gide_PRE_2019.Rdata")
data=Gide_PRE_data; name="Gide_PRE_2019"
clinical17 = kk_extract(data,name)
load("ICB_data/Gide et al/Gide_ON_2019.Rdata")
data=Gide_ON_data; name="Gide_ON_2019"
clinical18 = kk_extract(data,name)
load("ICB_data/Lee et al/Lee_PRE_2020.Rdata")
data=Lee_PRE_data; name="Lee_PRE_2020"
clinical19 = kk_extract(data,name)
load("ICB_data/Lee et al/Lee_ON_2020.Rdata")
data=Lee_ON_data; name="Lee_ON_2020"
clinical20 = kk_extract(data,name)
load("ICB_data/Riaz et al/Riaz_PRE_2017.Rdata")
data=Riaz_PRE_data; name="Riaz_PRE_2017"
clinical21 = kk_extract(data,name)
load("ICB_data/Riaz et al/Riaz_ON_2017.Rdata")
data=Riaz_ON_data; name="Riaz_ON_2017"
clinical22 = kk_extract(data,name)
load("ICB_data/Gide et al/Gide_MONO_2019.Rdata")
data = Gide_MONO_data; name = "Gide_MONO_2019"
clinical23 = kk_extract(data,name)
load("ICB_data/Gide et al/Gide_COMBINE_data.Rdata")
data = Gide_COMBINE_data; name = "Gide_COMBINE_2019"
clinical24 = kk_extract(data,name)
load("ICB_data/Riaz et al/Riaz_NAIVE_2017.Rdata")
data = Riaz_NAIVE_data; name = "Riaz_NAIVE_2017"
clinical25 = kk_extract(data,name)
load("ICB_data/Riaz et al/Riaz_EXPOSURE_2017.Rdata")
data = Riaz_EXPOSURE_data; name = "Riaz_EXPOSURE_2017"
clinical26 = kk_extract(data,name)
load("ICB_data/Liu et al/Liu_NAIVE_2019.Rdata")
data = Liu_NAIVE_data; name = "Liu_NAIVE_2019"
clinical27 = kk_extract(data,name)
load("ICB_data/Liu et al/Liu_EXPOSURE_2019.Rdata")
data = Liu_EXPOSURE_data; name = "Liu_EXPOSURE_2019"
clinical28 = kk_extract(data,name)
}
library(dplyr)
clinical = dplyr::bind_rows(clinical1,clinical2,clinical3,clinical4,clinical5,clinical6,clinical7,
                            clinical8,clinical9,clinical10,clinical11,clinical12,clinical13,clinical14,
                            clinical15,clinical16,clinical17,clinical18,clinical19,clinical20,clinical21,
                            clinical22,clinical23,clinical24,clinical25,clinical26,clinical27,clinical28)
colnames(clinical) = tolower(colnames(clinical))
dbWriteTable(mydb,"clinical",clinical, overwrite = F, append= T,
             row.names = F)


####
{
load("Results/benchmark_results.Rdata")
res1 = res
load("Results/average2ssgsea_results2.Rdata")
res2 = res
load("Results/Pre_On_benchmark_results.Rdata")
res3 = res
load("Results/treatment_benchmark_results.Rdata")
res4 = res

res = rbind(res1,res2,res3,res4)
##统一biomarker和dataset的名字
res$Biomarker = as.character(res$Biomarker)
res$Dataset = as.character(res$Dataset)
res$Dataset[which(res$Dataset=="ccRCC_124_Braun")] = "Braun_2020"
res$Dataset[which(res$Dataset=="Melanoma_19_MGH_PRE")] = "MGH_PRE_2021"
res$Dataset[which(res$Dataset=="Melanoma_31_MGH_ON")] = "MGH_ON_2021"
res$Dataset[which(res$Dataset=="Melanoma_90_Gide")] = "Gide_2019"
res$Dataset[which(res$Dataset=="Melanoma_26_Hugo")] = "Hugo_2016"
res$Dataset[which(res$Dataset=="NSCLC_27_Jung")] = "Jung_2019"
res$Dataset[which(res$Dataset=="Gastric-Cancer_45_Kim")] = "Kim_2018"
res$Dataset[which(res$Dataset=="Melanoma_79_Lee")] = "Lee_2020"
res$Dataset[which(res$Dataset=="Melanoma_121_Liu")] = "Liu_2019"
res$Dataset[which(res$Dataset=="Urothelial-Cancer_348_Mariathasan")] = "Mariathasan_2018"
res$Dataset[which(res$Dataset=="ccRCC_33_Miao")] = "Miao_2018"
res$Dataset[which(res$Dataset=="NSCLC_354_Motzer")] = "Motzer_2020"
res$Dataset[which(res$Dataset=="Melanoma_24_Nathanson")] = "Nathanson_2017"
res$Dataset[which(res$Dataset=="Melanoma_103_Riaz")] = "Riaz_2017"
res$Dataset[which(res$Dataset=="Urothelial-Cancer_26_Snyder")] = "Snyder_2017"
res$Dataset[which(res$Dataset=="Melanoma_42_VanAllen")] = "VanAllen_2015"

res$Dataset[which(res$Dataset=="Gide_Melanoma_90")] = "Gide_2019"
res$Dataset[which(res$Dataset=="Gide_PRE_Melanoma_72")] = "Gide_PRE_2019"
res$Dataset[which(res$Dataset=="Gide_ON_Melanoma_18")] = "Gide_ON_2019"
res$Dataset[which(res$Dataset=="Riaz_Melanoma_103")] = "Riaz_2017"
res$Dataset[which(res$Dataset=="Riaz_PRE_Melanoma_49")] = "Riaz_PRE_2017"
res$Dataset[which(res$Dataset=="Riaz_ON_Melanoma_54")] = "Riaz_ON_2017"
res$Dataset[which(res$Dataset=="Lee_Melanoma_79")] = "Lee_2020"
res$Dataset[which(res$Dataset=="Lee_PRE_Melanoma_44")] = "Lee_PRE_2020"
res$Dataset[which(res$Dataset=="Lee_ON_Melanoma_35")] = "Lee_ON_2020"

res$Dataset[which(res$Dataset=="MGH_PRE_Melanoma_19")] = "MGH_PRE_2021"
res$Dataset[which(res$Dataset=="MGH_ON_Melanoma_31")] = "MGH_ON_2021"
res$Dataset[which(res$Dataset=="Gide_MONO_Melanoma_50")] = "Gide_MONO_2019"
res$Dataset[which(res$Dataset=="Gide_COMBINE_Melanoma_40")] = "Gide_COMBINE_2019"
res$Dataset[which(res$Dataset=="Riaz_NAIVE_Melanoma_45")] = "Riaz_NAIVE_2017"
res$Dataset[which(res$Dataset=="Riaz_EXPOSURE_Melanoma_58")] = "Riaz_EXPOSURE_2017"
res$Dataset[which(res$Dataset=="Liu_Melanoma_121")] = "Liu_2019"
res$Dataset[which(res$Dataset=="Liu_NAIVE_Melanoma_74")] = "Liu_NAIVE_2019"
res$Dataset[which(res$Dataset=="Liu_EXPOSURE_Melanoma_47")] = "Liu_EXPOSURE_2019"
}
res = res[which(!duplicated.data.frame(res)),]
dim(res)
colnames(res) = c("pvalue","auc","dataset","biomarker")
res$wilcoxpath = paste0("Results/Results_wilcox_test/",res$dataset,"/",res$biomarker,".png")
res$aucpath = paste0("Results/Results_AUC/",res$dataset,"/",res$biomarker,".png")

dbWriteTable(mydb,"benchmark",res, overwrite = F, append= T,
             row.names = F)

####
load("Results/OS_results.Rdata")
load("Results/PFS_results.Rdata")
all_os = OS_results
all_pfs = PFS_results
load("Results/OS_results_subdatasets.Rdata")
load("Results/PFS_results_subdatasets.Rdata")
OS = rbind.data.frame(all_os,OS_results)
PFS = rbind.data.frame(all_pfs,PFS_results)

colnames(OS) = c("hr","hr_upper","hr_lower","pvalue","dataset","biomarker")
OS$kmpath = paste0("Results/Results_survival_analysis/OS/",OS$dataset,"/",OS$biomarker,".png")
OS$forestpath = paste0("Results/Results_survival_analysis/OS/forest/",OS$dataset,"/",OS$biomarker,".png")

colnames(PFS) = c("hr","hr_upper","hr_lower","pvalue","dataset","biomarker")
PFS$kmpath = paste0("Results/Results_survival_analysis/PFS/",PFS$dataset,"/",PFS$biomarker,".png")
PFS$forestpath = paste0("Results/Results_survival_analysis/PFS/forest/",PFS$dataset,"/",PFS$biomarker,".png")

dbWriteTable(mydb,"os",OS, overwrite = F, append= T,
             row.names = F)

dbWriteTable(mydb,"pfs",PFS, overwrite = F, append= T,
             row.names = F)


################biomarker_genes
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
gene = select(org.Hs.eg.db, keys=keys(org.Hs.eg.db), columns=c("SYMBOL","GENETYPE","ENSEMBL","ENTREZID","ALIAS","GENENAME","MAP"))
##integrate the alias to one row
new_genes = c()
for (i in gene$ENTREZID) {
  sub = gene[which(gene$ENTREZID==i),]
  row = sub[1,]
  row$ALIAS = paste0(sub$ALIAS,collapse = ', ')
  new_genes = rbind(new_genes,row)
}
save(new_genes,file = "Results/all_genes.Rdata")

load("Results/all_genes.Rdata")
load("marker/all_biomarker_genes.Rdata")
all_biomarker_genes = merge.data.frame(biomarker_genes,new_genes,all.x = T,by.x = "gene",by.y = "SYMBOL")
all_biomarker_genes = all_biomarker_genes[which(!duplicated.data.frame(all_biomarker_genes[,c(1,2)])),]
# non = all_biomarker_genes[which(is.na(all_biomarker_genes$ENTREZID)),]

colnames(all_biomarker_genes) = c("symbol","biomarker","entrezid","genetype","ensembl","alias","genename","map")

dbWriteTable(mydb,"biomarker_genes",all_biomarker_genes, overwrite = F, append= T,
             row.names = F)

#########TCGA
load("Results/TCGA_OS_results.Rdata")
colnames(OS_results) = c("hr","hr_upper","hr_lower","pvalue","dataset","biomarker")
OS_results$kmpath = paste0("Results/Results_survival_analysis/OS/",OS_results$dataset,"/",OS_results$biomarker,".png")
OS_results$forestpath = paste0("Results/Results_survival_analysis/OS/forest/",OS_results$dataset,"/",OS_results$biomarker,".png")

dbWriteTable(mydb,"tcga",OS_results, overwrite = F, append= T,
             row.names = F)


###landscape
{
  library(tidyr)
  library(dplyr)
  construct_landscape = function(data,name){
    landscape = data$Landscape
    landscape = landscape %>% dplyr::select(-c("MFP","Ecotype"))
    
    a = apply(landscape[,-1], 2, scale)
    b = as.data.frame(a)
    b$sample = landscape$Sample
    
    tmp = gather(b,biomarker,value,PD_L1,PD_1,PD_L2,CX3CL1,CTLA4,CYT_score,HLA_DRA,IFN_gamma,Expanded_immune_gene_signature,
                 T_cell_inflamed_GEP_score,Immunophenoscore,IMPRES_score,CRMA_score,The_immune_resistance_program,EMT_Stroma_core_signature,
                 F_TBRS,TMEscore,RiskScore,TLS_score,CXCL9,MPS_score,Renal_101_Immuno_signature,HRH1,TIDE,IIS_score,TIS_score,
                 APM_score,IPRES_score,C_ECM_score,IMS_score,PASS_PRE,PASS_ON,MIAS_score,CD8T_xCell,CD8T_MCPcounter,CD8T_CIBERSORTx,
                 Immunoscore_CIBERSORTx,IFN_gamma_ssGSEA,Expanded_immune_gene_ssGSEA,T_cell_inflamed_GEP_ssGSEA,CRMA_ssGSEA,
                 EMT_Stroma_core_ssGSEA,F_TBRS_ssGSEA,RiskScore_ssGSEA,TLS_score_ssGSEA,Renal_101_Immuno_ssGSEA)
    tmp$dataset = name
    colnames(tmp) = tolower(colnames(tmp))
    return(tmp)
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
  
  ###sub-datasets
  load("./ICB_data/Gide et al/Gide_PRE_2019.Rdata")
  data17=Gide_PRE_data
  name17="Gide_PRE_2019"
  
  load("./ICB_data/Gide et al/Gide_ON_2019.Rdata")
  data18=Gide_ON_data
  name18="Gide_ON_2019"
  
  load("./ICB_data/Riaz et al/Riaz_PRE_2017.Rdata")
  data19=Riaz_PRE_data
  name19="Riaz_PRE_2017"
  
  load("./ICB_data/Riaz et al/Riaz_ON_2017.Rdata")
  data20=Riaz_ON_data
  name20="Riaz_ON_2017"
  
  load("./ICB_data/Lee et al/Lee_PRE_2020.Rdata")
  data21=Lee_PRE_data
  name21="Lee_PRE_2020"
  
  load("./ICB_data/Lee et al/Lee_ON_2020.Rdata")
  data22=Lee_ON_data
  name22="Lee_ON_2020"
  
  load("./ICB_data/Gide et al/Gide_MONO_2019.Rdata")
  data23=Gide_MONO_data
  name23="Gide_MONO_2019"
  
  load("./ICB_data/Gide et al/Gide_COMBINE_2019.Rdata")
  data24=Gide_COMBINE_data
  name24="Gide_COMBINE_2019"
  
  load("./ICB_data/Riaz et al/Riaz_NAIVE_2017.Rdata")
  data25=Riaz_NAIVE_data
  name25="Riaz_NAIVE_2017"
  
  load("./ICB_data/Riaz et al/Riaz_EXPOSURE_2017.Rdata")
  data26=Riaz_EXPOSURE_data
  name26="Riaz_EXPOSURE_2017"
  
  load("./ICB_data/Liu et al/Liu_NAIVE_2019.Rdata")
  data27=Liu_NAIVE_data
  name27="Liu_NAIVE_2019"
  
  load("./ICB_data/Liu et al/Liu_EXPOSURE_2019.Rdata")
  data28=Liu_EXPOSURE_data
  name28="Liu_EXPOSURE_2019"
  
  all_landscape = rbind.data.frame(construct_landscape(data1,name1),
                                   construct_landscape(data2,name2),
                                   construct_landscape(data3,name3),
                                   construct_landscape(data4,name4),
                                   construct_landscape(data5,name5),
                                   construct_landscape(data6,name6),
                                   construct_landscape(data7,name7),
                                   construct_landscape(data8,name8),
                                   construct_landscape(data9,name9),
                                   construct_landscape(data10,name10),
                                   construct_landscape(data11,name11),
                                   construct_landscape(data12,name12),
                                   construct_landscape(data13,name13),
                                   construct_landscape(data14,name14),
                                   construct_landscape(data15,name15),
                                   construct_landscape(data16,name16),
                                   construct_landscape(data17,name17),
                                   construct_landscape(data18,name18),
                                   construct_landscape(data19,name19),
                                   construct_landscape(data20,name20),
                                   construct_landscape(data21,name21),
                                   construct_landscape(data22,name22),
                                   construct_landscape(data23,name23),
                                   construct_landscape(data24,name24),
                                   construct_landscape(data25,name25),
                                   construct_landscape(data26,name26),
                                   construct_landscape(data27,name27),
                                   construct_landscape(data28,name28))
}


dbWriteTable(mydb,"landscape",all_landscape, overwrite = F, append= T,
             row.names = F)


#####cibersort
{
  library(tidyr)
  construct_cibersort = function(data,name){
    file_path = paste0("CIBERSORTx_output/CIBERSORTx_",name,"_Results.txt")
    library(readr)
    CIBERSORTx_result <- read_delim(file_path,
                                    delim = "\t",
                                    escape_double = FALSE, 
                                    trim_ws = TRUE)
    
    tmp = gather(CIBERSORTx_result[,-c(24:26)],celltype,value,`B cells naive`,`B cells memory`,`Plasma cells`,`T cells CD8`,
                 `T cells CD4 naive`,`T cells CD4 memory resting`,`T cells CD4 memory activated`,`T cells follicular helper`,
                 `T cells regulatory (Tregs)`,`T cells gamma delta`,`NK cells resting`,`NK cells activated`,`Monocytes`,`Macrophages M0`,
                 `Macrophages M1`,`Macrophages M2`,`Dendritic cells resting`,`Dendritic cells activated`,`Mast cells resting`,
                 `Mast cells activated`,`Eosinophils`,`Neutrophils`)
    tmp$dataset = name
    colnames(tmp) = c("sample","celltype","value","dataset")
    return(tmp)
  }
  
  all_cibersort = rbind.data.frame(construct_cibersort(data1,name1),
                                   construct_cibersort(data2,name2),
                                   construct_cibersort(data3,name3),
                                   construct_cibersort(data4,name4),
                                   construct_cibersort(data5,name5),
                                   construct_cibersort(data6,name6),
                                   construct_cibersort(data7,name7),
                                   construct_cibersort(data8,name8),
                                   construct_cibersort(data9,name9),
                                   construct_cibersort(data10,name10),
                                   construct_cibersort(data11,name11),
                                   construct_cibersort(data12,name12),
                                   construct_cibersort(data13,name13),
                                   construct_cibersort(data14,name14),
                                   construct_cibersort(data15,name15),
                                   construct_cibersort(data16,name16),
                                   construct_cibersort(data17,name17),
                                   construct_cibersort(data18,name18),
                                   construct_cibersort(data19,name19),
                                   construct_cibersort(data20,name20),
                                   construct_cibersort(data21,name21),
                                   construct_cibersort(data22,name22),
                                   construct_cibersort(data23,name23),
                                   construct_cibersort(data24,name24),
                                   construct_cibersort(data25,name25),
                                   construct_cibersort(data26,name26),
                                   construct_cibersort(data27,name27),
                                   construct_cibersort(data28,name28))
}

dbWriteTable(mydb,"cibersort",all_cibersort, overwrite = F, append= T,
             row.names = F)

###xCell
{
  library(tidyr)
  library(xCell)
  construct_xCell = function(data,name){
    xCell = xCellAnalysis(data$TPM)
    tmp1 = t(xCell)
    tmp2 = as.data.frame(tmp1[,-c(65:67)])
    tmp2$sample = rownames(tmp2)
    tmp3 = gather(tmp2,celltype,value,colnames(tmp2)[-65])
    tmp3$dataset = name
    return(tmp3)
  }
  all_xcell = rbind.data.frame(construct_xCell(data1,name1),
                               construct_xCell(data2,name2),
                               construct_xCell(data3,name3),
                               construct_xCell(data4,name4),
                               construct_xCell(data5,name5),
                               construct_xCell(data6,name6),
                               construct_xCell(data7,name7),
                               construct_xCell(data8,name8),
                               construct_xCell(data9,name9),
                               construct_xCell(data10,name10),
                               construct_xCell(data11,name11),
                               construct_xCell(data12,name12),
                               construct_xCell(data13,name13),
                               construct_xCell(data14,name14),
                               construct_xCell(data15,name15),
                               construct_xCell(data16,name16),
                               construct_xCell(data17,name17),
                               construct_xCell(data18,name18),
                               construct_xCell(data19,name19),
                               construct_xCell(data20,name20),
                               construct_xCell(data21,name21),
                               construct_xCell(data22,name22),
                               construct_xCell(data23,name23),
                               construct_xCell(data24,name24),
                               construct_xCell(data25,name25),
                               construct_xCell(data26,name26),
                               construct_xCell(data27,name27),
                               construct_xCell(data28,name28))
  
}

dbWriteTable(mydb,"xcell",all_xcell, overwrite = F, append= T,
             row.names = F)

##MCPCounter
{
  library(tidyr)
  library(MCPcounter)
  construct_mcpcounter = function(data,name){
    
    MCPcounter_score = MCPcounter.estimate(data$TPM,featuresType="HUGO_symbols")
    tmp1 = as.data.frame(t(MCPcounter_score))
    tmp1$sample = rownames(tmp1)
    tmp2 = gather(tmp1,celltype,value,colnames(tmp1)[-11])
    tmp2$dataset = name
    return(tmp2)
  }
  all_mcpcounter = rbind.data.frame(construct_mcpcounter(data1,name1),
                                    construct_mcpcounter(data2,name2),
                                    construct_mcpcounter(data3,name3),
                                    construct_mcpcounter(data4,name4),
                                    construct_mcpcounter(data5,name5),
                                    construct_mcpcounter(data6,name6),
                                    construct_mcpcounter(data7,name7),
                                    construct_mcpcounter(data8,name8),
                                    construct_mcpcounter(data9,name9),
                                    construct_mcpcounter(data10,name10),
                                    construct_mcpcounter(data11,name11),
                                    construct_mcpcounter(data12,name12),
                                    construct_mcpcounter(data13,name13),
                                    construct_mcpcounter(data14,name14),
                                    construct_mcpcounter(data15,name15),
                                    construct_mcpcounter(data16,name16),
                                    construct_mcpcounter(data17,name17),
                                    construct_mcpcounter(data18,name18),
                                    construct_mcpcounter(data19,name19),
                                    construct_mcpcounter(data20,name20),
                                    construct_mcpcounter(data21,name21),
                                    construct_mcpcounter(data22,name22),
                                    construct_mcpcounter(data23,name23),
                                    construct_mcpcounter(data24,name24),
                                    construct_mcpcounter(data25,name25),
                                    construct_mcpcounter(data26,name26),
                                    construct_mcpcounter(data27,name27),
                                    construct_mcpcounter(data28,name28))
}

dbWriteTable(mydb,"mcpcounter",all_mcpcounter, overwrite = F, append= T,
             row.names = F)



dbDisconnect(mydb)




