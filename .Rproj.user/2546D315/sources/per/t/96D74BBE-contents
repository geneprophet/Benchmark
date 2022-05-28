##TCGA
# source("./src/utils.R")
# TCGA = c("ACC","CHOL","BLCA","BRCA","CESC","COAD","UCEC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LIHC",
#          "LGG","LUAD","LUSC","SKCM","MESO","UVM","OV","PAAD","PCPG","PRAD","READ","SARC","STAD","TGCT","THYM", "THCA","UCS")
# ## generate TIDE input files
# tide_input = function(data,name){
#   data = log2(data$TPM+1)
#   data = scale(data)
#   Gene = rownames(data)
#   data = cbind(Gene,data)
#   write.table(data,file = paste0("./TIDE_input/",name,".txt"),quote = F,sep = "\t",row.names = F,col.names = T)
# }
# ## generate CIBERSORTx input files
# CIBERSORTx_input = function(data,name){
#   data = as.data.frame(data$TPM)
#   Gene = rownames(data)
#   data = cbind(Gene,data)
#   write.table(data,file = paste0("./CIBERSORTx_input/",name,".txt"),quote = F,sep = "\t",row.names = F,col.names = T)
# }
# 
# 
# for(CANCER in TCGA){
#   fpkm_file_path = paste0("D:/TCGA/TCGA-",CANCER,".htseq_fpkm.tsv")
#   survival_file_path = paste0("D:/TCGA/TCGA-",CANCER,".survival.tsv")
#   library(readr)
#   htseq_fpkm <- read_delim(fpkm_file_path, 
#                            delim = "\t", escape_double = FALSE, 
#                            trim_ws = TRUE)
#   survival <- read_delim(survival_file_path, 
#                          delim = "\t", escape_double = FALSE, 
#                          trim_ws = TRUE)
#   ##log2(fpkm+1) to fpkm
#   2^htseq_fpkm[,-1] -1 -> fpkm
#   ##fpkm to tpm
#   apply(fpkm, 2, fpkmToTpm) -> tpm
#   rownames(tpm) = htseq_fpkm$Ensembl_ID
#   #colSums(tpm)
#   tpm[gtf_data$gene_id,] -> pcg
#   rownames(pcg) = gtf_data$gene_name
#   
#   dataset = list()
#   dataset[["TPM"]] = pcg
#   colnames(survival) = c("Sample","OS_CNSR","_PATIENT","OS")
#   dataset[["Clinical"]] = survival
#   assign(paste0("TCGA_",CANCER,"_data"),dataset)
#   save(list = paste0("TCGA_",CANCER,"_data"),file = paste0("TCGA/TCGA_",CANCER,"_data.Rdata"))
#   tide_input(get(paste0("TCGA_",CANCER,"_data")),paste0("TCGA_",CANCER,"_data"))
#   CIBERSORTx_input(get(paste0("TCGA_",CANCER,"_data")),paste0("TCGA_",CANCER,"_data"))
#   rm(list = paste0("TCGA_",CANCER,"_data"))
# }
# 
# 
# # library(readr)
# # TCGA_BRCA_data <- read_delim("CIBERSORTx_input/TCGA_BRCA_data.txt", 
# #                              delim = "\t", escape_double = FALSE, 
# #                              trim_ws = TRUE)
# # TCGA_BRCA_data_1 = TCGA_BRCA_data[,c(1:600)]
# # TCGA_BRCA_data_2 = TCGA_BRCA_data[,c(1,601:1218)]
# # write.table(TCGA_BRCA_data_1,file = "CIBERSORTx_input/TCGA_BRCA_data1.txt",sep = "\t",quote = F,row.names = F)
# # write.table(TCGA_BRCA_data_2,file = "CIBERSORTx_input/TCGA_BRCA_data2.txt",sep = "\t",quote = F,row.names = F)
# # 
# # MFP_input(TCGA_ACC_data,"TCGA_ACC_data")
# # MFP_input(TCGA_BLCA_data,"TCGA_BLCA_data")
# # MFP_input(TCGA_BRCA_data,"TCGA_BRCA_data")
# # MFP_input(TCGA_CESC_data,"TCGA_CESC_data")
# # MFP_input(TCGA_CHOL_data,"TCGA_CHOL_data")
# # MFP_input(TCGA_COAD_data,"TCGA_COAD_data")
# # MFP_input(TCGA_ESCA_data,"TCGA_ESCA_data")
# # MFP_input(TCGA_GBM_data,"TCGA_GBM_data")
# # MFP_input(TCGA_HNSC_data,"TCGA_HNSC_data")
# # MFP_input(TCGA_KICH_data,"TCGA_KICH_data")
# # MFP_input(TCGA_KIRC_data,"TCGA_KIRC_data")
# # MFP_input(TCGA_KIRP_data,"TCGA_KIRP_data")
# # MFP_input(TCGA_LGG_data,"TCGA_LGG_data")
# # MFP_input(TCGA_LIHC_data,"TCGA_LIHC_data")
# # MFP_input(TCGA_LUAD_data,"TCGA_LUAD_data")
# # MFP_input(TCGA_LUSC_data,"TCGA_LUSC_data")
# # MFP_input(TCGA_MESO_data,"TCGA_MESO_data")
# # MFP_input(TCGA_OV_data,"TCGA_OV_data")
# # MFP_input(TCGA_PAAD_data,"TCGA_PAAD_data")
# # MFP_input(TCGA_PCPG_data,"TCGA_PCPG_data")
# # MFP_input(TCGA_PRAD_data,"TCGA_PRAD_data")
# # MFP_input(TCGA_READ_data,"TCGA_READ_data")
# # MFP_input(TCGA_SARC_data,"TCGA_SARC_data")
# # MFP_input(TCGA_SKCM_data,"TCGA_SKCM_data")
# # MFP_input(TCGA_STAD_data,"TCGA_STAD_data")
# # MFP_input(TCGA_TGCT_data,"TCGA_TGCT_data")
# # MFP_input(TCGA_THCA_data,"TCGA_THCA_data")
# # MFP_input(TCGA_THYM_data,"TCGA_THYM_data")
# # MFP_input(TCGA_UCEC_data,"TCGA_UCEC_data")
# # MFP_input(TCGA_UCS_data,"TCGA_UCS_data")
# # MFP_input(TCGA_UVM_data,"TCGA_UVM_data")



##TCGA landscape
load("TCGA/TCGA_ACC_data.Rdata")
Landscape = dataset_landscape(TCGA_ACC_data,"TCGA_ACC_data")
TCGA_ACC_data[["Landscape"]] = Landscape
save(TCGA_ACC_data,file = 'TCGA/TCGA_ACC_data.Rdata')

load("TCGA/TCGA_BLCA_data.Rdata")
Landscape = dataset_landscape(TCGA_BLCA_data,"TCGA_BLCA_data")
TCGA_BLCA_data[["Landscape"]] = Landscape
save(TCGA_BLCA_data,file = 'TCGA/TCGA_BLCA_data.Rdata')

load("TCGA/TCGA_BRCA_data.Rdata")
Landscape = dataset_landscape(TCGA_BRCA_data,"TCGA_BRCA_data")
TCGA_BRCA_data[["Landscape"]] = Landscape
save(TCGA_BRCA_data,file = 'TCGA/TCGA_BRCA_data.Rdata')

load("TCGA/TCGA_CESC_data.Rdata")
Landscape = dataset_landscape(TCGA_CESC_data,"TCGA_CESC_data")
TCGA_CESC_data[["Landscape"]] = Landscape
save(TCGA_CESC_data,file = 'TCGA/TCGA_CESC_data.Rdata')

load("TCGA/TCGA_CHOL_data.Rdata")
Landscape = dataset_landscape(TCGA_CHOL_data,"TCGA_CHOL_data")
TCGA_CHOL_data[["Landscape"]] = Landscape
save(TCGA_CHOL_data,file = 'TCGA/TCGA_CHOL_data.Rdata')

load("TCGA/TCGA_COAD_data.Rdata")
Landscape = dataset_landscape(TCGA_COAD_data,"TCGA_COAD_data")
TCGA_COAD_data[["Landscape"]] = Landscape
save(TCGA_COAD_data,file = 'TCGA/TCGA_COAD_data.Rdata')


load("TCGA/TCGA_ESCA_data.Rdata")
Landscape = dataset_landscape(TCGA_ESCA_data,"TCGA_ESCA_data")
TCGA_ESCA_data[["Landscape"]] = Landscape
save(TCGA_ESCA_data,file = 'TCGA/TCGA_ESCA_data.Rdata')


load("TCGA/TCGA_GBM_data.Rdata")
Landscape = dataset_landscape(TCGA_GBM_data,"TCGA_GBM_data")
TCGA_GBM_data[["Landscape"]] = Landscape
save(TCGA_GBM_data,file = 'TCGA/TCGA_GBM_data.Rdata')


load("TCGA/TCGA_HNSC_data.Rdata")
Landscape = dataset_landscape(TCGA_HNSC_data,"TCGA_HNSC_data")
TCGA_HNSC_data[["Landscape"]] = Landscape
save(TCGA_HNSC_data,file = 'TCGA/TCGA_HNSC_data.Rdata')


load("TCGA/TCGA_KICH_data.Rdata")
Landscape = dataset_landscape(TCGA_KICH_data,"TCGA_KICH_data")
TCGA_KICH_data[["Landscape"]] = Landscape
save(TCGA_KICH_data,file = 'TCGA/TCGA_KICH_data.Rdata')


load("TCGA/TCGA_KIRC_data.Rdata")
Landscape = dataset_landscape(TCGA_KIRC_data,"TCGA_KIRC_data")
TCGA_KIRC_data[["Landscape"]] = Landscape
save(TCGA_KIRC_data,file = 'TCGA/TCGA_KIRC_data.Rdata')


load("TCGA/TCGA_KIRP_data.Rdata")
Landscape = dataset_landscape(TCGA_KIRP_data,"TCGA_KIRP_data")
TCGA_KIRP_data[["Landscape"]] = Landscape
save(TCGA_KIRP_data,file = 'TCGA/TCGA_KIRP_data.Rdata')


load("TCGA/TCGA_LGG_data.Rdata")
Landscape = dataset_landscape(TCGA_LGG_data,"TCGA_LGG_data")
TCGA_LGG_data[["Landscape"]] = Landscape
save(TCGA_LGG_data,file = 'TCGA/TCGA_LGG_data.Rdata')


load("TCGA/TCGA_LIHC_data.Rdata")
Landscape = dataset_landscape(TCGA_LIHC_data,"TCGA_LIHC_data")
TCGA_LIHC_data[["Landscape"]] = Landscape
save(TCGA_LIHC_data,file = 'TCGA/TCGA_LIHC_data.Rdata')


load("TCGA/TCGA_LUAD_data.Rdata")
Landscape = dataset_landscape(TCGA_LUAD_data,"TCGA_LUAD_data")
TCGA_LUAD_data[["Landscape"]] = Landscape
save(TCGA_LUAD_data,file = 'TCGA/TCGA_LUAD_data.Rdata')


load("TCGA/TCGA_LUSC_data.Rdata")
Landscape = dataset_landscape(TCGA_LUSC_data,"TCGA_LUSC_data")
TCGA_LUSC_data[["Landscape"]] = Landscape
save(TCGA_LUSC_data,file = 'TCGA/TCGA_LUSC_data.Rdata')


load("TCGA/TCGA_MESO_data.Rdata")
Landscape = dataset_landscape(TCGA_MESO_data,"TCGA_MESO_data")
TCGA_MESO_data[["Landscape"]] = Landscape
save(TCGA_MESO_data,file = 'TCGA/TCGA_MESO_data.Rdata')


load("TCGA/TCGA_OV_data.Rdata")
Landscape = dataset_landscape(TCGA_OV_data,"TCGA_OV_data")
TCGA_OV_data[["Landscape"]] = Landscape
save(TCGA_OV_data,file = 'TCGA/TCGA_OV_data.Rdata')


load("TCGA/TCGA_PAAD_data.Rdata")
Landscape = dataset_landscape(TCGA_PAAD_data,"TCGA_PAAD_data")
TCGA_PAAD_data[["Landscape"]] = Landscape
save(TCGA_PAAD_data,file = 'TCGA/TCGA_PAAD_data.Rdata')


load("TCGA/TCGA_PCPG_data.Rdata")
Landscape = dataset_landscape(TCGA_PCPG_data,"TCGA_PCPG_data")
TCGA_PCPG_data[["Landscape"]] = Landscape
save(TCGA_PCPG_data,file = 'TCGA/TCGA_PCPG_data.Rdata')


load("TCGA/TCGA_PRAD_data.Rdata")
Landscape = dataset_landscape(TCGA_PRAD_data,"TCGA_PRAD_data")
TCGA_PRAD_data[["Landscape"]] = Landscape
save(TCGA_PRAD_data,file = 'TCGA/TCGA_PRAD_data.Rdata')


load("TCGA/TCGA_READ_data.Rdata")
Landscape = dataset_landscape(TCGA_READ_data,"TCGA_READ_data")
TCGA_READ_data[["Landscape"]] = Landscape
save(TCGA_READ_data,file = 'TCGA/TCGA_READ_data.Rdata')


load("TCGA/TCGA_SARC_data.Rdata")
Landscape = dataset_landscape(TCGA_SARC_data,"TCGA_SARC_data")
TCGA_SARC_data[["Landscape"]] = Landscape
save(TCGA_SARC_data,file = 'TCGA/TCGA_SARC_data.Rdata')


load("TCGA/TCGA_SKCM_data.Rdata")
Landscape = dataset_landscape(TCGA_SKCM_data,"TCGA_SKCM_data")
TCGA_SKCM_data[["Landscape"]] = Landscape
save(TCGA_SKCM_data,file = 'TCGA/TCGA_SKCM_data.Rdata')


load("TCGA/TCGA_STAD_data.Rdata")
Landscape = dataset_landscape(TCGA_STAD_data,"TCGA_STAD_data")
TCGA_STAD_data[["Landscape"]] = Landscape
save(TCGA_STAD_data,file = 'TCGA/TCGA_STAD_data.Rdata')


load("TCGA/TCGA_TGCT_data.Rdata")
Landscape = dataset_landscape(TCGA_TGCT_data,"TCGA_TGCT_data")
TCGA_TGCT_data[["Landscape"]] = Landscape
save(TCGA_TGCT_data,file = 'TCGA/TCGA_TGCT_data.Rdata')


load("TCGA/TCGA_THCA_data.Rdata")
Landscape = dataset_landscape(TCGA_THCA_data,"TCGA_THCA_data")
TCGA_THCA_data[["Landscape"]] = Landscape
save(TCGA_THCA_data,file = 'TCGA/TCGA_THCA_data.Rdata')


load("TCGA/TCGA_THYM_data.Rdata")
Landscape = dataset_landscape(TCGA_THYM_data,"TCGA_THYM_data")
TCGA_THYM_data[["Landscape"]] = Landscape
save(TCGA_THYM_data,file = 'TCGA/TCGA_THYM_data.Rdata')


load("TCGA/TCGA_UCEC_data.Rdata")
Landscape = dataset_landscape(TCGA_UCEC_data,"TCGA_UCEC_data")
TCGA_UCEC_data[["Landscape"]] = Landscape
save(TCGA_UCEC_data,file = 'TCGA/TCGA_UCEC_data.Rdata')


load("TCGA/TCGA_UCS_data.Rdata")
Landscape = dataset_landscape(TCGA_UCS_data,"TCGA_UCS_data")
TCGA_UCS_data[["Landscape"]] = Landscape
save(TCGA_UCS_data,file = 'TCGA/TCGA_UCS_data.Rdata')


load("TCGA/TCGA_UVM_data.Rdata")
Landscape = dataset_landscape(TCGA_UVM_data,"TCGA_UVM_data")
TCGA_UVM_data[["Landscape"]] = Landscape
save(TCGA_UVM_data,file = 'TCGA/TCGA_UVM_data.Rdata')


##survival analysis 
OS_results = rbind(Biomarker_OS(TCGA_ACC_data,"TCGA_ACC_data"),
                   Biomarker_OS(TCGA_BLCA_data,"TCGA_BLCA_data"),
                   Biomarker_OS(TCGA_BRCA_data,"TCGA_BRCA_data"),
                   Biomarker_OS(TCGA_CESC_data,"TCGA_CESC_data"),
                   Biomarker_OS(TCGA_CHOL_data,"TCGA_CHOL_data"),
                   Biomarker_OS(TCGA_COAD_data,"TCGA_COAD_data"),
                   Biomarker_OS(TCGA_ESCA_data,"TCGA_ESCA_data"),
                   Biomarker_OS(TCGA_GBM_data,"TCGA_GBM_data"),
                   Biomarker_OS(TCGA_HNSC_data,"TCGA_HNSC_data"),
                   Biomarker_OS(TCGA_KICH_data,"TCGA_KICH_data"),
                   Biomarker_OS(TCGA_KIRC_data,"TCGA_KIRC_data"),
                   Biomarker_OS(TCGA_KIRP_data,"TCGA_KIRP_data"),
                   Biomarker_OS(TCGA_LGG_data,"TCGA_LGG_data"),
                   Biomarker_OS(TCGA_LIHC_data,"TCGA_LIHC_data"),
                   Biomarker_OS(TCGA_LUAD_data,"TCGA_LUAD_data"),
                   Biomarker_OS(TCGA_LUSC_data,"TCGA_LUSC_data"),
                   Biomarker_OS(TCGA_MESO_data,"TCGA_MESO_data"),
                   Biomarker_OS(TCGA_OV_data,"TCGA_OV_data"),
                   Biomarker_OS(TCGA_PAAD_data,"TCGA_PAAD_data"),
                   Biomarker_OS(TCGA_PCPG_data,"TCGA_PCPG_data"),
                   Biomarker_OS(TCGA_PRAD_data,"TCGA_PRAD_data"),
                   Biomarker_OS(TCGA_READ_data,"TCGA_READ_data"),
                   Biomarker_OS(TCGA_SARC_data,"TCGA_SARC_data"),
                   Biomarker_OS(TCGA_SKCM_data,"TCGA_SKCM_data"),
                   Biomarker_OS(TCGA_STAD_data,"TCGA_STAD_data"),
                   Biomarker_OS(TCGA_TGCT_data,"TCGA_TGCT_data"),
                   Biomarker_OS(TCGA_THCA_data,"TCGA_THCA_data"),
                   Biomarker_OS(TCGA_THYM_data,"TCGA_THYM_data"),
                   Biomarker_OS(TCGA_UCEC_data,"TCGA_UCEC_data"),
                   Biomarker_OS(TCGA_UCS_data,"TCGA_UCS_data"),
                   Biomarker_OS(TCGA_UVM_data,"TCGA_UVM_data"))

OS_results = as.data.frame(OS_results)

save(OS_results,file = "Results/TCGA_OS_results.Rdata")


