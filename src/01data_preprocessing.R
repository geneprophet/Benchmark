#######################################
#data curation
##Mariathasan:IMvigor210
Mariathasan_data = list() 
library(readr)
clinical <- read_delim("ICB_data/IMvigor210/pData.txt", 
                         delim = "\t", escape_double = FALSE, 
                         trim_ws = TRUE)
clinical$Resp_NoResp = clinical$binaryResponse
clinical$Resp_NoResp[which(clinical$Resp_NoResp=="CR/PR")] = "Response"
clinical$Resp_NoResp[which(clinical$Resp_NoResp=="SD/PD")] = "No_Response"
Mariathasan_data[["Samples"]] = clinical
counts <- read_delim("ICB_data/Mariathasan et al/IMvigor210/counts.txt", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)
fData <- read_delim("ICB_data/Mariathasan et al/IMvigor210/fData.txt", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)
countToTpm(counts = counts[,2],effLen = fData$length)
tpm = apply(counts[,-1], 2, function(x){
  return(countToTpm(x,effLen = fData$length))
})
##only retain protein coding genes
tpm = tpm[which(is.element(fData$symbol,pcg)),]
rownames(tpm) = fData$symbol[which(is.element(fData$symbol,pcg))]
tpm = tpm[which(rowSums(tpm)>0),]
Mariathasan_data[["TPM"]] = tpm
save(Mariathasan_data,file = 'ICB_data/Mariathasan et al/Mariathasan_data.Rdata')


##Motzer
Motzer = list()
library(readxl)
clinical <- read_excel("ICB_data/Motzer et al/41591_2020_1044_MOESM3_ESM.xlsx", 
                       sheet = "S11_Clinical_data", skip = 1)
TPM <- read_excel("ICB_data/Motzer et al/41591_2020_1044_MOESM3_ESM.xlsx", 
                  sheet = "S13_Gene_expression_TPM", skip = 1)

tpm = as.matrix(TPM[,-1])
rownames(tpm) = TPM$HUGO
Motzer[["TPM"]] = tpm
library(dplyr)
clinical =  clinical %>% filter(TRT01P == "Avelumab+Axitinib")
clinical$Sample = clinical$ID
clinical$Resp_NoResp = clinical$PFS_P
clinical$Resp_NoResp[which(clinical$Resp_NoResp >= 6)] = "Response"
clinical$Resp_NoResp[which(clinical$Resp_NoResp < 6)] = "No_Response"
Motzer[["Samples"]] = clinical
tpm = Motzer$TPM
tpm = tpm[which(is.element(rownames(tpm),pcg)),]
dim(tpm)
Motzer$TPM = tpm
Motzer_data = Motzer
save(Motzer_data,file = 'ICB_data/Motzer et al/Motzer_data.Rdata')


####J-Y-Kim = Jung data
J_Y_Kim_data = list()
library(readr)
TPM <- read_delim("ICB_data/Jeong Yeon Kim et al/GSE135222_GEO_RNA-seq_omicslab_exp.tsv", 
                                                 delim = "\t", escape_double = FALSE, 
                                                 trim_ws = TRUE)
tpm = as.matrix(TPM[,-1])
genes = ENSEMBLIDtoSYMBOL(substr(TPM$gene_id,1,15))
rownames(tpm) = as.character(genes)
tpm = tpm[is.element(rownames(tpm),pcg),]
dim(tpm)
J_Y_Kim_data[["TPM"]] = tpm
clinical <- read_excel("ICB_data/Jeong Yeon Kim et al/13148_2020_907_MOESM6_ESM.xlsx", 
                       sheet = "Table1", skip = 2)
clinical$Sample = paste0("NSCLC",clinical$`Patient ID`)
clinical$Resp_NoResp = clinical$benefit
clinical$Resp_NoResp[which(clinical$Resp_NoResp=="Y")] = "Response"
clinical$Resp_NoResp[which(clinical$Resp_NoResp=="N")] = "No_Response"
J_Y_Kim_data[["Samples"]] = clinical
save(J_Y_Kim_data,file = 'ICB_data/Jeong Yeon Kim et al/J_Y_Kim_data.Rdata')



##MGH: collect from "Pathway signatures derived from on-treatment tumor specimens predict response to anti-PD1 blockade in metastatic melanoma"
MGH_ON_data$Batch14
rownames(MGH_ON_data$Batch17)[which(rownames(MGH_ON_data$Batch17) %in% rownames(MGH_ON_data$Batch14))] -> common
TPM = cbind(MGH_ON_data$Batch14[common,],MGH_ON_data$Batch17[common,],MGH_ON_data$SN0119610[common,],MGH_ON_data$SN0123099[common,],MGH_ON_data$SN0131794[common,])
MGH_ON_data[["TPM"]] = TPM
save(MGH_ON_data,file = "ICB_data/Du et al/MGH_ON_data.RData")

MGH_PRE_data$Batch14
rownames(MGH_PRE_data$Batch17)[which(rownames(MGH_PRE_data$Batch17) %in% rownames(MGH_PRE_data$Batch14))] -> common
TPM = cbind(MGH_PRE_data$Batch14[common,],MGH_PRE_data$Batch17[common,],MGH_PRE_data$SN0119610[common,],MGH_PRE_data$SN0123099[common,],MGH_PRE_data$SN0131794[common,])
MGH_PRE_data[["TPM"]] = TPM
save(MGH_PRE_data,file = "ICB_data/Du et al/MGH_PRE_data.RData")


## Gide_data,Lee_data, Riaz_data were collected from "Pathway signatures derived from on-treatment tumor specimens predict response to anti-PD1 blockade in metastatic melanoma"
##
load('./ICB_data/Gide et al/Gide_data.RData')
library(readxl)
clinical1 <- read_excel("ICB_data/Gide et al/1-s2.0-S1535610819300376-mmc2.xlsx", 
                                             skip = 1)
clinical1_1 = clinical1[which(clinical1$`RNA Sequencing`!="-" & clinical1$`RNA Sequencing`!="PRE and EDT"),]
clinical1_2 = clinical1[which(clinical1$`RNA Sequencing`=="PRE and EDT"),]
clinical1_3 = clinical1_2
clinical1_2$`RNA Sequencing` = "PRE"
clinical1_3$`RNA Sequencing` = "EDT"
clinical1 = rbind(clinical1_1,clinical1_2,clinical1_3)
clinical1["Sample"] = paste0("PD1_",clinical1$`Patient no.`,"_",clinical1$`RNA Sequencing`)

clinical2 <- read_excel("ICB_data/Gide et al/1-s2.0-S1535610819300376-mmc3.xlsx", 
                                             skip = 1)
clinical2_1 = clinical2[which(clinical2$`RNA Sequencing`!="-" & clinical2$`RNA Sequencing`!="PRE and EDT"),]
clinical2_2 = clinical2[which(clinical2$`RNA Sequencing`=="PRE and EDT"),]
clinical2_3 = clinical2_2
clinical2_2$`RNA Sequencing` = "PRE"
clinical2_3$`RNA Sequencing` = "EDT"
clinical2 = rbind(clinical2_1,clinical2_2,clinical2_3)
clinical2["Sample"] = paste0("ipiPD1_",clinical2$`Patient no.`,"_",clinical2$`RNA Sequencing`)

clinical = rbind.data.frame(clinical1,clinical2)
merged_clinical = merge.data.frame(Gide_data$Samples,clinical,by="Sample")
Gide_data[["Samples"]] = merged_clinical
save(Gide_data,file = "ICB_data/Gide et al/Gide_data.RData")



##Kim
library(readr)
ICB_Kim2018_Pembrolizumab_Gastric <- read_delim("ICB_data/Kim et al/ICB.Kim2018_Pembrolizumab_Gastric.clinical", 
                                                delim = "\t", escape_double = FALSE, 
                                                trim_ws = TRUE)
ICB_Kim2018_Pembrolizumab_Gastric$Response[which(ICB_Kim2018_Pembrolizumab_Gastric$Response==1)] = "Response"
ICB_Kim2018_Pembrolizumab_Gastric$Response[which(ICB_Kim2018_Pembrolizumab_Gastric$Response==0)] = "No_Response"
Kim_data = list()
Kim_data[["Sample"]] = data.frame(Sample=ICB_Kim2018_Pembrolizumab_Gastric$patient,Resp_NoResp=ICB_Kim2018_Pembrolizumab_Gastric$Response)
library(readr)
fpkm_normalized <- read_delim("ICB_data/Kim et al/ICB.Kim2018_Pembrolizumab_Gastric.self_subtract", 
                                                delim = "\t", escape_double = FALSE, 
                                                trim_ws = TRUE)
fpkm = as.matrix(fpkm_normalized[,-1])
rownames(fpkm) = ENTREZtoSYMBOL(as.character(fpkm_normalized$Entrez))
fpkm = fpkm[intersect(pcg,rownames(fpkm)),]
dim(fpkm)
Kim_data[["TPM"]] = fpkm
save(Kim_data,file = "ICB_data/Kim et al/Kim_data.RData")


##Jung data
library(readr)
GSE135222 <- read_delim("ICB_data/Jung et al/GSE135222_GEO_RNA-seq_omicslab_exp.tsv", 
                                                 delim = "\t", escape_double = FALSE, 
                                                 trim_ws = TRUE)
gene = c()
for( i in GSE135222$gene_id) {
  gene = c(gene,unlist(strsplit(i,split = ".",fixed = T))[1])
}
GSE135222$gene_id = gene
GSE135222$gene_id = ENSEMBLIDtoSYMBOL(GSE135222$gene_id)
GSE135222 = GSE135222[which(is.element(GSE135222$gene_id,pcg)),]
dim(GSE135222)
tpm = as.matrix(GSE135222[,-1])
rownames(tpm) = GSE135222$gene_id
Jung_data = list()
Jung_data[["TPM"]] = tpm
library(readr)
clinical <- read_delim("ICB_data/Jung et al/clinical.txt", 
                       delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE)

clinical$Sample = paste0("NSCLC",clinical$`Sample ID`)
clinical$Resp_NoResp = clinical$`Clinical benefit`
clinical$Resp_NoResp[which(clinical$Resp_NoResp=="DCB")] = "Response"
clinical$Resp_NoResp[which(clinical$Resp_NoResp=="NDB")] = "No_Response"
Jung_data[["Samples"]] = clinical
save(Jung_data,file = 'ICB_data/Jung et al/Jung_data.Rdata')


##Braun
library(readxl)
sample <- read_excel("ICB_data/Braun et al/41591_2020_839_MOESM2_ESM.xlsx", 
                     sheet = "S1_Clinical_and_Immune_Data", 
                     skip = 1)
tpm <- read_excel("ICB_data/Braun et al/41591_2020_839_MOESM2_ESM.xlsx", 
                  sheet = "S4A_RNA_Expression", skip = 1)
Braun_data = list()
sample["Sample"] = sample$RNA_ID
sample["Resp_NoResp"] = sample$Benefit
sample$Resp_NoResp[which(sample$Resp_NoResp=="CB")] = "Response"
sample$Resp_NoResp[which(sample$Resp_NoResp=="NCB")] = "No_Response"
sample$Resp_NoResp[which(sample$Resp_NoResp=="ICB")] = "NA"
sample = sample[which(is.element(sample$RNA_ID,colnames(tpm))),]
sample = sample[which(sample$Resp_NoResp != "NA"),]
sample = sample[which(sample$Arm == "NIVOLUMAB"),]
Braun_data[["Samples"]] = sample
dim(tpm)
TPM = as.matrix(tpm[,-1])
rownames(TPM) = tpm$gene_name
TPM = TPM[which(is.element(rownames(TPM),pcg)),]
dim(TPM)
Braun_data[["TPM"]] = TPM
save(Braun_data,file = 'ICB_data/Braun et al/Braun_data.Rdata')


##Hugo: collect from "Pathway signatures derived from on-treatment tumor specimens predict response to anti-PD1 blockade in metastatic melanoma"
Hugo_data[["Samples"]]['Sample'] = Hugo_data[["Samples"]]$Patient.ID
Hugo_data[["Samples"]]['Resp_NoResp'] = Hugo_data[["Samples"]]$Response
Hugo_data[["Samples"]]['Resp_NoResp'] -> a
a$Resp_NoResp[which(a$Resp_NoResp=="NR")] = "No_Response"
a$Resp_NoResp[which(a$Resp_NoResp=="R")] = "Response"
Hugo_data[["Samples"]]['Resp_NoResp'] = a
## integrate sample clinical information from 
library(readxl)
sample <- read_excel("ICB_data/Hugo et al/sample.xls", 
                     sheet = "S1A", skip = 2, n_max = 40)
sample = sample[which(sample$RNAseq==1 & !is.na(sample$`Patient ID`)),]
merged_samples = merge.data.frame(Hugo_data$Samples,sample,by.x = "Patient.ID",by.y = "Patient ID",all.x = T)
Hugo_data$Samples = merged_samples
save(Hugo_data,file = 'ICB_data/Hugo et al/Hugo_data.Rdata')


#RECIST (v.1.1) criteria: progressive disease (PD), stable disease (SD), mixed response (MR), partial response (PR), complete response (CD) 
#Response: MR,PR,CR
#NoResponse: PD,SD
Liu_data=list()
library(readr)
clinical <- read_delim("ICB_data/Liu et al/clinical.txt", 
                       delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE)
Liu_data[["Samples"]] = clinical
tpm <- read_delim("ICB_data/Liu et al/tpm.txt", 
                  delim = "\t", escape_double = FALSE, 
                  trim_ws = TRUE)
tpm2 = as.matrix(tpm[,-1])
rownames(tpm2) = tpm$PatientID
tpm2 = t(tpm2)
tpm2 = tpm2[intersect(rownames(tpm2),pcg),]
Liu_data[["TPM"]] = tpm2
Liu_data[["Samples"]]['Sample'] = Liu_data[["Samples"]]$PatientID
Liu_data[["Samples"]]['Resp_NoResp'] = Liu_data[["Samples"]]$BR
Liu_data[["Samples"]]['Resp_NoResp'] -> a
a$Resp_NoResp[which(a$Resp_NoResp=="PD")] = "No_Response"
a$Resp_NoResp[which(a$Resp_NoResp=="SD")] = "Response"
a$Resp_NoResp[which(a$Resp_NoResp=="MR")] = "Response"
a$Resp_NoResp[which(a$Resp_NoResp=="PR")] = "Response"
a$Resp_NoResp[which(a$Resp_NoResp=="CR")] = "Response"
Liu_data[["Samples"]]['Resp_NoResp'] = a
save(Liu_data,file = 'ICB_data/Liu et al/Liu_data.Rdata')


##Miao data
Miao_data=list()
library(readr)
clinical1 <- read_delim("ICB_data/Miao et al/clinical1.txt", 
                        delim = "\t", escape_double = FALSE, 
                        trim_ws = TRUE)
clinical2 <- read_delim("ICB_data/Miao et al/clinical2.txt", 
                        delim = "\t", escape_double = FALSE, 
                        trim_ws = TRUE)
common = intersect(colnames(clinical1),colnames(clinical2))
clinical = rbind.data.frame(clinical1[,common],clinical2[,common])
Miao_data[["Samples"]] = clinical
Miao_data[["Samples"]]['Sample'] = Miao_data[["Samples"]]$patient_id
Miao_data[["Samples"]]['Resp_NoResp'] = Miao_data[["Samples"]]$response_category
Miao_data[["Samples"]]['Resp_NoResp'] -> a
a$Resp_NoResp[which(a$Resp_NoResp=="no clinical benefit")] = "No_Response"
a$Resp_NoResp[which(a$Resp_NoResp=="intermediate benefit")] = "Response"
a$Resp_NoResp[which(a$Resp_NoResp=="clinical benefit")] = "Response"
Miao_data[["Samples"]]['Resp_NoResp'] = a
library(readr)
tpm <- read_delim("ICB_data/Miao et al/tpm.txt", 
                  delim = "\t", escape_double = FALSE, 
                  trim_ws = TRUE)
gene = c()
for( i in tpm$gene_id) {
  gene = c(gene,unlist(strsplit(i,split = ".",fixed = T))[1])
}
gene = ENSEMBLIDtoSYMBOL(gene)
tpm = as.matrix(tpm[,-1])
rownames(tpm) = gene
tpm = tpm[is.element(rownames(tpm),pcg),]
colnames(tpm)
cnames = c()
for (i in colnames(tpm)) {
  if(startsWith(i,"RCC")){
    j = paste(unlist(strsplit(i,split = "_"))[1:2],collapse = "_")
    cnames = c(cnames,j)
  }else{
    j = paste0("RCC-",substr(i,1,6))
    cnames = c(cnames,j)
  }
}
colnames(tpm) = cnames
Miao_data[["TPM"]] = tpm
save(Miao_data,file = 'ICB_data/Miao et al/Miao_data.Rdata')


#RECIST (v.1.1) criteria: progressive disease (PD), stable disease (SD), mixed response (MR), partial response (PR), complete response (CD) 
#Response: MR,PR,CR
#NoResponse: PD,SD
Snyder_data = list()
library(readr)
sample <- read_csv("ICB_data/Snyder et al/sample.txt")
sample = sample[,-c(38:43)]
Snyder_data[["Samples"]] = sample
Snyder_data[["Samples"]]['Sample'] = Snyder_data[["Samples"]]$ID
Snyder_data[["Samples"]]['Resp_NoResp'] = Snyder_data[["Samples"]]$`Best Response RECIST 1.1`
Snyder_data[["Samples"]]['Resp_NoResp'] -> a
a$Resp_NoResp[which(a$Resp_NoResp=="PD")] = "No_Response"
a$Resp_NoResp[which(a$Resp_NoResp=="SD")] = "No_Response"
a$Resp_NoResp[which(a$Resp_NoResp=="MR")] = "Response"
a$Resp_NoResp[which(a$Resp_NoResp=="PR")] = "Response"
a$Resp_NoResp[which(a$Resp_NoResp=="CR")] = "Response"
a$Resp_NoResp[which(a$Resp_NoResp=="Only scanned baseline")] = NA
a$Resp_NoResp[which(a$Resp_NoResp=="Only scanned at baseline")] = NA
Snyder_data[["Samples"]]['Resp_NoResp'] = a
Snyder_data$Samples$Resp_NoResp[which(Snyder_data$Samples$Resp_NoResp == "No_Response" & Snyder_data$Samples$`PFS (RECIST 1.1) in days` > 180)] = "Response"

tpm <- read_csv("ICB_data/Snyder et al/tpm.txt")
gene = ENSEMBLIDtoSYMBOL(tpm$...1)
tpm = as.matrix(tpm[,-1])
rownames(tpm) = gene
tpm = tpm[intersect(rownames(tpm),pcg),]
Snyder_data [["TPM"]] = tpm
save(Snyder_data,file = 'ICB_data/Snyder et al/Snyder_data.Rdata')




##Nathanson
library(readr)
cohort <- read_csv("ICB_data/Nathanson et al/cohort.txt")
library(readr)
cufflinks <- read_csv("ICB_data/Nathanson et al/cufflinks.csv")
unique(cufflinks$sample)
unique(cohort$`Study ID`)
cohort$Resp_NoResp = cohort$Benefit
cohort$Benefit
cohort$Resp_NoResp[which(cohort$Benefit==T)] = 'Response'
cohort$Resp_NoResp[which(cohort$Benefit==F)] = 'No_Response'

##########
#转换基因ID,只取protein-coding
gene = as.character(ENSEMBLIDtoSYMBOL(cufflinks$gene_id))
cufflinks$Symbol = gene
gene = intersect(gene,pcg)
expMatrix_all = vector(mode = "numeric",length = length(gene))
for (i in unique(cufflinks$sample)) {
  print(i)
  expMatrix = data.frame(Symbol = cufflinks$Symbol[which(cufflinks$sample==i)],i=cufflinks$FPKM[which(cufflinks$sample==i)])
  colnames(expMatrix) = c("Symbol",i)
  expMatrix = expMatrix[which(expMatrix$Symbol %in% gene),]
  expMatrix = expMatrix[which(!duplicated(expMatrix$Symbol)),]
  rownames(expMatrix) = expMatrix$Symbol
  expMatrix = expMatrix[gene,]
  print(dim(expMatrix))
  expMatrix_all = cbind(expMatrix_all,expMatrix[,2])
}

expMatrix_all = expMatrix_all[,-1]
rownames(expMatrix_all) = gene
colnames(expMatrix_all) = unique(cufflinks$sample)

Nathanson_data = list()
Nathanson_data[["Samples"]] = cohort
Nathanson_data[["TPM"]] = expMatrix_all
Nathanson_data$TPM = Nathanson_data$TPM[which(rowSums(Nathanson_data[["TPM"]])>0),]
save(Nathanson_data,file = "ICB_data/Nathanson et al/Nathanson_data.Rdata")

##VanAllen_data collect from "Pathway signatures derived from on-treatment tumor specimens predict response to anti-PD1 blockade in metastatic melanoma"

library(readxl)
tables2_revised <- read_excel("ICB_data/VanAllen et al/tables2_revised.xlsx", 
                              sheet = "transcriptome analysis (n=42)")
merged_samples = merge.data.frame(tables2_revised,VanAllen_data$Samples,by = "patient",all.x = T)
merged_samples$Sample = merged_samples$patient
merged_samples$Resp_NoResp[which(is.na(merged_samples$Resp_NoResp))] = "Response"
VanAllen_data$Samples = merged_samples

tpm = VanAllen_data$TPM
tpm = tpm[intersect(rownames(tpm),pcg),]
VanAllen_data$TPM = tpm
save(VanAllen_data,file = "ICB_data/VanAllen et al/VanAllen_data.Rdata")



###################################
###Manual correction and harmonization of relevant clinical information from different datasets
##Braun
Samples = Braun_data$Samples
Samples$PFS = round(30*Samples$PFS)
Samples$OS = round(30*Samples$OS)
Samples$Sex[which(Samples$Sex=="F")] = "Female"
Samples$Sex[which(Samples$Sex=="FEMALE")] = "Female"
Samples$Sex[which(Samples$Sex=="M")] = "Male"
Samples$Sex[which(Samples$Sex=="MALE")] = "Male"
Samples$Treatment = "anti-PD1"
Braun_data[["Clinical"]] = Samples
save(Braun_data,file = "ICB_data/Braun et al/Braun_data.Rdata")

##Gide
Samples = Gide_data$Samples
Samples$Age = Samples$`Age (Years)`
Samples$Sex[which(Samples$Sex=="F")] = "Female"
Samples$Sex[which(Samples$Sex=="M")] = "Male"
Samples$OS = Samples$`Overall Survival (Days)`
Samples$OS_CNSR = Samples$`Last Followup Status`
Samples$OS_CNSR[which(Samples$OS_CNSR=="Alive")] = 0
Samples$OS_CNSR[which(Samples$OS_CNSR=="Dead, melanoma")] = 1
Samples$OS_CNSR[which(Samples$OS_CNSR=="Dead")] = 1
Samples$PFS = Samples$`Progression Free Survival (Days)`
Samples$PFS_CNSR = Samples$OS_CNSR
Samples$Pre_On = Samples$PREEDT
Samples$Pre_On[which(Samples$Pre_On=="PRE")] = "Pre"
Samples$Pre_On[which(Samples$Pre_On=="EDT")] = "On"
Samples$Treatment = Samples$Treatment.x
Samples$Treatment[which(Samples$Treatment=="ipiPD1")] = "anti-CTLA4 + anti-PD1"
Samples$Treatment[which(Samples$Treatment=="PD1")] = "anti-PD1"
Gide_data[["Clinical"]] = Samples
save(Gide_data,file = "ICB_data/Gide et al/Gide_data.Rdata")


##Hugo
Samples = Hugo_data$Samples
Samples$Sex = Samples$Gender
Samples$Sex[which(Samples$Sex=="F")] = "Female"
Samples$Sex[which(Samples$Sex=="M")] = "Male"
Samples$OS = Samples$Overall.Survival
Samples$OS_CNSR = Samples$Vital.Status
Samples$OS_CNSR[which(Samples$OS_CNSR=="Dead")] = 1
Samples$OS_CNSR[which(Samples$OS_CNSR=="Alive")] = 0
Samples$Treatment[which(Samples$Treatment=="Pembrolizumab")] = "anti-PD1"
Samples$Pre_On = Samples$`Biopsy Time`
Samples$Pre_On[which(Samples$Pre_On=="pre-treatment")] = "Pre"
Samples$Pre_On[which(Samples$Pre_On=="on-treatment")] = "On"
Hugo_data[["Clinical"]] = Samples
save(Hugo_data,file = "ICB_data/Hugo et al/Hugo_data.Rdata")


##Jung
Samples = Jung_data$Samples
Samples$PFS_CNSR = Samples$PD_Event.1_Censoring.0
Samples$Treatment = "anti-PD1/PDL1"
Samples$Pre_On = "Pre"
Jung_data[["Clinical"]] = Samples
save(Jung_data,file = "ICB_data/Jung et al/Jung_data.Rdata")


##Kim
Samples = Kim_data$Samples
Samples$Treatment = "anti-PD1"
Samples$Pre_On = "Pre"
Kim_data[["Clinical"]] = Samples
save(Kim_data,file = "ICB_data/Kim et al/Kim_data.Rdata")


##Lee
Samples_1 = Lee_data$Pre_Samples
Samples_1$OS = round(30*as.numeric(Samples_1$Overall.survival..months.))
Samples_1$OS_CNSR = Samples_1$Vital.status
Samples_1$OS_CNSR[which(Samples_1$OS_CNSR=="Dead")] = 1
Samples_1$OS_CNSR[which(Samples_1$OS_CNSR=="Alive")] = 0
Samples_1$Pre_On =  "Pre"
Samples_1$Treatment = "anti-PD1"
Samples_1 = Samples_1[c("Sample","Resp_NoResp","OS","OS_CNSR","Pre_On","Treatment")]

Samples_2 = Lee_data$On_Samples
Samples_2$OS = round(30*as.numeric(Samples_2$Overall.survival..months.))
Samples_2$OS_CNSR = Samples_2$Vital.status
Samples_2$OS_CNSR[which(Samples_2$OS_CNSR=="Dead")] = 1
Samples_2$OS_CNSR[which(Samples_2$OS_CNSR=="Alive")] = 0
Samples_2$Pre_On =  "On"
Samples_2$Treatment = "anti-PD1"
Samples_2 = Samples_2[c("Sample","Resp_NoResp","OS","OS_CNSR","Pre_On","Treatment")]

Samples = rbind(Samples_1,Samples_2)
Lee_data[["Clinical"]] = Samples
save(Lee_data,file = "ICB_data/Lee et al/Lee_data.Rdata")


##Liu, prior anti-CTLA4
Samples = Liu_data$Samples
Samples$Sex = Samples$`gender (Male=1, Female=0)`
Samples$Sex[which(Samples$Sex==1)] = "Male"
Samples$Sex[which(Samples$Sex==0)] = "Female"
Samples$Treatment = "anti-PD1"
Samples$Treatment[which(Samples$priorCTLA4=="1")] = "anti-PD1 with previous anti-CTLA4 exposure"
Samples$Pre_On =  "Pre"
Samples$PFS_CNSR = Samples$dead
Samples$OS_CNSR = Samples$dead
Liu_data[["Clinical"]] = Samples
save(Liu_data,file = "ICB_data/Liu et al/Liu_data.Rdata")


##Mariathasan
Samples = Mariathasan_data$Samples
Samples$Sex[which(Samples$Sex=="F")] = "Female"
Samples$Sex[which(Samples$Sex=="M")] = "Male"
Samples$OS = Samples$os
Samples$OS_CNSR = Samples$censOS
Samples$Treatment = "anti-PDL1"
Samples$Pre_On =  "Pre"
Mariathasan_data[["Clinical"]] = Samples
save(Mariathasan_data,file = "ICB_data/Mariathasan et al/Mariathasan_data.Rdata")


##MGH_ON
Samples = MGH_ON_data$Samples
Samples$OS = Samples$os
Samples$OS_CNSR = Samples$censure_os
Samples$PFS = Samples$pfs
Samples$PFS_CNSR = Samples$censure_pfs
Samples$Pre_On =  "On"
Samples$Treatment = Samples$treatment
Samples$Treatment[which(Samples$Treatment=="IPIPD1")] = "anti-CTLA4 + anti-PD1"
Samples$Treatment[which(Samples$Treatment=="Nivolumab")] = "anti-PD1"
Samples$Treatment[which(Samples$Treatment=="PD1")] = "anti-PD1"
Samples$Treatment[which(Samples$Treatment=="PDL1")] = "anti-PDL1"
MGH_ON_data[["Clinical"]] = Samples
save(MGH_ON_data,file = "ICB_data/Du et al/MGH_ON_data.Rdata")


##MGH_PRE
Samples = MGH_PRE_data$Samples
Samples$OS = Samples$os
Samples$OS_CNSR = Samples$censure_os
Samples$PFS = Samples$pfs
Samples$PFS_CNSR = Samples$censure_pfs
Samples$Pre_On =  "On"
Samples$Treatment = Samples$treatment
Samples$Treatment[which(Samples$Treatment=="IPIPD1")] = "anti-CTLA4 + anti-PD1"
Samples$Treatment[which(Samples$Treatment=="Nivolumab")] = "anti-PD1"
Samples$Treatment[which(Samples$Treatment=="PD1")] = "anti-PD1"
Samples$Treatment[which(Samples$Treatment=="PDL1")] = "anti-PDL1"
MGH_PRE_data[["Clinical"]] = Samples
save(MGH_PRE_data,file = "ICB_data/Du et al/MGH_PRE_data.Rdata")


##Miao
Samples = Miao_data$Samples
Samples$Sex = Samples$sex
Samples$Sex[which(Samples$Sex=="F")] = "Female"
Samples$Sex[which(Samples$Sex=="M")] = "Male"
Samples$Sex[which(Samples$Sex=="FEMALE")] = "Female"
Samples$Sex[which(Samples$Sex=="MALE")] = "Male"
Samples$Age = Samples$age
Samples$OS = Samples$os_days
Samples$OS_CNSR = Samples$os_censor
Samples$PFS = Samples$pfs_days
Samples$PFS_CNSR = Samples$pfs_censor
Samples$Treatment = Samples$drug
Samples$Treatment[which(Samples$Treatment=="atezolizumab")] = "anti-PDL1"
Samples$Treatment[which(Samples$Treatment=="nivolumab")] = "anti-PD1"
Samples$Treatment[which(Samples$Treatment=="nivolumab + ipilimumab")] = "anti-CTLA4 + anti-PD1"
Samples$Pre_On = "Pre"
Miao_data[["Clinical"]] = Samples
save(Miao_data,file = "ICB_data/Miao et al/Miao_data.Rdata")


##Motzer
Samples = Motzer_data$Samples
Samples$Sex = Samples$SEX
Samples$Sex[which(Samples$Sex=="F")] = "Female"
Samples$Sex[which(Samples$Sex=="M")] = "Male"
Samples$Age = Samples$AGE
Samples$PFS = round(30*as.numeric(Samples$PFS_P))
Samples$PFS_CNSR = Samples$PFS_P_CNSR
Samples$Treatment = "anti-PDL1 + TKI"
Motzer_data[["Clinical"]] = Samples
save(Motzer_data,file = "ICB_data/Motzer et al/Motzer_data.Rdata")


##Nathanson
Samples = Nathanson_data$Samples
Samples$Sex = Samples$Gender
Samples$Sex[which(Samples$Sex=="F")] = "Female"
Samples$Sex[which(Samples$Sex=="M")] = "Male"
Samples$OS = round(30*Samples$OS)
Samples$OS_CNSR = Samples$Alive
Samples$OS_CNSR[which(Samples$Alive==0)] = 1
Samples$OS_CNSR[which(Samples$Alive==1)] = 0
Samples$Treatment = "anti-CTLA4"
Samples$Pre_On = Samples$`Biopsy pre or post ipi`
Samples$Pre_On[which(Samples$Pre_On=="pre")] = "Pre"
Samples$Pre_On[which(Samples$Pre_On=="post")] = "On"
Nathanson_data[["Clinical"]] = Samples
save(Nathanson_data,file = "ICB_data/Nathanson et al/Nathanson_data.Rdata")


##Riaz
Samples = Riaz_data$Samples
Samples$Pre_On = Samples$PreOn
Samples$Treatment = Samples$PopCateg
Samples$Treatment[Samples$Treatment=="NaiveOn"] = "anti-PD1"
Samples$Treatment[Samples$Treatment=="NaivePre"] = "anti-PD1"
Samples$Treatment[Samples$Treatment=="ProgOn"] = "anti-PD1 with previous anti-CTLA4 exposure"
Samples$Treatment[Samples$Treatment=="ProgPre"] = "anti-PD1 with previous anti-CTLA4 exposure"
merged_samples = merge.data.frame(Samples,Riaz_data$Surv.info,by = "PatientID",all.x = T)
merged_samples$OS_CNSR = merged_samples$OS_SOR
merged_samples$PFS_CNSR = merged_samples$PFS_SOR
Riaz_data[["Clinical"]] = merged_samples
save(Riaz_data,file = "ICB_data/Riaz et al/Riaz_data.Rdata")


##Snyder
Samples = Snyder_data$Samples
Samples$Sex[which(Samples$Sex=="F")] = "Female"
Samples$Sex[which(Samples$Sex=="M")] = "Male"
Samples$OS = Samples$`OS in days`
Samples$OS_CNSR = Samples$`Alive Status`
Samples$OS_CNSR[which(Samples$OS_CNSR=="N")] = 1
Samples$OS_CNSR[which(Samples$OS_CNSR=="Y")] = 0
Samples$PFS = Samples$`PFS (RECIST 1.1) in days`
Samples$PFS_CNSR = Samples$OS_CNSR
Samples$Treatment = "anti-PDL1"
Samples$Pre_On = "Pre"
Snyder_data[["Clinical"]] = Samples
save(Snyder_data,file = "ICB_data/Snyder et al/Snyder_data.Rdata")


##VanAllen_data
Samples = VanAllen_data$Samples
Samples$Age = Samples$age_start
Samples$Sex = Samples$gender
Samples$Sex[which(Samples$Sex=="female")] = "Female"
Samples$Sex[which(Samples$Sex=="male")] = "Male"
Samples$OS = Samples$overall_survival.y
Samples$OS_CNSR = Samples$dead.y
Samples$PFS = Samples$progression_free.y
Samples$PFS_CNSR = Samples$progression.y
Samples$Treatment = "anto-CTLA4"
Samples$Pre_On = "Pre"
VanAllen_data[["Clinical"]] = Samples
save(VanAllen_data,file = "ICB_data/VanAllen et al/VanAllen_data.Rdata")




##############################################################
## split datasets to sub-datasets based on the biopsy time and treatment

##biopsy time

load("ICB_data/Gide et al/Gide_data.RData")
load("ICB_data/Riaz et al/Riaz_data.RData")
load("ICB_data/Lee et al/Lee_data.RData")
load("ICB_data/Nathanson et al/Nathanson_data.Rdata")

table(Gide_data$Clinical$Pre_On)
table(Riaz_data$Clinical$Pre_On)
table(Lee_data$Clinical$Pre_On)
table(Nathanson_data$Clinical$Pre_On)


Gide_PRE_data = list()
PRE_samples = Gide_data$Clinical$Sample[which(Gide_data$Clinical$Pre_On == "Pre")]
Gide_PRE_data[["Samples"]] = Gide_data$Samples[which(is.element(Gide_data$Samples$Sample,PRE_samples)),]
Gide_PRE_data[["Clinical"]] = Gide_data$Clinical[which(is.element(Gide_data$Clinical$Sample,PRE_samples)),]
Gide_PRE_data[["TPM"]] = Gide_data$TPM[,which(is.element(colnames(Gide_data$TPM),PRE_samples))]
save(Gide_PRE_data,file = "ICB_data/Gide et al/Gide_PRE_data.Rdata")

Gide_ON_data = list()
ON_samples = Gide_data$Clinical$Sample[which(Gide_data$Clinical$Pre_On == "On")]
Gide_ON_data[["Samples"]] = Gide_data$Samples[which(is.element(Gide_data$Samples$Sample,ON_samples)),]
Gide_ON_data[["Clinical"]] = Gide_data$Clinical[which(is.element(Gide_data$Clinical$Sample,ON_samples)),]
Gide_ON_data[["TPM"]] = Gide_data$TPM[,which(is.element(colnames(Gide_data$TPM),ON_samples))]
save(Gide_ON_data,file = "ICB_data/Gide et al/Gide_ON_data.Rdata")

Riaz_PRE_data = list()
PRE_samples = Riaz_data$Clinical$Sample[which(Riaz_data$Clinical$Pre_On == "Pre")]
Riaz_PRE_data[["Samples"]] = Riaz_data$Samples[which(is.element(Riaz_data$Samples$Sample,PRE_samples)),]
Riaz_PRE_data[["Clinical"]] = Riaz_data$Clinical[which(is.element(Riaz_data$Clinical$Sample,PRE_samples)),]
Riaz_PRE_data[["TPM"]] = Riaz_data$TPM[,which(is.element(colnames(Riaz_data$TPM),PRE_samples))]
save(Riaz_PRE_data,file = "ICB_data/Riaz et al/Riaz_PRE_data.RData")

Riaz_ON_data = list()
ON_samples = Riaz_data$Clinical$Sample[which(Riaz_data$Clinical$Pre_On == "On")]
Riaz_ON_data[["Samples"]] = Riaz_data$Samples[which(is.element(Riaz_data$Samples$Sample,ON_samples)),]
Riaz_ON_data[["Clinical"]] = Riaz_data$Clinical[which(is.element(Riaz_data$Clinical$Sample,ON_samples)),]
Riaz_ON_data[["TPM"]] = Riaz_data$TPM[,which(is.element(colnames(Riaz_data$TPM),ON_samples))]
save(Riaz_ON_data,file = "ICB_data/Riaz et al/Riaz_ON_data.Rdata")

Lee_PRE_data = list()
PRE_samples = Lee_data$Clinical$Sample[which(Lee_data$Clinical$Pre_On == "Pre")]
Lee_PRE_data[["Samples"]] = Lee_data$Samples[which(is.element(Lee_data$Samples$Sample,PRE_samples)),]
Lee_PRE_data[["Clinical"]] = Lee_data$Clinical[which(is.element(Lee_data$Clinical$Sample,PRE_samples)),]
Lee_PRE_data[["TPM"]] = Lee_data$TPM[,which(is.element(colnames(Lee_data$TPM),PRE_samples))]
save(Lee_PRE_data,file = "ICB_data/Lee et al/Lee_PRE_data.Rdata")

Lee_ON_data = list()
ON_samples = Lee_data$Clinical$Sample[which(Lee_data$Clinical$Pre_On == "On")]
Lee_ON_data[["Samples"]] = Lee_data$Samples[which(is.element(Lee_data$Samples$Sample,ON_samples)),]
Lee_ON_data[["Clinical"]] = Lee_data$Clinical[which(is.element(Lee_data$Clinical$Sample,ON_samples)),]
Lee_ON_data[["TPM"]] = Lee_data$TPM[,which(is.element(colnames(Lee_data$TPM),ON_samples))]
save(Lee_ON_data,file = "ICB_data/Lee et al/Lee_ON_data.Rdata")


###treatment

load("ICB_data/Gide et al/Gide_data.RData")
load("ICB_data/Riaz et al/Riaz_data.RData")
load("ICB_data/Liu et al/Liu_data.Rdata")

table(Gide_data$Clinical$Treatment)
table(Riaz_data$Clinical$Treatment)
table(Liu_data$Clinical$Treatment)

Gide_MONO_data = list()
MONO_samples = Gide_data$Clinical$Sample[which(Gide_data$Clinical$Treatment == "anti-PD1")]
Gide_MONO_data[["Samples"]] = Gide_data$Samples[which(is.element(Gide_data$Samples$Sample,MONO_samples)),]
Gide_MONO_data[["Clinical"]] = Gide_data$Clinical[which(is.element(Gide_data$Clinical$Sample,MONO_samples)),]
Gide_MONO_data[["TPM"]] = Gide_data$TPM[,which(is.element(colnames(Gide_data$TPM),MONO_samples))]
save(Gide_MONO_data,file = "ICB_data/Gide et al/Gide_MONO_data.Rdata")

Gide_COMBINE_data = list()
COMBINE_samples = Gide_data$Clinical$Sample[which(Gide_data$Clinical$Treatment == "anti-CTLA4 + anti-PD1")]
Gide_COMBINE_data[["Samples"]] = Gide_data$Samples[which(is.element(Gide_data$Samples$Sample,COMBINE_samples)),]
Gide_COMBINE_data[["Clinical"]] = Gide_data$Clinical[which(is.element(Gide_data$Clinical$Sample,COMBINE_samples)),]
Gide_COMBINE_data[["TPM"]] = Gide_data$TPM[,which(is.element(colnames(Gide_data$TPM),COMBINE_samples))]
save(Gide_COMBINE_data,file = "ICB_data/Gide et al/Gide_COMBINE_data.Rdata")

Riaz_NAIVE_data = list()
NAIVE_samples = Riaz_data$Clinical$Sample[which(Riaz_data$Clinical$Treatment == "anti-PD1")]
Riaz_NAIVE_data[["Samples"]] = Riaz_data$Samples[which(is.element(Riaz_data$Samples$Sample,NAIVE_samples)),]
Riaz_NAIVE_data[["Clinical"]] = Riaz_data$Clinical[which(is.element(Riaz_data$Clinical$Sample,NAIVE_samples)),]
Riaz_NAIVE_data[["TPM"]] = Riaz_data$TPM[,which(is.element(colnames(Riaz_data$TPM),NAIVE_samples))]
save(Riaz_NAIVE_data,file = "ICB_data/Riaz et al/Riaz_NAIVE_data.Rdata")

Riaz_EXPOSURE_data = list()
EXPOSURE_samples = Riaz_data$Clinical$Sample[which(Riaz_data$Clinical$Treatment == "anti-PD1 with previous anti-CTLA4 exposure")]
Riaz_EXPOSURE_data[["Samples"]] = Riaz_data$Samples[which(is.element(Riaz_data$Samples$Sample,EXPOSURE_samples)),]
Riaz_EXPOSURE_data[["Clinical"]] = Riaz_data$Clinical[which(is.element(Riaz_data$Clinical$Sample,EXPOSURE_samples)),]
Riaz_EXPOSURE_data[["TPM"]] = Riaz_data$TPM[,which(is.element(colnames(Riaz_data$TPM),EXPOSURE_samples))]
save(Riaz_EXPOSURE_data,file = "ICB_data/Riaz et al/Riaz_EXPOSURE_data.Rdata")

Liu_NAIVE_data = list()
NAIVE_samples = Liu_data$Clinical$Sample[which(Liu_data$Clinical$Treatment == "anti-PD1")]
Liu_NAIVE_data[["Samples"]] = Liu_data$Samples[which(is.element(Liu_data$Samples$Sample,NAIVE_samples)),]
Liu_NAIVE_data[["Clinical"]] = Liu_data$Clinical[which(is.element(Liu_data$Clinical$Sample,NAIVE_samples)),]
Liu_NAIVE_data[["TPM"]] = Liu_data$TPM[,which(is.element(colnames(Liu_data$TPM),NAIVE_samples))]
save(Liu_NAIVE_data,file = "ICB_data/Liu et al/Liu_NAIVE_data.Rdata")

Liu_EXPOSURE_data = list()
EXPOSURE_samples = Liu_data$Clinical$Sample[which(Liu_data$Clinical$Treatment == "anti-PD1 with previous anti-CTLA4 exposure")]
Liu_EXPOSURE_data[["Samples"]] = Liu_data$Samples[which(is.element(Liu_data$Samples$Sample,EXPOSURE_samples)),]
Liu_EXPOSURE_data[["Clinical"]] = Liu_data$Clinical[which(is.element(Liu_data$Clinical$Sample,EXPOSURE_samples)),]
Liu_EXPOSURE_data[["TPM"]] = Liu_data$TPM[,which(is.element(colnames(Liu_data$TPM),EXPOSURE_samples))]
save(Liu_EXPOSURE_data,file = "ICB_data/Liu et al/Liu_EXPOSURE_data.Rdata")



