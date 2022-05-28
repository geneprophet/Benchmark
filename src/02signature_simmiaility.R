##Biomarker Similarity analysis
rm(list = ls())
CYT = c("GZMA","PRF1")
IFN_gamma = c("IDO1","CXCL10","CXCL9","HLA-DRA","STAT1","IFNG")
Expanded_immune_gene_signature = c("CD3D","IDO1","CIITA","CD3E","CCL5","GZMK","CD2","HLA-DRA","CXCL13","IL2RG","NKG7","HLA-E","CXCR6","LAG3","TAGAP","CXCL10","STAT1","GZMB")
T_cell_inflamed_GEP = c("TIGIT","PDCD1LG2","CD27","CD8A","LAG3","CD274","CXCR6","CMKLR1","NKG7","CCL5","PSMB10","IDO1","CXCL9","HLA-DQA1","CD276","STAT1","HLA-DRB1","HLA-E")
library(readr)
IPS_genes <- read_delim("marker/Immunophenogram/IPS_genes.txt",
                        delim = "\t", escape_double = FALSE,
                        trim_ws = TRUE)
for( i in unique(IPS_genes$CLASS)){
  b = IPS_genes$GENE[which(IPS_genes$CLASS==i)]
  assign(paste0("IPS_",i),b)
}
rm("b","i","IPS_genes")
IMPRES_Gene1 = c("PDCD1","CD27","CTLA4","CD40","CD86","CD28","CD80","CD274","CD86","CD40","CD86","CD40","CD28","CD40","TNFRSF14")
IMPRES_Gene2 = c("TNFSF4","PDCD1","TNFSF4","CD28","TNFSF4","CD86","TNFSF9","VSIR","HAVCR2","PDCD1","CD200","CD80","CD276","CD274","CD86")
CRMA_genes = c("MAGEA3", "CSAG3", "CSAG2","MAGEA2", "MAGEA2B", "CSAG1", "MAGEA12", "MAGEA6")
load("marker/resistance.program.RData")
The_immune_resistance_program_up = res.sig$resu.up
The_immune_resistance_program_down = res.sig$resu.down
rm("res.sig")
EMT_Stroma_core_signature = c("FLNA","EMP3","CALD1","FN1","FOXC2","LOX","FBN1","TNC")
F_TBRS = c("ACTA2", "ACTG2", "ADAM12", "ADAM19", "CNN1", "COL4A1", "CCN2", "CTPS1",
                 "RFLNB", "FSTL3", "HSPB1", "IGFBP3", "PXDC1", "SEMA7A", "SH3PXD2A", "TAGLN", 
                 "TGFBI", "TNS1", "TPM1")
library(readr)
TME_signature <- read_delim("marker/TMEScore.txt", 
                            delim = "\t", escape_double = FALSE, 
                            trim_ws = TRUE)
TMEscore_immune = TME_signature$Symbol[which(TME_signature$`TME-signature-group`=="TME-gene-A")]
TMEscore_stroma = TME_signature$Symbol[which(TME_signature$`TME-signature-group`=="TME-gene-B")]
rm("TME_signature")

Risk_score_genes = c("LEPR","PRLHR","NR2F2","PRL","NRP1","ANGPTL5","IGF1","TNFRSF10B","TNFRSF10A","PLAU","IFI30")
#Tertiary lymphoid structures (TLS) gene signature score
TLS_score_genes = c("CD79B","CD1D","CCR6","LAT","SKAP1","CETP","EIF1AY","RBP5","PTGDS")

MPS_genes <- read_delim("marker/MPS_gene_list.txt", 
                            delim = "\t", escape_double = FALSE, 
                            trim_ws = TRUE)
MPS_genes = MPS_genes$`Gene Symbol`
Renal_101_Immuno_signature = c("CD3G","CD3E","CD8B","THEMIS","TRAT1","GRAP2","CD247",
                               "CD2","CD96","PRF1","CD6","IL7R","ITK","GPR18","EOMES",
                               "SIT1","NLRC3","CD244","KLRD1","SH2D1A","CCL5","XCL2",
                               "CST7","GFI1","KCNA3","PSTPIP1")
IIS_TIS_signature <- read_delim("marker/IIS_TIS_signature.txt",
                                delim = "\t", escape_double = FALSE,
                                trim_ws = TRUE)
for (i in unique(IIS_TIS_signature$`Cell type`)) {
  b = IIS_TIS_signature$Symbol[which(IIS_TIS_signature$`Cell type`==i)]
  assign(paste0("IIS_TIS_signature_",i),b)
}
rm("b","i","IIS_TIS_signature")

APM_signature = c("HLA-A","HLA-B","HLA-C","B2M","TAP1","TAP2","TAPBP")

library("qusage")
IPRES_signatures = read.gmt("marker/IPRES_signatures.gmt")
for (i in names(IPRES_signatures)) {
  b = IPRES_signatures[[i]]
  assign(paste0("IPRES_signature_",i),b)
}
rm("b","i","IPRES_signatures")

C_ECM_genes <- read_delim("marker/C_ECM.txt", delim = ";", 
                          escape_double = FALSE, col_names = FALSE, 
                          trim_ws = TRUE)
C_ECM_genes = C_ECM_genes$X1

###
IMS_signature <- read_delim("marker/IMS_signature.txt",
                            delim = "\t", escape_double = FALSE,
                            trim_ws = TRUE)
for (i in unique(IMS_signature$`Immune cell type`)) {
  b = IMS_signature$Gene[which(IMS_signature$`Immune cell type` == i)]
  assign(paste0("IMS_signature_",i),b)
}
rm("b","i","IMS_signature")
##
load('marker/Pathway_Singatures.Rdata')
#Prepare signaure to ssGSEA
prepare_sig <- function(sig){
  leadingEdge.list <- sig$leadingEdge
  names(leadingEdge.list) <- sig$pathway
  return(leadingEdge.list)
}
PASS.PRE.Sigs <- prepare_sig(Pathway.Sigs$PASS_PRE)
for (i in names(PASS.PRE.Sigs)) {
  b = PASS.PRE.Sigs[[i]]
  assign(paste0("PASS_PRE_Sigs_",i),b)
}
PASS.ON.Sigs <- prepare_sig(Pathway.Sigs$PASS_ON)
for (i in names(PASS.ON.Sigs)) {
  b = PASS.ON.Sigs[[i]]
  assign(paste0("PASS_ON_Sigs_",i),b)
}
rm("b","i","PASS.PRE.Sigs","PASS.ON.Sigs","prepare_sig","Pathway.Sigs")
##
MIAS_signature = read_csv("marker/MIAS.txt")
MIAS_signature = MIAS_signature$x
## 
library("qusage")
MFP_signatures = read.gmt("marker/MFP.gmt")
for (i in names(MFP_signatures)) {
  b = MFP_signatures[[i]]
  assign(paste0("MFP_signatures_",i),b)
}
rm("b","i","MFP_signatures")


all_signature = ls()

similarity = matrix(nrow = length(all_signature),ncol = length(all_signature))
for (i in 1:length(all_signature)) {
  for (j in 1:length(all_signature)) {
    signature_a = get(all_signature[i])
    signature_b = get(all_signature[j])
    jaccard_index = length(intersect(signature_a,signature_b)) / length(union(signature_a,signature_b))
    similarity[i,j] = jaccard_index
  }
}
rownames(similarity) = all_signature
colnames(similarity) = all_signature
library(pheatmap)
pheatmap(similarity,cluster_cols = T,show_colnames = F)


#################################################################
##每个biomarker所包含的基因
biomarker = c("PD_L1","PD_1","PD_L2","CX3CL1","CTLA4","HLA_DRA","CXCL9","HRH1")
gene = c("CD274","PDCD1","PDCD1LG2","CX3CL1","CTLA4","HLA-DRA","CXCL9","HRH1")

biomarker = c(biomarker,rep("CYT_score",length(CYT)))
gene = c(gene,CYT)

biomarker = c(biomarker,rep("IFN_gamma",length(IFN_gamma)))
gene = c(gene,IFN_gamma)
biomarker = c(biomarker,rep("IFN_gamma_ssGSEA",length(IFN_gamma)))
gene = c(gene,IFN_gamma)

biomarker = c(biomarker,rep("Expanded_immune_gene_signature",length(Expanded_immune_gene_signature)))
gene = c(gene,Expanded_immune_gene_signature)
biomarker = c(biomarker,rep("Expanded_immune_gene_ssGSEA",length(Expanded_immune_gene_signature)))
gene = c(gene,Expanded_immune_gene_signature)

biomarker = c(biomarker,rep("T_cell_inflamed_GEP_score",length(T_cell_inflamed_GEP)))
gene = c(gene,T_cell_inflamed_GEP)
biomarker = c(biomarker,rep("T_cell_inflamed_GEP_ssGSEA",length(T_cell_inflamed_GEP)))
gene = c(gene,T_cell_inflamed_GEP)

IPS_genes <- read_delim("marker/Immunophenogram/IPS_genes.txt",
                        delim = "\t", escape_double = FALSE,
                        trim_ws = TRUE)
biomarker = c(biomarker,rep("Immunophenoscore",length(IPS_genes$GENE)))
gene = c(gene,IPS_genes$GENE)

biomarker = c(biomarker,rep("IMPRES_score",length(IMPRES_Gene1)))
gene = c(gene,IMPRES_Gene1)
biomarker = c(biomarker,rep("IMPRES_score",length(IMPRES_Gene2)))
gene = c(gene,IMPRES_Gene2)

biomarker = c(biomarker,rep("CRMA_score",length(CRMA_genes)))
gene = c(gene,CRMA_genes)
biomarker = c(biomarker,rep("CRMA_ssGSEA",length(CRMA_genes)))
gene = c(gene,CRMA_genes)

biomarker = c(biomarker,rep("EMT_Stroma_core_signature",length(EMT_Stroma_core_signature)))
gene = c(gene,EMT_Stroma_core_signature)
biomarker = c(biomarker,rep("EMT_Stroma_core_ssGSEA",length(EMT_Stroma_core_signature)))
gene = c(gene,EMT_Stroma_core_signature)

biomarker = c(biomarker,rep("F_TBRS",length(F_TBRS)))
gene = c(gene,F_TBRS)
biomarker = c(biomarker,rep("F_TBRS_ssGSEA",length(F_TBRS)))
gene = c(gene,F_TBRS)

TME_signature <- read_excel("marker/TMEScore.xlsx",
                            sheet = "S8", skip = 1)
biomarker = c(biomarker,rep("TMEscore",length(TME_signature$Symbol)))
gene = c(gene,TME_signature$Symbol)

load("marker/resistance.program.RData")
biomarker = c(biomarker,rep("The_immune_resistance_program",length(res.sig$resu.up)))
gene = c(gene,res.sig$resu.up)
biomarker = c(biomarker,rep("The_immune_resistance_program",length(res.sig$resu.down)))
gene = c(gene,res.sig$resu.down)

biomarker = c(biomarker,rep("RiskScore",length(Risk_score_genes)))
gene = c(gene,Risk_score_genes)
biomarker = c(biomarker,rep("RiskScore_ssGSEA",length(Risk_score_genes)))
gene = c(gene,Risk_score_genes)

biomarker = c(biomarker,rep("TLS_score",length(TLS_score_genes)))
gene = c(gene,TLS_score_genes)
biomarker = c(biomarker,rep("TLS_score_ssGSEA",length(TLS_score_genes)))
gene = c(gene,TLS_score_genes)

biomarker = c(biomarker,rep("MPS_score",length(MPS_genes)))
gene = c(gene,MPS_genes)

CTL_genes = c('CD8A', 'CD8B', 'GZMA', 'GZMB', 'PRF1')
biomarker = c(biomarker,rep("TIDE",length(CTL_genes)))
gene = c(gene,CTL_genes)

biomarker = c(biomarker,rep("Renal_101_Immuno_signature",length(Renal_101_Immuno_signature)))
gene = c(gene,Renal_101_Immuno_signature)
biomarker = c(biomarker,rep("Renal_101_Immuno_ssGSEA",length(Renal_101_Immuno_signature)))
gene = c(gene,Renal_101_Immuno_signature)

IIS_score_genes = c(IIS_TIS_signature_aDC,IIS_TIS_signature_Angiogenesis,`IIS_TIS_signature_Antigen presenting machinery`,`IIS_TIS_signature_B cells`,`IIS_TIS_signature_CD8 T cells`,`IIS_TIS_signature_Cytotoxic cells`,
                    IIS_TIS_signature_DC,IIS_TIS_signature_Eosinophils,IIS_TIS_signature_iDC,IIS_TIS_signature_Macrophages,`IIS_TIS_signature_Mast cells`,
                    IIS_TIS_signature_Neutrophils,`IIS_TIS_signature_NK CD56bright cells`,`IIS_TIS_signature_NK CD56dim cells`,`IIS_TIS_signature_NK cells`,IIS_TIS_signature_pDC,
                    `IIS_TIS_signature_T cells`,`IIS_TIS_signature_T helper cells`,`IIS_TIS_signature_Tcm cells`,`IIS_TIS_signature_Tem cells`,`IIS_TIS_signature_Tfh cells`,
                    `IIS_TIS_signature_Tgd cells`,`IIS_TIS_signature_Th1 cells`,`IIS_TIS_signature_Th17 cells`,`IIS_TIS_signature_Th2 cells`,`IIS_TIS_signature_Treg cells`)
biomarker = c(biomarker,rep("IIS_score",length(IIS_score_genes)))
gene = c(gene,IIS_score_genes)

TIS_score_genes = c(`IIS_TIS_signature_CD8 T cells`,`IIS_TIS_signature_T helper cells`,`IIS_TIS_signature_Tcm cells`,`IIS_TIS_signature_Tem cells`,`IIS_TIS_signature_Th1 cells`,`IIS_TIS_signature_Th17 cells`,`IIS_TIS_signature_Th2 cells`,`IIS_TIS_signature_Treg cells`)
biomarker = c(biomarker,rep("TIS_score",length(TIS_score_genes)))
gene = c(gene,TIS_score_genes)

biomarker = c(biomarker,rep("APM_score",length(APM_signature)))
gene = c(gene,APM_signature)

IPRES_signatures = read.gmt("marker/IPRES_signatures.gmt")
biomarker = c(biomarker,rep("IPRES_score",length(as.character(unlist(IPRES_signatures)))))
gene = c(gene,as.character(unlist(IPRES_signatures)))

biomarker = c(biomarker,rep("C_ECM_score",length(C_ECM_genes)))
gene = c(gene,C_ECM_genes)

MFP_signatures = read.gmt("marker/MFP.gmt")
biomarker = c(biomarker,rep("MFP",length(as.character(unlist(MFP_signatures)))))
gene = c(gene,as.character(unlist(MFP_signatures)))

load('marker/Pathway_Singatures.Rdata')
prepare_sig <- function(sig){
  leadingEdge.list <- sig$leadingEdge
  names(leadingEdge.list) <- sig$pathway
  return(leadingEdge.list)
}
PASS.PRE.Sigs <- prepare_sig(Pathway.Sigs$PASS_PRE)
PASS.ON.Sigs <- prepare_sig(Pathway.Sigs$PASS_ON)
biomarker = c(biomarker,rep("PASS_PRE",length(as.character(unlist(PASS.PRE.Sigs)))))
gene = c(gene,as.character(unlist(PASS.PRE.Sigs)))
biomarker = c(biomarker,rep("PASS_ON",length(as.character(unlist(PASS.ON.Sigs)))))
gene = c(gene,as.character(unlist(PASS.ON.Sigs)))

IMS_signature <- read_delim("marker/IMS_signature.txt",
                            delim = "\t", escape_double = FALSE,
                            trim_ws = TRUE)
biomarker = c(biomarker,rep("IMS_score",length(IMS_signature$Gene)))
gene = c(gene,IMS_signature$Gene)

biomarker = c(biomarker,rep("MIAS_score",length(MIAS_signature)))
gene = c(gene,MIAS_signature)

library(readr)
LM22 <- read_delim("marker/LM22.txt", delim = "\t", 
                   escape_double = FALSE, trim_ws = TRUE)
biomarker = c(biomarker,rep("CD8T_CIBERSORTx",length(LM22$`Gene symbol`)))
gene = c(gene,LM22$`Gene symbol`)

library(readr)
MCPcounter <- read_delim("marker/MCPcounter.txt", 
                         delim = "\t", escape_double = FALSE, 
                         trim_ws = TRUE)
biomarker = c(biomarker,rep("CD8T_MCPcounter",length(MCPcounter$`HUGO symbols`)))
gene = c(gene,MCPcounter$`HUGO symbols`)

library(qusage)
xCell_signature = read.gmt("marker/xCell_signature.gmt")
xCell_signature = as.character(unlist(xCell_signature))
biomarker = c(biomarker,rep("CD8T_xCell",length(xCell_signature)))
gene = c(gene,xCell_signature)

biomarker = c(biomarker,rep("Immunoscore_CIBERSORTx",length(LM22$`Gene symbol`)))
gene = c(gene,LM22$`Gene symbol`)


library(readr)
TR4 <- read_delim("marker/TR4.txt", delim = "\t", 
                  escape_double = FALSE, trim_ws = TRUE)
biomarker = c(biomarker,rep("Ecotype",length(TR4$`Gene Symbol`)))
gene = c(gene,TR4$`Gene Symbol`)
biomarker = c(biomarker,rep("Ecotype",length(LM22$`Gene symbol`)))
gene = c(gene,LM22$`Gene symbol`)

biomarker_genes = data.frame(biomarker = biomarker,gene = gene)
biomarker_genes = biomarker_genes[which(!duplicated.data.frame(biomarker_genes)),]
save(biomarker_genes,file = "marker/all_biomarker_genes.Rdata")



