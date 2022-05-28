##analysis the associations between biomarkers and survival time
library(survminer)
library(survival)
library(png)
HR_95CI <- function(x){ 
  x <- summary(x)
  HR <-signif(x$coef[2], digits=2)
  HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
  HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
  res<- c(HR,HR.confint.lower,HR.confint.upper)
  names(res)<-c("Hazard Ratio","95% Upper CI","95% lower CI")
  return(res)
}
dir.create("Results/Results_survival_analysis")
dir.create("Results/Results_survival_analysis/OS")
dir.create("Results/Results_survival_analysis/OS/forest/")
dir.create("Results/Results_survival_analysis/PFS")
dir.create("Results/Results_survival_analysis/PFS/forest/")
Biomarker_OS <- function(data,name){
  dir.create(paste0("Results/Results_survival_analysis/OS/",name))
  dir.create(paste0("Results/Results_survival_analysis/OS/forest/",name))
  OS_results <- c()
  expMarker <- merge.data.frame(data$Clinical, data$Landscape,by="Sample",all.y = T)
  expMarker <- expMarker[which(!is.na(expMarker$OS)),]
  expMarker$OS <- as.numeric(expMarker$OS)
  expMarker$OS_CNSR <- as.numeric(expMarker$OS_CNSR)
  expMarker$Class <- ifelse(expMarker$PD_L1 > mean(expMarker$PD_L1),"PD_L1 High","PD_L1 Low")
  fit = survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/OS/",name,"/PD_L1.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/OS/forest/",name,"/PD_L1.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  OS_results = rbind(OS_results,PD_L1=result)
  
  expMarker$Class = ifelse(expMarker$PD_1 > mean(expMarker$PD_1),"PD_1 High","PD_1 Low")
  fit <- survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/OS/",name,"/PD_1.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/OS/forest/",name,"/PD_1.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  OS_results = rbind(OS_results,PD_1=result)
  
  expMarker$Class = ifelse(expMarker$PD_L2 > mean(expMarker$PD_L2),"PD_L2 High","PD_L2 Low")
  fit <- survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/OS/",name,"/PD_L2.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/OS/forest/",name,"/PD_L2.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  OS_results = rbind(OS_results,PD_L2=result)
  
  expMarker$Class = ifelse(expMarker$CX3CL1 > mean(expMarker$CX3CL1),"CX3CL1 High","CX3CL1 Low")
  fit <- survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/OS/",name,"/CX3CL1.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/OS/forest/",name,"/CX3CL1.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  OS_results = rbind(OS_results,CX3CL1=result)
  
  expMarker$Class = ifelse(expMarker$CTLA4 > mean(expMarker$CTLA4),"CTLA4 High","CTLA4 Low")
  fit <- survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/OS/",name,"/CTLA4.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/OS/forest/",name,"/CTLA4.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  OS_results = rbind(OS_results,CTLA4=result)
  
  expMarker$Class = ifelse(expMarker$CYT_score > mean(expMarker$CYT_score),"CYT_score High","CYT_score Low")
  fit <- survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/OS/",name,"/CYT_score.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/OS/forest/",name,"/CYT_score.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  OS_results = rbind(OS_results,CYT_score=result)
  
  expMarker$Class = ifelse(expMarker$HLA_DRA > mean(expMarker$HLA_DRA),"HLA_DRA High","HLA_DRA Low")
  fit <- survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/OS/",name,"/HLA_DRA.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/OS/forest/",name,"/HLA_DRA.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  OS_results = rbind(OS_results,HLA_DRA=result)
  
  expMarker$Class = ifelse(expMarker$IFN_gamma > mean(expMarker$IFN_gamma),"IFN_gamma High","IFN_gamma Low")
  fit <- survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/OS/",name,"/IFN_gamma.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/OS/forest/",name,"/IFN_gamma.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  OS_results = rbind(OS_results,IFN_gamma=result)
  
  expMarker$Class = ifelse(expMarker$Expanded_immune_gene_signature > mean(expMarker$Expanded_immune_gene_signature),"Expanded_immune_gene_signature High","Expanded_immune_gene_signature Low")
  fit <- survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/OS/",name,"/Expanded_immune_gene_signature.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/OS/forest/",name,"/Expanded_immune_gene_signature.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  OS_results = rbind(OS_results,Expanded_immune_gene_signature=result)
  
  expMarker$Class = ifelse(expMarker$T_cell_inflamed_GEP_score > mean(expMarker$T_cell_inflamed_GEP_score),"T_cell_inflamed_GEP_score High","T_cell_inflamed_GEP_score Low")
  fit <- survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/OS/",name,"/T_cell_inflamed_GEP_score.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/OS/forest/",name,"/T_cell_inflamed_GEP_score.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  OS_results = rbind(OS_results,T_cell_inflamed_GEP_score=result)
  
  expMarker$Class = ifelse(expMarker$Immunophenoscore > mean(expMarker$Immunophenoscore),"Immunophenoscore High","Immunophenoscore Low")
  fit <- survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/OS/",name,"/Immunophenoscore.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/OS/forest/",name,"/Immunophenoscore.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  OS_results = rbind(OS_results,Immunophenoscore=result)
  
  expMarker$Class = ifelse(expMarker$IMPRES_score > mean(expMarker$IMPRES_score),"IMPRES_score High","IMPRES_score Low")
  fit <- survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/OS/",name,"/IMPRES_score.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/OS/forest/",name,"/IMPRES_score.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  OS_results = rbind(OS_results,IMPRES_score=result)
  
  expMarker$Class = ifelse(expMarker$CRMA_score > mean(expMarker$CRMA_score),"CRMA_score High","CRMA_score Low")
  fit <- survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/OS/",name,"/CRMA_score.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/OS/forest/",name,"/CRMA_score.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  OS_results = rbind(OS_results,CRMA_score=result)
  
  expMarker$Class = ifelse(expMarker$The_immune_resistance_program > mean(expMarker$The_immune_resistance_program),"The_immune_resistance_program High","The_immune_resistance_program Low")
  fit <- survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/OS/",name,"/The_immune_resistance_program.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/OS/forest/",name,"/The_immune_resistance_program.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  OS_results = rbind(OS_results,The_immune_resistance_program=result)
  
  expMarker$Class = ifelse(expMarker$EMT_Stroma_core_signature > mean(expMarker$EMT_Stroma_core_signature),"EMT_Stroma_core_signature High","EMT_Stroma_core_signature Low")
  fit <- survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/OS/",name,"/EMT_Stroma_core_signature.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/OS/forest/",name,"/EMT_Stroma_core_signature.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  OS_results = rbind(OS_results,EMT_Stroma_core_signature=result)
  
  
  expMarker$Class = ifelse(expMarker$F_TBRS > mean(expMarker$F_TBRS),"F_TBRS High","F_TBRS Low")
  fit <- survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/OS/",name,"/F_TBRS.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/OS/forest/",name,"/F_TBRS.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  OS_results = rbind(OS_results,F_TBRS=result)
  
  expMarker$Class = ifelse(expMarker$TMEscore > mean(expMarker$TMEscore),"TMEscore High","TMEscore Low")
  fit <- survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/OS/",name,"/TMEscore.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/OS/forest/",name,"/TMEscore.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  OS_results = rbind(OS_results,TMEscore=result)
  
  expMarker$Class = ifelse(expMarker$RiskScore > mean(expMarker$RiskScore),"RiskScore High","RiskScore Low")
  fit <- survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/OS/",name,"/RiskScore.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/OS/forest/",name,"/RiskScore.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  OS_results = rbind(OS_results,RiskScore=result)
  
  expMarker$Class = ifelse(expMarker$TLS_score > mean(expMarker$TLS_score),"TLS_score High","TLS_score Low")
  fit <- survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/OS/",name,"/TLS_score.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class, data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/OS/forest/",name,"/TLS_score.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  OS_results = rbind(OS_results,TLS_score=result)
  
  expMarker$Class = ifelse(expMarker$CXCL9 > mean(expMarker$CXCL9),"CXCL9 High","CXCL9 Low")
  fit <- survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/OS/",name,"/CXCL9.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/OS/forest/",name,"/CXCL9.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  OS_results = rbind(OS_results,CXCL9=result)
  
  expMarker$Class = ifelse(expMarker$MPS_score > mean(expMarker$MPS_score),"MPS_score High","MPS_score Low")
  fit <- survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/OS/",name,"/MPS_score.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/OS/forest/",name,"/MPS_score.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  OS_results = rbind(OS_results,MPS_score=result)
  
  expMarker$Class = ifelse(expMarker$Renal_101_Immuno_signature > mean(expMarker$Renal_101_Immuno_signature),"Renal_101_Immuno_signature High","Renal_101_Immuno_signature Low")
  fit <- survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/OS/",name,"/Renal_101_Immuno_signature.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/OS/forest/",name,"/Renal_101_Immuno_signature.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  OS_results = rbind(OS_results,Renal_101_Immuno_signature=result)
  
  expMarker$Class = ifelse(expMarker$HRH1 > mean(expMarker$HRH1),"HRH1 High","HRH1 Low")
  fit <- survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/OS/",name,"/HRH1.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/OS/forest/",name,"/HRH1.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  OS_results = rbind(OS_results,HRH1=result)
  
  expMarker$Class = ifelse(expMarker$TIDE > mean(expMarker$TIDE),"TIDE High","TIDE Low")
  fit <- survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/OS/",name,"/TIDE.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/OS/forest/",name,"/TIDE.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  OS_results = rbind(OS_results,TIDE=result)
  
  expMarker$Class = ifelse(expMarker$IIS_score > mean(expMarker$IIS_score),"IIS_score High","IIS_score Low")
  fit <- survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/OS/",name,"/IIS_score.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/OS/forest/",name,"/IIS_score.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  OS_results = rbind(OS_results,IIS_score=result)
  
  expMarker$Class = ifelse(expMarker$TIS_score > mean(expMarker$TIS_score),"TIS_score High","TIS_score Low")
  fit <- survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/OS/",name,"/TIS_score.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/OS/forest/",name,"/TIS_score.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  OS_results = rbind(OS_results,TIS_score=result)
  
  expMarker$Class = ifelse(expMarker$APM_score > mean(expMarker$APM_score),"APM_score High","APM_score Low")
  fit <- survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/OS/",name,"/APM_score.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/OS/forest/",name,"/APM_score.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  OS_results = rbind(OS_results,APM_score=result)
  
  expMarker$Class = ifelse(expMarker$IPRES_score > mean(expMarker$IPRES_score),"IPRES_score High","IPRES_score Low")
  fit <- survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/OS/",name,"/IPRES_score.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/OS/forest/",name,"/IPRES_score.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  OS_results = rbind(OS_results,IPRES_score=result)
  
  expMarker$Class = ifelse(expMarker$C_ECM_score > mean(expMarker$C_ECM_score),"C_ECM_score High","C_ECM_score Low")
  fit <- survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/OS/",name,"/C_ECM_score.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/OS/forest/",name,"/C_ECM_score.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  OS_results = rbind(OS_results,C_ECM_score=result)
  
  expMarker$Class = ifelse(expMarker$IMS_score > mean(expMarker$IMS_score),"IMS_score High","IMS_score Low")
  fit <- survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/OS/",name,"/IMS_score.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/OS/forest/",name,"/IMS_score.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  OS_results = rbind(OS_results,IMS_score=result)
  
  expMarker$Class = ifelse(expMarker$PASS_PRE > mean(expMarker$PASS_PRE),"PASS_PRE High","PASS_PRE Low")
  fit <- survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/OS/",name,"/PASS_PRE.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/OS/forest/",name,"/PASS_PRE.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  OS_results = rbind(OS_results,PASS_PRE=result)
  
  expMarker$Class = ifelse(expMarker$PASS_ON > mean(expMarker$PASS_ON),"PASS_ON High","PASS_ON Low")
  fit <- survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/OS/",name,"/PASS_ON.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/OS/forest/",name,"/PASS_ON.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  OS_results = rbind(OS_results,PASS_ON=result)
  
  expMarker$Class = ifelse(expMarker$MIAS_score > mean(expMarker$MIAS_score),"MIAS_score High","MIAS_score Low")
  fit <- survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/OS/",name,"/MIAS_score.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/OS/forest/",name,"/MIAS_score.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  OS_results = rbind(OS_results,MIAS_score=result)
  
  expMarker$Class = ifelse(expMarker$CD8T_xCell > mean(expMarker$CD8T_xCell),"CD8T_xCell High","CD8T_xCell Low")
  fit <- survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/OS/",name,"/CD8T_xCell.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/OS/forest/",name,"/CD8T_xCell.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  OS_results = rbind(OS_results,CD8T_xCell=result)
  
  expMarker$Class = ifelse(expMarker$CD8T_MCPcounter > mean(expMarker$CD8T_MCPcounter),"CD8T_MCPcounter High","CD8T_MCPcounter Low")
  fit <- survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/OS/",name,"/CD8T_MCPcounter.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/OS/forest/",name,"/CD8T_MCPcounter.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  OS_results = rbind(OS_results,CD8T_MCPcounter=result)
  
  expMarker$Class = ifelse(expMarker$CD8T_CIBERSORTx > mean(expMarker$CD8T_CIBERSORTx),"CD8T_CIBERSORTx High","CD8T_CIBERSORTx Low")
  fit <- survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/OS/",name,"/CD8T_CIBERSORTx.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/OS/forest/",name,"/CD8T_CIBERSORTx.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  OS_results = rbind(OS_results,CD8T_CIBERSORTx=result)
  
  expMarker$Class = ifelse(expMarker$Immunoscore_CIBERSORTx > mean(expMarker$Immunoscore_CIBERSORTx),"Immunoscore_CIBERSORTx High","Immunoscore_CIBERSORTx Low")
  fit <- survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/OS/",name,"/Immunoscore_CIBERSORTx.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/OS/forest/",name,"/Immunoscore_CIBERSORTx.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  OS_results = rbind(OS_results,Immunoscore_CIBERSORTx=result)
  
  expMarker$Class = ifelse(expMarker$Ecotype == "CE9","CE9","NOT_CE9")
  # expMarker$Class[which(is.na(expMarker$Class))] = "NOT_CE9"
  fit <- survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/OS/",name,"/Ecotype.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  ## if there is no samples assign as CE9 EcoTyper, return NA
  if(length(levels(as.factor(expMarker$Class))) >1 ){
    res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
    hr = HR_95CI(res.cox)
    # if Inf exists, don't plot the forest
    if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
    ggsave(paste0("Results/Results_survival_analysis/OS/forest/",name,"/Ecotype.png"),width = 12,height = 6, bg="white") }
    result = c(hr,p$pval)
    names(result) = c(names(hr),"p_value")
    OS_results = rbind(OS_results,Ecotype=result)
  }else{
    result = c(NA,NA,NA,NA)
    names(result) = c("Hazard Ratio","95% Upper CI","95% lower CI","p_value")
    OS_results = rbind(OS_results,Ecotype=result)
  }
  
  
  
  expMarker$Class = ifelse(expMarker$MFP == "IE","IE","NOT_IE")
  # expMarker$Class[which(is.na(expMarker$Class))] = "NOT_IE"
  fit <- survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/OS/",name,"/MFP.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/OS/forest/",name,"/MFP.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  OS_results = rbind(OS_results,MFP=result)
  
  ####averager2ssGSEA: 9
  expMarker$Class = ifelse(expMarker$IFN_gamma_ssGSEA > mean(expMarker$IFN_gamma_ssGSEA),"IFN_gamma_ssGSEA High","IFN_gamma_ssGSEA Low")
  fit <- survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/OS/",name,"/IFN_gamma_ssGSEA.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/OS/forest/",name,"/IFN_gamma_ssGSEA.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  OS_results = rbind(OS_results,IFN_gamma_ssGSEA=result)
  
  expMarker$Class = ifelse(expMarker$Expanded_immune_gene_ssGSEA > mean(expMarker$Expanded_immune_gene_ssGSEA),"Expanded_immune_gene_ssGSEA High","Expanded_immune_gene_ssGSEA Low")
  fit <- survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/OS/",name,"/Expanded_immune_gene_ssGSEA.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/OS/forest/",name,"/Expanded_immune_gene_ssGSEA.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  OS_results = rbind(OS_results,Expanded_immune_gene_ssGSEA=result)
  
  expMarker$Class = ifelse(expMarker$T_cell_inflamed_GEP_ssGSEA > mean(expMarker$T_cell_inflamed_GEP_ssGSEA),"T_cell_inflamed_GEP_ssGSEA High","T_cell_inflamed_GEP_ssGSEA Low")
  fit <- survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/OS/",name,"/T_cell_inflamed_GEP_ssGSEA.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/OS/forest/",name,"/T_cell_inflamed_GEP_ssGSEA.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  OS_results = rbind(OS_results,T_cell_inflamed_GEP_ssGSEA=result)
  
  expMarker$Class = ifelse(expMarker$CRMA_ssGSEA > mean(expMarker$CRMA_ssGSEA),"CRMA_ssGSEA High","CRMA_ssGSEA Low")
  fit <- survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/OS/",name,"/CRMA_ssGSEA.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/OS/forest/",name,"/CRMA_ssGSEA.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  OS_results = rbind(OS_results,CRMA_ssGSEA=result)
  
  expMarker$Class = ifelse(expMarker$EMT_Stroma_core_ssGSEA > mean(expMarker$EMT_Stroma_core_ssGSEA),"EMT_Stroma_core_ssGSEA High","EMT_Stroma_core_ssGSEA Low")
  fit <- survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/OS/",name,"/EMT_Stroma_core_ssGSEA.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/OS/forest/",name,"/EMT_Stroma_core_ssGSEA.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  OS_results = rbind(OS_results,EMT_Stroma_core_ssGSEA=result)
  
  expMarker$Class = ifelse(expMarker$F_TBRS_ssGSEA > mean(expMarker$F_TBRS_ssGSEA),"F_TBRS_ssGSEA High","F_TBRS_ssGSEA Low")
  fit <- survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/OS/",name,"/F_TBRS_ssGSEA.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/OS/forest/",name,"/F_TBRS_ssGSEA.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  OS_results = rbind(OS_results,F_TBRS_ssGSEA=result)
  
  expMarker$Class = ifelse(expMarker$RiskScore_ssGSEA > mean(expMarker$RiskScore_ssGSEA),"RiskScore_ssGSEA High","RiskScore_ssGSEA Low")
  fit <- survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/OS/",name,"/RiskScore_ssGSEA.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/OS/forest/",name,"/RiskScore_ssGSEA.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  OS_results = rbind(OS_results,RiskScore_ssGSEA=result)
  
  expMarker$Class = ifelse(expMarker$TLS_score_ssGSEA > mean(expMarker$TLS_score_ssGSEA),"TLS_score_ssGSEA High","TLS_score_ssGSEA Low")
  fit <- survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/OS/",name,"/TLS_score_ssGSEA.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/OS/forest/",name,"/TLS_score_ssGSEA.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  OS_results = rbind(OS_results,TLS_score_ssGSEA=result)
  
  expMarker$Class = ifelse(expMarker$Renal_101_Immuno_ssGSEA > mean(expMarker$Renal_101_Immuno_ssGSEA),"Renal_101_Immuno_ssGSEA High","Renal_101_Immuno_ssGSEA Low")
  fit <- survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/OS/",name,"/Renal_101_Immuno_ssGSEA.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/OS/forest/",name,"/Renal_101_Immuno_ssGSEA.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  OS_results = rbind(OS_results,Renal_101_Immuno_ssGSEA=result)
  
  
  Dataset = rep(name,39+9) 
  Biomarker = rownames(OS_results)
  OS_results = cbind(OS_results,Dataset,Biomarker)
  
  return(OS_results)
}

Biomarker_PFS <- function(data,name){
  dir.create(paste0("Results/Results_survival_analysis/PFS/",name))
  dir.create(paste0("Results/Results_survival_analysis/PFS/forest/",name))
  PFS_results <- c()
  expMarker <- merge.data.frame(data$Clinical, data$Landscape,by="Sample",all.y = T)
  expMarker <- expMarker[which(!is.na(expMarker$PFS)),]
  expMarker$PFS <- as.numeric(expMarker$PFS)
  expMarker$PFS_CNSR <- as.numeric(expMarker$PFS_CNSR)
  expMarker$Class <- ifelse(expMarker$PD_L1 > mean(expMarker$PD_L1),"PD_L1 High","PD_L1 Low")
  fit = survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/PFS/",name,"/PD_L1.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/PFS/forest/",name,"/PD_L1.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  PFS_results = rbind(PFS_results,PD_L1=result)
  
  expMarker$Class = ifelse(expMarker$PD_1 > mean(expMarker$PD_1),"PD_1 High","PD_1 Low")
  fit <- survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/PFS/",name,"/PD_1.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/PFS/forest/",name,"/PD_1.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  PFS_results = rbind(PFS_results,PD_1=result)
  
  expMarker$Class = ifelse(expMarker$PD_L2 > mean(expMarker$PD_L2),"PD_L2 High","PD_L2 Low")
  fit <- survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/PFS/",name,"/PD_L2.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/PFS/forest/",name,"/PD_L2.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  PFS_results = rbind(PFS_results,PD_L2=result)
  
  expMarker$Class = ifelse(expMarker$CX3CL1 > mean(expMarker$CX3CL1),"CX3CL1 High","CX3CL1 Low")
  fit <- survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/PFS/",name,"/CX3CL1.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/PFS/forest/",name,"/CX3CL1.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  PFS_results = rbind(PFS_results,CX3CL1=result)
  
  expMarker$Class = ifelse(expMarker$CTLA4 > mean(expMarker$CTLA4),"CTLA4 High","CTLA4 Low")
  fit <- survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/PFS/",name,"/CTLA4.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/PFS/forest/",name,"/CTLA4.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  PFS_results = rbind(PFS_results,CTLA4=result)
  
  expMarker$Class = ifelse(expMarker$CYT_score > mean(expMarker$CYT_score),"CYT_score High","CYT_score Low")
  fit <- survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/PFS/",name,"/CYT_score.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/PFS/forest/",name,"/CYT_score.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  PFS_results = rbind(PFS_results,CYT_score=result)
  
  expMarker$Class = ifelse(expMarker$HLA_DRA > mean(expMarker$HLA_DRA),"HLA_DRA High","HLA_DRA Low")
  fit <- survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/PFS/",name,"/HLA_DRA.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/PFS/forest/",name,"/HLA_DRA.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  PFS_results = rbind(PFS_results,HLA_DRA=result)
  
  expMarker$Class = ifelse(expMarker$IFN_gamma > mean(expMarker$IFN_gamma),"IFN_gamma High","IFN_gamma Low")
  fit <- survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/PFS/",name,"/IFN_gamma.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/PFS/forest/",name,"/IFN_gamma.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  PFS_results = rbind(PFS_results,IFN_gamma=result)
  
  expMarker$Class = ifelse(expMarker$Expanded_immune_gene_signature > mean(expMarker$Expanded_immune_gene_signature),"Expanded_immune_gene_signature High","Expanded_immune_gene_signature Low")
  fit <- survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/PFS/",name,"/Expanded_immune_gene_signature.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/PFS/forest/",name,"/Expanded_immune_gene_signature.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  PFS_results = rbind(PFS_results,Expanded_immune_gene_signature=result)
  
  expMarker$Class = ifelse(expMarker$T_cell_inflamed_GEP_score > mean(expMarker$T_cell_inflamed_GEP_score),"T_cell_inflamed_GEP_score High","T_cell_inflamed_GEP_score Low")
  fit <- survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/PFS/",name,"/T_cell_inflamed_GEP_score.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/PFS/forest/",name,"/T_cell_inflamed_GEP_score.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  PFS_results = rbind(PFS_results,T_cell_inflamed_GEP_score=result)
  
  expMarker$Class = ifelse(expMarker$Immunophenoscore > mean(expMarker$Immunophenoscore),"Immunophenoscore High","Immunophenoscore Low")
  fit <- survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/PFS/",name,"/Immunophenoscore.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  PFS_results = rbind(PFS_results,Immunophenoscore=result)
  # if Inf exists, don't plot the forest
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/PFS/forest/",name,"/Immunophenoscore.png"),width = 12,height = 6, bg="white") }

  
  expMarker$Class = ifelse(expMarker$IMPRES_score > mean(expMarker$IMPRES_score),"IMPRES_score High","IMPRES_score Low")
  fit <- survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/PFS/",name,"/IMPRES_score.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/PFS/forest/",name,"/IMPRES_score.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  PFS_results = rbind(PFS_results,IMPRES_score=result)
  
  expMarker$Class = ifelse(expMarker$CRMA_score > mean(expMarker$CRMA_score),"CRMA_score High","CRMA_score Low")
  fit <- survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/PFS/",name,"/CRMA_score.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/PFS/forest/",name,"/CRMA_score.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  PFS_results = rbind(PFS_results,CRMA_score=result)
  
  expMarker$Class = ifelse(expMarker$The_immune_resistance_program > mean(expMarker$The_immune_resistance_program),"The_immune_resistance_program High","The_immune_resistance_program Low")
  fit <- survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/PFS/",name,"/The_immune_resistance_program.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/PFS/forest/",name,"/The_immune_resistance_program.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  PFS_results = rbind(PFS_results,The_immune_resistance_program=result)
  
  expMarker$Class = ifelse(expMarker$EMT_Stroma_core_signature > mean(expMarker$EMT_Stroma_core_signature),"EMT_Stroma_core_signature High","EMT_Stroma_core_signature Low")
  fit <- survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/PFS/",name,"/EMT_Stroma_core_signature.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/PFS/forest/",name,"/EMT_Stroma_core_signature.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  PFS_results = rbind(PFS_results,EMT_Stroma_core_signature=result)
  
  
  expMarker$Class = ifelse(expMarker$F_TBRS > mean(expMarker$F_TBRS),"F_TBRS High","F_TBRS Low")
  fit <- survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/PFS/",name,"/F_TBRS.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/PFS/forest/",name,"/F_TBRS.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  PFS_results = rbind(PFS_results,F_TBRS=result)
  
  expMarker$Class = ifelse(expMarker$TMEscore > mean(expMarker$TMEscore),"TMEscore High","TMEscore Low")
  fit <- survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/PFS/",name,"/TMEscore.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/PFS/forest/",name,"/TMEscore.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  PFS_results = rbind(PFS_results,TMEscore=result)
  
  expMarker$Class = ifelse(expMarker$RiskScore > mean(expMarker$RiskScore),"RiskScore High","RiskScore Low")
  fit <- survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/PFS/",name,"/RiskScore.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/PFS/forest/",name,"/RiskScore.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  PFS_results = rbind(PFS_results,RiskScore=result)
  
  expMarker$Class = ifelse(expMarker$TLS_score > mean(expMarker$TLS_score),"TLS_score High","TLS_score Low")
  fit <- survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/PFS/",name,"/TLS_score.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class, data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/PFS/forest/",name,"/TLS_score.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  PFS_results = rbind(PFS_results,TLS_score=result)
  
  expMarker$Class = ifelse(expMarker$CXCL9 > mean(expMarker$CXCL9),"CXCL9 High","CXCL9 Low")
  fit <- survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/PFS/",name,"/CXCL9.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/PFS/forest/",name,"/CXCL9.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  PFS_results = rbind(PFS_results,CXCL9=result)
  
  expMarker$Class = ifelse(expMarker$MPS_score > mean(expMarker$MPS_score),"MPS_score High","MPS_score Low")
  fit <- survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/PFS/",name,"/MPS_score.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/PFS/forest/",name,"/MPS_score.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  PFS_results = rbind(PFS_results,MPS_score=result)
  
  expMarker$Class = ifelse(expMarker$Renal_101_Immuno_signature > mean(expMarker$Renal_101_Immuno_signature),"Renal_101_Immuno_signature High","Renal_101_Immuno_signature Low")
  fit <- survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/PFS/",name,"/Renal_101_Immuno_signature.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/PFS/forest/",name,"/Renal_101_Immuno_signature.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  PFS_results = rbind(PFS_results,Renal_101_Immuno_signature=result)
  
  expMarker$Class = ifelse(expMarker$HRH1 > mean(expMarker$HRH1),"HRH1 High","HRH1 Low")
  fit <- survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/PFS/",name,"/HRH1.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/PFS/forest/",name,"/HRH1.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  PFS_results = rbind(PFS_results,HRH1=result)
  
  expMarker$Class = ifelse(expMarker$TIDE > mean(expMarker$TIDE),"TIDE High","TIDE Low")
  fit <- survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/PFS/",name,"/TIDE.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/PFS/forest/",name,"/TIDE.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  PFS_results = rbind(PFS_results,TIDE=result)
  
  expMarker$Class = ifelse(expMarker$IIS_score > mean(expMarker$IIS_score),"IIS_score High","IIS_score Low")
  fit <- survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/PFS/",name,"/IIS_score.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/PFS/forest/",name,"/IIS_score.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  PFS_results = rbind(PFS_results,IIS_score=result)
  
  expMarker$Class = ifelse(expMarker$TIS_score > mean(expMarker$TIS_score),"TIS_score High","TIS_score Low")
  fit <- survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/PFS/",name,"/TIS_score.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/PFS/forest/",name,"/TIS_score.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  PFS_results = rbind(PFS_results,TIS_score=result)
  
  expMarker$Class = ifelse(expMarker$APM_score > mean(expMarker$APM_score),"APM_score High","APM_score Low")
  fit <- survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/PFS/",name,"/APM_score.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/PFS/forest/",name,"/APM_score.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  PFS_results = rbind(PFS_results,APM_score=result)
  
  expMarker$Class = ifelse(expMarker$IPRES_score > mean(expMarker$IPRES_score),"IPRES_score High","IPRES_score Low")
  fit <- survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/PFS/",name,"/IPRES_score.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/PFS/forest/",name,"/IPRES_score.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  PFS_results = rbind(PFS_results,IPRES_score=result)
  
  expMarker$Class = ifelse(expMarker$C_ECM_score > mean(expMarker$C_ECM_score),"C_ECM_score High","C_ECM_score Low")
  fit <- survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/PFS/",name,"/C_ECM_score.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/PFS/forest/",name,"/C_ECM_score.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  PFS_results = rbind(PFS_results,C_ECM_score=result)
  
  expMarker$Class = ifelse(expMarker$IMS_score > mean(expMarker$IMS_score),"IMS_score High","IMS_score Low")
  fit <- survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/PFS/",name,"/IMS_score.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/PFS/forest/",name,"/IMS_score.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  PFS_results = rbind(PFS_results,IMS_score=result)
  
  expMarker$Class = ifelse(expMarker$PASS_PRE > mean(expMarker$PASS_PRE),"PASS_PRE High","PASS_PRE Low")
  fit <- survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/PFS/",name,"/PASS_PRE.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/PFS/forest/",name,"/PASS_PRE.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  PFS_results = rbind(PFS_results,PASS_PRE=result)
  
  expMarker$Class = ifelse(expMarker$PASS_ON > mean(expMarker$PASS_ON),"PASS_ON High","PASS_ON Low")
  fit <- survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/PFS/",name,"/PASS_ON.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/PFS/forest/",name,"/PASS_ON.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  PFS_results = rbind(PFS_results,PASS_ON=result)
  
  expMarker$Class = ifelse(expMarker$MIAS_score > mean(expMarker$MIAS_score),"MIAS_score High","MIAS_score Low")
  fit <- survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/PFS/",name,"/MIAS_score.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/PFS/forest/",name,"/MIAS_score.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  PFS_results = rbind(PFS_results,MIAS_score=result)
  
  expMarker$Class = ifelse(expMarker$CD8T_xCell > mean(expMarker$CD8T_xCell),"CD8T_xCell High","CD8T_xCell Low")
  fit <- survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/PFS/",name,"/CD8T_xCell.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/PFS/forest/",name,"/CD8T_xCell.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  PFS_results = rbind(PFS_results,CD8T_xCell=result)
  
  expMarker$Class = ifelse(expMarker$CD8T_MCPcounter > mean(expMarker$CD8T_MCPcounter),"CD8T_MCPcounter High","CD8T_MCPcounter Low")
  fit <- survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/PFS/",name,"/CD8T_MCPcounter.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/PFS/forest/",name,"/CD8T_MCPcounter.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  PFS_results = rbind(PFS_results,CD8T_MCPcounter=result)
  
  expMarker$Class = ifelse(expMarker$CD8T_CIBERSORTx > mean(expMarker$CD8T_CIBERSORTx),"CD8T_CIBERSORTx High","CD8T_CIBERSORTx Low")
  fit <- survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/PFS/",name,"/CD8T_CIBERSORTx.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/PFS/forest/",name,"/CD8T_CIBERSORTx.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  PFS_results = rbind(PFS_results,CD8T_CIBERSORTx=result)
  
  expMarker$Class = ifelse(expMarker$Immunoscore_CIBERSORTx > mean(expMarker$Immunoscore_CIBERSORTx),"Immunoscore_CIBERSORTx High","Immunoscore_CIBERSORTx Low")
  fit <- survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/PFS/",name,"/Immunoscore_CIBERSORTx.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/PFS/forest/",name,"/Immunoscore_CIBERSORTx.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  PFS_results = rbind(PFS_results,Immunoscore_CIBERSORTx=result)
  
  expMarker$Class = ifelse(expMarker$Ecotype == "CE9","CE9","NOT_CE9")
  fit <- survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/PFS/",name,"/Ecotype.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  ## if there is no samples assign as CE9 EcoTyper, return NA
  if(length(levels(as.factor(expMarker$Class))) >1 ){
    res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
    hr = HR_95CI(res.cox)
    # if Inf exists, don't plot the forest
    if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
    ggsave(paste0("Results/Results_survival_analysis/PFS/forest/",name,"/Ecotype.png"),width = 12,height = 6, bg="white") }
    result = c(hr,p$pval)
    names(result) = c(names(hr),"p_value")
    PFS_results = rbind(PFS_results,Ecotype=result)
  }else{
    result = c(NA,NA,NA,NA)
    names(result) = c("Hazard Ratio","95% Upper CI","95% lower CI","p_value")
    PFS_results = rbind(PFS_results,Ecotype=result)
  }

  expMarker$Class = ifelse(expMarker$MFP == "IE","IE","NOT_IE")
  fit <- survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/PFS/",name,"/MFP.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/PFS/forest/",name,"/MFP.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  PFS_results = rbind(PFS_results,MFP=result)
  
  ####averager2ssGSEA: 9
  expMarker$Class = ifelse(expMarker$IFN_gamma_ssGSEA > mean(expMarker$IFN_gamma_ssGSEA),"IFN_gamma_ssGSEA High","IFN_gamma_ssGSEA Low")
  fit <- survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/PFS/",name,"/IFN_gamma_ssGSEA.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/PFS/forest/",name,"/IFN_gamma_ssGSEA.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  PFS_results = rbind(PFS_results,IFN_gamma_ssGSEA=result)
  
  expMarker$Class = ifelse(expMarker$Expanded_immune_gene_ssGSEA > mean(expMarker$Expanded_immune_gene_ssGSEA),"Expanded_immune_gene_ssGSEA High","Expanded_immune_gene_ssGSEA Low")
  fit <- survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/PFS/",name,"/Expanded_immune_gene_ssGSEA.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/PFS/forest/",name,"/Expanded_immune_gene_ssGSEA.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  PFS_results = rbind(PFS_results,Expanded_immune_gene_ssGSEA=result)
  
  expMarker$Class = ifelse(expMarker$T_cell_inflamed_GEP_ssGSEA > mean(expMarker$T_cell_inflamed_GEP_ssGSEA),"T_cell_inflamed_GEP_ssGSEA High","T_cell_inflamed_GEP_ssGSEA Low")
  fit <- survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/PFS/",name,"/T_cell_inflamed_GEP_ssGSEA.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/PFS/forest/",name,"/T_cell_inflamed_GEP_ssGSEA.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  PFS_results = rbind(PFS_results,T_cell_inflamed_GEP_ssGSEA=result)
  
  expMarker$Class = ifelse(expMarker$CRMA_ssGSEA > mean(expMarker$CRMA_ssGSEA),"CRMA_ssGSEA High","CRMA_ssGSEA Low")
  fit <- survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/PFS/",name,"/CRMA_ssGSEA.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/PFS/forest/",name,"/CRMA_ssGSEA.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  PFS_results = rbind(PFS_results,CRMA_ssGSEA=result)
  
  expMarker$Class = ifelse(expMarker$EMT_Stroma_core_ssGSEA > mean(expMarker$EMT_Stroma_core_ssGSEA),"EMT_Stroma_core_ssGSEA High","EMT_Stroma_core_ssGSEA Low")
  fit <- survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/PFS/",name,"/EMT_Stroma_core_ssGSEA.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/PFS/forest/",name,"/EMT_Stroma_core_ssGSEA.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  PFS_results = rbind(PFS_results,EMT_Stroma_core_ssGSEA=result)
  
  expMarker$Class = ifelse(expMarker$F_TBRS_ssGSEA > mean(expMarker$F_TBRS_ssGSEA),"F_TBRS_ssGSEA High","F_TBRS_ssGSEA Low")
  fit <- survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/PFS/",name,"/F_TBRS_ssGSEA.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/PFS/forest/",name,"/F_TBRS_ssGSEA.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  PFS_results = rbind(PFS_results,F_TBRS_ssGSEA=result)
  
  expMarker$Class = ifelse(expMarker$RiskScore_ssGSEA > mean(expMarker$RiskScore_ssGSEA),"RiskScore_ssGSEA High","RiskScore_ssGSEA Low")
  fit <- survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/PFS/",name,"/RiskScore_ssGSEA.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/PFS/forest/",name,"/RiskScore_ssGSEA.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  PFS_results = rbind(PFS_results,RiskScore_ssGSEA=result)
  
  expMarker$Class = ifelse(expMarker$TLS_score_ssGSEA > mean(expMarker$TLS_score_ssGSEA),"TLS_score_ssGSEA High","TLS_score_ssGSEA Low")
  fit <- survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/PFS/",name,"/TLS_score_ssGSEA.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  
    ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
    ggsave(paste0("Results/Results_survival_analysis/PFS/forest/",name,"/TLS_score_ssGSEA.png"),width = 12,height = 6, bg="white") 
  }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  PFS_results = rbind(PFS_results,TLS_score_ssGSEA=result)
  
  expMarker$Class = ifelse(expMarker$Renal_101_Immuno_ssGSEA > mean(expMarker$Renal_101_Immuno_ssGSEA),"Renal_101_Immuno_ssGSEA High","Renal_101_Immuno_ssGSEA Low")
  fit <- survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  png(paste0("Results/Results_survival_analysis/PFS/",name,"/Renal_101_Immuno_ssGSEA.png"),width = 600,height = 600, bg="white")
  plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
  print(plot); dev.off()
  p = surv_pvalue(fit,data=expMarker)
  res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
  hr = HR_95CI(res.cox)
  if(!is.element(Inf,hr) & !is.element(NA,hr)){  ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
  ggsave(paste0("Results/Results_survival_analysis/PFS/forest/",name,"/Renal_101_Immuno_ssGSEA.png"),width = 12,height = 6, bg="white") }
  result = c(hr,p$pval)
  names(result) = c(names(hr),"p_value")
  PFS_results = rbind(PFS_results,Renal_101_Immuno_ssGSEA=result)
  
  Dataset = rep(name,39+9)
  Biomarker = rownames(PFS_results)
  PFS_results = cbind(PFS_results,Dataset,Biomarker)
  
  return(PFS_results)
}

load("./ICB_data/Braun et al/Braun_2020.Rdata")
data1=Braun_data
name1="Braun_2020"
load("./ICB_data/Du et al/MGH_PRE_2021.RData")
data2=MGH_PRE_data
name2="MGH_PRE_2021"
load("./ICB_data/Du et al/MGH_ON_2021.RData")
data3=MGH_ON_data
name3="MGH_ON_2021"
load("./ICB_data/Gide et al/Gide_2019.RData")
data4=Gide_data
name4="Gide_2019"
load("./ICB_data/Hugo et al/Hugo_2016.RData")
data5=Hugo_data
name5="Hugo_2016"
load("./ICB_data/Lee et al/Lee_2020.RData")
data6=Lee_data
name6="Lee_2020"
load("./ICB_data/Liu et al/Liu_2019.Rdata")
data7=Liu_data
name7="Liu_2019"
load("./ICB_data/Mariathasan et al/Mariathasan_2018.Rdata")
data8=Mariathasan_data
name8="Mariathasan_2018"
load("./ICB_data/Miao et al/Miao_2018.Rdata")
data9=Miao_data
name9="Miao_2018"
load("./ICB_data/Nathanson et al/Nathanson_2017.Rdata")
data10=Nathanson_data
name10="Nathanson_2017"
load("./ICB_data/Riaz et al/Riaz_2017.RData")
data11=Riaz_data
name11="Riaz_2017"
load("./ICB_data/Snyder et al/Snyder_2017.Rdata")
data12=Snyder_data
name12="Snyder_2017"
load("./ICB_data/VanAllen et al/VanAllen_2015.RData")
data13=VanAllen_data
name13="VanAllen_2015"

OS_results = rbind(Biomarker_OS(data1,name1),
                   Biomarker_OS(data2,name2),
                   Biomarker_OS(data3,name3),
                   Biomarker_OS(data4,name4),
                   Biomarker_OS(data5,name5),
                   Biomarker_OS(data6,name6),
                   Biomarker_OS(data7,name7),
                   Biomarker_OS(data8,name8),
                   Biomarker_OS(data9,name9),
                   Biomarker_OS(data10,name10),
                   Biomarker_OS(data11,name11),
                   Biomarker_OS(data12,name12),
                   Biomarker_OS(data13,name13))

OS_results = as.data.frame(OS_results)

save(OS_results,file = "Results/OS_results.Rdata")

load("./ICB_data/Jung et al/Jung_2019.Rdata")
data14=Jung_data
name14="Jung_2019"
load("./ICB_data/Motzer et al/Motzer_2020.Rdata")
data15=Motzer_data
name15="Motzer_2020"

PFS_results = rbind(Biomarker_PFS(data1,name1),
                    Biomarker_PFS(data2,name2),
                    Biomarker_PFS(data3,name3),
                    Biomarker_PFS(data4,name4),
                    Biomarker_PFS(data7,name7),
                    Biomarker_PFS(data9,name9),
                    Biomarker_PFS(data11,name11),
                    Biomarker_PFS(data12,name12),
                    Biomarker_PFS(data13,name13),
                    Biomarker_PFS(data14,name14),
                    Biomarker_PFS(data15,name15))

PFS_results = as.data.frame(PFS_results)

save(PFS_results,file = "Results/PFS_results.Rdata")


###visualization of the results
rm(list = ls())

load("Results/OS_results.Rdata")
OS_results$HR = as.numeric(OS_results$`Hazard Ratio`)
OS_results$p_value = as.numeric(OS_results$p_value)
res = OS_results
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
source("./src/visualization.R")
p = dotplot(res, size="HR")
p
ggsave(p,filename = "./figures/os_survival_benchmark.pdf",width = 16,height = 16)


load("Results/PFS_results.Rdata")
PFS_results$HR = as.numeric(PFS_results$`Hazard Ratio`)
PFS_results$p_value = as.numeric(PFS_results$p_value)
res = PFS_results
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
source("./src/visualization.R")
p = dotplot(res, size="HR")
p
ggsave(p,filename = "./figures/pfs_survival_benchmark.pdf",width = 16,height = 16)


#############################################
##perform survival analysis on sub-datasets
load("ICB_data/Gide et al/Gide_COMBINE_2019.Rdata")
load("ICB_data/Gide et al/Gide_MONO_2019.Rdata")
load("ICB_data/Gide et al/Gide_ON_2019.Rdata")
load("ICB_data/Gide et al/Gide_PRE_2019.Rdata")

load("ICB_data/Riaz et al/Riaz_EXPOSURE_2017.Rdata")
load("ICB_data/Riaz et al/Riaz_NAIVE_2017.Rdata")
load("ICB_data/Riaz et al/Riaz_ON_2017.Rdata")
load("ICB_data/Riaz et al/Riaz_PRE_2017.RData")

load("ICB_data/Lee et al/Lee_ON_2020.Rdata")
load("ICB_data/Lee et al/Lee_PRE_2020.Rdata")
load("ICB_data/Liu et al/Liu_EXPOSURE_2019.Rdata")
load("ICB_data/Liu et al/Liu_NAIVE_2019.Rdata")

OS_results = rbind(Biomarker_OS(Gide_COMBINE_data,"Gide_COMBINE_2019"),
                   Biomarker_OS(Gide_MONO_data,"Gide_MONO_2019"),
                   Biomarker_OS(Gide_ON_data,"Gide_ON_2019"),
                   Biomarker_OS(Gide_PRE_data,"Gide_PRE_2019"),
                   Biomarker_OS(Riaz_EXPOSURE_data,"Riaz_EXPOSURE_2017"),
                   Biomarker_OS(Riaz_NAIVE_data,"Riaz_NAIVE_2017"),
                   Biomarker_OS(Riaz_ON_data,"Riaz_ON_2017"),
                   Biomarker_OS(Riaz_PRE_data,"Riaz_PRE_2017"),
                   Biomarker_OS(Lee_ON_data,"Lee_ON_2020"),
                   Biomarker_OS(Lee_PRE_data,"Lee_PRE_2020"),
                   Biomarker_OS(Liu_EXPOSURE_data,"Liu_EXPOSURE_2019"),
                   Biomarker_OS(Liu_NAIVE_data,"Liu_NAIVE_2019"))

OS_results = as.data.frame(OS_results)
save(OS_results,file = "Results/OS_results_subdatasets.Rdata")

PFS_results = rbind(Biomarker_PFS(Gide_COMBINE_data,"Gide_COMBINE_2019"),
                    Biomarker_PFS(Gide_MONO_data,"Gide_MONO_2019"),
                    Biomarker_PFS(Gide_ON_data,"Gide_ON_2019"),
                    Biomarker_PFS(Gide_PRE_data,"Gide_PRE_2019"),
                    Biomarker_PFS(Riaz_EXPOSURE_data,"Riaz_EXPOSURE_2017"),
                    # Biomarker_PFS(Riaz_NAIVE_data,"Riaz_NAIVE_2017"),
                    Biomarker_PFS(Riaz_ON_data,"Riaz_ON_2017"),
                    Biomarker_PFS(Riaz_PRE_data,"Riaz_PRE_2017"),
                    Biomarker_PFS(Liu_EXPOSURE_data,"Liu_EXPOSURE_2017"),
                    Biomarker_PFS(Liu_NAIVE_data,"Liu_NAIVE_2017"))

PFS_results = as.data.frame(PFS_results)

save(PFS_results,file = "Results/PFS_results_subdatasets.Rdata")


