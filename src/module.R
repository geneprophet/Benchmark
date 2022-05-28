

##the usage of Rserve
# 
# library(Rserve)
# Rserve(port=6311,args="--no-save")
# library(RSclient)
# rsc <- RSconnect(port = 6311)
# RSshutdown(rsc)

library(ggplot2)
library(gridExtra)
library(grid)
library(ggpubr)
library(pROC)
library(survminer)
library(survival)
library(png)

HR_95CI <- function(x){ 
  x <- summary(x)
  HR <-signif(x$coef[2], digits=2)
  HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
  HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
  res<- c(HR,HR.confint.lower,HR.confint.upper)
  #names(res)<-c("Hazard Ratio","95% Upper CI","95% lower CI")
  names(res)<-c("HR","HR_Upper","HR_Lower")
  return(res)
}

kk_single = function(dataset,gene,uuid){
  my_comparisons = list(c("No_Response","Response"))
  file_path_dataset = paste0("./datasets/",dataset,".Rdata")
  load(file_path_dataset)
  data = get(dataset)
  result = c()
  if(is.element(gene,rownames(data$TPM))){
    expression = data$TPM[gene,]
    expMarker = merge(data$Samples,data.frame(Sample=names(expression),Marker=expression),by="Sample")
    expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
    g = ggboxplot(expMarker,x="Resp_NoResp",y="Marker",color = "Resp_NoResp",palette = "npg",add = "jitter") + 
      stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "two.sided"))
    ggsave(plot = g,filename = paste0("tmp/Results_wilcox_test/",uuid,".png"),width = 6,height = 6, bg="white")
    res = wilcox.test(Marker ~ Resp_NoResp, data = expMarker,alternative = "two.sided")
    auc = roc(Resp_NoResp ~ Marker, data = expMarker,auc = T)
    g = ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
    ggsave(plot = g,filename = paste0("tmp/Results_AUC/",uuid,".png"),width = 6,height = 6,bg = "white")
    result=c(res$p.value,as.numeric(auc$auc))
    names(result) = c("p_value","AUC")
    if(is.element("OS",colnames(data$Clinical))){
      OS_results <- c()
      expression = data$TPM[gene,]
      expMarker <- merge.data.frame(data$Clinical, data.frame(Sample=names(expression),Marker=expression),by="Sample",all.y = T)
      expMarker <- expMarker[which(!is.na(expMarker$OS)),]
      expMarker$OS <- as.numeric(expMarker$OS)
      expMarker$OS_CNSR <- as.numeric(expMarker$OS_CNSR)
      expMarker$Class <- ifelse(expMarker$Marker > mean(expMarker$Marker),"High","Low")
      fit = survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
      png(paste0("tmp/Results_os/",uuid,".png"),width = 600,height = 600, bg="white")
      plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
      print(plot);dev.off()
      p = surv_pvalue(fit,data=expMarker)
      res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
      hr = HR_95CI(res.cox)
      if(!is.element(Inf,hr) & !is.element(NA,hr)){ 
        plot = ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
        ggsave(plot=plot,filename = paste0("tmp/Results_os_forest/",uuid,".png"),width = 12,height = 6, bg="white") 
      }
      OS_results = c(hr,p$pval)
      names(OS_results) = c(paste0("OS_",names(hr)),"OS_p_value")
      result = c(result,OS_results)
    }
    if(is.element("PFS",colnames(data$Clinical))){
      PFS_results <- c()
      expression = data$TPM[gene,]
      expMarker <- merge.data.frame(data$Clinical, data.frame(Sample=names(expression),Marker=expression),by="Sample",all.y = T)
      expMarker <- expMarker[which(!is.na(expMarker$PFS)),]
      expMarker$PFS <- as.numeric(expMarker$PFS)
      expMarker$PFS_CNSR <- as.numeric(expMarker$PFS_CNSR)
      expMarker$Class <- ifelse(expMarker$Marker > mean(expMarker$Marker),"High","Low")
      fit = survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
      png(paste0("tmp/Results_pfs/",uuid,".png"),width = 600,height = 600, bg="white")
      plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
      print(plot);dev.off()
      p = surv_pvalue(fit,data=expMarker)
      res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
      hr = HR_95CI(res.cox)
      if(!is.element(Inf,hr) & !is.element(NA,hr)){ 
        plot = ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
        ggsave(plot=plot,filename = paste0("tmp/Results_pfs_forest/",uuid,".png"),width = 12,height = 6, bg="white") 
      }
      PFS_results = c(hr,p$pval)
      names(PFS_results) = c(paste0("PFS_",names(hr)),"PFS_p_value")
      result = c(result,PFS_results)
    }
    return(as.data.frame(t(result)))
  }else{
    return("The input gene is not included in the dataset.")
  }
}

# kk_single(dataset = "Gide_2019",gene = "CD274",uuid="test")

kk_set = function(dataset,geneset,uuid,type){
  
  my_comparisons = list(c("No_Response","Response"))
  file_path_dataset = paste0("./datasets/",dataset,".Rdata")
  load(file_path_dataset)
  data = get(dataset)
  result = c()
  
 # geneset = "CD3G,CD3E,CD8B,THEMIS,TRAT1,GRAP2,CD247,CD2,CD96,PRF1,CD6,IL7R,ITK,GPR18,EOMES,SIT1,NLRC3,CD244,KLRD1,SH2D1A,CCL5,XCL2,CST7,GFI1,KCNA3,PSTPIP1"
  
  signature = unlist(strsplit(geneset,split = ",",fixed = T))
  
  signature = intersect(signature,rownames(data$TPM))
  if(length(signature) >= 2){
    expression = data$TPM[signature,]
    if(type=="average"){
      expression = apply(expression, 2, mean)
    }else if(type=="sum"){
      expression = apply(expression, 2, sum)
    }else if(type=="ssGSEA"){
      library("GSVA")
      geneSet = list()
      geneSet[["signature"]] = signature
      gsva.es <- gsva(data$TPM, geneSet, method="ssgsea", verbose=T)
      expression = as.numeric(gsva.es)
      names(expression) = colnames(gsva.es)
    }
    expMarker = merge(data$Clinical,data.frame(Sample=names(expression),Signature=as.numeric(expression)),by="Sample")
    expMarker = expMarker[which(!is.na(expMarker$Resp_NoResp)),]
    g = ggboxplot(expMarker,x="Resp_NoResp",y="Signature",color = "Resp_NoResp",palette = "npg",add = "jitter") + 
      stat_compare_means(comparisons = my_comparisons,label = "p.signif",method.args = list(alternative = "two.sided"))
    ggsave(plot = g,filename = paste0("tmp/Results_wilcox_test/",uuid,".png"),width = 6,height = 6, bg="white")
    res = wilcox.test(Signature ~ Resp_NoResp, data = expMarker,alternative = "two.sided")
    auc = roc(Resp_NoResp ~ Signature, data = expMarker,auc = T)
    g = ggroc(auc,color="#4D96FF",linetype = 1.2,size= 1) + theme_minimal() + annotate("text", x=0.9, y=0.97, label=paste0("AUC = ",round(auc$auc,digits = 2)),size=6, fontface="bold") + xlab("Specificity") + ylab("Sensitivity")
    ggsave(plot = g,filename = paste0("tmp/Results_AUC/",uuid,".png"),width = 6,height = 6,bg = "white")
    result=c(res$p.value,as.numeric(auc$auc))
    names(result) = c("p_value","AUC")
    
    if(is.element("OS",colnames(data$Clinical))){
      OS_results <- c()
      expMarker = merge(data$Clinical,data.frame(Sample=names(expression),Signature=as.numeric(expression)),by="Sample",all.y = T)
      expMarker <- expMarker[which(!is.na(expMarker$OS)),]
      expMarker$OS <- as.numeric(expMarker$OS)
      expMarker$OS_CNSR <- as.numeric(expMarker$OS_CNSR)
      expMarker$Class <- ifelse(expMarker$Signature > mean(expMarker$Signature),"High","Low")
      fit = survfit(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
      png(paste0("tmp/Results_os/",uuid,".png"),width = 600,height = 600, bg="white")
      plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
      print(plot);dev.off()
      p = surv_pvalue(fit,data=expMarker)
      res.cox <- coxph(Surv(OS,OS_CNSR) ~ Class,data=expMarker)
      hr = HR_95CI(res.cox)
      if(!is.element(Inf,hr) & !is.element(NA,hr)){ 
        plot = ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
        ggsave(plot=plot,filename = paste0("tmp/Results_os_forest/",uuid,".png"),width = 12,height = 6, bg="white") 
      }
      OS_results = c(hr,p$pval)
      names(OS_results) = c(paste0("OS_",names(hr)),"OS_p_value")
      result = c(result,OS_results)
    }
    if(is.element("PFS",colnames(data$Clinical))){
      PFS_results <- c()
      expMarker = merge(data$Clinical,data.frame(Sample=names(expression),Signature=as.numeric(expression)),by="Sample",all.y = T)
      expMarker <- expMarker[which(!is.na(expMarker$PFS)),]
      expMarker$PFS <- as.numeric(expMarker$PFS)
      expMarker$PFS_CNSR <- as.numeric(expMarker$PFS_CNSR)
      expMarker$Class <- ifelse(expMarker$Signature > mean(expMarker$Signature),"High","Low")
      fit = survfit(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
      png(paste0("tmp/Results_pfs/",uuid,".png"),width = 600,height = 600, bg="white")
      plot = ggsurvplot(fit, data=expMarker,pval = TRUE,ggtheme = theme_minimal(),pval.method = TRUE,surv.median.line="hv",conf.int=T,risk.table=T)
      print(plot);dev.off()
      p = surv_pvalue(fit,data=expMarker)
      res.cox <- coxph(Surv(PFS,PFS_CNSR) ~ Class,data=expMarker)
      hr = HR_95CI(res.cox)
      if(!is.element(Inf,hr) & !is.element(NA,hr)){ 
        plot = ggforest(res.cox,data = expMarker,refLabel = 1,fontsize = 1, cpositions = c(0.00,0.06,0.4))
        ggsave(plot=plot,filename = paste0("tmp/Results_pfs_forest/",uuid,".png"),width = 12,height = 6, bg="white") 
      }
      PFS_results = c(hr,p$pval)
      names(PFS_results) = c(paste0("PFS_",names(hr)),"PFS_p_value")
      result = c(result,PFS_results)
    }
    return(as.data.frame(t(result)))
  }else{
    return("The number of input genes contained in the dataset is less than 2.")
  }
}

#kk_set(dataset = "Gide_2019",geneset = "CD3G,CD3E,CD8B,THEMIS,TRAT1,GRAP2,CD247,CD2,CD96,PRF1,CD6,IL7R,ITK,GPR18,EOMES,SIT1,NLRC3,CD244,KLRD1,SH2D1A,CCL5,XCL2,CST7,GFI1,KCNA3,PSTPIP1",uuid = "test",type = "sum")
