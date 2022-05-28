##Treatment classification:
## anti-PD-1/PD-L1 monotherapy 
## anti-CTLA4
## anti-PD-1/PD-L1 monotherapy with prior anti-CTLA4 monotherapy
## the combination of anti-PD-1 plus anti-CTLA-4 therapies
## Datasets used: Gide_data, Riaz_data, Liu_data

load("ICB_data/Gide et al/Gide_2019.RData")
load("ICB_data/Gide et al/Gide_MONO_2019.Rdata")
load("ICB_data/Gide et al/Gide_COMBINE_2019.Rdata")
load("ICB_data/Riaz et al/Riaz_2017.RData")
load("ICB_data/Riaz et al/Riaz_NAIVE_2017.Rdata")
load("ICB_data/Riaz et al/Riaz_EXPOSURE_2017.Rdata")
load("ICB_data/Liu et al/Liu_2019.Rdata")
load("ICB_data/Liu et al/Liu_NAIVE_2019.Rdata")
load("ICB_data/Liu et al/Liu_EXPOSURE_2019.Rdata")

data1 = Gide_data; name1 = "Gide_2019"
data2 = Gide_MONO_data; name2 = "Gide_MONO_2019"
data3 = Gide_COMBINE_data; name3 = "Gide_COMBINE_2019"
data4 = Riaz_data; name4 = "Riaz_2017"
data5 = Riaz_NAIVE_data; name5 = "Riaz_NAIVE_2017"
data6 = Riaz_EXPOSURE_data; name6 = "Riaz_EXPOSURE_2017"
data7 = Liu_data; name7 = "Liu_2019"
data8 = Liu_NAIVE_data; name8 = "Liu_NAIVE_2019"
data9 = Liu_EXPOSURE_data; name9 = "Liu_EXPOSURE_2019"


res = rbind(benchmark_result(data1,name1),
            benchmark_result(data2,name2),
            benchmark_result(data3,name3),
            benchmark_result(data4,name4),
            benchmark_result(data5,name5),
            benchmark_result(data6,name6),
            benchmark_result(data7,name7),
            benchmark_result(data8,name8),
            benchmark_result(data9,name9))


res$Dataset[which(res$Dataset=="Gide_2019")] = "Gide_Melanoma_90"
res$Dataset[which(res$Dataset=="Gide_MONO_2019")] = "Gide_MONO_Melanoma_50"
res$Dataset[which(res$Dataset=="Gide_COMBINE_2019")] = "Gide_COMBINE_Melanoma_40"
res$Dataset[which(res$Dataset=="Riaz_2017")] = "Riaz_Melanoma_103"
res$Dataset[which(res$Dataset=="Riaz_NAIVE_2017")] = "Riaz_NAIVE_Melanoma_45"
res$Dataset[which(res$Dataset=="Riaz_EXPOSURE_2017")] = "Riaz_EXPOSURE_Melanoma_58"
res$Dataset[which(res$Dataset=="Liu_2019")] = "Liu_Melanoma_121"
res$Dataset[which(res$Dataset=="Liu_NAIVE_2019")] = "Liu_NAIVE_Melanoma_74"
res$Dataset[which(res$Dataset=="Liu_EXPOSURE_2019")] = "Liu_EXPOSURE_Melanoma_47"


res$Dataset = as.factor(res$Dataset)
res$Biomarker = as.factor(res$Biomarker)
res$p_value = as.numeric(res$p_value)
res$AUC = as.numeric(res$AUC)

source("./src/visualization.R")
p = dotplot(res)
ggsave(p,filename = "./figures/treatment_benchmark.pdf",width = 10,height = 12)

save(list = ls(),file = "Results/treatment_benchmark_results.Rdata")

