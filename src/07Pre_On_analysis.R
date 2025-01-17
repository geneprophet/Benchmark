##analysis the samples from pre-treatment and on-treatment
## Datasets used: NGH_ON_data, MGH_PRE_data, Gide_data, Riaz_data,Lee_data

load("ICB_data/Gide et al/Gide_2019.RData")
load("ICB_data/Gide et al/Gide_PRE_2019.Rdata")
load("ICB_data/Gide et al/Gide_ON_2019.Rdata")
load("ICB_data/Lee et al/Lee_2020.RData")
load("ICB_data/Lee et al/Lee_PRE_2020.Rdata")
load("ICB_data/Lee et al/Lee_ON_2020.Rdata")
load("ICB_data/Du et al/MGH_ON_2021.RData")
load("ICB_data/Du et al/MGH_PRE_2021.RData")
load("ICB_data/Riaz et al/Riaz_2017.RData")
load("ICB_data/Riaz et al/Riaz_PRE_2017.RData")
load("ICB_data/Riaz et al/Riaz_ON_2017.Rdata")

data1=Gide_data; name1="Gide_2019"
data2=Gide_PRE_data; name2="Gide_PRE_2019"
data3=Gide_ON_data; name3="Gide_ON_2019"
data4=Lee_data; name4="Lee_2020"
data5=Lee_PRE_data; name5="Lee_PRE_2020"
data6=Lee_ON_data; name6="Lee_ON_2020"
data7=MGH_PRE_data; name7="MGH_PRE_2021"
data8=MGH_ON_data; name8="MGH_ON_2021"
data9=Riaz_data; name9="Riaz_2017"
data10=Riaz_PRE_data; name10="Riaz_PRE_2017"
data11=Riaz_ON_data; name11="Riaz_ON_2017"



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
            benchmark_result(data11,name11))

res$Dataset[which(res$Dataset=="MGH_PRE_2021")] = "MGH_PRE_Melanoma_19"
res$Dataset[which(res$Dataset=="MGH_ON_2021")] = "MGH_ON_Melanoma_31"
res$Dataset[which(res$Dataset=="Gide_2019")] = "Gide_Melanoma_90"
res$Dataset[which(res$Dataset=="Gide_PRE_2019")] = "Gide_PRE_Melanoma_72"
res$Dataset[which(res$Dataset=="Gide_ON_2019")] = "Gide_ON_Melanoma_18"
res$Dataset[which(res$Dataset=="Riaz_2017")] = "Riaz_Melanoma_103"
res$Dataset[which(res$Dataset=="Riaz_PRE_2017")] = "Riaz_PRE_Melanoma_49"
res$Dataset[which(res$Dataset=="Riaz_ON_2017")] = "Riaz_ON_Melanoma_54"
res$Dataset[which(res$Dataset=="Lee_2020")] = "Lee_Melanoma_79"
res$Dataset[which(res$Dataset=="Lee_PRE_2020")] = "Lee_PRE_Melanoma_44"
res$Dataset[which(res$Dataset=="Lee_ON_2020")] = "Lee_ON_Melanoma_35"


res$Dataset = as.factor(res$Dataset)
res$Biomarker = as.factor(res$Biomarker)
res$p_value = as.numeric(res$p_value)
res$AUC = as.numeric(res$AUC)

source("./src/visualization.R")
p = dotplot(res)
ggsave(p,filename = "./figures/Pre_On_benchmark.pdf",width = 10,height = 12)

save(list = ls(),file = "Results/Pre_On_benchmark_results.Rdata")

