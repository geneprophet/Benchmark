

save_file = function(name){
  load(paste0("datasets/",name,'.Rdata'))
  data = get(name)
  write.table(data$Clinical,file = paste0("datasets/",name,"_clinical.txt"),sep = "\t", quote = F, row.names = F)
  write.table(data$TPM,file = paste0("datasets/",name,"_expression.txt"),sep = "\t", quote = F, row.names = F)
  write.table(data$Landscape,file = paste0("datasets/",name,"_landscape.txt"),sep = "\t", quote = F, row.names = F)
}
library(readxl)
dataset <- read_excel("statistical.xlsx")
dataset$name
for (i in dataset$name) {
  print(i)
  save_file(i)
}

