library(ggplot2)
library(tidyverse)
## Importation des données en R et présentation du dataset

data <- readRDS("Datasets/data_vitC.Rds")
print(dim(data$prot))
print(dim(data$RNA))
print(dim(data$sample_info))
## 36 mices, 2392 protéines, 36066 gènes, 4 blocs (vitC+sex)

## Preprocessing

## Retirer les gènes avec des données manquantes
sum(is.na(data$RNA)) ## 0 donc pas de NA
## Garder uniquement les gènes avec de la variabilité 
threshold = 100
RNA_sd = apply(data$RNA, 2, sd)%>%
  sort()
genes_to_keep = RNA_sd[which(RNA_sd>threshold)]%>%
  names()
RNA = data$RNA[,genes_to_keep]
RNA[,colSums(RNA != 0) > 0]
RNA_log = apply(RNA, c(1,2),log)
RNA_scale = scale(RNA, center = TRUE, scale = TRUE)
dim(RNA_scale)
print(colnames(RNA_scale)[1:10])

threshold = 400000
prot_sd = apply(data$prot, 2, sd)%>%
  sort()
prot_to_keep = prot_sd[which(prot_sd>threshold)]%>%
  names()
prot = data$prot[,prot_to_keep]
prot[,colSums(prot != 0) > 0]
prot_log = apply(prot, c(1,2),log)
prot_scale = scale(prot, center = TRUE, scale = TRUE)
dim(prot_scale)
print(colnames(prot_scale)[1:10])
## D'abord log ou d'abord sd ?
## base de donnée pour associer gène à prot

## ACP
## gènes
pr = prcomp(RNA_scale)
pca = pr$x%>%
  as.data.frame()%>%
  mutate(id = row.names(.))%>%
  left_join(data$sample_info, by = join_by(id==sample))
ggplot(pca)+
  geom_point(aes(x = PC1, y = PC2, col = Sex))




# boxplot(head(data$RNA),col = data$sample_info, main = "boxplot mRNA")
# boxplot(head(t(data$RNA)), col = data$sample_info, main = "boxplot sample")

# datanett <- data$RNA[,colSums(data$RNA) > 10]

# datanett[datanett==0]<-0.1
# 
# coef.var <- function(x){
#   c.var = sd(x)/mean(x)
# }
# coef.mRNA <- as.numeric(apply(datanett,2,coef.var))
# hist(coef.mRNA, xlim = range(0,3))
# 
# data.filtered <- datanett[,abs(coef.mRNA) > 0.6]
# 
# data.logged <- log(data.filtered)
# data.scaled <- scale(data.filtered, center = TRUE, scale = TRUE)
