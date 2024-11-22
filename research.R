whead(pc_scores)
ls(pc_scores)
colnames(pc_scores)

T1_NCR.NET = 9
T1_ED = 7


colnames(pc_scores)[grepl("T1_NCR.NET",colnames(pc_scores))]


c(colnames(pc_scores)[(which(endsWith(colnames(pc_scores), ".1"))-1)], colnames(pc_scores)[ncol(pc_scores)])

list <- data.frame("T1_NCR.NET" = 9,  
                 "T1_ED"=7  ,      
                 "T1_ET"=12 ,        
                 "T1Gd_NCR.NET"=9,
                 "T1Gd_ED"=8     ,
                 "T1Gd_ET"=9     ,  
                 "T2_NCR.NET"=16  ,  
                 "T2_ED"=12    ,    
                 "T2_ET"=13    ,    
                 "FLAIR_NCR.NET"=18,
                 "FLAIR_ED"=13   ,  
                 "FLAIR_ET"=17 )

library(dplyr)
library(tidyverse)

pc_scores
decrease <- pc_scores[1:61, 71:82] %>%
  apply(2, sd)
decrease

scores <- pc_scores[,(which(endsWith(colnames(pc_scores), ".1")))]
together <- cbind(scores, pathway.scores) 

cor <- vector("list", 12)
for(1 in 1:12)

  
laspply
sapply
setNames
grepl
apply
gsub
unique
names
colnames

colnames(pathway.scores)



library(tidyverse)
colSums(is.na(pathway.scores))

pca <- prcomp(pathway.scores, center=TRUE, scale.=TRUE) 

colnames(pca)
results <- -1*pca
install.packages("ggbiplot")
library(ggbiplot)

biplot(pca, scale=0)

normal <- scale(pathway.scores)
head(normal)

corr_matrix <- cor(normal)
install.packages("ggcorrplot")
library('ggcorrplot')
ggcorrplot(corr_matrix)


data.pca <- summary(princomp(corr_matrix))
install.packages("FactoMineR")
library("FactoMineR")

install.packages("factoextra")
library('factoextra')


fviz_eig(data.pca, addlabels = TRUE)
?fviz_eig

colnames(pc_scores)
scaled <- scale(pc_scores)
imaging_matrix <- cor(scaled)
ggcorrplot(imaging_matrix)

final.pca <- princomp(imaging_matrix)
summary(final.pca)

get_eig(final.pca)

fviz_eig(final.pca, addlabels = TRUE)

imagingpca <- prcomp(pc_scores, center=TRUE, scale.=TRUE)
scale()

X <- pathway.scores
pathways <- as.data.frame(pathway.scores)
pathways$survival <- Y
library(MASS)

#LASSO on the pathway scores data
install.packages("glmnet")
library(glmnet)
model1 <- cv.glmnet(x = X, y = Y, alpha = 1, nfolds = 10, family = gaussian, 
                    type.measure = "mse")
summary(model1)
options(repr.plot.width = 5, repr.plot.height = 5)
plot(model1)
which(coef(model1, s = 'lambda.min')[,"s1"]!=0)

