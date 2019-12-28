library(FactoMineR)
library(fitdistrplus)
library(Hmisc)
library(MASS)
setwd('/Users/pedrogoncalves/Documents/OneDrive/MECD/AM/ProjectAM')
#Creating the Target variable
Target=c(1,2,2,2,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,
          1,2,2,2,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,
          1,2,2,2,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,
          1,2,2,2,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,
          1,2,2,2,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,
          1,2,2,2,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,
          1,2,2,2,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,
          1,2,2,2,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,
          1,2,2,2,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,
          1,2,2,2,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,
          1,2,2,2,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,
          1,2,2,2,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4)

PremierLeague <- read.csv("stats.csv")
PremierLeague$Target <- Target
PremierLeague <- PremierLeague[,-2]
PremierLeague <- PremierLeague[,-2]
PremierLeague$Target <- as.factor(PremierLeague$Target)

##Preliminary Data Analysis---------------------------------------------------------
#Summary Measures
lt2 <- c(1,40,41)
colMeans(PremierLeague[,-lt2])
table(PremierLeague$Target) #Few variables in Class 1 
describe(PremierLeague) #There are 6 variables with missing values
cor(PremierLeague$dispossessed, as.numeric(PremierLeague$Target), use = "complete.obs")
cor(PremierLeague$big_chance_missed, as.numeric(PremierLeague$Target), use = "complete.obs")
cor(PremierLeague$backward_pass, as.numeric(PremierLeague$Target), use = "complete.obs")
cor(PremierLeague$total_through_ball, as.numeric(PremierLeague$Target), use = "complete.obs")
cor(PremierLeague$head_clearance, as.numeric(PremierLeague$Target), use = "complete.obs")
cor(PremierLeague$saves, as.numeric(PremierLeague$Target), use = "complete.obs")
#Remove NA tries (variable dispossessed)
descdist(as.numeric(na.omit(PremierLeague$dispossessed)), boot=200) #Identify distr of complete obs
fitdistr(as.numeric(na.omit(PremierLeague$dispossessed)), "gamma")#Get the distr parameters
generate_obs <- rgamma(20,shape=32.457274300,scale=1/0.068624649)#Generate 20 obs from the above distr
generate_obs <- round(generate_obs)
#Plots
hist(PremierLeague$Target)
descdist(PremierLeague$dispossessed, boot=200)
pairs(PremierLeague,col=(PremierLeague[,41]))
#Correlation matrix
str(PremierLeague)
lt <- c(1,40)
corMatrix <- cor(PremierLeague[,-lt])
corMatrix <- as.data.frame(corMatrix)
corMatrix[,39]
write.csv(corMatrix,"/Users/pedrogoncalves/Documents/OneDrive/MECD/AM/ProjectAM\\corMatrix.csv", row.names = TRUE)

##Experimental PCA analysis--------------------------------------------------------
PCA(PremierLeague[,-lt], scale.unit = TRUE, ncp = 5, ind.sup = NULL, 
    quanti.sup = NULL, quali.sup = NULL, row.w = NULL, 
    col.w = NULL, graph = TRUE, axes = c(1,2))

#See how many variables explain more than 80% of the total variance
lt1 <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)
sum(pca$eig[lt1,1])/sum(pca$eig[,1])




#cena copiada so pa ver komo ficca
library(MASS)
library(e1071)
library(caret)

smp_size <- floor(0.75 * nrow(PremierLeague))
train_ind <- sample(seq_len(nrow(PremierLeague)), size = smp_size)

train <- PremierLeague[train_ind, ]
test <- PremierLeague[-train_ind, ]
#Build the model
model2<-lda(Target~goals+total_scoring_att+clean_sheet+interception+last_man_tackle+total_pass,data=train)
#Summarize the model
summary(model2)
#Predict using the model
test$pred_lda<-predict(model2,test)$class
#Accuracy of the model
mtab<-table(test$pred_lda,test$Target)
confusionMatrix(mtab)


#agora SVM 

library(kernlab)
library(SparseM)
smp_size <- floor(0.75 * nrow(PremierLeague))
train_ind <- sample(seq_len(nrow(PremierLeague)), size = smp_size)

train <- PremierLeague[train_ind, ]

test <- PremierLeague[-train_ind, ]
for(i in 1:ncol(test)){
  test[is.na(test[,i]), i] <- mean(test[,i], na.rm = TRUE)
}
#Build the model
model9<-svm(Target~goals+total_scoring_att+clean_sheet+interception+last_man_tackle+total_pass+touches,data=train,type = "C-classification")
#Summarize the model
summary(model9)
#Predict using the model
test$pred_svm<-predict(model9,test)
#Accuracy of the model
mtab2<-table(test$pred_svm,test$Target)
confusionMatrix(mtab2)
