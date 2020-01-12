#install.packages(c("mice","GGally","fastDummies","missMDA", "caret", "corrr", "dplyr", "e1071", "FactoMineR", "fitdistrplus", "Hmisc", "lsr", "naniar", "rcompanion", "tidyverse", "devtools", "fBasics", "e1071"))
library(tidyverse)
library(dplyr)
library(rpart)
library(naniar)
library(FactoMineR)
library(fitdistrplus)
library(Hmisc)
library(MASS)
library(rcompanion)
library(lsr)
require(corrr)
library(stats)
library(e1071)
library(gghighlight)
library(GGally)
ggpairs(seeds, aes(colour = Species, alpha = 0.4))
library(caret)
library(missMDA)
library(fastDummies)
library(mice)
library(reshape2)
library(ggplot2)
library(devtools)
library(ggbiplot)
library(fBasics)
library(corrplot)
library(numbers)
library(e1071)
library(psych)
#install_github("vqv/ggbiplot")


seeds = read.delim("~/Github/ProjectAM/Data/seeds_dataset.txt")
str(seeds)
seeds$type=as.factor(seeds$type)

#secondary function
mixed_assoc = function(df, cor_method="spearman", adjust_cramersv_bias=TRUE){
  df_comb = expand.grid(names(df), names(df),  stringsAsFactors = F) %>% set_names("X1", "X2")

  is_nominal = function(x) class(x) %in% c("factor", "character")
  # https://community.rstudio.com/t/why-is-purr-is-numeric-deprecated/3559
  # https://github.com/r-lib/rlang/issues/781
  is_numeric <- function(x) { is.integer(x) || is_double(x)}

  f = function(xName,yName) {
    x =  pull(df, xName)
    y =  pull(df, yName)

    result = if(is_nominal(x) && is_nominal(y)){
      # use bias corrected cramersV as described in https://rdrr.io/cran/rcompanion/man/cramerV.html
      cv = cramerV(as.character(x), as.character(y), bias.correct = adjust_cramersv_bias)
      data.frame(xName, yName, assoc=cv, type="cramersV")

    }else if(is_numeric(x) && is_numeric(y)){
      correlation = cor(x, y, method=cor_method, use="complete.obs")
      data.frame(xName, yName, assoc=correlation, type="correlation")

    }else if(is_numeric(x) && is_nominal(y)){
      # from https://stats.stackexchange.com/questions/119835/correlation-between-a-nominal-iv-and-a-continuous-dv-variable/124618#124618
      r_squared = summary(lm(x ~ y))$r.squared
      data.frame(xName, yName, assoc=sqrt(r_squared), type="anova")

    }else if(is_nominal(x) && is_numeric(y)){
      r_squared = summary(lm(y ~x))$r.squared
      data.frame(xName, yName, assoc=sqrt(r_squared), type="anova")

    }else {
      warning(paste("unmatched column type combination: ", class(x), class(y)))
    }

    # finally add complete obs number and ratio to table
    result %>% mutate(complete_obs_pairs=sum(!is.na(x) & !is.na(y)), complete_obs_ratio=complete_obs_pairs/length(x)) %>% rename(x=xName, y=yName)
  }

  # apply function to each variable combination
  map2_df(df_comb$X1, df_comb$X2, f)
}

##Preliminary Analysis (Frank) --- Standardization, Plots...
#summary of the data
summary(seeds)

#Scaling variables
seedsscale = cbind(scale(seeds[1:7]),seeds[8])

#conditional density estimate
ggplot(data = seeds, mapping = aes(x = area, fill= type,colour = type)) +
  geom_density(alpha = 0.5,position = "fill")


d <- melt(seeds[])


ggplot(d,aes(x = value,fill= type,colour = type)) +
     facet_wrap(~variable,scales = "free_x") +
    geom_freqpoly()


ggplot(d,aes(x = value,fill= type,colour = type)) +
  facet_wrap(~variable,scales = "free_x") +
  geom_area(stat = "bin")


#conditional density estimate
ggplot(d,aes(x = value, fill = type,colour = type)) +
  facet_wrap(~variable,scales = "free_x") +
  geom_density(alpha = 0.5,position = "fill")



ggplot(data=seeds,mapping = aes(x=type)) +
  geom_bar()        # we have the same amount of each type of plant

pairs(seeds,col=seeds[,8])
ggpairs(seeds, aes(colour = type, alpha = 0.4))
#correlation
df_res = mixed_assoc(seeds)

# plot results
#Correlation plots

corMatrix = df_res %>%
  ggplot(aes(x,y,fill=assoc))+
  geom_tile()+
  # geom_text(aes(x,y,label=assoc))+
  scale_fill_gradient(low="red", high="yellow")+
  theme_classic()
corMatrix

corrplot(cor(seeds[1:7]), type = "upper", order = "hclust", tl.col = "black", tl.srt = 45)

seeds %>%
  mixed_assoc() %>%
  dplyr::select(x, y, assoc) %>%
  spread(y, assoc) %>%
  column_to_rownames("x") %>%
  as.matrix %>%
  as_cordf %>%
  network_plot()


##LDA---------------------------------------------------------------------------------------

#lda
df=seeds
smp_size <- floor(0.75 * nrow(df))
train_ind <- sample(seq_len(nrow(df)), size = smp_size)

train <- df[train_ind, ]
test <- df[-train_ind, ]


#LDA
#Build the model
model2<-lda(x=train[,-c(8)],grouping = train[,8],prior = c(1/3,1/3,1/3),data=train,CV= FALSE)

#Summarize the model
summary(model2)

#Predict using the model
test1=test
predseeds= predict(model2,test[,-c(8)])
test1$pred_lda<-predict(model2,test[,-c(8)])$class

#Accuracy of the model
mtab<-table(test1$pred_lda,test[,8])

confusionMatrix(mtab)
plot(model2)
newdata <-data.frame(type = test[,8], lda = predseeds$x)
newdata2 <- data.frame(type = test1$pred_lda, lda = predseeds$x)
highlight_df = newdata[,1] != newdata2[,1]

ggplot(newdata) + geom_point(aes(lda.LD1, lda.LD2, colour = type), size = 2.5) +
  geom_point(data=newdata[highlight_df,], aes(x=lda.LD1,y=lda.LD2), color='tomato1',size=3,alpha=0.9) +
  geom_label(
    label="3", 
    data=newdata[highlight_df,], aes(x=lda.LD1,y=lda.LD2) ,
    label.padding = unit(0.1,"lines"), # Rectangle size around label
    label.size = 0.31,
    color = "black",
    fill="royalblue2",
    nudge_x = 0.18,alpha = 0.7)

ggbiplot(model2,groups = train$type, ellipse = TRUE, circle = TRUE)
#  gghighlight(newdata[,1] != newdata2[,1],use_group_by = FALSE)



##PCA (Pedro)----------------------------------------------------------------------------------
seeds.pca <- prcomp(seeds[,1:7], center = TRUE, scale = TRUE)
summary(seeds.pca)
seeds.pca_df <- as.data.frame(seeds.pca$x[,c(1,2,3)]) #Dataset with only the first 3 principal components
seeds.pca_df$type <- seeds$type

ggbiplot(seeds.pca)#lets you see how the data points relate to the axes

##Decision Tree--------------------------------------------------------------------------------
df = seeds
smp_size <- floor(0.75 * nrow(df))
train_ind <- sample(seq_len(nrow(df)), size = smp_size)

train <- df[train_ind, ]
test <- df[-train_ind, ]

#Building the model
dtfit <- rpart(type ~ area + perimeter + compactness + length + width + asymmetry + length_groove, method="class", data=train)
plotcp(dtfit) #plot of the cross validation step to choose the complexity parameter (cp)
summary(dtfit)
plot(dtfit, uniform=TRUE,
     main="DT without pca data")
text(dtfit, use.n=TRUE, all=TRUE, cex=.8) #plot of the tree


pred_dt <- predict(dtfit, test, type="class") #predictions
mtab<-table(pred_dt, test[,8])
confusionMatrix(mtab)

#DT with pca
df.pca = seeds.pca_df
smp_size <- floor(0.75 * nrow(df.pca))
train_ind <- sample(seq_len(nrow(df.pca)), size = smp_size)

train <- df.pca[train_ind, ]
test <- df.pca[-train_ind, ]

#Building the model
dtfit.pca <- rpart(type ~ PC1 + PC2 + PC3, method="class", data=train)
plotcp(dtfit.pca) #plot of the cross validation step to choose the complexity parameter (cp)
summary(dtfit.pca)
plot(dtfit.pca, uniform=TRUE,
     main="DT with pca data")
text(dtfit.pca, use.n=TRUE, all=TRUE, cex=.8) #plot of the tree


pred_dt.pca <- predict(dtfit.pca,test,type="class") #predictions
mtab<-table(pred_dt.pca, test[,4])
confusionMatrix(mtab)


###SUPPORT VECTOR MACHINES-------------------------------------------------------------------
##SVM without PCA
#Train/Test Split
df = seeds
smp_size <- floor(0.75 * nrow(df))
train_ind <- sample(seq_len(nrow(df)), size = smp_size)

train <- df[train_ind, ]
test <- df[-train_ind, ]

#Building the model
svmfit <- svm(type ~ area + perimeter + compactness + length + width + asymmetry + length_groove, data=train, scale=TRUE, kernel = "polynomial") #Polynomial kernel
svmfit.r <- svm(type ~ area + perimeter + compactness + length + width + asymmetry + length_groove, data=train, scale=TRUE, kernel = "radial") #Radial kernel
plot(cmdscale(dist(train[,-8])),col = as.integer(train[,8]),pch = c("o","+")[1:150 %in% svmfit$index + 1])


#Evaluation
pred_svm <- predict(svmfit,test) 
mtab<-table(pred_svm, test[,8])
confusionMatrix(mtab)

pred_svm.r <- predict(svmfit.r,test) 
mtab.r<-table(pred_svm.r, test[,8])
confusionMatrix(mtab.r)

##SVM with PCA
#Train/Test Split
df = seeds.pca_df
smp_size <- floor(0.75 * nrow(df))
train_ind <- sample(seq_len(nrow(df)), size = smp_size)

train <- df[train_ind, ]
test <- df[-train_ind, ]

#Building the model
svmfit <- svm(type ~ PC1 + PC2 + PC3, data=train, scale=TRUE, kernel = "polynomial") #Polynomial kernel
svmfit.r <- svm(type ~ PC1 + PC2 + PC3, data=train, scale=TRUE, kernel = "radial", method="C-classification") #Radial kernel
plot(cmdscale(dist(train[,-4])),col = as.integer(train[,4]),pch = c("o","+")[1:150 %in% svmfit$index + 1])

#Evaluation
pred_svm <- predict(svmfit,test) 
mtab<-table(pred_svm, test[,4])
confusionMatrix(mtab)

pred_svm.r <- predict(svmfit.r,test) 
mtab.r<-table(pred_svm.r, test[,4])
confusionMatrix(mtab.r)

###Clustering without PCA----------------------------------------------------------------------

#Function to correct clusters labels
label_correct = function(pred){
  count1 = count(pred[1:70])
  count2 = count(pred[71:140])
  count3 = count(pred[141:210])
  bad_labels = c(count1$x[which.max(count1$freq)], count2$x[which.max(count2$freq)], count3$x[which.max(count3$freq)])
  for (i in 1:210){
    if (pred[i]==bad_labels[1]){
      pred[i] = 1
    }
    else if(pred[i]==bad_labels[2]){
      pred[i] = 2
    }
    else{
      pred[i] = 3
    }
  }
  result =  pred
}

#Second Function to correct clusters labels

correct=function(cl){
  a=cl
  conf_matrix =confusionMatrix(table(true=seeds$type,pred=cl))
  h=conf_matrix$overall[1]
  for (i in 1:6){
    for(j in 1:length(cl)){
      if(cl[j]== mod(i-1,3)+1){
        cl[j]=mod(i,3)+1}
      else if(cl[j] == mod(i,3)+1){
        cl[j]=mod(i-1,3)+1}
    }
    
    conf_matrix =confusionMatrix(table(true=seeds$type,pred=cl))
    
    if (conf_matrix$overall[1]>h)
    { h= conf_matrix$overall[1]
    a=cl}
  }
  return(a)
} 




X = seeds[1:7]
X_scaled = scale(X)

#Hierarchical clustering with ward's distance
dist_X = dist(X)
ward_clust = hclust(dist_X, method="ward.D2")
cut_ward = cutree(ward_clust, k = 3)
cut_ward = label_correct(cut_ward)

#Hierarchical clustering with ward's distance and scaled variables
dist_X_scaled = dist(X_scaled)
ward_clust_scaled = hclust(dist_X_scaled, method="ward.D2")
cut_ward_scaled = cutree(ward_clust_scaled, k = 3)
cut_ward_scaled = label_correct(cut_ward_scaled)

#Confusion matrix - ward's clustering
conf_matrix_ward = confusionMatrix(as.factor(cut_ward), seeds$type)
conf_matrix_ward

conf_matrix_ward_scaled = confusionMatrix(as.factor(cut_ward_scaled), seeds$type)
conf_matrix_ward_scaled

#KMeans
kmeans_clust = kmeans(X, 3, iter.max=50, nstart=10)
classes_kmeans = fitted(kmeans_clust, method="classes")
classes_kmeans = label_correct(classes_kmeans)

#KMeans with scaled variables
kmeans_clust_scaled = kmeans(X_scaled, 3, iter.max=50, nstart=10)
classes_kmeans_scaled = fitted(kmeans_clust_scaled, method="classes")
classes_kmeans_scaled = label_correct(classes_kmeans_scaled)

#Confusion matrix - kmeans clustering
conf_matrix_kmeans = confusionMatrix(as.factor(classes_kmeans), seeds$type)
conf_matrix_kmeans

conf_matrix_kmeans_scaled = confusionMatrix(as.factor(classes_kmeans_scaled), seeds$type)
conf_matrix_kmeans_scaled


#Clustering for the first 3 Principal Components-------------------------------------------------------------------
X_pca3 = seeds.pca_df[1:3]
X_pca3_scaled = scale(X)

#Hierarchical clustering with ward's distance
dist_X_pca3 = dist(X_pca3)
ward_clust_pca3 = hclust(dist_X_pca3, method="ward.D2")
cut_ward_pca3 = cutree(ward_clust_pca3, k = 3)
cut_ward_pca3 = label_correct(cut_ward_pca3)

#Hierarchical clustering with ward's distance and scaled variables
dist_X_pca3_scaled = dist(X_pca3_scaled)
ward_clust_pca3_scaled = hclust(dist_X_pca3_scaled, method="ward.D2")
cut_ward_pca3_scaled = cutree(ward_clust_pca3_scaled, k = 3)
cut_ward_pca3_scaled = label_correct(cut_ward_pca3_scaled)

#Confusion matrix - ward's clustering
conf_matrix_ward_pca3 = confusionMatrix(as.factor(cut_ward_pca3), seeds.pca_df$type)
conf_matrix_ward_pca3

conf_matrix_ward_pca3_scaled = confusionMatrix(as.factor(cut_ward_pca3_scaled), seeds.pca_df$type)
conf_matrix_ward_pca3_scaled

#KMeans
kmeans_clust_pca3 = kmeans(X_pca3, 3, iter.max=50, nstart=10)
classes_kmeans_pca3 = fitted(kmeans_clust_pca3, method="classes")
classes_kmeans_pca3 = label_correct(classes_kmeans_pca3)

#KMeans with scaled variables
kmeans_clust_pca3_scaled = kmeans(X_pca3_scaled, 3, iter.max=50, nstart=10)
classes_kmeans_pca3_scaled = fitted(kmeans_clust_pca3_scaled, method="classes")
classes_kmeans_pca3_scaled = label_correct(classes_kmeans_pca3_scaled)

#Confusion matrix - kmeans clustering
conf_matrix_kmeans_pca3 = confusionMatrix(as.factor(classes_kmeans_pca3), seeds.pca_df$type)
conf_matrix_kmeans_pca3

conf_matrix_kmeans_pca3_scaled = confusionMatrix(as.factor(classes_kmeans_pca3_scaled),seeds.pca_df$type)
conf_matrix_kmeans_pca3_scaled


#Repeat the classification Study now using the partition clusters

# for the clustering we will use the hierachical clustering using the ward method and Euclidean distance


#1.Decision Tree

seedsnewlabel = seeds
seedsnewlabel$type = cut_ward_scaled
seedsnewlabel$type=as.factor(seedsnewlabel$type)

df = seedsnewlabel

##PCA (Pedro)----------------------------------------------------------------------------------

seeds.pca_df$type <- seedsnewlabel$type

##___________________________________
train <- df[train_ind, ]
test <- df[-train_ind, ]

#Building the model
dtfit <- rpart(type ~ area + perimeter + compactness + length + width + asymmetry + length_groove, method="class", data=train)
plotcp(dtfit) #plot of the cross validation step to choose the complexity parameter (cp)
summary(dtfit)
plot(dtfit, uniform=TRUE,
     main="DT without pca data")
text(dtfit, use.n=TRUE, all=TRUE, cex=.8) #plot of the tree


pred_dt <- predict(dtfit, test, type="class") #predictions
mtab<-table(pred_dt, test[,8])
confusionMatrix(mtab)

#DT with pca
df.pca = seeds.pca_df
smp_size <- floor(0.75 * nrow(df.pca))
train_ind <- sample(seq_len(nrow(df.pca)), size = smp_size)

train <- df.pca[train_ind, ]
test <- df.pca[-train_ind, ]
#Build the model
dtfit.pca <- rpart(type ~ PC1 + PC2 + PC3, method="class", data=train)
plotcp(dtfit.pca) #plot of the cross validation step to choose the complexity parameter (cp)
summary(dtfit.pca)
plot(dtfit.pca, uniform=TRUE,
     main="DT with pca data")
text(dtfit.pca, use.n=TRUE, all=TRUE, cex=.8) #plot of the tree


pred_dt.pca <- predict(dtfit.pca,test,type="class") #predictions
mtab<-table(pred_dt.pca, test[,4])
confusionMatrix(mtab)

# we see worst values for DT with PCA

##2.LDA

smp_size <- floor(0.75 * nrow(df))
train_ind <- sample(seq_len(nrow(df)), size = smp_size)

train <- df[train_ind, ]
test <- df[-train_ind, ]





#LDA
#Build the model
model2<-lda(x=train[,-c(8)],grouping = train[,8],data=train,CV= FALSE)

#Summarize the model
summary(model2)

#Predict using the model
test1=test
predseeds= predict(model2,test[,-c(8)])
test1$pred_lda<-predict(model2,test[,-c(8)])$class

#Accuracy of the model
mtab<-table(test1$pred_lda,test[,8])
confusionMatrix(mtab)
plot(model2)
newdata <- data.frame(type = test[,8], lda = predseeds$x)
library(ggplot2)
ggplot(newdata) + geom_point(aes(lda.LD1, lda.LD2, colour = type), size = 2.5)

#worst values with LDA

#QDA
smp_size <- floor(0.75 * nrow(df))
train_ind <- sample(seq_len(nrow(df)), size = smp_size)

train <- df[train_ind, ]
test <- df[-train_ind, ]





#LDA
#Build the model
model2<-qda(x=train[,-c(8)],grouping = train[,8],data=train,CV= FALSE)

#Summarize the model
summary(model2)

#Predict using the model
test1=test
predseeds= predict(model2,test[,-c(8)])
test1$pred_lda<-predict(model2,test[,-c(8)])$class

#Accuracy of the model
mtab<-table(test1$pred_lda,test[,8])
confusionMatrix(mtab)
#plot(model2)



#3.Support Vector Machines

seedsnewlabel = seeds
seedsnewlabel$type = cut_ward_scaled
seedsnewlabel$type=as.factor(seedsnewlabel$type)

df = seedsnewlabel
smp_size <- floor(0.75 * nrow(df))
train_ind <- sample(seq_len(nrow(df)), size = smp_size)

train <- df[train_ind, ]
test <- df[-train_ind, ]


#Building the model
svmfit <- svm(type ~ area + perimeter + compactness + length + width + asymmetry + length_groove, data=train, scale=TRUE, kernel = "polynomial") #Polynomial kernel
svmfit.r <- svm(type ~ area + perimeter + compactness + length + width + asymmetry + length_groove, data=train, scale=TRUE, kernel = "radial") #Radial kernel
plot(cmdscale(dist(train[,-8])),col = as.integer(train[,8]),pch = c("o","+")[1:150 %in% svmfit$index + 1])


#Evaluation
pred_svm <- predict(svmfit,test)
mtab<-table(pred_svm, test[,8])
confusionMatrix(mtab)

pred_svm.r <- predict(svmfit.r,test) 
mtab.r<-table(pred_svm.r, test[,8])
confusionMatrix(mtab.r)
