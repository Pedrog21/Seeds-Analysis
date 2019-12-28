#install.packages(c("mice","GGally","fastDummies","missMDA", "caret", "corrr", "dplyr", "e1071", "FactoMineR", "fitdistrplus", "Hmisc", "lsr", "naniar", "rcompanion", "tidyverse", "devtools"))
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
install_github("vqv/ggbiplot")


seeds = read.delim("~/GitHub/ProjectAM/Data/seeds_dataset.txt")
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


corMatrix = df_res %>%
  ggplot(aes(x,y,fill=assoc))+
  geom_tile()+
  # geom_text(aes(x,y,label=assoc))+
  scale_fill_gradient(low="red", high="yellow")+
  theme_classic()
corMatrix

seeds %>%
  mixed_assoc() %>%
  dplyr::select(x, y, assoc) %>%
  spread(y, assoc) %>%
  column_to_rownames("x") %>%
  as.matrix %>%
  as_cordf %>%
  network_plot()

##Decision Tree--------------------------------------------------------------------------------
df = seeds
smp_size <- floor(0.75 * nrow(df))
train_ind <- sample(seq_len(nrow(df)), size = smp_size)

train <- df[train_ind, ]
test <- df[-train_ind, ]
#Build the model
dtfit <- rpart(type ~ area + perimeter + compactnes + length + width + asymmetry + length_groove, method="class", data=train)
plotcp(dtfit)
summary(dtfit)
plot(dtfit, uniform=TRUE, 
     main="DT with pca data")
text(dtfit, use.n=TRUE, all=TRUE, cex=.8)


pred_dt <- predict(dtfit, test, type="class")
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
plotcp(dtfit.pca)
summary(dtfit.pca)
plot(dtfit.pca, uniform=TRUE, 
     main="DT with pca data")
text(dtfit.pca, use.n=TRUE, all=TRUE, cex=.8)


pred_dt.pca <- predict(dtfit.pca,test,type="class")
mtab<-table(pred_dt.pca, test[,4])
confusionMatrix(mtab)

##SVM/LDA---------------------------------------------------------------------------------------

#lda
df=seeds
smp_size <- floor(0.75 * nrow(df))
train_ind <- sample(seq_len(nrow(df)), size = smp_size)

train <- df[train_ind, ]
test <- df[-train_ind, ]
#Build the model
model2<-lda(x=train[,-c(8)],grouping = train[,8],prior = c(1/3,1/3,1/3),data=train,CV= FALSE)
#Summarize the model
summary(model2)
#Predict using the model
test1=test
predseeds= predict(model2,test[,-c(8)])
test1$pred_lda<-predict(model2,test[,-c(8)])$class
S#Accuracy of the model
mtab<-table(test1$pred_lda,test[,8])
confusionMatrix(mtab)
plot(model2)
newdata <- data.frame(type = test[,8], lda = predseeds$x)
library(ggplot2)
ggplot(newdata) + geom_point(aes(lda.LD1, lda.LD2, colour = type), size = 2.5)


##KNN


##PCA (Pedro)----------------------------------------------------------------------------------
seeds.pca <- prcomp(seeds[,1:7], center = TRUE, scale = TRUE)
summary(seeds.pca)

seeds.pca_df <- as.data.frame(seeds.pca$x[,c(1,2,3)]) #Dataset with only the first 3 principal components
seeds.pca_df$type <- seeds$type #Add the target variable in order to perform classification with this new dataset.

ggbiplot(seeds.pca)#lets you see how the data points relate to the axes



##Clustering (Ricardo) 3 clusters----------------------------------------------------------------