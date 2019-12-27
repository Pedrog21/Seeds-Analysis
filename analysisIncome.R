#install.packages(c("mice","fastDummies","missMDA", "caret", "corrr", "dplyr", "e1071", "FactoMineR", "fitdistrplus", "Hmisc", "lsr", "naniar", "rcompanion", "tidyverse"))
library(tidyverse)
library(dplyr)
library(naniar)
library(FactoMineR)
library(fitdistrplus)
library(Hmisc)
library(MASS)
library(rcompanion)
library(lsr)
require(corrr)
library(MASS)
library(e1071)
library(caret)
library(missMDA)
library(fastDummies)
library(mice)
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



income_evaluation <- read.csv("~/GitHub/ProjectAM/Data/adult.csv")



#1 Data treatment 

plot(education.num ~ education, data=income_evaluation) 
# since higher education level can be completely explained by education years 
# the variable higher education level is removed
income_evaluation = income_evaluation[,-c(4)]
# replace all the ? with NA for simplicity
for(i in 1:nrow(income_evaluation)){
  for (j in 1:ncol(income_evaluation)) {
    
  if(income_evaluation[i,j]== "?"){income_evaluation[i,j]="NA"}
}}
income_evaluation = na.omit(income_evaluation)

#we decided to remove the finalweight variable, because the finalweight is a value of the numbers of people each
#row represents. Since the dataset is big enough we can assume that the ratio of people from different races is present


income_evaluation = income_evaluation[,-c(3)]


#we decided to join the two variables capital gains and losses into one variable that takes negative (loss) or positive values 
income_evaluation$capital.gain = income_evaluation$capital.gain-income_evaluation$capital.loss
income_evaluation = income_evaluation[,-c(10)]

describe(income_evaluation) # some variables have missing data

#Now we look at the correlation matrix
df=income_evaluation 

df_res = mixed_assoc(df)

# plot results


corMatrix = df_res %>%
  ggplot(aes(x,y,fill=assoc))+
  geom_tile()+
 # geom_text(aes(x,y,label=assoc))+
  scale_fill_gradient(low="red", high="yellow")+
  theme_classic()
corMatrix

df %>%
  
  mixed_assoc() %>%
  dplyr::select(x, y, assoc) %>%
  spread(y, assoc) %>%
  column_to_rownames("x") %>%
  as.matrix %>%
  as_cordf %>%
  network_plot()


#the Higher correlation is between sex and relationship, and marital_status and age


smp_size <- floor(0.75 * nrow(df))
train_ind <- sample(seq_len(nrow(df)), size = smp_size)

train <- df[train_ind, ]
test <- df[-train_ind, ]
#Build the model
model2<-lda(income~workclass + education.num + marital.status + race + capital.gain + hours.per.week,data=train)
#Summarize the model
summary(model2)
#Predict using the model
test$pred_lda<-predict(model2,test)$class
#Accuracy of the model
mtab<-table(test$pred_lda,test$income)
confusionMatrix(mtab)

#Creating dummy variables for categorical data
dworkclass = dummy_cols(df$workclass)
dmarital_Status = dummy_cols(df$marital.status)
doccupation = dummy_cols(df$occupation)
drelationship = dummy_cols(df$relationship)
drace = dummy_cols(df$race)
dsex = dummy_cols(df$sex)
dnative_country = dummy_cols(df$native.country)
dincome = dummy_cols(df$income)

#Data Frame with dummy variables
dumdf = cbind(df[])

estn= estim_ncpFAMD(df,ncp.min=3, ncp.max= 12 )
res.famd <- FAMD(df)

## Gower dissimilarity matrix
dmm=df
dmm=dmm[complete.cases(dmm),]
dmm[,1]=scale(dmm$age)
dmm$education.num = scale(dmm$education.num)
dmm$capital.gain = scale(dmm$capital.gain)
dmm$hours.per.week = scale(dmm$hours.per.week)
dmn = std(df,include.fac = FALSE )
daisy(df)
###klustering
h2o.init()
h2o.df <- as.h2o(df[,-c(12)])
h2o_kmeans <- h2o.kmeans(training_frame = h2o.df,k=7, estimate_k = TRUE,nfolds = 5,keep_cross_validation_predictions = TRUE)
summary(h2o_kmeans)
dmn = std(df,include.fac = FALSE )
daisy(df)