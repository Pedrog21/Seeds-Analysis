library(tidyverse)
library(dplyr)




income_evaluation <- read.csv("~/GitHub/ProjectAM/Data/adult.csv")


#1 Data treatment 

plot(education.num ~ education, data=income_evaluation) 
# since higher education level can be completely explained by education years 
# the variable higher education level is removed
income_evaluation = income_evaluation[,-c(4)]



for(i in 1:nrow(income_evaluation)){
  if(income_evaluation$native.country[i]== "?"){income_evaluation$native.country[i]="NA"}
}


