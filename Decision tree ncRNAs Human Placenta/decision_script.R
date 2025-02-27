require(dplyr)
require(ggplot2)
require(tibble)
require(rattle)

#reading the database:
library(readxl)
edgeR_DE_miRNAs <- read_excel("significant_survival_all_tumour_miRNAs_expression.xlsx")
View(edgeR_DE_miRNAs)

#Check if the class (response) you want to classify is as a factor, if not, fix it
str(edgeR_DE_miRNAs)
edgeR_DE_miRNAs$Sex<-as.factor(edgeR_DE_miRNAs$Sex)

######################################################################
#dividing the database into 70% for training and 30% for testing -> most common
######################################################################
 count(edgeR_DE_miRNAs$Sex=="Male")

# Randomly sample 
df_filtrado <- edgeR_DE_miRNAs[edgeR_DE_miRNAs$Sex=="Male", ]

# 106 samples
df_novo <- df_filtrado %>% sample_n(106, replace = FALSE)

## just female
female <- edgeR_DE_miRNAs[edgeR_DE_miRNAs$Sex=="Female", ]

## conc dataframes
df_conc <- rbind(df_novo, female)


######################################################################
#dividing the database into 70% for training and 30% for testing -> most common
######################################################################

#Sample function randomly sorts the database
#sample(1: number of samples)
#size: number of items to be chosen (in this case, 70%)
#replace: Is sampling with replacement?

#in this case size is *0.7, because I'm training 70%

tr_idx1 <- sample(1:nrow(df_conc), size = round(0.7*nrow(df_conc)), replace = FALSE)

treino1 <- df_conc[tr_idx1,]

teste1 <- df_conc[-tr_idx1,]

######################################################################
# Now we make a classification tree from the training data
######################################################################

require(rpart)
#my dataset here is training1 because the other part of the dataset I will use as test data in the final tree
arv_treino1 <- rpart(Sex~ . ,  
                     method="class",
                     data = treino1,
                     parms = list(split="gini"), 
                     cp=0.0001,
                     control = rpart.control(
                       minsplit = 1, # you can adjust according to your interest
                       minbucket = 1, # you can adjust according to your interest
                       maxdepth = 25 # maximum tree depth
                     ))
summary(arv_treino1) # detailed summary of splits

#minsplit e minbucket = 1 force the tree to grow, rpart has a default value for these values


require(rpart.plot)
require(partykit)
require(plotmo)

#two ways to plot the trees
rpart.plot(arv_treino1, 
           type = 0, 
           extra = 101, 
           box.palette = "GnBu",
           branch.lty=3,
           shadow.col = "gray",
           nn=TRUE,
           under = TRUE)

plot(as.party(arv_treino1),
     tp_args = list(fill = c("blue", "lightgray")),
     ip_args = list(fill = c("lightgreen")))

#detailed information about the decision tree
summary(arv_treino1)
print(arv_treino1)


#model summary table
printcp(arv_treino1)
#the final tree is on the last row
#we look for a tree with a lower x error to know the best cp value to prune the tree

plotcp(arv_treino1, upper = "splits")
#this graph also helps to make a decision regarding the value of cp


######################################################################
#If we need to prune the decision tree
######################################################################

#Gets the CP index with the lowest xerror:
opt1 <- which.min (arv_treino1$cptable [, "xerror"])

#get CP
cp1 <- arv_treino1$cptable[opt1, "CP"]
cp1
# this cp value has the lowest x error (you can double-check in the model summary table)

# We then prune the tree based on that CP value:
pruned1<- prune(arv_treino1,cp1)

# pruned tree

rpart.plot(pruned1, 
           type = 0, 
           extra = 101, 
           box.palette = "GnBu",
           branch.lty=3,
           shadow.col = "gray",
           nn=TRUE,
           under = TRUE)

plot(as.party(pruned1),
     tp_args = list(fill = c("blue", "lightgray")),
     ip_args = list(fill = c("lightgreen")))

#The table with the model summary
printcp(pruned1)

#note that the summary table indicates that there is no other tree with lesser error,
#so there is no need to prune

#important variables for the final model of the tree decision
df <- data.frame(imp = pruned1$variable.importance)
df2 <- df %>%
  tibble::rownames_to_column() %>%
  dplyr::rename("variable" = rowname) %>%
  dplyr::arrange(imp) %>%
  dplyr::mutate(variable = forcats::fct_inorder(variable))
ggplot2::ggplot(df2) +
  geom_col(aes(x = variable, y = imp),
           col = "black", show.legend = F) +
  coord_flip() +
  scale_fill_grey() +
  theme_bw()

######################################################################
#Evaluate the model
######################################################################

#Confusion matrix:

pred2<-predict(pruned1, teste1, type = "class")

conf2 <- table(pred2, teste1$Sex, dnn = c("Predicted", "Actual"))
print(conf2)

sum(diag(conf2))/sum(conf2)
# --% of the tissues in the test database were classified correctly.

#you will notice that if you run this code again, you will find different results
#that's why I do simulations of 10-100k trees to find the best set of 70 and 30%
# and thus have a better tree
#depending on what you find, we can do this step to find out how your trees behave
