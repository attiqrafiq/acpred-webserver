#######set directory
#setwd('D:\\Peptide prediction\\Anticancer peptides\\Web server')
#######Load package
library(protr)
library(seqinr)
#library(RWeka)
#library(Interpol)
#library(caret)

library(kernlab)
library(e1071)
library(data.table)

####### Feature extraction for Training set
x <- read.fasta('acpred_training.fasta', seqtype="AA", as.string = TRUE)
label = read.csv("acpred_label.csv", header = TRUE) 
A <- x[(sapply(x, protcheck))]
m = length(A)

AAC <- t(sapply(A, extractAAC))

col = 20+ 2*4
APAAC  <- matrix(nrow = length(A), ncol = col)
for (i in 1:length(A)){
APAAC[i,] = extractAPAAC(A[[i]][1],lambda = 4, w = 0.01, customprops = NULL)
}

internal = data.frame(AAC,APAAC,Class = label)

Msvm =ksvm(Class ~ ., data = internal,kernel="rbfdot", cost = 2,gamma = 0.25,prob.model=TRUE)


####
saveRDS(Msvm, "Model.rds")

mod <- readRDS(file.path("data","Model.rds"))

#######################################################

#######Feature extraction for Test set
test <- read.fasta('acpred_example.fasta', seqtype="AA", as.string = TRUE)###read data
test <- test[(sapply(test, protcheck))]###check special symbol
AACtest <- t(sapply(test, extractAAC))
col = 20+ 2*4
APAACtest  <- matrix(nrow = length(test), ncol = col)
for (i in 1:length(test)){
APAACtest[i,] = extractAPAAC(test[[i]][1],lambda = 4, w = 0.01, customprops = NULL)
}


Dtest = data.frame(AACtest,APAACtest)
#######Predicting unknown sequences
data.frame(Prediction= predict(mod,Dtest),round(predict(mod,Dtest,type="prob"),3) )






