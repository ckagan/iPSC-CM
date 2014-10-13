#Analysis of lactate purification

lac = read.table('Lactate2.txt', as.is=T, header=T)
key2 = c(rep("Lactate Purification", 13), rep("No Lactate Purification", 7))
key = lac$lact
pdf("Lactate.pdf")
library(gplots)
#Do not re-order anything!!
#boxplot2(lac$post.lactate~key, bottom=T, xlab = "Lactate", ylab = "% cTnT+", main = "Purity  With/Without Lactate")
boxplot2(lac$post.lactate~key2, bottom=T, ylab = "% cTnT+", main = "Purity  With/Without Lactate Purification")
dev.off()
laconly = lac[1:13,]
boxplot2(laconly$post.lactate~ laconly$days, bottom=T, xlab = "Days in Lactate Media", ylab = "% cTnT+", main = "Purity After Lactate Purification")
boxplot2(laconly$post.lactate~ laconly$Indiv, bottom=T, xlab = "Individual", ylab = "% cTnT+", main = "Purity After Lactate Purification")


summary(laconly)
#    lactate       days in lactate  
#Min.   :57.50   Min.   :2.000  
#1st Qu.:65.90   1st Qu.:4.000   
#Median :69.10   Median :6.000   
#Mean   :74.61   Mean   :5.077  
#3rd Qu.:87.00   3rd Qu.:6.000      
#Max.   :92.90   Max.   :8.000  

summary(lac[14:20,])
# no lactate
#Min.   : 1.60  
#1st Qu.: 4.85
#Median :16.40 
#Mean   :34.64 
#3rd Qu.:63.40 
#Max.   :88.00 

t.test(lac$post.lactate~key2)

#Welch Two Sample t-test

#data:  lac$post.lactate by key2
#t = 2.7273, df = 6.785, p-value = 0.03036
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  5.090353 74.839317
#sample estimates:
#  mean in group Lactate Purification mean in group No Lactate Purification 
#74.60769                              34.64286 