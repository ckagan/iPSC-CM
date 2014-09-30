key = lac$lact
lac = read.table('Lactate2.txt', as.is=T, header=T)
pdf("Lactate.pdf")
boxplot2(lac$post.lactate~key, bottom=T, xlab = "Lactate", ylab = "% cTnT+", main = "Purity Pre/Post-Lactate")
laconly = lac[1:13,]
boxplot2(laconly$post.lactate~ laconly$days, bottom=T, xlab = "Days in Lactate Media", ylab = "% cTnT+", main = "Purity After Lactate")
boxplot2(laconly$post.lactate~ laconly$Indiv, bottom=T, xlab = "Individual", ylab = "% cTnT+", main = "Purity After Lactate")
dev.off()
