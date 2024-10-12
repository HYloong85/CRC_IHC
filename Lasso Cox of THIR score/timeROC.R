##Calculate the p-values of different positions and makers, and plot the timeROC curve
rm(list = ls())
library(ggpubr)
library(survminer)
library(survival)
library(glmnet)

library(survivalROC)
library(timeROC)
set.seed(1)
setwd("./")
list.files = list.files(pattern="all fea data*")

method = 'OS';
#1:OS,2:RFS;


data = read.csv(list.files[1])

data = data.matrix(data)
s1 = nrow(data)
s2 = ncol(data)
if (method == 'OS'){
  mySurv = Surv(data[, 1], data[, 2]);
}else{
  mySurv = Surv(data[, 3], data[, 4]);
}
x = data[, 5:s2]

# 10-fold CV
set.seed(1)
cvfit = cv.glmnet(x, mySurv, family = "cox", type.measure = "C")

beta = coef(cvfit, s = "lambda.min")

pred = predict(cvfit, newx = x, s = "lambda.min", type="response")
data = data[,-(6:19)]
data[,5] = pred

#Take the median and calculate the p-value;
if (method == 'OS'){
  mySurv = Surv(data[,1], data[,2]);
}else{
  mySurv = Surv(data[,3], data[,4]);
}
group = as.numeric(data[,5])
mv = median(group);
group = ifelse(group >= mv, '2', '1');
log1 = survdiff(mySurv ~ group);
p1 = pchisq(log1$chisq, 1, lower.tail=FALSE);
print(p1)


# plot ROC curve
if (method == 'OS'){
  quantile = quantile(data[,1])
  sROC <-timeROC(T=data[,1],delta = data[,2],marker = data[,5],
                 cause = 1,times = c(quantile[2],quantile[3],quantile[4]),ROC=T)
  png(filename = "figure/OS lasso_timeroc.png", width = 8, height = 8,
      units = "cm", res = 300, pointsize = 7)
}else{
  quantile = quantile(data[,3])
  sROC <-timeROC(T=data[,3],delta = data[,4],marker = data[,5],
                 cause = 1,times = c(quantile[2],quantile[3],quantile[4]),ROC=T)
  png(filename = "figure/RFS lasso_timeroc.png", width = 8, height = 8,
      units = "cm", res = 300, pointsize = 7)
}
par(mar= c(5,5,5,5),cex.lab=1.5,cex.axis=1.2) #Set graphic boundaries
plot(sROC,time=quantile[2],col="cyan3",title=F,lwd=2) #1 year ROC
plot(sROC,time=quantile[3],col="green",add=T,title=F,lwd=2) #3 year ROC
plot(sROC,time=quantile[4],col="salmon",add=T,title=F,lwd=2) #5 year ROC

if (method == 'OS'){
  
  legend("bottomright", # Set the position coordinates of the legend here: horizontal axis 0, vertical axis 1. 
         c(paste0("AUC 25%-OS  ", round(sROC$AUC[1], 2)),
           paste0("AUC 50%-OS  ", round(sROC$AUC[2], 2)),
           paste0("AUC 75%-OS  ", round(sROC$AUC[3], 2))),
         col=c("cyan3","green","salmon"),lwd=2,cex=1.2,bty="n")
  title("Time-ROC (OS)",cex.main = 1.5)
  
}else{
  
  legend("bottomright", 
         c(paste0("AUC 25%-RFS  ", round(sROC$AUC[1], 2)),
           paste0("AUC 50%-RFS  ", round(sROC$AUC[2], 2)),
           paste0("AUC 75%-RFS  ", round(sROC$AUC[3], 2))),
         col=c("cyan3","green","salmon"),lwd=2,cex=1.2,bty="n")
  title("Time-ROC (RFS)",cex.main = 1.5)
}
dev.off()
print(sROC$AUC)


if (method == 'OS'){
  ROC3 <- timeROC(T=data[,1], delta=data[,2],  marker=data[,5], cause=1, weighting="marginal", times = c(quantile[2],quantile[3],quantile[4]), iid=TRUE)
  png(filename = "OS timeROC all fea lasso coef.png", width = 8, height = 8,
      units = "cm", res = 300, pointsize = 7)
}else{
  ROC3 <- timeROC(T=data[,3], delta=data[,4],  marker=data[,5], cause=1, weighting="marginal", times = c(quantile[2],quantile[3],quantile[4]), iid=TRUE)
  png(filename = "RFS timeROC all fea lasso coef.png", width = 8, height = 8,
      units = "cm", res = 300, pointsize = 7)
}
plotAUCcurve(ROC3, conf.int=TRUE, conf.band=TRUE,col = "black")
dev.off()


