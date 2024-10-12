#Calculate the p-values of different makers at different positions
rm(list = ls())
library(ggpubr)
library(survminer)
library(survival)
library(glmnet)
setwd("F:/CRC-IHC/lasso cox/")
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

#Model prediction;
pred = predict(cvfit, newx = x, s = "lambda.min", type="response")
group = cbind(numeric(s1))
mv = median(pred)
group[pred<mv] = 1
group[pred>=mv] = 2


log1 = survdiff(mySurv ~ group)
p1 = pchisq(log1$chisq, 1, lower.tail=FALSE)
#print(p1)
res.cox = coxph(mySurv ~ group);
cox = summary(res.cox);

# plot KM curve
if (method == 'OS'){
  fit = survfit(mySurv ~ group)
  n1 = sum(group==1)
  leg1 = paste("Low risk(", n1, ")", sep = "")
  n2 = sum(group==2)
  leg2 = paste("High risk(", n2, ")", sep = "")
  
  png(filename = "figure/OS lasso_KM.png", width = 5.5, height = 5.5,
      units = "cm", res = 300, pointsize = 7)
  plot(fit, mark.time=TRUE, xlab = "Months", ylab = "Survival rate", lty = 1:2,
       col = 1:2, cex = 0.5)
  grid()
  legend(x = 94,y = 0.35, legend = c(leg1, leg2), lty = 1:2,
         col = 1:2, cex = 0.65)
  text(10, 0.1, paste("p=", formatC(p1, format="g", digits = 3), sep = ""),
       pos = 4, cex = 1)
  title("IHC Progression Score (OS)",cex.main = 1)
  dev.off()
}else{
  fit = survfit(mySurv ~ group)
  n1 = sum(group==1)
  leg1 = paste("Low risk(", n1, ")", sep = "")
  n2 = sum(group==2)
  leg2 = paste("High risk(", n2, ")", sep = "")
  
  png(filename = "figure/RFS lasso_KM.png", width = 5.5, height = 5.5,
      units = "cm", res = 300, pointsize = 7)
  plot(fit, mark.time=TRUE, xlab = "Months", ylab = "Relapse rate", lty = 1:2,
       col = 1:2, cex = 0.5)
  grid()
  legend(x = 94,y = 0.35, legend = c(leg1, leg2), lty = 1:2,
         col = 1:2, cex = 0.65)
  text(10, 0.1, paste("p=", formatC(p1, format="g", digits = 3), sep = ""),
       pos = 4, cex = 1)
  title("IHC Progression Score (RFS)",cex.main = 1)
  dev.off()
}



