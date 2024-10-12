## Calculate p-values for different positions and makers£»
rm(list = ls())
library(ggpubr)
library(survminer)
library(survival)
library(glmnet)
set.seed(1)
setwd("F:/CRC-IHC/univariate logranktest completed xtile")
list.files = list.files(pattern="completed_maker*")


method = 'OS';
#1:OS,2:RFS;

maker<-list('CD3', 'CD8', 'CD45RO','Granzyme_B','CD4','CD57','CD68','S100','Tryptase',
            'Foxp3','CD20','HLA_DR','FasL','IL_17','Fas')

matrix = matrix(nrow=15,ncol=8);
matrix1 = matrix(nrow=15,ncol=8);

#Read the cutoff value for xtile;
if (method == 'OS'){
  xtile = read.csv("OS xtile cutoff.csv")
}else{
  xtile = read.csv("RFS xtile cutoff.csv")
}
xtile = xtile[,-1]


#Group according to the threshold of xtile;
for (num in 1:length(list.files)){
  
  data = read.csv(list.files[num])
  m = nrow(data);
  n = ncol(data);
  
  for (i in 5:n){
    threshold = xtile[i-4, num]
    data[,i] = ifelse(data[,i] > threshold, '2', '1');
  }
  
#Calculate p-value using logrank test£»
  if (method == 'OS'){
    mySurv = Surv(data[, 1], data[, 2]);
  }else{
    mySurv = Surv(data[, 3], data[, 4]);
  }

  for (i in 5:n){
    group = data[, i];
    log1 = survdiff(mySurv ~ group);
    p = pchisq(log1$chisq, 1, lower.tail=FALSE);
    matrix[i-4, num] = p;
    
    # plot KM curve
    fit = survfit(mySurv ~ group)
    n1 = sum(group==1)
    leg1 = paste("Low group(", n1, ")", sep = "")
    n2 = sum(group==2)
    leg2 = paste("High group(", n2, ")", sep = "")
    
    if (method == 'OS'){
      name=sprintf('./KM_OS/KMCurve_%s_%s.png',num,maker[i-4])}
    else{
      name=sprintf('./KM_RFS/KMCurve_%s_%s.png',num,maker[i-4])
      }
    png(filename = name, width = 5.5, height = 5.5,
        units = "cm", res = 300, pointsize = 7)
    plot(fit, mark.time=TRUE, xlab = "Months", ylab = "Survival", lty = 1:2,
         col = 1:2, cex = 0.5, main=maker[i-4])
    grid()
    legend(x = "topright", legend = c(leg1, leg2), lty = 1:2,
           col = 1:2, cex = 0.65)
    text(10, 0.1, paste("p=", formatC(p, format="g", digits = 3), sep = ""),
         pos = 4, cex = 1)
    dev.off()
    
    # Determine whether high-value group has better survival outcome than low-value group
    coxph.fit <- coxph(mySurv ~ group, method="breslow")
    # hr <- exp(coef(coxph.fit))
    mysum = summary(coxph.fit)
    hr = mysum$conf.int[1]
    if(hr < 1){
      matrix1[i-4, num] = 1
    }else{
      matrix1[i-4, num] = 0
    }
  }
}
##Display the number and proportion of significant makers;
# matrix = ifelse(matrix < 0.05, '1', '0');
# sum(matrix==1)
# print(sum(matrix==1)/120)


if (method == 'OS'){
  write.csv(matrix, "OS p.csv", row.names = FALSE)
}else{
  write.csv(matrix, "RFS p.csv", row.names = FALSE)
}
if (method == 'OS'){
  write.csv(matrix1, "OS_cor.csv", row.names = FALSE)
}else{
  write.csv(matrix1, "RFS_cor.csv", row.names = FALSE)
}


