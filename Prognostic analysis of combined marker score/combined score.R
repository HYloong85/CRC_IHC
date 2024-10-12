##Calculate the score, with+1 point for those with good prognosis and high expression, 
##no points for those with good prognosis and low expression, no points for those with poor prognosis and high expression, 
##and 1 point for those with poor prognosis and low expression. Select a significant maker;
rm(list = ls())
library(ggpubr)
library(survminer)
library(survival)
setwd("F:/CRC-IHC/combined marker score/")
#traverse folder;
list.files = list.files(pattern="completed_maker*")

name = c("Normal tissue£ºGland","Normal tissue£ºStroma","Paracancerous tissue£ºGland","Paracancerous tissue£ºStroma","Invasive margin£ºStroma","Invasive margin£ºTumor"," Tumor center£ºStroma","Tumor center£ºTumor")

#1:OS,2:RFS;
method = 'OS';

#Read the cutoff value of xtile;
if (method == 'OS'){
  xtile = read.csv("OS xtile cutoff.csv")
}else{
  xtile = read.csv("RFS xtile cutoff.csv")
}
xtile = xtile[,-1]

#Select the log rank test P value to read based on the current mode;
if (method == 'OS'){
  lrtP = read.csv("OS p.csv")
}else{
  lrtP = read.csv("RFS p.csv")
}

lrtP = ifelse(lrtP < 0.05, '1', '0')

#Grouping based on the cutoff of xtile
for (num in 1:length(list.files)){
  
  data = read.csv(list.files[num])
  m = nrow(data);
  n = ncol(data);
  
  for (i in 5:n){
    threshold = xtile[i-4, num]
    data[,i] = ifelse(data[,i] > threshold, '2', '1');
  }
  
  #Select significant makers;
  select = cbind();
  for (z in 1:15){
    if (lrtP[z, num] == 1)
      select = c(select, z)
  }
  data = data[, c(1:4,select+4)]
  m = nrow(data);
  n = ncol(data);
  
  #Univariate Cox analysis;
  hr = matrix(nrow=n-4 ,ncol=1);
  if (method == 'OS'){
    mySurv = Surv(data[,1], data[,2]);
  }else{
    mySurv = Surv(data[,3], data[,4]);
  }
  
  for (i in 5:n){
    group = data[,i];
    res.cox = coxph(mySurv ~ group);
    cox = summary(res.cox);
    coef = cox$conf.int[1];
    
    hr[i-4] = coef;
  }
  
  #Score each maker based on HR and expression;
  for (i in 1:m){
    for (j in 5:n){
      if (hr[j-4] > 1){ #Markers with poor prognosis;
        if (data[i,j] == 2){
          data[i,j] = 0;}#Poor prognosis, high expression, 0 score;
        else{
          data[i,j] = 1;}#Poor prognosis, low expression, 1 score;
      }
      else{ # Markers with good prognosis;
        if (data[i,j] == 2){
          data[i,j] = 1;}#Good prognosis, high expression, 1 score;
        else{
          data[i,j] = 0;}#Good prognosis, low expression, 0 score;
      }
    }
  }
  
  #Sum up the score data and only retain the accumulated scores;
  for (i in 1:m){
  data[i,5] = sum(data[i, c(5:n)]==1)
  }
  data = data[,-c(6:n)];

  #Take the median and calculate the p-value;
  group = as.numeric(data[,5])
  mv = median(group);
  group = ifelse(group > mv, '2', '1');
  if(length(unique(group)) == 1){
    group = as.numeric(data[,5])
    mv = median(group) - 1;
    group = ifelse(group > mv, '2', '1');
  }
  log1 = survdiff(mySurv ~ group);
  p1 = pchisq(log1$chisq, 1, lower.tail=FALSE);
  
  # plot KM curve
  if (method == 'OS'){
    fit = survfit(mySurv ~ group)
    n1 = sum(group==1)
    leg1 = paste("Low score(", n1, ")", sep = "")
    n2 = sum(group==2)
    leg2 = paste("High score(", n2, ")", sep = "")
    
    png(filename = paste("combined score KM Selected/OS/OS Combined score_", num, ".png", seq=""), width = 5.5, height = 5.5,
        units = "cm", res = 300, pointsize = 7)
    plot(fit, mark.time=TRUE, xlab = "Months", ylab = "Survival rate", lty = 1:2,
         col = 1:2, cex = 0.5)
    grid()
    legend(x = 87,y = 0.35, legend = c(leg1, leg2), lty = 1:2,
           col = 1:2, cex = 0.65)
    text(10, 0.1, paste("p=", formatC(p1, format="g", digits = 3), sep = ""),
         pos = 4, cex = 1)
    title(paste(name[num]),cex.main = 1)
    dev.off()
  }else{
    fit = survfit(mySurv ~ group)
    n1 = sum(group==1)
    leg1 = paste("Low score(", n1, ")", sep = "")
    n2 = sum(group==2)
    leg2 = paste("High score(", n2, ")", sep = "")
    
    png(filename = paste("combined score KM Selected/RFS/RFS Combined score_", num, ".png", seq=""), width = 5.5, height = 5.5,
        units = "cm", res = 300, pointsize = 7)
    plot(fit, mark.time=TRUE, xlab = "Months", ylab = "Relapse rate", lty = 1:2,
         col = 1:2, cex = 0.5)
    grid()
    legend(x = 84,y = 0.35, legend = c(leg1, leg2), lty = 1:2,
           col = 1:2, cex = 0.65)
    text(10, 0.1, paste("p=", formatC(p1, format="g", digits = 3), sep = ""),
         pos = 4, cex = 1)
    title(paste(name[num]),cex.main = 1)
    dev.off()
  }
}
