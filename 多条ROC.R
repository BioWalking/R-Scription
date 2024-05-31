par(mar= c(5,5,1,1),cex.lab=1.2,cex.axis= 1.2)

#先画一个1年的图
sROC=survivalROC(Stime=overdata$month, status=overdata$OS, marker = overdata$score, predict.time =12, method="KM")
plot(sROC$FP, sROC$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="red", 
     xlab="False positive rate", ylab="True positive rate",
     lwd = 2, cex.main=1.3, cex.lab=1.5, cex.axis=1.2, font=1.2)
abline(0,1)
aucText=paste0("1 years"," (AUC=",sprintf("%.3f",sROC$AUC),")") #这个后面添加legend用

#再加一个3年的线
sROC3=survivalROC(Stime=overdata$month, status=overdata$OS, marker = 2-overdata$score, predict.time =36, method="KM")
lines(sROC3$FP, sROC3$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="green",lwd = 2)
aucText3=paste0("3 years"," (AUC=",sprintf("%.3f",sROC3$AUC),")") #这个后面添加legend用

#再加一个5年的线
sROC5=survivalROC(Stime=overdata$month, status=overdata$OS, marker = 2-overdata$score, predict.time =60, method="KM")
lines(sROC5$FP, sROC5$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="blue",lwd = 2)
aucText5=paste0("5 years"," (AUC=",sprintf("%.3f",sROC5$AUC),")") #这个后面添加legend用

#添加legend
legend("bottomright", c(aucText,aucText3,aucText5),
       lwd=2,bty="n",col=c("red","green","blue"),cex=1.2)