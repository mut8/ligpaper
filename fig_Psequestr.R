<<figphos, fig=T, echo=F, results=hide>>

phos <-  read.csv("raw_data/phos.csv")

phosplot.data<-data.frame(matrix(ncol=8, nrow=11))
  
colnames(phosplot.data)<-c("type","days", "P_lit_mean", "P_lit_se", "P_mic_mean", "P_mic_se", "P_mic_TDP_mean", "P_mic_TDP_se")

phosplot.data$type <- c("AK", "AK", "","KL", "KL","", "OS", "OS","", "SW", "SW")
phosplot.data$days <-   c(rep(c(14, 181, 99),3),14, 181)

for (i in 1:nrow(phosplot.data)) {
  if(is.element(i, c(3,6,9))) {phosplot.data[i,3:8]<-0} else {
  phosplot.data$P_lit_mean[i] <- mean(phos$P_lit[phos$type==phosplot.data$type[i] & phos$days==phosplot.data$days[i]])
  phosplot.data$P_lit_se[i] <- stderr(phos$P_lit[phos$type==phosplot.data$type[i] & phos$days==phosplot.data$days[i]])
  phosplot.data$P_mic_mean[i] <- mean(phos$P_mic[phos$type==phosplot.data$type[i] & phos$days==phosplot.data$days[i]])
  phosplot.data$P_mic_se[i] <- stderr(phos$P_mic[phos$type==phosplot.data$type[i] & phos$days==phosplot.data$days[i]])
  phosplot.data$P_mic_TDP_mean[i] <- mean(phos$P_mic_TDP[phos$type==phosplot.data$type[i] & phos$days==phosplot.data$days[i]])
  phosplot.data$P_mic_TDP_se[i] <- stderr(phos$P_mic_TDP[phos$type==phosplot.data$type[i] & phos$days==phosplot.data$days[i]])
  }
}

ylim<-c(0, max(phosplot.data$P_lit_mean+phosplot.data$P_lit_se)*1.2)
height<-phosplot.data$P_lit_mean
main<-"main"
ylab<-"µg P g -1 d.w."
xlab<-""
error<-  phosplot.data$P_lit_se

barplot2(height, ylim=c(ylim[1], ylim[2]), plot.ci=TRUE, ci.u=height+error,ci.l=height, #names.arg=names.arg,
           main=main, ylab=ylab, col=col, tck=0.01, xlab=xlab)

height<-phosplot.data$P_mic_TDP_mean
main<-""
ylab<-""
xlab<-""
error<-phosplot.data$P_mic_TDP_se
width2=0.5
space2=1.4

barplot2(height, ylim=c(ylim[1], ylim[2]), plot.ci=TRUE, ci.u=height+error,ci.l=height, #names.arg=names.arg,
         main=main, ylab=ylab, col="grey", tck=0.01, xlab=xlab, add=T, width=width2, space=space2)

height<-phosplot.data$P_mic_mean
main<-""
ylab<-""
xlab<-""
error<-phosplot.data$P_mic_se


barplot2(height, ylim=c(ylim[1], ylim[2]), plot.ci=TRUE, ci.u=height+error,ci.l=height, #names.arg=names.arg,
         main=main, ylab=ylab, col="black", tck=0.01, xlab=xlab, add=T, width=width2, space=space2)



@