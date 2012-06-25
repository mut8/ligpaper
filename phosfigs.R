pdf("phosphorus.pdf", height=10, width=5)

phos<-read.csv("raw_data/pmicdif.csv")
head(phos)
phos$TDP

var<-alldata$Gross.P.immobilization
cond<-is.na(var)==F
timeseries(var[cond], alldata$days[cond], alldata$Litter[cond], pch=pch, pt.bg=colscale)

par(mfrow=c(3,1))
plot(alldata$P_mic/10000/alldata$P_lit, alldata$days,type="n", xlim=c(0,500),ylim=c(0,1))
timeseries(alldata$P_mic/100/alldata$P_lit, alldata$days, alldata$Litter, pch=pch, pt.bg=colscale, add=F, xlab="Litter incubation (days)", ylab="microbial P (% of total P)",ylim=c(0,100))

par(mfrow=c(3,1))
plot(phos$TDP/10000/alldata$P_lit, alldata$days,type="n", xlim=c(0,500),ylim=c(0,1))

timeseries(alldata$P_lit*10000, alldata$days, alldata$Litter, pch=pch, pt.bg=colscale, add=F, xlab="Litter incubation (days)", ylab="total P (myg g-1)")

timeseries(
  phos$TDP
  , alldata$days, alldata$Litter, pch=pch, pt.bg=colscale, add=F, xlab="Litter incubation (days)", ylab="TDP P (myg g-1)")


timeseries(phos$TDP#/100/alldata$P_lit
           , alldata$days, alldata$Litter, pch=pch, pt.bg="white", xlab="Litter incubation (days)", ylab="% of total P",ylim=c(-1,1000), topsig=F, , lwd=.5)

timeseries(alldata$P_mic/100/alldata$P_lit, alldata$days, alldata$Litter, pch=pch, pt.bg=colscale, xlab="Litter incubation (days)", ylab="microbial P (% of total P)",ylim=c(-1,100), endsig=F, lwd=2, add=T, topsig=F)

timeseries(
  (alldata$P_mic+phos$TDP)#/100/alldata$P_lit
  , alldata$days, alldata$Litter, pch=pch, pt.bg=colscale, xlab="Litter incubation (days)", ylab="microbial P (% of total P)",ylim=c(-1,150), topsig=F, add=T, lwd=2)

abline(h=100, lty=3)

var<-10*alldata$N_.mic/alldata$N_lit
cond<-is.na(var)==F
timeseries(var[cond], alldata$days[cond], alldata$Litter[cond], pch=pch, pt.bg=colscale,ylim=c(0,10))


timeseries(, alldata$days, alldata$Litter, pch=pch, pt.bg=colscale, add=T, lty=2)

colnames(alldata)

plot(alldata$P_mic/10000/alldata$P_lit, alldata$days,type="n", xlim=c(0,500),ylim=c(0,1))
timeseries((alldata$P_lit-alldata$PO4/10000-alldata$P_mic/10000)/alldata$P_lit, alldata$days, alldata$Litter, pch=pch, pt.bg=colscale, add=T, lty=2)


write.csv(  cbind(alldata[,c("Litter", "days", "P_lit", "P_mic")], phos), "phosout.csv"
            )

dev.off()
