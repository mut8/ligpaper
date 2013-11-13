
timeseries(alldata$CN_inbal, alldata$days, alldata$Litter, pt.bg=colscale, pch=21:24, log="y")

timeseries(alldata$TDP/alldata$P_lit/10000, alldata$days, alldata$Litter, pt.bg=colscale, pch=21:24)
timeseries((alldata$TDP  +  alldata$P_mic*1000)/alldata$P_lit/10000
  , alldata$days, alldata$Litter, pt.bg=colscale, pch=21:24)

timeseries(alldata$Chitinase, alldata$days, alldata$Litter, pt.bg=colscale, pch=21:24, log="y")
timeseries(alldata$Chitinase/alldata$Phenoloxidase, alldata$days, alldata$Litter, pt.bg=colscale, pch=21:24, log="y")

timeseries(alldata$Chitinase/alldata$Cellulase, alldata$days, alldata$Litter, pt.bg=colscale, pch=21:24)

timeseries(alldata$DOC_k2SO4, alldata$days, alldata$Litter, pt.bg=colscale, pch=21:24)

timeseries(alldata$TN_K2SO4, alldata$days, alldata$Litter, pt.bg=colscale, pch=21:24)

timeseries(alldata$DOC_k2SO4/alldata$TN_K2SO4, alldata$days, alldata$Litter, pt.bg=colscale, pch=21:24)



timeseries(alldata$TN, alldata$days, alldata$Litter, pt.bg=colscale, pch=21:24)


var<-alldata$Respiration
cond <- is.na(var)==F
timeseries(var[cond],alldata$days[cond], alldata$Litter[cond])

var<-alldata$DOC
cond <- is.na(var)==F
timeseries(var[cond],alldata$days[cond], alldata$Litter[cond], pch=pch, pt.bg=colscale,endsig=T)

cor.test(alldata$DOC_k2SO4, alldata$N_lit)

var<-alldata$DOC_k2SO4
cond <-is.na(var)==F
timeseries(var[cond],alldata$days[cond], alldata$Litter[cond], pch=pch, pt.bg=colscale,endsig=T)

var<-alldata$DOC_k2SO4/alldata$TN_K2SO4
cond <-is.na(var)==F
timeseries(var[cond],alldata$days[cond], alldata$Litter[cond], pch=pch, pt.bg=colscale,endsig=T)
cor.test(alldata$DOC_k2SO4/alldata$TN_K2SO4, alldata$N_lit)

cond<-alldata$days==475
cor.test(alldata$DOC_k2SO4[cond], alldata$N_lit[cond])

 var<-alldata$TDP
cond <-is.na(var)==F
timeseries(var[cond],alldata$days[cond], alldata$Litter[cond], pch=pch, pt.bg=colscale,endsig=T)

var<-alldata$TN_K2SO4
cond <-is.na(var)==F
timeseries(var[cond],alldata$days[cond], alldata$Litter[cond], pch=pch, pt.bg=colscale,endsig=T)

cor.test(alldata$DOC_k2SO4[alldata$days==475], alldata$N_lit[alldata$days==475])
cor.test(alldata$DOC_k2SO4[alldata$days==181], alldata$N_lit[alldata$days==181])
cor.test(alldata$DOC_k2SO4[alldata$days==97], alldata$N_lit[alldata$days==97])
cor.test(alldata$DOC_k2SO4[alldata$days==14], alldata$N_lit[alldata$days==14])

cor.test(alldata$DOC_k2SO4[alldata$days==475], alldata$Respiration[alldata$days==475])
cor.test(alldata$DOC_k2SO4[alldata$days==181], alldata$Respiration[alldata$days==181])
cor.test(alldata$DOC_k2SO4[alldata$days==97], alldata$Respiration[alldata$days==97])
cor.test(alldata$DOC_k2SO4[alldata$days==14], alldata$Respiration[alldata$days==14])

cor.test(alldata$DOC_k2SO4[alldata$days==475], alldata$P_lit[alldata$days==475])
cor.test(alldata$DOC_k2SO4[alldata$days==181], alldata$P_lit[alldata$days==181])
cor.test(alldata$DOC_k2SO4[alldata$days==97], alldata$P_lit[alldata$days==97])
cor.test(alldata$DOC_k2SO4[alldata$days==14], alldata$P_lit[alldata$days==14])

cor.test(alldata$DOC_k2SO4, alldata$TN_K2SO4)

var<-alldata$DOC_k2SO4/alldata$TDP
cond <-is.na(var)==F
timeseries(var[cond],alldata$days[cond], alldata$Litter[cond], pch=pch, pt.bg=colscale,endsig=T)

var<-alldata$TN_K2SO4
cond <-is.na(var)==F
timeseries(var[cond],alldata$days[cond], alldata$Litter[cond], pch=pch, pt.bg=colscale,endsig=T)

var<-alldata$DOC_k2SO4/alldata$TN_K2SO4
cond <-is.na(var)==F
timeseries(var[cond],alldata$days[cond], alldata$Litter[cond], pch=pch, pt.bg=colscale,endsig=T)

cor.test(alldata$DOC_k2SO4/alldata$TN_K2SO4, alldata$C.N_lit)


pdf("fig_pca.pdf")
ord <- rda(rsim, scale=T)
ord.plot(ord, samples$days, samples$type, pch=pch, pt.bg=colscale, spe.label.type="text", spe.labels=peaks$origin, spe.mult=2)
ord.plot(ord, samples$days, samples$type, pch=pch, pt.bg=colscale, spe.label.type="text", spe.labels=peaks$origin, spe.mult=2, choices=c(1,3))
ord.plot(ord, samples$days, samples$type, pch=pch, pt.bg=colscale, spe.label.type="text", spe.labels=peaks$origin, spe.mult=2, choices=c(1,4))
ord.plot(ord, samples$days, samples$type, pch=pch, pt.bg=colscale, spe.label.type="text", spe.labels=peaks$origin, spe.mult=2, choices=c(3,4))
ord.plot(ord, samples$days, samples$type, pch=pch, pt.bg=colscale, spe.label.type="text", spe.labels=peaks$origin, spe.mult=2, choices=c(3,4))
plot(envfit(ord, orig_rsim, choices=c(3,4)))
plot(envfit(ord, samples[,c("days", "cons_acc_resp_litC")], choices=c(3,4)), col="green", labels=c("incubation(days)", "accumulated respiration (%C)"))
dev.off()

ord.plot(ord, samples$days, samples$type, pch=pch, pt.bg=colscale, spe.label.type="text", spe.labels=peaks$class, spe.mult=2, choices=c(3,4))
plot(envfit(ord, orig_rsim, choices=c(3,4)))

plot(envfit(ord, samples[,c("days", "cons_acc_resp_litC")], choices=c(3,4)), col="green", labels=c("incubation(days)", "accumulated respiration (%C)"))
dev.off()

var<-scores(ord, choices=3)[[2]]
cond <-is.na(var)==F
timeseries(var[cond],samples$days[cond], samples$type[cond], pch=pch, pt.bg=colscale,endsig=T)

var<-scores(ord, choices=4)[[2]]
cond <-is.na(var)==F
timeseries(var[cond],samples$days[cond], samples$type[cond], pch=pch, pt.bg=colscale,endsig=T)


corr.ab(class_rsim, scores(ord, choices=3:4)[[2]])

peaks

var
plot(
samples$days,
scores(ord, choices=1:3)[[2]]
, bg=colscale.all, pch=pch.all
     )
     


pdf("fig_mds.pdf", height=12, width=6)
ord <- metaMDS(rsim, scale=T, k=3)
par(mfrow=c(2,1))
ord.plot(ord, samples$days, samples$type, pch=pch, pt.bg=colscale, spe.label.type="text", spe.labels=peaks$class, spe.mult=1)
ord.plot(ord, samples$days, samples$type, pch=pch, pt.bg=colscale, spe.label.type="text", spe.labels=peaks$class, spe.mult=1, choices=c(1,3))
dev.off()

?metaMDS

var<-  scores(ord, choices=3)
cond <-is.na(var)==F
timeseries(var[cond],samples$days[cond], samples$type[cond], pch=pch, pt.bg=colscale,endsig=T)



corr.ab(class_rsim, scores(ord, choices=1:3))
corr.ab(orig_rsim, scores(ord, choices=1:3))

var<-alldata$Respiration/alldata$C_mic
cond <-is.na(var)==F
timeseries(var[cond],alldata$days[cond], alldata$Litter[cond], pch=pch, pt.bg=colscale,endsig=T)
samples$days[cond]

var<-alldata$Cellulase/alldata$C_mic
cond <-is.na(var)==F
timeseries(var[cond],alldata$days[cond], alldata$Litter[cond], pch=pch, pt.bg=colscale,endsig=T, log="y")

var<-alldata$Phenoloxidase/alldata$C_mic
cond <-is.na(var)==F
timeseries(var[cond],alldata$days[cond], alldata$Litter[cond], pch=pch, pt.bg=colscale,endsig=T, log="y")

samples$days[cond]

var<-alldata$Cellulase/alldata$Respiration
cond <-is.na(var)==F
timeseries(var[cond],alldata$days[cond], alldata$Litter[cond], pch=pch, pt.bg=colscale,endsig=T)

var<-alldata$P_mic+alldata$TN_K2SO4
cond <-is.na(var)==F
timeseries(var[cond],alldata$days[cond], alldata$Litter[cond], pch=pch, pt.bg=colscale,endsig=T)

var<-(alldata$P_mic+alldata$TDP)/alldata$P_lit/100
cond <-is.na(var)==F
timeseries(var[cond],alldata$days[cond], alldata$Litter[cond], pch=pch, pt.bg=colscale,endsig=T)

var<-(alldata$TDP)/alldata$P_lit/100
cond <-is.na(var)==F
timeseries(var[cond],alldata$days[cond], alldata$Litter[cond], pch=pch, pt.bg=colscale,endsig=T)

var<-(alldata$P_mic)/alldata$P_lit
cond <-is.na(var)==F
timeseries(var[cond],alldata$days[cond], alldata$Litter[cond], pch=pch, pt.bg=colscale,endsig=T)


var<-(alldata$P_mic)/alldata$TN_K2SO4
cond <-is.na(var)==F
timeseries(var[cond],alldata$days[cond], alldata$Litter[cond], pch=pch, pt.bg=colscale,endsig=T)

pcor.test(var4,var1, var2)

accresp<-read.csv("raw_data/accresp_weekly.csv")

plot(accresp$days, accresp$SW.mean, type="l", pch=24, col=colscale[4], xlim=c(0,500), ylim=c(0,.32), lwd=2)
lines(accresp$days, accresp$AK.mean, col=colscale[1], lwd=2)
lines(accresp$days, accresp$KL.mean, col=colscale[2], lwd=2)
lines(accresp$days, accresp$OS.mean, col=colscale[3], lwd=2)

lines(c(181, 475), c(mean(samples$cons_acc_resp_litC[samples$days==181& samples$type=="AK"]), mean(samples$cons_acc_resp_litC[samples$days==475& samples$type=="AK"])), lty=2, col=colscale[1], lwd=2)
lines(c(181, 475), c(mean(samples$cons_acc_resp_litC[samples$days==181& samples$type=="KL"]), mean(samples$cons_acc_resp_litC[samples$days==475& samples$type=="KL"])), lty=2, col=colscale[2], lwd=2)
lines(c(181, 475), c(mean(samples$cons_acc_resp_litC[samples$days==181& samples$type=="OS"]), mean(samples$cons_acc_resp_litC[samples$days==475& samples$type=="OS"])), lty=2, col=colscale[3], lwd=2)
lines(c(181, 475), c(mean(samples$cons_acc_resp_litC[samples$days==181& samples$type=="SW"]), mean(samples$cons_acc_resp_litC[samples$days==475& samples$type=="SW"])), lty=2, col=colscale[4], lwd=2)

plotCI(samples$days, samples$cons_acc_resp_litC, liw=samples$cons_acc_resp_litC_se, uiw=samples$cons_acc_resp_litC_se, cex=.7, add=T, gap=0, pt.bg= colscale[pch.all-20], pch=pch.all)

  accresp$days, accresp$KL.mean, uiw=accresp$LK.se, liw=accresp$SW.se, type="o", pch=22, pt.bg=colscale[2], xlim=c(0,500), ylim=c(0,.32), cex=.5, add=T)


plotCI(accresp$days, accresp$OS.mean, uiw=accresp$OS.se, liw=accresp$SW.se, type="o", pch=23, pt.bg=colscale[3], xlim=c(0,500), ylim=c(0,.32), cex=.5, add=T)


warnings()

?plotCI
       , accresp$AK.mean))