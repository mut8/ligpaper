require("fields")

metaprot.red<-metaprot[which(is.na(metaprot[6])==F),T]
metaprot.red[6:ncol(metaprot.red)]<-metaprot.red[6:ncol(metaprot.red)]/rowSums(metaprot.red[6:ncol(metaprot.red)])

cond<-which(is.na(metaprot[6])==F)


rda.metaprot<-rda(metaprot.red[,6:ncol(metaprot)])

pch.met<-c(rep(21,4), rep(22,4), rep(23,4),rep(24,4))
col.met<-rep(colscale, 4)
labels.metaprot<-c("Doth","Euro","Leot","Sacc","Sord","Agar","Trem", "Usti", "Ther","Bact", "Acti", "Cyan", "Firm", "Fuso","Verr", "Dict", "Alph", "Beta", "Gamm", "Delt", "Epsi")

constrains1<-data.frame(alldata$C.N_lit[cond], alldata$C.P_lit[cond], alldata$days[cond])
colnames(constrains)<-c("C:N lit", "C:P lit", "Incubation time")

constrains1a<-data.frame(alldata$C.N_lit[cond], alldata$C.P_lit[cond])
colnames(constrains)<-c("C:N lit", "C:P lit")

constrains2<-data.frame(alldata$N_lit[cond], alldata$P_lit[cond], log(alldata$days[cond]))
colnames(constrains2)<-c("N lit", "P lit", "log(Incubation time)")

constrains2a<-data.frame(alldata$N_lit[cond], alldata$P_lit[cond])
colnames(constrains2)<-c("N lit", "P lit")



rda.metaprot<-rda(metaprot.red[,6:ncol(metaprot)])
ord1<-cca(metaprot.red[,6:ncol(metaprot)], constrains1)
ord1a<-cca(metaprot.red[,6:ncol(metaprot)], constrains1a)
ord2<-cca(metaprot.red[,6:ncol(metaprot)], constrains2)
ord2a<-cca(metaprot.red[,6:ncol(metaprot)], constrains2a)

cor.test(alldata$N_lit[cond], alldata$P_lit[cond])

plot.meta(ord1a, enz=F, arrows=T)
plot.meta(ord1a, enz=F, choices=c(1,3), arrows=T)
plot.meta(ord2a, enz=F, arrows=T)
plot.meta(ord2a, enz=F, choices=c(1,3), arrows=T)


dataframe

plot(ord, type="n")
points(scores(ord2, display="sites", choices=1:2), pch=pch.met, bg=col.met, cex=1)
text(scores(ord2, display="species", choices=1:2), labels=labels.metaprot, cex=0.5)

text(scores(ord2, display="bp", 
            labels=
              colnames(constrains2)
            ))

pdf("metaprot_ords2.pdf")


plot.meta(ord, arrows=T, ylim=c(-3,4), enz=F)
plot.meta(ord2, arrows=T, enz=F)
plot.meta(ord, arrows=T, envfit=F, enz=F)
plot.meta(ord2, arrows=T, envfit=F, enz=F)
dev.off()

arrows(rep(0,3), rep(0,3), 
           scores(ord2, display="bp", choices=1),  
           scores(ord2, display="bp", choices=2), length=0.1, angle=10)
?arrows
plot.cca
,
         scaling = 2, type, xlim, ylim, const, ...)
plot(ord2, display=c("sp", "lc"))         
, scale=F)


pdf("metraprot_ords3.pdf")
plot.meta(cca(metaprot.red[,6:ncol(metaprot)], scale=T), enz=F, lab="CCA")
plot.meta(cca(metaprot.red[,6:ncol(metaprot)], scale=T), enz=F, lab="CCA", choices=c(1,3))
dev.off()

ord<-cca(metaprot.red[,6:ncol(metaprot)])

plot.meta<-function(ord, arrows=F, envfit=T, enz=F, spe.mult=1, site.mult=1, lab="PCA", choices=1:2, ...) {

xvar<-eigenvals(ord)/sum(eigenvals(ord))
plot(ord, choices=choices, type="n", tck=.01,
 xlab=paste(lab, choices[1], formatC(xvar[choices[1]]*100, digits=3), "% variance"),
 ylab=paste(lab, choices[2], formatC(xvar[choices[2]]*100, digits=3), "% variance"))  

points(scores(ord, display="sites", choices=choices)*site.mult, pch=pch.met, bg=col.met, cex=1)
text(scores(ord, display="species", choices=choices)*spe.mult, labels=labels.metaprot, cex=0.5)

if (arrows==T) {
arrows(rep(0,length(scores(ord2, display="bp", choices=choices))), rep(0,length(scores(ord, display="bp"))), 
           scores(ord, display="bp", choices=choices[1])*2,  
           scores(ord, display="bp", choices=choices[2])*2, length=0.1, angle=10)
text(scores(ord, display="bp", choices=choices)*2.4, labels=rownames(scores(ord, display="bp", choices=choices)), cex=0.5)
}

lims<-par("usr")
xspa<-(lims[2]-lims[1])/50
yspa<-(lims[4]-lims[3])/50

xmin<-tapply(scores(ord, display="sites", choices=choices[1]), metaprot.red$Harvest, min)
xmax<-tapply(scores(ord, display="sites", choices=choices[1]), metaprot.red$Harvest, max)
ymin<-tapply(scores(ord, display="sites", choices=choices[2]), metaprot.red$Harvest, min)
ymax<-tapply(scores(ord, display="sites", choices=choices[2]), metaprot.red$Harvest, max)

for(i in 1:4)
  rect(xmin[i]*site.mult-xspa, ymin[i]*site.mult-yspa, xmax[i]*site.mult+xspa,ymax[i]*site.mult+yspa)

if (envfit==T) {
fit.var<-data.frame(alldata$days[cond])
colnames(fit.var)<-"Incubation time"
plot(envfit(ord, fit.var, choices=choices), cex=.5, p.max=0.05)

fit.var2<-alldata[cond,38:43]
colnames(fit.var2)<-c("C lit", "N lit", "P lit", "C:N lit", "C:P lit", "N:P lit")
plot(envfit(ord, fit.var2, choices=choices), cex=.5, p.max=0.05)


fit.var3<-alldata[cond, 50:55]
colnames(fit.var3)<-c("C mic", "N mic", "P mic", "C:N mic", "C:P mic", "N:P mic")
plot(envfit(ord, fit.var3, choices=choices), cex=.5, p.max=0.05)


fit.var4<-alldata[cond, 75:77]
colnames(fit.var4)<-c("C:N imb", "C:P imb", "N:P imb")
plot(envfit(ord, fit.var4, choices=choices), cex=.5, p.max=0.05)


fit.var5<-data.frame(metaprot.red[,4])
colnames(fit.var5)<-"F:B"
plot(envfit(ord, fit.var5, choices=choices), cex=.5, p.max=0.05)
}

if (enz==T) {
fit.enz<-data.frame(alldata$phen2cell[cond], alldata[cond, 11:16]/alldata$C_mic[cond])
colnames(fit.enz)<-c("Phen:Cell", "Cell", "Chit","Phos", "Prot", "Pero", "Phen")
plot(envfit(ord, fit.enz, choices=choices), cex=.5, p.max=0.05)


}
}

