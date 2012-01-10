setwd("~/ligpaper")
require(vegan)
source("~/R/functions.R")
source("startup.R")


mikr<-read.csv("raw_data/metaprot_fung.csv")

mikr2<-mikr
mikr2[,22:26]<-100*mikr[,22:26]/mikr[,14]
mikr2<-mikr2[,-14]

fung<-100*mikr[c(2,6,10,14,1,5,9,13,3,7,11,15,4,8,12,16),5:12]/rowSums(mikr[c(2,6,10,14,1,5,9,13,3,7,11,15,4,8,12,16),5:12])
bak<-100*mikr[c(2,6,10,14,1,5,9,13,3,7,11,15,4,8,12,16),13:21]/rowSums(mikr[c(2,6,10,14,1,5,9,13,3,7,11,15,4,8,12,16),13:21])
prot<-100*mikr[c(2,6,10,14,1,5,9,13,3,7,11,15,4,8,12,16),22:26]/rowSums(mikr[c(2,6,10,14,1,5,9,13,3,7,11,15,4,8,12,16),22:26])

bak<-rbind(bak[1:4,T], rep(0, 9), bak[5:8,T], rep(0, 9), bak[9:12,T], rep(0, 9),bak[13:16,T])
prot<-rbind(prot[1:4,T], rep(0, 5), prot[5:8,T], rep(0, 5), prot[9:12,T], rep(0, 5), prot[13:16,T])
fung<-rbind(fung[1:4,T], rep(0, 8), fung[5:8,T], rep(0, 8), fung[9:12,T], rep(0, 8), fung[13:16,T])

par()

bak

pdf("mikr.pdf", width=6, height=10)
par(mfrow=c(3,1), mar=c( 2.1, 2.1, 3.1, 12.1), cex.axis=.5)
barplot(as.matrix(t(fung)), names.arg=c("AK1", "AK2", "AK3", "AK4","","KL1", "KL2", "KL3", "KL4","", "OS1", "OS2", "OS3","OS4","", "SW1","SW2", "SW3", "SW4"), main="fungi",  legend.text=colnames(fung), args.legend=list(x=32, y=100))
barplot(as.matrix(t(bak)), names.arg=c("AK1", "AK2", "AK3", "AK4","","KL1", "KL2", "KL3", "KL4","", "OS1", "OS2", "OS3","OS4","", "SW1","SW2", "SW3", "SW4"), main="bacteria",  legend.text=colnames(bak), args.legend=list(x=31, y=105))
barplot(as.matrix(t(prot)), names.arg=c("AK1", "AK2", "AK3", "AK4","","KL1", "KL2", "KL3", "KL4","", "OS1", "OS2", "OS3","OS4","", "SW1","SW2", "SW3", "SW4"), main="proteobacteria",  legend.text=colnames(prot), args.legend=list(x=28, y=90))
dev.off()

mikr2

pdf("mikr.pdf", width=10, height=10)
par(mfrow=c(2,1), mar=c( 2.1, 4.1, 3.1, 12.1), cex.axis=.5, tck=0.01, las=1)

bak<-100*mikr2[c(2,6,10,14,1,5,9,13,3,7,11,15,4,8,12,16),13:25]/rowSums(mikr2[c(2,6,10,14,1,5,9,13,3,7,11,15,4,8,12,16),13:25])
bak<-rbind(bak[1:4,T], rep(0, 13), bak[5:8,T], rep(0, 13), bak[9:12,T], rep(0, 13),bak[13:16,T])
barplot(as.matrix(t(fung)), names.arg=c("AK1", "AK2", "AK3", "AK4","","KL1", "KL2", "KL3", "KL4","", "OS1", "OS2", "OS3","OS4","", "SW1","SW2", "SW3", "SW4"), main="fungi",  legend.text=colnames(fung), args.legend=list(x=31, y=100), col=1:8, yaxt="n", ylab="protein abundance (% of fungi)")
axis(2)
axis(4, label=F)
barplot(as.matrix(t(bak)), names.arg=c("AK1", "AK2", "AK3", "AK4","","KL1", "KL2", "KL3", "KL4","", "OS1", "OS2", "OS3","OS4","", "SW1","SW2", "SW3", "SW4"), main="bacteria",  legend.text=colnames(bak), args.legend=list(x=31, y=100, cex=1), col=1:13, yaxt="n", ylab="protein abundance (% of bacteria)")
axis(2)
axis(4, label=F)

dev.off()

pca.bak<-rda(mikr2[,13:25], scale=F)
plot(pca.bak, type="n")
text(scores(pca.bak, display="sites")[,1:2], labels=mikr2$X, cex=.5)
text(scores(pca.bak, display="species")[,1:2], colnames(mikr2[,13:25]), cex=.5)

mds.bak<-metaMDS(mikr2[,13:25])
mds.fung<-metaMDS(mikr2[,5:12])

plot(mds.fung, type="n")
text(scores(mds.fung), labels=mikr$X, cex=.5)
text(scores(mds.fung, display="species") , colnames(mikr2[,5:12]), cex=.5)

plot(mds.bak, type="n")
text(scores(mds.bak), labels=mikr$X, cex=.5)
text(scores(mds.bak, display="species") , colnames(mikr2[,13:25]), cex=.5)



pca.fung<-rda(mikr2[,5:12,2:ncol(fung)], scale=F)
plot(pca.fung, type="n")
text(scores(pca.fung, display="sites")[,1:2], labels=mikr$X, cex=.5)
text(scores(pca.fung, display="species")[,1:2], colnames(mikr2[,5:12]), cex=.5)

H<-diversity(fung[,2:ncol(fung)])
J<-H/log(specnumber(fung[,2:ncol(fung)]))

plot(c(14,97,181,475), H[c(1,5,9,13)]  , type="o", bg=colscale[2], pch=pch[2], ylim=c(1, 2), ylab="shannon index (fungal classes)", xlab="litter incubation (days)")
lines(c(14,97,181,475), H[c(2,6,10,14)])
points(c(14,97,181,475), H[c(2,6,10,14)]  ,  bg=colscale[1], pch=pch[1])
lines(c(14,97,181,475), H[c(3,7,11,15)])
points(c(14,97,181,475), H[c(3,7,11,15)]  ,  bg=colscale[3], pch=pch[3])
lines(c(14,97,181,475), H[c(4,8,12,16)])
points(c(14,97,181,475), H[c(4,8,12,16)]  ,  bg=colscale[4], pch=pch[4])
legend("topleft", c("AK", "KL","OS", "SW"), lty=1, pch=pch, pt.bg=colscale)

plot(c(14,97,181,475), J[c(1,5,9,13)]  , type="o", bg=colscale[2], pch=pch[2], ylim=c(0.5, 1), ylab="shannon index (fungal classes)", xlab="litter incubation (days)")
lines(c(14,97,181,475), J[c(2,6,10,14)])
points(c(14,97,181,475), J[c(2,6,10,14)]  ,  bg=colscale[1], pch=pch[1])
lines(c(14,97,181,475), J[c(3,7,11,15)])
points(c(14,97,181,475), J[c(3,7,11,15)]  ,  bg=colscale[3], pch=pch[3])
lines(c(14,97,181,475), J[c(4,8,12,16)])
points(c(14,97,181,475), J[c(4,8,12,16)]  ,  bg=colscale[4], pch=pch[4])
#legend("topleft", c("AK", "KL","OS", "SW"), lty=1, pch=pch, pt.bg=colscale)


, c(rep(14,4), rep(97,4), rep(181,4), rep(475,4)), rep(c("KL", "AK", "OS","SW")))

