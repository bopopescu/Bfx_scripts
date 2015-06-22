setwd("C:/MyDoc/Projects/HIV/GenoRules/PBM/DRV")
tmp <- read.delim("DRV_finalCoef.txt",as.is=T)
library(RColorBrewer)

mcol <- colorRampPalette(brewer.pal(10,"Spectral"))(256)

isIAS <- tmp[,"isIAS"]=="y"

png("ColorCoef_DRV.png", width=200,height=800)
par(mar=c(2,2,2,8))
y <- 1:nrow(tmp)
x=1
op=par(mar=c(2,2,2,10))
image(x,y,as.matrix(-t(tmp[,2])),col=mcol,axes=F, main="DRV/r")
abline(h=c(1:nrow(tmp))+0.5)
box()
axis(4,at=1:nrow(tmp),labels=tmp[,1],las=2)
mtext(rep("*",length(isIAS)),side=4,at=c(1:nrow(tmp))[isIAS],line=4,col=2,cex=2)
mtext(tmp[,4],side=4,at=1:nrow(tmp),outer=F,line=6,las=2)
mtext("Frequency",side=3,adj=0.95,outer=T,line=-1.5,cex=1.2)
op
dev.off()



setwd("C:/MyDoc/Projects/HIV/GenoRules/PBM/ATV")
tmp <- read.delim("ATV.txt",as.is=T)[-48,]
mcol <- colorRampPalette(brewer.pal(10,"Spectral"))(256)

isIAS <- tmp[,"isIAS"]=="y"

png("ColorCoef_ATV.png", width=200,height=800)
par(mar=c(2,2,2,8))
y <- 1:nrow(tmp)
x=1
op=par(mar=c(2,2,2,10))
image(x,y,as.matrix(-t(tmp[,"Weight"])),col=mcol,axes=F, main="ATV")
abline(h=c(1:nrow(tmp))+0.5)
box()
axis(4,at=1:nrow(tmp),labels=tmp[,1],las=2)
mtext(rep("*",length(isIAS)),side=4,at=c(1:nrow(tmp))[isIAS],line=4,col=2,cex=2)
mtext(round(tmp[,"Frequency"],2),side=4,at=1:nrow(tmp),outer=F,line=6,las=2)
mtext("Frequency",side=3,adj=0.95,outer=T,line=-1.5,cex=1.2)
op
dev.off()



