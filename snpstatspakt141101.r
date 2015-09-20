getwd()
setwd("~/Dropbox/Rworkspace/projects/phlpp2snp")
rm(list = ls())
load("gbmrppaset.rda")
load("gbmsnp6.rda")
load("gbm133a.rda")
gbm133a = gbm133a[,-1]
gbmsnplibrary(GGtools)
ls()
xy = intersect(rownames(gbmsnp6),colnames(gbm133a))
length(xy)
ex = gbm133a[,xy]
sml = gbmsnp6[xy,]
sml = sml[,col.summary(sml)$MAF>0.05]
sml = sml[,!duplicated(colnames(sml))]
dim(sml)
exp = exprs(ex)
pda = pData(ex)
pakt = exprs(ex)[grep("Akt_p",rownames(exp)),]
pakt = as.data.frame(t(pakt))
colnames(pakt) = c("pakt473","pakt308")
paktp = data.frame(pakt,pda[,2:5])
paktp[which(paktp$gender == "male"),"gender"] = "1"
paktp[which(paktp$gender == "female"),"gender"] = "2"
for( i in 1:4)paktp[,i] = as.numeric(as.character(paktp[,i]))
all(rownames(paktp) == rownames(sml))

x = as(sml, "numeric")
x = x+1
all(rownames(paktp) == rownames(x))
all(colnames(ex) == rownames(x))
paktsnp = data.frame(paktp,x)
with(paktsnp,boxplot(pakt473~rs11079662,main = "pakt473 level of rs4791285"))
with(paktsnp,boxplot(pakt473~rs2007530,main = "pakt473 level of rs7207499"))
genePITPsnp = data.frame(t(ex),x)
library(plyr)
mod = lapply(1:nrow(ex),function(i){lm(genePITPsnp[,i]~genePITPsnp[,"rs4791285"])})
lp = ldply(mod,function(x){summary(x)$coefficients[2,]})

rownames(lp) = rownames(ex)[1:nrow(ex)]
View(lp)
lp["rs10500560",]
lp=lp[order(lp[,4]),]
colnames(lp)=c("estimate","stderr","tvalue","pvalue")
write.csv(lp, file = "rs4791285genecor.csv")
#pakthapresult$padj = p.adjust(pakthapresult$pvalue, method = "BH")
save(pakthapresult,file = "paktsnp6linearresult.rda")
i = 10
lm(genePITPsnp[,i]~genePITPsnp[,"rs4791285"])

boxplot(genePITPsnp[,"PRKCA"]~genePITPsnp[,"rs2007530"],main = "PRKCA level of rs4791285")

dim(paktsnp)
library(plyr)
mod = lapply(7:ncol(paktsnp),function(i){lm(paktsnp$pakt473~paktsnp[,i]+paktsnp$gender)})
pakthapresult = ldply(mod,function(x){summary(x)$coefficients[2,]})
rownames(pakthapresult) = colnames(paktsnp)[7:ncol(paktsnp)]
View(pakthapresult)
pakthapresult["rs10500560",]
pakthapresult=pakthapresult[order(pakthapresult[,4]),]
colnames(pakthapresult)=c("estimate","stderr","tvalue","pvalue")
pakthapresult$padj = p.adjust(pakthapresult$pvalue, method = "BH")
xy = intersect(rownames(gbmsnp6),colnames(gbmRppaSet))
length(xy)
ex = gbmRppaSet[,xy]
sml = gbmsnp6[xy,]
sml = sml[,col.summary(sml)$MAF>0.05]
sml = sml[,!duplicated(colnames(sml))]
exp = exprs(ex)
pda = pData(ex)
pakt = exprs(ex)[grep("Akt_p",rownames(exp)),]
pakt = as.data.frame(t(pakt))
colnames(pakt) = c("pakt473","pakt308")
paktp = data.frame(pakt,pda[,2:5])
pda[which(pda$gender == "male"),"gender"] = "1"
pda[which(pda$gender == "female"),"gender"] = "2"
for( i in 1:4)pda[,i] = as.numeric(as.character(pda[,i]))
head(paktp)
ls()
rm(xy)
rm(list = ls()[1:6])
rm(sml1)
rm(pda)
gc(reset = T)
all(rownames(paktp) == rownames(sml))

load("snp6ch19info.rda")

load("gbmhapsnp.rda")
ls()
gbmhapsnp[1:10,1:10]



setwd("/media/liu/My Book/gbmhap550")
file = dir()[-1]
x = read.table(file = file[1], header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE, check.names = FALSE)
for (i in 2:length(file)){
  x0 = read.table(file = file[i], header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE, check.names = FALSE)
  if (all(x0[,1]==x[,1])){
    x0 = x0[,-c(1:3)]
    x = cbind(x,x0)
    writeLines(paste("Calculating ", i,",",round(i/length(file)*100, digits = 2), "% done.", sep = ""))
  }
}
gbmhapsnp = x
gbmhapsnp = gbmhapsnp[-1,]
rownames(gbmhapsnp) = gbmhapsnp[,1]
save(gbmhapsnp, file = "gbmhapsnp.rda")
table(substr(colnames(gbmhapsnp),14,15))

gbmhap = gbmhapsnp[,!duplicated(substr(colnames(gbmhapsnp),1,12))]
haploc = gbmhapsnp[,1:3]
gbmhap = gbmhap[,-(1:3)]
#dim(gbmbldhap)
colnames(gbmhap) = substr(colnames(gbmhap),1,12)
#gbmbldhap[1:10,1:10]
for(i in 1:ncol(gbmhap)){
	gbmhap[which(gbmhap[,i] =="AA"),i] = "1"
	gbmhap[which(gbmhap[,i] =="AB"),i] = "2"
	gbmhap[which(gbmhap[,i] =="BB"),i] = "3"
}
#gbmhap = apply(gbmhap,2,function(x){
#	x[which(x=="AA")] = "1"
#	x[which(x=="AB")] = "1"
#	x[which(x=="BB")] = "3"
#})
for(i in 1:ncol(gbmhap))gbmhap[,i] = as.numeric(as.character(gbmhap[,i]))
for(i in 1:ncol(gbmhap))gbmhap[,i] = as.raw(gbmhap[,i])
gbmhap = new("SnpMatrix",as.matrix(t(gbmhap)))
save(gbmhap,file = "gbmallhap.rda")
load("gbmallhap.rda")
snpsum = col.summary(gbmhap)
xy = intersect(rownames(gbmhap),colnames(gbmRppaSet))
xy = intersect(rownames(gbmhap),colnames(gbm133a))
length(xy)
ex = gbmRppaSet[,xy]
sml = gbmhap[xy,]
sml = sml[,col.summary(sml)$MAF>0.05]
sml = sml[,!duplicated(colnames(sml))]
exp = exprs(ex)
pda = pData(ex)
pakt = exprs(ex)[grep("Akt_p",rownames(exp)),]
pakt = as.data.frame(t(pakt))
colnames(pakt) = c("pakt473","pakt308")
paktp = data.frame(pakt,pda[,2:5])
paktp[which(paktp$gender == "male"),"gender"] = "1"
paktp[which(paktp$gender == "female"),"gender"] = "2"
for( i in 1:4)paktp[,i] = as.numeric(as.character(paktp[,i]))
all(rownames(paktp) == rownames(sml))

x = as(sml, "numeric")
all(rownames(paktp) == rownames(x))
paktsnp = data.frame(paktp,x)

ex = gbm133a[,xy]
genePITPsnphap = data.frame(t(ex),x)
lm(genePITPsnphap[,"PRKCA"]~genePITPsnphap[,"rs7207499"])

boxplot(genePITPsnphap[,"PRKCA"]~genePITPsnphap[,"rs2007530"],main = "PRKCA level of rs7207499")


library(plyr)
mod = lapply(7:10,function(i){lm(paktsnp$pakt473~paktsnp[,i]+paktsnp$gender)})
lp = ldply(mod,function(x){summary(x)$coefficients[2,]})
lp
mod = lapply(7:ncol(paktsnp),function(i){lm(paktsnp$pakt473~paktsnp[,i]+paktsnp$gender)})
pakthapresult = ldply(mod,function(x){summary(x)$coefficients[2,]})
rownames(pakthapresult) = colnames(paktsnp)[7:ncol(paktsnp)]
View(pakthapresult)
pakthapresult["rs10500560",]
pakthapresult=pakthapresult[order(pakthapresult[,4]),]
colnames(pakthapresult)=c("estimate","stderr","tvalue","pvalue")
pakthapresult$padj = p.adjust(pakthapresult$pvalue, method = "BH")
save(pakthapresult,file = "paktsnp6linearresult.rda")
load("pakthapresult.rda")
load("hap550snpinfo.rda")
View(pakthapresult)
dim(pakthapresult)
pakthapresult$pvalue[1000]
ls()
View(hapsnpinfo)

nrow(pakthapresult[pakthapresult$pvalue<0.001,])
sig = x[,rownames(pakthapresult)[1:500]]
paktsnpsig = data.frame(paktp,sig)

dim(sig)
View(sig)
cor.test(sig[,2],sig[,3])
sig[2,2]
install.packages("psych")
install.packages("GPArotation")
library(psych)
sigCor = cor(sig,use ="pairwise.complete.obs")
x = which(abs(sigCor)>0.8 & sigCor<1)%%ncol(sigCor)
y = which(abs(sigCor)>0.8&sigCor<1)%/%ncol(sigCor)+1
xy = data.frame(x,y,xy = x-y)
xy = xy[xy$xy>0,]
sigx = colnames(sig[,-c(unique(xy$x))])
x = data.frame(hapsnpinfo[match(rownames(pakthapresult)[1:500],hapsnpinfo$marker),c(1,2,4:6)],pakthapresult[1:500,])
x1 = x[which(x$marker %in% sigx),]
gc(reset = T)
View(sigCor)
fa.parallel(sig,fa = "pc", n.iter = 500, show.legend=F,main = "scree plot with parrallel analysis")
pc = principal(sig, nfactors = 10,score = T, rotate = "varimax")