getwd()
setwd("~/Dropbox/Rworkspace/projects/phlpp2snp")
rm(list = ls())
library(GGtools)
library(snpStats)
#generate SnpMatrix set
load("gbmbldsnp6.rda")
gbmbldsnp6[1:5,1:5]
gbmbldsnp6 = gbmbldsnp6+1
gbmbldsnp6 = as.data.frame(gbmbldsnp6)
for (i in 1:ncol(gbmbldsnp6))gbmbldsnp6[,i] = as.raw(gbmbldsnp6[,i])
gbmsnp6 = new("SnpMatrix",as.matrix(t(gbmbldsnp6)))
snpsum = summary(gbmsnp6)
str(gbmsnp6)
save(gbmsnp6,file = "gbmsnp6.rda")
# intersect rppa and snp data
load("gbmrppaset.rda")
load("gbmsnp6.rda")
ls()
xy = intersect(rownames(gbmsnp6),colnames(gbmRppaSet))
length(xy)
ex = gbmRppaSet[,xy]
sml = gbmsnp6[xy,]
# shrink SNP numbers by MAF
x0 = col.summary(sml)
View(x0)
nrow(x0[which(x0$MAF>0.05),])
#result 697446
as(sml[1:5,1:5], "matrix")
str(x0)
sml = sml[,col.summary(sml)$MAF>0.05]

str(sml1)
snp = as(sml, "matrix")
exp = exprs(ex)
pda = pData(ex)
rm(sex)
head(exp)
pakt = exp[grep("Akt_p",rownames(exp)),]
head(pakt)
all(colnames(pakt) == rownames(snp))
aktsnp = data.frame(gender = pda$gender,t(pakt),snp)
aktsnp[which(aktsnp$gender == "male"),"gender"] = "1"
aktsnp[which(aktsnp$gender == "female"),"gender"] = "2"
head(aktsnp)
colnames(aktsnp)[2:3] = c("pakt473","pakt308")
aktsnp = as.data.frame(t(aktsnp))
for (i in 1:ncol(aktsnp)) aktsnp[,i] = as.numeric(as.character(aktsnp[,i]))
aktsnp = as.data.frame(t(aktsnp))
akt473snpresult = data.frame(snp = colnames(aktsnp)[4:ncol(aktsnp)], pvalue = 1, diff = 0,fdr = 1)
dim(akt473snpresult)
for (i in 1:nrow(akt473snpresult)){
	snpn = akt473snpresult[i,1]
	x = summary(lm(aktsnp$pakt473~aktsnp[,snpn]+aktsnp$gender))
	akt473snpresult[i,"pvalue"] = x$coefficients[2,4]
	akt473snpresult[i,"diff"] = x$coefficients[2,1]
	if (i %% 5000 == 0)writeLines(paste("Calculating ", i,",",round(i/nrow(akt473snpresult)*100, digits = 2), "% done.", sep = ""))
}
save(akt473snpresult,file = "akt473snplinear.rda")
akt473snpresult = akt473snpresult[order(akt473snpresult$pvalue),]
View(akt473snpresult)
load("akt473snplinear.rda")






clin = read.csv(file = "../glioma clin/tcgagbmclinfinal.csv")

load("~/Dropbox/Rworkspace/projects/glioma clin/gbm133expSet.rda")



View(clin)
dim(clin)
rownames(clin)= clin[,1]
xy = intersect(rownames(clin),colnames(rppa))
length(xy)
rppa = as.matrix(rppa[,xy])
clin = clin[xy,]