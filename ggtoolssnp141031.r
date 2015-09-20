getwd()
setwd("~/Dropbox/Rworkspace/projects/phlpp2snp")
rm(list = ls())
load("~/Dropbox/Rworkspace/projects/glioma clin/gbm133expSet.rda")
load("gbm_rppa.rda")
ls()
#generate SnpMatrix set
load("gbmbldsnp6.rda")
View(gbmbldsnp6)
class(gbmbldsnp6)
gbmbldsnp6 = as.data.frame(gbmbldsnp6)
for (i in 1:ncol(gbmbldsnp6))gbmbldsnp6[,i] = as.raw(gbmbldsnp6[,i])
gbmsnp6 = new("SnpMatrix",as.matrix(t(gbmbldsnp6)))
snpsum = summary(gbmsnp6)
str(gbmsnp6)
save(gbmsnp6,file = "gbmsnp6.rda")
load("~/annotation_data/snp6annot.rda")
xy = intersect(rownames(gbmsnp6),colnames(gbm133expSet))
length(xy)
#generate ExpressionSet for rppa data
ls()
View(rppa)
dim(rppa)
clin = read.csv(file = "../glioma clin/tcgagbmclinfinal.csv")
View(clin)
dim(clin)
rownames(clin)= clin[,1]
xy = intersect(rownames(clin),colnames(rppa))
length(xy)
rppa = as.matrix(rppa[,xy])
clin = clin[xy,]
all(rownames(clin)==colnames(rppa))
metadata <- data.frame(labelDescription=colnames(clin), row.names=colnames(clin))
adf <- new("AnnotatedDataFrame", data=clin, varMetadata=metadata)
rownames(clin)
sampleNames(adf)
gbmRppaSet = new("ExpressionSet", exprs=rppa,phenoData=adf)
save(gbmRppaSet,file = "gbmrppaset.rda")

xy = intersect(rownames(gbmsnp6),colnames(gbmRppaSet))
length(xy)
xy[1:10]
ex = gbmRppaSet[,xy]
sml = gbmsnp6[xy,]
all(colnames(ex)==rownames(sml))
gbmRppaSnpSet = make_smlSet(ex,sml)

head(gbmsnp6)