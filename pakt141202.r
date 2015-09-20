setwd("~/Dropbox/Rworkspace/projects/phlpp2snp")
load("pakthapresult.rda")
load("gbmallhap.rda")
load("gbmrppaset.rda")

View(head(gbmhap))
View(pakthapresult)
pakthapresult$padj = p.adjust(pakthapresult$pvalue, method = "BH")
pakthapsig = rownames(pakthapresult[pakthapresult$pvalue<0.01,])
pakthapsig[1:10]
hapsig = t(as(gbmhap[,pakthapsig],"numeric"))
#library(NCBI2R)
#snplist<-c("rs12345","rs333","rs624662")
#GetSNPInfo("rs1234567")
#library(SNPlocs.Hsapiens.dbSNP.20120608)
#rsid2loc("rs1234567")
#myrsids <- c("rs2639606", "rs75264089", "rs73396229", "rs55871206", "rs10932221", "rs56219727", "rs73709730", "rs55838886", "rs3734153", "rs79381275", "rs1516535")
#x = rsid2loc(myrsids)
#load("haphg19info.rda")
#sigsnphap = haploc[match(pakthapsig,haploc[,1]),]
#library(NCBI2R)
#snplist<-c("rs12345","rs333","rs624662")
#x = GetSNPInfo(pakthapsig)
#load("hap550snpinfo.rda")
#sigsnphap = hapsnpinfo[match(pakthapsig,hapsnpinfo[,1]),c(1,4,5)]
#all(sigsnphap[,1] == rownames(hapsig))
#hapsig = data.frame(sigsnphap,hapsig)
#colnames(hapsig) = gsub("\\.","-",colnames(hapsig))
#gr <- rsidsToGRanges(pakthapsig)
library("biomaRt")
snpmart = useMart("snp", dataset="hsapiens_snp")
filters = listFilters(snpmart)
attributes = listAttributes(snpmart)
View(attributes)
View(filters)
x = getBM(c('refsnp_id','allele','chr_name','chrom_start'), filters = "snp_filter", values = pakthapsig, mart = snpmart)
dim(x)
dim(hapsig)
length(unique(x[,1]))
library(NCBI2R)
rownames(hapsig)[1:10]
GetSNPInfo(setdiff(pakthapsig,unique(x[,1])))
rownames(hapsig)[grep("rs12020398",rownames(hapsig))] = GetSNPInfo(setdiff(pakthapsig,unique(x[,1])))$"current.rsid"
x = getBM(c('refsnp_id','allele','chr_name','chrom_start'), filters = "snp_filter", values = rownames(hapsig), mart = snpmart)
x1 = x[!duplicated(x[,1]),]
rownames(x1) = x1$refsnp_id
all(x1[rownames(hapsig),1] == rownames(hapsig))
hapsig = data.frame(x1[rownames(hapsig),],hapsig)
colnames(hapsig) = gsub("\\.","-",colnames(hapsig))
for(i in 1:ncol(hapsig))hapsig[,i] = as.character(hapsig[,i])
hapsig$chr_name = paste("chr",hapsig$chr_name,sep = "")
hapsig[which(hapsig$allele ==""),1:5]
GetSNPInfo(hapsig[which(hapsig$allele ==""),1])
grep(hapsig[which(hapsig$allele ==""),1][1],hapsig[,1])
grep(hapsig[which(hapsig$allele ==""),1][2],hapsig[,1])
hapsig[1114,2:4] = c("A/C","chr7","14280477")
hapsig= hapsig[-4347,]

View(hapsig)
h = hapsig
for(i in 1:nrow(h)){
	g = strsplit(as.character(h[i,2]),"/")[[1]]
	h[i,which(h[i,] == "0")] = paste(g[1],g[1],sep = "")
	h[i,which(h[i,] == "1")] = paste(g[1],g[2],sep = "")
    h[i,which(h[i,] == "2")] = paste(g[2],g[2],sep = "")
}
colnames(h) = gsub("\\.","-",colnames(h))

hapsig = h
hapsig$chr_name = sub("chr","",hapsig$chr_name)

write.csv(hapsig,file = "pakthapsiggeno.csv")
h = hapsig
rownames(h) = NULL
h = h[,-2]
h[is.na(h[,3]),3]
colnames(h)[1:3] = c("Name","Chr","Pos")
for (i in 4:ncol(h))h[is.na(h[,i]),i] = "00"
xy = intersect(colnames(h[,4:ncol(h)]),colnames(gbmRppaSet))
length(xy)
ex = gbmRppaSet[,xy]

exp = exprs(ex)
pda = pData(ex)
pakt = exprs(ex)[grep("Akt_p",rownames(exp)),]
pakt = as.data.frame(t(pakt))
colnames(pakt) = c("pakt473","pakt308")
paktp = data.frame(pda[,2:5],pakt)
paktp[which(paktp$gender == "male"),"gender"] = "1"
paktp[which(paktp$gender == "female"),"gender"] = "0"
for( i in 4:6)paktp[,i] = as.numeric(as.character(paktp[,i]))
paktp = data.frame(rownames(paktp),paktp)
colnames(paktp)[1:3] = c("id","sex","age")

write.table(paktp, file = "pakrhapsiggenopheno.dat",row.names = F,quote = F)


h = h[,c(colnames(h)[1:3],xy)]

write.table(h, file = "pakrhapsiggeno.txt",row.names = F,quote = F)
hapsig = read.csv(file = "pakthapsiggeno.csv")
View(hapsig)
tes = hapsig[1,]
tes[6] == "2"
GetSNPInfo("rs6947359")
tes[2]
library(GenABEL)
convert.snp.illumina(inf = "pakrhapsiggeno.txt",out = "genos.raw")
df <- load.gwaa.data(phe="pakrhapsiggenopheno.dat", gen="genos.raw", force=TRUE)
str(df)
export.merlin(df)
View(sigsnphap)
dim(hapsig)
xy = intersect(rownames(hapsig),colnames(gbmRppaSet))
#xy = intersect(rownames(hapsig),colnames(gbm133a))
length(xy)
ex = gbmRppaSet[,xy]
sml = hapsig[xy,]
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
dim(paktsnp)
colnames(paktsnp)[6]
detach("package:glmnet",unload=T)
library(glmnet)
View(paktp)
y = paktsnp[,1]
cvfit = cv.glmnet(x,y)

class(paktsnp[1,7])

class(x[1,1])




