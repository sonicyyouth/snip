#load("gbmhapsnp.rda")
#load("gbmallhap.rda")
#load("gbmsnp6.rda")
load("GBM_SNP_L2_FileStructure.rda")
filepath = "/media/liu/My Book/gbmsnpgenotype"
fmap = GbmSnpL2File
default = getwd()
setwd(filepath)
getwd()
x0 = read.table(file = fmap[1,2], header = TRUE, sep = "\t", quote = "", skip = 1, stringsAsFactors = FALSE, check.names = FALSE)
#x1 = x0
#x0 = x0[order(x0[,1]),]
#snp = x0[,1]
#snp1 = x1[,1]
colnames(x0) = c("snpprobe",fmap[1,1],"cof")
rownames(x0) =x0[,1]

#x0[which(x0$cof >0.05),2] = NA
x01 = x0[,1:2]
x02 = x0[,c(1,3)]
colnames(x01) = c("snpprobe",fmap[1,1])
colnames(x02) = c("snpprobe",fmap[1,1])
#m = matrix(c(rep(snp,each =2),rep(c("Call","Confidence"),length(snp))),nrow = 2*length(snp),ncol =2,byrow = F)
#m1 = matrix(c(rep(snp1,each =2),rep(c("Call","Confidence"),length(snp1))),nrow = 2*length(snp1),ncol =2,byrow = F)
#snpdata = as.data.frame(m)
#snpdata1 = as.data.frame(m1)
unfile = NA
for ( i in 2:nrow(fmap)){
  fname = fmap[i,2]
  sname = fmap[i,1]
  x = read.table(file = fname, header = TRUE, sep = "\t", quote = "", skip = 1, stringsAsFactors = FALSE, check.names = FALSE)
  #x1 = matrix(t(x[,2:3]),ncol = 1,byrow = T)
  #x[which(x$Confidence >0.1),2] = NA
  if(all(x[,1] == rownames(x0))){
    x01 = cbind(x01, sample = x[,2])
    x02 = cbind(x02, sample = x[,3])
    colnames(x01)[ncol(x01)] = sname
    colnames(x02)[ncol(x02)] = sname
    writeLines(paste("Calculating ", i,",",round(i/nrow(fmap)*100, digits = 2), "% done.", sep = ""))
  }else unfile = c(unfile,fname)
  #x = x[order(x[,1]),]
  #snpname = x[,1]
  #x1 = matrix(t(x[,2:3]),ncol = 1,byrow = T)
  #if (all(snpname == snp)){
  #  snpdata = cbind(snpdata, sample = x1)
  #  colnames(snpdata)[ncol(snpdata)] = sname
  #  writeLines(paste("Calculating ", i,",",round(i/nrow(fmap)*100, digits = 2), "% done.", sep = ""))
  #}
}


####for (i in 1:nrow(x)){
##  name = x[i,1]
##  ca = x[i,2]
##  co = x[i,3]
##  nu = grep(name,snpdata[,1])
##  snpdata[nu[1],ncol(snpdata)] = ca
##  snpdata[nu[2],ncol(snpdata)] = co
##}

gbmsnp6 = x01
gbmsnp6cof = x02
colnames(gbmsnp6cof) = substr(colnames(gbmsnp6cof),1,15)
colnames(gbmsnp6) = substr(colnames(gbmsnp6),1,15)
save(gbmsnp6,file = "gbmsnp6.rda")
save(gbmsnp6cof,file = "gbmsnp6cof.rda")

load("/media/liu/My Book/gbmsnpgenotype/gbmsnp6.rda")
load("/media/liu/My Book/gbmsnpgenotype/gbmsnp6cof.rda")
all(gbmsnp6[,1] == snp6anno[,2])
gbmsnp6 = data.frame(gbmsnp6[,1],snp6anno[,c(4,5,6,9,10)],gbmsnp6[3:ncol(gbmsnp6)],check.names = F)
gbmsnp6$rsid = snp6anno[match(gbmsnp6[,1],snp6anno[,"man_fsetid"]),"dbsnp_rs_id"]
all(gbmsnp6[,2] == snp6anno[,4])
gbmsnp6 = gbmsnp6[!duplicated(gbmsnp6[,2]),]
rownames(gbmsnp6) = gbmsnp6$dbsnp_rs_id

dim(gbmsnp6)
View(gbmsnp6[which(gbmsnp6$rsid == "rs12635398"),])

load("~/tcgadata/snp6annot.rda")

#txgasnp6 = read.table(file = "~/tcgadata/TCGA_Genome_Wide_SNP_6.ADF", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE, check.names = FALSE)

ncol(gbmsnp6)
gbmsnp6[1,7]
save(gbmsnp6,file = "gbmsnp6.rda")
load("gbmsnp6.rda")

#library("plyr")
#library("dplyr")
#library("data.table")
gbmsnp = gbmsnp6[,7:ncol(gbmsnp6)]
#rownames(gbmsnp) = colnames(gbmsnp6)[7:ncol(gbmsnp6)]
#snpdt = data.table(testset,keep.rownames = T)
#setkey(snpdt,rn)
rm(gbmsnp6)
#snp = sapply(gbmsnp)
#snp = gbmsnp %>%rowwise() %>% do(as.numeric)

gbmsnp = as.matrix(t(gbmsnp))
gbmsnp = t(gbmsnp)
#snpdt = t(snpdt)
colnames(gbmsnp) = rownames(gbmsnp6)

load("/media/liu/My Book/gbmsnpgenotype/gbmsnp6ad.rda")
dim(snp6)
load("/media/liu/My Book/gbmhap550/gbmhapsnp.rda")
View(head(gbmhapsnp))

colnames(gbmhapsnp) = substr(colnames(gbmhapsnp),1,15)
rownames(gbmhapsnp) = gbmhapsnp[,1]
gbmhap = lapply(4:ncol(gbmhapsnp),function(i)gbmhapsnp[which(gbmhapsnp[,i] =="AA"),i] = "0")
gbmhap = lapply(4:ncol(gbmhap),function(i)gbmhap[which(gbmhap[,i] =="AB"),i] = "1")
gbmhap = lapply(4:ncol(gbmhap),function(i)gbmhap[which(gbmhap[,i] =="BB"),i] = "2")

for(i in 4:ncol(gbmhapsnp)){
	gbmhapsnp[which(gbmhapsnp[,i] =="AA"),i] = "0"
	gbmhapsnp[which(gbmhapsnp[,i] =="AB"),i] = "1"
	gbmhapsnp[which(gbmhapsnp[,i] =="BB"),i] = "2"
	if (i %% 1000 == 0)writeLines(paste("Calculating ", i,",",round(i/nrow(gbmhapsnp)*100, digits = 2), "% done.", sep = ""))
}

for(i in 4:ncol(gbmhapsnp)){
	gbmhapsnp[,i] = as.numeric(as.character(gbmhapsnp[,i]))
}
x= gbmhapsnp[,4:ncol(gbmhapsnp)]
x = t(x)
cmean = colMeans(x,na.rm = T)
cmeany = colMeans(y,na.rm = T)

condy = which(cmeany>1)
condx = which(cmean>1)
length(cond)
length(condx)

x [,condx] = abs(x[,condx] - 2)
for(i in 7:ncol(snp6)){
	snp6[,i] = as.numeric(as.character(snp6[,i]))
}
save(snp6, file = "/media/liu/My Book/gbmsnpgenotype/gbmsnp6ad.rda")
save(gbmhapsnp,file = "/media/liu/My Book/gbmhap550/gbmhapsnpad.rda")
#snp = gbmsnp[,1:10000]
#genocount = lapply(1:ncol(gbmsnp),function(x)table(gbmsnp[,x]))
#save(genocount,file = "genocount.rda")
#len = lapply(genocount,function(x)length(x))
#len = unlist(len)
#clas = lapply(genocount,function(x)x[1]<x[length(x)])
#clas = unlist(clas)
#n = grep(T,clas)
#clas[1:1000]
#table(clas)
#x = gbmsnp[,clas]
#y = gbmsnp[,n]
#gbmsnp[,n] = abs(gbmsnp[,n] - 2)
#snpdt = data.table(gbmsnp,keep.rownames = T,key = rn)
#setkey(snpdt,rn)
#
#length(len[len==2])
#for(i in 1:nrow(gbmsnp6)){
#	x = table(as.character(gbmsnp6[i,7:1170]))
#	if ((length(x) == 3 & x[1]<x[3]) | (length(x) == 2 & x[1]<x[2])){
#		gbmsnp6[i,7:1170] = gbmsnp6[i,7:1170] - 2
#	}
#    if (i %% 1000 == 0)writeLines(paste("Calculating ", i,",",round(i/nrow(gbmsnp6)*100, digits = 2), "% done.", sep = ""))
#}
#save(gbmsnp6 , file = "gbmsnp6temp.rda")
#x1 = as.matrix(gbmsnp6[1:300000,7:1170])
#x2 = as.matrix(gbmsnp6[300001:600000,7:1170])
#x3 = as.matrix(gbmsnp6[600001:906598,7:1170])
#rm(gbmsnp6)
#gc(reset=T)
#for(i in 1:nrow(x1)){
#	x = table(as.character(x1[i,]))
#	if (length(x) = 3 & x[1]<x[3]){
#		x1[i,] = x1[i,] - 2
#	}
#}
#for(i in 1:nrow(x2)){
#	x = table(as.character(x2[i,]))
#	if (length(x) = 3 & x[1]<x[3]){
#		x2[i,] = x2[i,] - 2
#	}
#}
#for(i in 1:nrow(x3)){
#	x = table(as.character(x3[i,]))
#	if (length(x) = 3 & x[1]<x[3]){
#		x3[i,] = x3[i,] - 2
#	}
#}
#x0 = rbind(x1,x2)
#x0 = rbind(x0,x3)
#save(x0,file = "gbmsnp6temp.rda")
#
#i = 11
colnames(gbmhapsnp)[1:3] = c("dbsnp_rs_id", "chrom", "physical_pos")
colnames(snp6)[1] = "probeid"
colnames(snp6) = gsub("\\.","-",colnames(snp6))

u = union(colnames(snp6),colnames(gbmhapsnp))
s = intersect(colnames(snp6),colnames(gbmhapsnp))
snp = intersect(rownames(snp6),rownames(gbmhapsnp))
gsnp6 = snp6[snp,s]
hap = gbmhapsnp[snp,s]
View(gsnp6)
View(hap)
table(as.character(gsnp6[4,4:ncol(gsnp6)]))


table(as.numeric(as.character(snp6[1,])))
load("hapgbmhg19corrwith1000g.rda")
load("gbmrppaset.rda")
load("/media/liu/My Book/gbmsnpgenotype/gbmsnp6.rda")

table(substr(colnames(gbmhapsnp),14,15))

gbmhap = gbmhapsnp[,!duplicated(substr(colnames(gbmhapsnp),1,12))]
haploc = gbmhapsnp[,1:3]
gbmhap = gbmhap[,-(1:3)]

u = union(rownames(gbmsnp6),rownames(gbmhapsnp))
s = intersect(rownames(gbmsnp6),rownames(gbmhapsnp))
snp = intersect(colnames(gbmsnp6),colnames(gbmhapsnp))
snp6 = gbmsnp6[s,snp]
hap = gbmhap[s,snp]
snp6 = as(snp6, "numeric")
hap = as(hap,"numeric")
View(hap)
View(snp6)
length(snp)
length(u)
ls()
class(gbmRppaSet)
show(gbmRppaSet)
summary(gbmRppaSet)
rppa = exprs(gbmRppaSet)

x = intersect(u,colnames(rppa))