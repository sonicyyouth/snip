filepath = "/Volumes/My Passport/gbmsnpgenotype"
load("GBM_SNP_L2_FileStructure.rda")
fmap = GbmSnpL2File
default = getwd()
setwd(filepath)
x0 = read.table(file = "STAIR_p_TCGA_Batch7_Affx_N_GenomeWideSNP_6_F12_238496.birdseed.data.txt", header = TRUE, sep = "\t", quote = "", skip = 1, stringsAsFactors = FALSE, check.names = FALSE)
x1 = x0
x0 = x0[order(x0[,1]),]
snp = x0[,1]
snp1 = x1[,1]
m = matrix(c(rep(snp,each =2),rep(c("Call","Confidence"),length(snp))),nrow = 2*length(snp),ncol =2,byrow = F)
m1 = matrix(c(rep(snp1,each =2),rep(c("Call","Confidence"),length(snp1))),nrow = 2*length(snp1),ncol =2,byrow = F)
snpdata = as.data.frame(m)
snpdata1 = as.data.frame(m1)
for ( i in 1:nrow(fmap)){
  fname = fmap[i,2]
  sname = fmap[i,1]
  x = read.table(file = fname, header = TRUE, sep = "\t", quote = "", skip = 1, stringsAsFactors = FALSE, check.names = FALSE)
  x1 = matrix(t(x[,2:3]),ncol = 1,byrow = T)
  if(all(x[,1] == snp1)){
    snpdata1 = cbind(snpdata1, sample = x1)
    colnames(snpdata1)[ncol(snpdata1)] = sname
    writeLines(paste("Calculating ", i,",",round(i/nrow(fmap)*100, digits = 2), "% done.", sep = ""))
  }
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

gbmsnp6 = snpdata1
save(gbmsnp,file = "gbmsnp6.rda")
rm(list = ls())
filepath = "/Volumes/My Passport/gbmsnpgenotype"
setwd(filepath)
load("gbmsnp.rda")
library(pd.genomewidesnp.6)
con <- pd.genomewidesnp.6@getdb()
dbListTables(con)
anno <- dbGetQuery(con, "select * from featureSet")
dim(anno)
snpanno = anno
save(snpanno,file = "snp6annotation.rda")
phlpp2 = anno[grepl("PHLPP2",anno[,"gene_assoc"],ignore.case=T),]

setwd("/Volumes/My Passport/gbmhap550")
file = dir()
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
save(gbmhapsnp, file = "gbmhapsnp.rda")
grep("rs10500560",x[,1])
tar = x[c(1,469793),]
y1 = as.character(tar[2,4:ncol(tar)])
table(y1)
aa = colnames(tar[,which(tar[2,] == "AA")])
bb = colnames(tar[,which(tar[1,] == "BB")])
ab = colnames(tar[,which(tar[1,] == "AB")])

aa = substr(aa,1,15)
ab = substr(ab,1,15)
bb = substr(bb,1,15)

aas = substr(aa,1,12)
abs = substr(ab,1,12)
bbs = substr(bb,1,12)

aaless = setdiff(aas[!duplicated(aas)],aas[duplicated(aas)])
abless = setdiff(abs[!duplicated(abs)],abs[duplicated(abs)])
bbless = setdiff(bbs[!duplicated(bbs)],bbs[duplicated(bbs)])

aaname = setdiff(aas[!duplicated(aas)],intersect(aaless,abless))
bbname = setdiff(bbs[!duplicated(bbs)],intersect(bbless,abless))
abname = setdiff(abs[!duplicated(abs)],intersect(aaless,abless))
abname = setdiff(abname,intersect(bbless,abless))

aa = data.frame(sample = aaname, genotype = "AA")
ab = data.frame(sample = abname, genotype = "AB")
bb = data.frame(sample = bbname, genotype = "BB")
snpgenotype = rbind(aa,ab,bb)
rownames(snpgenotype) = snpgenotype[,1]
snpgenotype = snpgenotype[order(snpgenotype[,1]),]
save(snpgenotype, file = "ressnpsample.rda")

setwd("/Users/liuqun/cloud360/R workspace")
load("~/cloud360/R workspace/GBM_microarray_expression.rda")
load("~/cloud360/R workspace/projects/phlpp2snp/phlppmir.rda")

x = phlppmirgbm[,which(substr(colnames(phlppmirgbm),14,15) == "01")]
colnames(x) = substr(colnames(x),1,12)

xy = intersect(colnames(x),snpgenotype[,1])
y = t(snpgenotype[xy,])
x0 = as.data.frame(x[,xy])
all(colnames(x0) == colnames(y))
phlpp = data.frame(t(y),t(x0[535:536,]),t(x0[1:534,]),stringsAsFactors = F)
save(phlpp,file = "phlpp_SNP_MIRNA.rda")
phlpp = phlpp[,-4]
boxplot(as.numeric(phlpp[,3])~factor(phlpp[,2]))

x = phlpp[which(phlpp[,2] != "AB"),]
boxplot(as.numeric(x[,3])~factor(x[,2]))
t.test(as.numeric(x[,3])~factor(x[,2]))

x = phlpp[which(phlpp[,2] != "AA"),]
boxplot(as.numeric(x[,3])~factor(x[,2]))
t.test(as.numeric(x[,3])~factor(x[,2]))

x = phlpp[which(phlpp[,2] != "BB"),]
boxplot(as.numeric(x[,3])~factor(x[,2]))
t.test(as.numeric(x[,3])~factor(x[,2]))

aa = phlpp[which(phlpp[,2] == "AA"),]
ab = phlpp[which(phlpp[,2] == "AB"),]
bb = phlpp[which(phlpp[,2] == "BB"),]

data = bb

x = as.numeric(data[,3])
data = rbind(cor = rep(0,ncol(data)),corpvalue = rep(0,ncol(data)),data)
data = as.data.frame(data)
for (i in 3:ncol(data)) data[,i] = as.numeric(as.character(data[,i]))

for (i in 4:ncol(data)){
  y = as.numeric(data[3:nrow(data),i])
  c = cor.test(x,y,method = "spearman")
  data[2,i] = c$p.value
  data[1,i] = c$estimate
}
mircor = data[1:2,]
mircor = mircor[,-1]
mircor[1:2,2] = mean(x)

aamircor = mircor
aamircor[1:2,1] = "AA"

bbmircor = mircor
bbmircor[1:2,1] = "BB"

abmircor = mircor
abmircor[1:2,1] = "AB"

if(all(colnames(aamircor) ==colnames(bbmircor))&all(colnames(aamircor) ==colnames(abmircor))) snpmir = data.frame(t(aamircor),t(abmircor),t(bbmircor),stringsAsFactors = F)
snpmir[1,1] = "AA_cor"
snpmir[1,2] = "AA_cor_pvalue"
snpmir[1,3] = "AB_cor"
snpmir[1,4] = "AB_cor_pvalue"
snpmir[1,5] = "BB_cor"
snpmir[1,6] = "BB_cor_pvalue"

colnames(snpmir) = snpmir[1,]
snpmir = snpmir[-1,]
for (i in 1:ncol(snpmir)){ snpmir[,i] <- as.numeric(as.character(snpmir[,i]))}
snpmir = snpmir[order(snpmir[,1]),]
write.csv(snpmir,file = "gbmsnpmir.csv")

load("snp6annot.rda")
x = snp6anno[,4]
library(NCBI2R)
for(i in 1:40){
  n = i*10000
  snpinfo = GetSNPInfo(x[(490001+n):(500000+n)])
  var = paste("snpinfo",i,sep = "")
  assign(var,snpinfo)
  file = paste(var,".rda",sep = "")
  save()
  f
  
}


