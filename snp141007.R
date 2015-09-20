getwd()
setwd("~/boxsync/R workspace/projects/phlpp2snp")
rm(list = ls())
gc(reset = T)
load("gbmblsnp6.rda")
View(gbmblsnp6)
gbmbldsnp6 = gbmblsnp6[-(seq(2,nrow(gbmblsnp6),2)),]
View(head(gbmbldsnp6))
rm(gbmblsnp6)
rownames(gbmbldsnp6)=gbmbldsnp6[,1]
gbmbldsnp6 = gbmbldsnp6[,-(1:2)]
gbmbldsnp6 = as.matrix(gbmbldsnp6)
biocLite("pd.genomewidesnp.6")
biocLite("oligo")
library(pd.genomewidesnp.6)
help(pd.genomewidesnp.6)
con <- pd.genomewidesnp.6@getdb()
dbListTables(con)
snp6anno = dbGetQuery(con, "select * from featureSet")
dim(snp6anno)
View(head(snp6anno))
snp6an = data.frame(probe= rownames(gbmbldsnp6),snpid = snp6anno[match(rownames(gbmbldsnp6),snp6anno[,2]),4])
View(head(snp6an))
sum(is.na(snp6an[,2]))
all(rownames(gbmbldsnp6)==snp6an[,1])
load(snp6annot.rda)
x = rownames(gbmbldsnp6)[10001]
snp6anaffy[grep(x,snp6anaffy[,1]),2] == snp6an[10001,2]
snp6an[10001,]
rownames(gbmbldsnp6) = snp6an[,2]
colnames(gbmbldsnp6) = substr(colnames(gbmbldsnp6),1,12)
save(gbmbldsnp6, file = "gbmbldsnp6.rda")
load("gbmbldsnp6.rda")
load("snp6ch19info.rda")
load("F:/annotation_data/snp6annot.rda")
#ls()
rm(x)
#View(head(snp6annot))
#View(head(snp6anaffy))
#nrow(snp6anaffy)
#nrow(snp6annot)
#nrow(gbmbldsnp6)
ls()
rm(snp6annot)
rm(snp6anaffy)
rm(snp6antcga)
rm(snp6an)
snpgene = snp6annor[which(snp6annor[,7]!="intergenic"),]
table(snpgene[,7])
length(unique(snpgene[,1]))
#x = unique(snpgene[,1])
genebldsnp = gbmbldsnp6[x,]
colnames(genebldsnp) = substr(colnames(genebldsnp),1,12)
save(genebldsnp,file = "genebldsnp.rda")
load("genebldsnp.rda")
gene = unique(snpgene$GENEID)
library(org.Hs.eg.db)
Symbol2id <- as.list(org.Hs.egSYMBOL2EG)
id2Symbol <- rep( names(Symbol2id),sapply(Symbol2id,length))
names(id2Symbol) <- unlist(Symbol2id)
#x <- unique( with(out, c(levels(GENEID), levels(PRECEDEID), levels(FOLLOWID))) )
table(gene %in% names(id2Symbol) ) # good, all found
snpgene$GENESYMBOL <- id2Symbol[ as.character(snpgene$GENEID) ]
n = snpgene[is.na(snpgene[,8]),]
View(snpgene)
annox = snp6anno[match(n[,1],snp6anno[,4]),]
n$enst = substr(annox[,11],1,15)
x = select(org.Hs.eg.db,n$enst,"SYMBOL","ENSEMBLTRANS")
x2 = x[!duplicated(x[,1]),]
for (i in 1:nrow(n)) if (n[i,12] %in% x2[,1])n[i,11] = x2[which(x2[,1] == n[i,12]),2]
for (i in 1:nrow(n))snpgene[which(snpgene[,1] == n[i,1]),11] = n[i,11]
snpgene$affid = snp6anno[match(snpgene[,1],snp6anno[,4]),2]
gbm133a = read.table("GBM.transcriptome__ht_hg_u133a__broad_mit_edu__Level_3__gene_rma__data.data.txt",head = T, sep = "\t",stringsAsFactor = F,check.names = F)
gbm133a = gbm133a[-1,]
rownames(gbm133a) = gbm133a[,1]

source("D:/boxsync/R workspace/Module_B.r")

#gbm133a[,1] = checkgenesymbol(gbm133a[,1])
data(hgnc.table)

adf = read.table(file = "d:/boxsync/R workspace/HT_HG-U133A.tcga.adf",head = T, row.names = NULL,sep = "\t", as.is = T,stringsAsFactors = F, check.names = F)
View(adf)
dim(adf)
x = adf[,9]
x = x[!duplicated(x)]
x1 = x[-(grep(";",x))]

length(x1)
dim(gbm133a)
xy = intersect(gbm1[,1],snpgene$GENESYMBOL)
length(xy)
length(x1[x1 %in% snpgene$GENESYMBOL])
View(snpgene)
# #require function:CheckGeneSmybol
#gbm133[,1] = CheckGeneSymbol(gbm133[,1])
#gbm133[gbm133[,1] %in% hgnc.table[,1],1] =  hgnc.table[match(gbm133[gbm133[i,1] %in% hgnc.table[,1],1],hgnc.table[,1]),2]
# }
con = hgnc.table[1:2,]
for (i in 1:nrow(gbm133a)){
  if (gbm133a[i,1] %in% hgnc.table[,1]){
    if(length(unique(hgnc.table[which(hgnc.table[,1]==gbm133a[i,1]),2]))>1){
      m = hgnc.table[which(hgnc.table[,1]==gbm133a[i,1]),]
      con = rbind(con,m)
    }
    gbm133a[i,1] = hgnc.table[match(gbm133a[i,1],hgnc.table[,1]),2]
  }
}
dim(con)
length(table(con[,1])==2)
length(con[!duplicated(con[,1]),1])
con[!is.element(con[,2],gbm133a[,1]),2]
x15 = gbm133a[duplicated(gbm133a[,1]),1]
xy = intersect(x15,con[,2])
length(x15)
grep(x15[1],hgnc.table[,2])
for (i in 2:ncol(gbm133a)) gbm133a[,i] = as.numeric(as.character(gbm133a[,i]))
gbm13 = gbm133a[order(rowMeans(gbm133a[,2:539],na.rm=T),decreasing = T),]
View(gbm13)
gbm1 = gbm13[!duplicated(gbm13[,1]),]
nrow(gbm1)
#mean(as.numeric(gbm133[grep(x15[1],gbm133[,1])[1],2:539]))
#mean(as.numeric(gbm133[grep(x15[1],gbm133[,1])[2],2:539]))
#gbm133 = gbm133[-(grep(x15[4],gbm133[,1])[2]),]
#gbm133 = gbm133[-(grep(x15[3],gbm133[,1])[1]),]
#gbm133 = gbm133[-(grep(x15[2],gbm133[,1])[1]),]
#gbm133 = gbm133[-(grep(x15[1],gbm133[,1])[1]),]
table(substr(colnames(gbm133a),14,15))
gbm133a = gbm1[,c(1,which(substr(colnames(gbm1),14,15) == "01"))]
colnames(gbm133a) = substr(colnames(gbm133a),1,12)
#rownames(gbm133a) = gbm133a[,1]
gbm133a[is.na(gbm133a[,1]),1] = rownames(gbm133a)[is.na(gbm133a[,1])]
save(gbm133a,file = "gbm133a.rda")
x = colnames(gbm133a)[2:ncol(gbm133a)]
ls()
y = colnames(genebldsnp)
x[1]
y[1]

gbm133a = gbm133a[,-1]
xy = intersect(x,y)
length(xy)
datax = gbm133a[,xy]
datay = genebldsnp[,xy]
dim(datay)
View(snpgene)
snpsum = data.frame(rownames(datay),aa= 0,ab = 0,bb =0)
for(i in 1:nrow(snpsum)){
  snpsum[i,2:4] = table(factor(datay[i,],levels = c("0","1","2")))
  if (i %% 1000 == 0)writeLines(paste("Calculating ", i,",",round(i/nrow(snpsum)*100, digits = 2), "% done.", sep = ""))
}

table(factor(datay[50,],levels = c("0","1","2")))
su = snpsum[which(snpsum[,1]>6&snpsum[,3]>6),]
dim(su)
View(snpsum)
sin = snpgene[which(snpgene[,11] %in% rownames(datax)),]
sin1 = sin[which(sin[,1] %in% su[,1]),]
su1 = su[which(su[,1] %in% sin1[,1]),]
datay1 = datay[sin1[,1],]
all(colnames(datax) == colnames(datay1))  
tmx = rbind(datax,datay1)
snpresult = data.frame(sin1[,12],sin1[,1],sin1[,11],aa = 0,bb = 0, diff = 99, p.value = 99,mean = 99,stddiff = 99)
sult = snpresult
all(snpresult[,3] %in% rownames(tmx))
sum(is.na(snpresult[,3]))     
# snpresult = snpresult[!is.na(snpresult[,3]),]
# snpresult = snpresult[snpresult[,3] %in% rownames(tmx),]
# library(stringr)
# colnames(datay) = str_replace_all(colnames(datay),"\\.","-")
View(sin1)
View(su)
save(snpresult,tmx, file = "snpexp.rda")
# tmx[,1] = as.character(tmx[,1]

for (i in 1:nrow(snpresult)){
    n = which(rownames(tmx)==snpresult[i,2])
    gene = tmx[which(rownames(tmx) == snpresult[i,3]),which(tmx[n,] == "0" | tmx[n,] == "2")]
    x = tmx[n,which(tmx[n,] == "0" | tmx[n,] == "2")]
    w = wilcox.test(as.numeric(as.character(gene))~factor(as.character(x)))
    t = t.test(as.numeric(as.character(gene))~factor(as.character(x)))
#   snpresult[i,4] = length(tmx[n,which(tmx[n,] == "0")])
#   snpresult[i,5] = length(tmx[n,which(tmx[n,] == "2")])
    snpresult[i,7] = w$p.value
    snpresult[i,6] = t$estimate[1]-t$estimate[2]
#   snpresult[i,8] = mean(as.numeric(tmx[which(tmx[,1] == snpresult[i,3]),2:ncol(tmx)]),na.rm = T)
#   snpresult[i,9] = snpresult[i,6]/snpresult[i,8]
    if (i %% 1000 == 0)writeLines(paste("Calculating ", i,",",round(i/nrow(snpresult)*100, digits = 2), "% done.", sep = ""))
}
load("snpexp.rda")
save(snpresult,file = "snp6bldresult.rda")
getwd()
load("snp6bldresult.rda")
View(snpresult)
dim(snpresult1)
snpresult = snpresult[order(snpresult[,7]),]
ls()
