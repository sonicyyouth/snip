load("~/cloud360/R workspace/gbm_hap_3utrgeno.rda")
x = gbmhap3utrdata
colnames(x) = substr(colnames(x),1,15)
colnames(x) = gsub("\\.","-",colnames(x))
nam = colnames(x)[5:849]
x01 = nam[substr(nam,14,15) == "01"]
x10 = nam[substr(nam,14,15) == "10"]

x010 = unique(substr(x01,1,12))
x100 = unique(substr(x10,1,12))
x00 = intersect(x010,x100)
x001 = paste(x00,"-10",sep = "")
nam1 = setdiff(nam,x001)
x11 = nam[substr(nam,14,15) == "11"]
x110 = unique(substr(x11,1,12))
x111 = paste(x110,"-11",sep = "")
nam2 = setdiff(nam1,x111)

x0 = data.frame(x[,1:4],x[,nam2])
colnames(x0) = substr(colnames(x0),1,12)
colnames(x0) = gsub("\\.","-",colnames(x0))
geno = x0[,5:447]

save(gbmhap3utrdata,gbmhap3utrdatacur,file = "gbm_hap_3utrgeno.rda")
load("~/cloud360/R workspace/GBM_microarray_expression.rda")
x = colnames(gbm133a)
table(substr(x,14,15))
x = x[which(substr(x,14,15) == "01")]
gbm133a = gbm133a[,x]
colnames(gbm133a) = substr(colnames(gbm133a),1,12)
gbm133a = gbm133a[,unique(colnames(gbm133a))]
colnames(gbm133a) = substr(colnames(gbm133a),1,12)

gbm4502a = gbm4502a[,unique(colnames(gbm4502a))]
colnames(gbm4502a) = substr(colnames(gbm4502a),1,12)
namex = colnames(geno)
y = colnames(gbm133a)
y = colnames(gbm4502a)
xy = intersect(namex,y)
snpgene = rbind(geno[,xy],gbm133a[,xy])
snpgene = rbind(geno[,xy],gbm4502a[,xy])
snpname = data.frame(x0[,1:4],pvalue = 99,diff = 99,AA = 0, BB = 0, mean = 0, stringsAsFactors = F)
sum(is.na(snpname[,4]))
all(rownames(snpname) == rownames(snpgene)[1:6176])

library(HGNChelper);
x133 = rownames(gbm133a)
x4502 = rownames(gbm4502a)
res <- checkGeneSymbols(x133)
res <- checkGeneSymbols(x4502)

for ( i in 1:nrow(snpname)){
  gen = snpname[i,4]
  snpname[i,7] = table(as.character(snpgene[i,]))["AA"]
  snpname[i,8] = table(as.character(snpgene[i,]))["BB"]
  n = which(rownames(snpgene) == gen)
  if(length(n) == 0){
    gene = res[which(res[,3] == gen),1]
    if(length(gene) == 1)
    n = which(rownames(snpgene) == gene)
  }
  if(length(n) >0){
    x = snpgene[,which(snpgene[i,] == "AA" | snpgene[i,] == "BB")]
    snpname[i,9] = mean(as.numeric(as.character(x[n,])))
    if (length(levels(factor(as.character(x[i,])))) == 2){
      if (table(factor(as.character(x[i,])))[1] > 6 & table(factor(as.character(x[i,])))[2] > 6 ){
        w = wilcox.test(as.numeric(as.character(x[n,]))~factor(as.character(x[i,])))
        t = t.test(as.numeric(as.character(x[n,]))~factor(as.character(x[i,])))
        snpname[i,5] = w$p.value
        snpname[i,6] = t$estimate[1]-t$estimate[2]
      }else  snpname[i,5] = 999
    }
  }
}
a22 = snpname[which(snpname[,5]==999),]
a33 = snpname[which(snpname[,5]==99),]
snpname = snpname[order(snpname[,5]),]

snp4502result = snpname

snp133result = snpname

y = snp4502result
y = snp133result
y = data.frame(y[,1:5],padjust = p.adjust(y[,5],"BH"),y[,6:9],check.names = F)
snp4502result = y
snp133result = y

save(snp133result,snp4502result,file = "gbm_3utrsnpexpfinal.rda")

x = phlpp[which(phlpp[,2] != "AB"),]
boxplot(as.numeric(x[,3])~factor(x[,2]))
t.test(as.numeric(x[,3])~factor(x[,2]))







