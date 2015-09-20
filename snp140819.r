x = select(org.Hs.eg.db, "GO:0003723", cols, keytype="GO")
load("gbm_3utrsnpexpfinal.rda")
load("~/cloud360/R workspace/gbm_hap_3utrgeno.rda")
load("~/cloud360/R workspace/GBM_microarray_expression.rda")
load("~/cloud360/R workspace/tcga/ProcessedData.Firehose/GBM_mirna_array.rda")
mir = gbm_mirna[,-1]
x = colnames(mir)
table(substr(x,14,15))
x = x[which(substr(x,14,15) == "01")]
mir = mir[,x]
colnames(mir) = substr(colnames(mir),1,12)

x = colnames(gbm133a)
table(substr(x,14,15))
x = x[which(substr(x,14,15) == "01")]
gbm133a = gbm133a[,x]
colnames(gbm133a) = substr(colnames(gbm133a),1,12)
gbm133a = gbm133a[,unique(colnames(gbm133a))]



gbm4502a = gbm4502a[,unique(colnames(gbm4502a))]
colnames(gbm4502a) = substr(colnames(gbm4502a),1,12)

geneexp = gbm133a
y = snp133result
geneexp = gbm4502a
y = snp4502result

geneexp = as.data.frame(geneexp)
for (i in 1:ncol(geneexp)){ geneexp[,i] <- as.numeric(as.character(geneexp[,i]))}
mir = as.data.frame(mir)
for (i in 1:ncol(mir)){ mir[,i] <- as.numeric(as.character(mir[,i]))}
genx = colnames(geneexp)
mirn = colnames(mir)
geno = gbmhap3utrdatacur[,5:447]
snpy = colnames(geno)
xy = intersect(genx,snpy)
xy1 = intersect(xy,mirn)
snpgene = data.frame(t(geno[,xy]),t(geneexp[,xy]),check.names = F,stringsAsFactors = F)
snpgenemir = data.frame(t(geno[,xy1]),t(geneexp[,xy1]),t(mir[,xy1]),check.names = F,stringsAsFactors = F)
snpgene = snpgenemir
n1 = nrow(y)+1
num = length(which(y$padjust <0.05))
rs = rownames(y)[1:num]
gen = y[1:num,4]
na = paste(rs,gen,sep = "_")
x = paste(rep(na,each = 8),rep(c("AAcorr","AApval","ABcorr","ABpval","BBcorr","BBpval","AABBdiff","pval"),100),sep = "-")
result = matrix(99,nrow = nrow(geneexp)+nrow(mir),ncol = length(x))
colnames(result) = x
rownames(result) = c(rownames(geneexp),rownames(mir))
all(rownames(result) == colnames(snpgene)[n1:ncol(snpgene)])

library(HGNChelper);
xname = rownames(geneexp)
res <- checkGeneSymbols(xname)

library(psych)
options(warn = -3)

for (i in 1:num){
  rs0 = rs[i]
  gen0 = gen[i]
  n = which(colnames(snpgene) == gen0)
  if(length(n) == 0){
    gene = res[which(res[,3] == gen0),1]
    if(length(gene) == 1)
      gen0 = colnames(snpgene)[which(colnames(snpgene) == gene)]
  }
  aa = snpgene[which(snpgene[,rs0]=="AA"),]
  ab = snpgene[which(snpgene[,rs0]=="AB"),]
  bb = snpgene[which(snpgene[,rs0]=="BB"),]
  n = (i-1)*8

  aamean = mean(aa[,gen0],na.rm = T)
  abmean = mean(ab[,gen0],na.rm = T)
  bbmean = mean(bb[,gen0],na.rm = T)
  result[,(n+1)] = cor(aa[,n1:ncol(snpgene)],method = "spearman")[,gen0]
  result[,(n+2)] = r.test(result[,(n+1)],n = nrow(aa))$p
  result[,(n+3)] = cor(ab[,n1:ncol(snpgene)],method = "spearman")[,gen0]
  result[,(n+4)] = r.test(result[,(n+3)],n = nrow(ab))$p
  result[,(n+5)] = cor(bb[,n1:ncol(snpgene)],method = "spearman")[,gen0]
  result[,(n+6)] = r.test(result[,(n+5)],n = nrow(bb))$p
  if(aamean > bbmean) result[,(n+7)] = result[,(n+1)]-result[,(n+5)]
  if(aamean < bbmean) result[,(n+7)] = result[,(n+5)]-result[,(n+1)]
  result[,(n+8)] = r.test(n = nrow(aa),result[,(n+1)],result[,(n+5)],n2 = nrow(bb))$p
  for(j in 1:4) gc()

#   for(j in 1:nrow(result)){
#     tes = rownames(result)[j]
#     aas = cor.test(aa[,gen0],aa[,tes],method = "spearman")
#     abs = cor.test(ab[,gen0],ab[,tes],method = "spearman")
#     bbs = cor.test(bb[,gen0],bb[,tes],method = "spearman")
#     result[j,(n+1)] = aas$estimate
#     result[j,(n+2)] = aas$p.value
#     result[j,(n+3)] = abs$estimate
#     result[j,(n+4)] = abs$p.value
#     result[j,(n+5)] = bbs$estimate
#     result[j,(n+6)] = bbs$p.value
#     result[j,(n+7)] = result[j,(n+1)]-result[j,(n+5)]
#     result[j,(n+8)] = r.test(n = nrow(aa),result[j,(n+1)],result[j,(n+5)],n2 = nrow(bb))$p
#   }
  writeLines(paste("Calculating ", i,",",round(i/num*100, digits = 2), "% done.", sep = ""))
}
options(warn = 0)



a = seq(7,8*num-1,8)
b = seq(8,8*num,8)
ab = c(a,b)
ab = ab[order(ab)]
x = result[,ab]
x1 = x[,seq(1,2*num-1,2)]
x2 = data.frame(sum = rowSums(x1,na.rm = T),x1)

x3 = data.frame(x[,1:2], p.adjust(x[,2],"BH"),check.names = F)
colnames(x3)[3] = paste(strsplit(colnames(x)[1],split = "-")[[1]][1],"padjust",sep = "-")
for (i in 1:c(ncol(x2)-2)){
  n = i*2+1
  x3 = data.frame(x3,x[,n:(n+1)], p.adjust(x[,(n+1)],"BH"),check.names = F)
  colnames(x3)[ncol(x3)] = paste(strsplit(colnames(x)[n],split = "-")[[1]][1],"padjust",sep = "-")
}
gbmhap4502corpadj = x3
gbmhap4502diff = x2
gbmhap4502cor = result
save(gbmhap4502cor, gbmhap4502diff,gbmhap4502corpadj,file = "gbmhapsnp4502cor.rda")

gbmhap133corpadj = x3
gbmhap133diff = x2
gbmhap133cor = result
save(gbmhap133cor, gbmhap133diff,gbmhap133corpadj,file = "gbmhapsnp133cor.rda")

# result1 = result
#
# for (i in 1:num){
#   rs0 = rs[i]
#   gen0 = gen[i]
#   n = which(colnames(snpgene) == gen0)
#   if(length(n) == 0){
#     gene = res[which(res[,3] == gen0),1]
#     if(length(gene) == 1)
#       gen0 = colnames(snpgene)[which(colnames(snpgene) == gene)]
#   }
#   aa = snpgene[which(snpgene[,rs0]=="AA"),]
#   ab = snpgene[which(snpgene[,rs0]=="AB"),]
#   bb = snpgene[which(snpgene[,rs0]=="BB"),]
#   n = (i-1)*8
#
#
#   aamean = mean(aa[,gen0],na.rm = T)
#   abmean = mean(ab[,gen0],na.rm = T)
#   bbmean = mean(bb[,gen0],na.rm = T)
#   if(aamean > bbmean) result1[,(n+7)] = result1[,(n+1)]-result1[,(n+5)]
#   if(aamean < bbmean) result1[,(n+7)] = result1[,(n+5)]-result1[,(n+1)]
#   result1[,(n+8)] = r.test(n = nrow(aa),result1[,(n+1)],result1[,(n+5)],n2 = nrow(bb))$p
#   writeLines(paste("Calculating ", i,",",round(i/num*100, digits = 2), "% done.", sep = ""))
# }