rm(list = ls())
gc(reset = T)
load("gbmblsnp6.rda")
rownames(gbmsnp6)=gbmsnp6[,1]
gbmsnp6 = gbmsnp6[,-(1:2)]
gbmsnp6 = as.matrix(gbmsnp6)
save(gbmsnp6, file = "gbmsmpsnp6.rda")
load("snp6ch19info.rda")
snpgene = snp6annor[which(snp6annor[,7]!="intergenic"),]
table(snpgene[,7])
length(unique(snpgene[,1]))
x = unique(snpgene[,1])
y = gbmsnp6[x,]
xy = snp6anno[,c(2,4)]
xy1 = xy[match(x,xy[,2]),]
y1 = xy1[,1]
genesnp= gbmsnp6[y1,]
save(genesnp,file = "genesnp.rda")
gene = unique(snpgene[,11])
library(org.Hs.eg.db)
Symbol2id <- as.list(org.Hs.egSYMBOL2EG)
id2Symbol <- rep( names(Symbol2id),sapply(Symbol2id,length))
names(id2Symbol) <- unlist(Symbol2id)
#x <- unique( with(out, c(levels(GENEID), levels(PRECEDEID), levels(FOLLOWID))) )
table(gene %in% names(id2Symbol) ) # good, all found
snpgene$GENESYMBOL <- id2Symbol[ as.character(snpgene$GENEID)]
n = snpgene[is.na(snpgene[,8]),]
annox = snp6anno[match(n[,1],snp6anno[,4]),]
n$enst = substr(annox[,11],1,15)
x = select(org.Hs.eg.db,n$enst,"SYMBOL","ENSEMBLTRANS")
x2 = x[!duplicated(x[,1]),]
for (i in 1:nrow(n)) if (n[i,12] %in% x2[,1])n[i,11] = x2[which(x2[,1] == n[i,12]),2]
for (i in 1:nrow(n))snpgene[which(snpgene[,1] == n[i,1]),11] = n[i,11]
snpgene$affid = snp6anno[match(snpgene[,1],snp6anno[,4]),2]
gbm133 = read.table("GBM.transcriptome__ht_hg_u133a__broad_mit_edu__Level_3__gene_rma__data.data.txt",head = T, sep = "\t",stringsAsFactor = F,check.names = F)
gbm133 = gbm133[-1,]
rownames(gbm133) = gbm133[,1]
data(hgnc.table)
# #require function:CheckGeneSmybol
gbm133[,1] = CheckGeneSymbol(gbm133[,1])
# gbm133[gbm133[i,1] %in% hgnc.table[,1],1] =  hgnc.table[match(gbm133[gbm133[i,1] %in% hgnc.table[,1],1],hgnc.table[,1]),2]
# }

# for (i in 1:nrow(gbm133)){
#   if (gbm133[i,1] %in% hgnc.table[,1]){
#     gbm133[i,1] = hgnc.table[match(gbm133[i,1],hgnc.table[,1]),2]
#   }
# }
x15 = gbm133[duplicated(gbm133[,1]),1]
mean(as.numeric(gbm133[grep(x15[1],gbm133[,1])[1],2:539]))
mean(as.numeric(gbm133[grep(x15[1],gbm133[,1])[2],2:539]))
gbm133 = gbm133[-(grep(x15[4],gbm133[,1])[2]),]
gbm133 = gbm133[-(grep(x15[3],gbm133[,1])[1]),]
gbm133 = gbm133[-(grep(x15[2],gbm133[,1])[1]),]
gbm133 = gbm133[-(grep(x15[1],gbm133[,1])[1]),]

gbm133a = gbm133[,c(1,which(substr(colnames(gbm133),14,15) == "01"))]
colnames(gbm133a) = substr(colnames(gbm133a),1,12)
rownames(gbm133a) = gbm133a[,1]

x = colnames(gbm133a)[2:ncol(gbm133a)]
y = colnames(genesnp)
rownames(gbm133a) = gbm133a[,1]
gbm133a = gbm133a[,-1]
xy = intersect(x,y)
datax = gbm133a[,xy]
datay = genesnp[,xy]

snpsum = data.frame(rownames(datay),aa= 0,ab = 0,bb =0)
for(i in 1:nrow(snpsum)){
  snpsum[i,] = table(factor(datay[i,],levels = c("0","1","2")))
  if (i %% 500 == 0)writeLines(paste("Calculating ", i,",",round(i/nrow(snpsum)*100, digits = 2), "% done.", sep = ""))
}
su = snpsum[which(snpsum[,1]>6&snpsum[,3]>6),]
sin = snpgene[which(snpgene[,11] %in% rownames(datax)),]
sin1 = sin[which(sin[,12] %in% rownames(su)),]
su1 = su[which(rownames(su) %in% sin1[,12]),]
datay1 = datay[rownames(su1),]
all(colnames(datax) == colnames(datay1))  
tmx = rbind(datax,datay1)

snpresult = data.frame(sin1[,12],sin1[,1],sin1[,11],aa = 0,bb = 0, diff = 99, p.value = 99,mean = 99,stddiff = 99)
sult = snpresult
snpresult = snpresult[!is.na(snpresult[,3]),]
snpresult = snpresult[(snpresult[,3] %in% tmx[,1]),]
# library(stringr)
# colnames(datay) = str_replace_all(colnames(datay),"\\.","-")

# tmx[,1] = as.character(tmx[,1]

for (i in 1:nrow(snpresult)){
    n = grep(snpresult[i,1],rownames(tmx))
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
save(snpresult,file = "snp6result.rda")
   