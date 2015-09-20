load("~/cloud360/R workspace/projects/phlpp2snp/gbmsmpsnp6.rda")
snp = gbmsmpsnp[which(gbmsmpsnp[,2] == "Call"),]
gbmsnp6 = snp
save(gbmsnp6,file = "gbmsmpsnp6.rda")
cof = gbmsmpsnp[which(gbmsmpsnp[,2] == "Confidence"),]
rownames(cof) = cof[,1]
cof = cof[,-1]
cof = cof[,-1]
# for(i in 1:nrow(cof)){
#   cof[i,1] = length(which(gbmsmcof[i,]>0.05))
#   if (i %% 1000 == 0) writeLines(paste("Calculating ", i,",",round(i/nrow(cof)*100, digits = 2), "% done.", sep = ""))
# }
ncount = function(x){
  n = length(which(x>0.05))
  return(n)
}
cof = data.frame(outer = apply(cof,1,ncount),cof)
cof = t(cof)
cof1 =  data.frame(count = apply(cof,1,ncount),cof)
ex = rownames(cof1[which(cof1[,1]>300000),])
smp = rownames(cof1)

save(cof1,file = "gbmsmsnp6cof.rda")

# n = rownames(cof[which(cof[,1]<150),])
cof1 = as.matrix(cof1)
n = colnames(cof1[,which(cof1[1,]<150)])
rownames(snp) = snp[,1]
snp = snp[,-1]
snp = as.mtrix(snp[,-1])
snp1 = snp[n,]
gc(reset=TRUE)

snp11 = data.frame(maf = 0, snp1)
snp11 = as.matrix(snp11)
snp11[,1] = rowSums(snp11)
snp11[,1] = snp11[,1]/(2*611)


maf5 = rownames(snp11[which(snp11[,1]>0.05&snp11[,1]<0.95),])

snp2 = snp11[maf5,]

hwepval = function(x){
  n = table(as.factor(as.character(x)))
  if(length(n)==3){
    m = data.frame(c(n[1],n[2]/2),c(n[2]/2,n[3]))
    p = chisq.test(m,correct = T)$p.value
  }else p = 0
  return(p)
}

b = hwepval(snp2[1,2:612])
snp3 = snp2[,-1]
hwe = apply(snp3,1,hwepval)
snp4 = data.frame(hwe,snp3)
snp5 = snp4[which(snp4[,1]>0.0001),]
library(stringr)
snp5 = snp5[,-1]
colnames(snp5) =  str_replace(colnames(snp5),"\\.","-")
colnames(snp5) =  str_replace(colnames(snp5),"\\.","-")
gbmsnp6final = snp5
save(gbmsnp6final,file = "gbmsnp6final.rda")
rppa = read.table(file = "GBM.protein_exp__mda_rppa_core__mdanderson_org__Level_3__protein_normalization__data.data.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE, check.names = FALSE)
colnames(rppa) = substr(colnames(rppa),1,12)
x = colnames(rppa)
y = colnames(snp5)
xy = intersect(x,y)
rppa1 = rppa[,xy]
snp = snp5[,xy]
hwe = apply(snp,1,hwepval)
all(colnames(rppa1) == colnames(snp))
sr = rbind(rppa1,snp)
grep("Akt_pS473",rownames(sr))
snpname = data.frame(rownames(sr)[172:nrow(sr)],pvalue = 99,diff = 99,AA = 0, BB = 0, stringsAsFactors = F)
for (i in 172:nrow(sr)){
  akt = sr[11,which(sr[i,] == "0" | sr[i,] == "2")]
  x = sr[i,which(sr[i,] == "0" | sr[i,] == "2")]
  n = i-171
  snpname[n,4] = table(as.character(x))["0"]
  snpname[n,5] = table(as.character(x))["2"]
  if (length(table(factor(as.character(x)))) ==3){
    if (table(factor(as.character(x)))[1] > 4 & table(factor(as.character(x)))[2] > 4 ){
      w = wilcox.test(as.numeric(as.character(akt))~factor(as.character(x)))
      t = t.test(as.numeric(as.character(akt))~factor(as.character(x)))
      snpname[n,2] = w$p.value
      snpname[n,3] = t$estimate[1]-t$estimate[2]
    }
  }
  if (i %% 1000 == 0) writeLines(paste("Calculating ", i,",",round(i/nrow(snpname)*100, digits = 2), "% done.", sep = ""))
}
sn = snpname[order(snpname[,2]),]
save(sn,file = "gbmsnp6result.rda")
x = sn[1:4,1]
load("~/cloud360/R workspace/projects/phlpp2snp/snp6annot.rda")
pn = vector(mode = "character",length= 4)
for ( i in 1:4){
  n[i] = snp6anno[which(snp6anno[,2] == rownames(x)[i]),4]
}
snpname = data.frame(rownames(sr)[1:171],pvalue = 9,diff = 9)
for (i in 1:171){
  gene = sr[i,which(sr[226,] == "0" | sr[226,] == "2")]
  x = sr[226,which(sr[226,] == "0" | sr[226,] == "2")]
  #w = wilcox.test(as.numeric(as.character(akt))~factor(as.character(x)))
  t = t.test(as.numeric(as.character(gene))~factor(as.character(x)))
  snpname[i,2] = t$p.value
  snpname[i,3] = t$estimate[1]-t$estimate[2]
  writeLines(paste("Calculating ", i,",",round(i/nrow(snpname)*100, digits = 2), "% done.", sep = ""))
}
snpname = data.frame(rownames(g133), mean = 99,pvalue = 99,diff = 99,diffnorm = 99)
tmx = gbm133snp
n = 2
for (i in 1:nrow(snpname)){
  gene = tmx[(10+i),which(tmx[n,] == "0" | tmx[n,] == "2")]
  x = tmx[n,which(tmx[n,] == "0" | tmx[n,] == "2")]
  w = wilcox.test(as.numeric(as.character(gene))~factor(as.character(x)))
  t = t.test(as.numeric(as.character(gene))~factor(as.character(x)))
  snpname[i,3] = w$p.value
  snpname[i,4] = t$estimate[1]-t$estimate[2]
  snpname[i,2] = mean(as.numeric(tmx[(10+i),]),na.rm = T)
  snpname[i,5] = snpname[i,4]/snpname[i,2]
  if (i %% 500 == 0)writeLines(paste("Calculating ", i,",",round(i/nrow(snpname)*100, digits = 2), "% done.", sep = ""))
}
snpname = snpname[order(snpname[,3]),]
snp5 = data.frame(snpname,diffev = snpname[,4]/snpname[,2])
genex = snp5[which(snp5[,3]<0.01),1]

library(biomaRt)
library(pd.genomewidesnp.6)
con <- pd.genomewidesnp.6@getdb()
snp6anno <- dbGetQuery(con, "select * from featureSet")
data <- data.frame(rownames(nx), snp6anno[match(rownames(nx), snp6anno[,2]),])

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
snpmart = useMart("snp", dataset="hsapiens_snp")
input1 = getBM(c("refsnp_id","chr_name","chrom_start"), filters = "snp_filter",values = data[,5],mart = snpmart)
getBM(attributes=c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "band"), filters ="hgnc_symbol", values = genex, mart = ensembl)
