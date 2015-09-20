hapannot = read.csv(file = "~/Dropbox/Rworkspace/TCGA_Illumina_HumanHap550K.csv",head = T, sep = ",", as.is = T,stringsAsFactors = F, check.names = F)
hapannot$chromosome = sapply(strsplit(hapannot$Aliases,split = ":"),"[[",1)
hapannot$location = sapply(strsplit(hapannot$Aliases,split = ":"),"[[",2)
hapannot$ABgenotype = sapply(strsplit(hapannot$genotype,split = "="),"[[",2)
hapannot$ABgenotype = gsub("\\[", "",hapannot$ABgenotype)
hapannot$ABgenotype = gsub("\\]", "",hapannot$ABgenotype)
hapannot = hapannot[,-3]
hapannot = hapannot[,-6]
hapannot[which(hapannot$chromosome == "ChrM"),"chromosome"] = "chrM"
rownames(hapannot) = hapannot$Reporter
save(hapannot,file = "hap550annotTCGA.rda")
load("hap550annotTCGA.rda")
library(rtracklayer)
hap550gr <- with(hapannot, GRanges(chromosome, IRanges(as.numeric(as.character(location)), as.numeric(as.character(location)),names = Reporter), strand = "+", gene, Reporter, Sequence = hapannot$"Reporter Sequence"))
chain <- import.chain("~/annotation_data/hg18ToHg19.over.chain")
hap550grhg19 <- liftOver(hap550gr, chain)
hapgrhg19 = unlist(hap550grhg19)
haphg19 <- data.frame(names=mcols(hapgrhg19)$Reporter, chromosome=seqnames(hapgrhg19), location=start(hapgrhg19), strands=strand(hapgrhg19), sequence = mcols(hapgrhg19)$Sequence,stringsAsFactors = F)
save(haphg19,file = "hap550annothg19.rda")
#use ncbi2r to get hg38 annotation from NCBI website
load("~/tcgadata/hap550annot.rda")
library(NCBI2R)
xinfo = GetSNPInfo(rownames(hapannot)[grep("rs",rownames(hapannot))])
load("xinfo.rda")
hapinfo = xinfo[,c("marker","chr","chrpos","current.rsid")]
View(hapinfo)
x = hap550annot[,1:3]
x = unique(hap550annot[,1:3])
table(hapinfo$chr)

load("gbmallhap.rda")
hapsnp = t(as(gbmhap,"numeric"))
xdy = setdiff(rownames(hapsnp),hapinfo$marker)
xdyinfo = hapannot[xdy,] 
xdyinfo = xdyinfo[,c("Reporter","chromosome","location")]
colnames(xdyinfo) = c("marker","chr","chrpos")
xdyinfo$current.rsid = xdyinfo$marker
hap550gr <- with(xdyinfo, GRanges(chr, IRanges(as.numeric(as.character(chrpos)), as.numeric(as.character(chrpos)),names = marker), strand = "+", marker))
chain <- import.chain("~/annotation_data/hg18ToHg19.over.chain")
hap550grhg19 <- liftOver(hap550gr, chain)
hapgrhg19 = unlist(hap550grhg19)
haphg19 <- data.frame(marker=mcols(hapgrhg19)$marker, chr=seqnames(hapgrhg19),chrpos=start(hapgrhg19), current.rsid=mcols(hapgrhg19)$marker, stringsAsFactors = F)

hapinfo = rbind(xdyinfo,hapinfo)
rownames(hapinfo) = hapinfo$marker
dim(hapinfo)

hapinfo$ABgenotype = hapannot[match(rownames(hapinfo),rownames(hapannot)),"ABgenotype"]

xy = intersect(rownames(hapsnp),hapinfo$marker)
hapinfo = hapinfo[xy,]
hapannot = hapannot[xy,]
table(hapinfo$chr)
hapinfo[which(hapinfo$chr == "chrM"),"chr"] = "mt"
hapinfo1 = hapinfo[which(hapinfo$chr != ""),]
hapinfo1$chr = toupper(hapinfo1$chr)
hapinfo1$chr = paste("chr",hapinfo1$chr,sep = "")
hapinfo1$chr = sub("chrMT","chrM",hapinfo1$chr)

all(hapinfo1$ABgenotype == hapannot$ABgenotype)
hapinfo1$ABgenotype = chartr("ATCG","TAGC",hapinfo1$ABgenotype)

View(head(hapinfo1))
chain <- import.chain("~/annotation_data/hg38ToHg19.over.chain")
hap550gr <- with(hapinfo1, GRanges(chr, IRanges(as.numeric(as.character(chrpos)), as.numeric(as.character(chrpos)),names = marker), strand = "+", marker,ABgenotype))
hap550grhg19 <- liftOver(hap550gr, chain)
hapgrhg19 = unlist(hap550grhg19)
haphg19 <- data.frame(marker=mcols(hapgrhg19)$marker, chr=seqnames(hapgrhg19),chrpos=start(hapgrhg19), current.rsid=mcols(hapgrhg19)$marker, ABgenotype = mcols(hapgrhg19)$ABgenotype, stringsAsFactors = F)

haphg19$Agenotype = sapply(strsplit(haphg19$ABgenotype,split = "/"),"[[",1)
haphg19$Bgenotype = sapply(strsplit(haphg19$ABgenotype,split = "/"),"[[",2)
xy = intersect(rownames(hapsnp),haphg19$marker)
hapgbm19 = data.frame(haphg19[xy,c("marker","chr","chrpos","ABgenotype","Agenotype","Bgenotype")],hapsnp[xy,])
colnames(hapgbm19) = gsub("\\.","-",colnames(hapgbm19))
hapgbm[which(hapgbm$chr == "chrM"),"chr"] = "mt"
table(hapgbm19$chr)
hapgbm = hapgbm[which(hapgbm$chr != ""),]
dim(hapgbm)
View(head(h))
h = hapgbm
h$chr = paste("chr",h$chr)
for(i in 8:ncol(h)){
	h[which(h[,i] == 0),i] = paste(h[which(h[,i] == 0),"Agenotype"],h[which(h[,i] == 0),"Agenotype"],sep = "")
	h[which(h[,i] == 1),i] = paste(h[which(h[,i] == 1),"Agenotype"],h[which(h[,i] == 1),"Bgenotype"],sep = "")
	h[which(h[,i] == 2),i] = paste(h[which(h[,i] == 2),"Bgenotype"],h[which(h[,i] == 2),"Bgenotype"],sep = "")
	if(i%%10 == 0) writeLines(paste("Calculating ", i,",",round(i/ncol(h)*100, digits = 2), "% done.", sep = ""))
}
i=7

library(GenABEL)

sapply(4:ncol(h),function(i)h[is.na(h[,i]),i]== "--")
h[4,8] == NA
h= h1
h1[is.na(h1)] = "--"

colnames(h)[1:3] = c("Name","Chr","Pos")
clingbm = read.csv(file = "~/tcgadata/Clinical and Molecular Subclass Data Table.csv")
xny = intersect(colnames(h)[4:ncol(h)],clingbm$Case.ID)
length(xny)
hapmeta = clingbm[which(clingbm[,1] %in% xny),c(1,4,3)]
dim(hapmeta)
colnames(hapmeta) =  c("id","sex","age")
hapmeta[which(hapmeta$sex == "MALE"),"sex"] = "1"
hapmeta[which(hapmeta$sex == "FEMALE"),"sex"] = "0"
rownames(hapmeta) = hapmeta$id
hapmeta = hapmeta[xny,]

#h = hapsig
#rownames(h) = NULL
#h = h[,-2]
#h[is.na(h[,3]),3]
#colnames(h)[1:3] = c("Name","Chr","Pos")
#for (i in 4:ncol(h))h[is.na(h[,i]),i] = "00"
#xy = intersect(colnames(h[,4:ncol(h)]),colnames(gbmRppaSet))
#length(xy)
#ex = gbmRppaSet[,xy]
#
#exp = exprs(ex)
#pda = pData(ex)
#pakt = exprs(ex)[grep("Akt_p",rownames(exp)),]
#pakt = as.data.frame(t(pakt))
#colnames(pakt) = c("pakt473","pakt308")
#paktp = data.frame(pda[,2:5],pakt)
#paktp[which(paktp$gender == "male"),"gender"] = "1"
#paktp[which(paktp$gender == "female"),"gender"] = "0"
#for( i in 4:6)paktp[,i] = as.numeric(as.character(paktp[,i]))
#paktp = data.frame(rownames(paktp),paktp)
#colnames(paktp)[1:3] = c("id","sex","age")
#hapmeta$id = paste(paste('\"',hapmeta$id,sep = ''),'\"',sep = "")
#hapmeta$id = gsub('\\\\',"",hapmeta$id)
#hapmeta$id = gsub('\"',"",hapmeta$id)
write.table(hapmeta, file = "happheno.txt",row.names = F,quote = F)

View(h[1:10,])
h = h[,c(colnames(h)[1:3],xny)]
rownames(h) = NULL
#colnames(h)[4:ncol(h)] = paste(paste('\\"',colnames(h)[4:ncol(h)],sep = ''),'\\"',sep = "")
#colnames(h)[4:ncol(h)] = gsub('\\\\',"",colnames(h)[4:ncol(h)])
#colnames(h)[4:ncol(h)] = gsub('\"',"",colnames(h)[4:ncol(h)])
load("haphg38gbm.rda")
all(colnames(h)[4:ncol(h)] == hapmeta$id)
for (i in names(table(h$Chr))){
	x = h[which(h$Chr == i),]
	name = paste("hapgbmhg38chr",i,sep = "")
	filename = paste(name,".txt",sep = "")
	raw = paste(name,".raw",sep = "")
	ped = paste(paste("merlin",name,sep = ""),".ped",sep = "")
	dat = paste(paste("merlin",name,sep = ""),".dat",sep = "")
	map = paste(paste("merlin",name,sep = ""),".map",sep = "")
	write.table(x, file = filename,row.names = F,quote = F)
	library(GenABEL)
    convert.snp.illumina(inf = filename,out = raw,strand = "+")
    df <- load.gwaa.data(phe="happheno.txt", gen=raw, force=TRUE)
    export.merlin(df,pedfile = ped, datafile = dat, mapfile = map)
}

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

#shell commands
mach1 -d merlinhapgbmhg38chr1.dat -p merlinhapgbmhg38chr1.ped -h ~/snp/CEU+YRI/hm3_r2_b36_fwd.CEU+YRI.chr1.hap -s hm3_r2_b36_fwd.CEU+YRI.chr1.snps --rounds 50 --greedy --geno
mach1 -d merlinhapgbmhg38chr7.dat -p merlinhapgbmhg38chr7.ped --rounds 50 --states 200 --phase --interim 5 --sample 5  --prefix sample7.pp | tee machchr7.log

for ()
for i in $(ls *.bam); do (samtools view -bh -t hg19.fa.fai $m 'chr16:71677000-71759000' > phlpp${m%.sorted*}.bam); done
for i in $(seq 20 22); do (mach1 -d merlinhapgbmhg38chr$i.dat -p merlinhapgbmhg38chr$i.ped --rounds 50 --states 200 --phase --interim 5 --sample 5  --prefix merlinhapgbmhg38chr$i.pp | tee machchr$i.log); done


x = read.table(file = "hudsonalpha.org_GBM.HumanHap550.10.3.0.sdrf.txt",head = T, sep = "\t", as.is = T,stringsAsFactors = F, check.names = F)