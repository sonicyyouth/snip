## build lgg snp file including all availible data
getwd()
setwd("/home/liuqun/Dropbox/Rworkspace/projects/snp")
setwd("~/Downloads/snp")
fil <- dir()
 filename <- read.table(file = "../broad.mit.edu_LGG.Genome_Wide_SNP_6.sdrf.txt",head = T,stringsAsFactors = F, sep = "\t")
snpall <- read.table(file = fil[1] ,head = F,skip = 2, stringsAsFactors = F,quote = "")
## snpall1 = snpall[,1:2]
## snpall2 <- snpall[,c(1,3)]
## colnames(snpall1)[2] = filename[match(fil[1],filename[,31]),2]
## colnames(snpall2)[2] = filename[match(fil[1],filename[,31]),2]
snpaall <- lapply(1:length(fil),function(x){
    rppa <- read.table(fil[x],head = F,skip = 2,quote = ""h, stringsAsFactors = F,sep = "\t")
    colnames(rppa)[2] = filename[match(fil[x],filename[,31]),2]
    return(rppa)
})
x =filename[match(fil,filename[,31]),2]
for(i in 1:length(snpaall))colnames(snpaall[[i]])[2] = x[i]
len <- sapply(snpaall,nrow)
all(len == 906600)
nam <- sapply(snpaall,function(x)all(x[,1] == snpaall[[1]][,1]))
length(nam) == sum(nam)
head(snpaall[[1]])
snpa = do.call(cbind,snpaall)
head(snpa)
snpa1 = snpa[,c(1,seq(2,ncol(snpa),3))]
snpa2 = snpa[,c(1,seq(3,ncol(snpa),3))]
colnames(snpa2) = colnames(snpa1)
lggsnp = snpa1
lggsnpcof = snpa2
save(lggsnp,lggsnpcof,file = "lggsnpdata.rda")
load("lggsnpsimple.rda")

load("~/Documents/annotation_data/snp6annot.rda")

rownames(lggsnp) = lggsnp[,1]
lggsnp$rsid = snp6anaffy[match(lggsnp[,1], snp6anaffy[,1]),2]
lggsnp$rsid[1:10]
ncol(lggsnp)
lggsnp[1,1000]
lggsnp = lggsnp[,c(1000,1:999)]
View(lggsnp)
rownames(lggsnp) = lggsnp$rsid
table(substr(colnames(lggsnp),14,15))
save(lggsnp,file = "lggsnpsimple.rda")

load("lggsnpsimple.rda")


mkdir snp

for FIR in ./*.0;do mv ./${FIR}/*.birdseed.data.txt ./snp; done

