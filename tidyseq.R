setwd("~/Dropbox/Rworkspace/projects/snp")
library(survival)
library(tidyr)
library(dplyr)
library(limma)
genesurv = function(expset,entz,gene,set = "TCGA GBM"){
    source("~/Dropbox/Rworkspace/r scripts/ggsurv.R")
    gbmgene = exprs(expset)[entz,]
    gbmgene1 = median(gbmgene)
    gbmgene[gbmgene <= gbmgene1] = 1
    gbmgene[gbmgene>gbmgene1] = 2
    gbmgene = ifelse(gbmgene == 1, "Low expression","High expression")

    sv = Surv(expset$time,as.numeric(expset$vitalstatus))
   surv = survfit(sv~gbmgene)
    logr = survdiff(sv ~ gbmgene)
    cox = coxph(sv ~ gbmgene)
    hr = hazardr(logr)
    hr = hr[order(hr)]
    pval = round(1- pchisq(logr$chisq,1),2)
    lb = paste("Log-rank P = ",pval,"\nHR = ", round(hr[2],2),"(",round(hr[1],2),"-", round(hr[3],2), ")",sep = "")
    tl = paste("Survival curve of ", gene," in ", set,sep = "")
    f0 = ggsurv(surv)+theme_classic() + labs(list(title = tl,x = "Time (days)",y = "proportion survival/surviving"))+ annotate("text", label = lb ,x = 0.72*max(expset$time), y = 0.6,size = 6) + theme(legend.position = c(.7,.8),text = element_text(size = 16),legend.title = element_blank(),legend.text = element_text(size = 20))
    return(f0)
}
load("../glioma_clin/tcgagbm133entrezmr.rda")
RNAexp <- read.table(file = "../tcgafile/LGG.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt", head = T, stringsAsFactors = F,check.names = F, sep = "\t")
View(RNAexp)

setwd("../tcgafile")
fil <- dir("./Level_3")
## fil <- paste("./Level_3/",fil, sep = "")
filename <- read.table(file = "FILE_SAMPLE_MAP.txt",head = T,stringsAsFactors = F, sep = "\t")
rppaall <- read.table(file = paste("./Level_3/",fil[1], sep = ""),head = T, stringsAsFactors = F,sep = "\t")
xd = setdiff(x[,1],x400[,1])

colnames(rppaall)[2] = filename[match(fil[1],filename[,1]),2]
ppaall <- lapply(1:length(fil),function(x){
    rppa <- read.table( paste("./Level_3/",fil[x], sep = ""),head = T, stringsAsFactors = F,sep = "\t")
    dir = "/home/liuqun/Dropbox/Rworkspace/projects/"
    colnames(rppa)[2] = filename[match(fil[x],filename[,1]),2]
    return(rppa)
})
pakt = lapply(ppaall,function(x)x[grep("Akt",x[,1]),])
pakta = do.call(cbind,pakt)
chk = sapply(seq(1,ncol(pakta),2),function(x)all(pakta[,x] == pakta[,1]))
length(chk) == sum(chk)

pakta = pakta[,c(1,seq(2,ncol(pakta),2))]

View(pakta)
save(pakta, file = "lggpakt.rda")
ls()
colnames(pakta) = substr(colnames(pakta), 1,15)
colnames(lggsnp) = substr(colnames(lggsnp), 1,15)
xy = intersect(colnames(pakta),colnames(lggsnp))
length(xy)
grep("rs10500560",lggsnp[,1])
lggsnp[grep("rs1050560",lggsnp[,1]),1]

length(ppaall)
chk = lapply(1:length(ppaall),function(x)all(ppaall[[x]][,1] == ppaall[[1]][,1]))
rppall <- do.call(rbind,ppaall)
dim(rppall)
head(ppaall[[1]])
dim(rppaall)
xd = setdiff(x[,1],x400[,1])
xd1 = setdiff(x400[,1],x[,1])
xd = setdiff(x[,1],x400[,1])
x = read.csv(file = "lggclin1.csv", sep = ",", check.names = F,stringsAsFactors = F)
x[,171] = gsub(",", ";",x[,171])
write.csv(x,file = "lggclin1.csv")
View(x)
rownames(x) = toupper(x[,6])
x = x[,-1]
x[1,11]
table(x[,6])
x[1,15:20]
x1 = x[,c(1,6,9,12,15,18,41,62,63,81,83:85,92,95:98,102,104,105,109,111:117,148,158,163:166,170:172)]
write.csv(x1,file = "lggclinsimple.csv")
View(x1)
getwd()
glmclin = read.csv(file = "gbmclin1.csv",head = T, stringsAsFactors = F, check.names = F)
View(glmclin)
glmrna = read.csv(file = "primyglmrna.csv",head = T, stringsAsFactors = F, quote = "",check.names = F)
View(glmrna)
colnames(glmrna) = substr(colnames(glmrna),1,15)
glmclin = glmclin[which(glmclin[,3] != ""),]
glmclin[,3] = toupper(glmclin[,3])
rownames(glmclin) = glmclin[,3]
write.csv(glmclin,file = "gbmclin1.csv")
rownames(glmrna) = glmrna[,1]
xy = intersect(colnames(glmrna), rownames(glmclin))
length(xy)
grep('CDC20',rownames(glmrna))
rownames(glmrna)[grep('MYC',rownames(glmrna))]
glmclin[which(glmclin[,"disease"] == "gbm"),"grade"] = "g4"
table(glmrna[3399,xy],glmclin[xy,"histology"])
myccdc = data.frame(glmclin[xy,],myc = log2(as.numeric(glmrna[11402,xy])),cdc20 = log2(as.numeric(glmrna[3399,xy])))
myccdc = myccdc[which(myccdc[,"grade"] != ""),]

## colnames(glmclin)[2] = "ages"
## colnames(glmclin)[4] = "cqcfhistology"
## colnames(glmclin)[20] = "gender"
## colnames(glmclin)[30] = "grade"
## colnames(glmclin)[22] = "histology"
View(myccdc)
library(ggplot2)
f1 = ggplot(myccdc,aes(x = grade,y = myc,colors = grade))+geom_boxplot(width = 0.7,outlier.shape = NA)  + geom_jitter(aes(color = grade),position = position_jitter(width = 0.2))+theme_classic()+theme(legend.position="none",text = element_text(size = 16),axis.text.x = element_text(size = 11))+ylab("MYC expression")
f2 = ggplot(myccdc,aes(x = histology,y = myc,colors = histology))+geom_boxplot(width = 0.7,outlier.shape = NA)  + geom_jitter(aes(color = histology),position = position_jitter(width = 0.2))+theme_classic()+theme(legend.position="none",text = element_text(size = 16),axis.text.x = element_text(size = 11))+ylab("MYC expression")
f3 = ggplot(myccdc,aes(x = grade,y = cdc20,colors = grade))+geom_boxplot(width = 0.7,outlier.shape = NA)  + geom_jitter(aes(color = grade),position = position_jitter(width = 0.2))+theme_classic()+theme(legend.position="none",text = element_text(size = 16),axis.text.x = element_text(size = 11))+ylab("CDC20 expression")
f4 = ggplot(myccdc,aes(x = histology,y = cdc20,colors = histology))+geom_boxplot(width = 0.7,outlier.shape = NA)  + geom_jitter(aes(color = histology),position = position_jitter(width = 0.2))+theme_classic()+theme(legend.position="none",text = element_text(size = 16),axis.text.x = element_text(size = 11))+ylab("CDC20 expression")
f1
f2
f3
f4
pdf("f1.pdf")
f1
dev.off()
pdf("f2.pdf")
f2
dev.off()
pdf("f3.pdf")
f3
dev.off()
pdf("f4.pdf")
f4
dev.off()
t.test(myccdc$cdc20~myccdc$grade)

hazardr = function(x){
        hr = (x$obs[1]/x$exp[1])/(x$obs[2]/x$exp[2])
        up95 = exp(log(hr) + qnorm(0.975)*sqrt(1/x$exp[2]+1/x$exp[1]))
        low95 = exp(log(hr) - qnorm(0.975)*sqrt(1/x$exp[2]+1/x$exp[1]))
        return(c(up95,hr,low95))
}
source("~/Dropbox/Rworkspace/r scripts/ggsurv.R")
myc = myccdc[,"myc"]
myc1 = median(myc)
myc[myc <= myc1] = 1
myc[myc>myc1] = 2
myc = ifelse(myc == 1, "Low expression","High expression")

sv = Surv(myccdc$time,as.numeric(myccdc$vitual_status))
surv = survfit(sv~myc)
logr = survdiff(sv ~ myc)
cox = coxph(sv ~ myc)
hr = hazardr(logr)

f5 = ggsurv(surv)+theme_classic() + labs(list(title = "Survival curve of MYC in TCGA Glioma",x = "Time (days)",y = "proportion survival/surviving"))+ annotate("text", label = "Log-rank P = 0.00006\nHR = 0.55(0.41-0.74)",x = 4500, y = 0.6,size = 6) + theme(legend.position = c(.7,.8),text = element_text(size = 20),legend.title = element_blank(),legend.text = element_text(size = 20))
pdf("f5.pdf")
f5
dev.off()

cdc = myccdc[,"cdc20"]
cdc1 = median(cdc)
cdc[cdc <= cdc1] = 1
cdc[cdc>cdc1] = 2
cdc = ifelse(cdc == 1, "Low expression","High expression")

sv = Surv(myccdc$time,as.numeric(myccdc$vitual_status))
surv = survfit(sv~cdc)
logr = survdiff(sv ~ cdc)
cox = coxph(sv ~ cdc)
hr = hazardr(logr)
f6 = ggsurv(surv)+theme_classic() + labs(list(title = "Survival curve of CDC20 in TCGA Glioma",x = "Time (days)",y = "proportion survival/surviving"))+ annotate("text", label = "Log-rank P = 0.00006\nHR = 0.55(0.41-0.74)",x = 4500, y = 0.6,size = 6) + theme(legend.position = c(.7,.8),text = element_text(size = 20),legend.title = element_blank(),legend.text = element_text(size = 20))
pdf("f6.pdf")
f6
dev.off()
load("../glioma_clin/tcgagbm133entrezmr.rda")

gbmmyc = exprs(gbm133)["4609",]
gbmmyc1 = median(gbmmyc)
gbmmyc[gbmmyc <= gbmmyc1] = 1
gbmmyc[gbmmyc>gbmmyc1] = 2
gbmmyc = ifelse(gbmmyc == 1, "Low expression","High expression")

sv = Surv(gbm133$time,as.numeric(gbm133$vitalstatus))
surv = survfit(sv~gbmmyc)
logr = survdiff(sv ~ gbmmyc)
cox = coxph(sv ~ gbmmyc)
hr = hazardr(logr)

f7 = ggsurv(surv)+theme_classic() + labs(list(title = "Survival curve of MYC in TCGA GBM",x = "Time (days)",y = "proportion survival/surviving"))+ annotate("text", label = "Log-rank P = 0.36\nHR = 0.91(0.75-1.10)",x = 2800, y = 0.6,size = 6) + theme(legend.position = c(.7,.8),text = element_text(size = 20),legend.title = element_blank(),legend.text = element_text(size = 20))
pdf("f7.pdf")
f7
dev.off()

gbmcdc = exprs(gbm133)["991",]
gbmcdc1 = median(gbmcdc)
gbmcdc[gbmcdc <= gbmcdc1] = 1
gbmcdc[gbmcdc>gbmcdc1] = 2
gbmcdc = ifelse(gbmcdc == 1, "Low expression","High expression")

sv = Surv(gbm133$time,as.numeric(gbm133$vitalstatus))
surv = survfit(sv~gbmcdc)
logr = survdiff(sv ~ gbmcdc)
cox = coxph(sv ~ gbmcdc)
hr = hazardr(logr)

f8 = ggsurv(surv)+theme_classic() + labs(list(title = "Survival curve of CDC20 in TCGA GBM",x = "Time (days)",y = "proportion survival/surviving"))+ annotate("text", label = "Log-rank P = 0.97\nHR = 1.003(0.83-1.22)",x = 2800, y = 0.6,size = 6) + theme(legend.position = c(.7,.8),text = element_text(size = 20),legend.title = element_blank(),legend.text = element_text(size = 20))
pdf("f8.pdf")
f8
dev.off()

gbmmyccdc = data.frame(sampleNames(gbm133),cdc = gbmcdc,myc = gbmmyc)

View(gbmmyccdc)
gbmmyccdc$both = 1
gbmmyccdc[which(gbmmyccdc$cdc == "High expression" & gbmmyccdc$myc == "Low expression"),"both"] = 3
gbmmyccdc[which(gbmmyccdc$myc == "High expression" & gbmmyccdc$cdc == "Low expression"),"both"] = 2
gbmmyccdc[which(gbmmyccdc$myc == "High expression" & gbmmyccdc$cdc == "High expression"),"both"] = 4

table(gbmmyccdc$both)
sv = Surv(gbm133$time,as.numeric(gbm133$vitalstatus))
surv = survfit(sv~gbmmyccdc$both)
logr = survdiff(sv ~ gbmmyccdc$both)
cox = coxph(sv ~ gbmmyccdc$both)
hr = hazardr(logr)

f10 = ggsurv(surv)+theme_classic() + labs(list(title = "Survival curve of CDC20 in TCGA GBM",x = "Time (days)",y = "proportion survival/surviving"))+ annotate("text", label = "Log-rank P = 0.97\nHR = 1.003(0.83-1.22)",x = 2800, y = 0.6,size = 6) + theme(legend.position = c(.7,.8),text = element_text(size = 20),legend.title = element_blank(),legend.text = element_text(size = 20))
pdf("f10.pdf")
f10
dev.off()

sampleNames(gbm133)[1:10]


xy = intersect(sampleNames(gbm133), rownames(myccdc))

length(xy)
pdf("f11.pdf")

plot(exprs(gbm133)["4609",xy],myccdc[xy,"myc"])

dev.off()

pdf("f12.pdf")
plot(exprs(gbm133)["4609",xy],exprs(gbm133)["991",xy])
dev.off()
cor.test(exprs(gbm133)["4609",xy],exprs(gbm133)["991",xy])

xy = intersect(sampleNames(gbm133),colnames(gbmmh27))

cyb5mhexp = data.frame(sample = xy,cg03826976 = gbmmh27["cg03826976",xy],CYB5R2 = exprs(gbm133)["51700",xy])

m = lm(CYB5R2~cg03826976,data = cyb5mhexp)

f14 = ggplot(cyb5mhexp,aes(x = cg03826976,y = CYB5R2))+geom_point()  + geom_abline(aes(intercept = coef(m)[1],slope = coef(m)[2]))+theme_classic() +annotate("text", label = "Linear regression P=0.0004",x = 0.6, y = 9,size = 6.5)+theme(text = element_text(size = 18))+labs(x = "cg03826976 methylation",y = "CYB5R2 expression")
pdf("figure14.pdf",width = 5,height = 5)
f14
dev.off()

rembanno <- read.csv(file = "rembanno.csv", head = T, stringsAsFactors = F,check.names = F)
library(org.Hs.eg.db)
genesynb = AnnotationDbi::select(org.Hs.eg.db, keys=c('CDC20','MYC'), columns=c("SYMBOL","ENTREZID"), keytype="SYMBOL")
genesynb$probe = rembanno[match(genesynb$ENTREZID, rembanno$Entrez),"Composite Element Name"]
rownames(genesynb) = genesynb$SYMBOL

load("/home/liu/Documents/bigfile/rembrandt_caArray_eset.rda")
rmbpda = pData(rembcarma)
View(pda)
cybclass = rembcarma[,rembcarma$disease != ""]
cybclass = cybclass[,cybclass$disease != "MIXED"]
cybclass = cybclass[,cybclass$disease != "UNCLASSIFIED"]
cybclass = cybclass[,cybclass$disease != "UNKNOWN"]

table(rmbpda[which(rmbpda$disease == "GBM"),"grading"])
table(rmbpda$grading, rmbpda$disease)
table(rmbpda$disease)
rmbpda[which(rmbpda$disease == "GBM"),"grading"] = "Grade 4"
rmbpda[which(rmbpda$disease == "NON_TUMOR"),"grading"] = "NON_TUMOR"
table(rmbpda[which(rmbpda$disease == "NON_TUMOR"),"grading"])
pData(rembcarma) = rmbpda
cybgrading = rembcarma[,rembcarma$grading != ""]
mycdc = data.frame(sample = sampleNames(cybclass),class =cybclass$disease, "MYC"= exprs(cybclass)[genesynb["MYC","probe"],], "CDC20" = exprs(cybclass)[genesynb["CDC20","probe"],])
mycdcg = data.frame(sample = sampleNames(cybgrading),grading =cybgrading$grading, "MYC"= exprs(cybgrading)[genesynb["MYC","probe"],], "CDC20" = exprs(cybgrading)[genesynb["CDC20","probe"],])
table(cyb5g$grading)

#cyb5c = cyb5c %>% filter(class %in% c("NON_TUMOR","ASTROCYTOMA","OLIGODENDROGLIOMA","GBM"))
f14 = ggplot(mycdc,aes(x = class,y = MYC,colors = class))+geom_boxplot(width = 0.7,outlier.shape = NA)  + geom_jitter(aes(color = class),position = position_jitter(width = 0.2))+scale_x_discrete(limits = c("NON_TUMOR","ASTROCYTOMA","OLIGODENDROGLIOMA","GBM"),breaks = c("NON_TUMOR","ASTROCYTOMA","OLIGODENDROGLIOMA","GBM"),labels = c("NON TUMOR","A/AA","O/AO","GBM"))+theme_classic()+theme(legend.position="none",text = element_text(size = 16),axis.text.x = element_text(size = 11))+ylab("MYC expression")
f15 = ggplot(mycdcg,aes(x = grading,y = MYC,colors = grading))+geom_boxplot(width = 0.7,outlier.shape = NA)  + geom_jitter(aes(color = grading),position = position_jitter(width = 0.2))+scale_x_discrete(limits = c("NON_TUMOR","Grade 2","Grade 3","Grade 4"),breaks = c("NON_TUMOR","Grade 2","Grade 3","Grade 4"),labels = c("NON TUMOR","Grade II","Grade III","Grade IV"))+theme_classic()+theme(legend.position="none",text = element_text(size = 16),axis.text.x = element_text(size = 10))+ylab("MYC expression")+xlab(NULL)
pdf("f14.pdf")
f14
dev.off()
pdf("f15.pdf")
f15
dev.off()
f16 = ggplot(mycdc,aes(x = class,y = CDC20,colors = class))+geom_boxplot(width = 0.7,outlier.shape = NA)  + geom_jitter(aes(color = class),position = position_jitter(width = 0.2))+scale_x_discrete(limits = c("NON_TUMOR","ASTROCYTOMA","OLIGODENDROGLIOMA","GBM"),breaks = c("NON_TUMOR","ASTROCYTOMA","OLIGODENDROGLIOMA","GBM"),labels = c("NON TUMOR","A/AA","O/AO","GBM"))+theme_classic()+theme(legend.position="none",text = element_text(size = 16),axis.text.x = element_text(size = 11))+ylab("CDC20 expression")
f17 = ggplot(mycdcg,aes(x = grading,y = CDC20,colors = grading))+geom_boxplot(width = 0.7,outlier.shape = NA)  + geom_jitter(aes(color = grading),position = position_jitter(width = 0.2))+scale_x_discrete(limits = c("NON_TUMOR","Grade 2","Grade 3","Grade 4"),breaks = c("NON_TUMOR","Grade 2","Grade 3","Grade 4"),labels = c("NON TUMOR","Grade II","Grade III","Grade IV"))+theme_classic()+theme(legend.position="none",text = element_text(size = 16),axis.text.x = element_text(size = 10))+ylab("CDC20 expression")+xlab(NULL)
pdf("f16.pdf")
f16
dev.off()
pdf("f17.pdf")
f17
Gdev.off()

x = read.table(file = "../tcgafile/gdac.broadinstitute.org_GBMLGG.Mutation_Packager_Raw_Calls.Level_3.2015060100.0.0/TCGA-02-0003-01.maf.txt", header = T, sep = "\t")
View(x)
dim(x)
write.csv(x, file = "muttest.csv")

x1 = read.csv(file = "../tcgafile/TCGA_GBM_mutation_broad_gene-2015-02-24/genomicMatrix.csv", header = T)
write.csv(x1,file = "muttest1.csv")

glmcna = read.table(file = "../tcgafile/GBMLGG.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_hg19__seg.seg.txt", head = T, sep = "\t",quote = "")
dim(glmcna)
head(glmcna)

load("glioma_gene_cna.rda")
load("glioma_gene_cnaall.rda")
p16 = glmgenecnall %>% filter(gene == "1029")
sum(is.na(p16$segmean))
hist(p16$segmean)
p16 = p16[which(substr(p16$sample,14,15) == "01"),]
sum(p16$segmean < -1)
p16in = glmgenecn %>% filter(gene == "1029")

p16in = p16in[which(substr(p16in$sample,14,15) == "01" |substr(p16in$sample,14,15) == "02"),]
p16in1 = p16in[which(substr(p16in$sample,14,15) == "01"),]

rownames(p16) = substr(p16$sample,1,12)
length(unique(substr(p16in$sample,1,12)))
gbmrn = exprs(gbm133)
xy = intersect(rownames(p16in1),colnames(gbmrn))
xydel = rownames(p16in1[which(p16in1[,"segmean"] <= -1 & rownames(p16in1) %in% xy),])
xyndel = rownames(p16in1[which(p16in1[,"segmean"] > -1 & rownames(p16in1) %in% xy),])
expset = gbm133[,xydel]
entz = "4609"
gene= "MYC"
set= "TCGA p16 del GBM"
file = "f18"
f18 = genesurv(expset = expset,entz= "4609",gene = "MYC",set = "TCGA p16 del GBM")
pdf("f18.pdf")
f18
dev.off()
f19 = genesurv(expset = expset,entz= "991",gene = "CDC20",set = "TCGA p16 del GBM")
pdf("f19.pdf")
f19
dev.off()
gbmcdnk = gbm133[,xy]
gbmcdnk$cdkn = ifelse(sampleNames(gbmcdnk) %in% xydel,"del","nondel")

sv = Surv(gbmcdnk$time,as.numeric(gbmcdnk$vitalstatus))
surv = survfit(sv~gbmcdnk$cdkn)
logr = survdiff(sv ~ gbmcdnk$cdkn)
cox = coxph(sv ~ gbmcdnk$cdkn)
hr = hazardr(logr)
    pval = round(1- pchisq(logr$chisq,1),3)
    lb = paste("Log-rank P = ",pval,"\nHR = ", round(hr[2],2),"(",round(hr[1],2),"-", round(hr[3],2), ")",sep = "")
    tl = "Survival curve of p16 del in TCGA GBM"
    f0 = ggsurv(surv)+theme_classic() + labs(list(title = tl,x = "Time (days)",y = "proportion survival/surviving"))+ annotate("text", label = lb ,x = 0.72*max(expset$time), y = 0.6,size = 6) + theme(legend.position = c(.7,.8),text = element_text(size = 16),legend.title = element_blank(),legend.text = element_text(size = 20))
pdf("f20.pdf")
f0
dev.off()

gbmcdkn = data.frame(p16 = gbmcdnk$cdkn,cdc20 = exprs(gbmcdnk)["991",],myc = exprs(gbmcdnk)["4609",])
t.test(gbmcdkn$myc ~ gbmcdkn$p16)
f21 = ggplot(gbmcdkn,aes(x = p16,y = myc,colors = p16))+geom_boxplot(width = 0.7,outlier.shape = NA)  + geom_jitter(aes(color = p16),position = position_jitter(width = 0.2))+annotate("text", label = "P=0.006",x = 1.5, y = 12,size = 6.5)+theme_classic()+theme(legend.position="none",text = element_text(size = 16),axis.text.x = element_text(size = 11))+ylab("MYC expression")

pdf("f21.pdf")

f21
dev.off()
t.test(gbmcdkn$cdc20 ~ gbmcdkn$p16)
f22 = ggplot(gbmcdkn,aes(x = p16,y = cdc20,colors = p16))+geom_boxplot(width = 0.7,outlier.shape = NA)  + geom_jitter(aes(color = p16),position = position_jitter(width = 0.2))+annotate("text", label = "P=0.57",x = 1.5, y = 12,size = 6.5)+theme_classic()+theme(legend.position="none",text = element_text(size = 16),axis.text.x = element_text(size = 11))+ylab("CDC20 expression")
pdf("f22.pdf")
f22
dev.off()
du = substr(p16[duplicated(substr(p16$sample,1,12)),"sample"],1,12)
pdu = p16[which(substr(p16$sample,1,12) %in% du),]
p16 = p16[order(abs(p16$segmean),decreasing = T),]
p16 = p16[!duplicated(p16$sample),]
p16$sample = substr(p16$sample,1,12)

rownames(p16) = p16$sample

rownames(glmclin) = glmclin[,1]
p16g2del = intersect(rownames(p16[which(p16$segmean < -1),]),rownames(glmclin[which(glmclin$grade == "g2"),]))
p16g3del = intersect(rownames(p16[which(p16$segmean < -1),]),rownames(glmclin[which(glmclin$grade == "g3"),]))
p16g4del = intersect(rownames(p16[which(p16$segmean < -1),]),rownames(glmclin[which(glmclin$grade == "g4"),]))

## make the figure ji asked

load("glioma_gene_cnaall.rda")
ji = glmgenecnall %>% filter(gene == "1029" | gene == "5728" | gene == "10000")
ji = ji[which(substr(ji$sample,14,15) == "01"),]
ji = ji[order(abs(ji$segmean),decreasing = T),]

ji = ji %>% distinct(sample,gene) %>% spread(gene,segmean)
ji$sample = substr(ji$sample,1,12)
rownames(ji) = ji$sample
gbmrn = exprs(gbm133)
xy = intersect(rownames(ji),colnames(gbmrn))
jiall = data.frame(ji[xy,], cdc20 = gbmrn["991",xy], myc = gbmrn["4609",xy])
colnames(jiall)[2:4]  = c("akt3","p16","pten")

x1 = read.csv(file = "../tcgafile/TCGA_GBM_mutation_broad_gene-2015-02-24/genomicMatrix.csv", header = T,check.names = F)
colnames(x1) = substr(colnames(x1),1,12)
jimu = x1[which(x1[,1] == "AKT3" | x1[,1] == "PTEN"),]
jimu = t(jimu)
colnames(jimu) = jimu[1,]
jimu = jimu[-1,]
colnames(jimu) = paste(colnames(jimu),"mu",sep = "")
jimu = data.frame(sample = rownames(jimu),jimu,stringsAsFactors = F)
jif = left_join(jiall,jimu)
jif = jif[,-7]
save(jif, file = "jifigure.rda")

## make figure

jifmu = jif[which(!is.na(jif$PTENmu)),]
jifmu[which(jif$akt3 <= -1),"akt3"] = -1
jifmu[which(jif$akt3 >= 1),"akt3"] = 1
jifmu[which(jif$akt3 < 1 $ jif$akt3 > -1),"akt3"] = 0
jifmu[which(jif$p16 <= -1),"p16"] = -1
jifmu[which(jif$p16 >= 1),"p16"] = 1
jifmu[which(jif$p16 < 1 $ jif$p16 > -1),"p16"] = 0
jifmu[which(jif$pten <= -1),"pten"] = -1
jifmu[which(jif$pten >= 1),"pten"] = 1
jifmu[which(jif$pten < 1 $ jif$pten > -1),"pten"] = 0

sum(jif$akt3 > 1)
