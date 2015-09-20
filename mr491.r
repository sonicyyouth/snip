#library(mirbase.db)
#ls("package:mirbase.db")
#qcdata = capture.output(mirbase())
#class(qcdata)
#qcdata[1:50]
#mirbaseMAPCOUNTS
#mir = ls(mirbaseCHRLOC)
#mir[grep("mir-491",mir)]
#mir491loc = unlist(mget("hsa-mir-491",mirbaseCHRLOC))
#class(mir491loc)

library(biomaRt)
ensembl = useMart("ensembl",dataset = "hsapiens_gene_ensembl")
#listDatasets(ensembl)
#pages = attributePages(ensembl)
#filters = listFilters(ensembl)
#attrib = listAttributes(ensembl)
#dim(attrib)
#filters[grep("mir",filters[,2]),]
x = getBM(attributes = c("chromosome_name","start_position","end_position"),filters = "mirbase_id", values = "hsa-mir-491",mart= ensembl)
gene = getBM(attributes = c("hgnc_symbol", "chromosome_name","start_position","end_position"),filters = c("chromosome_name","start","end") ,values = list(x[1,1],x[1,2],x[1,3]), mart = ensembl)
gene
#hgnc_symbol chromosome_name start_position end_position
#1      MIR491               9       20716100     20716209
#2       FOCAD               9       20658309     20995955
getwd()
setwd("/home/liu/Dropbox/Rworkspace")
rm(list=ls())
source("./TCGA-Assembler/Module_A.r")
source("./TCGA-Assembler/Module_B.r")
Methylation450Data = ProcessMethylation450Data(inputFilePath ="/home/liu/tcgaorginfile/LGG.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt", outputFileName ="LGG_humanmethylation450_assembler", outputFileFolder = "/home/liu/tcgadata",fileSource = "Firehose");

Methylation450DataGBM = ProcessMethylation450Data(inputFilePath ="/home/liu/tcgaorginfile/GBM.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt", outputFileName ="GBM_humanmethylation450_assembler", outputFileFolder = "/home/liu/tcgadata",fileSource = "Firehose");
Methylation27DataGBM = ProcessMethylation450Data(inputFilePath ="/home/liu/tcgaorginfile/GBM.methylation__humanmethylation27__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt", outputFileName ="GBM_humanmethylation27_assembler", outputFileFolder = "/home/liu/tcgadata",fileSource = "Firehose");
rowgdac2rds = function(orginfile){
	outfile = strsplit(orginfile,split = "__")[[1]][1]
	inloc = "~/tcgaorginfile"
	outloc = "~/tcgadata/"
	fulloutname = paste(outloc,paste(outfile,".rds",sep= ""),sep ="")
	fullinname = paste(inloc,orginfile,sep = "/")
	rowdata = read.table(file = fullinname,head = T, sep = "\t", as.is = T,quote = "",stringsAsFactors = F, check.names = F)
	if(length(grep("mirnaseq", outfile))==1){
		rowdata = rowdata[,c(1,grep("reads_per_million",rowdata[1,]))]
	}
	rowdata = rowdata[-1,]
	rownames(rowdata) = rowdata[,1]
	rowdata = rowdata[,-1]
	colnames(rowdata) = substr(colnames(rowdata),1,15)
	saveRDS(rowdata,file = fulloutname)
}
rowgdac2rds("GBM.transcriptome__agilentg4502a_07_1__unc_edu__Level_3__unc_lowess_normalization_gene_level__data.data.txt")
rowgdac2rds("GBM.transcriptome__agilentg4502a_07_2__unc_edu__Level_3__unc_lowess_normalization_gene_level__data.data.txt")
rowgdac2rds("LGG.mirnaseq__illuminahiseq_mirnaseq__bcgsc_ca__Level_3__miR_gene_expression__data.data.txt")
rowgdac2rds("LGG.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt")
rowgdac2rds("GBM.transcriptome__ht_hg_u133a__broad_mit_edu__Level_3__gene_rma__data.data.txt")
rowgdac2rds("GBM.mirna__h_mirna_8x15k__unc_edu__Level_3__unc_DWD_Batch_adjusted__data.data.txt")
gbmrna4502a1 = readRDS("~/tcgadata/GBM.transcriptome4502a1.rds")
gbmrna4502a2 = readRDS("~/tcgadata/GBM.transcriptome4502a2.rds")
x = intersect(rownames(gbmrna4502a1),rownames(gbmrna4502a2))
gbmrna4502a = data.frame(gbmrna4502a1[x,],gbmrna4502a2[x,])
colnames(gbmrna4502a) = gsub("\\.","-",colnames(gbmrna4502a))

saveRDS(gbmrna4502a,file = "~/tcgadata/GBM.transcriptome4502a.rds")

lggmrna = readRDS("~/tcgadata/LGG.rnaseqv2.rds")
lggmirna = readRDS("~/tcgadata/LGG.mirnaseq.rds")
gbmmrna = readRDS("~/tcgadata/GBM.transcriptome133a.rds")
gbmmrna4502 = readRDS("~/tcgadata/GBM.transcriptome4502a.rds")
y = setdiff(rownames(gbmrna4502a2),rownames(gbmrna4502a1))

dim(gbmrna4502a1)
dim(gbmrna4502a2)
gbmmirna = readRDS("~/tcgadata/GBM.mirna.rds")
load("~/tcgadata/LGG_humanmethylation450_assembler.rda")
lggmethy= cbind(as.data.frame(Des),as.data.frame(Data))
lggmethy[,"CoordinateID"] = as.numeric(as.character(lggmethy[,"CoordinateID"]))

ink4a = lggmethy[which(lggmethy[,2] == "CDKN2A"),]
ink4agbm = gbmmethy450[which(gbmmethy450[,2] == "CDKN2A"),]

hist(as.numeric(ink4agbm[,5:ncol(ink4agbm)]))

hist(as.numeric(ink4a[1,5:ncol(ink4a)]))
focad = lggmethy[which(lggmethy$ChromosomeID == "9" & lggmethy$CoordinateID > (gene[2,3] - 5000) & lggmethy$CoordinateID < gene[2,4]),] 

genemethyhist = function(gene){
	mgenelgg = lggmethy[which(lggmethy[,2] == gene),]
	mgenegbm = gbmmethy450[which(gbmmethy450[,2] == gene),]
	mainlgg = paste(paste("Methylation status of",gene),"in LGG")
	maingbm = paste(paste("Methylation status of",gene),"in GBM")
	outlgg = paste(gene,"lgg.pdf",sep = "")
	outgbm = paste(gene,"gbm.pdf",sep = "")

	pdf(outlgg)
    hist(as.numeric(mgenelgg[4,5:ncol(mgenelgg)]),main = mainlgg, xlab = "Methylation values",ylab = "Sample Frequency")
    dev.off()
    pdf(outgbm)
    hist(as.numeric(mgenegbm[4,5:ncol(mgenegbm)]),main = maingbm, xlab = "Methylation values",ylab = "Sample Frequency")
    dev.off()
    return(list(mgenelgg,mgenegbm))
}
x = genemethyhist("AJAP1")
y = genemethyhist("MIIP")
gene = "MIIP"
genemethyhist("AJAP1")
dim(mgenelgg)
ink4a = lggmethy[which(lggmethy[,2] == "CDKN2A"),]
save(gbmmethy450,file = "GBM_humanmethylation450.rda")
save(lggmethy,file = "LGG_methylation450.rda")
ink4agbm = gbmmethy450[which(gbmmethy450[,2] == "CDKN2A"),]
mainlgg
ls()
pdf("lggccdk.pdf")
hist(as.numeric(ink4a[1,5:ncol(ink4a)]),main = "Methylation status of CDKN2A in LGG", xlab = "Methylation values",ylab = "Sample Frequency")
dev.off()

pdf("gbmcdk.pdf")
hist(as.numeric(ink4agbm[1,5:ncol(ink4agbm)]),main = "Methylation status of CDKN2A in GBM", xlab = "Methylation values",ylab = "Sample Frequency")
dev.off()

pdf("lggcg01947138.pdf")
hist(as.numeric(focad[1,5:ncol(focad)]),main = "Methylation status of cg01947138 in LGG", xlab = "Methylation values",ylab = "Sample Frequency")
dev.off()

pdf("lggcg13517516.pdf")
hist(as.numeric(focad[3,5:ncol(focad)]),main = "Methylation status of cg13517516 in LGG", xlab = "Methylation values",ylab = "Sample Frequency")
dev.off()

hist(as.numeric(focad[2,5:ncol(focad)]),main = "methylation status of ")
focadgbm[,1:4]

dim(focad)

mir491[,1:10]
load("~/tcgadata/GBM_humanmethylation450_assembler.rda")
gbmmethy450 = cbind(as.data.frame(Des),as.data.frame(Data))
focadgbm = gbmmethy450[which(gbmmethy450$ChromosomeID == "1" & gbmmethy450$CoordinateID > (gene[2,3] - 5000) & gbmmethy450$CoordinateID < gene[2,4]),] 

hist(as.numeric(focadgbm[3,5:ncol(focadgbm)]))

pdf("GBMcg13134607.pdf")
hist(as.numeric(focadgbm[2,5:ncol(focadgbm)]),main = "Methylation status of cg13134607 in GBM", xlab = "Methylation values",ylab = "Sample Frequency")
dev.off()

pdf("GBMcg01947138.pdf")
hist(as.numeric(focadgbm[1,5:ncol(focadgbm)]),main = "Methylation status of cg01947138 in GBM", xlab = "Methylation values",ylab = "Sample Frequency")
dev.off()

pdf("GBMcg13517516.pdf")
hist(as.numeric(focadgbm[3,5:ncol(focadgbm)]),main = "Methylation status of cg13517516 in GBM", xlab = "Methylation values",ylab = "Sample Frequency")
dev.off()

ink4a = lggmethy[grep("IIP",lggmethy[,2]),]

lggmrna = readRDS("~/tcgadata/LGG.rnaseqv2.rds")
lggmirna = readRDS("~/tcgadata/LGG.mirnaseq.rds")
gbmmrna = readRDS("~/tcgadata/GBM.transcriptome133a.rds")
gbmmirna = readRDS("~/tcgadata/GBM.mirna.rds")
gbmmrna4502 = readRDS("~/tcgadata/GBM.transcriptome4502a.rds")

xylgg = intersect(colnames(lggmrna),colnames(lggmirna))
xygbm = intersect(colnames(gbmrna4502a),colnames(gbmmirna))
colnames(gbmrna4502a)[1:10]
cor.test(as.numeric(lggmirna["hsa-mir-491",xylgg]),as.numeric(lggmrna["KIAA1797|54914",xylgg]),method= "spearman")
plot(lggmirna["hsa-mir-491",xylgg],lggmrna["KIAA1797|54914",xylgg])
rownames(gbmmrna)[grep("KIAA",rownames(gbmmrna))]

cor.test(log2(as.numeric(lggmirna["hsa-mir-491",xylgg])),log2(as.numeric(lggmrna["KIAA1797|54914",xylgg])))

plot(log2(as.numeric(lggmirna["hsa-mir-491",xylgg])),log2(as.numeric(lggmrna["KIAA1797|54914",xylgg])))

hist(log2(as.numeric(lggmirna["hsa-mir-491",xylgg])))
hist(log2(as.numeric(lggmrna["KIAA1797|54914",xylgg])))

hist(as.numeric(lggmirna["hsa-mir-491",xylgg]))
hist(as.numeric(lggmrna["KIAA1797|54914",xylgg]))
lgg491 = log2(as.numeric(lggmirna["hsa-mir-491",xylgg]))
lggfocad = log2(as.numeric(lggmrna["KIAA1797|54914",xylgg]))
pdf("lgg491focad.pdf")
plot(lgg491~lggfocad,main = "expression correlation of miR-491 and FOCAD in LGG", xlab = "expression of miR-491",ylab = "expression of FOCAD")
abline(lm(lgg491~lggfocad))
cor.test(lgg491,lggfocad)
dev.off()
gbmfocad = as.numeric(gbmrna4502a[grep("KIAA1797",rownames(gbmrna4502a)),xygbm])
gbm491 = log2(as.numeric(gbmmirna["hsa-miR-491",xygbm]))
pdf("gbm491focad.pdf")
plot(gbm491~gbmfocad,main = "expression correlation of miR-491 and FOCAD in GBM", xlab = "expression of miR-491",ylab = "expression of FOCAD")
)
abline(lm(gbm491~gbmfocad))

cor.test(gbm491,gbmfocad)
dev.off()

x[x == 0]

gbm4502 = read.table(file = "~/tcgaorginfile/GBM.transcriptome__agilentg4502a_07_2__unc_edu__Level_3__unc_lowess_normalization_gene_level__data.data.txt",head = T, sep = "\t", as.is = T,stringsAsFactors = F, check.names = F)