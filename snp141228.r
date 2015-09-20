load("gbmallhap.rda")
load("hap550annothg19.rda")
load("hap550annotTCGA.rda")
load("xinfo.rda")
library(rtracklayer)


hapinfo = xinfo[,c("marker","chr","chrpos","current.rsid")]

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


hapinfo1$ABgenotype = chartr("ATCG","TAGC",hapinfo1$ABgenotype)

View(head(hapinfo1))
chain <- import.chain("~/annotation_data/hg38ToHg19.over.chain")
hap550gr <- with(hapinfo1, GRanges(chr, IRanges(as.numeric(as.character(chrpos)), as.numeric(as.character(chrpos)),names = marker), strand = "+", marker,ABgenotype))
hap550grhg19 <- liftOver(hap550gr, chain)
hapgrhg19 = unlist(hap550grhg19)
haphg19 <- data.frame(marker=mcols(hapgrhg19)$marker, chr=seqnames(hapgrhg19),chrpos=start(hapgrhg19), current.rsid=mcols(hapgrhg19)$marker, ABgenotype = mcols(hapgrhg19)$ABgenotype, stringsAsFactors = F)
save(haphg19,file = "haphg19annot.rda")
haphg19[]
haphg19$Agenotype = sapply(strsplit(haphg19$ABgenotype,split = "/"),"[[",1)
haphg19$Bgenotype = sapply(strsplit(haphg19$ABgenotype,split = "/"),"[[",2)

xy = intersect(rownames(hapsnp),haphg19$marker)
rownames(haphg19) = haphg19$marker
hapsnp19 = hapsnp[rownames(haphg19),]
all(rownames(hapsnp19) == rownames(haphg19))
hapgbm19 = data.frame(haphg19[,c("marker","chr","chrpos","ABgenotype","Agenotype","Bgenotype")],hapsnp19)
colnames(hapgbm19) = gsub("\\.","-",colnames(hapgbm19))
table(hapgbm19$chr)

h = hapgbm19
#h$chr = paste("chr",h$chr)
for(i in 7:ncol(h)){
	h[which(h[,i] == 0),i] = paste(h[which(h[,i] == 0),"Agenotype"],h[which(h[,i] == 0),"Agenotype"],sep = "")
	h[which(h[,i] == 1),i] = paste(h[which(h[,i] == 1),"Agenotype"],h[which(h[,i] == 1),"Bgenotype"],sep = "")
	h[which(h[,i] == 2),i] = paste(h[which(h[,i] == 2),"Bgenotype"],h[which(h[,i] == 2),"Bgenotype"],sep = "")
	if(i%%10 == 0) writeLines(paste("Calculating ", i,",",round(i/ncol(h)*100, digits = 2), "% done.", sep = ""))
}

h[is.na(h)] = "--"
h = h[,-(4:6)]
View(head(h))
hapgbm19 = h
save(hapgbm19, file = "hapgbmhg19correctwith1000g.rda")
colnames(h)[1:3] = c("Name","Chr","Pos")

clingbm = read.csv(file = "~/tcgadata/Clinical and Molecular Subclass Data Table.csv")
xny = intersect(colnames(h)[4:ncol(h)],clingbm$Case.ID)
length(xny)
xcy = setdiff(colnames(h)[4:ncol(h)],clingbm$Case.ID)
length(xcy)


hapmeta = clingbm[which(clingbm[,1] %in% xny),c(1,4,3)]
dim(hapmeta)
colnames(hapmeta) =  c("id","sex","age")
hapmeta[which(hapmeta$sex == "MALE"),"sex"] = "1"
hapmeta[which(hapmeta$sex == "FEMALE"),"sex"] = "0"
rownames(hapmeta) = hapmeta$id
hapmeta = hapmeta[xny,]

hapmismeta = data.frame(id = xcy,sex = hapmeta$sex[1:16], age = hapmeta$age[1:16])
hapallmeta = rbind(hapmeta, hapmismeta)
View(hapallmeta)
rownames(hapallmeta) = hapallmeta$id
hapallmeta = hapallmeta[colnames(h)[4:ncol(h)],]
all(colnames(h)[4:ncol(h)] == hapallmeta$id)
hapallmeta = read.table(file = "happheno.txt",head = T,stringsAsFactors = F)
write.table(hapallmeta, file = "happheno.txt",row.names = F,quote = F)

for (i in names(table(h$Chr))){
	x = h[which(h$Chr == i),]
	x = x[order(x$Pos),]
	name = paste("hapgbmhg19chr",i,sep = "")
	filename = paste(name,".txt",sep = "")
	raw = paste(name,".raw",sep = "")
	#ped = paste(paste("merlin",name,sep = ""),".ped",sep = "")
	#dat = paste(paste("merlin",name,sep = ""),".dat",sep = "")
	#map = paste(paste("merlin",name,sep = ""),".map",sep = "")
	sample = paste(paste("impute",name,sep = ""),".sample",sep = "")
	strand = paste(paste("impute",name,sep = ""),".strand",sep = "")
	geno = paste(paste("impute",name,sep = ""),".gen",sep = "")
	write.table(x, file = filename,row.names = F,quote = F)
	library(GenABEL)
    convert.snp.illumina(inf = filename,out = raw,strand = "+")
    df <- load.gwaa.data(phe="happheno.txt", gen=raw, force=TRUE)
    export.impute(df,genofile=geno,samplefile=sample,strandfile=strand,cachesizeMb=4096)
    #export.merlin(df,pedfile = ped, datafile = dat, mapfile = map)
}
i = "chr1"
x = x[order(x$Pos),]

for i in $(seq 1 22); do (shapeit -check -G imputehapgbmhg19chrchr$i\.gen imputehapgbmhg19chrchr$i\.sample -M ./phase1nosing/genetic_map_chr$i\_combined_b37.txt --input-ref ./phase1nosing/ALL.chr$i\.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nosing.haplotypes.gz ./phase1nosing/ALL.chr$i\.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nosing.legend.gz ./phase1nosing/ALL.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.sample --output-log gbmhapchr$i\.alignments); done

setwd("~/1000GP_Phase3")
dir()
for (i in 2:22){
	#name = paste("hapgbmhg19chrchr",i,sep = "")
	#filename = paste(name,".txt",sep = "")
	strand = paste("gbmhapchr",i,sep = "")
	strandfile = paste(strand,".alignments.snp.strand",sep = "")
	wrong = read.table(strandfile,head = T,row.names = NULL,stringsAsFactors=F,check.names= F)
	wrong = wrong[which(wrong$type == "Strand"),]
	wrong$chr = i
	#x= wrong[which(wrong$main_id == wrong$ref_id),]
	wrong1 = rbind(wrong1,wrong)
	#rownames(x) = x[,"main_id"]
    #file = read.table(filename,head= T,stringsAsFactors=F,check.names= F)
    #rownames(file) = file$Name
    #row = x[,"main_id"]
    #for(i in 1:length(row)){
    #	file[row[i],(4:ncol(file))] = gsub(x[row[i],"main_A"],x[row[i],"ref_A"],file[row[i],(4:ncol(file))])
    #    file[row[i],(4:ncol(file))] = gsub(x[row[i],"main_B"],x[row[i],"ref_B"],file[row[i],(4:ncol(file))])
    #    if(i%%500 == 0) writeLines(paste("Calculating ", i,",",round(i/length(row)*100, digits = 2), "% done.", sep = "")) 
    #}
    

}
i=1
wrong = wrong1
x= wrong[which(wrong$main_id == wrong$ref_id),]
dim(x)
table(nchar(x$ref_B))
x = nchar(x[,"ref_B"])

setwd("~/Dropbox/Rworkspace/projects/phlpp2snp")

save(wrong,file = "genouncorrectwith1000g.rda")
x= wrong[which(wrong$main_id == wrong$ref_id),]
x$ABgenotype = paste(x$ref_A,x$ref_B,sep = "/")
rownames(x) = x$main_id
length(unique(x$main_id))
dim(x)
load("haphg19annot.rda")
rownames(haphg19) = haphg19$marker
haphg19[rownames(x),"ABgenotype"] = x[rownames(x),"ABgenotype"]
save(hapgbm19, file = "hapgbmhg19corrwith1000g.rda")


