load("hap550annotTCGA.rda")
library(rtracklayer)
hap550gr <- with(hapannot, GRanges(chromosome, IRanges(as.numeric(as.character(location)), as.numeric(as.character(location)),names = Reporter), strand = "+", gene, Reporter, Sequence = hapannot$"Reporter Sequence"))
chain <- import.chain("~/annotation_data/hg18ToHg19.over.chain")
hap550grhg19 <- liftOver(hap550gr, chain)
hapgrhg19 = unlist(hap550grhg19)
haphg19 <- data.frame(names=mcols(hapgrhg19)$Reporter, chromosome=seqnames(hapgrhg19), location=start(hapgrhg19), strands=strand(hapgrhg19), sequence = mcols(hapgrhg19)$Sequence,stringsAsFactors = F)
save(haphg19,file = "hap550annothg19.rda")
haphg19$chromosome = gsub("chr","",haphg19$chromosome)

nam = c(1:22,"mt","X","Y")
nam[1]
i = 1
for (i in 2:length(nam)){
	name = paste("hapgbmhg38chr",nam[i],sep = "")
	filename = paste(name,".txt",sep = "")
	data = read.table(file = filename,head = T,check.name = F)
	datax = rbind(datax,data)
	#data$Chr = haphg19[match(data$Name, haphg19$names),"chromosome"]
	#data$Pos = haphg19[match(data$Name, haphg19$names),"location"]
	#raw = paste(name,".raw",sep = "")
	#ped = paste(paste("merlin",name,sep = ""),".ped",sep = "")
	#dat = paste(paste("merlin",name,sep = ""),".dat",sep = "")
	#map = paste(paste("merlin",name,sep = ""),".map",sep = "")
	#write.table(x, file = filename,row.names = F,quote = F)
	#library(GenABEL)
    #convert.snp.illumina(inf = filename,out = raw,strand = "+")
    #df <- load.gwaa.data(phe="happheno.txt", gen=raw, force=TRUE)
    #export.merlin(df,pedfile = ped, datafile = dat, mapfile = map)
}

gbmhapgh19 = datax
hapgr = datax[,1:3]
hapgr$Chr = toupper(hapgr$Chr)
hapgr$Chr = paste("chr",hapgr$Chr,sep = "")
hapgr$Chr = sub("chrMT","chrM",hapgr$Chr)

hap550gr38 <- with(hapgr, GRanges(Chr, IRanges(as.numeric(as.character(Pos)), as.numeric(as.character(Pos)),names = Name), strand = "+", Name))
#hap550gr38
chain <- import.chain("~/annotation_data/hg38ToHg19.over.chain")

hap550gr19 <- liftOver(hap550gr38, chain)
hap19 = unlist(hap550grhg19)
hapgrhg19 = unlist(hap550gr19)
haphg19 <- data.frame(names=mcols(hapgrhg19)$Name, chromosome=seqnames(hapgrhg19), location=start(hapgrhg19), strands=strand(hapgrhg19), stringsAsFactors = F)
hap1819 = data.frame(names=mcols(hap19)$Reporter, chromosome=seqnames(hap19), location=start(hap19), strands=strand(hap19), sequence = mcols(hap19)$Sequence,stringsAsFactors = F)
#View(haphg19)
#View(hap1819)
#xy = intersect(hap1819$names,haphg19$names)
#length(xy)
#x = haphg19[xy,]
#y = hap1819[xy,]
#all(haphg19$names == hapgr$Name)
haphg19$pos18 = hap1819[match(haphg19$names, hap1819$names),"location"]
View(haphg19)
all(haphg19$location == haphg19$pos18)
datax$Chr = haphg19[match(datax$Name, haphg19$names),"chromosome"]
datax$Pos = haphg19[match(datax$Name, haphg19$names),"location"]
gbmhaphg19 = datax
save(gbmhaphg19,file = "GBMhaphg19.rda")
happheno = read.table(file = "happheno.txt",stringsAsFactors = F, head= T)
all(colnames(datax)[4:ncol(datax)] == happheno$id)
for (i in names(table(h$Chr))){
	x = h[which(h$Chr == i),]
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
    export.impute(df,genofile=geno,samplefile=sample,strandfile=strand,cachesizeMb=2048)
    #export.merlin(df,pedfile = ped, datafile = dat, mapfile = map)
}
for m in $(ls *.gz); do (tar -zxvf $m); done
hapannot["rs2342694",]