getwd()
setwd("/home/liu/Dropbox/Rworkspace/projects/phlpp2snp/tcgarnabam")
library(VariantAnnotation)
source("http://bioconductor.org/biocLite.R")
biocLite("SNPlocs.Hsapiens.dbSNP.20120608")
biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")

library(org.Hs.eg.db)
enid = select(org.Hs.eg.db, keys="PHLPP2", columns="ENTREZID", keytype="SYMBOL")$ENTREZID
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
ex = exons(txdb)
phlpp = select(txdb, keys=enid, columns=columns(txdb), keytype="GENEID")
phlppexon = unique(phlpp$EXONID)
phlppexonGrange = ex[which(ex$exon_id %in% phlppexon),]
phlppGrange = reduce(phlppexonGrange)
ranges(phlppGrange)
str(phlppexonGrange)
phlpp = phlpp[order(phlpp$EXONSTART),]
?plotRanges

library(biomaRt)
listDatasets(useMart("snp"))
snpmart = useMart("snp", dataset="hsapiens_snp")
listFilters(snpmart)
listAttributes(snpmart)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb = Txls(Db.Hsapiens.UCSC.hg19.knownGene
keytypes(txdb)
columns(txdb)

input1 = getBM(c("refsnp_id","chr_name","chrom_start","ensembl_gene_stable_id"), filters = "chromosomal_region",values = list(16,71677000,71759000),mart = snpmart)

library(SNPlocs.Hsapiens.dbSNP.20120608)

ch16snp = getSNPlocs("ch16",as.GRanges = T)
phlppsnp = ch16snp[findOverlaps(ranges(phlppGrange),ranges(ch16snp))@subjectHits]
strsave(phlppsnp,file = "phlpp2snp.rda")
x@subjectHits
71677000-71759000
phlpp2snp = ch16snp[which(ch16snp$loc > 71677000 & ch16snp$loc < 71759000),]
head(phlpp2snp)
nrow(phlpp2snp)
gr = rsidsToGRanges("rs10500560")
loc = start(gr)


setwd("/media/liu/My Passport/bamfile")

phlpp = dir("./tcgarnabam")[grep("phlpp",dir())]
bamfile = dir("./tcgarnabam")[grep("bam",dir("./tcgarnabam"))]
bamfullnam = file.path("./tcgarnabam",bamfile)
phlpp = intersect(phlpp,bamfile)
bai = dir()[grep("bai",dir())]
phlpp = setdiff(phlpp,bai)
phlpp = phlpp[-1]
phlpp2snp = data.frame(file = phlpp, sample = substr(phlpp,6,15),genotype= "x",freq = "0")
all(phlpp2snp[,1] == phlpp)
for (i in 1:length(phlpp)){
	require(Rsamtools)
	bam = scanBam(phlpp[i])
	y = do.call(DataFrame,bam[[1]])
	y1 = y[which(y$pos >= (loc-100) & y$pos <= (loc)),]

	n = loc-y1[1,"pos"]
	snp = substr(y1[1,"seq"],n+1,n+1)
	if(nrow(y1) >1){
		for (j in 2:nrow(y1)){
		    n = loc-y1[j,"pos"]
		    snp1 = substr(y1[j,"seq"],n+1,n+1)
		    snp = c(snp,snp1)
	    }
	}

	phlpp2snp[i,"freq"] = paste(as.vector(table(snp)),collapse = "_")
	phlpp2snp[i,"genotype"] = paste(names(table(snp)),collapse = "_")
}

View(phlpp2snp)
snpgfile = phlpp2snp[phlpp2snp$genotype == "G",1]
snpafile = phlpp2snp[phlpp2snp$genotype == "A",1]
snpagfile = phlpp2snp[phlpp2snp$genotype == "A_G",1]
x = do.call(snpgfile,function(x)paste(x,".bai",sep=""))
snpgfile = c(paste(snpgfile,".bai",sep = ""),snpgfile)
snpafile = c(paste(snpafile,".bai",sep = ""),snpafile)
snpagfile = c(paste(snpagfile,".bai",sep = ""),snpagfile)

getwd()
dir.create("snpg")
file.copy(snpgfile,"snpg")
dir.create("snpa")
file.copy(snpafile,"snpa")
dir.create("snpag")
file.copy(snpagfile,"snpag")
system("samtools merge snpg.bam ./snpg/*.bam")
system("samtools merge snpa.bam ./snpa/*.bam")
system("samtools merge snpag.bam ./snpag/*.bam")

Sys.setenv(RTEST = "testit", "A+C" = 123)
Sys.getenv("R_TEST")
Sys.getenv()

write.table(snpgfile, file = "snpgfile", quote = F,row.names = F)

system("samtools merge snpgbam.bam ${snpgfile}")
x = as.vector(table(snp))
paste(as.vector(table(snp)),collapse = "_")
table(snp)
showMethods(table)
i = 4

x = as.data.frame.matrix(table(snp))
table(snp)
names(x)
View(y)
View(y1)
vcffile = dir()[grep("vcf",dir())]
vcf = readVcf(vcffile[1],"hg.19")




gr <-
  GRanges(seqnames =
          Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
          ranges =
          IRanges(1:10, width = 10:1, names = head(letters,10)),
          strand =
          Rle(strand(c("-", "+", "*", "+", "-")),
              c(1, 2, 2, 3, 2)),
          score = 1:10,
          GC = seq(1, 0, length=10))
library(XML)

gbmrnaseq = xmlParse("gbmrnaseq.xml")
gbmseqxml = xmlToList(gbmrnaseq)
str(gbmseqxml)
dim(gbmseqxml)

x = gbmseqxml
x[[1]] = NULL
x[[1]] = NULL
x[[520]] = NULL
x[[520]] = NULL
library(plyr)
x1 = ldply(x,unlist())
unlist(x[[1]])
str(x)
cgquery "disease_abbr=GBM&library_strategy=RNA-Seq&study=phs000178&state=live&refassem_short_name=HG19" -o gbmrnaseqhg19.xml
gtdownload -vv -c /usr/bin/cghub.key -d gbmrnaseqhg19.xml -p /media/liu/HD-LCU3/tcgaGBMRnaSeq
cgquery "disease_abbr=GBM&library_strategy=RNA-Seq&study=phs000178&state=live&refassem_short_name=HG19" -i -c cghub.key -v -p /media/liu/HD-LCU3/tcga

gtdownload -vv -d 89324e86-5b7a-4f69-92c1-3b67293f8748 -c https://cghub.ucsc.edu/software/downloads/cghub_public.key

gtdownload -vv -d 88109f88-2e2b-4872-b872-e7e3a8d79b24 -c ~/Downloads/cghub.key

gtdownload -vv -d 89324e86-5b7a-4f69-92c1-3b67293f8748 -c ~/Downloads/cghub.key


cgquery -o t1.xml "state=live&start=0&rows=100"
