setwd("~/Dropbox/Rworkspace/projects/phlpp2snp")
library(SNPlocs.Hsapiens.dbSNP.20120608)
gr = rsidsToGRanges("rs10500560")
loc = start(gr)
phlpp = dir("./tcgarnabam")[grep("phlpp",dir())]
bamfile = dir("./tcgarnabam")[grep("bam",dir("./tcgarnabam"))]

bai = dir("./tcgarnabam")[grep("bai",dir("./tcgarnabam"))]
phlpp = setdiff(bamfile,bai)
bamfullnam = file.path("./tcgarnabam",phlpp)

phlpp2snp = data.frame(file = phlpp, sample = substr(phlpp,1,15),genotype= "x",freq = "0")
all(phlpp2snp[,1] == phlpp)
for (i in 1:length(bamfullnam)){
	require(Rsamtools)
	bam = scanBam(bamfullnam[i])
	y = do.call(DataFrame,bam[[1]])
	y1 = y[which(y$pos >= (loc-y$qwidth[1]+1) & y$pos <= (loc)),]
	if(nrow(y1) >1){
		n = loc-y1[1,"pos"]
	    snp = substr(y1[1,"seq"],n+1,n+1)
		for (j in 2:nrow(y1)){
		    n = loc-y1[j,"pos"]
		    snp1 = substr(y1[j,"seq"],n+1,n+1)
		    snp = c(snp,snp1)
	    }
	    phlpp2snp[i,"freq"] = paste(as.vector(table(snp)),collapse = "_")
	    phlpp2snp[i,"genotype"] = paste(names(table(snp)),collapse = "_")
	}
}
View(phlpp2snp)
phlpp2snp = phlpp2snp[order(phlpp2snp$sample),]
i = 79
bamfullnam[145]
phlpp2snp[79,1]
table(phlpp2snp$genotype)
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


load("exonhg19ref.rda")
l1 = 71678852 - 1000
l2 = l1 + 1000000
gro <- GRanges(c("chr16"), IRanges(l1, width = 90000))
txgrchr = subsetByOverlaps(txgr,gro)
#xgro = subsetByOverlaps(x16,gro)
#mcols(xgro) = mcols(xgro)$logp

grtrack = GeneRegionTrack(txgrchr,genome = "hg19", chromosome = "chr16",name = "Gene Model", transcriptAnnotation="symbol",background.title="brown",collapseTranscripts = F)
bmt = BiomartGeneRegionTrack(genome = "hg19", chromosome = "chr16",start = l1, end = l2, filter = list(with_ox_refseq_mrna = T),starking = "dense")
refgene <- UcscTrack(genome = "hg19", chromosome = "chr16", track = "refGene", from = l1, to = l2, trackType = "GeneRegionTrack", rstarts = "exonStarts", rends = "exonEnds", gene = "name", symbol = "name2", transcript = "name", strand = "strand",fill = "#8282d2", name = "Ref Genes",transcriptAnnotation="symbol")
lst = min(start(refgene[which(symbol(refgene) == "PHLPP2")]))
led = max(end(refgene[which(symbol(refgene) == "PHLPP2")]))

gtrack <- GenomeAxisTrack(littleTicks = TRUE)
itrack <- IdeogramTrack(genome="hg19", chromosome="chr16",showId = F)
library(Rsamtools)
bam = scanBam(bamfullnam[145])
y = do.call(cbind.data.frame,bam[[1]])
bam1 = scanBam("phlppSRR1521352.bam")
sapply(bam[[1]],class)
indexBam("snpa.bam")
indexBam("snpg.bam")
indexBam("snpag.bam")
altrack = AlignmentsTrack(bamfullnam[145],isPaired = T)
alatrack = AlignmentsTrack("snpa.bam",isPaired = F)
algtrack = AlignmentsTrack("snpg.bam",isPaired = F)
alagtrack = AlignmentsTrack("snpag.bam",isPaired = F)
#atrack <- AnnotationTrack(xgro, name="SNP")
#dtrack = DataTrack(xgro,name = "log(1/pvalue)")
pdf("snpgenotypeseq.pdf")
plotTracks(list(itrack,gtrack,grtrack,altrack))
plotTracks(list(itrack,gtrack,refgene,alatrack,algtrack,alagtrack),from = lst,to = led-2000,chromosome = "chr16")
dev.off()
bmt

plotTracks(list(itrack,gtrack,refgene,altrack),from = lst,to = lst+3000,chromosome = "chr16")
