setwd("~/Dropbox/Rworkspace/projects/phlpp2snp")
library(Rsamtools)
bamfile = dir("./tcgarnabam")[grep("bam",dir("./tcgarnabam"))]
file = read.csv("~/Dropbox/Rworkspace/gbmrnaseqhg19.csv")

filename = bamfile
renam = gsub("phlpp","",filename)
setwd("~/Dropbox/Rworkspace/projects/phlpp2snp/gbmrnabam")
file.rename(filename,renam)
filename = dir()[grep("bam",dir())]
file$uucid =sapply(file$filename,function(x)strsplit(x,"\\.")[[1]][1])
uucidname = sapply(filename,function(x)strsplit(x,"\\.")[[1]][1])
sample = file[match(uucidname,file$uucid),"legacy_sample_id"]
samplname = paste(substr(sample,1,15),".bam",sep = "")
file.rename(filename,samplname)
samplname =samplname[-113]
indexBam(samplname)
setwd("~/Dropbox/Rworkspace/projects/phlpp2snp")
bamfullnam = file.path("./tcgarnabam",samplname)
load("gbmallhap.rda")
seqsampl = sub(".bam","",samplname)
seqsampl[1:10]
dim(gbmallhap)
load("gbmhapsnp.rda")
colnames(gbmhapsnp) = substr(colnames(gbmhapsnp),1,15)
xy = intersect(colnames(gbmhapsnp), seqsampl)
ls()
library(Gviz)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb19 <- TxDb.Hsapiens.UCSC.hg19.knownGene
#columns(txdb)
#bmTrack <- BiomartGeneRegionTrack(start=26682683, end=26711643,chromosome=7, genome="GRCh38")
txgr = exons(txdb19,columns = c("tx_id","tx_name","gene_id"))
mcols(txgr)$gene_id = unlist(mcols(txgr)$gene_id)
txgr = txgr[unlist(lapply(mcols(txgr)$gene_id,function(x) !all(is.na(x))))]
library(org.Hs.eg.db)
ref = select(org.Hs.eg.db, keys=unlist(mcols(txgr)$gene_id), columns="SYMBOL", keytype="ENTREZID")
#id = as.character(mcols(txgr)$gene_id)
mcols(txgr)$symbol = lapply(mcols(txgr)$gene_id,function(x)ref[match(x,ref[,1]),2])
mcols(txgr)$gene = mcols(txgr)$gene_id
txgr1 = txgr
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
bam = scanBam("phlpp2bam.bam")
y = do.call(cbind.data.frame,bam[[1]])
bam1 = scanBam("phlppSRR1521352.bam")
sapply(bam[[1]],class)
indexBam("snpa.bam")
indexBam("snpg.bam")
indexBam("snpag.bam")
altrack = AlignmentsTrack(bamfullnam[1],isPaired = F)
alatrack = AlignmentsTrack("snpa.bam",isPaired = F)
algtrack = AlignmentsTrack("snpg.bam",isPaired = F)
alagtrack = AlignmentsTrack("snpag.bam",isPaired = F)
#atrack <- AnnotationTrack(xgro, name="SNP")
#dtrack = DataTrack(xgro,name = "log(1/pvalue)")
pdf("snpgenotypeseq.pdf")
plotTracks(list(itrack,gtrack,grtrack,altrack))
plotTracks(list(itrack,gtrack,refgene,altrack,alatrack,algtrack,alagtrack),from = lst,to = led-2000,chromosome = "chr16")
dev.off()
bmt

plotTracks(list(itrack,gtrack,refgene,altrack),from = lst,to = lst+3000,chromosome = "chr16")

