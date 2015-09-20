#library(SNPlocs.Hsapiens.dbSNP.20120608)
#x = read.table(file = "/media/liu/My Passport/gbmhap550/hudsonalpha.org_GBM.HumanHap550.3.5.0.Genotypes.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE, check.names = FALSE)
#x = x[-1,]
#haploc =x[,1:3]
#save(haploc,file = "haphg19info.rda")
#load("/home/liu/annotation_data/snp6annotsimple.rda")
#lsxy = intersect(snp6id,snp6infoid)
#unsnp = setdiff(snp6id,snp6infoid)
#length(xy)
#length(unsnp)
#library(biomaRt)
#snpmart = useMart("snp", dataset="hsapiens_snp")
#unsnpid = getBM(c('refsnp_id','allele','chr_name','chrom_start','chrom_strand'), filters = 'snp_filter', values = unsnp, mart = snpmart)
#
#listFilters(snpmart)
#snp6info = getBM(c('refsnp_id','allele','chr_name','chrom_start','chrom_strand'), filters = 'snp_filter', values = snpid, mart = snpmart)
#
#head(pakthapresult)
#snp6id = rownames(pakthapresult)
#snp6infoid = snp6annot[,2]
#length(intersect(snp6id,snp6infoid))
#length(snp6id)
#snpid[1:10]
#snp6gr <- rsidsToGRanges(snpid[1:100000])
getwd()
setwd("~/Dropbox/Rworkspace/projects/phlpp2snp")
load("/home/liu/annotation_data/hap550annot.rda")
load("paktsnp6linearresult.rda")
paktsnp6result = pakthapresult
load("pakthapresult.rda")
snpname = unique(c(rownames(pakthapresult),rownames(paktsnp6result)))

length(snpname)

library(biomaRt)

snpmart = useMart("snp", dataset="hsapiens_snp")
unsnphap = setdiff(rownames(pakthapresult),hapsnpinfo$marker)
table(substr(rownames(pakthapresult),1,2))
pakthapresult = pakthapresult[which(substr(rownames(pakthapresult),1,2) == "rs"),]
colnames(hapsnpinfo)[1] = colnames(unsnp6id)[1]
colnames(hapsnpinfo)[4] = colnames(unsnp6id)[3]
colnames(hapsnpinfo)[5] = colnames(unsnp6id)[4]
hapsnpinfo$chr_name = paste("chr",hapsnpinfo$chr_name,sep = "")
snphapgrange = hapsnpinfo[,c(1,4,5)]
snphapgrange$chr_name = paste("chr",snphapgrange$chr_name,sep = "")
rownames(snphapgrange) = snphapgrange$refsnp_id
unhap = setdiff(rownames(pakthapresult),snphapgrange$refsnp_id)

unsnp6 = setdiff(rownames(paktsnp6result),snp6annor$names)
unsnp6id = getBM(c('refsnp_id','allele','chr_name','chrom_start','chrom_strand'), filters = 'snp_filter', values = unsnp6, mart = snpmart)
unsnp6id$chr_name = paste("chr",unsnp6id$chr_name,sep = "")
colnames(snp6annor)[1] = colnames(unsnp6id)[1]
colnames(snp6annor)[2] = colnames(unsnp6id)[3]
colnames(snp6annor)[3] = colnames(unsnp6id)[4]
snp6grange = rbind(snp6annor[1:3],unsnp6id[,c(1,3,4)])
snp6grange = snp6grange[!duplicated(snp6grange$refsnp_id),]
rownames(snp6grange) = snp6grange$refsnp_id
save(snphapgrange,file = "hapsnpgrange.rda")
save(snp6grange,file = "snp6grange.rda")
load("snp6grange.rda")
load("hapsnpgrange.rda")
gbmsnphap = data.frame(snphapgrange[rownames(pakthapresult),],pakthapresult)
gbmsnp6pakt = data.frame(snp6grange[rownames(paktsnp6result),],paktsnp6result)
dim(snp6grange)
snp6grange[97000,]

for (i in 1:105){
	snpid1 = snpid = getBM(c('refsnp_id','allele','chr_name','chrom_start','chrom_strand'), filters = 'snp_filter', values = snpname[(10000*i+1),(10000*(i+1))], mart = snpmart)
	snpid = rbind(snpid,snpid1)

}
snpid = getBM(c('refsnp_id','allele','chr_name','chrom_start','chrom_strand'), filters = 'snp_filter', values = snpname, mart = snpmart)
#unsnp = setdiff(rownames(paktsnp6result),snp6annor$names)

load("snp6ch19info.rda")

load("hap550snpinfo.rda")

x = data.frame(hapsnpinfo[match(rownames(paktsnp6result),snp6annor$names),c(1:3,6)],paktsnp6result)
length(unsnp)
View(unsnpid)
x = data.frame(hapsnpinfo[match(rownames(pakthapresult)[1:500],hapsnpinfo$marker),c(1,2,4:6)],pakthapresult)

library(BSgenome)
genome = getBSgenome("GRCh38")
x = seqinfo(genome)
snpGrangechr <- function (input) {
  require(GenomicRanges)
  require(BSgenome)
  genome = getBSgenome("GRCh38")
  seqname = paste("chr",c(as.character(seq(1,22)),"X","Y"),sep = "")
  seqinfo = seqinfo(genome)[c(as.character(seq(1,22)),"X","Y")]
  seqlevels(seqinfo) = paste("chr",seqlevels(seqinfo),sep = "")
  #library(TxDb.Hsapiens.UCSC.hg19.knownGene) # for annotation
  colnames(input)[1:3] = c("rsid","chr","pos")
  #input$chr = sub("chr","",input$chr)
  input$pos = as.numeric(as.character(input$pos))
  input = input[input$chr %in% seqname,]
  target <- with(input, GRanges( seqnames = Rle(chr), ranges = IRanges(pos, end=pos, names=rsid), strand   = Rle(strand("*")), seqinfo = seqinfo) )
  mcols(target) <- input[,c(1,4:ncol(input))]
  #target = keepSeqlevels(target, seqname)
  return(target)
}

gbmhapgr = snpGrange(gbmsnphap)
gbmsnp6pakt[,2] = paste("chr",gbmsnp6pakt[,2],sep = "")
gbmsnp6gr = snpGrange(gbmsnp6pakt)
mcols(gbmhapgr)$platform = "hapman550"
mcols(gbmsnp6gr)$platform = "affysnp6"
gbmhapgr[1:10,]
mcols(gbmhapgr)$logp = - log10(mcols(gbmhapgr)$pvalue)
mcols(gbmsnp6gr)$logp = - log10(mcols(gbmsnp6gr)$pvalue)
mcols(gbmhapgr) = mcols(gbmhapgr[,-6])

save(gbmhapgr,file = "gbmhapgr.rda")
save(gbmsnp6gr, file = "gbmsnp6gr.rda")
gbmsnpcom = c(gbmhapgr,gbmsnp6gr)
save(gbmsnpcom, file = "gbmsnpcom.rda")
load("gbmhapgr.rda")
load("gbmsnp6gr.rda")
load("gbmsnpcom.rda")
library("ggbio")


plotGrandLinear(gbmhapgr, geom = "point", coord = "genome", aes(y = logp), cutoff = 4, cutoff.color = "blue", cutoff.size=0.2,xlab = "SNP in hapman550 platform")
plotGrandLinear(gbmsnp6gr, geom = "point", coord = "genome", aes(y = logp), cutoff = 4, cutoff.color = "blue", cutoff.size=0.2,xlab = "SNP in Affy snp6 platform")
plotGrandLinear(gbmsnpcom, geom = "point", coord = "genome", aes(y = logp), cutoff = 4, cutoff.color = "blue", cutoff.size=0.2,xlab = "SNP in both platforms")

autoplot(x16, geom = "point", coord = "genome", aes(y = logp), cutoff = 3, cutoff.color = "blue", cutoff.size=0.2)
l1 = 71300000
l2 = l1 + 2e6
gro <- GRanges(c("chr16"), IRanges(l1, width = 8e5))
x17 = gbmsnpcom[seqnames(gbmsnpcom) == "chr17"]
x16 = gbmsnpcom[seqnames(gbmsnpcom) == "chr16"]
#seqlevels(x17) = paste("chr",seqlevels(x17),sep = "")
#seqlevels(x16) = paste("chr",seqlevels(x16),sep = "")
plotGrandLinear(x16, aes(y = logp), highlight.gr = gro, cutoff = 4, cutoff.color = "blue", cutoff.size=0.2)
x = keepSeqlevels(gbmhapgr, seqname[1:22])
mcols(x2) = mcols(x2)$logp
library(Gviz)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#columns(txdb)
#bmTrack <- BiomartGeneRegionTrack(start=26682683, end=26711643,chromosome=7, genome="GRCh38")
txgr = transcripts(txdb,columns = c("tx_id","tx_name","gene_id"))
mcols(txgr)$gene_id = as.character(mcols(txgr)$gene_id)
txgr = txgr[!is.na(mcols(txgr)$gene_id)]
library(org.Hs.eg.db)
x = select(org.Hs.eg.db, keys=as.character(mcols(txgr)$gene_id), columns="SYMBOL", keytype="ENTREZID")
id = as.character(mcols(txgr)$gene_id)
sym = x[match(id,x[,1]),2]
mcols(txgr)$symbol = x[match(id,x[,1]),2]
mcols(txgr)$gene = mcols(txgr)$gene_id

l1 = 71300000
gro <- GRanges(c("chr16"), IRanges(l1, width = 2e6))
txgrchr = subsetByOverlaps(txgr,gro)
xgro = subsetByOverlaps(x16,gro)
mcols(xgro) = mcols(xgro)$logp
grtrack = GeneRegionTrack(txgrchr,genome = "hg19", chromosome = "chr16",name = "Gene Model", transcriptAnnotation="symbol",background.title="brown",collapseTranscripts = T)
gtrack <- GenomeAxisTrack(littleTicks = TRUE)
itrack <- IdeogramTrack(genome="hg38", chromosome="chr16",showId = F)
atrack <- AnnotationTrack(xgro, name="SNP")
dtrack = DataTrack(xgro,name = "log(1/pvalue)")
plotTracks(list(itrack,gtrack,grtrack,dtrack))
#bmTrack <- BiomartGeneRegionTrack(start=26682683, end=26711643,chromosome=7, genome="hg19")

grtrack <- GeneRegionTrack(geneModels, name = "foo")

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb19 <- TxDb.Hsapiens.UCSC.hg19.knownGene
#columns(txdb)
#bmTrack <- BiomartGeneRegionTrack(start=26682683, end=26711643,chromosome=7, genome="GRCh38")
txgr = exons(txdb19,columns = c("tx_id","tx_name","gene_id"))
mcols(txgr)$gene_id = unlist(mcols(txgr)$gene_id)
x = mcols(txgr)
txgr = txgr[unlist(lapply(mcols(txgr)$gene_id,function(x) !all(is.na(x))))]
library(org.Hs.eg.db)
ref = select(org.Hs.eg.db, keys=unlist(mcols(txgr)$gene_id), columns="SYMBOL", keytype="ENTREZID")
#id = as.character(mcols(txgr)$gene_id)
sym = x[match(id,x[,1]),2]
mcols(txgr)$symbol = lapply(mcols(txgr)$gene_id,function(x)ref[match(x,ref[,1]),2])
mcols(txgr)$gene = mcols(txgr)$gene_id
l1 = 71678852 - 1000
l2 = l1 + 100000
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
altrack = AlignmentsTrack("phlpp2bam.bam",isPaired = F)
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

plotTracks(list(itrack,gtrack,refgene,alatrack),from = lst,to = led-2000,chromosome = "chr16")
A <- 2
f <- function(x) print(x^2)
env <- new.env()
assign("A", 10, envir = env)
assign("f", f, envir = env)
f <- function(x) print(x)
f(A)                                      # 2
do.call("f", list(A))                     # 2
do.call("f", list(A), envir = env)        # 4
do.call(f, list(A), envir = env)          # 2
do.call("f", list(quote(A)), envir = env) # 100
do.call(f, list(quote(A)), envir = env)   # 10
do.call("f", list(as.name("A")), envir = env) # 100

gc(reset = T)