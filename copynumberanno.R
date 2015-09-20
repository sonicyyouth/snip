glmcna = read.table(file = "../tcgafile/GBMLGG.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_hg19__seg.seg.txt", head = T, sep = "\t")
glmcna = glmcna[which(glmcna$Chromosome < 23),]
glmcna$Chromosome = paste("chr",glmcna$Chromosome,sep = "")
library(GenomicRanges)
library(dplyr)
glmcn <- GRanges(seqnames =Rle(glmcna$Chromosome),ranges =IRanges(glmcna$Start, end = glmcna$End, names =glmcna$Sample),strand =Rle(strand(c("*")),),sample = glmcna$Sample,segmean = glmcna$Segment_Mean)

table(glmcna$Chromosome)

chrinfo <- read.table(file = "../tcgafile/chromInfo.txt", head = F)
chrinfo = chrinfo %>% filter(grepl("^chr[0-9]+$", V1)) %>%select(V1,V2) %>% arrange(V1)
View(chrinfo)
maxl = glmcna %>% group_by(Chromosome) %>% summarise(max = max(End))
chrinfo$V2 > maxl$max

seqlengths(glmcn) <- chrinfo[match(names(seqlengths(glmcn)),chrinfo$V1),"V2"]

library(tidyr)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
hg19genes = genes(txdb)

xin = findOverlaps(hg19genes, glmcn,type = "within")
x0 = findOverlaps(hg19genes, glmcn)
xin = as.data.frame(xin)
x0 = as.data.frame(x0)
x = glmcn[xin$subjectHits]
x1 = hg19genes[xin$queryHits]
x01 = glmcn[x0$subjectHits]
x02 = hg19genes[x0$queryHits]
sample = x$sample
segmean = x$segmean
gene = x1$gene_id
glmgenecn = data.frame(sample = sample,gene = gene, segmean = segmean)
glmgenecnall = data.frame(sample = x01$sample,gene = x02$gene_id, segmean = x01$segmean)
glmgenecn$sample = substr(glmgenecn$sample,1,15)
glmgenecna = spread(glmgenecn,gene,segmean)
glmgenecnaall= spread(glmgenecnall,gene,segmean)
dim(glmgenecna)
genesynb = AnnotationDbi::select(org.Hs.eg.db, keys=unique(as.character(glmgenecn$gene)), columns=c("SYMBOL","ENTREZID"), keytype="ENTREZID")
View(glmgenecna)
save(glmgenecn,file = "glioma_gene_cna.rda")
save(glmgenecnall,file = "glioma_gene_cnaall.rda")
p16 = glmgenecnall %>% filter(gene == "1029")
sum(is.na(p16$segmean))
hist(p16$segmean)
p16 = p16[which(substr(p16$sample,14,15) == "01" |substr(p16$sample,14,15) == "02"),]
sum(p16$segmean < -1)
p16in = glmgenecn %>% filter(gene == "1029")
p16in = p16in[which(substr(p16in$sample,14,15) == "01" |substr(p16in$sample,14,15) == "02"),]
sum(is.na(p16in$segmean))
hist(p16$segmean)
p16 = p16[which(substr(p16$sample,14,15) == "01" |substr(p16$sample,14,15) == "02"),]



