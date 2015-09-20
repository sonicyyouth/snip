library(org.Hs.eg.db)                      # to convert from Entrez Gene ID to Symbol
library(biomaRt)
snpannot <- function (input) {
  library(VariantAnnotation)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene) # for annotation
  colnames(input) = c("rsid","chr","pos")
  input$pos = as.numeric(as.character(input$pos))
  target <- with(input, GRanges( seqnames = Rle(chr), ranges = IRanges(pos, end=pos, names=rsid), strand   = Rle(strand("*")) ) )
  loc <- locateVariants(target, TxDb.Hsapiens.UCSC.hg19.knownGene, AllVariants())
  names(loc) <- NULL
  out <- as.data.frame(loc)
  out$names <- names(target)[ out$QUERYID ]
  out <- out[ , c("names",colnames(out)[1:12])]
  out <- out[,-c(8:10)]
  out <- unique(out)
}

probe = read.table(file = "hudsonalpha.org_TCGA_Illumina_HumanHap550K.adf", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE, check.names = FALSE)
x = probe[,9:10]
pos = strsplit(x[,2],split = ",")
x0 = matrix(unlist(pos), ncol=3, byrow=TRUE)
x2 = data.frame(x,matrix(unlist(strsplit(x0[,1],split = ":")),ncol = 2, byrow = T),x0[,2:3])
input = x2[,c(1,3,4)]
rsid = input[which(substr(input[,1],1,2) == "rs"),1]
gr <- rsidsToGRanges(rsid)

hap550annot1 = snpannot(input)
save(hap550annot, file = "hap550annot.rda")

con <- pd.genomewidesnp.6@getdb()
dbListTables(con)
snp6anno <- dbGetQuery(con, "select * from featureSet")
""
input = snpanno[,4:6]
sum(is.na(input[,2]))
sum(is.na(input[,3]))
na = input[which(is.na(input[,2]) | is.na(input[,3])),]

input = input[which(!is.na(input[,2]) & !is.na(input[,3])),]

snpmart = useMart("snp", dataset="hsapiens_snp")
input1 = getBM(c("refsnp_id","chr_name","chrom_start"), filters = "snp_filter",values = rsid,mart = snpmart)

input1 = input1[order(input1[,2]),]
View(input1)
input1 = input1[-(1:116),]
input1 = input1[order(input1[,2],decreasing = T),]
View(input1)
input1 = input1[-(452:611),]
View(input1)
input1 = input1[order(input1[,2]),]
View(input1)

colnames(input1) = colnames(input)
input = rbind(input,input1)
input[,2] = paste("chr",input[,2],sep = "")
snp6ann = snpannot(input)
save(snp6annot, file = "snp6annot.rda")

snp6antcga = read.table(file = "TCGA_Genome_Wide_SNP_6.adf", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE, check.names = FALSE)

snp6anaffy = read.csv(file = "GenomeWideSNP_6.na34.annot.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE, check.names = FALSE,comment.char = "#")

save(snp6annot, snp6anaffy,snp6antcga,file = "snp6annot.rda")

#allannot = snp6ann
allannot = hap550annot1
table(allannot[,7])
snp = allannot[which(allannot[,7] == "threeUTR"),]

snpname = unique(snp[,1])

out = snp

library(org.Hs.eg.db)
Symbol2id <- as.list( org.Hs.egSYMBOL2EG )
id2Symbol <- rep( names(Symbol2id), sapply(Symbol2id, length) )
names(id2Symbol) <- unlist(Symbol2id)

x <- unique(c(out$GENEID, out$PRECEDEID, out$FOLLOWID))
table( x %in% names(id2Symbol) ) # good, all found

out$GENESYMBOL <- id2Symbol[ as.character(out$GENEID) ]
out$PRECEDESYMBOL <- id2Symbol[ as.character(out$PRECEDEID) ]
out$FOLLOWSYMBOL <- id2Symbol[ as.character(out$FOLLOWID) ]
snp = out

hap5503utr = snp
#snp63utr = snp
save(snp63utr, hap5503utr, file = "snp3utrlist.rda")

# ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")



# library(NCBI2R)
# library(plyr)
# scott = AnnotateSNPList(na[,1])
# scottsobject = AnnotateSNPList(na[1:100,1])
# s = scottsobject
# t = nrow(na)%/%100-1
# for (i in 6:t){
#   n =1+ 100*i
#   scot = AnnotateSNPList(na[n:(99+n),1])
#   s = rbind.fill(s,scot)
#  }
# s = s[order(s[,3]),]
# s0 = s[735:1123,]
# mis = setdiff(na[,1],s0[,1])




# hap550annot = s
# save(hap550annot, file = "hap550annot.rda")
