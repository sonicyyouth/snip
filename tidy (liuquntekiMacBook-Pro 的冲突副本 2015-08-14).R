getwd()

setwd("~/Dropbox/Rworkspace/projects/snp")
setwd("~/Downloads/snp")
x = dir()[1]
snp1 = read.table(file = x, head = F, skip =2,stringsAsFactors = F, quote = "")
dim(snp1)

head(snp1)
