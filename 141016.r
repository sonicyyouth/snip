setwd("~/Dropbox/Rworkspace/projects/phlpp2snp")
rm(list = ls())
load("snpexp.rda")
ls()
View(tmx)
View(snpresult)
dim(tmx)
dim(snpresult)
for (i in 1:nrow(snpresult)){
    n = which(rownames(tmx)==snpresult[i,2])
    gene = tmx[which(rownames(tmx) == snpresult[i,3]),which(tmx[n,] == "0" | tmx[n,] == "2")]
    x = tmx[n,which(tmx[n,] == "0" | tmx[n,] == "2")]
    w = wilcox.test(as.numeric(as.character(gene))~factor(as.character(x)))
    t = t.test(as.numeric(as.character(gene))~factor(as.character(x)))
    snpresult[i,4] = length(tmx[n,which(tmx[n,] == "0")])
    snpresult[i,5] = length(tmx[n,which(tmx[n,] == "2")])
    snpresult[i,7] = w$p.value
    snpresult[i,6] = t$estimate[1]-t$estimate[2]
    snpresult[i,8] = mean(as.numeric(tmx[which(tmx[,1] == snpresult[i,3]),2:ncol(tmx)]),na.rm = T)
#   snpresult[i,9] = snpresult[i,6]/snpresult[i,8]
    if (i %% 1000 == 0)writeLines(paste("Calculating ", i,",",round(i/nrow(snpresult)*100, digits = 2), "% done.", sep = ""))
}
save(snpresult,file = "snp6bldresult.rda")
View(snpresult)
snpresult = snpresult[order(snpresult$p.value,decreasing = F),]