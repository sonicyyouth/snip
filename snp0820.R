load("~/cloud360/R workspace/GBM_microarray_expression.rda")
load("~/cloud360/R workspace/tcga/ProcessedData.Firehose/GBM_mirna_array.rda")
gbm4502a = gbm4502a[,unique(colnames(gbm4502a))]
colnames(gbm4502a) = substr(colnames(gbm4502a),1,12)

mir = gbm_mirna[,-1]
x = colnames(mir)
table(substr(x,14,15))
x = x[which(substr(x,14,15) == "01")]
mir = mir[,x]
colnames(mir) = substr(colnames(mir),1,12)

x = colnames(gbm4502a)
y = colnames(mir)
xy = intersect(x,y)
genmir = data.frame(t(gbm4502a[,xy]),t(mir[,xy]))

library(org.Hs.eg.db)
cols <- c("PFAM","GO", "SYMBOL")
x = select(org.Hs.eg.db, "GO:0010507", cols, keytype="GO")
rname = unique(x[,5] )
rn = intersect(rownames(gbm4502a),rname)
genmir = data.frame(t(mir[,xy]),t(gbm4502a[rn,xy]))
grep("hcmv.mir.us5.2",colnames(genmir))
hsagenmir = genmir[,-(1:46)]
x = cor(hsagenmir)
n = nrow(mir)
n=n-46
x = x[-(1:n),]
x = t(x[,1:n])
x = data.frame(rowSums(x),x)
x = x[order(x[,1]),]
xall =x
