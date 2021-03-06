source("http://bioconductor.org/biocLite.R")
biocLite("GGtools")
biocLite("GGdata")
install.packages("MatrixEQTL")
library(MatrixEQTL)
library(GGdata)
library(GGtools)

cc = new("CisConfig")
chrnames(cc) = "16"
estimates(cc) = FALSE
f1 <- All.cis(cc) 
length(f1)
grep("rs1057147",x$snp)
View(x)
x=mcols(f1)
base.dir = find.package('MatrixEQTL')
base.dir
useModel = "modelLINEAR"
SNP_file_name = paste(base.dir, "/data/SNP.txt", sep="");
expression_file_name = paste(base.dir, "/data/GE.txt", sep="");
covariates_file_name = paste(base.dir, "/data/Covariates.txt", sep="");
output_file_name = tempfile();
pvOutputThreshold = 1e-2;
errorCovariance = numeric();
SNP_file_name
snps = SlicedData$new();
snps$fileDelimiter = "\t"; 
## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
cvrt$LoadFile(covariates_file_name);
}

## Run the analysis

me = Matrix_eQTL_engine(
snps = snps,
gene = gene,
cvrt = cvrt,
output_file_name = output_file_name,
pvOutputThreshold = pvOutputThreshold,
useModel = useModel,
errorCovariance = errorCovariance,
verbose = TRUE,
pvalue.hist = TRUE,
min.pv.by.genesnp = FALSE,
noFDRsaveMemory = FALSE);

unlink(output_file_name);

## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected eQTLs:', '\n');
show(me$all$eqtls)

## Plot the histogram of all p-values

plot(me)

suppressPackageStartupMessages(library(GGtools))
library(parallel)
g16 = getSS("GGdata", "16")
t1 = gwSnpTests(genesym("PHLPP2")~male, g16)
plot(t1, snplocsDefault())
x = topSnps(t1, n=109692)[[1]][1]
grep("rs1050560",rownames(x))
head(x)
x[1:100,]
smList(g16)
as(smList(g16)[[1]][1:5,1:5], "matrix")
gwSnpTests
ls()
rm(list =ls())

