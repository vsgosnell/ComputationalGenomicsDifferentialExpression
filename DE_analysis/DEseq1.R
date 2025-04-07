#remove all variables from environment
rm(list = ls())

#reset plot margins -> so 4 plots don't display at same time
par(mfrow=c(1,1))

#set wd
setwd('/Users/veronicagosnell/Desktop/R/ComputationalGenomicsDifferentialExpression/DE_analysis')

#load the data
df <- read.table('vehicle_drug_feature_counts.txt',
                 header = T, sep = '\t', row.names = 1)

#check the df
head(df)

#take only the gene read counts from the file -> columns 6-9
df <- df[,6:9]

#check again
head(df)


#check the distribution of counts per sample
#set parameters for 4 pane histogram plots
par(mfrow = c(2,2), mar = c(4,4,1,1))
print( apply(df, 2, function(x) quantile(as.numeric(x))) )


#create a 4 pane histogram plot
hist.plots <- apply(df, 2, function(x) {hist(log2(x), breaks = 100)})

#if any number is bigger than 20 in read counts, it is kept
abline(v = log2(20))

expID <- apply(df, 1, function(x){any(x > 20)})

dfexpID <- df[expID, ]

apply(dfexpID, 2, function(x){hist(log2(x),
                                   breaks = 100)})

#initial filter
#Low-expressed genes can make for noisy data. We'll filter out the lowest-expressed genes (bottom 25%), which corresponds to less than ~20 (or 2^ ~4.3) reads.
#Keep genes with expression bigger than 20
expressed.ids <- apply(df, 1, function(x) any(x > 20))
dfExp <- df[expressed.ids, ]


#Check distribution again after filtering
par(mfrow = c(2,2), mar = c(4,4,1,1))
print( apply(dfExp, 2, function(x) quantile(as.numeric(x))) )

#plot the new histograms
hist.plots <- apply(dfExp, 2, function(x) {hist(log2(x), breaks = 100)})


#We have filtered data ready for normalization and expression-analysis with DESeq2.



#install bioconductor for DEseq2
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomeInfoDbData", force = TRUE)

BiocManager::install("DESeq2")

library(DESeq2)


#Normalized reads count
norCounts <- counts(dds, normalized=TRUE)

#set wd
setwd('/Users/veronicagosnell/Desktop/R/ComputationalGenomicsDifferentialExpression/DE_analysis')

res <- read.table('/Users/veronicagosnell/Desktop/R/ComputationalGenomicsDifferentialExpression/deseq.txt')

head(res)

# extract the genes with adjusted P-value < 0.01
resSig <- res[ which(res$padj < 0.01), ]
dim(resSig)

# sort by log2FoldChange to get the significant genes with the strongest down-regulation
head( resSig[ order( resSig$log2FoldChange ), ] )
#order it by the log2fold change, and rank the matrix



# and strongest up-regulation
tail( resSig[ order( resSig$log2FoldChange ), ] )



library("RColorBrewer")
install.packages("gplots")
library("gplots")


#want to see gene expression between control and treatment
#have to use normalized counts
#we only want the significant genes to be represented
#plot heatmap of normalized read-counts

sigNorData <- norCounts[rownames(norCounts) %in% rownames(resSig),]
hmcol <-  colorRampPalette(brewer.pal(9, "GnBu"))(100)
## expression heatmap
heatMapDF_nor <- t( apply(sigNorData, 1, function(x){(x-mean(x))/sd(x)}) )

colnames(heatMapDF_nor) <- c('control1','control2',
                             'treat1'  , 'treat2')

## comment out lines starting with "pdf" & "dev.off" to view your plot in your "plots" window, or leave as-is to save your plot to .pdf
#pdf('heatmap.pdf', height = 10, width = 10)
heatmap.2(heatMapDF_nor, col = hmcol, trace="none", margin=c(10, 10),labRow = F)


