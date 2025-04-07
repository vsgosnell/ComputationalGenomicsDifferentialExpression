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


table2 <- data.frame(name = colnames(dfExp),
                     condition = c('control','control', 'treatment', 'treatment'))
dds <- DESeqDataSetFromMatrix(dfExp, 
                              colData=table2, design= ~ condition)


dds <- DESeq(dds)


#Normalized reads count
norCounts <- counts(dds, normalized=TRUE)


res <- results(dds)


#Getting differentially expressed genes 
#Take a look at our results object ("res")
res


# res is a dataframe object, so we can check out metadata for what the columns mean
def <- mcols(res, use.names=TRUE)
def@listData[["description"]]


# extract the genes with adjusted P-value < 0.01
resSig <- res[ which(res$padj < 0.01), ]
dim(resSig)


# sort by log2FoldChange to get the significant genes with the strongest down-regulation
head( resSig[ order( resSig$log2FoldChange ), ] )


# and strongest up-regulation
tail( resSig[ order( resSig$log2FoldChange ), ] )


resSig.sorted.df <- as.data.frame(resSig[ order( resSig$log2FoldChange ), ])
# write the full list to a sorted, tab-delimited table for your records. The file will be created in your Routput folder.
# write.table(resSig.sorted.df,
#             file = "Routput/resSig.sorted.tab.txt",
#             sep = "\t",
#             row.names = TRUE,
#             quote = FALSE)


