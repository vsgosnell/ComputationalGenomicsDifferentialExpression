install.packages("gplots")

library("RColorBrewer")
library("gplots")


#plot heatmap of normalized read-counts

sigNorData <- norCounts[rownames(norCounts) %in% rownames(resSig),]
hmcol <-  colorRampPalette(brewer.pal(9, "GnBu"))(100)
## expression heatmap
heatMapDF_nor <- t( apply(sigNorData, 1, function(x){(x-mean(x))/sd(x)}) )

colnames(heatMapDF_nor) <- c('control1','control2',
                             'treat1'  , 'treat2')


## comment out lines starting with "pdf" & "dev.off" to view your plot in your "plots" window, or leave as-is to save your plot to .pdf

#save the plot
pdf("heatmap.pdf", height = 10, width = 10)

#create the plot 
heatmap.2(heatMapDF_nor, col = hmcol, trace = "none", margin = c(10, 10), labRow = FALSE)

#make sure to run this in order to save the pdf properly
dev.off()
