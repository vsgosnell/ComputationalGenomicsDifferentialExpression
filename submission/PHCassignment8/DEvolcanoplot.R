# generate a data frame with a column assigning color to significant differentially expressed genes
res_plot      <- data.frame( res )
res_plot$col  <- 'gray40'

# setting cutoffs for significantly up (red) and down (blue) genes; everything else is grey
res_plot$col[res_plot$log2FoldChange > 1 & res_plot$padj < 0.01] <- 'red'
res_plot$col[res_plot$log2FoldChange < -1 & res_plot$padj < 0.01] <- 'cornflowerblue'


## comment out lines starting with "pdf" & "dev.off" to view your plot in your "plots" window, or leave as-is to save your plot to .pdf

#save the plot
pdf('volcano.pdf', height = 5, width = 5)

#create margins for plot
par(mar = c(5,5,1,1))

#create plot
plot(res_plot$log2FoldChange,
     -log10(res_plot$padj),
     col = res_plot$col, pch = 19, xlab = 'log2(fold change)',
     ylab = '-log10(p-adj)', xlim = c(-8,8),
     
)


#make sure to run this so the pdf saves correctly
dev.off()

