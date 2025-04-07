#load function 'go_enrichment' for go enrichment analysis
source("source.R")

# load GO gaf file (GO annotation file)
geneGO <- read.delim2('PlasmoDB-48_Pfalciparum3D7_GO.gaf.txt', header = F, sep = '\t')
geneGO <- geneGO[,c(2,5)]


# load ontology
library(GO.db)
go_db <- Term(GOTERM)
go_On <- Ontology(GOTERM)
go_inf<- data.frame(ontology = go_On,
                    term     = go_db)


# now we'll run CW's go_enrichment function perform a Fisher test to find any GO terms that appear in our differentially expressed genes more often than would be expected by chance (see "RSource/source.R" for details)
goEnrichment <- go_enrichment(rownames(resSig)[resSig$log2FoldChange > 1],
                              rownames(df),
                              geneGO)

idMatch <- match(rownames(goEnrichment), rownames(go_inf))
goEnrichment <- data.frame(goEnrichment,
                           go_inf[idMatch, ])

# List the enriched GO terms
goEnrichment$term[goEnrichment$padj < 0.1]



#bar plot for GO enrichment

goSig   <- goEnrichment[goEnrichment$padj < 0.01,]
goSig$expection <- goSig$X.BackgroundWithGOterm/(goSig$X.BackgroundWithGOterm +
                                                   goSig$X.BackgroundWithoutGOterm) * (goSig$X.QueryWithGOterm + 
                                                                                         goSig$X.QueryWithoutGOterm)

goSig <- goSig[order(goSig$X.QueryWithGOterm), ]
goSigDraw <- t( goSig[,c(1,9)] )
colnames(goSigDraw) <- goSig$term


## comment out lines starting with "pdf" & "dev.off" to view your plot in your "plots" window, or leave as-is to save your plot to .pdf

#save the plot
pdf('GO_barplot.pdf', height = 10, width = 10)

#create margins for plot
par(mar = c(4,20,1,1))

#create the plot
barplot(goSigDraw, horiz = T, las = 1)

#make sure to run this so pdf saves correctly
dev.off()
