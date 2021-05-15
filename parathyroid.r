
BiocManager::install("GenomicFeatures")
BiocManager::install( "DESeq2" )
BiocManager::install( "Rsamtools" )
BiocManager::install("GenomicAlignments")
BiocManager::install("rtracklayer")
BiocManager::install("parathyroidSE")
BiocManager::install("pasillaBamSubset")
BiocManager::install("genefilter")

library( "GenomicFeatures" )
library( "DESeq2" )
library( "Rsamtools" )
library( "GenomicAlignments" )
library( "rtracklayer" )
library( "parathyroidSE" )
library( "pasillaBamSubset" )
library( "gplots" )
library( "RColorBrewer" )
library( "genefilter" )


data( "parathyroidGenesSE" )
se <- parathyroidGenesSE
colnames(se) <- se$run


colnames(se)


countdata <- assay( parathyroidGenesSE )
head( countdata )

coldata <- colData( parathyroidGenesSE )
rownames( coldata ) <- coldata$run
colnames( countdata ) <- coldata$run
head( coldata[ , c( "patient", "treatment", "time" ) ] )

ddsFull <- DESeqDataSetFromMatrix(
  countData = countdata,
  colData = coldata,
  design = ~ patient + treatment)

ddsCollapsed <- collapseReplicates( ddsFull,
                                    groupby = ddsFull$sample,
                                    run = ddsFull$run )
head( as.data.frame( colData(ddsCollapsed)[ ,c("sample","runsCollapsed") ] ), 12 )


dds <- ddsCollapsed[ , ddsCollapsed$time == "48h" ]
dds$time <- droplevels( dds$time )

dds$treatment <- relevel( dds$treatment, "Control" )
as.data.frame( colData(dds) )

dds <- DESeq(dds)
res <- results( dds )

res

mcols(res, use.names=TRUE)

plotMA( res, ylim = c(-1, 1) )
plotDispEsts( dds, ylim = c(1e-6, 1e1) )
hist( res$pvalue, breaks=20, col="grey" )


rld <- rlog( dds )
par( mfrow = c( 1, 2 ) )
plot( log2( 1+counts(dds, normalized=TRUE)[, 1:2] ), col="#00000020", pch=20, cex=0.3 )
plot( assay(rld)[, 1:2], col="#00000020", pch=20, cex=0.3 )

sampleDists <- dist( t( assay(rld) ) )
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$treatment,
                                     rld$patient, sep="-" )
colnames(sampleDistMatrix) <- NULL

colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrix, trace="none", col=colours)


topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 35 )
heatmap.2( assay(rld)[ topVarGenes, ], scale="row",
           trace="none", dendrogram="column",
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))