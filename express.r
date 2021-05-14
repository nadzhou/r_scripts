if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")


BiocManager::install("GenomicFeatures")
BiocManager::install( "DESeq2" )
BiocManager::install( "Rsamtools" )
BiocManager::install("GenomicAlignments")

library( "GenomicFeatures" )
library( "DESeq2" )
library( "Rsamtools" )
library( "GenomicAlignments" )



hse <- makeTranscriptDbFromGFF( "annotation.gtf", format="gtf" )
exonsByGene <- exonsBy( hse, by="gene" )

bamLst <- BamFileList( fls, yieldSize=100000 )
fls <- list.files( "expression.bam", pattern="bam$", full=TRUE )


library( "GenomicAlignments" )
se <- summarizeOverlaps( exonsByGene, bamLst,
                         mode="Union",
                         singleEnd=FALSE,
                         ignore.strand=TRUE,
                         fragments=TRUE )


print(se)

