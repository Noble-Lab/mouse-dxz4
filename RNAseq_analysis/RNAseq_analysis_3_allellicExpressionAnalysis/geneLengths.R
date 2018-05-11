###########################################################
#
# Gene lengths & GC content
# Environment: Server
# Author: GAB
# 01/14/2015
# Adapted from version by Irsan
# Ref.2: https://www.biostars.org/p/83901/
#
###########################################################


refdataDir="<refdataDirectory>"


###########################################################
# *********************************************************
# Change input & output files. 
# *********************************************************
###########################################################

# GTF file.
GTFfile = file.path(refdataDir, "iGenomes/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf")

# Output File
geneLengthsFile = file.path(refdataDir, "iGenomes/Mus_musculus/UCSC/mm10/Annotation/Genes/GeneLengths.tsv")

###########################################################

# install.packages("GenomicFeatures")
# Didn't work, so did this:
# source("http://bioconductor.org/biocLite.R")
# biocLite("GenomicFeatures")


##############
##############
# First, import the GTF-file that you have also used as input for htseq-count

require(GenomicFeatures)
# txdb <- makeTranscriptDbFromGFF(GTFfile, format="gtf")
# # Error: 'makeTranscriptDbFromGFF' is defunct.
# # Use 'makeTxDbFromGFF' instead.
# # See help("Defunct")
txdb <- makeTxDbFromGFF(GTFfile, format="gtf")

# then collect the exons per gene id
exons.list.per.gene <- exonsBy(txdb,by="gene")

# then for each gene, reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum then
sumReducedWidths <- function(x){sum(width(reduce(x)))}

# Test
exons.list.per.gene[[1]]
sumReducedWidths(exons.list.per.gene[[1]])

exonic.gene.sizes <- lapply(exons.list.per.gene, sumReducedWidths)
length(exonic.gene.sizes)
head(exonic.gene.sizes)

output2 = as.data.frame(unlist(exonic.gene.sizes))
output2 = cbind(rownames(output2), output2)
colnames(output2) <- c("geneId", "Length")
head(output2)
nrow(output2)

write.table(output2, file=geneLengthsFile, sep="\t", row.names=FALSE, col.names=TRUE)
