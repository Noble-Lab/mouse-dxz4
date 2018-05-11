###########################################################
###########################################################
###########################################################
#
# DEseq Analysis from HTseq-count output
# 
# Author: GB
# 20161010
# Ref.1: http://www.bioconductor.org/packages/2.13/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf
# Ref.2: http://www.dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/ 
#
###########################################################
###########################################################
###########################################################




###########################################################
###########################################################
###########################################################
# This is analysis continues from where 'allelicExpressionAnalysis.R' left off.

projectDir = "<projectDirectory>"
subDir = "20170620_RNAseq_InvDxz4abAndWTab"
workDir = "Patski.combo.AllData"

outputDir = file.path(projectDir, subDir, "expressionAnalysis")
if ( ! file.exists( outputDir ) ) {
  if ( ! dir.create( outputDir, showWarnings=FALSE, recursive=TRUE )  ) outputDir = getwd()
}
setwd(outputDir)
getwd()

#save.image()
#load(file = ".RData") # From .RData





###########################################################
###########################################################
###########################################################
# Load DESeq2 and other code

require(DESeq2)
require(openxlsx)

sourceDir = "<scriptDirectory>"
# version of vioplot with colour
source(file.path(sourceDir, "colVioplot.R"))





###########################################################
###########################################################
###########################################################
# Determine allelic counts based on adjusted allelic proportion

expressionTables[["counts"]] = countData
expressionTables[["altCounts"]] = expressionTables[["counts"]] * geneLevelASEdata.uniqued.adjAltProp 
expressionTables[["refCounts"]] = expressionTables[["counts"]] * geneLevelASEdata.uniqued.adjRefProp

###########################################################
# Save expression tables with expression for all samples
tableDir = file.path(outputDir, "expressionTables")
if ( ! file.exists( tableDir ) ) {
  if ( ! dir.create( tableDir, showWarnings=FALSE, recursive=TRUE )  ) tableDir = outputDir
}

datasetSelections = colnames(countData)
datasetRenames = colnames(countData)
names(datasetRenames) = datasetSelections

for (expDataType in grep("counts$", names(expressionTables), value=TRUE, ignore.case=TRUE)) {
  bigOldTable = cbind( as.character(geneLevelASEdata.uniqued[,4]),
                       as.character(geneLevelASEdata.uniqued[,1]),
                       format(geneLevelASEdata.uniqued[,2], scientific=FALSE, trim=TRUE, justify="none"),
                       format(geneLevelASEdata.uniqued[,3], scientific=FALSE, trim=TRUE, justify="none"),      
                       expressionTables[[expDataType]][,datasetSelections]
  )
  colnames(bigOldTable) = c("geneId", "chr", "start", "stop", datasetRenames)
  bigOldFilename = file.path(tableDir, paste(expDataType, ".tsv", sep = ""))
  write.table(bigOldTable, bigOldFilename, 
              append = FALSE, sep = "\t", 
              col.names = TRUE, row.names = FALSE, quote = FALSE)
  
  try(system(paste("gzip -f \"", sub("~/", "$HOME/", bigOldFilename), "\"", sep = "")))
}


#save.image() # to .Rdata
#load(file = ".RData")











###########################################################
###########################################################
###########################################################
# DESeq2 work

whichComp = "WT-vs-Del1a"
# whichComp = "WT-vs-Del1ab"
whichComp = "newWT-vs-Del1a"
whichComp = "newWT-vs-InvDxz4"
whichComp = "WT-vs-InvDxz4"
whichComp = "WT-vs-newWT"

DESeq2Dir = file.path(outputDir, "DESeq2", whichComp )
if ( ! file.exists( DESeq2Dir ) ) {
  if ( ! dir.create( DESeq2Dir, showWarnings=FALSE, recursive=TRUE )  ) DESeq2Dir = outputDir
}

setwd(DESeq2Dir)
getwd()

#save.image()
#load(file = ".RData") # From .RData




###########################################################
###########################################################
###########################################################
# Subset relevant count data 
# And pool where necessary

# colnames(expressionTables[[1]])

NArowSums = function(x){
  ifelse(all(is.na(x)), NaN, sum(x, na.rm=TRUE))
} 

NArowMeans = function(x){
  ifelse(all(is.na(x)), NaN, mean(x, na.rm=TRUE))
} 

expressionTables.subset = list()

if (whichComp == "WT-vs-Del1a") {

  for (genomeType in names(expressionTables)) {
    # genomeType="counts"
    cat(genomeType, "\n", sep="")
    if (grepl("counts", genomeType, ignore.case=TRUE)) {
      expressionTables.subset[[genomeType]] = data.frame(expressionTables[[genomeType]][,c("PatskiWT")])
      colnames(expressionTables.subset[[genomeType]]) = "PatskiWT"
      expressionTables.subset[[genomeType]]$PatskiWT_0uMaza = apply(expressionTables[[genomeType]][, c("PatskiWT_0uMaza_rep1", "PatskiWT_0uMaza_rep2")], 1, NArowSums)
      expressionTables.subset[[genomeType]]$PatskiDel1 = expressionTables[[genomeType]][,c("PatskiDel1")]
      expressionTables.subset[[genomeType]]$PatskiDel1_0uMaza = apply(expressionTables[[genomeType]][, c("PatskiDel1_0uMaza_rep1", "PatskiDel1_0uMaza_rep2")], 1, NArowSums)
    } else {
      expressionTables.subset[[genomeType]] = data.frame(expressionTables[[genomeType]][,c("PatskiWT")])
      colnames(expressionTables.subset[[genomeType]]) = "PatskiWT"
      expressionTables.subset[[genomeType]]$PatskiWT_0uMaza = apply(expressionTables[[genomeType]][, c("PatskiWT_0uMaza_rep1", "PatskiWT_0uMaza_rep2")], 1, NArowMeans)
      expressionTables.subset[[genomeType]]$PatskiDel1 = expressionTables[[genomeType]][,c("PatskiDel1")]
      expressionTables.subset[[genomeType]]$PatskiDel1_0uMaza = apply(expressionTables[[genomeType]][, c("PatskiDel1_0uMaza_rep1", "PatskiDel1_0uMaza_rep2")], 1, NArowMeans)
    }
  }

  geneLevelASEdata.uniqued.totCov.subset = data.frame(geneLevelASEdata.uniqued.totCov[,c("PatskiWT")])
  colnames(geneLevelASEdata.uniqued.totCov.subset) = "PatskiWT"
  geneLevelASEdata.uniqued.totCov.subset$PatskiWT_0uMaza = apply(geneLevelASEdata.uniqued.totCov[, c("PatskiWT_0uMaza_rep1", "PatskiWT_0uMaza_rep2")], 1, NArowSums)
  geneLevelASEdata.uniqued.totCov.subset$PatskiDel1 = geneLevelASEdata.uniqued.totCov[,c("PatskiDel1")]
  geneLevelASEdata.uniqued.totCov.subset$PatskiDel1_0uMaza = apply(geneLevelASEdata.uniqued.totCov[, c("PatskiDel1_0uMaza_rep1", "PatskiDel1_0uMaza_rep2")], 1, NArowSums)
  
  head(geneLevelASEdata.uniqued.totCov.subset)
  head(geneLevelASEdata.uniqued.totCov[,c("PatskiWT", "PatskiWT_0uMaza_rep1", "PatskiWT_0uMaza_rep2",
                                          "PatskiDel1", "PatskiDel1_0uMaza_rep1", "PatskiDel1_0uMaza_rep2")])

} else if (whichComp == "WT-vs-Del1ab") {
 
  for (genomeType in names(expressionTables)) {
    # genomeType="counts"
    cat(genomeType, "\n", sep="")
    if (grepl("counts", genomeType, ignore.case=TRUE)) {
      expressionTables.subset[[genomeType]] = data.frame(expressionTables[[genomeType]][,c("PatskiWT")])
      colnames(expressionTables.subset[[genomeType]]) = "PatskiWT"
      expressionTables.subset[[genomeType]]$PatskiWT_0uMaza = apply(expressionTables[[genomeType]][, c("PatskiWT_0uMaza_rep1", "PatskiWT_0uMaza_rep2")], 1, NArowSums)
      expressionTables.subset[[genomeType]]$PatskiDel1 = expressionTables[[genomeType]][,c("PatskiDel1")]
      expressionTables.subset[[genomeType]]$PatskiDel1b = expressionTables[[genomeType]][,c("PatskiDel1b")]
    } else {
      expressionTables.subset[[genomeType]] = data.frame(expressionTables[[genomeType]][,c("PatskiWT")])
      colnames(expressionTables.subset[[genomeType]]) = "PatskiWT"
      expressionTables.subset[[genomeType]]$PatskiWT_0uMaza = apply(expressionTables[[genomeType]][, c("PatskiWT_0uMaza_rep1", "PatskiWT_0uMaza_rep2")], 1, NArowMeans)
      expressionTables.subset[[genomeType]]$PatskiDel1 = expressionTables[[genomeType]][,c("PatskiDel1")]
      expressionTables.subset[[genomeType]]$PatskiDel1b = expressionTables[[genomeType]][,c("PatskiDel1b")]
    }
  }
  
  geneLevelASEdata.uniqued.totCov.subset = data.frame(geneLevelASEdata.uniqued.totCov[,c("PatskiWT")])
  colnames(geneLevelASEdata.uniqued.totCov.subset) = "PatskiWT"
  geneLevelASEdata.uniqued.totCov.subset$PatskiWT_0uMaza = apply(geneLevelASEdata.uniqued.totCov[, c("PatskiWT_0uMaza_rep1", "PatskiWT_0uMaza_rep2")], 1, NArowSums)
  geneLevelASEdata.uniqued.totCov.subset$PatskiDel1 = geneLevelASEdata.uniqued.totCov[,c("PatskiDel1")]
  geneLevelASEdata.uniqued.totCov.subset$PatskiDel1 = geneLevelASEdata.uniqued.totCov[,c("PatskiDel1")]
  geneLevelASEdata.uniqued.totCov.subset$PatskiDel1b = geneLevelASEdata.uniqued.totCov[,c("PatskiDel1b")]  
  
  head(geneLevelASEdata.uniqued.totCov.subset)
  head(geneLevelASEdata.uniqued.totCov[,c("PatskiWT", "PatskiWT_0uMaza_rep1", "PatskiWT_0uMaza_rep2",
                                          "PatskiDel1", "PatskiDel1b")])
     
} else if (whichComp == "newWT-vs-Del1a") {
  
  for (genomeType in names(expressionTables)) {
    # genomeType="counts"
    cat(genomeType, "\n", sep="")
    if (grepl("counts", genomeType, ignore.case=TRUE)) {
      expressionTables.subset[[genomeType]] = data.frame(expressionTables[[genomeType]][,c("PatskiWTa")])
      colnames(expressionTables.subset[[genomeType]]) = "PatskiWTa"
      expressionTables.subset[[genomeType]]$PatskiWTb = expressionTables[[genomeType]][,c("PatskiWTb")]
      expressionTables.subset[[genomeType]]$PatskiDel1 = expressionTables[[genomeType]][,c("PatskiDel1")]
      expressionTables.subset[[genomeType]]$PatskiDel1_0uMaza = apply(expressionTables[[genomeType]][, c("PatskiDel1_0uMaza_rep1", "PatskiDel1_0uMaza_rep2")], 1, NArowSums)
    } else {
      expressionTables.subset[[genomeType]] = data.frame(expressionTables[[genomeType]][,c("PatskiWTa")])
      colnames(expressionTables.subset[[genomeType]]) = "PatskiWTa"
      expressionTables.subset[[genomeType]]$PatskiWTb = expressionTables[[genomeType]][,c("PatskiWTb")]
      expressionTables.subset[[genomeType]]$PatskiDel1 = expressionTables[[genomeType]][,c("PatskiDel1")]
      expressionTables.subset[[genomeType]]$PatskiDel1_0uMaza = apply(expressionTables[[genomeType]][, c("PatskiDel1_0uMaza_rep1", "PatskiDel1_0uMaza_rep2")], 1, NArowMeans)
    }
  }
  
  geneLevelASEdata.uniqued.totCov.subset = data.frame(geneLevelASEdata.uniqued.totCov[,c("PatskiWTa")])
  colnames(geneLevelASEdata.uniqued.totCov.subset) = "PatskiWTa"
  geneLevelASEdata.uniqued.totCov.subset$PatskiWTb = geneLevelASEdata.uniqued.totCov[,c("PatskiWTb")]
  geneLevelASEdata.uniqued.totCov.subset$PatskiDel1 = geneLevelASEdata.uniqued.totCov[,c("PatskiDel1")]
  geneLevelASEdata.uniqued.totCov.subset$PatskiDel1_0uMaza = apply(geneLevelASEdata.uniqued.totCov[, c("PatskiDel1_0uMaza_rep1", "PatskiDel1_0uMaza_rep2")], 1, NArowSums)
  
  head(geneLevelASEdata.uniqued.totCov.subset)
  head(geneLevelASEdata.uniqued.totCov[,c("PatskiWTa", "PatskiWTb", "PatskiDel1", "PatskiDel1_0uMaza_rep1", "PatskiDel1_0uMaza_rep2")])
  
} else if (whichComp == "newWT-vs-InvDxz4") {
  
  for (genomeType in names(expressionTables)) {
    # genomeType="counts"
    cat(genomeType, "\n", sep="")
    expressionTables.subset[[genomeType]] = expressionTables[[genomeType]][, c("PatskiWTa", "PatskiWTb", "PatskiInvDxz4a", "PatskiInvDxz4b")]
  }
  geneLevelASEdata.uniqued.totCov.subset = geneLevelASEdata.uniqued.totCov[, c("PatskiWTa", "PatskiWTb", "PatskiInvDxz4a", "PatskiInvDxz4b")]
  
  head(geneLevelASEdata.uniqued.totCov.subset)
  head(geneLevelASEdata.uniqued.totCov[,c("PatskiWTa", "PatskiWTb", "PatskiInvDxz4a", "PatskiInvDxz4b")])

  
} else if (whichComp == "WT-vs-InvDxz4") {
  
  for (genomeType in names(expressionTables)) {
    # genomeType="counts"
    cat(genomeType, "\n", sep="")
    if (grepl("counts", genomeType, ignore.case=TRUE)) {
      expressionTables.subset[[genomeType]] = data.frame(expressionTables[[genomeType]][,c("PatskiWT")])
      colnames(expressionTables.subset[[genomeType]]) = "PatskiWT"
      expressionTables.subset[[genomeType]]$PatskiWT_0uMaza = apply(expressionTables[[genomeType]][, c("PatskiWT_0uMaza_rep1", "PatskiWT_0uMaza_rep2")], 1, NArowSums)
      expressionTables.subset[[genomeType]]$PatskiInvDxz4a = expressionTables[[genomeType]][,c("PatskiInvDxz4a")]
      expressionTables.subset[[genomeType]]$PatskiInvDxz4b = expressionTables[[genomeType]][,c("PatskiInvDxz4b")]
    } else {
      expressionTables.subset[[genomeType]] = data.frame(expressionTables[[genomeType]][,c("PatskiWT")])
      colnames(expressionTables.subset[[genomeType]]) = "PatskiWT"
      expressionTables.subset[[genomeType]]$PatskiWT_0uMaza = apply(expressionTables[[genomeType]][, c("PatskiWT_0uMaza_rep1", "PatskiWT_0uMaza_rep2")], 1, NArowMeans)
      expressionTables.subset[[genomeType]]$PatskiInvDxz4a = expressionTables[[genomeType]][,c("PatskiInvDxz4a")]
      expressionTables.subset[[genomeType]]$PatskiInvDxz4b = expressionTables[[genomeType]][,c("PatskiInvDxz4b")]
    }
  }
  
  geneLevelASEdata.uniqued.totCov.subset = data.frame(geneLevelASEdata.uniqued.totCov[,c("PatskiWT")])
  colnames(geneLevelASEdata.uniqued.totCov.subset) = "PatskiWT"
  geneLevelASEdata.uniqued.totCov.subset$PatskiWT_0uMaza = apply(geneLevelASEdata.uniqued.totCov[, c("PatskiWT_0uMaza_rep1", "PatskiWT_0uMaza_rep2")], 1, NArowSums)
  geneLevelASEdata.uniqued.totCov.subset$PatskiInvDxz4a = geneLevelASEdata.uniqued.totCov[,c("PatskiInvDxz4a")]
  geneLevelASEdata.uniqued.totCov.subset$PatskiInvDxz4b = geneLevelASEdata.uniqued.totCov[,c("PatskiInvDxz4b")]
  
  head(geneLevelASEdata.uniqued.totCov.subset)
  head(geneLevelASEdata.uniqued.totCov[,c("PatskiWT", "PatskiWT_0uMaza_rep1", "PatskiWT_0uMaza_rep2", "PatskiInvDxz4a", "PatskiInvDxz4b")])
  
} else if (whichComp == "WT-vs-newWT") {
  
  for (genomeType in names(expressionTables)) {
    # genomeType="counts"
    cat(genomeType, "\n", sep="")
    if (grepl("counts", genomeType, ignore.case=TRUE)) {
      expressionTables.subset[[genomeType]] = data.frame(expressionTables[[genomeType]][,c("PatskiWT")])
      colnames(expressionTables.subset[[genomeType]]) = "PatskiWT"
      expressionTables.subset[[genomeType]]$PatskiWT_0uMaza = apply(expressionTables[[genomeType]][, c("PatskiWT_0uMaza_rep1", "PatskiWT_0uMaza_rep2")], 1, NArowSums)
      expressionTables.subset[[genomeType]]$PatskiWTa = expressionTables[[genomeType]][,c("PatskiWTa")]
      expressionTables.subset[[genomeType]]$PatskiWTb = expressionTables[[genomeType]][,c("PatskiWTb")]
    } else {
      expressionTables.subset[[genomeType]] = data.frame(expressionTables[[genomeType]][,c("PatskiWT")])
      colnames(expressionTables.subset[[genomeType]]) = "PatskiWT"
      expressionTables.subset[[genomeType]]$PatskiWT_0uMaza = apply(expressionTables[[genomeType]][, c("PatskiWT_0uMaza_rep1", "PatskiWT_0uMaza_rep2")], 1, NArowMeans)
      expressionTables.subset[[genomeType]]$PatskiWTa = expressionTables[[genomeType]][,c("PatskiWTa")]
      expressionTables.subset[[genomeType]]$PatskiWTb = expressionTables[[genomeType]][,c("PatskiWTb")]
  }
  
  geneLevelASEdata.uniqued.totCov.subset = data.frame(geneLevelASEdata.uniqued.totCov[,c("PatskiWT")])
  colnames(geneLevelASEdata.uniqued.totCov.subset) = "PatskiWT"
  geneLevelASEdata.uniqued.totCov.subset$PatskiWT_0uMaza = apply(geneLevelASEdata.uniqued.totCov[, c("PatskiWT_0uMaza_rep1", "PatskiWT_0uMaza_rep2")], 1, NArowSums)
  geneLevelASEdata.uniqued.totCov.subset$PatskiWTa = geneLevelASEdata.uniqued.totCov[,c("PatskiWTa")]
  geneLevelASEdata.uniqued.totCov.subset$PatskiWTb = geneLevelASEdata.uniqued.totCov[,c("PatskiWTb")]
  
  head(geneLevelASEdata.uniqued.totCov.subset)
  head(geneLevelASEdata.uniqued.totCov[,c("PatskiWT", "PatskiWT_0uMaza_rep1", "PatskiWT_0uMaza_rep2", "PatskiWTa", "PatskiWTb")])
  
}


###########################################################
# Common meta data for design matrix 

cat(colnames(expressionTables.subset[["counts"]]), sep="\n")
# PatskiWT/ PatskiWTa
# PatskiWT_0uMaza / PatskiWTb
# PatskiDel1
# PatskiDel1_0uMaza

dds = dds.diffexp = list()

# Create meta-data dataframe to correspond with countData
if (grepl("WT-vs-InvDxz4", whichComp)) {
  countData.design  <- data.frame(row.names=colnames(expressionTables.subset[["counts"]]),
                                  condition=as.factor(c(rep("WT", 2), rep("InvDxz4", 2))),
                                  libtype=rep("single-end", 4))
} else if (grepl("WT-vs-newWT", whichComp)) {
  countData.design  <- data.frame(row.names=colnames(expressionTables.subset[["counts"]]),
                                  condition=as.factor(c(rep("WT", 2), rep("newWT", 2))),
                                  libtype=rep("single-end", 4))
} else {
  countData.design  <- data.frame(row.names=colnames(expressionTables.subset[["counts"]]),
                                  condition=as.factor(c(rep("WT", 2), rep("Del1", 2))),
                                  libtype=rep("single-end", 4))
}


###########################################################
###########################################################
# Diploid differential expression
# Note: no prefiltering performed.

# Construct a DESeqDataSet
dds[["diploid"]] <- DESeqDataSetFromMatrix(countData = expressionTables.subset[["counts"]],
                              colData   = countData.design,
                              design    = ~ condition)
head(assay(dds[["diploid"]]))

# Run DESeq's differenctial expression analysis
dds.diffexp[["diploid"]] <- DESeq(dds[["diploid"]])
# Same as DESeq.v1 built-in functions below:
# dds <- estimateSizeFactors(dds)
# dds <- estimateDispersions(dds)
# dds <- nbinomWaldTest(dds)

# Standard differential expression analysis processing:
#    - estimating size factors
#    - estimating dispersions
#    - gene-wise dispersion estimates
#    - mean-dispersion relationship
#    - final dispersion estimates
#    - fitting model and testing

resultsNames(dds.diffexp[["diploid"]])
attributes(dds.diffexp[["diploid"]])
# Compute size factors 
# size factors normalize each sample to the reference for one scaling factor for counts
sizeFactors(dds.diffexp[["diploid"]])




##########################################
#
# If one wants to compute size factors from scratch, use code below:
#
# geomeans <- exp( rowMeans( log( counts(dds.diffexp) ) ) )
# median( ( counts(dds.diffexp)[,12]/geomeans )[ geomeans>0 ] ) 
#
##########################################








###########################################################
###########################################################
# Genes with sufficient coverage to infer allelic expression data
# This is a form of pre-filtering really.

minSNPcov = 5

# Prefilter data:
# Remove any genes with NAs
naBoolIdx = !apply(apply(expressionTables.subset[["refCounts"]], 1, is.na), 2, any)
head(naBoolIdx)

# Note: this is the same for altCounts
all(naBoolIdx == !apply(apply(expressionTables.subset[["altCounts"]], 1, is.na), 2, any))
# TRUE

# Remove any genes with coverage less than threshold
minSNPcovBoolIdx = apply(geneLevelASEdata.uniqued.totCov.subset, 1, min, na.rm=TRUE) >= minSNPcov

# Consider both
minSNPcovBoolIdx = minSNPcovBoolIdx & naBoolIdx

head(geneLevelASEdata.uniqued.totCov.subset)
head(expressionTables.subset[["refCounts"]])
head(minSNPcovBoolIdx)

# %
sum(minSNPcovBoolIdx) / length(minSNPcovBoolIdx) * 100


###########################################################
###########################################################
# black6 differential expression

# Construct a DESeqDataSet
dds[["B6"]] <- DESeqDataSetFromMatrix(countData = round(data.matrix(expressionTables.subset[["refCounts"]][minSNPcovBoolIdx,])),
                                           colData   = countData.design,
                                           design    = ~ condition)
head(assay(dds[["B6"]]))

# Run DESeq's differenctial expression analysis
dds.diffexp[["B6"]] <- DESeq(dds[["B6"]])
# Same as DESeq.v1 built-in functions below:
# dds <- estimateSizeFactors(dds)
# dds <- estimateDispersions(dds)
# dds <- nbinomWaldTest(dds)

# Standard differential expression analysis processing:
#    - estimating size factors
#    - estimating dispersions
#    - gene-wise dispersion estimates
#    - mean-dispersion relationship
#    - final dispersion estimates
#    - fitting model and testing

attributes(dds.diffexp[["B6"]])
# Compute size factors 
# size factors normalize each sample to the reference for one scaling factor for counts
sizeFactors(dds.diffexp[["B6"]])






###########################################################
###########################################################
# spretus differential expression

# Construct a DESeqDataSet
dds[["sp"]] <- DESeqDataSetFromMatrix(countData = round(data.matrix(expressionTables.subset[["altCounts"]][minSNPcovBoolIdx,])),
                                      colData   = countData.design,
                                      design    = ~ condition)
head(assay(dds[["sp"]]))

# Run DESeq's differenctial expression analysis
dds.diffexp[["sp"]] <- DESeq(dds[["sp"]])
# Same as DESeq.v1 built-in functions below:
# dds <- estimateSizeFactors(dds)
# dds <- estimateDispersions(dds)
# dds <- nbinomWaldTest(dds)

# Standard differential expression analysis processing:
#    - estimating size factors
#    - estimating dispersions
#    - gene-wise dispersion estimates
#    - mean-dispersion relationship
#    - final dispersion estimates
#    - fitting model and testing

attributes(dds.diffexp[["sp"]])
# Compute size factors 
# size factors normalize each sample to the reference for one scaling factor for counts
sizeFactors(dds.diffexp[["sp"]])





###########################################################
###########################################################
# Tables and plots

DEadjPthresh = 0.05; DElogFCthresh = 0.5
# DEadjPthresh = 0.01; DElogFCthresh = 1
threshText = paste("logP=", round(abs(log10(DEadjPthresh)), 2), "_", "logFC=", DElogFCthresh, sep="")

# samples
DEsamples = unique(countData.design$condition)
pwDEsamples = combn(DEsamples, 2)

###########################################################
# DE tables

DEresults = list()
for (genomeType in names(dds.diffexp)) {
  cat(genomeType, "\n", sep="")
  DEresults[[genomeType]] = list()
  for (pwIdx in 1:ncol(pwDEsamples)) {
    sampl1 = as.character(pwDEsamples[,pwIdx][1])
    sampl2 = as.character(pwDEsamples[,pwIdx][2])
    DEdesc = paste(sampl1, "-vs-", sampl2, sep="")
    cat(DEdesc, "\n", sep="")
    DEresults[[genomeType]][[DEdesc]] = as.data.frame(results(dds.diffexp[[genomeType]], contrast=c('condition', sampl1, sampl2)))
  }
}


# Note: these parameters are implied:
# independentFiltering=TRUE, alpha=0.1, pAdjustMethod="BH", lfcThreshold=0
#
# From ?results:
#
# lfcThreshold 	a non-negative value, which specifies the test which should be applied to the log2 fold changes. 
#               The standard is a test that the log2 fold changes are not equal to zero. 
#               However, log2 fold changes greater or less than lfcThreshold can also be tested. 
#               Specify the alternative hypothesis using the altHypothesis argument. 
#               If lfcThreshold is specified, the results are Wald tests, and LRT p-values will be overwritten.
#
# By default, independent filtering is performed to select a set of genes for multiple test correction 
# which will optimize the number of adjusted p-values less than a given critical value alpha (by default 0.1). 
# The adjusted p-values for the genes which do not pass the filter threshold are set to NA. 
# By default, the mean of normalized counts is used to perform this filtering, though other statistics can be provided. 
# Several arguments from the filtered_p function of genefilter are provided here to control or turn off the independent filtering behavior.
#
# By default, results assigns a p-value of NA to genes containing count outliers, as identified using Cook's distance (by default 0.99). 
# See the cooksCutoff argument for control of this behavior. Cook's distances for each sample are accessible as a matrix "cooks" stored in the assays() list. 
# This measure is useful for identifying rows where the observed counts might not fit to a Negative Binomial distribution.


# Save DE results for all genes
for (genomeType in names(dds.diffexp)) {
  cat(genomeType, "\n", sep="")
  ouputfilename = file.path(DESeq2Dir, 
                            paste(workDir, "_DESeq2genes_", threshText, "_", genomeType, ".xlsx", sep=""))
  write.xlsx(DEresults[[genomeType]], file=ouputfilename, colNames=TRUE, rowNames=TRUE, asTable=TRUE, keepNA=TRUE)  

  ouputfilename = file.path(DESeq2Dir, 
                            paste(workDir, "_DESeq2genes_", threshText, "_", genomeType, ".tsv", sep=""))
  write.table(DEresults[[genomeType]], ouputfilename, row.names=TRUE, col.names=NA, sep="\t", quote=FALSE) 
}


# Save DE results for all chrX genes
DEresults.chrX = list()
whichChrom = "chrX"
for (genomeType in names(dds.diffexp)) {
  cat(genomeType, "\n", sep="")
  
  # Note: DEresults different genes for diploid and haploid analyses
  # i.e. subset chrBoolindices by minSNPcovBoolIdx
  if (genomeType != "diploid") {
    geneSubsetBoolIdx = chrBoolindices[[whichChrom]][minSNPcovBoolIdx]
  } else {
    geneSubsetBoolIdx = chrBoolindices[[whichChrom]]
  }
  
  for (DEdesc in names(DEresults[[genomeType]])) {
    DEresults.chrX[[genomeType]][[DEdesc]] = DEresults[[genomeType]][[DEdesc]][geneSubsetBoolIdx,]
    
    ouputfilename = file.path(DESeq2Dir, 
                              paste(workDir, "_DESeq2genes_", threshText, "_", genomeType, ".", DEdesc, ".chrX.tsv", sep=""))
    write.table(DEresults.chrX[[genomeType]], ouputfilename, row.names=TRUE, col.names=NA, sep="\t", quote=FALSE) 
  }
  
  ouputfilename = file.path(DESeq2Dir, 
                            paste(workDir, "_DESeq2genes_", threshText, "_", genomeType, ".chrX.xlsx", sep=""))
  write.xlsx(DEresults.chrX[[genomeType]], file=ouputfilename, colNames=TRUE, rowNames=TRUE, asTable=TRUE, keepNA=TRUE)  
}


# Save sig results for all chrX genes
DEresults.chrX.sig = list()
whichChrom = "chrX"
for (genomeType in names(dds.diffexp)) {
  cat(genomeType, "\n", sep="")
  
  # Note: DEresults different genes for diploid and haploid analyses
  # i.e. subset chrBoolindices by minSNPcovBoolIdx
  if (genomeType != "diploid") {
    geneSubsetBoolIdx = chrBoolindices[[whichChrom]][minSNPcovBoolIdx]
  } else {
    geneSubsetBoolIdx = chrBoolindices[[whichChrom]]
  }
  
  for (DEdesc in names(DEresults[[genomeType]])) {
    DEgeneBoolIdx = DEresults[[genomeType]][[DEdesc]][geneSubsetBoolIdx,]$padj <= DEadjPthresh &
      abs(DEresults[[genomeType]][[DEdesc]][geneSubsetBoolIdx,]$log2FoldChange) >= DElogFCthresh
    DEresults.chrX.sig[[genomeType]][[DEdesc]] = DEresults[[genomeType]][[DEdesc]][geneSubsetBoolIdx,][which(DEgeneBoolIdx),]
    
    ouputfilename = file.path(DESeq2Dir, 
                              paste(workDir, "_DESeq2genes_", threshText, "_", genomeType, ".", DEdesc, ".chrX.sig.tsv", sep=""))
    write.table(DEresults.chrX.sig[[genomeType]][[DEdesc]], ouputfilename, row.names=TRUE, col.names=NA, sep="\t", quote=FALSE) 
    
    DEidx = match(rownames(DEresults.chrX.sig[[genomeType]][[DEdesc]]), geneLevelASEdata.uniqued$geneID)
    DElist = geneLevelASEdata.uniqued[DEidx, c(4,1,2,3,6)]
    DElist = DElist[order(DElist$start),]
    
    ouputfilename = file.path(DESeq2Dir, 
                              paste(workDir, "_DESeq2genes_", threshText, "_", genomeType, ".", DEdesc, ".chrX.sig.geneList.tsv", sep=""))
    write.table(DElist, ouputfilename, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE) 
  }
  
  ouputfilename = file.path(DESeq2Dir, 
                            paste(workDir, "_DESeq2genes_", threshText, "_", genomeType, ".chrX.sig.xlsx", sep=""))
  write.xlsx(DEresults.chrX.sig[[genomeType]] , file=ouputfilename, colNames=TRUE, rowNames=TRUE, asTable=TRUE, keepNA=TRUE)  
}




###########################################################
# MA plots

for (pwIdx in 1:ncol(pwDEsamples)) {
  # pwIdx = 1
  sampl1 = as.character(pwDEsamples[,pwIdx][1])
  sampl2 = as.character(pwDEsamples[,pwIdx][2])
  DEdesc = paste(sampl1, "-vs-", sampl2, sep="")

  tsampl1 = sampl1
  tsampl2 = sampl2
  titleDesc = paste(tsampl1, " vs ", tsampl2, sep="")
  titleDescAlt = paste(tsampl2, " vs ", tsampl1, sep="")
  cat(DEdesc, "; ", titleDesc, "; ", titleDescAlt, "\n", sep="")
  
  for (whichChrom in c("autosomes", "chrX")) {
    
    # png(file.path(DESeq2Dir, paste(workDir, "_DESeq2_MAplot_", DEdesc, "_", whichChrom, ".png", sep="")), type='cairo', width=15, height=5, units="in", res=150)
    # # had to use type='cairo'
    # # https://stat.ethz.ch/R-manual/R-devel/library/grDevices/html/png.html
    # # getOption("bitmapType")
    pdf(file.path(DESeq2Dir, paste(workDir, "_DESeq2_MAplot_", threshText, "_", DEdesc, "_", whichChrom, ".pdf", sep="")), width=15, height=5)
    par(mfrow=c(1,3), oma=c(0,0,2,0))
    
    maxLog10mean = max(log10(DEresults[["diploid"]][[DEdesc]][chrBoolindices[[whichChrom]],]$baseMean + 1), na.rm=TRUE)
    maxLogFC = max(abs(DEresults[["diploid"]][[DEdesc]][chrBoolindices[[whichChrom]],]$log2FoldChange), na.rm=TRUE)
    
    for (genomeType in c("diploid", "B6", "sp")) {
      
      # Note: DEresults different genes for diploid and haploid analyses
      # i.e. subset chrBoolindices by minSNPcovBoolIdx
      if (genomeType != "diploid") {
        geneSubsetBoolIdx = chrBoolindices[[whichChrom]][minSNPcovBoolIdx]
      } else {
        geneSubsetBoolIdx = chrBoolindices[[whichChrom]]
      }
      DEgeneBoolIdx = DEresults[[genomeType]][[DEdesc]][geneSubsetBoolIdx,]$padj <= DEadjPthresh &
        abs(DEresults[[genomeType]][[DEdesc]][geneSubsetBoolIdx,]$log2FoldChange) >= DElogFCthresh

      #   plotMA(dds.diffexp[[genomeType]][chrBoolindices[[whichChrom]],], contrast=c('condition', sampl1, sampl2), 
      #          ylim=c(-5,5), main=paste(DEdesc, "diploid expression"), alpha=0)
      
      plot(log10(DEresults[[genomeType]][[DEdesc]][geneSubsetBoolIdx,]$baseMean + 1),
           DEresults[[genomeType]][[DEdesc]][geneSubsetBoolIdx,]$log2FoldChange,
           col="gray", pch=20, cex=0.4, ylim=c(-maxLogFC,maxLogFC), xlim=c(0, maxLog10mean),
           main=ifelse(whichChrom=="autosomes",
                       paste(titleDesc, genomeType, "expression"), 
                       paste(titleDesc, genomeType, "expression", 
                             ifelse(genomeType=="diploid", "(chrX)", 
                                    ifelse(genomeType=="B6", "(Xi)", "(Xa)")))),
           xlab="mean counts (log10)", ylab="log fold change",
           xaxt = "n")
      axis(1, at=log10(10^(1:5)), labels=format(10^(1:5), scientific=TRUE))
      abline(h=0, col="cyan", lwd=2)
      points(log10(DEresults[[genomeType]][[DEdesc]][geneSubsetBoolIdx,]$baseMean + 1)[DEgeneBoolIdx],
             DEresults[[genomeType]][[DEdesc]][geneSubsetBoolIdx,]$log2FoldChange[DEgeneBoolIdx],
             col="red", cex=0.5)
      #op <- par(family = "Courier New")
      legend("topright", 
             legend=c(paste("Genes:", sum(geneSubsetBoolIdx)),
                      paste("DE:   ", sum(DEgeneBoolIdx, na.rm=TRUE)),
                      paste("Up:   ", sum(DEgeneBoolIdx & DEresults[[genomeType]][[DEdesc]][geneSubsetBoolIdx,]$log2FoldChange > 0, na.rm=TRUE)),
                      paste("Down: ", sum(DEgeneBoolIdx & DEresults[[genomeType]][[DEdesc]][geneSubsetBoolIdx,]$log2FoldChange < 0, na.rm=TRUE)))
      )
      # par(op)
    }  
    
    title(main=whichChrom, outer=TRUE)
    dev.off()
    
    
    # Alternative inverted plots
    # png(file.path(DESeq2Dir, paste(workDir, "_DESeq2_MAplot_", DEdesc, "_", whichChrom, ".png", sep="")), type='cairo', width=15, height=5, units="in", res=150)
    # # had to use type='cairo'
    # # https://stat.ethz.ch/R-manual/R-devel/library/grDevices/html/png.html
    # # getOption("bitmapType")
    pdf(file.path(DESeq2Dir, paste(workDir, "_DESeq2_MAplot_", threshText, "_", DEdesc, "_", whichChrom, ".V2.pdf", sep="")), width=15, height=5)
    par(mfrow=c(1,3), oma=c(0,0,2,0))
    
    maxLog10mean = max(log10(DEresults[["diploid"]][[DEdesc]][chrBoolindices[[whichChrom]],]$baseMean + 1), na.rm=TRUE)
    maxLogFC = max(abs(DEresults[["diploid"]][[DEdesc]][chrBoolindices[[whichChrom]],]$log2FoldChange), na.rm=TRUE)
    
    for (genomeType in c("diploid", "B6", "sp")) {
      
      # Note: DEresults different genes for diploid and haploid analyses
      # i.e. subset chrBoolindices by minSNPcovBoolIdx
      if (genomeType != "diploid") {
        geneSubsetBoolIdx = chrBoolindices[[whichChrom]][minSNPcovBoolIdx]
      } else {
        geneSubsetBoolIdx = chrBoolindices[[whichChrom]]
      }
      DEgeneBoolIdx = DEresults[[genomeType]][[DEdesc]][geneSubsetBoolIdx,]$padj <= DEadjPthresh &
        abs(DEresults[[genomeType]][[DEdesc]][geneSubsetBoolIdx,]$log2FoldChange) >= DElogFCthresh

      #   plotMA(dds.diffexp[[genomeType]][chrBoolindices[[whichChrom]],], contrast=c('condition', sampl1, sampl2), 
      #          ylim=c(-5,5), main=paste(DEdesc, "diploid expression"), alpha=0)
      
      plot(log10(DEresults[[genomeType]][[DEdesc]][geneSubsetBoolIdx,]$baseMean + 1),
           -1*(DEresults[[genomeType]][[DEdesc]][geneSubsetBoolIdx,]$log2FoldChange),
           col="gray", pch=20, cex=0.4, ylim=c(-maxLogFC,maxLogFC), xlim=c(0, maxLog10mean),
           main=ifelse(whichChrom=="autosomes",
                       paste(titleDescAlt, genomeType, "expression"), 
                       paste(titleDescAlt, genomeType, "expression", 
                             ifelse(genomeType=="diploid", "(chrX)", 
                                    ifelse(genomeType=="B6", "(Xi)", "(Xa)")))),
           xlab="mean counts (log10)", ylab="log fold change",
           xaxt = "n")
      axis(1, at=log10(10^(1:5)), labels=format(10^(1:5), scientific=TRUE))
      abline(h=0, col="cyan", lwd=2)
      points(log10(DEresults[[genomeType]][[DEdesc]][geneSubsetBoolIdx,]$baseMean + 1)[DEgeneBoolIdx],
             -1*(DEresults[[genomeType]][[DEdesc]][geneSubsetBoolIdx,]$log2FoldChange[DEgeneBoolIdx]),
             col="red", cex=0.5)
      #op <- par(family = "Courier New")
      legend("topright", 
             legend=c(paste("Genes:", sum(geneSubsetBoolIdx)),
                      paste("DE:   ", sum(DEgeneBoolIdx, na.rm=TRUE)),
                      paste("Up:   ", sum(DEgeneBoolIdx & DEresults[[genomeType]][[DEdesc]][geneSubsetBoolIdx,]$log2FoldChange < 0, na.rm=TRUE)),
                      paste("Down: ", sum(DEgeneBoolIdx & DEresults[[genomeType]][[DEdesc]][geneSubsetBoolIdx,]$log2FoldChange > 0, na.rm=TRUE)))
      )
      # par(op)
    }  
    
    title(main=whichChrom, outer=TRUE)
    dev.off()    
    
  }
}


###########################################################
# -logP distributions filtered by significance

diffTests = list()

for (pwIdx in 1:ncol(pwDEsamples)) {
  # pwIdx = 1
  sampl1 = as.character(pwDEsamples[,pwIdx][1])
  sampl2 = as.character(pwDEsamples[,pwIdx][2])
  DEdesc = paste(sampl1, "-vs-", sampl2, sep="")
  cat(DEdesc, "\n", sep="")

  logPList = list()
  altlogPListLabel = list()
 
  #for (genomeType in c("diploid", "B6", "sp")) {
  for (whichChrom in c("autosomes", "chrX")) {
    for (genomeType in c("sp", "B6")) {
      
      # Note: DEresults different genes for diploid and haploid analyses
      # i.e. subset chrBoolindices by minSNPcovBoolIdx
      if (genomeType != "diploid") {
        geneSubsetBoolIdx = chrBoolindices[[whichChrom]][minSNPcovBoolIdx]
      } else {
        geneSubsetBoolIdx = chrBoolindices[[whichChrom]]
      }
      DEgeneBoolIdx = DEresults[[genomeType]][[DEdesc]][geneSubsetBoolIdx,]$padj <= DEadjPthresh &
        abs(DEresults[[genomeType]][[DEdesc]][geneSubsetBoolIdx,]$log2FoldChange) >= DElogFCthresh
      
      logPlistBoolIdx = DEgeneBoolIdx & DEresults[[genomeType]][[DEdesc]][geneSubsetBoolIdx,]$log2FoldChange > 0
      logPlistLabel = paste(sampl1, "_", whichChrom, "_", genomeType,  sep="")
      logPList[[logPlistLabel]] = as.vector(na.omit(abs(DEresults[[genomeType]][[DEdesc]][geneSubsetBoolIdx,]$log2FoldChange[logPlistBoolIdx])))
      altlogPListLabel[[logPlistLabel]] = paste("Up in ", sampl1, "\n", 
                                                whichChrom, " ",
                                                ifelse(genomeType=="autosomes", "chrA", genomeType), 
                                                " (n=", sum(logPlistBoolIdx, na.rm=TRUE), ")", sep="")
    }
  }

  #for (genomeType in c("diploid", "B6", "sp")) {
  for (whichChrom in c("autosomes", "chrX")) {
    for (genomeType in c("sp", "B6")) {
      
      # Note: DEresults different genes for diploid and haploid analyses
      # i.e. subset chrBoolindices by minSNPcovBoolIdx
      if (genomeType != "diploid") {
        geneSubsetBoolIdx = chrBoolindices[[whichChrom]][minSNPcovBoolIdx]
      } else {
        geneSubsetBoolIdx = chrBoolindices[[whichChrom]]
      }
      DEgeneBoolIdx = DEresults[[genomeType]][[DEdesc]][geneSubsetBoolIdx,]$padj <= DEadjPthresh &
        abs(DEresults[[genomeType]][[DEdesc]][geneSubsetBoolIdx,]$log2FoldChange) >= DElogFCthresh
      
      logPlistBoolIdx = DEgeneBoolIdx & DEresults[[genomeType]][[DEdesc]][geneSubsetBoolIdx,]$log2FoldChange < 0
      logPlistLabel = paste(sampl2, "_", whichChrom, "_", genomeType,  sep="")
      logPList[[logPlistLabel]] = as.vector(na.omit(abs(DEresults[[genomeType]][[DEdesc]][geneSubsetBoolIdx,]$log2FoldChange[logPlistBoolIdx])))
      altlogPListLabel[[logPlistLabel]] = paste("Up in ", sampl2, "\n", 
                                                whichChrom, " ",
                                                ifelse(genomeType=="autosomes", "chrA", genomeType), 
                                                " (n=", sum(logPlistBoolIdx, na.rm=TRUE), ")", sep="")
    }  
  }

  maxLogFC = max(unlist(logPList), na.rm=TRUE)
  
  # To add transparency to colors
  alphaLvl = 0.1
  txparentColsTable = cbind(t(col2rgb(c("red", "blue"), alpha=FALSE)), alphaLvl*255)/255
  colnames(txparentColsTable) = c("red", "green", "blue", "alpha")
  txparentCols4Dists = c()
  for (colRowIdx in 1:nrow(txparentColsTable)) {
    colRowAsList = as.list(txparentColsTable[colRowIdx,])
    txparentCols4Dists = c(txparentCols4Dists, do.call(rgb, colRowAsList))
  }
  alphaLvl = 0.5
  txparentColsTable = cbind(t(col2rgb(c("red", "blue"), alpha=FALSE)), alphaLvl*255)/255
  colnames(txparentColsTable) = c("red", "green", "blue", "alpha")
  txparentCols4Strip = c()
  for (colRowIdx in 1:nrow(txparentColsTable)) {
    colRowAsList = as.list(txparentColsTable[colRowIdx,])
    txparentCols4Strip = c(txparentCols4Strip, do.call(rgb, colRowAsList))
  }
  
  # png(file.path(DESeq2Dir, 
  #               paste(workDir, "_DESeq2_", threshText, "_Boxplots_", DEdesc, ".png", sep="")), 
  #     type='cairo',  width=16, height=10, units="in", res=150)
  # # had to use type='cairo'
  # # https://stat.ethz.ch/R-manual/R-devel/library/grDevices/html/png.html
  # # getOption("bitmapType")
  pdf(file.path(DESeq2Dir, 
                paste(workDir, "_DESeq2_", threshText, "_Boxplots_", DEdesc, ".pdf", sep="")), 
      width=16, height=10)  # had to use type='cairo'
  par(mfrow=c(1,1), oma=c(5,0,2,0))
  boxplot(logPList, names=altlogPListLabel, ylim=c(0,maxLogFC), outline=FALSE,
          main=paste(DEdesc, genomeType, "log2 fold expression distibutions"),
          ylab="log fold change of signicantly differentially expressed genes", las=2, xaxs="i",
          col=txparentCols4Dists)
  stripchart(logPList,  method="jitter", jitter=0.35, pch=20, cex=1, vertical=TRUE,
             # panel.first = grid(),
             col=txparentCols4Strip, add=TRUE)  
  boxplot(logPList, names=altlogPListLabel, ylim=c(0,maxLogFC), outline=FALSE, las=2, col=txparentCols4Dists, add=TRUE )
  abline(v=2.5, lwd=1, lty="dashed", col="black")
  abline(v=4.5, lwd=1, col="black")
  abline(v=6.5, lwd=1, lty="dashed", col="black")
  dev.off()

  logPList4Viol = list()
  for (dataCat in names(logPList)) {
    if (length(logPList[[dataCat]]) == 0) {
      logPList4Viol[[dataCat]] = 0
    } else {
      logPList4Viol[[dataCat]] = logPList[[dataCat]] 
    }
  }
  # png(file.path(DESeq2Dir, 
  #               paste(workDir, "_DESeq2_", threshText, "_Violins_", DEdesc, ".png", sep="")), 
  #    type='cairo', width=16, height=10, units="in", res=150)
  # # had to use type='cairo'
  # # https://stat.ethz.ch/R-manual/R-devel/library/grDevices/html/png.html
  # # getOption("bitmapType")
  pdf(file.path(DESeq2Dir, 
                paste(workDir, "_DESeq2_", threshText, "_Violins_", DEdesc, ".pdf", sep="")), 
      width=16, height=10)
  par(mfrow=c(1,1), oma=c(5,0,2,0))
  colVioplot2(datas=logPList4Viol, names=altlogPListLabel, ylim=c(0,maxLogFC),
          ylab="log fold change of signicantly differentially expressed genes", las=2,
          col=rep(txparentCols4Dists, 4))
  title(main=paste(DEdesc, genomeType, "log2 fold expression distibutions"))
  stripchart(logPList,  method="jitter", jitter=0.35, pch=20, cex=1, vertical=TRUE,
             # panel.first = grid(),
             col=txparentCols4Strip, add=TRUE)  
  colVioplot2(datas=logPList4Viol, col=txparentCols4Dists, add=TRUE)
  abline(v=2.5, lwd=1, lty="dashed", col="black")
  abline(v=4.5, lwd=1, col="black")
  abline(v=6.5, lwd=1, lty="dashed", col="black")
  dev.off()

  ####
  # Tests of difference
  diffTests[[DEdesc]] = list()
  diffTests[[DEdesc]][["wilcoxtest.pvals"]] = c()
  
  # Log & test data
  x = list()
  for(itemA in 1:length(logPList)) {
    for(itemB in 1:length(logPList)) {
      cat(names(logPList)[itemA], "-v-", names(logPList)[itemB], "; ", sep="")
      
      wval = try(wilcox.test(logPList[[itemA]], logPList[[itemB]])$p.value)
      diffTests[[DEdesc]][["wilcoxtest.pvals"]] = c(diffTests[[DEdesc]][["wilcoxtest.pvals"]], ifelse(class(wval)=="try-error",1,wval))
    } 
  }
  cat("\n", sep="")
  
  diffTests[[DEdesc]][["wilcoxtest.pvals"]] = matrix(data=diffTests[[DEdesc]][["wilcoxtest.pvals"]], 
                                                     nrow=length(logPList), ncol=length(logPList), byrow=TRUE)
  colnames(diffTests[[DEdesc]][["wilcoxtest.pvals"]]) = names(logPList)
  rownames(diffTests[[DEdesc]][["wilcoxtest.pvals"]]) = names(logPList)
  
  wilcoxTestFile = file.path(DESeq2Dir, 
                             paste(workDir, "_DESeq2_", threshText, "_wilcoxTest_", DEdesc, ".tsv", sep=""))
  
  write.table(diffTests[[DEdesc]][["wilcoxtest.pvals"]] , wilcoxTestFile, 
              sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
}  
