###########################################################
###########################################################
#
# 171204
# Comparison of CTCFvsATAC and other features
#
# 180129
# Regions
#
# 180130
# 25Mb bins
#
###########################################################
###########################################################

CTCFprojectDir="<CTCFprojDirName>"
ATACprojectDir="<ATACprojDirName>"
HiCprojectDir="<HiCprojDirName>"
RNAseqProjectDir="<RNAseqProjDirName>"

outputDir = file.path(CTCFprojectDir, "CTCFvsATACvsFeatures")
setwd(outputDir)
getwd()

#save.image()
#load(file = ".RData") # From .RData


#########################
# Code 
require(plotrix)


#########################
# Variables 

dataTypes=c("CTCF", "ATAC")
samples=c( "WT", "Del", "InvDxz4" )
alleles=c("Ambig", "Xi", "Xa")

peakCols =  list()
peakCols[["WT"]] = "blue"
peakCols[["Del"]] = "red"
peakCols[["InvDxz4"]] = "black"

binSize=500000

firrebin = (50563120+(abs(50563120-50635321)/2)) / binSize
dxz4bin = (75764754+(abs(75637518-75764754)/2)) / binSize
xistbin = (103460373+(abs(103460373-103483233)/2)) / binSize

regionBoundaries = c((50563120+(abs(50563120-50635321)/2)), # Firre midpoint
                    (75764754+(abs(75637518-75764754)/2)), # Dxz4 midpoint
                    (103460373+(abs(103460373-103483233)/2))) # Xist
names(regionBoundaries) = c("Firre", "Dxz4", "Xist")

regional25MbBoundaries = seq(0, 175 * 10^6, 25 * 10^6)
names(regional25MbBoundaries) = as.character(seq(1, length(regional25MbBoundaries)))

chrApeakTotals = list()
chrApeakTotals[["ATAC"]] = c(160042, 342409, 80489)
names(chrApeakTotals[["ATAC"]]) = samples
chrApeakTotals[["CTCF"]] = c(91044, 75191, 63613)
names(chrApeakTotals[["CTCF"]]) = samples

chrAcovPeakTotals = list()
chrAcovPeakTotals[["ATAC"]] = c(138221, 180343,	69292)
names(chrAcovPeakTotals[["ATAC"]]) = samples
chrAcovPeakTotals[["CTCF"]] = c(82694, 66652, 57827)
names(chrAcovPeakTotals[["CTCF"]]) = samples

chrXcovPeakTotals = list()
chrXcovPeakTotals[["ATAC"]] = c(3054, 3012, 1902)
names(chrXcovPeakTotals[["ATAC"]]) = samples
chrXcovPeakTotals[["CTCF"]] = c(2169, 2012, 1528)
names(chrXcovPeakTotals[["CTCF"]]) = samples


#########################
# Input files

BEDdir = list()
PEAKfiles=list()

BEDdir[["CTCF"]] = file.path(CTCFprojectDir, "diploidPeakSegregation/peakAssignment/BEDfiles")
PEAKfiles[["CTCF"]] = list()
PEAKfiles[["CTCF"]][["WT"]]=list()
PEAKfiles[["CTCF"]][["WT"]][["Ambig"]]="patskiWTvsDelvsInv.allelicCTCF.WT.ambiguous.minCov=5.adjPropThresh=0.7.chrX.bed.gz"
PEAKfiles[["CTCF"]][["WT"]][["Xi"]]="patskiWTvsDelvsInv.allelicCTCF.WT.B6-specific.minCov=5.adjPropThresh=0.7.chrX.bed.gz"
PEAKfiles[["CTCF"]][["WT"]][["Xa"]]="patskiWTvsDelvsInv.allelicCTCF.WT.Sp-specific.minCov=5.adjPropThresh=0.7.chrX.bed.gz"
PEAKfiles[["CTCF"]][["Del"]]=list()
PEAKfiles[["CTCF"]][["Del"]][["Ambig"]]="patskiWTvsDelvsInv.allelicCTCF.Del.ambiguous.minCov=5.adjPropThresh=0.7.chrX.bed.gz"
PEAKfiles[["CTCF"]][["Del"]][["Xi"]]="patskiWTvsDelvsInv.allelicCTCF.Del.B6-specific.minCov=5.adjPropThresh=0.7.chrX.bed.gz"
PEAKfiles[["CTCF"]][["Del"]][["Xa"]]="patskiWTvsDelvsInv.allelicCTCF.Del.Sp-specific.minCov=5.adjPropThresh=0.7.chrX.bed.gz"
PEAKfiles[["CTCF"]][["InvDxz4"]]=list()
PEAKfiles[["CTCF"]][["InvDxz4"]][["Ambig"]]="patskiWTvsDelvsInv.allelicCTCF.InvDxz4.ambiguous.minCov=5.adjPropThresh=0.7.chrX.bed.gz"
PEAKfiles[["CTCF"]][["InvDxz4"]][["Xi"]]="patskiWTvsDelvsInv.allelicCTCF.InvDxz4.B6-specific.minCov=5.adjPropThresh=0.7.chrX.bed.gz"
PEAKfiles[["CTCF"]][["InvDxz4"]][["Xa"]]="patskiWTvsDelvsInv.allelicCTCF.InvDxz4.Sp-specific.minCov=5.adjPropThresh=0.7.chrX.bed.gz"

BEDdir[["ATAC"]] = file.path(ATACprojectDir, "diploidPeakSegregation/peakAssignment/BEDfiles")
PEAKfiles[["ATAC"]] = list()
PEAKfiles[["ATAC"]][["WT"]]=list()
PEAKfiles[["ATAC"]][["WT"]][["Ambig"]]="patskiWTvsDelvsInvDxz4.allelicATAC.WT.ambiguous.minCov=5.adjPropThresh=0.7.chrX.bed.gz"
PEAKfiles[["ATAC"]][["WT"]][["Xi"]]="patskiWTvsDelvsInvDxz4.allelicATAC.WT.B6-specific.minCov=5.adjPropThresh=0.7.chrX.bed.gz"
PEAKfiles[["ATAC"]][["WT"]][["Xa"]]="patskiWTvsDelvsInvDxz4.allelicATAC.WT.Sp-specific.minCov=5.adjPropThresh=0.7.chrX.bed.gz"
PEAKfiles[["ATAC"]][["Del"]]=list()
PEAKfiles[["ATAC"]][["Del"]][["Ambig"]]="patskiWTvsDelvsInvDxz4.allelicATAC.Del.ambiguous.minCov=5.adjPropThresh=0.7.chrX.bed.gz"
PEAKfiles[["ATAC"]][["Del"]][["Xi"]]="patskiWTvsDelvsInvDxz4.allelicATAC.Del.B6-specific.minCov=5.adjPropThresh=0.7.chrX.bed.gz"
PEAKfiles[["ATAC"]][["Del"]][["Xa"]]="patskiWTvsDelvsInvDxz4.allelicATAC.Del.Sp-specific.minCov=5.adjPropThresh=0.7.chrX.bed.gz"
PEAKfiles[["ATAC"]][["InvDxz4"]]=list()
PEAKfiles[["ATAC"]][["InvDxz4"]][["Ambig"]]="patskiWTvsDelvsInvDxz4.allelicATAC.InvDxz4.ambiguous.minCov=5.adjPropThresh=0.7.chrX.bed.gz"
PEAKfiles[["ATAC"]][["InvDxz4"]][["Xi"]]="patskiWTvsDelvsInvDxz4.allelicATAC.InvDxz4.B6-specific.minCov=5.adjPropThresh=0.7.chrX.bed.gz"
PEAKfiles[["ATAC"]][["InvDxz4"]][["Xa"]]="patskiWTvsDelvsInvDxz4.allelicATAC.InvDxz4.Sp-specific.minCov=5.adjPropThresh=0.7.chrX.bed.gz"

DSCOREdir = list()
DSCOREfiles = list()

DSCOREdir[["CTCF"]] = file.path(CTCFprojectDir, "diploidPeakSegregation/allelicRatios/adjAllelicPropBedGraphs")
DSCOREfiles[["CTCF"]] = list()
DSCOREfiles[["CTCF"]][["WT"]]="patskiWTvsDelvsInv.allelicCTCF.WT.dScore.minCov=5.chrX.bedgraph.gz"
DSCOREfiles[["CTCF"]][["Del"]]="patskiWTvsDelvsInv.allelicCTCF.Del.dScore.minCov=5.chrX.bedgraph.gz"
DSCOREfiles[["CTCF"]][["InvDxz4"]]="patskiWTvsDelvsInv.allelicCTCF.InvDxz4.dScore.minCov=5.chrX.bedgraph.gz"

DSCOREdir[["ATAC"]] = file.path(ATACprojectDir, "diploidPeakSegregation/allelicRatios/adjAllelicPropBedGraphs")
DSCOREfiles[["ATAC"]]=list()
DSCOREfiles[["ATAC"]][["WT"]]="patskiWTvsDelvsInvDxz4.allelicATAC.WT.dScore.minCov=5.chrX.bedgraph.gz"
DSCOREfiles[["ATAC"]][["Del"]]="patskiWTvsDelvsInvDxz4.allelicATAC.Del.dScore.minCov=5.chrX.bedgraph.gz"
DSCOREfiles[["ATAC"]][["InvDxz4"]]="patskiWTvsDelvsInvDxz4.allelicATAC.InvDxz4.dScore.minCov=5.chrX.bedgraph.gz"

# PC scores
PCSCOREdir =  file.path(HiCprojectDir, "compartments.vSNPs.CIGARsed.spearman/scoreTables")
PCSCOREfiles = list()
PCSCOREfiles[["WT"]] = list()
PCSCOREfiles[["Del"]] = list()
PCSCOREfiles[["InvDxz4"]] = list()
PCSCOREfiles[["WT"]][["Xi"]]="compartments.vSNPs.CIGARsed.500000.PCA.scores.WT_ref.chrX.tsv"
PCSCOREfiles[["Del"]][["Xi"]]="compartments.vSNPs.CIGARsed.500000.PCA.scores.Del_ref.chrX.tsv"
PCSCOREfiles[["InvDxz4"]][["Xi"]]="compartments.vSNPs.CIGARsed.500000.PCA.scores.invDxz4_ref.chrX.tsv"
PCSCOREfiles[["WT"]][["Xa"]]="compartments.vSNPs.CIGARsed.500000.PCA.scores.WT_alt.chrX.tsv"
PCSCOREfiles[["Del"]][["Xa"]]="compartments.vSNPs.CIGARsed.500000.PCA.scores.Del_alt.chrX.tsv"
PCSCOREfiles[["InvDxz4"]][["Xa"]]="compartments.vSNPs.CIGARsed.500000.PCA.scores.invDxz4_alt.chrX.tsv"

# Expression data
RNASEQdir = file.path(RNAseqProjectDir, "expressionAnalysis")

# DESeq2 results tables
RNASEQaltdir = file.path(RNAseqProjectDir, "expressionAnalysis.correctedAltToRefRatio")

DEresultDirs = list()
DEresultDirs[["WT-vs-Del"]] = file.path(RNASEQaltdir, "DESeq2/WT-vs-Del1a")
DEresultDirs[["WT-vs-InvDxz4"]] = file.path(RNASEQdir, "DESeq2/WT-vs-InvDxz4")
DEresultDirs[["newWT-vs-InvDxz4"]] = file.path(RNASEQdir, "DESeq2/newWT-vs-InvDxz4")

DEresultFiles = list()
DEresultFiles[["WT-vs-Del"]] = list ()
DEresultFiles[["WT-vs-Del"]][["Xi"]] = file.path(DEresultDirs[["WT-vs-Del"]], 
                                                 "Patski.AllData_DESeq2genes_logP=1.3_logFC=0.5_B6.WT-vs-Del1.chrX.tsv")
DEresultFiles[["WT-vs-Del"]][["Xa"]] = file.path(DEresultDirs[["WT-vs-Del"]], 
                                                 "Patski.AllData_DESeq2genes_logP=1.3_logFC=0.5_sp.WT-vs-Del1.chrX.tsv")
DEresultFiles[["WT-vs-InvDxz4"]] = list ()
DEresultFiles[["WT-vs-InvDxz4"]][["Xi"]] = file.path(DEresultDirs[["WT-vs-InvDxz4"]], 
                                                     "Patski.combo.AllData_DESeq2genes_logP=1.3_logFC=0.5_B6.WT-vs-InvDxz4.chrX.tsv")
DEresultFiles[["WT-vs-InvDxz4"]][["Xa"]] = file.path(DEresultDirs[["WT-vs-InvDxz4"]], 
                                                     "Patski.combo.AllData_DESeq2genes_logP=1.3_logFC=0.5_sp.WT-vs-InvDxz4.chrX.tsv")

DEresultFiles[["newWT-vs-InvDxz4"]] = list ()
DEresultFiles[["newWT-vs-InvDxz4"]][["Xi"]] = file.path(DEresultDirs[["newWT-vs-InvDxz4"]], 
                                                        "Patski.combo.AllData_DESeq2genes_logP=1.3_logFC=0.5_B6.WT-vs-InvDxz4.chrX.tsv")
DEresultFiles[["newWT-vs-InvDxz4"]][["Xa"]] = file.path(DEresultDirs[["newWT-vs-InvDxz4"]], 
                                                        "Patski.combo.AllData_DESeq2genes_logP=1.3_logFC=0.5_sp.WT-vs-InvDxz4.chrX.tsv")


#########################
# Load data

# Load peak data
peakData = list()
for( theDataType in dataTypes) {
  peakData[[theDataType]] = list()
  for(theSample in names(PEAKfiles[[theDataType]])) {
    peakData[[theDataType]][[theSample]] = list()
    for(theAllele in names(PEAKfiles[[theDataType]][[theSample]])) {
      peakData[[theDataType]][[theSample]][[theAllele]] = as.data.frame(read.table(gzfile(file.path(BEDdir[[theDataType]], PEAKfiles[[theDataType]][[theSample]][[theAllele]])), header=FALSE, skip=1, sep="\t")) #, fill = TRUE)
      colnames(peakData[[theDataType]][[theSample]][[theAllele]]) = c("chr", "start", "end")
    }
  }
}


# Load dscore data
dscoreData = list()
for (theDataType in dataTypes) {
  dscoreData[[theDataType]] = list()
  for(theSample in names(DSCOREfiles[[theDataType]] )) {
    dscoreData[[theSample]][[theDataType]]  = as.data.frame(read.table(gzfile(file.path(DSCOREdir[[theDataType]], DSCOREfiles[[theDataType]][[theSample]])), header=FALSE, skip=1, sep="\t")) #, fill = TRUE)
    colnames(dscoreData[[theSample]][[theDataType]] ) = c("chr", "start", "end", "score")
  }
}


# Load PC score data
PCscoreData = list()
for(theSample in names(PCSCOREfiles)) {
  PCscoreData[[theSample]] = list()
  for(theAllele in names(PCSCOREfiles[[theSample]])) { 
    PCscoreData[[theSample]][[theAllele]] = as.data.frame(read.table(gzfile(file.path(PCSCOREdir, PCSCOREfiles[[theSample]][[theAllele]])), header=TRUE, skip=0, sep="\t")) #, fill = TRUE)
  }
}
PCscoreData.orig = PCscoreData


# Load expression tables
tableDir = file.path(RNASEQdir, "expressionTables.withStrand")
expressionTables.chrX = list()
for (expDataType in c("Counts", "TPM", "RPKM")) {
  for (genome in c("", "ref", "alt")) {
    if ( genome=="" && expDataType == "Counts") {
      expDataType = "counts"
    }
    tableName = paste(genome, expDataType, sep = "")
    print(tableName)
    bigOldFilename = file.path(tableDir, paste(tableName, ".chrX.tsv.gz", sep = ""))
    expressionTables.chrX[[tableName]] = read.table(gzfile(bigOldFilename), sep = "\t", header=TRUE)
  }
}


# Load DE result tables
DEtables.chrX = list()
for (comp in names(DEresultFiles)) {
  DEtables.chrX[[comp]] = list()
  for (allele in names(DEresultFiles[[comp]])) {
    tableName = DEresultFiles[[comp]][[allele]]
    print(tableName)
    DEtables.chrX[[comp]][[allele]] = read.table(tableName, sep = "\t", header=TRUE)
  }
}


#########################
# Bin data

binnedPeakCounts = list()
for (theDataType in dataTypes) {
  binnedPeakCounts[[theDataType]] = list()
  for(theSample in names(peakData[[theDataType]])) {
    binnedPeakCounts[[theDataType]][[theSample]] = list()
    for(theAllele in names(peakData[[theDataType]][[theSample]])) {
      thisPeakData = peakData[[theDataType]][[theSample]][[theAllele]]
      midPoints = as.numeric(as.character(thisPeakData$start)) + (as.numeric(as.character(thisPeakData$end)) - as.numeric(as.character(thisPeakData$start)))/2
      midPoints = midPoints/binSize
      histData = hist(midPoints, breaks=seq(0,(180 * (1000000/binSize)) + 1,1), plot=FALSE)
      binnedPeakCounts[[theDataType]][[theSample]][[theAllele]] = as.numeric(as.character(histData$counts))
    }
  }
}


#########################
# 180129
# Bin data into regions

regionalBinnedPeakCounts = list()
for (theDataType in dataTypes) {
  regionalBinnedPeakCounts[[theDataType]] = list()
  for(theSample in names(peakData[[theDataType]])) {
    regionalBinnedPeakCounts[[theDataType]][[theSample]] = list()
    for(theAllele in names(peakData[[theDataType]][[theSample]])) {
      thisPeakData = peakData[[theDataType]][[theSample]][[theAllele]]
      midPoints = as.numeric(as.character(thisPeakData$start)) + (as.numeric(as.character(thisPeakData$end)) - as.numeric(as.character(thisPeakData$start)))/2
      histData = hist(midPoints, breaks=c(0, regionBoundaries, 180 * 1000000), plot=FALSE)
      regionalBinnedPeakCounts[[theDataType]][[theSample]][[theAllele]] = as.numeric(as.character(histData$counts))
    }
  }
}


#########################
# 180130
# Bin data into 25Mb regions

regional25MbBinnedPeakCounts = list()
for (theDataType in dataTypes) {
  regional25MbBinnedPeakCounts[[theDataType]] = list()
  for(theSample in names(peakData[[theDataType]])) {
    regional25MbBinnedPeakCounts[[theDataType]][[theSample]] = list()
    for(theAllele in names(peakData[[theDataType]][[theSample]])) {
      thisPeakData = peakData[[theDataType]][[theSample]][[theAllele]]
      midPoints = as.numeric(as.character(thisPeakData$start)) + (as.numeric(as.character(thisPeakData$end)) - as.numeric(as.character(thisPeakData$start)))/2
      histData = hist(midPoints, breaks=regional25MbBoundaries, plot=FALSE)
      regional25MbBinnedPeakCounts[[theDataType]][[theSample]][[theAllele]] = as.numeric(as.character(histData$counts))
    }
  }
}



#########################
# Pooled Xi-associated data

for (theDataType in dataTypes) {
  for(theSample in names(peakData[[theDataType]])) {
    thisPeakData.ambig = binnedPeakCounts[[theDataType]][[theSample]][["Ambig"]]
    thisPeakData.xi = binnedPeakCounts[[theDataType]][[theSample]][["Xi"]]
    binnedPeakCounts[[theDataType]][[theSample]][["Xi-assoc"]] = thisPeakData.ambig + thisPeakData.xi
  
    thisPeakData.ambig = regionalBinnedPeakCounts[[theDataType]][[theSample]][["Ambig"]]
    thisPeakData.xi = regionalBinnedPeakCounts[[theDataType]][[theSample]][["Xi"]]
    regionalBinnedPeakCounts[[theDataType]][[theSample]][["Xi-assoc"]] = thisPeakData.ambig + thisPeakData.xi
    
    thisPeakData.ambig = regional25MbBinnedPeakCounts[[theDataType]][[theSample]][["Ambig"]]
    thisPeakData.xi = regional25MbBinnedPeakCounts[[theDataType]][[theSample]][["Xi"]]
    regional25MbBinnedPeakCounts[[theDataType]][[theSample]][["Xi-assoc"]] = thisPeakData.ambig + thisPeakData.xi
  }
}


#########################
# Pooled Xa-associated data

for (theDataType in dataTypes) {
  for(theSample in names(peakData[[theDataType]])) {
    thisPeakData.ambig = binnedPeakCounts[[theDataType]][[theSample]][["Ambig"]]
    thisPeakData.xa = binnedPeakCounts[[theDataType]][[theSample]][["Xa"]]
    binnedPeakCounts[[theDataType]][[theSample]][["Xa-assoc"]] = thisPeakData.ambig + thisPeakData.xa
  
    thisPeakData.ambig = regionalBinnedPeakCounts[[theDataType]][[theSample]][["Ambig"]]
    thisPeakData.xa = regionalBinnedPeakCounts[[theDataType]][[theSample]][["Xa"]]
    regionalBinnedPeakCounts[[theDataType]][[theSample]][["Xa-assoc"]] = thisPeakData.ambig + thisPeakData.xa
    
    thisPeakData.ambig = regional25MbBinnedPeakCounts[[theDataType]][[theSample]][["Ambig"]]
    thisPeakData.xa = regional25MbBinnedPeakCounts[[theDataType]][[theSample]][["Xa"]]
    regional25MbBinnedPeakCounts[[theDataType]][[theSample]][["Xa-assoc"]] = thisPeakData.ambig + thisPeakData.xa
  }
}


#########################
# Calculate Xa bin count scaling factors 

XaPeakTotals = list()
XaPeakRatios = list()

for (theDataType in dataTypes) {
  XaPeakTotals[[theDataType]] = c()
  for(theSample in names(peakData[[theDataType]])) {
    theAllele = "Xa"
    XaPeakTotals[[theDataType]][theSample] = nrow(peakData[[theDataType]][[theSample]][[theAllele]])
  }
  minXaPeakIdx = which.min(XaPeakTotals[[theDataType]])
  XaPeakRatios[[theDataType]] = XaPeakTotals[[theDataType]][minXaPeakIdx]/XaPeakTotals[[theDataType]]
}



######
# Calculate scaling factors 
# using autosomal diploid count 

chrAPeakRatios = list()
for (theDataType in dataTypes) {
  minPeakIdx = which.min(chrApeakTotals[[theDataType]])
  chrAPeakRatios[[theDataType]] = chrApeakTotals[[theDataType]][minPeakIdx]/chrApeakTotals[[theDataType]]
}



######
# Calculate scaling factors 
# using autosomal diploid count 

chrXcovPeakRatios = list()
for (theDataType in dataTypes) {
  minPeakIdx = which.min(chrXcovPeakTotals[[theDataType]])
  chrXcovPeakRatios[[theDataType]] = chrXcovPeakTotals[[theDataType]][minPeakIdx]/chrXcovPeakTotals[[theDataType]]
}



#####
# Calculate scaling factors 
# using autosomal diploid count 

chrAcovPeakRatios = list()
for (theDataType in dataTypes) {
  minPeakIdx = which.min(chrAcovPeakTotals[[theDataType]])
  chrAcovPeakRatios[[theDataType]] = chrAcovPeakTotals[[theDataType]][minPeakIdx]/chrAcovPeakTotals[[theDataType]]
}



#####
# Use chrAcovPeakRatios:

# PeakRatioScalingFactors = XaPeakRatios; scalingText = "scaled by relative ratios of total Xa peaks"
# PeakRatioScalingFactors = chrXcovPeakRatios; scalingText = "scaled by relative ratios of total diploid X peaks"
PeakRatioScalingFactors = chrAcovPeakRatios; scalingText = "scaled by relative ratios of total diploid autosomal peaks"



#########################
# Calculate scaled counts

binnedPeakCounts.scaled = list()
for (theDataType in dataTypes) {
  binnedPeakCounts.scaled[[theDataType]] = list()
  for(theSample in names(binnedPeakCounts[[theDataType]])) {
    binnedPeakCounts.scaled[[theDataType]][[theSample]] = list()
    cat(theDataType, theSample, PeakRatioScalingFactors[[theDataType]][theSample], "\n", sep=" ")
    for(theAllele in names(binnedPeakCounts[[theDataType]][[theSample]])) {
      binnedPeakCounts.scaled[[theDataType]][[theSample]][[theAllele]] = binnedPeakCounts[[theDataType]][[theSample]][[theAllele]] * PeakRatioScalingFactors[[theDataType]][theSample]
    }
  }
}


regionalBinnedPeakCounts.scaled = list()
for (theDataType in dataTypes) {
  regionalBinnedPeakCounts.scaled[[theDataType]] = list()
  for(theSample in names(binnedPeakCounts[[theDataType]])) {
    regionalBinnedPeakCounts.scaled[[theDataType]][[theSample]] = list()
    cat(theDataType, theSample, PeakRatioScalingFactors[[theDataType]][theSample], "\n", sep=" ")
    for(theAllele in names(binnedPeakCounts[[theDataType]][[theSample]])) {
      regionalBinnedPeakCounts.scaled[[theDataType]][[theSample]][[theAllele]] = regionalBinnedPeakCounts[[theDataType]][[theSample]][[theAllele]] * PeakRatioScalingFactors[[theDataType]][theSample]
    }
  }
}


regional25MbBinnedPeakCounts.scaled = list()
for (theDataType in dataTypes) {
  regional25MbBinnedPeakCounts.scaled[[theDataType]] = list()
  for(theSample in names(binnedPeakCounts[[theDataType]])) {
    regional25MbBinnedPeakCounts.scaled[[theDataType]][[theSample]] = list()
    cat(theDataType, theSample, PeakRatioScalingFactors[[theDataType]][theSample], "\n", sep=" ")
    for(theAllele in names(binnedPeakCounts[[theDataType]][[theSample]])) {
      regional25MbBinnedPeakCounts.scaled[[theDataType]][[theSample]][[theAllele]] = regional25MbBinnedPeakCounts[[theDataType]][[theSample]][[theAllele]] * PeakRatioScalingFactors[[theDataType]][theSample]
    }
  }
}


#########################
# Get total counts 

totalCounts = list()
totalCounts.scaled = list()
for (theDataType in dataTypes) {
  totalCounts[[theDataType]] = list()
  totalCounts.scaled[[theDataType]] = list()
  for(theSample in names(binnedPeakCounts[[theDataType]])) {
    totalCounts[[theDataType]][[theSample]]  = list()
    totalCounts.scaled[[theDataType]][[theSample]]  = list()
    for(theAllele in names(binnedPeakCounts[[theDataType]][[theSample]])) {
      totalCounts[[theDataType]][[theSample]][[theAllele]] = sum(binnedPeakCounts[[theDataType]][[theSample]][[theAllele]])
      totalCounts.scaled[[theDataType]][[theSample]][[theAllele]]  = sum(binnedPeakCounts.scaled[[theDataType]][[theSample]][[theAllele]])
    }
  }
}

# Just 3 original alleles
totalCounts.orig = list()
totalCounts.orig.scaled = list()
for (theDataType in dataTypes) {
  totalCounts.orig[[theDataType]] = list()
  totalCounts.orig.scaled[[theDataType]] = list()
  for(theSample in names(binnedPeakCounts[[theDataType]])) {
    totalCounts.orig[[theDataType]][[theSample]]  = list()
    totalCounts.orig.scaled[[theDataType]][[theSample]]  = list()
    for(theAllele in alleles) {
      totalCounts.orig[[theDataType]][[theSample]][[theAllele]] = sum(binnedPeakCounts[[theDataType]][[theSample]][[theAllele]])
      totalCounts.orig.scaled[[theDataType]][[theSample]][[theAllele]]  = sum(binnedPeakCounts.scaled[[theDataType]][[theSample]][[theAllele]])
    }
  }
}



#########################
# Calculate max Y-axis-values

maxY = regionalMaxY = regional25MbMaxY = list()

for (theDataType in dataTypes) {
  maxY[[theDataType]] = -Inf
  regionalMaxY[[theDataType]] = -Inf
  regional25MbMaxY[[theDataType]] = -Inf
  for(theSample in names(binnedPeakCounts.scaled[[theDataType]])) {
    for(theAllele in names(binnedPeakCounts.scaled[[theDataType]][[theSample]])) {
      maxVal = max(binnedPeakCounts.scaled[[theDataType]][[theSample]][[theAllele]], na.rm=TRUE)
      if(maxY[[theDataType]] < maxVal) {
        maxY[[theDataType]] = maxVal
      }    
      maxVal = max(regionalBinnedPeakCounts.scaled[[theDataType]][[theSample]][[theAllele]], na.rm=TRUE)
      if(regionalMaxY[[theDataType]] < maxVal) {
        regionalMaxY[[theDataType]] = maxVal
      }    
      maxVal = max(regional25MbBinnedPeakCounts.scaled[[theDataType]][[theSample]][[theAllele]], na.rm=TRUE)
      if(regional25MbMaxY[[theDataType]] < maxVal) {
        regional25MbMaxY[[theDataType]] = maxVal
      }    
    }
  }
}


#########################
# Calculate max PC values

PCscoreData = PCscoreData.orig

# First correct some known incorrect values
PCscoreData[["Del"]][["Xi"]][,"PC2"] = PCscoreData[["Del"]][["Xi"]][,"PC2"] * -1
PCscoreData[["InvDxz4"]][["Xi"]][,"PC2"] = PCscoreData[["InvDxz4"]][["Xi"]][,"PC2"] * -1
PCscoreData[["WT"]][["Xi"]][,"PC3"] = PCscoreData[["WT"]][["Xi"]][,"PC3"] * -1
PCscoreData[["InvDxz4"]][["Xi"]][,"PC3"] = PCscoreData[["InvDxz4"]][["Xi"]][,"PC3"] * -1
PCscoreData[["Del"]][["Xi"]][,"PC4"] = PCscoreData[["Del"]][["Xi"]][,"PC4"] * -1

PCylims = list()
for(PCnum in paste("PC", 1:5, sep="")) {
  PCylims[[PCnum]] = list()
  for(theAllele in names(PCscoreData[[theSample]])) {
    PCylims[[PCnum]][[theAllele]] = c(Inf, -Inf)
    for(theSample in names(PCscoreData)) {
      rangeVal = range(PCscoreData[[theSample]][[theAllele]][,PCnum], na.rm=TRUE)
      if(PCylims[[PCnum]][[theAllele]][1] > rangeVal[1]) {
        PCylims[[PCnum]][[theAllele]][1] = rangeVal[1]
      }    
      if(PCylims[[PCnum]][[theAllele]][2] < rangeVal[2]) {
        PCylims[[PCnum]][[theAllele]][2] = rangeVal[2]
      }    
    }
  }
}


#########################
# Gene density

chrXTSS = expressionTables.chrX[["counts"]]$start
chrXTSS[expressionTables.chrX[["counts"]]$strand == "-"] = expressionTables.chrX[["counts"]]$stop[expressionTables.chrX[["counts"]]$strand == "-"]
chrXTSS = chrXTSS/binSize
histData = hist(chrXTSS, breaks=seq(0,(180 * (1000000/binSize)) + 1,1), plot=FALSE)
binnedGeneCounts = as.numeric(as.character(histData$counts))









###########################################################
# ATAC vs CTCF density Plots


#########################
# 180130
# 25Mb density plots

densityPlotDir = file.path(outputDir, "regional25MbDensityPlots")
if ( ! file.exists( densityPlotDir ) ) {
  if ( ! dir.create( densityPlotDir, showWarnings=FALSE, recursive=TRUE )  ) densityPlotDir = outputDir
}

# regionMidpoints = c(0, regionBoundaries, 180 * 1000000)[1:4]+ (c(0, regionBoundaries, 180 * 1000000)[2:5] - c(0, regionBoundaries, 180 * 1000000)[1:4])/2

for (theDataType in dataTypes) {
  theXlims = c(0, 175 * (1000000))
  theYlims = c(0, regional25MbMaxY[[theDataType]])
  
  for(theAllele in c(alleles, "Xi-assoc", "Xa-assoc")) {
    
    dataDesc = paste("regional25MbBinnedPeakCounts_scaled_", theDataType, "_", theAllele, sep="")
    pdf(file.path(densityPlotDir, paste(dataDesc,".pdf",sep="")), width=20, height=5)
    
    # Layout
    layout(matrix(c(1,1,1,1,2,3,4,5), nrow=4, ncol=2, byrow=FALSE), widths=c(1,0.2))
    # layout.show(4)
    par(mar=c(2,4,6,4), oma=c(3,0,2,0))
    
    # Plot WT counts
    # Plot an empty graph and legend to get the size of the legend
    # https://stackoverflow.com/questions/8929663/r-legend-placement-in-a-plot
    plot(regional25MbBoundaries,
         c(regional25MbBinnedPeakCounts.scaled[[theDataType]][["WT"]][[theAllele]], 
           rev(regional25MbBinnedPeakCounts.scaled[[theDataType]][["WT"]][[theAllele]])[1]), 
         type="n", xaxt="n", yaxt="n",  xlab=NA, ylab=NA, ylim=theYlims, xlim=theXlims)
    # # Custom ylim. Add the height of legend to upper bound of the range
    # ymax = max(c(daYlim, ymax), na.rm=TRUE)
    # my.range <- c(0, ymax)
    # ##my.range[2] <- 1.04*(my.range[2]+my.legend.size$rect$h)
    
    # Draw the WT plot with custom ylim
    par(new=TRUE)
    plot(regional25MbBoundaries,
         c(regional25MbBinnedPeakCounts.scaled[[theDataType]][["WT"]][[theAllele]], 
           rev(regional25MbBinnedPeakCounts.scaled[[theDataType]][["WT"]][[theAllele]])[1]), 
         type="s", lty="solid", lwd=2, col=peakCols[["WT"]], axes=F, main=NA, xlab=NA, ylab=NA, ylim=theYlims, xlim=theXlims)
    axis(side=2, line=0)
    mtext(side=2, line=2.5, text="Normalized counts")
    title(main=paste("Counts (", scalingText, ")", sep=""), line=3)
    
    # Plot other samples
    for (theSample in samples) {
      if (theSample != "WT") {
        lines(regional25MbBoundaries,
              c(regional25MbBinnedPeakCounts.scaled[[theDataType]][[theSample]][[theAllele]], 
                rev(regional25MbBinnedPeakCounts.scaled[[theDataType]][[theSample]][[theAllele]])[1]),
              type="s", lty="solid", lwd=2, col=peakCols[[theSample]])
        # abline(h=expectedCounts[[theSample]], lty="dashed", col=peakCols[[theSample]])
      }
    }
    
    # LOI
    abline(v=regionBoundaries["Firre"], lty="dashed", col="black")
    abline(v=regionBoundaries["Dxz4"], lty="dashed", col="black")
    abline(v=regionBoundaries["Xist"], lty="dashed", col="black")
    axis(3, at=regionBoundaries, labels=names(regionBoundaries), tick = FALSE)
    axis(1, at=regional25MbBoundaries, labels=regional25MbBoundaries/10^6)
    abline(v=regional25MbBoundaries, col="gray", lty="dotted", lw=0.75)    
    
    title(main=paste(dataDesc, sep=""), outer=TRUE)
    title(xlab="chrX (Mb)", outer=TRUE, line=1)
    
    # Legend
    par(mar=c(0,0,0,0))
    plot(seq(0,10), type="n", axes=FALSE, xlab="", ylab="")
    legend("top", 
           legend=c(
             paste("WT", sep=""),
             paste("Del", sep=""),
             paste("InvDxz4", sep="")
           ),
           lty=c(1,1,1), col=c(peakCols[["WT"]], peakCols[["Del"]], peakCols[["InvDxz4"]]), bty='n', cex=1)
    
    # add total counts table
    theTotCounts = unlist(totalCounts[[theDataType]])[grep(paste(theAllele, "$", sep=""), names(unlist(totalCounts[[theDataType]])) )]
    theTotCounts.scaled = unlist(totalCounts.scaled[[theDataType]])[grep(paste(theAllele, "$", sep=""), names(unlist(totalCounts.scaled[[theDataType]])) )]
    statsTable = cbind(theTotCounts, round(PeakRatioScalingFactors[[theDataType]], 2), round(theTotCounts.scaled))
    colnames(statsTable) = c("total", "Xa SF", "scaled")
    rownames(statsTable) = samples 
    
    par(mar=c(0,0,2,0))
    plot(seq(1,10), type="n", axes=FALSE, xlab="", ylab="")
    # https://stackoverflow.com/questions/18663159/conditional-coloring-of-cells-in-table
    # require(plotrix)
    addtable2plot("top", table=statsTable, bg="white", title="Total peaks before and after scaling",
                  bty="o", display.rownames=TRUE, display.colnames=TRUE,
                  hlines=TRUE, vlines=TRUE,
                  xjust=2, yjust=1, cex=1)
  
    dev.off()
  }
}




#########################
# 180130
# Relative 25Mb density plots

densityPlotDir = file.path(outputDir, "log2RatioRegional25MbDensityPlots")
if ( ! file.exists( densityPlotDir ) ) {
  if ( ! dir.create( densityPlotDir, showWarnings=FALSE, recursive=TRUE )  ) densityPlotDir = outputDir
}

# regionMidpoints = c(0, regionBoundaries, 180 * 1000000)[1:4]+ (c(0, regionBoundaries, 180 * 1000000)[2:5] - c(0, regionBoundaries, 180 * 1000000)[1:4])/2

for (theDataType in dataTypes) {
  theXlims = c(0, 175 * (1000000))
  theYlims = c(-2, +6)
  
  for(theAllele in c(alleles, "Xi-assoc", "Xa-assoc")) {
    
    dataDesc = paste("log2RatioRegional25MbBinnedPeakCounts_scaled_", theDataType, "_", theAllele, sep="")
    pdf(file.path(densityPlotDir, paste(dataDesc,".pdf",sep="")), width=20, height=5)
    
    # Layout
    layout(matrix(c(1,1,1,1,2,3,4,5), nrow=4, ncol=2, byrow=FALSE), widths=c(1,0.2))
    # layout.show(4)
    par(mar=c(2,4,6,4), oma=c(3,0,2,0))
    
    # Plot WT counts
    # Plot an empty graph and legend to get the size of the legend
    # https://stackoverflow.com/questions/8929663/r-legend-placement-in-a-plot
    plot(regional25MbBoundaries, rep(0, length(regional25MbBoundaries)), 
         type="n", xaxt="n", yaxt="n",  xlab=NA, ylab=NA, ylim=theYlims, xlim=theXlims)
    # # Custom ylim. Add the height of legend to upper bound of the range
    # ymax = max(c(daYlim, ymax), na.rm=TRUE)
    # my.range <- c(0, ymax)
    # ##my.range[2] <- 1.04*(my.range[2]+my.legend.size$rect$h)
    
    # Draw the WT plot with custom ylim
    par(new=TRUE)
    plot(regional25MbBoundaries, 
         log2(c(regional25MbBinnedPeakCounts.scaled[[theDataType]][["Del"]][[theAllele]], 
                rev(regional25MbBinnedPeakCounts.scaled[[theDataType]][["Del"]][[theAllele]])[1]) / 
                c(regional25MbBinnedPeakCounts.scaled[[theDataType]][["WT"]][[theAllele]], 
                  rev(regional25MbBinnedPeakCounts.scaled[[theDataType]][["WT"]][[theAllele]])[1])), 
         type="s", lty="solid", lwd=2, col=peakCols[["Del"]], axes=F, main=NA, xlab=NA, ylab=NA, ylim=theYlims, xlim=theXlims)
    axis(side=2, line=0)
    mtext(side=2, line=2.5, text="Log ratio")
    axis(1, at=regional25MbBoundaries, labels=regional25MbBoundaries/10^6)
    title(main=paste("Log2Ratio vs. WT (", scalingText, ")", sep=""), line=3)
    
    # Plot other samples
    lines(c(0, regionBoundaries, 180 * 1000000), 
          log2(c(regionalBinnedPeakCounts.scaled[[theDataType]][["InvDxz4"]][[theAllele]], 
                 rev(regionalBinnedPeakCounts.scaled[[theDataType]][["InvDxz4"]][[theAllele]])[1]) / 
                 c(regionalBinnedPeakCounts.scaled[[theDataType]][["WT"]][[theAllele]], 
                   rev(regionalBinnedPeakCounts.scaled[[theDataType]][["WT"]][[theAllele]])[1])), 
          type="s", lty="solid", lwd=2, col=peakCols[["InvDxz4"]])
    
    # LOI
    abline(v=regionBoundaries["Firre"], lty="dashed", col="black")
    abline(v=regionBoundaries["Dxz4"], lty="dashed", col="black")
    abline(v=regionBoundaries["Xist"], lty="dashed", col="black")
    axis(3, at=regionBoundaries, labels=names(regionBoundaries), tick = FALSE)
    abline(h=0, col="gray")
    abline(v=regional25MbBoundaries, col="gray", lty="dotted", lw=0.75)    
    title(main=paste(dataDesc, sep=""), outer=TRUE)
    title(xlab="chrX (Mb)", outer=TRUE, line=1)
    
    # Legend
    par(mar=c(0,0,0,0))
    plot(seq(0,10), type="n", axes=FALSE, xlab="", ylab="")
    legend("top", 
           legend=c(
             paste("Del", sep=""),
             paste("InvDxz4", sep="")
           ),
           lty=c(1,1,1), col=c(peakCols[["Del"]], peakCols[["InvDxz4"]]), bty='n', cex=1)
    
    # add total counts table
    theTotCounts = unlist(totalCounts[[theDataType]])[grep(paste(theAllele, "$", sep=""), names(unlist(totalCounts[[theDataType]])) )]
    theTotCounts.scaled = unlist(totalCounts.scaled[[theDataType]])[grep(paste(theAllele, "$", sep=""), names(unlist(totalCounts.scaled[[theDataType]])) )]
    statsTable = cbind(theTotCounts, round(PeakRatioScalingFactors[[theDataType]], 2), round(theTotCounts.scaled))
    colnames(statsTable) = c("total", "Xa SF", "scaled")
    rownames(statsTable) = samples 
    
    par(mar=c(0,0,2,0))
    plot(seq(1,10), type="n", axes=FALSE, xlab="", ylab="")
    # https://stackoverflow.com/questions/18663159/conditional-coloring-of-cells-in-table
    # require(plotrix)
    addtable2plot("top", table=statsTable, bg="white", title="Total peaks before and after scaling",
                  bty="o", display.rownames=TRUE, display.colnames=TRUE,
                  hlines=TRUE, vlines=TRUE,
                  xjust=2, yjust=1, cex=1)
    
    dev.off()
  }
}

















#########################
# 180129
# Regional density plots

densityPlotDir = file.path(outputDir, "regionalDensityPlots")
if ( ! file.exists( densityPlotDir ) ) {
  if ( ! dir.create( densityPlotDir, showWarnings=FALSE, recursive=TRUE )  ) densityPlotDir = outputDir
}

# regionMidpoints = c(0, regionBoundaries, 180 * 1000000)[1:4]+ (c(0, regionBoundaries, 180 * 1000000)[2:5] - c(0, regionBoundaries, 180 * 1000000)[1:4])/2

for (theDataType in dataTypes) {
  theXlims = c(0, 180 * (1000000))
  theYlims = c(0, regionalMaxY[[theDataType]])
  
  for(theAllele in c(alleles, "Xi-assoc", "Xa-assoc")) {
    
    dataDesc = paste("regionalBinnedPeakCounts_scaled_", theDataType, "_", theAllele, sep="")
    pdf(file.path(densityPlotDir, paste(dataDesc,".pdf",sep="")), width=20, height=5)
    
    # Layout
    layout(matrix(c(1,1,1,1,2,3,4,5), nrow=4, ncol=2, byrow=FALSE), widths=c(1,0.2))
    # layout.show(4)
    par(mar=c(2,4,6,4), oma=c(3,0,2,0))
    
    # Plot WT counts
    # Plot an empty graph and legend to get the size of the legend
    # https://stackoverflow.com/questions/8929663/r-legend-placement-in-a-plot
    plot(c(0, regionBoundaries, 180 * 1000000), c(regionalBinnedPeakCounts.scaled[[theDataType]][["WT"]][[theAllele]], 
                                                  rev(regionalBinnedPeakCounts.scaled[[theDataType]][["WT"]][[theAllele]])[1]), 
         type="n", xaxt="n", yaxt="n",  xlab=NA, ylab=NA, ylim=theYlims, xlim=theXlims)
    # # Custom ylim. Add the height of legend to upper bound of the range
    # ymax = max(c(daYlim, ymax), na.rm=TRUE)
    # my.range <- c(0, ymax)
    # ##my.range[2] <- 1.04*(my.range[2]+my.legend.size$rect$h)
    
    # Draw the WT plot with custom ylim
    par(new=TRUE)
    plot(c(0, regionBoundaries, 180 * 1000000), c(regionalBinnedPeakCounts.scaled[[theDataType]][["WT"]][[theAllele]], 
                                                  rev(regionalBinnedPeakCounts.scaled[[theDataType]][["WT"]][[theAllele]])[1]), 
         type="s", lty="solid", lwd=2, col=peakCols[["WT"]], axes=F, main=NA, xlab=NA, ylab=NA, ylim=theYlims, xlim=theXlims)
    axis(side=2, line=0)
    mtext(side=2, line=2.5, text="Normalized counts")
    axis(1, at=c(0, regionBoundaries, 180 * 1000000), labels=c("0", "Firre", "Dxz4", "Xist", "180"))
    title(main=paste("Counts (", scalingText, ")", sep=""), line=3)
    
    # Plot other samples
    for (theSample in samples) {
      if (theSample != "WT") {
        lines(c(0, regionBoundaries, 180 * 1000000), c(regionalBinnedPeakCounts.scaled[[theDataType]][[theSample]][[theAllele]], 
                                                       rev(regionalBinnedPeakCounts.scaled[[theDataType]][[theSample]][[theAllele]])[1]),
              type="s", lty="solid", lwd=2, col=peakCols[[theSample]])
        # abline(h=expectedCounts[[theSample]], lty="dashed", col=peakCols[[theSample]])
      }
    }
    
    # LOI
    abline(v=regionBoundaries["Firre"], lty="dashed", lw=0.5, col="black")
    abline(v=regionBoundaries["Dxz4"], lty="dashed", lw=0.5, col="black")
    abline(v=regionBoundaries["Xist"], lty="dashed", lw=0.5, col="black")
    axis(3, at=regionBoundaries, labels=names(regionBoundaries), tick = FALSE)
    
    title(main=paste(dataDesc, sep=""), outer=TRUE)
    title(xlab="chrX (Mb)", outer=TRUE, line=1)
    
    # Legend
    par(mar=c(0,0,0,0))
    plot(seq(0,10), type="n", axes=FALSE, xlab="", ylab="")
    legend("top", 
           legend=c(
             paste("WT", sep=""),
             paste("Del", sep=""),
             paste("InvDxz4", sep="")
           ),
           lty=c(1,1,1), col=c(peakCols[["WT"]], peakCols[["Del"]], peakCols[["InvDxz4"]]), bty='n', cex=1)
    
    # add total counts table
    theTotCounts = unlist(totalCounts[[theDataType]])[grep(paste(theAllele, "$", sep=""), names(unlist(totalCounts[[theDataType]])) )]
    theTotCounts.scaled = unlist(totalCounts.scaled[[theDataType]])[grep(paste(theAllele, "$", sep=""), names(unlist(totalCounts.scaled[[theDataType]])) )]
    statsTable = cbind(theTotCounts, round(PeakRatioScalingFactors[[theDataType]], 2), round(theTotCounts.scaled))
    colnames(statsTable) = c("total", "Xa SF", "scaled")
    rownames(statsTable) = samples 
    
    par(mar=c(0,0,2,0))
    plot(seq(1,10), type="n", axes=FALSE, xlab="", ylab="")
    # https://stackoverflow.com/questions/18663159/conditional-coloring-of-cells-in-table
    # require(plotrix)
    addtable2plot("top", table=statsTable, bg="white", title="Total peaks before and after scaling",
                  bty="o", display.rownames=TRUE, display.colnames=TRUE,
                  hlines=TRUE, vlines=TRUE,
                  xjust=2, yjust=1, cex=1)
    
    # add regional counts table
    # theRegCounts = unlist(regionalBinnedPeakCounts[[theDataType]])[grep(paste(theAllele, "$", sep=""), names(unlist(regionalBinnedPeakCounts[[theDataType]])) )]
    statsTable = c()
    for (theSample in samples) {
      statsTable = rbind(statsTable,
                         regionalBinnedPeakCounts.scaled[[theDataType]][[theSample]][[theAllele]])
    }
    statsTable = round(statsTable)
    colnames(statsTable) = paste( c("start", names(regionBoundaries)), "-", c(names(regionBoundaries), "stop"), sep="")
    rownames(statsTable) = samples 
    
    par(mar=c(0,0,2,0))
    plot(seq(1,10), type="n", axes=FALSE, xlab="", ylab="")
    # https://stackoverflow.com/questions/18663159/conditional-coloring-of-cells-in-table
    # require(plotrix)
    addtable2plot("top", table=statsTable, bg="white", title="Regional peaks (scaled)",
                  bty="o", display.rownames=TRUE, display.colnames=TRUE,
                  hlines=TRUE, vlines=TRUE,
                  xjust=2, yjust=1, cex=1)
    
    dev.off()
  }
}




#########################
# 180129
# Relative regional density plots

densityPlotDir = file.path(outputDir, "log2RatioRegionalDensityPlots")
if ( ! file.exists( densityPlotDir ) ) {
  if ( ! dir.create( densityPlotDir, showWarnings=FALSE, recursive=TRUE )  ) densityPlotDir = outputDir
}

# regionMidpoints = c(0, regionBoundaries, 180 * 1000000)[1:4]+ (c(0, regionBoundaries, 180 * 1000000)[2:5] - c(0, regionBoundaries, 180 * 1000000)[1:4])/2

for (theDataType in dataTypes) {
  theXlims = c(0, 180 * (1000000))
  theYlims = c(-2, +6)
  
  for(theAllele in c(alleles, "Xi-assoc", "Xa-assoc")) {
    
    dataDesc = paste("log2RatioRegionalBinnedPeakCounts_scaled_", theDataType, "_", theAllele, sep="")
    pdf(file.path(densityPlotDir, paste(dataDesc,".pdf",sep="")), width=20, height=5)
    
    # Layout
    layout(matrix(c(1,1,1,1,2,3,4,5), nrow=4, ncol=2, byrow=FALSE), widths=c(1,0.2))
    # layout.show(4)
    par(mar=c(2,4,6,4), oma=c(3,0,2,0))
    
    # Plot WT counts
    # Plot an empty graph and legend to get the size of the legend
    # https://stackoverflow.com/questions/8929663/r-legend-placement-in-a-plot
    plot(c(0, regionBoundaries, 180 * 1000000), rep(0, length(regionBoundaries)+2), 
         type="n", xaxt="n", yaxt="n",  xlab=NA, ylab=NA, ylim=theYlims, xlim=theXlims)
    # # Custom ylim. Add the height of legend to upper bound of the range
    # ymax = max(c(daYlim, ymax), na.rm=TRUE)
    # my.range <- c(0, ymax)
    # ##my.range[2] <- 1.04*(my.range[2]+my.legend.size$rect$h)
    
    # Draw the WT plot with custom ylim
    par(new=TRUE)
    plot(c(0, regionBoundaries, 180 * 1000000), 
         log2(c(regionalBinnedPeakCounts.scaled[[theDataType]][["Del"]][[theAllele]], 
                rev(regionalBinnedPeakCounts.scaled[[theDataType]][["Del"]][[theAllele]])[1]) / 
                c(regionalBinnedPeakCounts.scaled[[theDataType]][["WT"]][[theAllele]], 
                rev(regionalBinnedPeakCounts.scaled[[theDataType]][["WT"]][[theAllele]])[1])), 
         type="s", lty="solid", lwd=2, col=peakCols[["Del"]], axes=F, main=NA, xlab=NA, ylab=NA, ylim=theYlims, xlim=theXlims)
    axis(side=2, line=0)
    mtext(side=2, line=2.5, text="Log ratio")
    axis(1, at=c(0, regionBoundaries, 180 * 1000000), labels=c("0", "Firre", "Dxz4", "Xist", "180"))
    title(main=paste("Log2Ratio vs. WT (", scalingText, ")", sep=""), line=3)
    
    # Plot other samples
    lines(c(0, regionBoundaries, 180 * 1000000), 
          log2(c(regionalBinnedPeakCounts.scaled[[theDataType]][["InvDxz4"]][[theAllele]], 
                 rev(regionalBinnedPeakCounts.scaled[[theDataType]][["InvDxz4"]][[theAllele]])[1]) / 
                 c(regionalBinnedPeakCounts.scaled[[theDataType]][["WT"]][[theAllele]], 
                   rev(regionalBinnedPeakCounts.scaled[[theDataType]][["WT"]][[theAllele]])[1])), 
          type="s", lty="solid", lwd=2, col=peakCols[["InvDxz4"]])

    # LOI
    abline(v=regionBoundaries["Firre"], lty="dashed", lw=0.5, col="black")
    abline(v=regionBoundaries["Dxz4"], lty="dashed", lw=0.5, col="black")
    abline(v=regionBoundaries["Xist"], lty="dashed", lw=0.5, col="black")
    axis(3, at=regionBoundaries, labels=names(regionBoundaries), tick = FALSE)
    abline(h=0, col="gray")
    
    title(main=paste(dataDesc, sep=""), outer=TRUE)
    title(xlab="chrX (Mb)", outer=TRUE, line=1)
    
    # Legend
    par(mar=c(0,0,0,0))
    plot(seq(0,10), type="n", axes=FALSE, xlab="", ylab="")
    legend("top", 
           legend=c(
             paste("Del", sep=""),
             paste("InvDxz4", sep="")
           ),
           lty=c(1,1,1), col=c(peakCols[["Del"]], peakCols[["InvDxz4"]]), bty='n', cex=1)
    
    # add total counts table
    theTotCounts = unlist(totalCounts[[theDataType]])[grep(paste(theAllele, "$", sep=""), names(unlist(totalCounts[[theDataType]])) )]
    theTotCounts.scaled = unlist(totalCounts.scaled[[theDataType]])[grep(paste(theAllele, "$", sep=""), names(unlist(totalCounts.scaled[[theDataType]])) )]
    statsTable = cbind(theTotCounts, round(PeakRatioScalingFactors[[theDataType]], 2), round(theTotCounts.scaled))
    colnames(statsTable) = c("total", "Xa SF", "scaled")
    rownames(statsTable) = samples 
    
    par(mar=c(0,0,2,0))
    plot(seq(1,10), type="n", axes=FALSE, xlab="", ylab="")
    # https://stackoverflow.com/questions/18663159/conditional-coloring-of-cells-in-table
    # require(plotrix)
    addtable2plot("top", table=statsTable, bg="white", title="Total peaks before and after scaling",
                  bty="o", display.rownames=TRUE, display.colnames=TRUE,
                  hlines=TRUE, vlines=TRUE,
                  xjust=2, yjust=1, cex=1)
    
    # add regional counts table
    # theRegCounts = unlist(regionalBinnedPeakCounts[[theDataType]])[grep(paste(theAllele, "$", sep=""), names(unlist(regionalBinnedPeakCounts[[theDataType]])) )]
    statsTable = c()
    for (theSample in samples) {
      statsTable = rbind(statsTable,
                         regionalBinnedPeakCounts.scaled[[theDataType]][[theSample]][[theAllele]])
    }
    statsTable = round(statsTable)
    colnames(statsTable) = paste( c("start", names(regionBoundaries)), "-", c(names(regionBoundaries), "stop"), sep="")
    rownames(statsTable) = samples 
    
    par(mar=c(0,0,2,0))
    plot(seq(1,10), type="n", axes=FALSE, xlab="", ylab="")
    # https://stackoverflow.com/questions/18663159/conditional-coloring-of-cells-in-table
    # require(plotrix)
    addtable2plot("top", table=statsTable, bg="white", title="Regional peaks (scaled)",
                  bty="o", display.rownames=TRUE, display.colnames=TRUE,
                  hlines=TRUE, vlines=TRUE,
                  xjust=2, yjust=1, cex=1)

    # add relative regional log2 table
    # theRegCounts = unlist(regionalBinnedPeakCounts[[theDataType]])[grep(paste(theAllele, "$", sep=""), names(unlist(regionalBinnedPeakCounts[[theDataType]])) )]
    statsTable = rbind(
      log2(regionalBinnedPeakCounts.scaled[[theDataType]][["Del"]][[theAllele]] / 
             regionalBinnedPeakCounts.scaled[[theDataType]][["WT"]][[theAllele]]),
      log2(regionalBinnedPeakCounts.scaled[[theDataType]][["Del"]][[theAllele]] / 
             regionalBinnedPeakCounts.scaled[[theDataType]][["InvDxz4"]][[theAllele]])    )
    
    statsTable = round(statsTable,1)
    colnames(statsTable) = paste( c("start", names(regionBoundaries)), "-", c(names(regionBoundaries), "stop"), sep="")
    rownames(statsTable) = samples[-1]
    
    par(mar=c(0,0,2,0))
    plot(seq(1,10), type="n", axes=FALSE, xlab="", ylab="")
    # https://stackoverflow.com/questions/18663159/conditional-coloring-of-cells-in-table
    # require(plotrix)
    addtable2plot("top", table=statsTable, bg="white", title="Regional log2 ratio vs. WT (scaled)",
                  bty="o", display.rownames=TRUE, display.colnames=TRUE,
                  hlines=TRUE, vlines=TRUE,
                  xjust=2, yjust=1, cex=1)
    
        
    dev.off()
  }
}







#########################
# Plots density plots

densityPlotDir = file.path(outputDir, "densityPlots")
if ( ! file.exists( densityPlotDir ) ) {
  if ( ! dir.create( densityPlotDir, showWarnings=FALSE, recursive=TRUE )  ) densityPlotDir = outputDir
}

for (theDataType in dataTypes) {
  theXlims = c(0, 180 * (1000000/binSize))
  theYlims = c(0, maxY[[theDataType]])
  
  for(theAllele in c(alleles, "Xi-assoc", "Xa-assoc")) {

    dataDesc = paste("binnedPeakCounts_", format(binSize, scientific=FALSE), "_scaled_", theDataType, "_", theAllele, sep="")
    pdf(file.path(densityPlotDir, paste(dataDesc,".pdf",sep="")), width=20, height=5)

    # Layout
    layout(matrix(c(1,1,1,1,2,3,4,5), nrow=4, ncol=2, byrow=FALSE), widths=c(1,0.1))
    # layout.show(4)
    par(mar=c(2,4,6,4), oma=c(3,0,2,0))
    
    # Plot WT counts
    # Plot an empty graph and legend to get the size of the legend
    # https://stackoverflow.com/questions/8929663/r-legend-placement-in-a-plot
    plot(binnedPeakCounts.scaled[[theDataType]][["WT"]][[theAllele]], type="n", xaxt="n", yaxt="n",  xlab=NA, ylab=NA, ylim=theYlims, xlim=theXlims)
    # # Custom ylim. Add the height of legend to upper bound of the range
    # ymax = max(c(daYlim, ymax), na.rm=TRUE)
    # my.range <- c(0, ymax)
    # ##my.range[2] <- 1.04*(my.range[2]+my.legend.size$rect$h)
    
    # Draw the WT plot with custom ylim
    par(new=TRUE)
    plot(binnedPeakCounts.scaled[[theDataType]][["WT"]][[theAllele]], type="l", lty="solid", lwd=2, col=peakCols[["WT"]], axes=F, main=NA, xlab=NA, ylab=NA, ylim=theYlims, xlim=theXlims)
    axis(side=2, line=0)
    mtext(side=2, line=2.5, text="Normalized counts")
    axis(1, at=seq(0, 180 * (1000000/binSize), 20 * (1000000/binSize)) + 1, labels=seq(0, 180, 20))
    title(main=paste("Counts (", scalingText, ")", sep=""), line=3)
    
    # Plot other samples
    for (theSample in samples) {
      if (theSample != "WT") {
        lines(binnedPeakCounts.scaled[[theDataType]][[theSample]][[theAllele]], type="l", lty="solid", lwd=2, col=peakCols[[theSample]])
        # abline(h=expectedCounts[[theSample]], lty="dashed", col=peakCols[[theSample]])
      }
    }
    
    # LOI
    abline(v=firrebin, lty="dashed", lw=0.5, col="black")
    abline(v=dxz4bin, lty="dashed", lw=0.5, col="black")
    abline(v=xistbin, lty="dashed", lw=0.5, col="black")
    axis(3, at=c(firrebin, dxz4bin, xistbin), labels=c("Firre", "Dxz4", "Xist"), tick = FALSE)
    
    title(main=paste(dataDesc, " (", format(binSize, scientific=FALSE), "bp bins)", sep=""), outer=TRUE)
    title(xlab="chrX (Mb)", outer=TRUE, line=1)
    
    # Legend
    par(mar=c(0,0,0,0))
    plot(seq(0,10), type="n", axes=FALSE, xlab="", ylab="")
    legend("top", 
    legend=c(
    paste("WT", sep=""),
    paste("Del", sep=""),
    paste("InvDxz4", sep="")
    ),
    lty=c(1,1,1), col=c(peakCols[["WT"]], peakCols[["Del"]], peakCols[["InvDxz4"]]), bty='n', cex=1)
    
    # add counts table
    theTotCounts = unlist(totalCounts[[theDataType]])[grep(paste(theAllele, "$", sep=""), names(unlist(totalCounts[[theDataType]])) )]
    theTotCounts.scaled = unlist(totalCounts.scaled[[theDataType]])[grep(paste(theAllele, "$", sep=""), names(unlist(totalCounts.scaled[[theDataType]])) )]
    statsTable = cbind(theTotCounts, round(PeakRatioScalingFactors[[theDataType]], 2), round(theTotCounts.scaled))
    colnames(statsTable) = c("total", "Xa SF", "scaled")
    rownames(statsTable) = samples 
    
    par(mar=c(0,0,2,0))
    plot(seq(1,10), type="n", axes=FALSE, xlab="", ylab="")
    # https://stackoverflow.com/questions/18663159/conditional-coloring-of-cells-in-table
    # require(plotrix)
    addtable2plot("top", table=statsTable, bg="white", title="Total peaks before and after scaling",
    bty="o", display.rownames=TRUE, display.colnames=TRUE,
    hlines=TRUE, vlines=TRUE,
    xjust=2, yjust=1, cex=1)
    
    # add correlation table
    Pcorr = c()
    Scorr = c()
    for (theSampleA in samples) {
    PcorrVect = c()
    ScorrVect = c()
    for (theSampleB in samples) {
    PcorrVect = c(PcorrVect, round(cor(binnedPeakCounts.scaled[[theDataType]][[theSampleA]][[theAllele]], binnedPeakCounts.scaled[[theDataType]][[theSampleB]][[theAllele]], method="pearson"),2))
    ScorrVect = c(ScorrVect, round(cor(binnedPeakCounts.scaled[[theDataType]][[theSampleA]][[theAllele]], binnedPeakCounts.scaled[[theDataType]][[theSampleB]][[theAllele]], method="spearman"),2))
    }
    Pcorr =rbind(Pcorr, PcorrVect)
    Scorr =rbind(Scorr, ScorrVect)
    }
    rownames(Pcorr) = colnames(Pcorr) = rownames(Scorr) = colnames(Scorr) = samples
    
    par(mar=c(0,0,2,0))
    plot(seq(1,10), type="n", axes=FALSE, xlab="", ylab="")
    # https://stackoverflow.com/questions/18663159/conditional-coloring-of-cells-in-table
    # require(plotrix)
    addtable2plot("top", table=Pcorr, bg="white", title="Pearson correlation matrix",
    bty="o", display.rownames=TRUE, display.colnames=TRUE,
    hlines=TRUE, vlines=TRUE,
    xjust=2, yjust=1, cex=1)
    
    par(mar=c(0,0,2,0))
    plot(seq(1,10), type="n", axes=FALSE, xlab="", ylab="")
    # https://stackoverflow.com/questions/18663159/conditional-coloring-of-cells-in-table
    # require(plotrix)
    addtable2plot("top", table=Scorr, bg="white", title="Spearman correlation matrix",
    bty="o", display.rownames=TRUE, display.colnames=TRUE,
    hlines=TRUE, vlines=TRUE,
    xjust=2, yjust=1, cex=1)
    
    dev.off()
  }
}



#########################
# 180129
# Relative density plots

densityPlotDir = file.path(outputDir, "log2RatioDensityPlots")
if ( ! file.exists( densityPlotDir ) ) {
  if ( ! dir.create( densityPlotDir, showWarnings=FALSE, recursive=TRUE )  ) densityPlotDir = outputDir
}

for (theDataType in dataTypes) {
  theXlims = c(0, 180 * (1000000/binSize))
  theYlims = c(-2, +6)
  
  for(theAllele in c(alleles, "Xi-assoc", "Xa-assoc")) {
    
    dataDesc = paste("binnedPeakCounts_", format(binSize, scientific=FALSE), "_scaled_", theDataType, "_", theAllele, sep="")
    pdf(file.path(densityPlotDir, paste(dataDesc,".pdf",sep="")), width=20, height=5)
    
    # Layout
    layout(matrix(c(1,1,1,1,2,3,4,5), nrow=4, ncol=2, byrow=FALSE), widths=c(1,0.1))
    # layout.show(4)
    par(mar=c(2,4,6,4), oma=c(3,0,2,0))
    
    # Plot WT counts
    # Plot an empty graph and legend to get the size of the legend
    # https://stackoverflow.com/questions/8929663/r-legend-placement-in-a-plot
    plot(rep(0, length(binnedPeakCounts.scaled[[theDataType]][["WT"]][[theAllele]])), type="n", xaxt="n", yaxt="n",  xlab=NA, ylab=NA, ylim=theYlims, xlim=theXlims)
    # # Custom ylim. Add the height of legend to upper bound of the range
    # ymax = max(c(daYlim, ymax), na.rm=TRUE)
    # my.range <- c(0, ymax)
    # ##my.range[2] <- 1.04*(my.range[2]+my.legend.size$rect$h)
    
    # Draw the WT plot with custom ylim
    par(new=TRUE)
    plot(log2((binnedPeakCounts.scaled[[theDataType]][["Del"]][[theAllele]] + 1) / 
                (binnedPeakCounts.scaled[[theDataType]][["WT"]][[theAllele]]) + 1), 
         type="h", lty="solid", lwd=2, col=peakCols[["Del"]], 
         axes=F, main=NA, xlab=NA, ylab=NA, ylim=theYlims, xlim=theXlims)
    axis(side=2, line=0)
    mtext(side=2, line=2.5, text="Normalized counts")
    axis(1, at=seq(0, 180 * (1000000/binSize), 20 * (1000000/binSize)) + 1, labels=seq(0, 180, 20))
    title(main=paste("Log2Ratio vs. WT (", scalingText, ")", sep=""), line=3)
    
    # Plot other samples
    lines(log2(binnedPeakCounts.scaled[[theDataType]][["InvDxz4"]][[theAllele]] / 
                 binnedPeakCounts.scaled[[theDataType]][["WT"]][[theAllele]]), 
          type="h", lty="solid", lwd=2, col=peakCols[["InvDxz4"]])

    # LOI
    abline(v=firrebin, lty="dashed", lw=0.5, col="black")
    abline(v=dxz4bin, lty="dashed", lw=0.5, col="black")
    abline(v=xistbin, lty="dashed", lw=0.5, col="black")
    axis(3, at=c(firrebin, dxz4bin, xistbin), labels=c("Firre", "Dxz4", "Xist"), tick = FALSE)
    abline(h=0, col="gray")
    
    title(main=paste(dataDesc, " (", format(binSize, scientific=FALSE), "bp bins)", sep=""), outer=TRUE)
    title(xlab="chrX (Mb)", outer=TRUE, line=1)
    
    # Legend
    par(mar=c(0,0,0,0))
    plot(seq(0,10), type="n", axes=FALSE, xlab="", ylab="")
    legend("top", 
           legend=c(
             paste("WT", sep=""),
             paste("Del", sep=""),
             paste("InvDxz4", sep="")
           ),
           lty=c(1,1,1), col=c(peakCols[["WT"]], peakCols[["Del"]], peakCols[["InvDxz4"]]), bty='n', cex=1)
    
    # add counts table
    theTotCounts = unlist(totalCounts[[theDataType]])[grep(paste(theAllele, "$", sep=""), names(unlist(totalCounts[[theDataType]])) )]
    theTotCounts.scaled = unlist(totalCounts.scaled[[theDataType]])[grep(paste(theAllele, "$", sep=""), names(unlist(totalCounts.scaled[[theDataType]])) )]
    statsTable = cbind(theTotCounts, round(PeakRatioScalingFactors[[theDataType]], 2), round(theTotCounts.scaled))
    colnames(statsTable) = c("total", "Xa SF", "scaled")
    rownames(statsTable) = samples 
    
    par(mar=c(0,0,2,0))
    plot(seq(1,10), type="n", axes=FALSE, xlab="", ylab="")
    # https://stackoverflow.com/questions/18663159/conditional-coloring-of-cells-in-table
    # require(plotrix)
    addtable2plot("top", table=statsTable, bg="white", title="Total peaks before and after scaling",
                  bty="o", display.rownames=TRUE, display.colnames=TRUE,
                  hlines=TRUE, vlines=TRUE,
                  xjust=2, yjust=1, cex=1)
    
    # add correlation table
    Pcorr = c()
    Scorr = c()
    for (theSampleA in samples) {
      PcorrVect = c()
      ScorrVect = c()
      for (theSampleB in samples) {
        PcorrVect = c(PcorrVect, round(cor(binnedPeakCounts.scaled[[theDataType]][[theSampleA]][[theAllele]], binnedPeakCounts.scaled[[theDataType]][[theSampleB]][[theAllele]], method="pearson"),2))
        ScorrVect = c(ScorrVect, round(cor(binnedPeakCounts.scaled[[theDataType]][[theSampleA]][[theAllele]], binnedPeakCounts.scaled[[theDataType]][[theSampleB]][[theAllele]], method="spearman"),2))
      }
      Pcorr =rbind(Pcorr, PcorrVect)
      Scorr =rbind(Scorr, ScorrVect)
    }
    rownames(Pcorr) = colnames(Pcorr) = rownames(Scorr) = colnames(Scorr) = samples
    
    par(mar=c(0,0,2,0))
    plot(seq(1,10), type="n", axes=FALSE, xlab="", ylab="")
    # https://stackoverflow.com/questions/18663159/conditional-coloring-of-cells-in-table
    # require(plotrix)
    addtable2plot("top", table=Pcorr, bg="white", title="Pearson correlation matrix",
                  bty="o", display.rownames=TRUE, display.colnames=TRUE,
                  hlines=TRUE, vlines=TRUE,
                  xjust=2, yjust=1, cex=1)
    
    par(mar=c(0,0,2,0))
    plot(seq(1,10), type="n", axes=FALSE, xlab="", ylab="")
    # https://stackoverflow.com/questions/18663159/conditional-coloring-of-cells-in-table
    # require(plotrix)
    addtable2plot("top", table=Scorr, bg="white", title="Spearman correlation matrix",
                  bty="o", display.rownames=TRUE, display.colnames=TRUE,
                  hlines=TRUE, vlines=TRUE,
                  xjust=2, yjust=1, cex=1)

    dev.off()
  }
}
      



#########################
# Plots special density plots

densityPlotDir = file.path(outputDir, "specialDensityPlots")
if ( ! file.exists( densityPlotDir ) ) {
  if ( ! dir.create( densityPlotDir, showWarnings=FALSE, recursive=TRUE )  ) densityPlotDir = outputDir
}

theAllele = "Xi-assoc"
dataDesc = paste("binnedPeakCounts_", format(binSize, scientific=FALSE), "_scaled_", theAllele, sep="")
pdf(file.path(densityPlotDir, paste(dataDesc,".pdf",sep="")), width=20, height=(2+1+1+2)*5)

par(mfrow=c(2+1+1+2, 1), mar=c(2,4,6,4), oma=c(3,0,2,0))


###
# ATAC and CTCF density
for (theDataType in rev(dataTypes)) {
  #theXlims = c(0, 180 * (1000000/binSize))
  theXlims = c(1, 343)
  #theYlims = c(0, maxY[[theDataType]])
  theYlims = c(0, max(binnedPeakCounts.scaled[[theDataType]][["Del"]][[theAllele]], na.rm=TRUE))

  plot(binnedPeakCounts.scaled[[theDataType]][["WT"]][[theAllele]], type="l", lty="solid", lwd=2, col=peakCols[["WT"]], axes=F, main=NA, xlab=NA, ylab=NA, ylim=theYlims, xlim=theXlims)
  axis(side=2, line=0)
  mtext(side=2, line=2.5, text=paste(theDataType, " normalized counts", sep=""))
  axis(1, at=seq(0, 180 * (1000000/binSize), 20 * (1000000/binSize)) + 1, labels=seq(0, 180, 20))
  title(main=paste(theDataType, " binned counts (scaled)", sep=""), line=3)
  
  # Plot other samples
  for (theSample in samples) {
    if (theSample != "WT") {
      lines(binnedPeakCounts.scaled[[theDataType]][[theSample]][[theAllele]], type="l", lty="solid", lwd=2, col=peakCols[[theSample]])
      # abline(h=expectedCounts[[theSample]], lty="dashed", col=peakCols[[theSample]])
    }
  }
  
  # LOI
  abline(v=firrebin, lty="dashed", lw=0.5, col="black")
  abline(v=dxz4bin, lty="dashed", lw=0.5, col="black")
  abline(v=xistbin, lty="dashed", lw=0.5, col="black")
  axis(3, at=c(firrebin, dxz4bin, xistbin), labels=c("Firre", "Dxz4", "Xist"), tick = FALSE)
}

# Correlate ATAC and CTCF scores
PcorrVect = c()
ScorrVect = c()
for (theSample in samples) {
    ATACvalues = binnedPeakCounts.scaled[["ATAC"]][[theSample]][[theAllele]][1:343]
    CTCFvalues = binnedPeakCounts.scaled[["CTCF"]][[theSample]][[theAllele]][1:343]
    PcorrVect = c(PcorrVect, round(cor(ATACvalues, CTCFvalues, method="pearson", use="pairwise.complete"),2))
    ScorrVect = c(ScorrVect, round(cor(ATACvalues, CTCFvalues, method="spearman", use="pairwise.complete"),2))
}
names(PcorrVect) = names(ScorrVect) = samples
print(paste("ATAC vs CTCF", sep=""))
print(ScorrVect)

# https://stackoverflow.com/questions/18663159/conditional-coloring-of-cells-in-table
# require(plotrix)
addtable2plot("topright", table=t(as.matrix(ScorrVect)), bg="white", title=paste("ATAC vs CTCF", sep=""),
              bty="o", display.rownames=FALSE, display.colnames=TRUE,
              hlines=TRUE, vlines=TRUE,
              xjust=2, yjust=1, cex=1.5)


###
# Gene density
#theXlims = c(0, 180 * (1000000/binSize))
theXlims = c(1, 343)

plot(binnedGeneCounts, type="l", lty="solid", lwd=2, col="green", axes=F, main=NA, xlab=NA, ylab=NA, xlim=theXlims)
axis(side=2, line=0)
mtext(side=2, line=2.5, text=paste("binnned gene counts", sep=""))
axis(1, at=seq(0, 180 * (1000000/binSize), 20 * (1000000/binSize)) + 1, labels=seq(0, 180, 20))
title(main=paste("Binned gene counts (", format(binSize, scientific=FALSE), "bp)", sep=""), line=3)

# LOI
abline(v=firrebin, lty="dashed", lw=0.5, col="black")
abline(v=dxz4bin, lty="dashed", lw=0.5, col="black")
abline(v=xistbin, lty="dashed", lw=0.5, col="black")
axis(3, at=c(firrebin, dxz4bin, xistbin), labels=c("Firre", "Dxz4", "Xist"), tick = FALSE)

# Correlate Gene density with features
Pcorr = c()
Scorr = c()
for (theDataType in rev(dataTypes)) {
  PcorrVect = c()
  ScorrVect = c()
  for (theSample in samples) {
    featureValues = binnedPeakCounts.scaled[[theDataType]][["WT"]][[theAllele]][1:343]
    PcorrVect = c(PcorrVect, round(cor(featureValues, binnedGeneCounts[1:343], method="pearson", use="pairwise.complete"),2))
    ScorrVect = c(ScorrVect, round(cor(featureValues, binnedGeneCounts[1:343], method="spearman", use="pairwise.complete"),2))
  }
  Pcorr =rbind(Pcorr, PcorrVect)
  Scorr =rbind(Scorr, ScorrVect)
}
rownames(Pcorr) = rownames(Scorr) = rev(dataTypes)
colnames(Pcorr) =  colnames(Scorr) = samples
print(paste("Gene density", sep=""))
print(Scorr)

# https://stackoverflow.com/questions/18663159/conditional-coloring-of-cells-in-table
# require(plotrix)
addtable2plot("topright", table=Scorr, bg="white", title=paste("Gene density", sep=""),
              bty="o", display.rownames=TRUE, display.colnames=TRUE,
              hlines=TRUE, vlines=TRUE,
              xjust=2, yjust=1, cex=1.5)


###
# Now plot Xa PC1 scores
for(PCnum in paste("PC", 1, sep="")) {
  theXlims = c(1, 343)
  theYlims = PCylims[[PCnum]][["Xa"]]
  
  plot(PCscoreData[["WT"]][["Xa"]][,PCnum], type="l", lty="solid", lwd=2, col=peakCols[["WT"]], axes=F, main=NA, xlab=NA, ylab=NA, ylim=theYlims, xlim=theXlims)
  axis(side=2, line=0)
  mtext(side=2, line=2.5, text=paste("Xa ", PCnum, " scores", sep=""))
  axis(1, at=seq(0, 180 * (1000000/binSize), 20 * (1000000/binSize)) + 1, labels=seq(0, 180, 20))
  title(main=paste("Xa ", PCnum, sep=""), line=3)
  
  # Plot other samples
  for (theSample in samples) {
    if (theSample != "WT") {
      lines(PCscoreData[[theSample]][["Xa"]][,PCnum], type="l", lty="solid", lwd=2, col=peakCols[[theSample]])
      # abline(h=expectedCounts[[theSample]], lty="dashed", col=peakCols[[theSample]])
    }
  }
  
  # LOI
  abline(v=firrebin, lty="dashed", lw=0.5, col="black")
  abline(v=dxz4bin, lty="dashed", lw=0.5, col="black")
  abline(v=xistbin, lty="dashed", lw=0.5, col="black")
  axis(3, at=c(firrebin, dxz4bin, xistbin), labels=c("Firre", "Dxz4", "Xist"), tick = FALSE)

  # Correlate PC scores with features
  Pcorr = c()
  Scorr = c()
  for (theDataType in rev(dataTypes)) {
    PcorrVect = c()
    ScorrVect = c()
    for (theSample in samples) {
      featureValues = binnedPeakCounts.scaled[[theDataType]][["WT"]][[theAllele]][1:343]
      PCvalues = PCscoreData[[theSample]][["Xa"]][,PCnum]
      PcorrVect = c(PcorrVect, round(cor(featureValues, PCvalues, method="pearson", use="pairwise.complete"),2))
      ScorrVect = c(ScorrVect, round(cor(featureValues, PCvalues, method="spearman", use="pairwise.complete"),2))
    }
    Pcorr =rbind(Pcorr, PcorrVect)
    Scorr =rbind(Scorr, ScorrVect)
  }
  rownames(Pcorr) = rownames(Scorr) = rev(dataTypes)
  colnames(Pcorr) =  colnames(Scorr) = samples
  print(paste("Xa ", PCnum, sep=""))
  print(Scorr)
  
  # https://stackoverflow.com/questions/18663159/conditional-coloring-of-cells-in-table
  # require(plotrix)
  addtable2plot("topright", table=Scorr, bg="white", title=paste("Xa ", PCnum, " Spearman", sep=""),
                bty="o", display.rownames=TRUE, display.colnames=TRUE,
                hlines=TRUE, vlines=TRUE,
                xjust=2, yjust=1, cex=1.5)
  
}  

###
# Now plot Xi PC2/3 scores
for(PCnum in paste("PC", 2:3, sep="")) {
  theXlims = c(1, 343)
  theYlims = PCylims[[PCnum]][["Xi"]]
  
  plot(PCscoreData[["WT"]][["Xi"]][,PCnum], type="l", lty="solid", lwd=2, col=peakCols[["WT"]], axes=F, main=NA, xlab=NA, ylab=NA, ylim=theYlims, xlim=theXlims)
  axis(side=2, line=0)
  mtext(side=2, line=2.5, text=paste("Xi ", PCnum, " scores", sep=""))
  axis(1, at=seq(0, 180 * (1000000/binSize), 20 * (1000000/binSize)) + 1, labels=seq(0, 180, 20))
  title(main=paste("Xi ", PCnum, sep=""), line=3)
  
  # Plot other samples
  for (theSample in samples) {
    if (theSample != "WT") {
      lines(PCscoreData[[theSample]][["Xi"]][,PCnum], type="l", lty="solid", lwd=2, col=peakCols[[theSample]])
      # abline(h=expectedCounts[[theSample]], lty="dashed", col=peakCols[[theSample]])
    }
  }
  
  # LOI
  abline(v=firrebin, lty="dashed", lw=0.5, col="black")
  abline(v=dxz4bin, lty="dashed", lw=0.5, col="black")
  abline(v=xistbin, lty="dashed", lw=0.5, col="black")
  axis(3, at=c(firrebin, dxz4bin, xistbin), labels=c("Firre", "Dxz4", "Xist"), tick = FALSE)
  
  # Correlate PC scores with features
  Pcorr = c()
  Scorr = c()
  for (theDataType in rev(dataTypes)) {
    PcorrVect = c()
    ScorrVect = c()
    for (theSample in samples) {
      featureValues = binnedPeakCounts.scaled[[theDataType]][["WT"]][[theAllele]][1:343]
      PCvalues = PCscoreData[[theSample]][["Xi"]][,PCnum]
      PcorrVect = c(PcorrVect, round(cor(featureValues, PCvalues, method="pearson", use="pairwise.complete"),2))
      ScorrVect = c(ScorrVect, round(cor(featureValues, PCvalues, method="spearman", use="pairwise.complete"),2))
    }
    Pcorr =rbind(Pcorr, PcorrVect)
    Scorr =rbind(Scorr, ScorrVect)
  }
  rownames(Pcorr) = rownames(Scorr) = rev(dataTypes)
  colnames(Pcorr) =  colnames(Scorr) = samples
  print(paste("Xi ", PCnum, sep=""))
  print(Scorr)
  
  # https://stackoverflow.com/questions/18663159/conditional-coloring-of-cells-in-table
  # require(plotrix)
  addtable2plot("topright", table=Scorr, bg="white", title=paste("Xi ", PCnum, " Spearman", sep=""),
                bty="o", display.rownames=TRUE, display.colnames=TRUE,
                hlines=TRUE, vlines=TRUE,
                xjust=2, yjust=1, cex=1.5)
}  

title(main=paste(dataDesc, " (", format(binSize, scientific=FALSE), "bp bins)", sep=""), outer=TRUE)
title(xlab="chrX (Mb)", outer=TRUE, line=1)

dev.off()










###########################################################
# Gene expression vs ATAC vs CTCF peak correlation

#########################
# Counts around genes
# TSS-5kb to TTS

# Get gene regions around genes
for ( upstreamTSSbp in c(5000, 10000, 50000, 100000) ) {
  # upstreamTSSbp=10000
  chrXupTSS = expressionTables.chrX[["counts"]]$start - upstreamTSSbp
  chrXupTSS[expressionTables.chrX[["counts"]]$strand == "-"] = expressionTables.chrX[["counts"]]$stop[expressionTables.chrX[["counts"]]$strand == "-"] + upstreamTSSbp
  chrXTTS = expressionTables.chrX[["counts"]]$stop
  chrXTTS[expressionTables.chrX[["counts"]]$strand == "-"] = expressionTables.chrX[["counts"]]$start[expressionTables.chrX[["counts"]]$strand == "-"]
  chrXgeneRanges = cbind(chrXupTSS, chrXTTS)
  chrXgeneRanges[expressionTables.chrX[["counts"]]$strand == "-"] = cbind(chrXTTS, chrXupTSS)[expressionTables.chrX[["counts"]]$strand == "-"]
  colnames(chrXgeneRanges) = c("start", "end")
  rownames(chrXgeneRanges) = expressionTables.chrX[["counts"]]$geneId
  chrXgeneRanges = as.data.frame(chrXgeneRanges)
  # all(chrXgeneRanges$stop - chrXgeneRanges$start > 0)
  
  # Tally peaks therein
  genePeakCounts = list()
  for (theDataType in dataTypes) {
    genePeakCounts[[theDataType]] = list()
    for(theSample in names(peakData[[theDataType]])) {
      genePeakCounts[[theDataType]][[theSample]] = list()
      for(theAllele in names(peakData[[theDataType]][[theSample]])) {
        print(c(theDataType, theSample, theAllele))
        thisPeakData = peakData[[theDataType]][[theSample]][[theAllele]]
        midPoints = as.numeric(as.character(thisPeakData$start)) + (as.numeric(as.character(thisPeakData$end)) - as.numeric(as.character(thisPeakData$start)))/2
        genePeakCountsPerGene = c()
        for( i in 1:nrow(chrXgeneRanges)) {
          peaksInGeneBoolIdx = midPoints >= chrXgeneRanges$start[i] & midPoints <= chrXgeneRanges$end[i]
          genePeakCountsPerGene[i] = sum(peaksInGeneBoolIdx, na.rm=TRUE) / ((chrXgeneRanges$end[i]-chrXgeneRanges$start[i]+1)/1000)
        }
        names(genePeakCountsPerGene) = rownames(chrXgeneRanges)
        genePeakCounts[[theDataType]][[theSample]][[theAllele]] = genePeakCountsPerGene   
      }
    }
  }
  
  # pool Xi-associated data
  for (theDataType in dataTypes) {
    for(theSample in names(genePeakCounts[[theDataType]])) {
      thisPeakData.ambig = genePeakCounts[[theDataType]][[theSample]][["Ambig"]]
      thisPeakData.xi = genePeakCounts[[theDataType]][[theSample]][["Xi"]]
      genePeakCounts[[theDataType]][[theSample]][["Xi-assoc"]] = thisPeakData.ambig + thisPeakData.xi
    }
  }

  # pool Xa-associated data
  for (theDataType in dataTypes) {
    for(theSample in names(genePeakCounts[[theDataType]])) {
      thisPeakData.ambig = genePeakCounts[[theDataType]][[theSample]][["Ambig"]]
      thisPeakData.xi = genePeakCounts[[theDataType]][[theSample]][["Xa"]]
      genePeakCounts[[theDataType]][[theSample]][["Xa-assoc"]] = thisPeakData.ambig + thisPeakData.xi
    }
  }
  
  # Scale the peak counts
  genePeakCounts.scaled = list()
  for (theDataType in dataTypes) {
    genePeakCounts.scaled[[theDataType]] = list()
    for(theSample in names(genePeakCounts[[theDataType]])) {
      genePeakCounts.scaled[[theDataType]][[theSample]] = list()
      for(theAllele in names(genePeakCounts[[theDataType]][[theSample]])) {
        genePeakCounts.scaled[[theDataType]][[theSample]][[theAllele]] = genePeakCounts[[theDataType]][[theSample]][[theAllele]] * PeakRatioScalingFactors[[theDataType]][theSample]
      }
    }
  }
  
  head(genePeakCounts.scaled[[theDataType]][[theSample]][[theAllele]])
  
  #########################
  # DE vs non-DE  genes
  comp = "WT-vs-Del"
  cellLines = strsplit(comp, "-vs-")[[1]]
  
  currTable = DEtables.chrX[[comp]][["Xi"]]
  names(currTable) = c("geneId", "basemean", "log2FoldChange", 
                       "lfcSE", "stat", "pvalue", "padj")
  geneIdXref = match(currTable$geneId, expressionTables.chrX[["refTPM"]]$geneId)
  # theLoci = (expressionTables.chrX[["refTPM"]]$start + (expressionTables.chrX[["refTPM"]]$stop - expressionTables.chrX[["refTPM"]]$start)/2)[geneIdXref]
  all(names(genePeakCounts.scaled[[theDataType]][["Del"]][[allele]][geneIdXref]) == currTable$geneId)
  
  theLogRatio = -1 * currTable$log2FoldChange
  names(theLogRatio) = currTable$geneId
  thePvalue = currTable$pvalue
  thePadj = currTable$padj
  DEboolIdx = thePadj<=0.05 & abs(theLogRatio) >= 0.5
  
  
  #####
  # Boxplots of density
  
  boxplotsDir = file.path(outputDir, "genebodyPeakDensityBoxplots")
  if ( ! file.exists( boxplotsDir ) ) {
    if ( ! dir.create( boxplotsDir, showWarnings=FALSE, recursive=TRUE )  ) boxplotsDir = outputDir
  }

  dataDesc = paste("genebodyPeakDensity_scaled.TSS-", upstreamTSSbp/1000, "kbToTTS.", comp, sep="")
  pdf(file.path(boxplotsDir, paste(dataDesc,".pdf",sep="")), width=2*5, height=2*5)
  
  par(mfrow=c(2, 2), mar=c(2,6,6,4), oma=c(3,0,2,0))
  
  for(theDataType in rev(dataTypes)) {
    for (theAllele in c("Xi-assoc", "Xa-assoc")) {
      bpYlim = max( c(genePeakCounts.scaled[[theDataType]][[cellLines[2]]][[theAllele]][geneIdXref],
                      genePeakCounts.scaled[[theDataType]][[cellLines[1]]][[theAllele]][geneIdXref]), na.rm=TRUE)
      bpx = boxplot(genePeakCounts.scaled[[theDataType]][[cellLines[2]]][[theAllele]][geneIdXref][DEboolIdx],
              genePeakCounts.scaled[[theDataType]][[cellLines[1]]][[theAllele]][geneIdXref][DEboolIdx],
              genePeakCounts.scaled[[theDataType]][[cellLines[2]]][[theAllele]][geneIdXref][!DEboolIdx],
              genePeakCounts.scaled[[theDataType]][[cellLines[1]]][[theAllele]][geneIdXref][!DEboolIdx],
              names=c(cellLines[2], cellLines[1], cellLines[2], cellLines[1]), col=c(peakCols[[cellLines[2]]], peakCols[[cellLines[1]]]),
              ylab="gene body peak density",
              ylim=c(0,bpYlim))
      abline(v=2.5, lwd=1)
  
      DEtest2 = wilcox.test(genePeakCounts.scaled[[theDataType]][[cellLines[2]]][[theAllele]][geneIdXref][DEboolIdx],
                  genePeakCounts.scaled[[theDataType]][[cellLines[1]]][[theAllele]][geneIdXref][DEboolIdx],
                  alternative="two")
      nonDEtest2 = wilcox.test(genePeakCounts.scaled[[theDataType]][[cellLines[2]]][[theAllele]][geneIdXref][!DEboolIdx],
                           genePeakCounts.scaled[[theDataType]][[cellLines[1]]][[theAllele]][geneIdXref][!DEboolIdx],
                           alternative="two")
      DEtest1 = wilcox.test(genePeakCounts.scaled[[theDataType]][[cellLines[2]]][[theAllele]][geneIdXref][DEboolIdx],
                           genePeakCounts.scaled[[theDataType]][[cellLines[1]]][[theAllele]][geneIdXref][DEboolIdx],
                           alternative="greater")
      nonDEtest1 = wilcox.test(genePeakCounts.scaled[[theDataType]][[cellLines[2]]][[theAllele]][geneIdXref][!DEboolIdx],
                              genePeakCounts.scaled[[theDataType]][[cellLines[1]]][[theAllele]][geneIdXref][!DEboolIdx],
                              alternative="greater")
      par(xpd=TRUE)    
      text(1.5, bpYlim*1.1, paste("DE genes\n(n=", sum(DEboolIdx, na.rm=TRUE), 
                                  "; p2=",  round(DEtest2$p.value, 2 ), 
                                  "; p1=",  round(DEtest1$p.value, 2 ), 
                                  ")", sep="" ), cex=0.95)
      text(3.5, bpYlim*1.1, paste("non DE genes\n(n=", sum(!DEboolIdx, na.rm=TRUE),
                                  "; p2=",  round(nonDEtest2$p.value, 2 ), 
                                  "; p1=",  round(nonDEtest1$p.value, 2 ), 
                                  ")", sep="" ), cex=0.95)
      par(xpd=FALSE)    
      
      title(main=paste(theAllele, " ", theDataType, sep=""), line=3, cex=1)
    }
  }
  
  title(main=paste("Peak counts per 1kb of gene bodies (TSS-", upstreamTSSbp/1000, "kb to TTS): ", comp, " DE vs non-DE genes\n",
                   "normalized by Xa peak ratios; unpaired Wilcoxon bidirectional (p2) an unidirectional testes (p1)" , sep=""), outer=TRUE, cex=2.5)
  dev.off()

  
  #####
  # Boxplots of Differential density
  
  boxplotsDir = file.path(outputDir, "genebodyPeakDifferentialDensityBoxplots")
  if ( ! file.exists( boxplotsDir ) ) {
    if ( ! dir.create( boxplotsDir, showWarnings=FALSE, recursive=TRUE )  ) boxplotsDir = outputDir
  }
  
  dataDesc = paste("genebodyPeakDensity_scaled.TSS-", upstreamTSSbp/1000, "kbToTTS.", comp, sep="")
  pdf(file.path(boxplotsDir, paste(dataDesc,".pdf",sep="")), width=2*5, height=2*5)
  
  par(mfrow=c(2, 2), mar=c(2,6,6,4), oma=c(3,0,2,0))
  
  for(theDataType in rev(dataTypes)) {
    for (theAllele in c("Xi-assoc", "Xa-assoc")) {
      bpYlims = range( c(genePeakCounts.scaled[[theDataType]][[cellLines[2]]][[theAllele]][geneIdXref][which(DEboolIdx)] - genePeakCounts.scaled[[theDataType]][[cellLines[1]]][[theAllele]][geneIdXref][which(DEboolIdx)],
                      genePeakCounts.scaled[[theDataType]][[cellLines[2]]][[theAllele]][geneIdXref][which(!DEboolIdx)] - genePeakCounts.scaled[[theDataType]][[cellLines[1]]][[theAllele]][geneIdXref][which(!DEboolIdx)]), na.rm=TRUE)
      bpYlims = c(-abs(bpYlims[which.max(abs(bpYlims))]),
                  abs(bpYlims[which.max(abs(bpYlims))]))
      bpx = boxplot(genePeakCounts.scaled[[theDataType]][[cellLines[2]]][[theAllele]][geneIdXref][which(DEboolIdx)] - genePeakCounts.scaled[[theDataType]][[cellLines[1]]][[theAllele]][geneIdXref][which(DEboolIdx)],
                    genePeakCounts.scaled[[theDataType]][[cellLines[2]]][[theAllele]][geneIdXref][which(!DEboolIdx)] - genePeakCounts.scaled[[theDataType]][[cellLines[1]]][[theAllele]][geneIdXref][which(!DEboolIdx)],
                    names=c("DE genes", "non-DE genes"), 
                    col=c("green", "gray"),
                    ylab="Differential gene body peak density",
                    ylim=bpYlims)
      abline(v=1.5, lwd=1)
      abline(h=0, lwd=1, lty="dashed", col="blue")
      
      DEtest2 = wilcox.test(genePeakCounts.scaled[[theDataType]][[cellLines[2]]][[theAllele]][geneIdXref][DEboolIdx],
                            genePeakCounts.scaled[[theDataType]][[cellLines[1]]][[theAllele]][geneIdXref][DEboolIdx],
                            alternative="two", paired=TRUE)
      nonDEtest2 = wilcox.test(genePeakCounts.scaled[[theDataType]][[cellLines[2]]][[theAllele]][geneIdXref][!DEboolIdx],
                               genePeakCounts.scaled[[theDataType]][[cellLines[1]]][[theAllele]][geneIdXref][!DEboolIdx],
                               alternative="two", paired=TRUE)
      DEtest1 = wilcox.test(genePeakCounts.scaled[[theDataType]][[cellLines[2]]][[theAllele]][geneIdXref][DEboolIdx],
                            genePeakCounts.scaled[[theDataType]][[cellLines[1]]][[theAllele]][geneIdXref][DEboolIdx],
                            alternative="greater", paired=TRUE)
      nonDEtest1 = wilcox.test(genePeakCounts.scaled[[theDataType]][[cellLines[2]]][[theAllele]][geneIdXref][!DEboolIdx],
                               genePeakCounts.scaled[[theDataType]][[cellLines[1]]][[theAllele]][geneIdXref][!DEboolIdx],
                               alternative="greater", paired=TRUE)
      par(xpd=TRUE)    
      text(1, bpYlims[2]*1.2, paste("n=", sum(DEboolIdx, na.rm=TRUE), 
                                  "; p2=",  round(DEtest2$p.value, 2 ), 
                                  "; p1=",  round(DEtest1$p.value, 2 ), 
                                  sep="" ), cex=0.95)
      text(2, bpYlims[2]*1.2, paste("n=", sum(!DEboolIdx, na.rm=TRUE),
                                  "; p2=",  round(nonDEtest2$p.value, 2 ), 
                                  "; p1=",  round(nonDEtest1$p.value, 2 ), 
                                  sep="" ), cex=0.95)
      par(xpd=FALSE)    
      
      title(main=paste(theAllele, " ", theDataType, sep=""), line=3, cex=1)
    }
  }
  
  title(main=paste("Peak counts per 1kb of gene bodies (TSS-", upstreamTSSbp/1000, "kb to TTS): ", comp, " DE vs non-DE genes\n",
                   "normalized by Xa peak ratios; paired Wilcoxon bidirectional (p2) an unidirectional testes (p1)" , sep=""), outer=TRUE, cex=2.5)
  dev.off()

  
  
  #####
  # Boxplots of Differential density - rescaled to exclude outliers
  
  boxplotsDir = file.path(outputDir, "genebodyPeakDifferentialDensityBoxplots.V2.outliersExcluded")
  if ( ! file.exists( boxplotsDir ) ) {
    if ( ! dir.create( boxplotsDir, showWarnings=FALSE, recursive=TRUE )  ) boxplotsDir = outputDir
  }
  
  dataDesc = paste("genebodyPeakDensity_scaled.TSS-", upstreamTSSbp/1000, "kbToTTS.", comp, sep="")
  pdf(file.path(boxplotsDir, paste(dataDesc,".pdf",sep="")), width=2*5, height=2*5)
  
  par(mfrow=c(2, 2), mar=c(2,6,6,4), oma=c(3,0,2,0))
  
  for(theDataType in rev(dataTypes)) {
    for (theAllele in c("Xi-assoc", "Xa-assoc")) {
      bpx = boxplot( genePeakCounts.scaled[[theDataType]][[cellLines[2]]][[theAllele]][geneIdXref][which(DEboolIdx)] - genePeakCounts.scaled[[theDataType]][[cellLines[1]]][[theAllele]][geneIdXref][which(DEboolIdx)],
                    genePeakCounts.scaled[[theDataType]][[cellLines[2]]][[theAllele]][geneIdXref][which(!DEboolIdx)] - genePeakCounts.scaled[[theDataType]][[cellLines[1]]][[theAllele]][geneIdXref][which(!DEboolIdx)],
                    plot=FALSE)
      bpYlims = max(abs(bpx$stats[c(1,5), ]), na.rm=TRUE)
      bpYlims = c(-bpYlims, +bpYlims)
      bpx = boxplot(genePeakCounts.scaled[[theDataType]][[cellLines[2]]][[theAllele]][geneIdXref][which(DEboolIdx)] - genePeakCounts.scaled[[theDataType]][[cellLines[1]]][[theAllele]][geneIdXref][which(DEboolIdx)],
                    genePeakCounts.scaled[[theDataType]][[cellLines[2]]][[theAllele]][geneIdXref][which(!DEboolIdx)] - genePeakCounts.scaled[[theDataType]][[cellLines[1]]][[theAllele]][geneIdXref][which(!DEboolIdx)],
                    names=c("DE genes", "non-DE genes"), 
                    col=c("green", "gray"),
                    ylab="Differential gene body peak density",
                    ylim=bpYlims)
      abline(v=1.5, lwd=1)
      abline(h=0, lwd=1, lty="dashed", col="blue")
      
      DEtest2 = wilcox.test(genePeakCounts.scaled[[theDataType]][[cellLines[2]]][[theAllele]][geneIdXref][DEboolIdx],
                            genePeakCounts.scaled[[theDataType]][[cellLines[1]]][[theAllele]][geneIdXref][DEboolIdx],
                            alternative="two", paired=TRUE)
      nonDEtest2 = wilcox.test(genePeakCounts.scaled[[theDataType]][[cellLines[2]]][[theAllele]][geneIdXref][!DEboolIdx],
                               genePeakCounts.scaled[[theDataType]][[cellLines[1]]][[theAllele]][geneIdXref][!DEboolIdx],
                               alternative="two", paired=TRUE)
      DEtest1 = wilcox.test(genePeakCounts.scaled[[theDataType]][[cellLines[2]]][[theAllele]][geneIdXref][DEboolIdx],
                            genePeakCounts.scaled[[theDataType]][[cellLines[1]]][[theAllele]][geneIdXref][DEboolIdx],
                            alternative="greater", paired=TRUE)
      nonDEtest1 = wilcox.test(genePeakCounts.scaled[[theDataType]][[cellLines[2]]][[theAllele]][geneIdXref][!DEboolIdx],
                               genePeakCounts.scaled[[theDataType]][[cellLines[1]]][[theAllele]][geneIdXref][!DEboolIdx],
                               alternative="greater", paired=TRUE)
      par(xpd=TRUE)    
      text(1, bpYlims[2]*1.2, paste("n=", sum(DEboolIdx, na.rm=TRUE), 
                                    "; p2=",  round(DEtest2$p.value, 2 ), 
                                    "; p1=",  round(DEtest1$p.value, 2 ), 
                                    sep="" ), cex=0.95)
      text(2, bpYlims[2]*1.2, paste("n=", sum(!DEboolIdx, na.rm=TRUE),
                                    "; p2=",  round(nonDEtest2$p.value, 2 ), 
                                    "; p1=",  round(nonDEtest1$p.value, 2 ), 
                                    sep="" ), cex=0.95)
      par(xpd=FALSE)    
      
      title(main=paste(theAllele, " ", theDataType, sep=""), line=3, cex=1)
    }
  }
  
  title(main=paste("Peak counts per 1kb of gene bodies (TSS-", upstreamTSSbp/1000, "kb to TTS): ", comp, " DE vs non-DE genes\n",
                   "normalized by Xa peak ratios; paired Wilcoxon bidirectional (p2) an unidirectional testes (p1)" , sep=""), outer=TRUE, cex=2.5)
  dev.off()
  
}

 
