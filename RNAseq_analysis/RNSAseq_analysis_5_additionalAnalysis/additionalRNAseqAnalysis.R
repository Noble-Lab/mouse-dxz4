###########################################################
###########################################################
#
# 170504
# Scatterplots of WT vs. Del-hinge1a, Del-hinge1b, Del-Dxz4, Del-Ds-TR TPMs
# Overlaid with WT vs Del DE genes
#
# Run on Macbook in R GUI (not studio) and 
# save from Quartz window to get arrows.
#
###########################################################

projectDir = "<projectDirectory>"
subDir = "20170620_RNAseq_InvDxz4abAndWTab"
workDir = "Patski.AllData"

outputDir = file.path(projectDir, subDir, "expressionAnalysis")
setwd(outputDir)
getwd()

#save.image()
#load(file = ".RData") # From .RData

scatterPlotIIdir = file.path(outputDir, "scatterplotsII")
if ( ! file.exists( scatterPlotIIdir ) ) {
  if ( ! dir.create( scatterPlotIIdir, showWarnings=FALSE, recursive=TRUE )  ) scatterPlotIIdir = outputDir
}

# Load expression tables
tableDir = file.path(outputDir, "expressionTables")
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

# WTvsDel DE genes
WTvsDelDEgenes = c(
"Shroom4",
"Clcn5",
"Tspan7",
"Ddx3x",
"Rhox5",
"Fhl1",
"Flna",
"Zfx",
"Eif2s3x",
"Msn",
"Eda2r",
"Pcdh19",
"Wbp5",
"Ngfrap1",
"Ammecr1",
"Mid1"
)
# WTvsDelDEgenesIdx = match(WTvsDelDEgenes, expressionTables.chrX[["TPM"]]$geneId)

# Plot variables
# colnames(expressionTables.chrX[["refTPM"]])
yLines2plot = c("PatskiDel1", "PatskiDel1b", "PatskiDel5", "PatskiDel2")
yLinesLabels = c("Del-hinge a", "Del-hinge b", "Del-Dxz4", "Del-Ds-TR")
yLinesColors = c('red', 'darkred', 'orange', 'green')
# DATASETS=(Del1_Xa Del1_Xi Del5_Xa Del5_Xi Del2_Xa Del2_Xi Inv1_Xa Inv1_Xi WT_Xa WT_Xi )
# DATACOLS=('red' 'red' 'orange' 'orange' 'green' 'green' 'violet' 'violet' 'blue' 'blue'  )
xLine="PatskiWT"

# Plot variables 2
# colnames(expressionTables.chrX[["refTPM"]])
yLines2plot = c("PatskiInvDxz4a", "PatskiInvDxz4b")
yLinesLabels = c("Inv-Dxz4 a", "Inv-Dxz4 b")
yLinesColors = c(gray(0.2), gray(0.2))
# DATASETS=(Del1_Xa Del1_Xi Del5_Xa Del5_Xi Del2_Xa Del2_Xi Inv1_Xa Inv1_Xi WT_Xa WT_Xi )
# DATACOLS=('red' 'red' 'orange' 'orange' 'green' 'green' 'violet' 'violet' 'blue' 'blue'  )
xLine="PatskiWTa"


# 2^-4 = 0.0625; 2^-3 = 0.125; 2^-2 = 0.25
expressionType = "TPM"
lowMeanThreshold = -4
xyMin = -4
xyMax = 10
  
dev.new(width=4*2, height=4*length(yLines2plot))  
par(mfrow=c(length(yLines2plot), 2))
par(mar=c(2,3,2,1), oma=c(2,2,2,2), pty="s")
# # Move labels in: http://www.programmingr.com/content/controlling-margins-and-axes-oma-and-mgp/
# par(mgp = c(0, 1, 0))

for (lineIdx in 1:length(yLines2plot)) {
  # lineIdx = 1; par( mfrow=c(1,2) )
  for (genome in c("ref", "alt")) {
    # lineIdx = 1; genome = "alt"
    
    genomeExpType = paste(genome, expressionType, sep="")
    
    meanTPM = rowMeans(cbind(expressionTables.chrX[[genomeExpType]][,xLine],
                   expressionTables.chrX[[genomeExpType]][, yLines2plot[lineIdx] ]), na.rm=TRUE)
    lowMeanIdex = log2(meanTPM) < lowMeanThreshold | is.nan(meanTPM)
  
    x = log2(expressionTables.chrX[[genomeExpType]][,xLine])[!lowMeanIdex]
    y = log2(expressionTables.chrX[[genomeExpType]][, yLines2plot[lineIdx]])[!lowMeanIdex]
    WTvsDelDEgenesIdx = match(WTvsDelDEgenes, expressionTables.chrX[[expressionType]]$geneId[!lowMeanIdex])
    
    xLow = x < xyMin
    x[xLow] = xyMin
    xHigh = x > xyMax
    x[xHigh] = xyMax
    
    yLow = y < xyMin
    y[yLow] = xyMin
    yHigh = y > xyMax
    y[yHigh] = xyMax
    
    pchSymbols = rep(18, length(x))
    pchSymbols[xLow] = -as.hexmode("25C4") # <
    pchSymbols[xHigh] = -as.hexmode("25BA") # >
    pchSymbols[yLow] = -as.hexmode("25BC") # v
    pchSymbols[yHigh] = -as.hexmode("25B2") # ^
    
    colSymbols = rep(yLinesColors[lineIdx], length(x))
    colSymbols[xLow] = "blue"
    colSymbols[xHigh] = "blue"
    colSymbols[yLow] = "blue"
    colSymbols[yHigh] = "blue"
    
    cexSymbols = rep(1.5, length(x))
    cexSymbols[xLow] = 0.5
    cexSymbols[xHigh] = 0.5
    cexSymbols[yLow] = 0.5
    cexSymbols[yHigh] = 0.5

    plot(x, y, pch=pchSymbols, col=colSymbols, cex=cexSymbols,
           xlim=c(xyMin, xyMax), ylim=c(xyMin, xyMax),
           xlab="", ylab="", main="", axes=FALSE)
           
    axis(1, pos=0)
    axis(2, pos=0,las=1) 
        
    if (lineIdx == 1) {
      title(main=ifelse(genome=="ref", "Xi", "Xa"), cex.main=2, line=0.75)
    }
    
    if (lineIdx == length(yLines2plot)) {
      title(xlab=paste("WT [log2(", expressionType, ")]", sep=""), cex.lab=1.5, line=1)
    }
    
    if (genome == "ref") {
      title(ylab=paste(yLinesLabels[lineIdx], " [log2(", expressionType, ")]", sep=""), cex.lab=1.5, line=1)
    }
    
    # Xi DE genes 
    if (genome == "ref") {
      points(x[WTvsDelDEgenesIdx], y[WTvsDelDEgenesIdx], 
           pch=pchSymbols[WTvsDelDEgenesIdx], 
           cex=cexSymbols[WTvsDelDEgenesIdx]*2.5,
           col="black")
      points(x[WTvsDelDEgenesIdx], y[WTvsDelDEgenesIdx], 
           pch=pchSymbols[WTvsDelDEgenesIdx], 
           cex=cexSymbols[WTvsDelDEgenesIdx]*2, 
           col="white")
      points(x[WTvsDelDEgenesIdx], y[WTvsDelDEgenesIdx], 
      	   pch=pchSymbols[WTvsDelDEgenesIdx], 
      	   cex=cexSymbols[WTvsDelDEgenesIdx], 
      	   col=colSymbols[WTvsDelDEgenesIdx])
    }
    
    # log2FC=0.5 lines;  2^0.5 = 1.414214
    #par(xpd=FALSE) 
    abline(a=0.5, b=1, lt="dotted", lw=0.5)
    abline(a=-0.5, b=1, lt="dotted", lw=0.5)

  }
}

# dev.off()












###########################################################
###########################################################
# 
# 170628
# Added WT vs. Inv-Dxz4a and b
# Xi and Xa printed seperately
#
# 170504
#
# Scatterplots of WT vs. Del-hinge1a, Del-hinge1b, Del-Dxz4, Del-Ds-TR TPMs
# Overlaid with WT vs Del DE genes
#
# Run on Macbook in R GUI (not studio) and 
# save from Quartz window to get arrows.
#
# Must first generate chrX-specific expression tables using code in script:
# 'RNAseq_HTSeqAndExpressionAnalysis_workflow.sh'
#
###########################################################

projectDir = "<projectDirectory>"
subDir = "20170620_RNAseq_InvDxz4abAndWTab"
workDir = "Patski.AllData"

outputDir = file.path(projectDir, subDir, "expressionAnalysis")
setwd(outputDir)
getwd()

#save.image()
#load(file = ".RData") # From .RData

scatterPlotIIdir = file.path(outputDir, "scatterplotsII")
if ( ! file.exists( scatterPlotIIdir ) ) {
  if ( ! dir.create( scatterPlotIIdir, showWarnings=FALSE, recursive=TRUE )  ) scatterPlotIIdir = outputDir
}

# Load expression tables
tableDir = file.path(outputDir, "expressionTables")
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

# WTvsDel DE genes
WTvsDelDEgenes = c(
  "Shroom4",
  "Clcn5",
  "Tspan7",
  "Ddx3x",
  "Rhox5",
  "Fhl1",
  "Flna",
  "Zfx",
  "Eif2s3x",
  "Msn",
  "Eda2r",
  "Pcdh19",
  "Wbp5",
  "Ngfrap1",
  "Ammecr1",
  "Mid1"
)
# WTvsDelDEgenesIdx = match(WTvsDelDEgenes, expressionTables.chrX[["TPM"]]$geneId)

# # Plot variables
# # colnames(expressionTables.chrX[["refTPM"]])
# yLines2plot = c("PatskiDel1", "PatskiDel1b", "PatskiDel5", "PatskiDel2", "PatskiInvDxz4a", "PatskiInvDxz4b")
# yLinesLabels = c("Del-hinge a", "Del-hinge b", "Del-Dxz4", "Del-Ds-TR", "Inv-Dxz4 a", "Inv-Dxz4 b")
# yLinesColors = c('red', 'darkred', 'orange', 'green', gray(0.2), gray(0.2))
# # DATASETS=(Del1_Xa Del1_Xi Del5_Xa Del5_Xi Del2_Xa Del2_Xi Inv1_Xa Inv1_Xi WT_Xa WT_Xi )
# # DATACOLS=('red' 'red' 'orange' 'orange' 'green' 'green' 'violet' 'violet' 'blue' 'blue'  )

# Plot variables
# colnames(expressionTables.chrX[["refTPM"]])
yLines2plot = c("PatskiInvDxz4a", "PatskiInvDxz4b")
yLinesLabels = c("Inv-Dxz4 a", "Inv-Dxz4 b")
yLinesColors = c(gray(0.2), gray(0.2))
# DATASETS=(Del1_Xa Del1_Xi Del5_Xa Del5_Xi Del2_Xa Del2_Xi Inv1_Xa Inv1_Xi WT_Xa WT_Xi )
# DATACOLS=('red' 'red' 'orange' 'orange' 'green' 'green' 'violet' 'violet' 'blue' 'blue'  )

# 2^-4 = 0.0625; 2^-3 = 0.125; 2^-2 = 0.25
expressionType = "TPM"
lowMeanThreshold = -4
xyMin = -4
xyMax = 10

for (WTline in grep("WT", colnames(expressionTables.chrX[["refTPM"]]), value=TRUE)) {
  # [1] "PatskiWT"             "PatskiWT_0uMaza_rep1" "PatskiWT_0uMaza_rep2" "PatskiWT_4uMaza_rep1"
  # [5] "PatskiWT_4uMaza_rep2" "PatskiWT_DelFirre"    "PatskiWTa"            "PatskiWTb" 
  # WTline = "PatskiWT"
  # WTline = "PatskiWTa"
   
  print(WTline)
  
  for (genome in c("ref", "alt")) {
    # lineIdx = 1; genome = "alt"
    #genome = "ref"  # Xi
    #genome = "alt" # Xa
    
    print(genome)
          
	dev.new(width=4*2, height=4*length(yLines2plot))  
	par(mfrow=c(length(yLines2plot), 2))
	par(mar=c(2,3,2,1), oma=c(2,2,2,2), pty="s")

    for (lineIdx in 1:length(yLines2plot)) {
      # lineIdx = 1; par( mfrow=c(1,2) )
        
      genomeExpType = paste(genome, expressionType, sep="")
      
      meanTPM = rowMeans(cbind(expressionTables.chrX[[genomeExpType]][,WTline],
                               expressionTables.chrX[[genomeExpType]][, yLines2plot[lineIdx] ]), na.rm=TRUE)
      lowMeanIdex = log2(meanTPM) < lowMeanThreshold | is.nan(meanTPM)
      
      x = log2(expressionTables.chrX[[genomeExpType]][,WTline])[!lowMeanIdex]
      y = log2(expressionTables.chrX[[genomeExpType]][, yLines2plot[lineIdx]])[!lowMeanIdex]
      WTvsDelDEgenesIdx = match(WTvsDelDEgenes, expressionTables.chrX[[expressionType]]$geneId[!lowMeanIdex])
      
      xLow = x < xyMin
      x[xLow] = xyMin
      xHigh = x > xyMax
      x[xHigh] = xyMax
      
      yLow = y < xyMin
      y[yLow] = xyMin
      yHigh = y > xyMax
      y[yHigh] = xyMax
      
      pchSymbols = rep(18, length(x))
      pchSymbols[xLow] = -as.hexmode("25C4") # <
      pchSymbols[xHigh] = -as.hexmode("25BA") # >
      pchSymbols[yLow] = -as.hexmode("25BC") # v
      pchSymbols[yHigh] = -as.hexmode("25B2") # ^
      
      colSymbols = rep(yLinesColors[lineIdx], length(x))
      colSymbols[xLow] = "blue"
      colSymbols[xHigh] = "blue"
      colSymbols[yLow] = "blue"
      colSymbols[yHigh] = "blue"
      
      cexSymbols = rep(1.5, length(x))
      cexSymbols[xLow] = 0.5
      cexSymbols[xHigh] = 0.5
      cexSymbols[yLow] = 0.5
      cexSymbols[yHigh] = 0.5
      
      plot(x, y, pch=pchSymbols, col=colSymbols, cex=cexSymbols,
           xlim=c(xyMin, xyMax), ylim=c(xyMin, xyMax),
           xlab="", ylab="", main="", axes=FALSE)
      
      axis(1, pos=0)
      axis(2, pos=0,las=1) 
      
      if (lineIdx == 1) {
        title(main=ifelse(genome=="ref", "Xi", "Xa"), cex.main=2, line=-2, outer=TRUE)
        title(main=WTline, cex.main=2, line=0, outer=TRUE)
      }
      
      if (lineIdx == length(yLines2plot)) {
        title(xlab=paste("WT [log2(", expressionType, ")]", sep=""), cex.lab=1.5, line=1)
      }
      
      #if (genome == "ref") {
        title(ylab=paste(yLinesLabels[lineIdx], " [log2(", expressionType, ")]", sep=""), cex.lab=1.5, line=1)
      #}
      
      # Xi DE genes 
      if (genome == "ref") {
        points(x[WTvsDelDEgenesIdx], y[WTvsDelDEgenesIdx], 
               pch=pchSymbols[WTvsDelDEgenesIdx], 
               cex=cexSymbols[WTvsDelDEgenesIdx]*2.5,
               col="black")
        points(x[WTvsDelDEgenesIdx], y[WTvsDelDEgenesIdx], 
               pch=pchSymbols[WTvsDelDEgenesIdx], 
               cex=cexSymbols[WTvsDelDEgenesIdx]*2, 
               col="white")
        points(x[WTvsDelDEgenesIdx], y[WTvsDelDEgenesIdx], 
               pch=pchSymbols[WTvsDelDEgenesIdx], 
               cex=cexSymbols[WTvsDelDEgenesIdx], 
               col=colSymbols[WTvsDelDEgenesIdx])
      }
      
      # log2FC=0.5 lines;  2^0.5 = 1.414214
      #par(xpd=FALSE) 
      abline(a=0.5, b=1, lt="dotted", lw=0.5)
      abline(a=-0.5, b=1, lt="dotted", lw=0.5)
    }
    
    # dev.off()
  }
}  










###########################################################
###########################################################
# 
# 171202
# Xi Proportion (i.e. Xi/(Xi+Xa)) and WT-vs-mutant log2FC along chrX
#
# Must first generate chrX-specific expression tables using code in script:
# 'RNAseq_HTSeqAndExpressionAnalysis_workflow.sh'
# 
# Xi Proportion Notes: 
# Based on TPM values in expression tables.
# Take the average TPM across the samples.
# Thought about using DESeq2 normalized counts, but those are normalized only within the allelic dataset (Sp or B6).
# Would've also had to save the normalized counts, which I didn't do.
#
# Log2FC Notes:
# Used the DESeq2 log2 fold change values, so...
# Must first generate chrX-specific DEseq result tables (which contain log fold change) using code in this script above
#
###########################################################

projectDir = "<projectDirectory>"
subDir = "20170620_RNAseq_InvDxz4abAndWTab"
workDir = "Patski.AllData"

outputDir = file.path(projectDir, subDir, "expressionAnalysis")
setwd(outputDir)
getwd()

#save.image()
#load(file = ".RData") # From .RData

mainChrXplotDir = file.path(outputDir, "chrXplots")
if ( ! file.exists( mainChrXplotDir ) ) {
  if ( ! dir.create( mainChrXplotDir, showWarnings=FALSE, recursive=TRUE )  ) mainChrXplotDir = outputDir
}


################################
# Functions

# Add an alpha value to a colour
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

# Useful functions
NArowSums = function(x){
  ifelse(all(is.na(x)), NaN, sum(x, na.rm=TRUE))
} 

NArowMeans = function(x){
  ifelse(all(is.na(x)), NaN, mean(x, na.rm=TRUE))
} 

################################
# Data

###############
# Load DE result tables

projectDir = "<projectDirectory>"
subDir = "20170220_RNAseq_reprocessAllLibs"

altOutputDir = file.path(projectDir, subDir, "expressionAnalysis")

DEresultDirs = list()
DEresultDirs[["WT-vs-Del"]] = file.path(altOutputDir, "DESeq2/WT-vs-Del1a")
DEresultDirs[["WT-vs-InvDxz4"]] = file.path(outputDir, "DESeq2/WT-vs-InvDxz4")
DEresultDirs[["newWT-vs-InvDxz4"]] = file.path(outputDir, "DESeq2/newWT-vs-InvDxz4")

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


DEtables.chrX = list()
for (comp in names(DEresultFiles)) {
  DEtables.chrX[[comp]] = list()
  for (allele in names(DEresultFiles[[comp]])) {
    tableName = DEresultFiles[[comp]][[allele]]
    print(tableName)
    DEtables.chrX[[comp]][[allele]] = read.table(tableName, sep = "\t", header=TRUE)
  }
}


###############
# Load expression tables
tableDir = file.path(outputDir, "expressionTables")
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

# Get avereage expression levels
expressionTables.chrX.average = list()
for (genomeType in names(expressionTables.chrX)) {
  # genomeType="counts"
  cat(genomeType, "\n", sep="")
  if (grepl("counts", genomeType, ignore.case=TRUE)) {
    expressionTables.chrX.average[[genomeType]] = data.frame(expressionTables.chrX[[genomeType]][,c("geneId", "chr", "start", "stop", "PatskiWT")])
    colnames(expressionTables.chrX.average[[genomeType]]) = c("geneId", "chr", "start", "stop", "PatskiWT")
    expressionTables.chrX.average[[genomeType]]$PatskiWT = apply(expressionTables.chrX[[genomeType]][, c("PatskiWT", "PatskiWT_0uMaza_rep1", "PatskiWT_0uMaza_rep2")], 1, NArowSums)
    expressionTables.chrX.average[[genomeType]]$PatskiDel = apply(expressionTables.chrX[[genomeType]][, c("PatskiDel1", "PatskiDel1_0uMaza_rep1", "PatskiDel1_0uMaza_rep2")], 1, NArowSums)
    expressionTables.chrX.average[[genomeType]]$PatskiNewWT = apply(expressionTables.chrX[[genomeType]][, c("PatskiWTa", "PatskiWTb")], 1, NArowSums)
    expressionTables.chrX.average[[genomeType]]$PatskiInvDxz4 = apply(expressionTables.chrX[[genomeType]][, c("PatskiInvDxz4a", "PatskiInvDxz4b")], 1, NArowSums)
  } else {
    expressionTables.chrX.average[[genomeType]] = data.frame(expressionTables.chrX[[genomeType]][,c("geneId", "chr", "start", "stop", "PatskiWT")])
    colnames(expressionTables.chrX.average[[genomeType]]) = c("geneId", "chr", "start", "stop", "PatskiWT")
    expressionTables.chrX.average[[genomeType]]$PatskiWT = apply(expressionTables.chrX[[genomeType]][, c("PatskiWT", "PatskiWT_0uMaza_rep1", "PatskiWT_0uMaza_rep2")], 1, NArowMeans)
    expressionTables.chrX.average[[genomeType]]$PatskiDel = apply(expressionTables.chrX[[genomeType]][, c("PatskiDel1", "PatskiDel1_0uMaza_rep1", "PatskiDel1_0uMaza_rep2")], 1, NArowMeans)
    expressionTables.chrX.average[[genomeType]]$PatskiNewWT = apply(expressionTables.chrX[[genomeType]][, c("PatskiWTa", "PatskiWTb")], 1, NArowSums)
    expressionTables.chrX.average[[genomeType]]$PatskiInvDxz4 = apply(expressionTables.chrX[[genomeType]][, c("PatskiInvDxz4a", "PatskiInvDxz4b")], 1, NArowMeans)
  }
}

# WTvsDel DE genes
WTvsDelDEgenes = c(
  "Shroom4",
  "Clcn5",
  "Tspan7",
  "Ddx3x",
  "Rhox5",
  "Fhl1",
  "Flna",
  "Zfx",
  "Eif2s3x",
  "Msn",
  "Eda2r",
  "Pcdh19",
  "Wbp5",
  "Ngfrap1",
  "Ammecr1",
  "Mid1"
)
# WTvsDelDEgenesIdx = match(WTvsDelDEgenes, expressionTables.chrX[["TPM"]]$geneId)

newEscapeGenes = c(
  "Xist",
  "Zfp300",
  "Pbdc1",
  "Ftx",
  "Eif2s3x",
  "Kdm6a",
  "Kdm5c",
  "5530601H04Rik",
  "Xkrx",
  "Atp6ap1",
  "Rp2h",
  "Ddx3x",
  "1810030O07Rik",
  "Mid1",
  "Shroom4",
  "Gm14827",
  "Firre",
  "Slc16a2",
  "Car5b",
  "Magt1",
  "Tmem29",
  "Efnb1",
  "Pim2",
  "Hmgn5",
  "Gnl3l",
  "Yipf6",
  "Asmt",
  "Sh3bgrl",
  "Rlim"
)


#################################################################
# Plot Xi/(Xi+Xa) Proportions

chrXplotDir = file.path(mainChrXplotDir, "XiTPMproportion")
if ( ! file.exists( chrXplotDir ) ) {
  if ( ! dir.create( chrXplotDir, showWarnings=FALSE, recursive=TRUE )  ) chrXplotDir = mainChrXplotDir
}

# Plot Xi/(Xi+Xa)
sampleNames = c("PatskiWT", "PatskiDel", "PatskiNewWT", "PatskiInvDxz4")

sampleCols = c(add.alpha("blue", alpha=1),
               add.alpha("red", alpha=1),
               add.alpha("blue", alpha=1),
               add.alpha("black", alpha=1))
names(sampleCols) = sampleNames

escapeeCol=add.alpha("yellow", alpha=1)

samplePch = c(21, 24, 21, 23)
names(samplePch) = sampleNames

samplePchEscapee = c(20, 17, 20, 18)
names(samplePchEscapee) = sampleNames

samplePchEscapeeCex = c(1, .75, 1, 1)
names(samplePchEscapeeCex) = sampleNames

theLoci = expressionTables.chrX.average[["refTPM"]]$start + (expressionTables.chrX.average[["refTPM"]]$stop - expressionTables.chrX.average[["refTPM"]]$start)/2
theWTXiTPMprop = expressionTables.chrX.average[["refTPM"]]$PatskiWT / ( expressionTables.chrX.average[["refTPM"]]$PatskiWT  + expressionTables.chrX.average[["altTPM"]]$PatskiWT ) 
theDelXiTPMprop = expressionTables.chrX.average[["refTPM"]]$PatskiDel / ( expressionTables.chrX.average[["refTPM"]]$PatskiDel + expressionTables.chrX.average[["altTPM"]]$PatskiDel ) 
theNewWTXiTPMprop = expressionTables.chrX.average[["refTPM"]]$PatskiNewWT / ( expressionTables.chrX.average[["refTPM"]]$PatskiNewWT  + expressionTables.chrX.average[["altTPM"]]$PatskiNewWT ) 
theInvDxz4XiTPMprop = expressionTables.chrX.average[["refTPM"]]$PatskiInvDxz4 / ( expressionTables.chrX.average[["refTPM"]]$PatskiInvDxz4 + expressionTables.chrX.average[["altTPM"]]$PatskiInvDxz4 ) 

names(theWTXiTPMprop) = names(theDelXiTPMprop) = names(theNewWTXiTPMprop) = names(theInvDxz4XiTPMprop) = expressionTables.chrX.average[["refTPM"]]$geneId

escapeBoolIdx = expressionTables.chrX.average[["refTPM"]]$geneId %in% newEscapeGenes

lowWTXiExprBoolIdx = expressionTables.chrX.average[["refTPM"]]$PatskiWT < 0.15 |
  (expressionTables.chrX.average[["refTPM"]]$PatskiWT + expressionTables.chrX.average[["altTPM"]]$PatskiWT) < 1 |
  theWTXiTPMprop < 0.015 |
  expressionTables.chrX.average[["refTPM"]]$geneId %in% c("Gm6568", "Xlr")
theWTXiTPMpropFixed = theWTXiTPMprop
theWTXiTPMpropFixed[lowWTXiExprBoolIdx] = 0

lowDelXiExprBoolIdx = expressionTables.chrX.average[["refTPM"]]$PatskiDel < 0.15 |
  (expressionTables.chrX.average[["refTPM"]]$PatskiDel + expressionTables.chrX.average[["altTPM"]]$PatskiDel) < 1 |
  theDelXiTPMprop < 0.015 |
  expressionTables.chrX.average[["refTPM"]]$geneId %in% c("Gm6568", "Xlr")
theDelXiTPMpropFixed = theDelXiTPMprop
theDelXiTPMpropFixed[lowDelXiExprBoolIdx] = 0

lowNewWTXiExprBoolIdx = expressionTables.chrX.average[["refTPM"]]$PatskiNewWT < 0.15 |
  (expressionTables.chrX.average[["refTPM"]]$PatskiNewWT + expressionTables.chrX.average[["altTPM"]]$PatskiNewWT) < 1 |
  theNewWTXiTPMprop < 0.015 |
  expressionTables.chrX.average[["refTPM"]]$geneId %in% c("Gm6568", "Xlr")
theNewWTXiTPMpropFixed = theWTXiTPMprop
theNewWTXiTPMpropFixed[lowWTXiExprBoolIdx] = 0

lowInvDxz4XiExprBoolIdx = expressionTables.chrX.average[["refTPM"]]$PatskiInvDxz4 < 0.15 |
  (expressionTables.chrX.average[["refTPM"]]$PatskiInvDxz4 + expressionTables.chrX.average[["altTPM"]]$PatskiInvDxz4) < 1 |
  theInvDxz4XiTPMprop < 0.015 |
  expressionTables.chrX.average[["refTPM"]]$geneId %in% c("Gm6568", "Xlr")
theInvDxz4XiTPMpropFixed = theInvDxz4XiTPMprop
theInvDxz4XiTPMpropFixed[lowInvDxz4XiExprBoolIdx] = 0

theYlims=c(0,1)
theXlims=c(0, 180*10^6)

firrePos = (50563120+(abs(50563120-50635321)/2))
dxz4Pos = (75764754+(abs(75637518-75764754)/2))
xistPos = (103460373+(abs(103460373-103483233)/2))


################################
# Plot Xi proportions for both WT And Del

dataDesc = "XiProportions_WT-and-Del"
pdf(file.path(chrXplotDir, paste(dataDesc,".pdf",sep="")), width=15, height=5)

plot(theLoci, theWTXiTPMpropFixed, 
     pch=samplePch["PatskiWT"], col=sampleCols["PatskiWT"], lwd=2,
     axes=F, main=NA, xlab=NA, ylab=NA, ylim=theYlims, xlim=theXlims)
points(theLoci[escapeBoolIdx], theWTXiTPMpropFixed[escapeBoolIdx], 
       pch=samplePchEscapee["PatskiWT"], cex=samplePchEscapeeCex["PatskiWT"], col="yellow")

points(theLoci, theDelXiTPMpropFixed,
       pch=samplePch["PatskiDel"], col=sampleCols["PatskiDel"], lwd=2)
points(theLoci[escapeBoolIdx], theDelXiTPMpropFixed[escapeBoolIdx], 
       pch=samplePchEscapee["PatskiDel"], cex=samplePchEscapeeCex["PatskiDel"], col="yellow")

# hiDelBoolIdx = theDelXiTPMpropFixed[!escapeBoolIdx] > 0.05
# names(theDelXiTPMpropFixed[!escapeBoolIdx])[which(hiDelBoolIdx)]

segments(theLoci, theWTXiTPMpropFixed, theLoci, theDelXiTPMpropFixed, col=sampleCols["PatskiDel"], lty="solid", lwd=0.5)

axis(side=2)
mtext(side=2, line=2.5, text="Xi/(Xi+Xa)")
axis(1, at=seq(0, 180*10^6, 20*10^6), labels=seq(0, 180, 20))
title(main="Xi expression proportion", line=3)

legend("topright", bty="n",
       legend=c(paste("WT genes (", sum(!is.na(theWTXiTPMpropFixed)), ")", sep=""),
                paste("Del genes (", sum(!is.na(theDelXiTPMpropFixed)), ")", sep=""),
                paste("Escape genes (", sum(escapeBoolIdx), ")", sep="" )),
       pch=c(samplePch["PatskiWT"], samplePch["PatskiDel"], 20),
       col=c(sampleCols["PatskiWT"], sampleCols["PatskiDel"], "yellow")
       )

# LOI
abline(v=firrePos, lty="dashed", lw=0.5, col="black")
abline(v=dxz4Pos, lty="dashed", lw=0.5, col="black")
abline(v=xistPos, lty="dashed", lw=0.5, col="black")
axis(3, at=c(firrePos, dxz4Pos, xistPos), labels=c("Firre", "Dxz4", "Xist"), tick = FALSE)

dev.off()

################################
# Plot Xi proportions for both WT and InvDxz4

dataDesc = "XiProportions_WT-and-InvDxz4"
pdf(file.path(chrXplotDir, paste(dataDesc,".pdf",sep="")), width=15, height=5)

plot(theLoci, theWTXiTPMpropFixed, 
     pch=samplePch["PatskiWT"], col=sampleCols["PatskiWT"], lwd=2,
     axes=F, main=NA, xlab=NA, ylab=NA, ylim=theYlims, xlim=theXlims)
points(theLoci[escapeBoolIdx], theWTXiTPMpropFixed[escapeBoolIdx], 
       pch=samplePchEscapee["PatskiWT"], cex=samplePchEscapeeCex["PatskiWT"], col="yellow")

points(theLoci, theInvDxz4XiTPMpropFixed, 
       pch=samplePch["PatskiInvDxz4"], col=sampleCols["PatskiInvDxz4"], lwd=2)
points(theLoci[escapeBoolIdx], theInvDxz4XiTPMpropFixed[escapeBoolIdx], 
       pch=samplePchEscapee["PatskiInvDxz4"], cex=samplePchEscapeeCex["PatskiInvDxz4"], col="yellow")

# hiDelBoolIdx = theInvDxz4XiTPMpropFixed[!escapeBoolIdx] > 0.05
# names(theInvDxz4XiTPMpropFixed[!escapeBoolIdx])[which(hiDelBoolIdx)]

segments(theLoci, theWTXiTPMpropFixed, theLoci, theInvDxz4XiTPMpropFixed, col=sampleCols["PatskiInvDxz4"], lty="solid", lwd=0.5)

axis(side=2)
mtext(side=2, line=2.5, text="Xi/(Xi+Xa)")
axis(1, at=seq(0, 180*10^6, 20*10^6), labels=seq(0, 180, 20))
title(main="Xi expression proportion", line=3)

legend("topright", bty="n",
       legend=c(paste("WT genes (", sum(!is.na(theWTXiTPMpropFixed)), ")", sep=""),
                paste("InvDxz4 genes (", sum(!is.na(theInvDxz4XiTPMpropFixed)), ")", sep=""),
                paste("Escape genes (", sum(escapeBoolIdx), ")", sep="" )),
       pch=c(samplePch["PatskiWT"], samplePch["PatskiInvDxz4"], 20),
       col=c(sampleCols["PatskiWT"], sampleCols["PatskiInvDxz4"], "yellow"))

# LOI
abline(v=firrePos, lty="dashed", lw=0.5, col="black")
abline(v=dxz4Pos, lty="dashed", lw=0.5, col="black")
abline(v=xistPos, lty="dashed", lw=0.5, col="black")
axis(3, at=c(firrePos, dxz4Pos, xistPos), labels=c("Firre", "Dxz4", "Xist"), tick = FALSE)

dev.off()


################################
# Plot Xi proportions for both newWT and InvDxz4

dataDesc = "XiProportions_newWT-and-InvDxz4"
pdf(file.path(chrXplotDir, paste(dataDesc,".pdf",sep="")), width=15, height=5)

plot(theLoci, theNewWTXiTPMpropFixed, 
     pch=samplePch["PatskiNewWT"], col=sampleCols["PatskiNewWT"], lwd=2,
     axes=F, main=NA, xlab=NA, ylab=NA, ylim=theYlims, xlim=theXlims)
points(theLoci[escapeBoolIdx], theNewWTXiTPMpropFixed[escapeBoolIdx], 
       pch=samplePchEscapee["PatskiNewWT"], cex=samplePchEscapeeCex["PatskiNewWT"], col="yellow")

points(theLoci, theInvDxz4XiTPMpropFixed, 
       pch=samplePch["PatskiInvDxz4"], col=sampleCols["PatskiInvDxz4"], lwd=2)
points(theLoci[escapeBoolIdx], theInvDxz4XiTPMpropFixed[escapeBoolIdx], 
       pch=samplePchEscapee["PatskiInvDxz4"], cex=samplePchEscapeeCex["PatskiInvDxz4"], col="yellow")

# hiDelBoolIdx = theInvDxz4XiTPMpropFixed[!escapeBoolIdx] > 0.05
# names(theInvDxz4XiTPMpropFixed[!escapeBoolIdx])[which(hiDelBoolIdx)]

segments(theLoci, theNewWTXiTPMpropFixed, theLoci, theInvDxz4XiTPMpropFixed, col=sampleCols["PatskiInvDxz4"], lty="solid", lwd=0.5)

axis(side=2)
mtext(side=2, line=2.5, text="Xi/(Xi+Xa)")
axis(1, at=seq(0, 180*10^6, 20*10^6), labels=seq(0, 180, 20))
title(main="Xi expression proportion", line=3)

legend("topright", bty="n",
       legend=c(paste("NewWT genes (", sum(!is.na(theNewWTXiTPMpropFixed)), ")", sep=""),
                paste("InvDxz4 genes (", sum(!is.na(theInvDxz4XiTPMpropFixed)), ")", sep=""),
                paste("Escape genes (", sum(escapeBoolIdx), ")", sep="" )),
       pch=c(samplePch["PatskiNewWT"], samplePch["PatskiInvDxz4"], 20),
       col=c(sampleCols["PatskiNewWT"], sampleCols["PatskiInvDxz4"], "yellow"))

# LOI
abline(v=firrePos, lty="dashed", lw=0.5, col="black")
abline(v=dxz4Pos, lty="dashed", lw=0.5, col="black")
abline(v=xistPos, lty="dashed", lw=0.5, col="black")
axis(3, at=c(firrePos, dxz4Pos, xistPos), labels=c("Firre", "Dxz4", "Xist"), tick = FALSE)

dev.off()


#################################################################
# Plot log2 FCs (from DESeq2 results)

chrXplotDir = file.path(mainChrXplotDir, "Log2FC")
if ( ! file.exists( chrXplotDir ) ) {
  if ( ! dir.create( chrXplotDir, showWarnings=FALSE, recursive=TRUE )  ) chrXplotDir = mainChrXplotDir
}

# Plot the figures
theYlims=c(-3,3)
theXlims=c(0, 180*10^6)

firrePos = (50563120+(abs(50563120-50635321)/2))
dxz4Pos = (75764754+(abs(75637518-75764754)/2))
xistPos = (103460373+(abs(103460373-103483233)/2))

for (comp in names(DEresultFiles)) {
  for (allele in names(DEresultFiles[[comp]])) {
    
    dataDesc = paste("log2FC", comp, allele, sep="_")
    pdf(file.path(chrXplotDir, paste(dataDesc,".pdf",sep="")), width=15, height=5)
    
    currTable = DEtables.chrX[[comp]][[allele]]
    names(currTable) = c("geneId", "basemean", "log2FoldChange", 
                         "lfcSE", "stat", "pvalue", "padj")
    geneIdXref = match(currTable$geneId, expressionTables.chrX.average[["refTPM"]]$geneId)
    
    theLoci = (expressionTables.chrX.average[["refTPM"]]$start + (expressionTables.chrX.average[["refTPM"]]$stop - expressionTables.chrX.average[["refTPM"]]$start)/2)[geneIdXref]
    theLogRatio = -1 * currTable$log2FoldChange
    thePvalue = currTable$pvalue
    thePadj = currTable$padj
    DEboolIdx = thePadj<=0.05 & abs(theLogRatio) >= 0.5
    
    sum(!is.na(thePvalue))
    plot(theLoci, theLogRatio, type="h", col="gray", lwd=4, 
         axes=F, main=NA, xlab=NA, ylab=NA, ylim=theYlims, xlim=theXlims)
    lines(theLoci[DEboolIdx], theLogRatio[DEboolIdx], type="h", col="red", lwd=4)
         
    axis(side=2)
    mtext(side=2, line=2.5, text=paste(allele, " log2FC", sep=""))
    axis(1, at=seq(0, 180*10^6, 20*10^6), labels=seq(0, 180, 20))
    title(main=comp, line=3)
    
    legend("bottomleft", bty="n",
           legend=c(paste("Genes with FC values (", sum(!is.na(thePvalue)), ")", sep=""),
                    paste("DE genes (", sum(DEboolIdx, na.rm=TRUE), ")", sep="")),
           lty=c(1,1), col=c("gray", "red"), lwd=c(2,2))
    
    # Log2FC thresholds
    abline(h=0.5, lty="dashed", lw=0.5, col="blue")
    abline(h=-0.5, lty="dashed", lw=0.5, col="blue")
    
    # LOI
    abline(v=firrePos, lty="dashed", lw=0.5, col="black")
    abline(v=dxz4Pos, lty="dashed", lw=0.5, col="black")
    abline(v=xistPos, lty="dashed", lw=0.5, col="black")
    axis(3, at=c(firrePos, dxz4Pos, xistPos), labels=c("Firre", "Dxz4", "Xist"), tick = FALSE)
    
    dev.off()
  }
}


#################################################################
# Plot log2 FCs (from DESeq2 results)
# Escape vs not escape

chrXplotDir = file.path(mainChrXplotDir, "Log2FC.escapees")
if ( ! file.exists( chrXplotDir ) ) {
  if ( ! dir.create( chrXplotDir, showWarnings=FALSE, recursive=TRUE )  ) chrXplotDir = mainChrXplotDir
}

# Plot the figures
theYlims=c(-3,3)
theXlims=c(0, 180*10^6)

firrePos = (50563120+(abs(50563120-50635321)/2))
dxz4Pos = (75764754+(abs(75637518-75764754)/2))
xistPos = (103460373+(abs(103460373-103483233)/2))

lwdA = 6
lwdB = 2

for (comp in names(DEresultFiles)) {
  for (allele in names(DEresultFiles[[comp]])) {
    
    dataDesc = paste("log2FC", comp, allele, sep="_")
    pdf(file.path(chrXplotDir, paste(dataDesc,".pdf",sep="")), width=15, height=5)
    
    currTable = DEtables.chrX[[comp]][[allele]]
    names(currTable) = c("geneId", "basemean", "log2FoldChange", 
                         "lfcSE", "stat", "pvalue", "padj")
    geneIdXref = match(currTable$geneId, expressionTables.chrX.average[["refTPM"]]$geneId)
    
    theLoci = (expressionTables.chrX.average[["refTPM"]]$start + (expressionTables.chrX.average[["refTPM"]]$stop - expressionTables.chrX.average[["refTPM"]]$start)/2)[geneIdXref]
    theLogRatio = -1 * currTable$log2FoldChange
    thePvalue = currTable$pvalue
    thePadj = currTable$padj
 
    DEboolIdx = thePadj<=0.05 & abs(theLogRatio) >= 0.5
    escapeBoolIdx = currTable$geneId %in% newEscapeGenes
    
    sum(!is.na(thePvalue))
    plot(theLoci, theLogRatio, type="h", col="gray", lwd=lwdA, 
         axes=F, main=NA, xlab=NA, ylab=NA, ylim=theYlims, xlim=theXlims)
    lines(theLoci[escapeBoolIdx], theLogRatio[escapeBoolIdx], type="h", col="green", lwd=lwdA)
    lines(theLoci[DEboolIdx], theLogRatio[DEboolIdx], type="h", col="red", lwd=lwdB)
    
    axis(side=2)
    mtext(side=2, line=2.5, text=paste(allele, " log2FC", sep=""))
    axis(1, at=seq(0, 180*10^6, 20*10^6), labels=seq(0, 180, 20))
    title(main=comp, line=3)
    
    legend("bottomleft", bty="n",
           legend=c(paste("Genes with FC values (", sum(!is.na(thePvalue)), ")", sep=""),
                    paste("Escape genes (", sum(escapeBoolIdx, na.rm=TRUE), ")", sep=""),
                    paste("DE genes (", sum(DEboolIdx, na.rm=TRUE), ")", sep="")),
           lty=c(1,1,1), col=c("gray", "green", "red"), lwd=c(lwdA,lwdA,lwdB))
    
    # Log2FC thresholds
    abline(h=0.5, lty="dashed", lw=0.5, col="blue")
    abline(h=-0.5, lty="dashed", lw=0.5, col="blue")
    
    # LOI
    abline(v=firrePos, lty="dashed", lw=0.5, col="black")
    abline(v=dxz4Pos, lty="dashed", lw=0.5, col="black")
    abline(v=xistPos, lty="dashed", lw=0.5, col="black")
    axis(3, at=c(firrePos, dxz4Pos, xistPos), labels=c("Firre", "Dxz4", "Xist"), tick = FALSE)
    
    dev.off()
  }
}










###########################################################
###########################################################
# 
# 171205
# Get TPM values for DESeq2 comparisons
#
###########################################################

dataOutputDir = file.path(outputDir, "tabulatedData4rebuttal")
if ( ! file.exists( dataOutputDir ) ) {
  if ( ! dir.create( dataOutputDir, showWarnings=FALSE, recursive=TRUE )  ) dataOutputDir = outputDir
}


################################
# Data

###############
# Load DE result tables
altOutputDir = "~/proj/2017_RNAseq_analysis/results/gbonora/20170220_RNAseq_reprocessAllLibs/expressionAnalysis.correctedAltToRefRatio"

DEresultDirs = list()
DEresultDirs[["WT-vs-Del"]] = file.path(altOutputDir, "DESeq2/WT-vs-Del1a")
DEresultDirs[["WT-vs-InvDxz4"]] = file.path(outputDir, "DESeq2/WT-vs-InvDxz4")
DEresultDirs[["newWT-vs-InvDxz4"]] = file.path(outputDir, "DESeq2/newWT-vs-InvDxz4")
DEresultDirs[["WT-vs-WT_5aza"]] = file.path(altOutputDir, "DESeq2/5-aza")
DEresultDirs[["Del1-vs-Del1_5aza"]] = file.path(altOutputDir, "DESeq2/5-aza")

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

DEresultFiles[["WT-vs-WT_5aza"]] = list ()
DEresultFiles[["WT-vs-WT_5aza"]][["Xi"]] = file.path(DEresultDirs[["WT-vs-WT_5aza"]], 
                                                     "Patski.AllData_DESeq2genes_logP=1.3_logFC=0.5_B6.WT-vs-WT_5aza.chrX.tsv")
DEresultFiles[["WT-vs-WT_5aza"]][["Xa"]] = file.path(DEresultDirs[["WT-vs-WT_5aza"]], 
                                                     "Patski.AllData_DESeq2genes_logP=1.3_logFC=0.5_sp.WT-vs-WT_5aza.chrX.tsv")

DEresultFiles[["Del1-vs-Del1_5aza"]] = list ()
DEresultFiles[["Del1-vs-Del1_5aza"]][["Xi"]] = file.path(DEresultDirs[["Del1-vs-Del1_5aza"]], 
                                                     "Patski.AllData_DESeq2genes_logP=1.3_logFC=0.5_B6.Del1-vs-Del1_5aza.chrX.tsv")
DEresultFiles[["Del1-vs-Del1_5aza"]][["Xa"]] = file.path(DEresultDirs[["Del1-vs-Del1_5aza"]], 
                                                     "Patski.AllData_DESeq2genes_logP=1.3_logFC=0.5_sp.Del1-vs-Del1_5aza.chrX.tsv")

DEtables.chrX = list()
for (comp in names(DEresultFiles)) {
  DEtables.chrX[[comp]] = list()
  for (allele in names(DEresultFiles[[comp]])) {
    tableName = DEresultFiles[[comp]][[allele]]
    print(tableName)
    DEtables.chrX[[comp]][[allele]] = read.table(tableName, sep = "\t", header=TRUE)
  }
}


###############
# Load expression tables
tableDir = file.path(outputDir, "expressionTables")
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


###############
# Get WT-vs-Del TPM data for Table S8

fieldsOI = c("PatskiWT", "PatskiWT_0uMaza_rep1", "PatskiWT_0uMaza_rep2", 
             "PatskiDel1", "PatskiDel1_0uMaza_rep1", "PatskiDel1_0uMaza_rep2" )

head(DEtables.chrX[["WT-vs-Del"]][["Xi"]])

pvals = DEtables.chrX[["WT-vs-Del"]][["Xi"]][,7]
absLog2FC = abs(DEtables.chrX[["WT-vs-Del"]][["Xi"]][,3])
baseMean = abs(DEtables.chrX[["WT-vs-Del"]][["Xi"]][,2])

orderedGeneIds = DEtables.chrX[["WT-vs-Del"]][["Xi"]][ order(pvals, -absLog2FC, -baseMean), ][,1]

orderedGeneIdsIdx = match(orderedGeneIds, expressionTables.chrX[["refTPM"]]$geneId)
expressionTables.chrX[["refTPM"]][orderedGeneIdsIdx ,c("geneId")]

TPMdata2export = cbind(DEtables.chrX[["WT-vs-Del"]][["Xi"]][orderedGeneIds,],  # This sorting works b/c gene names are factors!
                       # Reverse comparison
                       -1* DEtables.chrX[["WT-vs-Del"]][["Xi"]][orderedGeneIds,3],
                       -1* DEtables.chrX[["WT-vs-Del"]][["Xi"]][orderedGeneIds,5],                       
                       expressionTables.chrX[["refTPM"]][orderedGeneIdsIdx ,c("geneId", fieldsOI)],
                       expressionTables.chrX[["altTPM"]][orderedGeneIdsIdx ,c("geneId", fieldsOI)])
nrow(TPMdata2export)

exportFilename = file.path(dataOutputDir, "WT-vs-Del_DEgenes_refAltTPM.tsv")

write.table(TPMdata2export, exportFilename, sep="\t", row.names=FALSE, col.names=TRUE)

# Slc16a2
Slc16a2Idx = match("Slc16a2", expressionTables.chrX[["refTPM"]]$geneId)
expressionTables.chrX[["refTPM"]][Slc16a2Idx, c("geneId",fieldsOI)]
expressionTables.chrX[["refcounts"]][Slc16a2Idx, c("geneId",fieldsOI)]
# geneId PatskiWT PatskiWT_0uMaza_rep1
# 872 Slc16a2 122.4019             32.84298
# PatskiWT_0uMaza_rep2 PatskiDel1 PatskiDel1b
# 872             34.77983        NaN    665.6986
# PatskiDel1_0uMaza_rep1
# 872                      0

# Escape genes
EscapeGeneIdsIdx = match(newEscapeGenes, expressionTables.chrX[["refTPM"]]$geneId)
# expressionTables.chrX[["refcounts"]][EscapeGeneIdsIdx, c("geneId",fieldsOI)]
TPMdata2export = cbind(expressionTables.chrX[["counts"]][EscapeGeneIdsIdx ,c("geneId", fieldsOI)],
                       expressionTables.chrX[["refcounts"]][EscapeGeneIdsIdx ,c("geneId", fieldsOI)],
                       expressionTables.chrX[["altcounts"]][EscapeGeneIdsIdx ,c("geneId", fieldsOI)])
exportFilename = file.path(dataOutputDir, "Escape_WT-and-Del_refAltCounts.tsv")
write.table(TPMdata2export, exportFilename, sep="\t", row.names=FALSE, col.names=TRUE)

# Genes with no allelic data
naIndex = apply(is.na(expressionTables.chrX[["refcounts"]][ ,fieldsOI]), 1, any)
sum(naIndex)
TPMdata2export = cbind(expressionTables.chrX[["counts"]][naIndex ,c("geneId", fieldsOI)],
                       expressionTables.chrX[["refcounts"]][naIndex ,c("geneId", fieldsOI)],
                       expressionTables.chrX[["altcounts"]][naIndex ,c("geneId", fieldsOI)])
exportFilename = file.path(dataOutputDir, "NoAllelicData_WT-and-Del_refAltCounts.tsv")
write.table(TPMdata2export, exportFilename, sep="\t", row.names=FALSE, col.names=TRUE)


###############
# Get Inv-Dxz4 vs Patski2-4 WT data for Table S10

fieldsOI = c("PatskiWTa", "PatskiWTb", "PatskiInvDxz4a", "PatskiInvDxz4b" )

head(DEtables.chrX[["newWT-vs-InvDxz4"]][["Xi"]])

pvals = DEtables.chrX[["newWT-vs-InvDxz4"]][["Xi"]][,7]
absLog2FC = abs(DEtables.chrX[["newWT-vs-InvDxz4"]][["Xi"]][,3])
baseMean = abs(DEtables.chrX[["newWT-vs-InvDxz4"]][["Xi"]][,2])

orderedGeneIds = DEtables.chrX[["newWT-vs-InvDxz4"]][["Xi"]][ order(pvals, -absLog2FC, -baseMean), ][,1]

orderedGeneIdsIdx = match(orderedGeneIds, expressionTables.chrX[["refTPM"]]$geneId)
expressionTables.chrX[["refTPM"]][orderedGeneIdsIdx ,c("geneId")]

TPMdata2export = cbind(DEtables.chrX[["newWT-vs-InvDxz4"]][["Xi"]][orderedGeneIds,],  # This sorting works b/c gene names are factors!
                       # Reverse comparison
                       -1* DEtables.chrX[["newWT-vs-InvDxz4"]][["Xi"]][orderedGeneIds,3],
                       -1* DEtables.chrX[["newWT-vs-InvDxz4"]][["Xi"]][orderedGeneIds,5],                       
                       expressionTables.chrX[["refTPM"]][orderedGeneIdsIdx ,c("geneId", fieldsOI)],
                       expressionTables.chrX[["altTPM"]][orderedGeneIdsIdx ,c("geneId", fieldsOI)])

exportFilename = file.path(dataOutputDir, "newWT-vs-InvDxz4_DEgenes_refAltTPM.tsv")

write.table(TPMdata2export, exportFilename, sep="\t", row.names=FALSE, col.names=TRUE)






###############
# Get  5-aza-2dC treated and untreated WT TPM data for Table S11

fieldsOI = c("PatskiWT_0uMaza_rep1", "PatskiWT_0uMaza_rep2", 
             "PatskiWT_4uMaza_rep1", "PatskiWT_4uMaza_rep2"  )

head(DEtables.chrX[["WT-vs-WT_5aza"]][["Xi"]])

pvals = DEtables.chrX[["WT-vs-WT_5aza"]][["Xi"]][,7]
absLog2FC = abs(DEtables.chrX[["WT-vs-WT_5aza"]][["Xi"]][,3])
baseMean = abs(DEtables.chrX[["WT-vs-WT_5aza"]][["Xi"]][,2])

orderedGeneIds = DEtables.chrX[["WT-vs-WT_5aza"]][["Xi"]][ order(pvals, -absLog2FC, -baseMean), ][,1]

orderedGeneIdsIdx = match(orderedGeneIds, expressionTables.chrX[["refTPM"]]$geneId)
expressionTables.chrX[["refTPM"]][orderedGeneIdsIdx ,c("geneId")]

TPMdata2export = cbind(DEtables.chrX[["WT-vs-WT_5aza"]][["Xi"]][orderedGeneIds,],  # This sorting works b/c gene names are factors!
                       # Reverse comparison
                       -1* DEtables.chrX[["WT-vs-WT_5aza"]][["Xi"]][orderedGeneIds,3],
                       -1* DEtables.chrX[["WT-vs-WT_5aza"]][["Xi"]][orderedGeneIds,5],                       
                       expressionTables.chrX[["refTPM"]][orderedGeneIdsIdx ,c("geneId", fieldsOI)],
                       expressionTables.chrX[["altTPM"]][orderedGeneIdsIdx ,c("geneId", fieldsOI)])

exportFilename = file.path(dataOutputDir, "WT-vs-WT_5aza_DEgenes_refAltTPM.tsv")

write.table(TPMdata2export, exportFilename, sep="\t", row.names=FALSE, col.names=TRUE)






###############
# Get 5-aza-2dC treated and untreated Del1 TPM data for Table S12

fieldsOI = c("PatskiDel1_0uMaza_rep1", "PatskiDel1_0uMaza_rep2", 
             "PatskiDel1_4uMaza_rep1", "PatskiDel1_4uMaza_rep2"  )

head(DEtables.chrX[["Del1-vs-Del1_5aza"]][["Xi"]])

# For some reason the Del1-vs-Del1_5aza tables contain genome-wide data
temp = DEtables.chrX[["Del1-vs-Del1_5aza"]][["Xi"]]
temp = temp[,c(1, grep("Del1.vs.Del1_5aza", names(temp)))]
# Use WT-vs-WT_5aza DE genes to get x-linked genes
XlinkedGenesIdx = match(DEtables.chrX[["WT-vs-WT_5aza"]][["Xi"]][,1], temp[,1])
temp = temp[XlinkedGenesIdx,]
DEtables.chrX[["Del1-vs-Del1_5aza"]][["Xi"]] = temp

pvals = DEtables.chrX[["Del1-vs-Del1_5aza"]][["Xi"]][,7]
absLog2FC = abs(DEtables.chrX[["Del1-vs-Del1_5aza"]][["Xi"]][,3])
baseMean = abs(DEtables.chrX[["Del1-vs-Del1_5aza"]][["Xi"]][,2])

orderedGeneIds = DEtables.chrX[["Del1-vs-Del1_5aza"]][["Xi"]][ order(pvals, -absLog2FC, -baseMean), ][,1]

orderedGeneIdsIdx = match(orderedGeneIds, expressionTables.chrX[["refTPM"]]$geneId)
expressionTables.chrX[["refTPM"]][orderedGeneIdsIdx ,c("geneId")]

TPMdata2export = cbind(DEtables.chrX[["Del1-vs-Del1_5aza"]][["Xi"]][orderedGeneIds,],  # This sorting works b/c gene names are factors!
                       # Reverse comparison
                       -1* DEtables.chrX[["Del1-vs-Del1_5aza"]][["Xi"]][orderedGeneIds,3],
                       -1* DEtables.chrX[["Del1-vs-Del1_5aza"]][["Xi"]][orderedGeneIds,5],                       
                       expressionTables.chrX[["refTPM"]][orderedGeneIdsIdx ,c("geneId", fieldsOI)],
                       expressionTables.chrX[["altTPM"]][orderedGeneIdsIdx ,c("geneId", fieldsOI)])

exportFilename = file.path(dataOutputDir, "Del1-vs-Del1_5aza_DEgenes_refAltTPM.tsv")

write.table(TPMdata2export, exportFilename, sep="\t", row.names=FALSE, col.names=TRUE)









###########################################################
###########################################################
# 
# 180115
# Get loci of autosomal diploid DE genes
#
###########################################################

dataOutputDir = file.path(outputDir, "tabulatedData4rebuttal")
if ( ! file.exists( dataOutputDir ) ) {
  if ( ! dir.create( dataOutputDir, showWarnings=FALSE, recursive=TRUE )  ) dataOutputDir = outputDir
}

DEadjPthresh = 0.05; DElogFCthresh = 0.5

################################
# Data

###############
# Load DE result tables
altOutputDir = "~/proj/2017_RNAseq_analysis/results/gbonora/20170220_RNAseq_reprocessAllLibs/expressionAnalysis.correctedAltToRefRatio"

DEresultDirs = list()
DEresultDirs[["WT-vs-Del"]] = file.path(altOutputDir, "DESeq2/WT-vs-Del1a")

DEresultFiles = list()
DEresultFiles[["WT-vs-Del"]] = list ()
DEresultFiles[["WT-vs-Del"]][["autosomalDiploid"]] = file.path(DEresultDirs[["WT-vs-Del"]], 
                                                  "Patski.AllData_DESeq2genes_logP=1.3_logFC=0.5_diploid.WT-vs-Del1.autosomes.tsv")
DEresultFiles[["WT-vs-Del"]][["chrXdiploid"]] = file.path(DEresultDirs[["WT-vs-Del"]], 
                                                               "Patski.AllData_DESeq2genes_logP=1.3_logFC=0.5_diploid.WT-vs-Del1.chrX.tsv")

DEtables = list()
for (comp in names(DEresultFiles)) {
  DEtables[[comp]] = list()
  for (allele in names(DEresultFiles[[comp]])) {
    tableName = DEresultFiles[[comp]][[allele]]
    print(tableName)
    DEtables[[comp]][[allele]] = read.table(tableName, sep = "\t", header=TRUE)
  }
}

# Must convert column to string
str(DEtables[["WT-vs-Del"]][["autosomalDiploid"]])
DEtables[["WT-vs-Del"]][["autosomalDiploid"]][,1] = as.character(DEtables[["WT-vs-Del"]][["autosomalDiploid"]][,1])
str(DEtables[["WT-vs-Del"]][["autosomalDiploid"]])

DEtables[["WT-vs-Del"]][["chrXdiploid"]][,1] = as.character(DEtables[["WT-vs-Del"]][["chrXdiploid"]][,1])


###############
# Load expression tables
tableDir = file.path(outputDir, "expressionTables")
expressionTables = list()
for (expDataType in c("Counts", "TPM", "RPKM")) {
  #for (genome in c("", "ref", "alt")) {
  for (genome in c("")) {
    if ( genome=="" && expDataType == "Counts") {
      expDataType = "counts"
    }
    tableName = paste(genome, expDataType, sep = "")
    print(tableName)
    #bigOldFilename = file.path(tableDir, paste(tableName, ".chrX.tsv.gz", sep = ""))
    bigOldFilename = file.path(tableDir, paste(tableName, ".tsv.gz", sep = ""))
    expressionTables[[tableName]] = read.table(gzfile(bigOldFilename), sep = "\t", header=TRUE)
  }
}

###############
# Get WT-vs-Del position and TPM data for Table S9

fieldsOI = c("PatskiWT", "PatskiWT_0uMaza_rep1", "PatskiWT_0uMaza_rep2", 
             "PatskiDel1", "PatskiDel1_0uMaza_rep1", "PatskiDel1_0uMaza_rep2" )

head(DEtables[["WT-vs-Del"]][["autosomalDiploid"]])

pvals = DEtables[["WT-vs-Del"]][["autosomalDiploid"]][,7]
absLog2FC = abs(DEtables[["WT-vs-Del"]][["autosomalDiploid"]][,3])
baseMean = abs(DEtables[["WT-vs-Del"]][["autosomalDiploid"]][,2])

# Select sig DE genes
DEgeneIdx = which(pvals<=DEadjPthresh & absLog2FC>=DElogFCthresh)
length(DEgeneIdx)
# 2689
DEtables[["WT-vs-Del"]][["autosomalDiploid"]] = DEtables[["WT-vs-Del"]][["autosomalDiploid"]][DEgeneIdx,]

pvals = DEtables[["WT-vs-Del"]][["autosomalDiploid"]][,7]
absLog2FC = abs(DEtables[["WT-vs-Del"]][["autosomalDiploid"]][,3])
baseMean = abs(DEtables[["WT-vs-Del"]][["autosomalDiploid"]][,2])

orderedGeneIds.auto = DEtables[["WT-vs-Del"]][["autosomalDiploid"]][ order(pvals, -absLog2FC, -baseMean), ][,1]
head(orderedGeneIds.auto)
orderedGeneIdsDEtableIdx = match(orderedGeneIds.auto, DEtables[["WT-vs-Del"]][["autosomalDiploid"]][,1])
head(DEtables[["WT-vs-Del"]][["autosomalDiploid"]][orderedGeneIdsDEtableIdx,1])

orderedGeneIdsIdx = match(orderedGeneIds, expressionTables[["TPM"]]$geneId)
head(expressionTables[["TPM"]][orderedGeneIdsIdx ,c("geneId")])

TPMdata2export.autosomes = cbind(DEtables[["WT-vs-Del"]][["autosomalDiploid"]][orderedGeneIdsDEtableIdx,],
                       # Reverse comparison
                       -1* DEtables[["WT-vs-Del"]][["autosomalDiploid"]][orderedGeneIdsDEtableIdx,3],
                       -1* DEtables[["WT-vs-Del"]][["autosomalDiploid"]][orderedGeneIdsDEtableIdx,5],                       
                       expressionTables[["TPM"]][orderedGeneIdsIdx ,c("geneId", "chr",	"start", "stop", 
                                                                        fieldsOI)])
nrow(TPMdata2export.autosomes)

# Reorder by chrom and position
chr = TPMdata2export.autosomes[,"chr"]
chr = as.numeric(sub("chr", "", chr))
pos = TPMdata2export.autosomes[,"start"]
TPMdata2export.autosomes = TPMdata2export.autosomes[order(chr, pos), ]

exportFilename = file.path(dataOutputDir, "WT-vs-Del_DEgenes_autosomalDiploidTPM.tsv")
write.table(TPMdata2export.autosomes, exportFilename, sep="\t", row.names=FALSE, col.names=TRUE)







###############
# Tabulate DE gene counts

###
# Add chrX data

# Doesn't need to be ordered:
chrXIdx = match(DEtables[["WT-vs-Del"]][["chrXdiploid"]][,1], expressionTables[["TPM"]]$geneId)

# Must run above code first to populate DEtables.chrX and expressionTables.chrX
TPMdata2export.chrX = cbind(DEtables[["WT-vs-Del"]][["chrXdiploid"]],  
                       # Reverse comparison
                       -1* DEtables[["WT-vs-Del"]][["chrXdiploid"]][,3],
                       -1* DEtables[["WT-vs-Del"]][["chrXdiploid"]][,5],                       
                       expressionTables.chrX[["refTPM"]][chrXIdx ,c("geneId", fieldsOI)],
                       expressionTables.chrX[["altTPM"]][chrXIdx ,c("geneId", fieldsOI)])
nrow(TPMdata2export.chrX)

pvals.chrX = DEtables[["WT-vs-Del"]][["chrXdiploid"]][,7]
absLog2FC.chrX = abs(DEtables[["WT-vs-Del"]][["chrXdiploid"]][,3])
baseMean.chrX = abs(DEtables[["WT-vs-Del"]][["chrXdiploid"]][,2])

# Select sig DE genes
DEgeneIdx.chrX = which(pvals.chrX<=DEadjPthresh & absLog2FC.chrX>=DElogFCthresh)
length(DEgeneIdx.chrX)


###
# tabulate

TPMdata2export.ALL = c(as.character(TPMdata2export.autosomes[,"chr"]),
                           rep("chrX", length(DEgeneIdx.chrX)))

chromOrder = c(paste("chr", seq(1:19), sep=""), "chrX")

tabulatedDEgenes = table(TPMdata2export.ALL)
tabulatedDEgenes[chromOrder]

# Get gene density per chromosome to normalize counts
tabulatedGenesPerChrom = table(expressionTables[[1]][,"chr"])
tabulatedGenesPerChrom[chromOrder]

# output
outputTable = cbind(tabulatedDEgenes[chromOrder], 
                    tabulatedGenesPerChrom[chromOrder],
                    tabulatedDEgenes[chromOrder]/tabulatedGenesPerChrom[chromOrder])
colnames(outputTable) = c("DEgenes", "totalGenes", "DEproportion")

exportFilename = file.path(dataOutputDir, "WT-vs-Del_DEgenes_autosomalDiploidTPM_tabulatedByChromosome.tsv")
write.table(outputTable, exportFilename, sep="\t", row.names=TRUE, col.names=NA)
