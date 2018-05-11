######################################################################################################################
######################################################################################################################
#
# CTCF Diploid peak Analysis
# Author: GB
# 20170331
#
######################################################################################################################
######################################################################################################################

sessionInfo()

require(RColorBrewer)
require(gplots)

require(ggplot2)
# Fonts: http://www.cookbook-r.com/Graphs/Fonts/
require(gridExtra)
require(grid)

# http://seananderson.ca/2013/10/19/reshape.html
require(reshape2)

#require(genefilter)
#require(calibrate)

require(TeachingDemos)
require(GGally)
require(colorspace)

# # install.packages("xlsx")
# require(xlsx)
require(openxlsx)

######################################################################################################################
######################################################################################################################
# Important functions for later

# version of vioplot with colour
sourceDir = "<scriptDirectory>"
source(file.path(sourceDir, "colVioplot.R"))

# The normalize.quantiles function in the Bolstad 'preprocessCore' R package 
# does not appear to retain the order of 'probes' of equal value.
# So wrote a custom version that does.
source(file.path(sourceDir, "quantNormGB.R"))

source(file.path(sourceDir, "corChart.R"))

# source(file.path(sourceDir, "multiplot.R"))

# For loggin RPKMs in a list
addThenLog = function(x, add=1, base=2) {
  return(log(x+add, base))
}

# read key
readkey <- function() {
  cat ("Press [enter] to continue")
  line <- readline()
}


# https://github.com/tidyverse/ggplot2/wiki/share-a-legend-between-two-ggplot2-graphs
grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right"), title = "") {
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight),
                                            top = title),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth),
                                           top = title))
  
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
}

# Just extract legend
# https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
# http://stackoverflow.com/questions/13649473/add-a-common-legend-for-combined-ggplots
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}





######################################################################################################################
######################################################################################################################
# Preliminaries

projectDir = "<projectDirectory>"
subDir = "20171005_ChIPseq_CTCF"
workDir = "patskiWTvsDelvsInv.allelicCTCF"

outputDir = file.path(projectDir, subDir, "diploidPeakSegregation")
if ( ! file.exists( outputDir ) ) {
  if ( ! dir.create( outputDir, showWarnings=FALSE, recursive=TRUE )  ) outputDir = getwd()
}
setwd(outputDir)
getwd()

#save.image() # To .RData
#load(file = ".RData")

refdataDir="<refdataDirectory>"

dataDir="<dataDirectory>"

dataSets=c("WT_CTCF", "2-4clone_CTCF", "Del1_CTCF", "InvDxz4_CTCF")

peakLevelASFolder= file.path(dataDir, "snpPileups")

sampleCols=c("blue", "lightBlue", "red", gray(0.5))
names(sampleCols) =c("WT", "2-4clone", "Del", "InvDxz4")

# http://stackoverflow.com/questions/30219738/is-there-a-way-to-programmatically-darken-the-color-given-rgb-values/30412528
# Lighter versions of colors
sampleColsLight <- readhex(file = textConnection(paste(rgb(t(col2rgb(sampleCols)), maxColorValue=255), collapse = "\n")),
                 class = "RGB")

#additive decrease of lightness
#transform to hue/lightness/saturation colorspace
sampleColsLightS25 = as(sampleColsLight, "HLS")
sampleColsLightS25@coords[, "L"] <- pmax(0, sampleColsLightS25@coords[, "L"] - 0.25)
#going via rgb seems to work better  
sampleColsLightS25 <- hex(as(sampleColsLightS25, "RGB"))

#multiplicative decrease of lightness
#transform to hue/lightness/saturation colorspace
sampleColsLightD2 <- as(sampleColsLight, "HLS")
sampleColsLightD4 <- as(sampleColsLight, "HLS")
sampleColsLightD8 <- as(sampleColsLight, "HLS")

sampleColsLightD2@coords[, "L"] <- sampleColsLightD2@coords[, "L"] * 0.50
sampleColsLightD4@coords[, "L"] <- sampleColsLightD4@coords[, "L"] * 0.25
sampleColsLightD8@coords[, "L"] <- sampleColsLightD8@coords[, "L"] * 0.125

#going via rgb seems to work better  
sampleColsLightD2 <- hex(as(sampleColsLightD2, "RGB"))
sampleColsLightD4 <- hex(as(sampleColsLightD4, "RGB"))
sampleColsLightD8 <- hex(as(sampleColsLightD8, "RGB"))

names(sampleColsLightD2) = names(sampleColsLightD4) =  names(sampleColsLightD8) = names(sampleCols)
names(sampleColsLightS25) = names(sampleCols)

# barplot(rep(10, 4), col=rbind(sampleCols, sampleColsLightD2, sampleColsLightD4, sampleColsLightD8) [,3])







######################################################################################################################
######################################################################################################################
# Allele-specific coverage data

# SNP level data not necessary here
# SNPlevelASEfile = file.path(peakLevelASFolder, "Patski.SNPpileup.annotated.tsv.gz")
# SNPlevelASdata = read.table(gzfile(SNPlevelASEfile), header=TRUE, sep="\t")
# head(SNPlevelASdata[,1:10])
# dim(SNPlevelASdata)

# peak level data
peakLevelASdata = list()
for (dataSet in dataSets) {
  peakLevelASfile = file.path(peakLevelASFolder, paste(dataSet, ".SNPcoverage.sumOverPeaks.tsv.gz", sep=""))
  peakLevelASdata[[dataSet]] = read.table(gzfile(peakLevelASfile), header=TRUE, sep="\t") #, fill = TRUE)
}

# Rename columns, removing "ChIPseq_Patski_"
for (dataSet in dataSets) {
  print( sub("ChIPseq_Patski_", "", colnames(peakLevelASdata[[dataSet]])) )
  colnames(peakLevelASdata[[dataSet]]) = sub("ChIPseq_Patski_", "", sub("2.4", "2-4", colnames(peakLevelASdata[[dataSet]])))
}

for (dataSet in dataSets) {
  # print(colnames(peakLevelASdata[[dataSet]]))
  head(peakLevelASdata[[dataSet]])
  cat (dataSet, "\t", paste(dim(peakLevelASdata[[dataSet]]), collapse="\t"), "\n", sep="")
}

# Chr indices
chrIndices = list()
for (dataSet in dataSets) {
  chrAidx = which(as.character(peakLevelASdata[[dataSet]]$chr) %in% paste("chr", 1:22, sep=""))
  chrDidx = which(as.character(peakLevelASdata[[dataSet]]$chr) %in% paste("chr", c(1,2,5:22), sep=""))
  chrXidx = which(as.character(peakLevelASdata[[dataSet]]$chr) == "chrX")
  chrIndices[[dataSet]] = list(chrAidx, chrDidx, chrXidx)
  names(chrIndices[[dataSet]]) = c("autosomes", "diplosomes", "chrX")
  for (currChr in paste("chr", seq(1,19,1), sep="")) {
    chrIndices[[dataSet]][[currChr]] = which(as.character(peakLevelASdata[[dataSet]]$chr) %in% currChr)
  }
}
names(chrIndices[[1]])

# peakLevelASdata[[dataSet]][chrIndices[[dataSet]][["chrX"]],][1:5,1:3]












######################################################################################################################
######################################################################################################################
# SNP counts per peak (coverage)
# Total SNP coverage per cell line
# allelic SNP proportion
# allelic prop correction factor

peakLevelAScovAndRatios = list()
obsAutosomeAltRefRatio = list()
obsXAltRefRatio = list()

for (dataSet in dataSets) {
  SNPcounts = peakLevelASdata[[dataSet]][,7] # Same for all peaks: peakLevelASdata[[dataSet]][1:10,c(7,11)]
  altCol = paste(dataSet, ".altCov.sumOverPeaks", sep="")
  refCol= paste(dataSet, ".refCov.sumOverPeaks", sep="")
  altCov = as.numeric(as.character(peakLevelASdata[[dataSet]][, altCol]))
  refCov = as.numeric(as.character(peakLevelASdata[[dataSet]][, refCol]))
  totCov = altCov + refCov
  altProp = altCov / totCov
  refProp = 1-altProp

  ######
  # Adjust Reference counts
  # Expect this to be 1:1
  obsAutosomeAltRefRatio[[dataSet]] = sum(altCov[ chrIndices[[dataSet]][["autosomes"]] ], na.rm=TRUE) / sum(refCov[ chrIndices[[dataSet]][["autosomes"]] ], na.rm=TRUE)
  obsXAltRefRatio[[dataSet]] = sum(altCov[ chrIndices[[dataSet]][["chrX"]] ], na.rm=TRUE) / sum(refCov[ chrIndices[[dataSet]][["chrX"]] ], na.rm=TRUE)

    # Correction of refcov
  if (obsAutosomeAltRefRatio[[dataSet]]  < 1) {
    adjRefCov = refCov * obsAutosomeAltRefRatio[[dataSet]]
    adjAltCov = altCov
  } else {
    adjAltCov = altCov * (1/obsAutosomeAltRefRatio[[dataSet]])
    adjRefCov = refCov 
  }
  adjTotCov = adjAltCov + adjRefCov
  adjAltProp = adjAltCov / adjTotCov
  adjRefProp = 1-adjAltProp
  
  peakLevelAScovAndRatios[[dataSet]] = cbind(peakLevelASdata[[dataSet]][,1:3],
                                             SNPcounts, altCov, refCov, totCov, altProp, refProp,
                                             adjAltCov, adjRefCov, adjTotCov, adjAltProp, adjRefProp)
  colnames(peakLevelAScovAndRatios[[dataSet]]) = c("chr", "start", "end", 
                                                   "SNPcounts", "altCov", "refCov", "totCov", "altProp", "refProp",
                                                   "adjAltCov", "adjRefCov", "adjTotCov", "adjAltProp", "adjRefProp")
  rownames(peakLevelAScovAndRatios[[dataSet]]) = rownames(peakLevelASdata[[dataSet]])
}







######################################################################################################################
######################################################################################################################
# SNP coverage

SNPcovDir = file.path(outputDir, "SNPcoverage")
if ( ! file.exists( SNPcovDir ) ) {
  if ( ! dir.create( SNPcovDir, showWarnings=FALSE, recursive=TRUE )  ) SNPcovDir = outputDir
}

###########################################################
# Overall averaage Sp:B6 coverage ratios

SNPcovRatioDir = file.path(SNPcovDir, "overallAverageCoverageRatios")
if ( ! file.exists( SNPcovRatioDir ) ) {
  if ( ! dir.create( SNPcovRatioDir, showWarnings=FALSE, recursive=TRUE )  ) SNPcovRatioDir = SNPcovDir
}

# print(obsAutosomeAltRefRatio)
# print(obsXAltRefRatio)

altRefRatioTable = c()
for (dataSet in dataSets) {
  altRefRatioTable = rbind(altRefRatioTable,
                             c(dataSet, obsAutosomeAltRefRatio[[dataSet]],  obsXAltRefRatio[[dataSet]]))
}
colnames(altRefRatioTable) = c("dataSet", "meanAutosomalSpB6coverageRatio", "meanChrXSpB6coverageRatio")

fileName = file.path(SNPcovRatioDir, paste(projectName, ".SpB6coverageRatios.tsv", sep=""))
write.table(altRefRatioTable, fileName, row.names=FALSE, col.names=TRUE, sep="\t")  


##############
# Bar plots

dat = as.data.frame(altRefRatioTable)
# sapply(dat, class)
# sapply(dat, typeof)
# sapply(dat, mode)

dat[, c(2,3)] <- sapply(dat[, c(2,3)], function(x) as.numeric(as.character(x)))
dat[, 1] <- sapply(dat[, 1], function(x) sub("_CTCF", "", sub("ChIPseq_Patski_", "", x)))

# Must factor datasets here so correctly ordered in plot
dat[,1] <- sapply(dat[,1], function(x) as.factor(x))

theCols2use = sampleCols
names(theCols2use) = dat[,1]

ggp = list()

ggp[["autosomes"]] = ggplot(dat, aes(x=dataSet, y=meanAutosomalSpB6coverageRatio, fill=dataSet)) + 
  geom_bar(stat="identity", color="black") +
  scale_y_continuous(limits=c(0, 2.5)) +
  scale_x_discrete(limits=c("WT", "2-4clone", "Del1", "InvDxz4")) +
  scale_fill_manual(values=theCols2use) +
  ggtitle(paste("Autosomes", sep="")) +
  labs(x="") + labs(y="mean Sp/B6 coverage ratio per SNP") +
  theme(legend.position='none')

ggp[["chrX"]] = ggplot(dat, aes(x=dataSet, y=meanChrXSpB6coverageRatio, fill=dataSet)) + 
  geom_bar(stat="identity", color="black") +
  scale_y_continuous(limits=c(0, 2.5)) +
  scale_x_discrete(limits=c("WT", "2-4clone", "Del1", "InvDxz4")) +
  scale_fill_manual(values=theCols2use) +
  ggtitle(paste("chrX", sep="")) +
  labs(x="") + labs(y="mean Xa/Xi coverage ratio per SNP") +   
  theme(legend.position='none')

PDFfilename = file.path(SNPcovRatioDir, 
                        paste(projectName, ".barplots.SpB6coverageRatios.pdf", sep=""))
pdf(PDFfilename, height=6, width=12)

grid.arrange(ggp[["autosomes"]], ggp[["chrX"]], nrow=1, ncol=2, 
             top=textGrob(paste("Overall SNP CTCF peak coverage ratio (Sp/B6)", sep=""),
                           gp=gpar(fontsize=16)))
dev.off()


###########################################################
# Allelic proportions wrt coverage

SNPallelicProportionsWRTcoverageDir = file.path(SNPcovDir, "allelicProportionsWRTcoverage")
if ( ! file.exists( SNPallelicProportionsWRTcoverageDir ) ) {
  if ( ! dir.create( SNPallelicProportionsWRTcoverageDir, showWarnings=FALSE, recursive=TRUE )  ) SNPallelicProportionsWRTcoverageDir = SNPcovDir
}

minCovTableResults = c()

for (dataSet in dataSets) {
  
  for (minCov in c(0, 5, 10)) {
    
    theseMinCovIndices = which(peakLevelAScovAndRatios[[dataSet]]$totCov >= minCov)
  
    cat( "\nminCov = ", minCov, "; n = ", length(theseMinCovIndices), "\n", sep="" )
  
#     cat( "Unadjusted:\n" )
#     cat(dataSet, "\tautosomes\t", mean(peakLevelAScovAndRatios[[dataSet]]$altProp[ intersect(chrIndices[[dataSet]][["autosomes"]], theseMinCovIndices) ], na.rm=TRUE), "\n", sep="")
#     cat(dataSet, "\tchrX\t", mean(peakLevelAScovAndRatios[[dataSet]]$altProp[ intersect(chrIndices[[dataSet]][["chrX"]], theseMinCovIndices) ], na.rm=TRUE), "\n", sep="")
# 
#     cat( "Adjusted:\n" )
#     cat(dataSet, "\tautosomes\t", mean(peakLevelAScovAndRatios[[dataSet]]$adjAltProp[intersect(chrIndices[[dataSet]][["autosomes"]], theseMinCovIndices) ], na.rm=TRUE), "\n", sep="")
#     cat(dataSet, "\tchrX\t", mean(peakLevelAScovAndRatios[[dataSet]]$adjAltProp[ intersect(chrIndices[[dataSet]][["chrX"]], theseMinCovIndices) ], na.rm=TRUE), "\n", sep="")

    minCovTableResults = rbind(minCovTableResults,
                               c(dataSet, minCov, length(theseMinCovIndices),  length(intersect(chrIndices[[dataSet]][["autosomes"]], theseMinCovIndices)),  length(intersect(chrIndices[[dataSet]][["chrX"]], theseMinCovIndices)), 
                                 mean(peakLevelAScovAndRatios[[dataSet]]$altProp[ intersect(chrIndices[[dataSet]][["autosomes"]], theseMinCovIndices) ], na.rm=TRUE),
                                 mean(peakLevelAScovAndRatios[[dataSet]]$altProp[ intersect(chrIndices[[dataSet]][["chrX"]], theseMinCovIndices) ], na.rm=TRUE),
                                 mean(peakLevelAScovAndRatios[[dataSet]]$adjAltProp[intersect(chrIndices[[dataSet]][["autosomes"]], theseMinCovIndices) ], na.rm=TRUE), 
                                 mean(peakLevelAScovAndRatios[[dataSet]]$adjAltProp[ intersect(chrIndices[[dataSet]][["chrX"]], theseMinCovIndices) ], na.rm=TRUE)
                                 ))
  }
}

colnames(minCovTableResults) = c("dataset", "minCoverage", "numberCovered", "autosomesNumberCovered", "chrXNumberCovered",
                                 "meanAutosomalSpProp", "meanXaProp", 
                                 "adjMeanAutosomalSpProp", "adjMeanXaProp" )

fileName = file.path(SNPallelicProportionsWRTcoverageDir, paste(projectName, ".meanSpPropWRTminCoverage.tsv", sep=""))
write.table(minCovTableResults, fileName, row.names=FALSE, col.names=TRUE, sep="\t")  


##############
# Bar plots: peaks vs coverage

dat = as.data.frame(minCovTableResults[,c(1,2,4,5)])
# sapply(dat, class)
# sapply(dat, typeof)
# sapply(dat, mode)

dat[, c(2,3,4)] <- sapply(dat[, c(2,3,4)], function(x) as.numeric(as.character(x)))
dat[, 1] <- sapply(dat[, 1], function(x) sub("_CTCF", "", sub("ChIPseq_Patski_", "", x)))

# Must factor datasets here so correctly ordered in plot
dat[,1] <- sapply(dat[,1], function(x) as.factor(x))

# Must sort by coverage levels here so correctly ordered in plot
dat = with(dat, dat[order(dat$minCoverage),])
dat[,2] <- sapply(dat[,2], function(x) as.factor(x))

theCols2use = sampleCols
names(theCols2use) = unique(dat[,1])

ggp = list()
ggp[["autosomes"]] = ggplot(dat, aes(x=dataset, y=autosomesNumberCovered)) + 
  geom_bar(aes(fill=dataset, alpha=factor(minCoverage)), stat="identity", position="dodge", color="black") +
  scale_x_discrete(limits=c("WT", "2-4clone", "Del1", "InvDxz4")) +
  scale_fill_manual(values=theCols2use) +
  ggtitle(paste("autosomes", sep="")) +
  labs(x="") + labs(y="# covered peaks") +
  theme_bw() + theme( plot.title=element_text(size=16),
                      strip.background  = element_blank(),
                      panel.grid.major = element_line(colour = "grey80"),
                      panel.border = element_blank(),
                      axis.ticks = element_blank(),
                      panel.grid.minor.x=element_blank(),
                      panel.grid.major.x=element_blank(),
                      legend.position='none')
# ggp[["autosomes"]]$labels$fill="Data set"
# ggp[["autosomes"]]$labels$alpha="Minimum coverage"
# ggp[["autosomes"]]

ggp[["chrX"]] = ggplot(dat, aes(x=dataset, y=chrXNumberCovered)) + 
  geom_bar(aes(fill=dataset, alpha=factor(minCoverage)), stat="identity", position="dodge", color="black") +
  scale_x_discrete(limits=c("WT", "2-4clone", "Del1", "InvDxz4")) +
  scale_fill_manual(values=theCols2use) +
  ggtitle(paste("chrX", sep="")) +
  labs(x="") + labs(y="# covered peaks") +
  theme_bw() + theme( plot.title=element_text(size=16),
                      strip.background  = element_blank(),
                      panel.grid.major = element_line(colour = "grey80"),
                      panel.border = element_blank(),
                      axis.ticks = element_blank(),
                      panel.grid.minor.x=element_blank(),
                      panel.grid.major.x=element_blank() )
ggp[["chrX"]]$labels$fill="Data set"
ggp[["chrX"]]$labels$alpha="Minimum coverage"
# ggp[["chrX"]]

# grid.arrange(ggp[["autosomes"]], ggp[["chrX"]], nrow=1, ncol=2, 
#              top=textGrob(paste("Peaks covered at different coverage thresholds", sep=""),
#                            gp=gpar(fontsize=16)))
# 
# grid_arrange_shared_legend(ggp[["autosomes"]], ggp[["chrX"]], ncol=2, nrow=1, position="right",
#                            title=textGrob(paste("Counts of CTCF peak types", sep="")))

mylegend<-g_legend(ggp[["chrX"]])

PDFfilename = file.path(SNPallelicProportionsWRTcoverageDir, 
                        paste(projectName, ".barplots.peaksVsCoverage.pdf", sep=""))
pdf(PDFfilename, height=6, width=10)

grid.arrange(ggp[["autosomes"]] + theme(legend.position="none"), 
             ggp[["chrX"]] + theme(legend.position="none"),
             mylegend, nrow=1, ncol=3, widths=c(10, 10, 5),
             top=textGrob(paste("Peaks covered at different coverage thresholds", sep=""), gp=gpar(fontsize=16)))

dev.off()


##############
# Bar plots: Allelic proportion vs coverage

dat = as.data.frame(minCovTableResults[,c(1,2,8,9)])
# sapply(dat, class)
# sapply(dat, typeof)
# sapply(dat, mode)

dat[, c(2,3,4)] <- sapply(dat[, c(2,3,4)], function(x) as.numeric(as.character(x)))
dat[, 1] <- sapply(dat[, 1], function(x) sub("_CTCF", "", sub("ChIPseq_Patski_", "", x)))

# Must factor datasets here so correctly ordered in plot
dat[,1] <- sapply(dat[,1], function(x) as.factor(x))

# Must sort by coverage levels here so correctly ordered in plot
dat = with(dat, dat[order(dat$minCoverage),])
dat[,2] <- sapply(dat[,2], function(x) as.factor(x))

theCols2use = sampleCols
names(theCols2use) = unique(dat[,1])

ggp = list()
ggp[["autosomes"]] = ggplot(dat, aes(x=dataset, y=adjMeanAutosomalSpProp)) + 
  geom_bar(aes(fill=dataset, alpha=factor(minCoverage)), stat="identity", position="dodge", color="black") +
  scale_x_discrete(limits=c("WT", "2-4clone", "Del1", "InvDxz4")) +
  scale_fill_manual(values=theCols2use) +
  ggtitle(paste("autosomes", sep="")) +
  labs(x="") + labs(y="Sp/[Sp+B6]") +
  theme_bw() + theme( plot.title=element_text(size=16),
                      strip.background  = element_blank(),
                      panel.grid.major = element_line(colour = "grey80"),
                      panel.border = element_blank(),
                      axis.ticks = element_blank(),
                      panel.grid.minor.x=element_blank(),
                      panel.grid.major.x=element_blank(),
                      legend.position='none')
# ggp[["autosomes"]]$labels$fill="Data set"
# ggp[["autosomes"]]$labels$alpha="Minimum coverage"
# ggp[["autosomes"]]

ggp[["chrX"]] = ggplot(dat, aes(x=dataset, y=adjMeanXaProp)) + 
  geom_bar(aes(fill=dataset, alpha=factor(minCoverage)), stat="identity", position="dodge", color="black") +
  scale_x_discrete(limits=c("WT", "2-4clone", "Del1", "InvDxz4")) +
  scale_fill_manual(values=theCols2use) +
  ggtitle(paste("chrX", sep="")) +
  labs(x="") + labs(y="Xa/[Xa+Xi]") +
  theme_bw() + theme( plot.title=element_text(size=16),
                      strip.background  = element_blank(),
                      panel.grid.major = element_line(colour = "grey80"),
                      panel.border = element_blank(),
                      axis.ticks = element_blank(),
                      panel.grid.minor.x=element_blank(),
                      panel.grid.major.x=element_blank() )
ggp[["chrX"]]$labels$fill="Data set"
ggp[["chrX"]]$labels$alpha="Minimum coverage"
# ggp[["chrX"]]

# grid.arrange(ggp[["autosomes"]], ggp[["chrX"]], nrow=1, ncol=2, 
#              top=textGrob(paste("Allelic proportions at different coverage levels", sep=""),
#                            gp=gpar(fontsize=16)))
# 
# grid_arrange_shared_legend(ggp[["autosomes"]], ggp[["chrX"]], ncol=2, nrow=1, position="right",
#                            title=textGrob(paste("Counts of CTCF peak types", sep="")))

mylegend<-g_legend(ggp[["chrX"]])

PDFfilename = file.path(SNPallelicProportionsWRTcoverageDir, 
                        paste(projectName, ".barplots.allelicProportionVsCoverage.pdf", sep=""))
pdf(PDFfilename, height=6, width=12)

grid.arrange(ggp[["autosomes"]] + theme(legend.position="none"), 
             ggp[["chrX"]] + theme(legend.position="none"),
             mylegend, nrow=1, ncol=3, widths=c(10, 10, 5),
             top=textGrob(paste("Allelic proportions at different coverage levels", sep=""), gp=gpar(fontsize=16)))

dev.off()






######################################################################################################################
######################################################################################################################
# Peak SNP count and coverage plots

###
# Peak SNP count distributions

SNPcountsPerPeakDir = file.path(SNPcovDir, "SNPsPerPeak")
if ( ! file.exists( SNPcountsPerPeakDir ) ) {
  if ( ! dir.create( SNPcountsPerPeakDir, showWarnings=FALSE, recursive=TRUE )  ) SNPcountsPerPeakDir = SNPcovDir
}

for (dataSet in dataSets) {

  cellLine = sub("ChIPseq_Patski_", "", dataSet)
  cellLine = sub("_CTCF", "", cellLine)
  cellLine = sub("1", "", cellLine)
  
  PDFfilename = file.path(SNPcountsPerPeakDir, 
                          paste(projectName, ".", dataSet, ".histogram.SNPcountsPerPeak.pdf", sep=""))
  pdf(PDFfilename, height=6, width=10)
  par(mfrow=c(1,2), oma=c(0,0,2.5,0))
  # layout(matrix(c(1,2,3,4),2,2,TRUE))
  # par(oma=c(0,0,2.5,0))

  #for (whichChrs in c("genome-wide", "autosomes", "chrX") ) {
  for (whichChrs in c("autosomes", "chrX") ) {
    
    if (whichChrs == "genome-wide") {
      SNPcountVector = peakLevelAScovAndRatios[[dataSet]]$SNPcounts
    } else {
      SNPcountVector = peakLevelAScovAndRatios[[dataSet]]$SNPcounts[chrIndices[[dataSet]][[whichChrs]]]
    }
    
    ret = hist(SNPcountVector, col=sampleCols[cellLine],
               freq=TRUE, xlim=c(0, ceiling(max(SNPcountVector, na.rm=TRUE)) ), 
               breaks=seq(-0.5, ceiling(max(SNPcountVector, na.rm=TRUE))+0.5, 1),
               main=paste(whichChrs, " (n=", length(SNPcountVector), ")", sep=""),
               xlab="SNPs per peak")
    
    pieValues = c(ret$counts[1], sum(ret$counts[-1]) )
    cols <- c("grey90", "grey50")
    percentlabels <- round(100*pieValues/sum(pieValues), 1)
    pielabels <- paste(percentlabels, "%", sep="")
    #     subPlotRef = subplot(pie(pieValues, main="Genes with SNPs", col=cols, labels=pielabels, cex=0.8), 
    #                          "top", size=c(3,3), inset=c(0, 0.2))
    #     op <- par(no.readonly=TRUE)
    #     par(subPlotRef)
    #     legend("topright", c("0 SNPs", "1+ SNPs"), 
    #            cex=0.8, fill=cols, inset=c(-0.1,0), xpd=TRUE)
    #     par(op)
    subplot( pie(pieValues, main="Peaks with SNPs", col=cols, labels=pielabels, cex=0.8, cex.main=0.8), 
             "top", size=c(2,2), inset=c(0, 0.2))
    legend("topright", c("0 SNPs", "1+ SNPs"), cex=0.8, fill=cols, inset=c(0.15, 0.15))
    
  }
  title(main=paste(dataSet, ": distribution of SNP counts within peaks", sep=""), outer=TRUE)
  dev.off()
}


###
# Peak SNP count distributions for peaks with at least one peak

SNPcountsPerPeakDir = file.path(SNPcovDir, "SNPsPerPeak_alop")
if ( ! file.exists( SNPcountsPerPeakDir ) ) {
  if ( ! dir.create( SNPcountsPerPeakDir, showWarnings=FALSE, recursive=TRUE )  ) SNPcountsPerPeakDir = SNPcovDir
}

for (dataSet in dataSets) {
  
  cellLine = sub("ChIPseq_Patski_", "", dataSet)
  cellLine = sub("_CTCF", "", cellLine)
  cellLine = sub("1", "", cellLine)
  
  PDFfilename = file.path(SNPcountsPerPeakDir, 
                          paste(projectName, ".", dataSet, ".histogram.SNPcountsPerPeak.pdf", sep=""))
  pdf(PDFfilename, height=6, width=10)
  par(mfrow=c(1,2), oma=c(0,0,2.5,0))
  # layout(matrix(c(1,2,3,4),2,2,TRUE))
  # par(oma=c(0,0,2.5,0))
  
  #for (whichChrs in c("genome-wide", "autosomes", "chrX") ) {
  for (whichChrs in c("autosomes", "chrX") ) {
    
    if (whichChrs == "genome-wide") {
      SNPcountVector = peakLevelAScovAndRatios[[dataSet]]$SNPcounts
    } else {
      SNPcountVector = peakLevelAScovAndRatios[[dataSet]]$SNPcounts[chrIndices[[dataSet]][[whichChrs]]]
    }
    # at least one peak
    SNPcountVector = SNPcountVector[SNPcountVector>0]
    
    ret = hist(SNPcountVector, col=sampleCols[cellLine],
               freq=TRUE, xlim=c(0, ceiling(max(SNPcountVector, na.rm=TRUE)) ), 
               breaks=seq(-0.5, ceiling(max(SNPcountVector, na.rm=TRUE))+0.5, 1),
               main=paste(whichChrs, " (n=", length(SNPcountVector), ")", sep=""),
               xlab="SNPs per peak")
    
    pieValues = c(sum(ret$counts[1:4]), sum(ret$counts[5:length(ret$counts)]) )
    cols <- c("grey90", "grey50")
    percentlabels <- round(100*pieValues/sum(pieValues), 1)
    pielabels <- paste(percentlabels, "%", sep="")
    #     subPlotRef = subplot(pie(pieValues, main="Genes with SNPs", col=cols, labels=pielabels, cex=0.8), 
    #                          "top", size=c(3,3), inset=c(0, 0.2))
    #     op <- par(no.readonly=TRUE)
    #     par(subPlotRef)
    #     legend("topright", c("0 SNPs", "1+ SNPs"), 
    #            cex=0.8, fill=cols, inset=c(-0.1,0), xpd=TRUE)
    #     par(op)
    subplot( pie(pieValues, main="Peaks with SNPs", col=cols, labels=pielabels, cex=0.8, cex.main=0.8), 
             "top", size=c(2,2), inset=c(0, 0.2))
    legend("topright", c("1-4 SNPs", "5+ SNPs"), cex=0.8, fill=cols, inset=c(0.15, 0.15))
    
  }
  title(main=paste(dataSet, ": distribution of SNP counts within peaks with at least one SNP", sep=""), outer=TRUE)
  dev.off()
}


###
# SNP coverage distributions
# for peaks with SNPs (SNPcountVector>0)

SNPcoveragePerPeakDir = file.path(SNPcovDir, "SNPcoveragePerPeak")
if ( ! file.exists( SNPcoveragePerPeakDir ) ) {
  if ( ! dir.create( SNPcoveragePerPeakDir, showWarnings=FALSE, recursive=TRUE )  ) SNPcoveragePerPeakDir = SNPcovDir
}

for (dataSet in dataSets) {
  
  cellLine = sub("ChIPseq_Patski_", "", dataSet)
  cellLine = sub("_CTCF", "", cellLine)
  cellLine = sub("1", "", cellLine)
  
  # plot distribution
  PDFfilename = file.path(SNPcoveragePerPeakDir, 
                          paste(projectName, ".", dataSet, ".histogram.SNPtotalReadCoveragePerPeak.pdf", sep=""))
  pdf(PDFfilename, height=6, width=10)
  par(mfrow=c(1,2), oma=c(0,0,2.5,0))
  # layout(matrix(c(1,2,3,4),2,2,TRUE))
  # par(oma=c(0,0,2.5,0))

  #for (whichChrs in c("genome-wide", "autosomes", "chrX") ) {
  for (whichChrs in c("autosomes", "chrX") ) {
    
    if (whichChrs == "genome-wide") {
      SNPcountVector = peakLevelAScovAndRatios[[dataSet]]$SNPcounts
      SNPcovVector = peakLevelAScovAndRatios[[dataSet]]$totCov[SNPcountVector>0]
    } else {
      SNPcountVector = peakLevelAScovAndRatios[[dataSet]]$SNPcounts[chrIndices[[dataSet]][[whichChrs]]]
      SNPcovVector = peakLevelAScovAndRatios[[dataSet]]$totCov[chrIndices[[dataSet]][[whichChrs]]][SNPcountVector>0]
    }
    # SNPcountVector[SNPcovVector==0] # NAs only i.e. Only have zero coverage where no SNPs
    
    ret = hist(SNPcovVector, col=sampleCols[cellLine],
               freq=TRUE, xlim=c(0, quantile(SNPcovVector, 0.99, na.rm=TRUE)), 
               breaks=c(0, 2^seq(0, ceiling(max(log2(SNPcovVector), na.rm=TRUE)), 1)),
               main=paste(whichChrs, " (n=", length(SNPcovVector), ")", sep=""),
               xlab="SNP coverage per peak (log2)")
    
    pieValues = c(sum(SNPcovVector<5 | is.na(SNPcovDir), na.rm=TRUE), sum(SNPcovVector>=5, na.rm=TRUE))
    cols <- c("grey90", "grey50")
    percentlabels <- round(100*pieValues/length(SNPcovVector), 1)
    pielabels <- paste(percentlabels, "%", sep="")
    subplot( pie(pieValues, main="Peaks with 5X", col=cols, labels=pielabels, cex=0.8, cex.main=0.8), 
             "top", size=c(2,2), inset=c(0, 0.2))
    legend("topright", c("<5X", ">=5X"), cex=0.8, fill=cols, inset=c(0.15, 0.15))
  }     
  title(main=paste(dataSet, ": distribution of SNP coverage within peaks", sep=""), outer=TRUE)
  dev.off()
}











######################################################################################################################
######################################################################################################################
# Allelic Ratios

ASEdir = file.path(outputDir, "allelicRatios")
if ( ! file.exists( ASEdir ) ) {
  if ( ! dir.create( ASEdir, showWarnings=FALSE, recursive=TRUE )  ) ASEdir = outputDir
}

minCov = 5

minCovIndices = list()
for (dataSet in dataSets) {
  minCovIndices[[dataSet]] = which(peakLevelAScovAndRatios[[dataSet]]$totCov >= minCov) 
}

# ###
# # Allelic proportion distributions
# for (dataSet in dataSets) {
#   
#   cellLine = sub("ChIPseq_Patski_", "", dataSet)
#   cellLine = sub("_CTCF", "", cellLine)
#   cellLine = sub("1", "", cellLine)
#   
#   # plot distribution
#   PDFfilename = file.path(ASEdir, 
#                           paste(projectName, ".", dataSet, ".histogram.AllelicRatios.minCov=", minCov, ".pdf", sep=""))
#   pdf(PDFfilename, height=6, width=10)
#   par(mfrow=c(1,2), oma=c(0,0,2.5,0))
#   # layout(matrix(c(1,2,3,4),2,2,TRUE))
#   # par(oma=c(0,0,2.5,0))
#   
#   #for (whichChrs in c("genome-wide", "autosomes", "chrX") ) {
#   for (whichChrs in c("autosomes", "chrX") ) {
# 
#     if (whichChrs == "genome-wide") {
#       allelicPropVector = peakLevelAScovAndRatios[[dataSet]]$altProp[ minCovIndices[[dataSet]] ]
#     } else {
#       allelicPropVector = peakLevelAScovAndRatios[[dataSet]]$altProp[ intersect(chrIndices[[dataSet]][[whichChrs]], minCovIndices[[dataSet]]) ]
#     }
#     
#     ret = hist(allelicPropVector, col=sampleCols[cellLine],
#                freq=TRUE, xlim =c(0, 1.05), 
#                breaks=seq(-0.05, 1.05, 0.1),
#                main=paste(whichChrs, " (n=", length(allelicPropVector), ")", sep=""),
#                xlab=ifelse(whichChrs=="chrX", "(Xa/[Xa+Xi])", "(Sp/[Sp+B6])"))
#   }     
#   title(main=paste(dataSet, ": distribution of allelic ratios (minCov=", minCov, ")", sep=""), outer=TRUE)
#   dev.off()  
# }
# 
# 
# ###
# # Allelic proportion distributions
# for (dataSet in dataSets) {
#   
#   cellLine = sub("ChIPseq_Patski_", "", dataSet)
#   cellLine = sub("_CTCF", "", cellLine)
#   cellLine = sub("1", "", cellLine)
#   
#   # plot distribution
#   PDFfilename = file.path(ASEdir, 
#                           paste(projectName, ".", dataSet, ".histogram.AdjustedAllelicRatios.minCov=", minCov, ".pdf", sep=""))
#   pdf(PDFfilename, height=6, width=10)
#   par(mfrow=c(1,2), oma=c(0,0,2.5,0))
#   # layout(matrix(c(1,2,3,4),2,2,TRUE))
#   # par(oma=c(0,0,2.5,0))
#   
#   #for (whichChrs in c("genome-wide", "autosomes", "chrX") ) {
#   for (whichChrs in c("autosomes", "chrX") ) {
#     # cellLines = c("PatskiDel1", "PatskiDel2", "PatskiDel5", "PatskiWT")
#     
#     if (whichChrs == "genome-wide") {
#       allelicPropVector = peakLevelAScovAndRatios[[dataSet]]$adjAltProp[ minCovIndices[[dataSet]] ]
#     } else {
#       allelicPropVector = peakLevelAScovAndRatios[[dataSet]]$adjAltProp[ intersect(chrIndices[[dataSet]][[whichChrs]], minCovIndices[[dataSet]]) ]
#     }
#     
#     ret = hist(allelicPropVector, col=sampleCols[cellLine],
#                freq=TRUE, xlim =c(0, 1.05), 
#                breaks=seq(-0.05, 1.05, 0.1),
#                main=paste(whichChrs, " (n=", length(allelicPropVector), ")", sep=""),
#                xlab=ifelse(whichChrs=="chrX", "(Xa/[Xa+Xi])", "(Sp/[Sp+B6])"))
#   }     
#   title(main=paste(dataSet, ": distribution of allelic ratios (minCov=", minCov, ")", sep=""), outer=TRUE)
#   dev.off()  
# }


###
# Allelic proportion distributions - Adjusted vs unadjusted - overlaid density plots

allelicPropDistDir = file.path(ASEdir, "adjVsUnadjDensityPlots")
if ( ! file.exists( allelicPropDistDir ) ) {
  if ( ! dir.create( allelicPropDistDir, showWarnings=FALSE, recursive=TRUE )  ) allelicPropDistDir = ASEdir
}

for (dataSet in dataSets) {
  
  cellLine = sub("ChIPseq_Patski_", "", dataSet)
  cellLine = sub("_CTCF", "", cellLine)
  cellLine = sub("1", "", cellLine)
  
  allelicPropVector=list()
  ggplotList = list()
  
  #for (whichChrs in c("genome-wide", "autosomes", "chrX") ) {
  for (whichChrs in c("autosomes", "chrX") ) {
    
    if (whichChrs == "genome-wide") {
      allelicPropVector[["unadjusted"]] = as.data.frame(peakLevelAScovAndRatios[[dataSet]]$altProp[ minCovIndices[[dataSet]] ])
      allelicPropVector[["adjusted"]] = as.data.frame(peakLevelAScovAndRatios[[dataSet]]$adjAltProp[ minCovIndices[[dataSet]] ])
    } else {
      allelicPropVector[["unadjusted"]] = as.data.frame(peakLevelAScovAndRatios[[dataSet]]$altProp[ intersect(chrIndices[[dataSet]][[whichChrs]], minCovIndices[[dataSet]]) ])
      allelicPropVector[["adjusted"]] = as.data.frame(peakLevelAScovAndRatios[[dataSet]]$adjAltProp[ intersect(chrIndices[[dataSet]][[whichChrs]], minCovIndices[[dataSet]]) ])
    }
    names(allelicPropVector[["unadjusted"]]) = "ratio"
    names(allelicPropVector[["adjusted"]]) = "ratio"
    
    list.names <- names(allelicPropVector)
    lns <- sapply(allelicPropVector, nrow)
    dat <- as.data.frame(do.call("rbind", allelicPropVector))
    colnames(dat)[1] = "ratio"
    dat$group <- rep(list.names, lns)

    theCols2use = c(sampleCols[cellLine], sampleColsLightD8[cellLine])
    names(theCols2use) = list.names

    ggplotList[[whichChrs]] <- ggplot(dat, aes(ratio, fill=group)) + 
      geom_density(alpha=0.6, size=0.25) +
      ggtitle(paste(whichChrs, "\n(n=", lns[1], ")", sep="")) +
      scale_fill_manual(values=theCols2use) +
      labs(x=ifelse(whichChrs=="chrX", "(Xa/[Xa+Xi])", "(Sp/[Sp+B6])")) + 
      theme(legend.position="top") +
      theme(legend.title=element_blank())
  }

  PDFfilename = file.path(allelicPropDistDir, 
                          paste(projectName, ".", dataSet, ".densityOverlaid.unadjustedVsAdjustedAllelicRatios.minCov=", minCov, ".pdf", sep=""))
  pdf(PDFfilename, height=6, width=10)
  #par(mfrow=c(1,2), oma=c(0,0,2.5,0))
  
  grid.arrange(ggplotList[[1]], ggplotList[[2]], nrow=1, ncol=2, 
               top=textGrob(paste(cellLine, " distributions of ", list.names[1], " and ", list.names[2], " allelic ratios (minCov=", minCov, ")", sep=""),
                             gp=gpar(fontsize=16)))
  dev.off()
}



###
# Allelic proportion distributions for all chromosomes - Adjusted vs unadjusted - overlaid density plots

ggplotList = list()

for (dataSet in dataSets) {
  
  cellLine = sub("ChIPseq_Patski_", "", dataSet)
  cellLine = sub("_CTCF", "", cellLine)
  cellLine = sub("1", "", cellLine)
  
  allelicPropVector=list()
  
  for (whichChrs in c(paste("chr", seq(1,19,1), sep=""), "chrX") ) {
    
    if (whichChrs == "genome-wide") {
      allelicPropVector[["unadjusted"]] = as.data.frame(peakLevelAScovAndRatios[[dataSet]]$altProp[ minCovIndices[[dataSet]] ])
      allelicPropVector[["adjusted"]] = as.data.frame(peakLevelAScovAndRatios[[dataSet]]$adjAltProp[ minCovIndices[[dataSet]] ])
    } else {
      allelicPropVector[["unadjusted"]] = as.data.frame(peakLevelAScovAndRatios[[dataSet]]$altProp[ intersect(chrIndices[[dataSet]][[whichChrs]], minCovIndices[[dataSet]]) ])
      allelicPropVector[["adjusted"]] = as.data.frame(peakLevelAScovAndRatios[[dataSet]]$adjAltProp[ intersect(chrIndices[[dataSet]][[whichChrs]], minCovIndices[[dataSet]]) ])
    }
    names(allelicPropVector[["unadjusted"]]) = "ratio"
    names(allelicPropVector[["adjusted"]]) = "ratio"
    
    list.names <- names(allelicPropVector)
    lns <- sapply(allelicPropVector, nrow)
    dat <- as.data.frame(do.call("rbind", allelicPropVector))
    colnames(dat)[1] = "ratio"
    dat$group <- rep(list.names, lns)
    
    theCols2use = c(sampleCols[cellLine], sampleColsLightD8[cellLine])
    names(theCols2use) = list.names
    
    ggplotList[[paste(dataSet, "_", whichChrs, sep="")]] <- ggplot(dat, aes(ratio, fill=group)) + 
      geom_density(alpha=0.6, size=0.25) +
      ggtitle(paste(whichChrs, " (n=", lns[1], ")", sep="")) +
      scale_fill_manual(values=theCols2use) +
      labs(x=ifelse(whichChrs=="chrX", "(Xa/[Xa+Xi])", "(Sp/[Sp+B6])")) + 
      guides(fill=FALSE)
  }
}

PDFfilename = file.path(allelicPropDistDir, 
                        paste(projectName, ".densityOverlaid.unadjustedVsAdjustedAllelicRatios.allChromosomes.minCov=", minCov, ".pdf", sep=""))
pdf(PDFfilename, height=4*3+1, width=20*3+1)
grid.arrange(do.call("arrangeGrob", c(ggplotList, nrow=4)))
dev.off()




###
# Allelic proportion distributions V2 - Adjusted vs unadjusted - overlaid density plots

allelicPropDistDir = file.path(ASEdir, "adjVsUnadjDensityPlotsV2")
if ( ! file.exists( allelicPropDistDir ) ) {
  if ( ! dir.create( allelicPropDistDir, showWarnings=FALSE, recursive=TRUE )  ) allelicPropDistDir = ASEdir
}

# vector of spretus minima
autosomalMinima = list()
autosomalMaxima = list()
autosomalProps = list()
chrXmaxima = list()
chrXSubProps = list()

for (dataSet in dataSets) {
  # dataSet = "ChIPseq_Patski_WT_CTCF"
  
  cellLine = sub("ChIPseq_Patski_", "", dataSet)
  cellLine = sub("_CTCF", "", cellLine)
  cellLine = sub("1", "", cellLine)
  
  allelicPropVector=list()
  ggplotList = list()
  
  autosomalMinima[[cellLine]] = list()
  autosomalMaxima[[cellLine]] = list()
  chrXmaxima[[cellLine]] = list()
  
  autosomalProps[[cellLine]] = list()
  chrXSubProps[[cellLine]] = list()
  
  # Density plots
  #for (whichChrs in c("genome-wide", "autosomes", "chrX") ) {
  for (whichChrs in c("autosomes", "diplosomes", "chrX") ) {
    # whichChrs = "autosomes"
    # whichChrs = "diplosomes"
    
    if (whichChrs == "genome-wide") {
      allelicPropVector[["unadjusted"]] = as.data.frame(peakLevelAScovAndRatios[[dataSet]]$altProp[ minCovIndices[[dataSet]] ])
      allelicPropVector[["adjusted"]] = as.data.frame(peakLevelAScovAndRatios[[dataSet]]$adjAltProp[ minCovIndices[[dataSet]] ])
    } else {
      allelicPropVector[["unadjusted"]] = as.data.frame(peakLevelAScovAndRatios[[dataSet]]$altProp[ intersect(chrIndices[[dataSet]][[whichChrs]], minCovIndices[[dataSet]]) ])
      allelicPropVector[["adjusted"]] = as.data.frame(peakLevelAScovAndRatios[[dataSet]]$adjAltProp[ intersect(chrIndices[[dataSet]][[whichChrs]], minCovIndices[[dataSet]]) ])
    }
    names(allelicPropVector[["unadjusted"]]) = "ratio"
    names(allelicPropVector[["adjusted"]]) = "ratio"
    
    list.names <- names(allelicPropVector)
    lns <- sapply(allelicPropVector, nrow)
    dat <- as.data.frame(do.call("rbind", allelicPropVector))
    colnames(dat)[1] = "ratio"
    dat$group <- rep(list.names, lns)
    
    theCols2use = c(sampleCols[cellLine], sampleColsLightD8[cellLine])
    names(theCols2use) = list.names
    
    ggplotList[[whichChrs]] <- ggplot(dat, aes(ratio, fill=group)) + 
      geom_density(alpha=0.6, size=0.25) +
      ggtitle(paste(whichChrs, "\n(n=", lns[1], ")", sep="")) +
      scale_fill_manual(values=theCols2use) +
      labs(x=ifelse(whichChrs=="chrX", "(Xa/[Xa+Xi])", "(Sp/[Sp+B6])")) + 
      theme(legend.position="top") +
      theme(legend.title=element_blank())
    
    # Get B6 and Sp-specifc autosomal minima & proportions below and above, resp.
    if ( whichChrs == "autosomes" || whichChrs == "diplosomes"  ) {
      
      # autosomal minima
      ggplotData = print(ggplotList[[whichChrs]])
      # names(ggplotData)
      # str(ggplotData[["data"]])
      # attributes(ggplotData[["data"]][[1]])
      yVals = ggplotData[["data"]][[1]]$y
      xVals = ggplotData[["data"]][[1]]$x
      idx1 = which.min(yVals[xVals<0.5])
      yMin1 = yVals[xVals<0.5][idx1]
      xMin1 = xVals[xVals<0.5][idx1]
      idx2 = which.min(yVals[xVals>0.5])
      yMin2 = yVals[xVals>0.5][idx2]
      xMin2 = xVals[xVals>0.5][idx2]
      idxMax = which.max(yVals)
      yMax = yVals[idxMax]
      xMax = xVals[idxMax]
      
      autosomalMinima[[cellLine]][[whichChrs]] = list()
      autosomalMinima[[cellLine]][[whichChrs]][["xMin1"]] = xMin1
      autosomalMinima[[cellLine]][[whichChrs]][["yMin1"]] = yMin1
      autosomalMinima[[cellLine]][[whichChrs]][["xMin2"]] = xMin2
      autosomalMinima[[cellLine]][[whichChrs]][["yMin1"]] = yMin2
      autosomalMaxima[[cellLine]][[whichChrs]][["xMax"]] = xMax
      autosomalMaxima[[cellLine]][[whichChrs]][["yMax"]] = yMax
      
      # B6 and Sp-specific proportions 
      B6autoProp = sum(ggplotData[["data"]][[1]]$y[ggplotData[["data"]][[1]]$x<xMin1]) / sum(ggplotData[["data"]][[1]]$y)
      SpAutoProp = sum(ggplotData[["data"]][[1]]$y[ggplotData[["data"]][[1]]$x>xMin2]) / sum(ggplotData[["data"]][[1]]$y)
      SubMaxAutoProp = sum(ggplotData[["data"]][[1]]$y[ggplotData[["data"]][[1]]$x<xMax]) / sum(ggplotData[["data"]][[1]]$y)
      
      autosomalProps[[cellLine]][[whichChrs]] = list()
      autosomalProps[[cellLine]][[whichChrs]][["B6autoProp"]] = B6autoProp
      autosomalProps[[cellLine]][[whichChrs]][["SpAutoProp"]] = SpAutoProp
      autosomalProps[[cellLine]][[whichChrs]][["SubMaxAutoProp"]] = SubMaxAutoProp
      
      yMid = ggplotData[["panel"]]$ranges[[1]]$y.range[1] +
        (ggplotData[["panel"]]$ranges[[1]]$y.range[2] - ggplotData[["panel"]]$ranges[[1]]$y.range[1])*0.65
      
      ggplotList[[whichChrs]] = ggplotList[[whichChrs]] +
        geom_vline(xintercept=xMin1, linetype="dashed") +
        geom_vline(xintercept=xMin2, linetype="dashed") +
        geom_vline(xintercept=xMax, linetype="dotted") +
        annotate(geom="text", label=paste("(", round(xMin1, 2), "; ", round(yMin1, 2), "); Pr(Y|X<", round(xMin1, 2), ") = ", round(B6autoProp, 2), sep=""), x=xMin1, y=yMid,  angle = 90, vjust=-1) +
        annotate(geom="text", label=paste("(", round(xMin2, 2), "; ", round(yMin2, 2), "); Pr(Y|X>", round(xMin2, 2), ") = ", round(SpAutoProp, 2), sep=""), x=xMin2, y=yMid,  angle = 90, vjust=+1.5) +
        annotate(geom="text", label=paste("(", round(xMax, 2), "; ", round(yMax, 2), "); Pr(Y|X<", round(xMax, 2), ") = ", round(SubMaxAutoProp, 2), sep=""), x=xMax, y=yMid,  angle = 90, vjust=-1)
    }    

    if ( whichChrs == "chrX" ) {
      xMin1 = autosomalMinima[[cellLine]][["diplosomes"]][["xMin1"]]
      yMin1 = autosomalMinima[[cellLine]][["diplosomes"]][["yMin1"]] 
      xMin2 = autosomalMinima[[cellLine]][["diplosomes"]][["xMin2"]] 
      yMin2 = autosomalMinima[[cellLine]][["diplosomes"]][["yMin1"]] 

      ggplotData = print(ggplotList[[whichChrs]])
      yVals = ggplotData[["data"]][[1]]$y
      xVals = ggplotData[["data"]][[1]]$x
      idxMax = which.max(yVals)
      yMax = yVals[idxMax]
      xMax = xVals[idxMax]

      chrXmaxima[[cellLine]][[whichChrs]][["xMax"]] = xMax
      chrXmaxima[[cellLine]][[whichChrs]][["yMax"]] = yMax
      
      yMid = ggplotData[["panel"]]$ranges[[1]]$y.range[1] +
        (ggplotData[["panel"]]$ranges[[1]]$y.range[2] - ggplotData[["panel"]]$ranges[[1]]$y.range[1])*0.65
      
      # B6 and Sp-sub proportions 
      B6subProp = sum(ggplotData[["data"]][[1]]$y[ggplotData[["data"]][[1]]$x<xMin1]) / sum(ggplotData[["data"]][[1]]$y)
      SpSubProp = sum(ggplotData[["data"]][[1]]$y[ggplotData[["data"]][[1]]$x>xMin2]) / sum(ggplotData[["data"]][[1]]$y)
      SubMaxAutoProp = sum(ggplotData[["data"]][[1]]$y[ggplotData[["data"]][[1]]$x<xMax]) / sum(ggplotData[["data"]][[1]]$y)
      
      chrXSubProps[[cellLine]][["B6subProp"]] = B6subProp
      chrXSubProps[[cellLine]][["SpSubProp"]] = SpSubProp
      chrXSubProps[[cellLine]][["SubMaxProp"]] = SubMaxAutoProp
      
      ggplotList[[whichChrs]] = ggplotList[[whichChrs]] +
        geom_vline(xintercept=xMin1, linetype="dashed") +
        geom_vline(xintercept=xMax, linetype="dotted") +
        geom_vline(xintercept=xMin2, linetype="dashed") +
        annotate(geom="text", label=paste("(", round(xMin1, 2), "; ", round(yMin1, 2), "); Pr(Y|X<", round(xMin1, 2), ") = ", round(B6subProp, 2), sep=""), x=xMin1, y=yMid,  angle = 90, vjust=-1) +
        annotate(geom="text", label=paste("(", round(xMin2, 2), "; ", round(yMin2, 2), "); Pr(Y|X>", round(xMin2, 2), ") = ", round(SpSubProp, 2), sep=""), x=xMin2, y=yMid,  angle = 90, vjust=+1.5) +
        annotate(geom="text", label=paste("(", round(xMax, 2), "; ", round(yMax, 2), "); Pr(Y|X<", round(xMax, 2), ") = ", round(SubMaxAutoProp, 2), sep=""), x=xMax, y=yMid,  angle = 90, vjust=-1)
    }    
  }
  
  # Histograms
  #for (whichChrs in c("genome-wide", "autosomes", "chrX") ) {
  for (whichChrs in c("autosomes", "diplosomes", "chrX") ) {
    # whichChrs = "chrX"

    if (whichChrs == "genome-wide") {
      allelicPropVector[["unadjusted"]] = as.data.frame(peakLevelAScovAndRatios[[dataSet]]$altProp[ minCovIndices[[dataSet]] ])
      allelicPropVector[["adjusted"]] = as.data.frame(peakLevelAScovAndRatios[[dataSet]]$adjAltProp[ minCovIndices[[dataSet]] ])
    } else {
      allelicPropVector[["unadjusted"]] = as.data.frame(peakLevelAScovAndRatios[[dataSet]]$altProp[ intersect(chrIndices[[dataSet]][[whichChrs]], minCovIndices[[dataSet]]) ])
      allelicPropVector[["adjusted"]] = as.data.frame(peakLevelAScovAndRatios[[dataSet]]$adjAltProp[ intersect(chrIndices[[dataSet]][[whichChrs]], minCovIndices[[dataSet]]) ])
    }
    names(allelicPropVector[["unadjusted"]]) = "ratio"
    names(allelicPropVector[["adjusted"]]) = "ratio"    
        
    list.names <- names(allelicPropVector)
    lns <- sapply(allelicPropVector, nrow)
    dat <- as.data.frame(do.call("rbind", allelicPropVector))
    colnames(dat)[1] = "ratio"
    dat$group <- rep(list.names, lns)
    
    theCols2use = c(sampleCols[cellLine], sampleColsLightD8[cellLine])
    names(theCols2use) = list.names
    
    ggplotList[[paste(whichChrs, "_hist", sep="")]] <- ggplot(dat, aes(ratio, fill=group)) + 
      geom_histogram(breaks=seq(0, 1, by=0.05), alpha=0.6, position='dodge',
                     aes(y=c(..count..[..group..==1]/sum(..count..[..group..==1]),
                             ..count..[..group..==2]/sum(..count..[..group..==2]))*100)) +
      ggtitle(paste("chrX zoomed\n(n=", lns[1], ")", sep="")) +
      scale_fill_manual(values=theCols2use) +
      labs(x=ifelse(whichChrs=="chrX", "(Xa/[Xa+Xi])", "(Sp/[Sp+B6])"), y="%") + 
      theme(legend.position="top") +
      theme(legend.title=element_blank())
    
    # Get B6 and Sp-specifc autosomal minima
    if ( whichChrs == "autosomes" || whichChrs == "diplosomes"  ) {
      xMin1 = autosomalMinima[[cellLine]][[whichChrs]][["xMin1"]]
      yMin1 = autosomalMinima[[cellLine]][[whichChrs]][["yMin1"]] 
      xMin2 = autosomalMinima[[cellLine]][[whichChrs]][["xMin2"]] 
      yMin2 = autosomalMinima[[cellLine]][[whichChrs]][["yMin1"]] 
      xMax = autosomalMaxima[[cellLine]][[whichChrs]][["xMax"]]
      yMax = autosomalMaxima[[cellLine]][[whichChrs]][["yMax"]]
      
      B6autoProp = autosomalProps[[cellLine]][[whichChrs]][["B6autoProp"]]
      SpAutoProp = autosomalProps[[cellLine]][[whichChrs]][["SpAutoProp"]]
      SubMaxAutoProp = autosomalProps[[cellLine]][[whichChrs]][["SubMaxAutoProp"]]
      
      ggplotData = print(ggplotList[[paste(whichChrs, "_hist", sep="")]])
      yMid = ggplotData[["panel"]]$ranges[[1]]$y.range[1] +
        (ggplotData[["panel"]]$ranges[[1]]$y.range[2] - ggplotData[["panel"]]$ranges[[1]]$y.range[1])*0.65
      
      ggplotList[[paste(whichChrs, "_hist", sep="")]] = ggplotList[[paste(whichChrs, "_hist", sep="")]] +
        geom_vline(xintercept=xMin1, linetype="dashed") +
        geom_vline(xintercept=xMin2, linetype="dashed") +
        geom_vline(xintercept=xMax, linetype="dotted") +
        annotate(geom="text", label=paste("(", round(xMin1, 2), "; ", round(yMin1, 2), "); Pr(Y|X<", round(xMin1, 2), ") = ", round(B6autoProp, 2), sep=""), x=xMin1, y=yMid,  angle = 90, vjust=-1) +
        annotate(geom="text", label=paste("(", round(xMin2, 2), "; ", round(yMin2, 2), "); Pr(Y|X>", round(xMin2, 2), ") = ", round(SpAutoProp, 2), sep=""), x=xMin2, y=yMid,  angle = 90, vjust=+1.5) +
        annotate(geom="text", label=paste("(", round(xMax, 2), "; ", round(yMax, 2), "); Pr(Y|X<", round(xMax, 2), ") = ", round(SubMaxAutoProp, 2), sep=""), x=xMax, y=yMid,  angle = 90, vjust=-1)
    }
    
    if ( whichChrs == "chrX" ) {
      xMin1 = autosomalMinima[[cellLine]][["diplosomes"]][["xMin1"]]
      yMin1 = autosomalMinima[[cellLine]][["diplosomes"]][["yMin1"]] 
      xMin2 = autosomalMinima[[cellLine]][["diplosomes"]][["xMin2"]] 
      yMin2 = autosomalMinima[[cellLine]][["diplosomes"]][["yMin1"]] 
      xMax = chrXmaxima[[cellLine]][[whichChrs]][["xMax"]]
      yMax = chrXmaxima[[cellLine]][[whichChrs]][["yMax"]]
      
      B6subProp = chrXSubProps[[cellLine]][["B6subProp"]]
      SpSubProp = chrXSubProps[[cellLine]][["SpSubProp"]]          
      SubMaxAutoProp = chrXSubProps[[cellLine]][["SubMaxProp"]]
      
      ggplotData = print(ggplotList[[paste(whichChrs, "_hist", sep="")]])
      yMid = ggplotData[["panel"]]$ranges[[1]]$y.range[1] +
        (ggplotData[["panel"]]$ranges[[1]]$y.range[2] - ggplotData[["panel"]]$ranges[[1]]$y.range[1])*0.65
      
      ggplotList[[paste(whichChrs, "_hist", sep="")]] = ggplotList[[paste(whichChrs, "_hist", sep="")]] +
        geom_vline(xintercept=xMin1, linetype="dashed") +
        geom_vline(xintercept=xMin2, linetype="dashed") +
        geom_vline(xintercept=xMax, linetype="dotted") +
        annotate(geom="text", label=paste("(", round(xMin1, 2), "; ", round(yMin1, 2), "); Pr(Y|X<", round(xMin1, 2), ") = ", round(B6autoProp, 2), sep=""), x=xMin1, y=yMid,  angle = 90, vjust=-1) +
        annotate(geom="text", label=paste("(", round(xMin2, 2), "; ", round(yMin2, 2), "); Pr(Y|X>", round(xMin2, 2), ") = ", round(SpAutoProp, 2), sep=""), x=xMin2, y=yMid,  angle = 90, vjust=+1.5) +
        annotate(geom="text", label=paste("(", round(xMax, 2), "; ", round(yMax, 2), "); Pr(Y|X<", round(xMax, 2), ") = ", round(SubMaxAutoProp, 2), sep=""), x=xMax, y=yMid,  angle = 90, vjust=-1)
    }    
  }
  
  PDFfilename = file.path(allelicPropDistDir, 
                          paste(projectName, ".", dataSet, ".densityOverlaid.unadjustedVsAdjustedAllelicRatios.V2.minCov=", minCov, ".pdf", sep=""))
  pdf(PDFfilename, height=12, width=15)
  #par(mfrow=c(1,2), oma=c(0,0,2.5,0))
  
  grid.arrange(ggplotList[[1]], ggplotList[[2]], ggplotList[[3]], 
               ggplotList[[4]], ggplotList[[5]], ggplotList[[6]], 
               nrow=2, ncol=3, 
               top=textGrob(paste(cellLine, " distributions of ", list.names[1], " and ", list.names[2], " allelic ratios (minCov=", minCov, ")", sep=""),
                            gp=gpar(fontsize=16)))
  dev.off()
}






###
# Allelic proportion distributions  - Adjusted histograms

allelicPropDistDir = file.path(ASEdir, "adjHistograms")
if ( ! file.exists( allelicPropDistDir ) ) {
  if ( ! dir.create( allelicPropDistDir, showWarnings=FALSE, recursive=TRUE )  ) allelicPropDistDir = ASEdir
}

ggplotList = list()

for (dataSet in dataSets) {
  # dataSet = "WT_ATAC"
  
  cellLine = sub("ChIPseq_Patski_", "", dataSet)
  cellLine = sub("_CTCF", "", cellLine)
  cellLine = sub("1", "", cellLine)
  
  autosomalMinima[[cellLine]] = list()
  autosomalMaxima[[cellLine]] = list()
  chrXmaxima[[cellLine]] = list()
  
  autosomalProps[[cellLine]] = list()
  chrXSubProps[[cellLine]] = list()
  
  # Histograms
  #for (whichChrs in c("genome-wide", "autosomes", "chrX") ) {
  for (whichChrs in c("autosomes", "diplosomes", "chrX") ) {
    # whichChrs = "chrX"
    
    if (whichChrs == "genome-wide") {
      allelicPropVector= as.data.frame(peakLevelAScovAndRatios[[dataSet]]$adjAltProp[ minCovIndices[[dataSet]] ])
    } else {
      allelicPropVector = as.data.frame(peakLevelAScovAndRatios[[dataSet]]$adjAltProp[ intersect(chrIndices[[dataSet]][[whichChrs]], minCovIndices[[dataSet]]) ])
    }
    names(allelicPropVector) = "ratio"
    
    ggplotList[[paste(dataSet, "_", whichChrs, "_hist", sep="")]] <- ggplot(allelicPropVector, aes(ratio)) + 
      geom_histogram(breaks=seq(0, 1, by=0.05), alpha=0.6, color="black", fill=sampleCols[cellLine],
                     aes(y=(..count../sum(..count..)*100))) +
      ggtitle(paste("chrX (n=", nrow(allelicPropVector), ")", sep="")) +
      labs(x=ifelse(whichChrs=="chrX", "(Xa/[Xa+Xi])", "(Sp/[Sp+B6])"), y="%") + 
      coord_cartesian(ylim=c(0, 15))
  }
}

PDFfilename = file.path(allelicPropDistDir, 
                        paste(projectName, ".histogram.adjustedAllelicRatios.minCov=", minCov, ".pdf", sep=""))
pdf(PDFfilename, height=25, width=15)
#par(mfrow=c(1,2), oma=c(0,0,2.5,0))

grid.arrange(ggplotList[[1]], ggplotList[[2]], ggplotList[[3]], 
             ggplotList[[4]], ggplotList[[5]], ggplotList[[6]], 
             ggplotList[[7]], ggplotList[[8]], ggplotList[[9]], 
             ggplotList[[10]], ggplotList[[11]], ggplotList[[12]], 
             nrow=4, ncol=3, 
             top=textGrob(paste("Distributions of adjusted allelic ratios (minCov=", minCov, ")", sep=""),
                          gp=gpar(fontsize=16)))
dev.off()




###
# Allelic proportion distributions for all chroms - Adjusted histograms

ggplotList = list()

for (dataSet in dataSets) {
  # dataSet = "WT_ATAC"
  
  cellLine = sub("ChIPseq_Patski_", "", dataSet)
  cellLine = sub("_CTCF", "", cellLine)
  cellLine = sub("1", "", cellLine)
  
  autosomalMinima[[cellLine]] = list()
  autosomalMaxima[[cellLine]] = list()
  chrXmaxima[[cellLine]] = list()
  
  autosomalProps[[cellLine]] = list()
  chrXSubProps[[cellLine]] = list()
  
  # Histograms
  for (whichChrs in c(paste("chr", seq(1,19,1), sep=""), "chrX") ) {
    # whichChrs = "chrX"
    
    if (whichChrs == "genome-wide") {
      allelicPropVector= as.data.frame(peakLevelAScovAndRatios[[dataSet]]$adjAltProp[ minCovIndices[[dataSet]] ])
    } else {
      allelicPropVector = as.data.frame(peakLevelAScovAndRatios[[dataSet]]$adjAltProp[ intersect(chrIndices[[dataSet]][[whichChrs]], minCovIndices[[dataSet]]) ])
    }
    names(allelicPropVector) = "ratio"
    
    ggplotList[[paste(dataSet, "_", whichChrs, "_hist", sep="")]] <- ggplot(allelicPropVector, aes(ratio)) + 
      geom_histogram(breaks=seq(0, 1, by=0.05), alpha=0.6, color="black", fill=sampleCols[cellLine],
                     aes(y=(..count../sum(..count..)*100))) +
      ggtitle(paste(whichChrs, " (n=", nrow(allelicPropVector), ")", sep="")) +
      labs(x=ifelse(whichChrs=="chrX", "(Xa/[Xa+Xi])", "(Sp/[Sp+B6])"), y="%") + 
      coord_cartesian(ylim=c(0, 20))
  }
}

PDFfilename = file.path(allelicPropDistDir, 
                        paste(projectName, ".histogram.adjustedAllelicRatios.allChromosomes.minCov=", minCov, ".pdf", sep=""))
pdf(PDFfilename, height=4*3+1, width=20*3+1)
grid.arrange(do.call("arrangeGrob", c(ggplotList, nrow=4)))
dev.off()





###
# Adjusted allelic proportion distributions - WT vs Del - overlaid density plots

allelicPropDistDir = file.path(ASEdir, "overlaidAdjAllelicPropDensityPlots")
if ( ! file.exists( allelicPropDistDir ) ) {
  if ( ! dir.create( allelicPropDistDir, showWarnings=FALSE, recursive=TRUE )  ) allelicPropDistDir = ASEdir
}

allelicPropVector=list()
ggplotList = list()

#for (whichChrs in c("genome-wide", "autosomes", "chrX") ) {
for (whichChrs in c("autosomes", "chrX") ) {
  
  for (dataSet in dataSets) {
    
    cellLine = sub("ChIPseq_Patski_", "", dataSet)
    cellLine = sub("_CTCF", "", cellLine)
    cellLine = sub("1", "", cellLine)
    
    if (whichChrs == "genome-wide") {
      allelicPropVector[[cellLine]] = as.data.frame(peakLevelAScovAndRatios[[dataSet]]$adjAltProp[ minCovIndices[[dataSet]] ])
    } else {
      allelicPropVector[[cellLine]] = as.data.frame(peakLevelAScovAndRatios[[dataSet]]$adjAltProp[ intersect(chrIndices[[dataSet]][[whichChrs]], minCovIndices[[dataSet]]) ])
    }
  }
  
  list.names <- names(allelicPropVector)
  lns <- sapply(allelicPropVector, nrow)
  dat <- as.data.frame(do.call("rbind", allelicPropVector))
  colnames(dat)[1] = "ratio"
  dat$group <- rep(list.names, lns)
  
  ggplotList[[whichChrs]] <- ggplot(dat, aes(ratio, fill=group)) + 
    geom_density(alpha=0.3, size=0.25) +
    ggtitle(paste(whichChrs, "\n(|WT|=", lns["WT"], "; |2-4|=", lns["2-4clone"], ";\n|Del|=", lns["Del"], "; |Inv|=", lns["InvDxz4"], ")", sep="")) +
    scale_fill_manual(values=sampleCols[list.names]) +
    labs(x=ifelse(whichChrs=="chrX", "(Xa/[Xa+Xi])", "(Sp/[Sp+B6])")) + 
    theme(legend.position="top") +
    theme(legend.title=element_blank())
}


PDFfilename = file.path(allelicPropDistDir, 
                        paste(projectName, ".densityOverlaid.WTvsDelAdjustedAllelicRatios.minCov=", minCov, ".pdf", sep=""))
pdf(PDFfilename, height=6, width=10)
#par(mfrow=c(1,2), oma=c(0,0,2.5,0))


grid.arrange(ggplotList[[1]], ggplotList[[2]], nrow=1, ncol=2, 
             top=textGrob(paste("Distributions of ", list.names[1], " and ", list.names[2], " allelic ratios (minCov=", minCov, ")", sep=""),
                       gp=gpar(fontsize=16)))
dev.off()






###
# Adjusted allelic proportion distributions - WT vs Del - overlaid density plots
# V2 autosomes + chrX

allelicPropVector=list()
ggplotList = list()

#for (whichChrs in c("genome-wide", "autosomes", "chrX") ) {
for (whichChrs in c("autosomes", "chrX") ) {
  
  for (dataSet in dataSets) {
    
    cellLine = sub("ChIPseq_Patski_", "", dataSet)
    cellLine = sub("_CTCF", "", cellLine)
    cellLine = sub("1", "", cellLine)
    
    if (whichChrs == "genome-wide") {
      allelicPropVector[[cellLine]] = as.data.frame(peakLevelAScovAndRatios[[dataSet]]$adjAltProp[ minCovIndices[[dataSet]] ])
    } else {
      allelicPropVector[[cellLine]] = as.data.frame(peakLevelAScovAndRatios[[dataSet]]$adjAltProp[ intersect(chrIndices[[dataSet]][[whichChrs]], minCovIndices[[dataSet]]) ])
    }
  }
  
  list.names <- names(allelicPropVector)
  lns <- sapply(allelicPropVector, nrow)
  dat <- as.data.frame(do.call("rbind", allelicPropVector))
  colnames(dat)[1] = "ratio"
  dat$group <- rep(list.names, lns)
  
  ggplotList[[whichChrs]] <- ggplot(dat, aes(ratio, fill=group)) + 
    geom_density(alpha=0.2, size=0.25) +
    ggtitle(paste(whichChrs, sep="")) +
    scale_fill_manual(values=sampleCols[list.names]) +
    labs(x=ifelse(whichChrs=="chrX", "(Xa/[Xa+Xi])", "(Sp/[Sp+B6])")) + 
    theme(legend.position="top") +
    theme(legend.title=element_blank())

  if ( whichChrs == "chrX" ) {
  
    ggplotData = print(ggplotList[[whichChrs]])

    for ( grpId in unique(ggplotData[["data"]][[1]]$group) ) { 
      yVals = ggplotData[["data"]][[1]]$y[ggplotData[["data"]][[1]]$group == grpId]
      xVals = ggplotData[["data"]][[1]]$x[ggplotData[["data"]][[1]]$group == grpId]
      idxMax = which.max(yVals)
      yMax = yVals[idxMax]
      xMax = xVals[idxMax]
      cat( "xMax: ", xMax,"; ", yMax, "\n", sep="")
      
      yMid = ggplotData[["panel"]]$ranges[[1]]$y.range[1] +
        (ggplotData[["panel"]]$ranges[[1]]$y.range[2] - ggplotData[["panel"]]$ranges[[1]]$y.range[1])*0.65
    
      ggplotList[[whichChrs]] = ggplotList[[whichChrs]] +
        geom_vline(xintercept=xMax, linetype="dotted") +
        annotate(geom="text", label=paste("(", round(xMax, 2), "; ", round(yMax, 2), ")", sep=""), x=xMax, y=yMid,  angle = 90, vjust=+1.5)
    }
  }
}

# # Get max y-limit
# yMax = 0
# for (theItem in names(ggplotList)) {
#   ggplotData = print(ggplotList[[theItem]])
#   # names(ggplotData)
#   # str(ggplotData[["data"]])
#   # attributes(ggplotData[["data"]][[1]])
#   yVals = ggplotData[["data"]][[1]]$y
#   thisYmax = max(yVals, na.rm=TRUE)
#   if (thisYmax > yMax) {
#     yMax = thisYmax
#   }
# }
# 
# # set max y-limit
# for (theItem in names(ggplotList)) {
#   ggplotList[[theItem]] = ggplotList[[theItem]] +
#     coord_cartesian(ylim = c(0, yMax))
# }

PDFfilename = file.path(allelicPropDistDir, 
                        paste(projectName, ".densityOverlaid.WTvsDelAdjustedAllelicRatios.V2.minCov=", minCov, ".pdf", sep=""))
pdf(PDFfilename, height=6, width=10)
#par(mfrow=c(1,2), oma=c(0,0,2.5,0))

grid.arrange(ggplotList[[1]], ggplotList[[2]], nrow=1, ncol=2, 
             top=textGrob(paste("Distributions of ", list.names[1], " and ", list.names[2], " allelic ratios (minCov=", minCov, ")", sep=""),
                          gp=gpar(fontsize=16)))
dev.off()






###
# Adjusted allelic proportion distributions - WT vs Del - overlaid density plots
# V3 diplosomes + chrX

allelicPropVector=list()
ggplotList = list()

#for (whichChrs in c("genome-wide", "autosomes", "chrX") ) {
for (whichChrs in c("diplosomes", "chrX") ) {
  
  for (dataSet in dataSets) {
    
    cellLine = sub("ChIPseq_Patski_", "", dataSet)
    cellLine = sub("_CTCF", "", cellLine)
    cellLine = sub("1", "", cellLine)
    
    if (whichChrs == "genome-wide") {
      allelicPropVector[[cellLine]] = as.data.frame(peakLevelAScovAndRatios[[dataSet]]$adjAltProp[ minCovIndices[[dataSet]] ])
    } else {
      allelicPropVector[[cellLine]] = as.data.frame(peakLevelAScovAndRatios[[dataSet]]$adjAltProp[ intersect(chrIndices[[dataSet]][[whichChrs]], minCovIndices[[dataSet]]) ])
    }
  }
  
  list.names <- names(allelicPropVector)
  lns <- sapply(allelicPropVector, nrow)
  dat <- as.data.frame(do.call("rbind", allelicPropVector))
  colnames(dat)[1] = "ratio"
  dat$group <- rep(list.names, lns)
  
  ggplotList[[whichChrs]] <- ggplot(dat, aes(ratio, fill=group)) + 
    geom_density(alpha=0.6, size=0.25) +
    ggtitle(paste(whichChrs, sep="")) +
    scale_fill_manual(values=sampleCols[list.names]) +
    labs(x=ifelse(whichChrs=="chrX", "(Xa/[Xa+Xi])", "(Sp/[Sp+B6])")) + 
    theme(legend.position="top") +
    theme(legend.title=element_blank())
  
  if ( whichChrs == "chrX" ) {
    
    ggplotData = print(ggplotList[[whichChrs]])
    
    for ( grpId in unique(ggplotData[["data"]][[1]]$group) ) { 
      yVals = ggplotData[["data"]][[1]]$y[ggplotData[["data"]][[1]]$group == grpId]
      xVals = ggplotData[["data"]][[1]]$x[ggplotData[["data"]][[1]]$group == grpId]
      idxMax = which.max(yVals)
      yMax = yVals[idxMax]
      xMax = xVals[idxMax]
      cat( "xMax: ", xMax,"; ", yMax, "\n", sep="")
      
      yMid = ggplotData[["panel"]]$ranges[[1]]$y.range[1] +
        (ggplotData[["panel"]]$ranges[[1]]$y.range[2] - ggplotData[["panel"]]$ranges[[1]]$y.range[1])*0.65
      
      ggplotList[[whichChrs]] = ggplotList[[whichChrs]] +
        geom_vline(xintercept=xMax, linetype="dotted") +
        annotate(geom="text", label=paste("(", round(xMax, 2), "; ", round(yMax, 2), ")", sep=""), x=xMax, y=yMid,  angle = 90, vjust=+1.5)
    }
  }
}

# # Get max y-limit
# yMax = 0
# for (theItem in names(ggplotList)) {
#   ggplotData = print(ggplotList[[theItem]])
#   # names(ggplotData)
#   # str(ggplotData[["data"]])
#   # attributes(ggplotData[["data"]][[1]])
#   yVals = ggplotData[["data"]][[1]]$y
#   thisYmax = max(yVals, na.rm=TRUE)
#   if (thisYmax > yMax) {
#     yMax = thisYmax
#   }
# }
# 
# # set max y-limit
# for (theItem in names(ggplotList)) {
#   ggplotList[[theItem]] = ggplotList[[theItem]] +
#     coord_cartesian(ylim = c(0, yMax))
# }

PDFfilename = file.path(allelicPropDistDir, 
                        paste(projectName, ".densityOverlaid.WTvsDelAdjustedAllelicRatios.V3.minCov=", minCov, ".pdf", sep=""))
pdf(PDFfilename, height=6, width=10)
#par(mfrow=c(1,2), oma=c(0,0,2.5,0))

grid.arrange(ggplotList[[1]], ggplotList[[2]], nrow=1, ncol=2, 
             top=textGrob(paste("Distributions of ", list.names[1], " and ", list.names[2], " allelic ratios (minCov=", minCov, ")", sep=""),
                          gp=gpar(fontsize=16)))
dev.off()






###
# Adjusted allelic proportion distributions - WT vs Del - overlaid density plots
# V4 autosomes, diplosomes + chrX

allelicPropVector=list()
ggplotList = list()

#for (whichChrs in c("genome-wide", "autosomes", "chrX") ) {
for (whichChrs in c("autosomes", "diplosomes", "chrX") ) {
  
  for (dataSet in dataSets) {
    
    cellLine = sub("ChIPseq_Patski_", "", dataSet)
    cellLine = sub("_CTCF", "", cellLine)
    cellLine = sub("1", "", cellLine)
    
    if (whichChrs == "genome-wide") {
      allelicPropVector[[cellLine]] = as.data.frame(peakLevelAScovAndRatios[[dataSet]]$adjAltProp[ minCovIndices[[dataSet]] ])
    } else {
      allelicPropVector[[cellLine]] = as.data.frame(peakLevelAScovAndRatios[[dataSet]]$adjAltProp[ intersect(chrIndices[[dataSet]][[whichChrs]], minCovIndices[[dataSet]]) ])
    }
  }
  
  list.names <- names(allelicPropVector)
  lns <- sapply(allelicPropVector, nrow)
  dat <- as.data.frame(do.call("rbind", allelicPropVector))
  colnames(dat)[1] = "ratio"
  dat$group <- rep(list.names, lns)
  
  ggplotList[[whichChrs]] <- ggplot(dat, aes(ratio, fill=group)) + 
    geom_density(alpha=0.2, size=0.25) +
    ggtitle(paste(whichChrs, sep="")) +
    scale_fill_manual(values=sampleCols[list.names]) +
    labs(x=ifelse(whichChrs=="chrX", "(Xa/[Xa+Xi])", "(Sp/[Sp+B6])")) + 
    theme(legend.position="top") +
    theme(legend.title=element_blank())
  
  if ( whichChrs == "chrX" ) {
    
    ggplotData = print(ggplotList[[whichChrs]])
    
    for ( grpId in unique(ggplotData[["data"]][[1]]$group) ) { 
      yVals = ggplotData[["data"]][[1]]$y[ggplotData[["data"]][[1]]$group == grpId]
      xVals = ggplotData[["data"]][[1]]$x[ggplotData[["data"]][[1]]$group == grpId]
      idxMax = which.max(yVals)
      yMax = yVals[idxMax]
      xMax = xVals[idxMax]
      cat( "xMax: ", xMax,"; ", yMax, "\n", sep="")
      
      yMid = ggplotData[["panel"]]$ranges[[1]]$y.range[1] +
        (ggplotData[["panel"]]$ranges[[1]]$y.range[2] - ggplotData[["panel"]]$ranges[[1]]$y.range[1])*0.65
      
      ggplotList[[whichChrs]] = ggplotList[[whichChrs]] +
        geom_vline(xintercept=xMax, linetype="dotted") +
        annotate(geom="text", label=paste("(", round(xMax, 2), "; ", round(yMax, 2), ")", sep=""), x=xMax, y=yMid,  angle = 90, vjust=+1.5)
    }
  }
}

# # Get max y-limit
# yMax = 0
# for (theItem in names(ggplotList)) {
#   ggplotData = print(ggplotList[[theItem]])
#   # names(ggplotData)
#   # str(ggplotData[["data"]])
#   # attributes(ggplotData[["data"]][[1]])
#   yVals = ggplotData[["data"]][[1]]$y
#   thisYmax = max(yVals, na.rm=TRUE)
#   if (thisYmax > yMax) {
#     yMax = thisYmax
#   }
# }
# 
# # set max y-limit
# for (theItem in names(ggplotList)) {
#   ggplotList[[theItem]] = ggplotList[[theItem]] +
#     coord_cartesian(ylim = c(0, yMax))
# }

PDFfilename = file.path(allelicPropDistDir, 
                        paste(projectName, ".densityOverlaid.WTvsDelAdjustedAllelicRatios.V4.minCov=", minCov, ".pdf", sep=""))
pdf(PDFfilename, height=6, width=15)
#par(mfrow=c(1,2), oma=c(0,0,2.5,0))

grid.arrange(ggplotList[[1]], ggplotList[[2]], ggplotList[[3]], nrow=1, ncol=3, 
             top=textGrob(paste("Distributions of ", list.names[1], ", and ", list.names[2], " allelic ratios (minCov=", minCov, ")", sep=""),
                          gp=gpar(fontsize=16)))
dev.off()




###
# Adjusted allelic proportion distributions - WT vs 2-4 only - overlaid density plots
# diplosomes + chrX

allelicPropVector=list()
ggplotList = list()

#for (whichChrs in c("genome-wide", "autosomes", "chrX") ) {
for (whichChrs in c("diplosomes", "chrX") ) {
  
  for (dataSet in dataSets[1:2]) {
    
    cellLine = sub("_rep1", "", dataSet)
    cellLine = sub("1", "", cellLine)
    cellLine = sub("_CTCF", "", cellLine) 
    
    if (whichChrs == "genome-wide") {
      allelicPropVector[[cellLine]] = as.data.frame(peakLevelAScovAndRatios[[dataSet]]$adjAltProp[ minCovIndices[[dataSet]] ])
    } else {
      allelicPropVector[[cellLine]] = as.data.frame(peakLevelAScovAndRatios[[dataSet]]$adjAltProp[ intersect(chrIndices[[dataSet]][[whichChrs]], minCovIndices[[dataSet]]) ])
    }
  }
  
  list.names <- names(allelicPropVector)
  lns <- sapply(allelicPropVector, nrow)
  dat <- as.data.frame(do.call("rbind", allelicPropVector))
  colnames(dat)[1] = "ratio"
  dat$group <- rep(list.names, lns)
  
  ggplotList[[whichChrs]] <- ggplot(dat, aes(ratio, fill=group)) + 
    geom_density(adjust=1, alpha=0.5, size=0.25) +
    ggtitle(paste(whichChrs, sep="")) +
    scale_fill_manual(values=sampleCols[list.names]) +
    labs(x=ifelse(whichChrs=="chrX", "(Xa/[Xa+Xi])", "(Sp/[Sp+B6])")) + 
    theme(legend.position="top") +
    theme(legend.title=element_blank())
  
  if ( whichChrs == "chrX" ) {
    
    ggplotData = print(ggplotList[[whichChrs]])
    
#     ggplotList[[whichChrs]] = ggplotList[[whichChrs]] +
#       coord_cartesian(ylim=c(0, 2))
          
          for ( grpId in unique(ggplotData[["data"]][[1]]$group) ) { 
            yVals = ggplotData[["data"]][[1]]$y[ggplotData[["data"]][[1]]$group == grpId]
            xVals = ggplotData[["data"]][[1]]$x[ggplotData[["data"]][[1]]$group == grpId]
            idxMax = which.max(yVals)
            yMax = yVals[idxMax]
            xMax = xVals[idxMax]
            
            yMid = ggplotData[["panel"]]$ranges[[1]]$y.range[1] +
              (ggplotData[["panel"]]$ranges[[1]]$y.range[2] - ggplotData[["panel"]]$ranges[[1]]$y.range[1])*0.65
            
            ggplotList[[whichChrs]] = ggplotList[[whichChrs]] +
              geom_vline(xintercept=xMax, linetype="dotted") +
              annotate(geom="text", label=paste("(", round(xMax, 2), "; ", round(yMax, 2), ")", sep=""), x=xMax, y=yMid,  angle = 90, vjust=+1.5)
          }
  }
}

# # Get max y-limit
# yMax = 0
# for (theItem in names(ggplotList)) {
#   ggplotData = print(ggplotList[[theItem]])
#   # names(ggplotData)
#   # str(ggplotData[["data"]])
#   # attributes(ggplotData[["data"]][[1]])
#   yVals = ggplotData[["data"]][[1]]$y
#   thisYmax = max(yVals, na.rm=TRUE)
#   if (thisYmax > yMax) {
#     yMax = thisYmax
#   }
# }
# 
# # set max y-limit
# for (theItem in names(ggplotList)) {
#   ggplotList[[theItem]] = ggplotList[[theItem]] +
#     coord_cartesian(ylim = c(0, yMax))
# }

PDFfilename = file.path(allelicPropDistDir, 
                        paste(projectName, ".densityOverlaid.WTvs2-4AdjustedAllelicRatios.minCov=", minCov, ".pdf", sep=""))
pdf(PDFfilename, height=6, width=10)
#par(mfrow=c(1,2), oma=c(0,0,2.5,0))

grid.arrange(ggplotList[[1]], ggplotList[[2]], nrow=1, ncol=2, 
             top=textGrob(paste("Distributions of ", list.names[1], " and ", list.names[2], " allelic ratios (minCov=", minCov, ")", sep=""),
                          gp=gpar(fontsize=16)))
dev.off()





###
# Adjusted allelic proportion distributions - all 4 samples - overlaid density plots
# diplosomes + chrX

allelicPropVector=list()
ggplotList = list()

#for (whichChrs in c("genome-wide", "autosomes", "chrX") ) {
for (whichChrs in c("diplosomes", "chrX") ) {
  
  for (dataSet in dataSets) {
    
    cellLine = sub("_rep1", "", dataSet)
    cellLine = sub("1", "", cellLine)
    cellLine = sub("_CTCF", "", cellLine) 
    
    if (whichChrs == "genome-wide") {
      allelicPropVector[[cellLine]] = as.data.frame(peakLevelAScovAndRatios[[dataSet]]$adjAltProp[ minCovIndices[[dataSet]] ])
    } else {
      allelicPropVector[[cellLine]] = as.data.frame(peakLevelAScovAndRatios[[dataSet]]$adjAltProp[ intersect(chrIndices[[dataSet]][[whichChrs]], minCovIndices[[dataSet]]) ])
    }
  }
  
  list.names <- names(allelicPropVector)
  lns <- sapply(allelicPropVector, nrow)
  dat <- as.data.frame(do.call("rbind", allelicPropVector))
  colnames(dat)[1] = "ratio"
  dat$group <- rep(list.names, lns)
  
  ggplotList[[whichChrs]] <- ggplot(dat, aes(ratio, fill=group)) + 
    geom_density(adjust=1.25, alpha=0.5) +
    ggtitle(paste(whichChrs, sep="")) +
    scale_fill_manual(values=sampleCols[list.names]) +
    labs(x=ifelse(whichChrs=="chrX", "(Xa/[Xa+Xi])", "(Sp/[Sp+B6])")) + 
    theme(legend.position="top") +
    theme(legend.title=element_blank())
  
  if ( whichChrs == "chrX" ) {
    
    ggplotData = print(ggplotList[[whichChrs]])
    
    ggplotList[[whichChrs]] = ggplotList[[whichChrs]] +
      coord_cartesian(ylim=c(0, 2.2))
    
    #       for ( grpId in unique(ggplotData[["data"]][[1]]$group) ) { 
    #         yVals = ggplotData[["data"]][[1]]$y[ggplotData[["data"]][[1]]$group == grpId]
    #         xVals = ggplotData[["data"]][[1]]$x[ggplotData[["data"]][[1]]$group == grpId]
    #         idxMax = which.max(yVals)
    #         yMax = yVals[idxMax]
    #         xMax = xVals[idxMax]
    #         
    #         yMid = ggplotData[["panel"]]$ranges[[1]]$y.range[1] +
    #           (ggplotData[["panel"]]$ranges[[1]]$y.range[2] - ggplotData[["panel"]]$ranges[[1]]$y.range[1])*0.65
    #         
    #         ggplotList[[whichChrs]] = ggplotList[[whichChrs]] +
    #           geom_vline(xintercept=xMax, linetype="dotted") +
    #           annotate(geom="text", label=paste("(", round(xMax, 2), "; ", round(yMax, 2), ")", sep=""), x=xMax, y=3.25,  angle = 90, vjust=+1.5)
    #       }
  }
}

# # Get max y-limit
# yMax = 0
# for (theItem in names(ggplotList)) {
#   ggplotData = print(ggplotList[[theItem]])
#   # names(ggplotData)
#   # str(ggplotData[["data"]])
#   # attributes(ggplotData[["data"]][[1]])
#   yVals = ggplotData[["data"]][[1]]$y
#   thisYmax = max(yVals, na.rm=TRUE)
#   if (thisYmax > yMax) {
#     yMax = thisYmax
#   }
# }
# 
# # set max y-limit
# for (theItem in names(ggplotList)) {
#   ggplotList[[theItem]] = ggplotList[[theItem]] +
#     coord_cartesian(ylim = c(0, yMax))
# }

PDFfilename = file.path(allelicPropDistDir, 
                        paste(projectName, ".densityOverlaid.4samples.AdjustedAllelicRatios.minCov=", minCov, ".pdf", sep=""))
pdf(PDFfilename, height=6, width=10)
#par(mfrow=c(1,2), oma=c(0,0,2.5,0))

grid.arrange(ggplotList[[1]], ggplotList[[2]], nrow=1, ncol=2, 
             top=textGrob(paste("Distributions of ", list.names[1], " and ", list.names[2], " allelic ratios (minCov=", minCov, ")", sep=""),
                          gp=gpar(fontsize=16)))
dev.off()







###
# Adjusted allelic proportion distributions - all 4 samples - overlaid density plots w/ lines
# diplosomes + chrX

allelicPropVector=list()
ggplotList = list()

#for (whichChrs in c("genome-wide", "autosomes", "chrX") ) {
for (whichChrs in c("diplosomes", "chrX") ) {
  
  for (dataSet in dataSets) {
    
    cellLine = sub("_rep1", "", dataSet)
    cellLine = sub("1", "", cellLine)
    cellLine = sub("_CTCF", "", cellLine) 
    
    if (whichChrs == "genome-wide") {
      allelicPropVector[[cellLine]] = as.data.frame(peakLevelAScovAndRatios[[dataSet]]$adjAltProp[ minCovIndices[[dataSet]] ])
    } else {
      allelicPropVector[[cellLine]] = as.data.frame(peakLevelAScovAndRatios[[dataSet]]$adjAltProp[ intersect(chrIndices[[dataSet]][[whichChrs]], minCovIndices[[dataSet]]) ])
    }
  }
  
  list.names <- names(allelicPropVector)
  lns <- sapply(allelicPropVector, nrow)
  dat <- as.data.frame(do.call("rbind", allelicPropVector))
  colnames(dat)[1] = "ratio"
  dat$group <- rep(list.names, lns)
  
  ggplotList[[whichChrs]] <- ggplot(dat, aes(ratio, fill=group, color=group)) + 
    geom_density(adjust=1.25, size=1.5, alpha=0.1) +
    ggtitle(paste(whichChrs, sep="")) +
    scale_fill_manual(values=sampleCols[list.names]) +
    scale_color_manual(values=sampleCols[list.names]) +
    labs(x=ifelse(whichChrs=="chrX", "(Xa/[Xa+Xi])", "(Sp/[Sp+B6])")) + 
    theme(legend.position="top") +
    theme(legend.title=element_blank())
  
  if ( whichChrs == "chrX" ) {
    
    ggplotData = print(ggplotList[[whichChrs]])
    
    ggplotList[[whichChrs]] = ggplotList[[whichChrs]] +
      coord_cartesian(ylim=c(0, 2.2))
    
    #       for ( grpId in unique(ggplotData[["data"]][[1]]$group) ) { 
    #         yVals = ggplotData[["data"]][[1]]$y[ggplotData[["data"]][[1]]$group == grpId]
    #         xVals = ggplotData[["data"]][[1]]$x[ggplotData[["data"]][[1]]$group == grpId]
    #         idxMax = which.max(yVals)
    #         yMax = yVals[idxMax]
    #         xMax = xVals[idxMax]
    #         
    #         yMid = ggplotData[["panel"]]$ranges[[1]]$y.range[1] +
    #           (ggplotData[["panel"]]$ranges[[1]]$y.range[2] - ggplotData[["panel"]]$ranges[[1]]$y.range[1])*0.65
    #         
    #         ggplotList[[whichChrs]] = ggplotList[[whichChrs]] +
    #           geom_vline(xintercept=xMax, linetype="dotted") +
    #           annotate(geom="text", label=paste("(", round(xMax, 2), "; ", round(yMax, 2), ")", sep=""), x=xMax, y=3.25,  angle = 90, vjust=+1.5)
    #       }
  }
}

# # Get max y-limit
# yMax = 0
# for (theItem in names(ggplotList)) {
#   ggplotData = print(ggplotList[[theItem]])
#   # names(ggplotData)
#   # str(ggplotData[["data"]])
#   # attributes(ggplotData[["data"]][[1]])
#   yVals = ggplotData[["data"]][[1]]$y
#   thisYmax = max(yVals, na.rm=TRUE)
#   if (thisYmax > yMax) {
#     yMax = thisYmax
#   }
# }
# 
# # set max y-limit
# for (theItem in names(ggplotList)) {
#   ggplotList[[theItem]] = ggplotList[[theItem]] +
#     coord_cartesian(ylim = c(0, yMax))
# }

PDFfilename = file.path(allelicPropDistDir, 
                        paste(projectName, ".densityOverlaid.4samples.lines.AdjustedAllelicRatios.minCov=", minCov, ".pdf", sep=""))
pdf(PDFfilename, height=6, width=10)
#par(mfrow=c(1,2), oma=c(0,0,2.5,0))

grid.arrange(ggplotList[[1]], ggplotList[[2]], nrow=1, ncol=2, 
             top=textGrob(paste("Distributions of ", list.names[1], " and ", list.names[2], " allelic ratios (minCov=", minCov, ")", sep=""),
                          gp=gpar(fontsize=16)))
dev.off()






















###
# Adjusted allelic proportion distributions - violin plots

allelicPropDistDir = file.path(ASEdir, "AdjAllelicPropViolinPlots")
if ( ! file.exists( allelicPropDistDir ) ) {
  if ( ! dir.create( allelicPropDistDir, showWarnings=FALSE, recursive=TRUE )  ) allelicPropDistDir = ASEdir
}

allelicPropVector=list()

#for (whichChrs in c("genome-wide", "autosomes", "chrX") ) {
for (whichChrs in c("autosomes", "chrX") ) {
  
  for (dataSet in dataSets) {
    
    cellLine = sub("ChIPseq_Patski_", "", dataSet)
    cellLine = sub("_CTCF", "", cellLine)
    cellLine = sub("1", "", cellLine)
    
    label2Use = paste(whichChrs, "_", cellLine, sep="")
    
    if (whichChrs == "genome-wide") {
      allelicPropVector[[label2Use]] = as.data.frame(peakLevelAScovAndRatios[[dataSet]]$adjAltProp[ minCovIndices[[dataSet]] ])
    } else {
      allelicPropVector[[label2Use]] = as.data.frame(peakLevelAScovAndRatios[[dataSet]]$adjAltProp[ intersect(chrIndices[[dataSet]][[whichChrs]], minCovIndices[[dataSet]]) ])
    }
  }
}

list.names <- names(allelicPropVector)
lns <- sapply(allelicPropVector, nrow)
dat <- as.data.frame(do.call("rbind", allelicPropVector))
colnames(dat)[1] = "ratio"
dat$group <- rep(list.names, lns)

theCols2use = rep(sampleCols, 2)
names(theCols2use) = list.names

PDFfilename = file.path(allelicPropDistDir, 
                        paste(projectName, ".violins.WTvsDelAdjustedAllelicRatios.minCov=", minCov, ".pdf", sep=""))
pdf(PDFfilename, height=6, width=10)
#par(mfrow=c(1,2), oma=c(0,0,2.5,0))

ggplot(dat, aes(x=group, y=ratio, fill=group)) + 
  #geom_boxplot() +
  #geom_violin(trim=FALSE) + 
  #geom_jitter(alpha=0.5, aes(color=theCols2use), position=position_jitter(width = .2)) +
  geom_violin(trim=FALSE, alpha=0.75, width=1) +  # scale = "count"
  scale_fill_manual(values=theCols2use) +
  geom_boxplot(width = 0.15) +
  scale_y_continuous(minor_breaks = seq(0 , 1, 0.05), breaks = seq(0, 1, 0.1)) +
  coord_cartesian(ylim=c(-0.1, 1.1)) + 
  scale_x_discrete(limits=list.names) +
  ggtitle(paste("Distributions of allelic ratio (minCov=", minCov, ")", sep="")) +
  theme(plot.title=element_text(size=16)) + 
  labs(x="") +
  labs(y="allelic ratios per peak (Sp/[Sp+B6])") + 
  theme(legend.position='none')

dev.off()



PDFfilename = file.path(allelicPropDistDir, 
                        paste(projectName, ".violins.WTvsDelAdjustedAllelicRatios.widthScaled.minCov=", minCov, ".pdf", sep=""))
pdf(PDFfilename, height=6, width=10)
#par(mfrow=c(1,2), oma=c(0,0,2.5,0))

ggplot(dat, aes(x=group, y=ratio, fill=group)) + 
  #geom_boxplot() +
  #geom_violin(trim=FALSE) + 
  #geom_jitter(alpha=0.5, aes(color=theCols2use), position=position_jitter(width = .2)) +
  geom_violin(trim=FALSE, alpha=0.75, width=1, scale="width") + 
  scale_fill_manual(values=theCols2use) +
  geom_boxplot(width = 0.15) +
  scale_y_continuous(minor_breaks = seq(0 , 1, 0.05), breaks = seq(0, 1, 0.1)) +
  coord_cartesian(ylim=c(-0.1, 1.1)) + 
  scale_x_discrete(limits=list.names) +
  ggtitle(paste("Distributions of allelic ratio (minCov=", minCov, ")", sep="")) +
  theme(plot.title=element_text(size=16)) + 
  labs(x="") +
  labs(y="allelic ratios per peak (Sp/[Sp+B6])") + 
  theme(legend.position='none')

dev.off()






###########################################################
# peak allelic ratio bedgraphs

peakAssignmentBedDir = file.path(ASEdir, "adjAllelicPropBedGraphs")
if ( ! file.exists( peakAssignmentBedDir ) ) {
  if ( ! dir.create( peakAssignmentBedDir, showWarnings=FALSE, recursive=TRUE )  ) peakAssignmentBedDir = ASEdir
}

for (dataSet in dataSets) {
  
  cellLine = sub("ChIPseq_Patski_", "", dataSet)
  cellLine = sub("_CTCF", "", cellLine)
  cellLine = sub("1", "", cellLine)
  
  theseMinCovIndices = which(peakLevelAScovAndRatios[[dataSet]]$totCov >= minCov) 
  
  bedTable = cbind( as.character(peakLevelAScovAndRatios[[dataSet]][theseMinCovIndices,1]),
                    format(peakLevelAScovAndRatios[[dataSet]][theseMinCovIndices,2], scientific=FALSE, trim=TRUE, justify="none"),
                    format(peakLevelAScovAndRatios[[dataSet]][theseMinCovIndices,3], scientific=FALSE, trim=TRUE, justify="none"), 
                    format(peakLevelAScovAndRatios[[dataSet]][theseMinCovIndices,'adjAltProp'] - 0.5, scientific=FALSE, trim=TRUE, justify="none") )
  
  bedFilename = file.path(peakAssignmentBedDir, paste(projectName, ".", cellLine, ".dScore.minCov=", minCov, ".bedgraph", sep = ""))
  
  cat("track type=bed ",
      "name=\"", paste(projectName, ".", cellLine, ".dScore.minCov=", minCov, sep = ""), "\"", 
      "description=\"", paste(projectName, ".", cellLine, ".dScore.minCov=", minCov, sep = ""), "\"", 
      "visibility=full color=200,100,0 altColor=0,100,200 priority=20 viewLimits=0:1000 autoScale=off\n", 
      file=bedFilename, sep="")
  
  write.table(bedTable, bedFilename, 
              append = TRUE, sep = "\t", 
              col.names = FALSE, row.names = FALSE, quote = FALSE)
  try(system(paste("gzip -f \"", sub("~/", "$HOME/", bedFilename), "\"", sep = "")))
  
  # chrX
  chrXminCovIndices = intersect(theseMinCovIndices, chrIndices[[dataSet]][["chrX"]])
  
  bedTable = cbind( as.character(peakLevelAScovAndRatios[[dataSet]][chrXminCovIndices,1]),
                    format(peakLevelAScovAndRatios[[dataSet]][chrXminCovIndices,2], scientific=FALSE, trim=TRUE, justify="none"),
                    format(peakLevelAScovAndRatios[[dataSet]][chrXminCovIndices,3], scientific=FALSE, trim=TRUE, justify="none"),
                    format(peakLevelAScovAndRatios[[dataSet]][chrXminCovIndices,'adjAltProp'] - 0.5, scientific=FALSE, trim=TRUE, justify="none") )
  
  bedFilename = file.path(peakAssignmentBedDir, paste(projectName, ".", cellLine, ".dScore.minCov=", minCov, ".chrX.bedgraph", sep = ""))
  
  cat("track type=bed ",
      "name=\"", paste(projectName, ".", cellLine, ".dScore.minCov=", minCov,  ".chrX", sep = ""), "\"", 
      "description=\"", paste(projectName, ".", cellLine, ".dScore.minCov=", minCov,  ".chrX", sep = ""), "\"", 
      "visibility=full color=200,100,0 altColor=0,100,200 priority=20 viewLimits=0:1000 autoScale=off\n", 
      file=bedFilename, sep="")
  
  write.table(bedTable, bedFilename, 
              append = TRUE, sep = "\t", 
              col.names = FALSE, row.names = FALSE, quote = FALSE)
  try(system(paste("gzip -f \"", sub("~/", "$HOME/", bedFilename), "\"", sep = "")))
  
}






# ###########################################################
# # Binned d-scores
# 
# binSize = 1000000
# 
# midPoints = list()
# theIntervals = list()
# dscores = list()
# smoothedDscores = list()
# meanDscores = list()
# meanBinnedDscores = list()
# smoothedBinnedDscores = list()
# smoothedBinnedDscoresMse = list()
# smoothedBinnedDscoresPse = list()
# 
# for (dataSet in dataSets) {
#   
#   cellLine = sub("ChIPseq_Patski_", "", dataSet)
#   cellLine = sub("_CTCF", "", cellLine)
#   cellLine = sub("1", "", cellLine)
#   
#   theseMinCovIndices = which(peakLevelAScovAndRatios[[dataSet]]$totCov >= minCov) 
#   # chrX
#   chrXminCovIndices = intersect(theseMinCovIndices, chrIndices[[dataSet]][["chrX"]])
# 
#   midPoints[[dataSet]] = as.numeric(as.character(peakLevelAScovAndRatios[[dataSet]][chrXminCovIndices, 2])) + 
#     (peakLevelAScovAndRatios[[dataSet]][chrXminCovIndices, 3] - peakLevelAScovAndRatios[[dataSet]][chrXminCovIndices, 2])/2
#   midPoints[[dataSet]] = midPoints[[dataSet]]/binSize
#   
#   theIntervals[[dataSet]] = cut(midPoints[[dataSet]], breaks=seq(0,(180 * (1000000/binSize)) + 1,1))
# 
#   dscores[[dataSet]] = peakLevelAScovAndRatios[[dataSet]][chrXminCovIndices, 'adjAltProp'] - 0.5  
#   lo = loess(dscores[[dataSet]]~midPoints[[dataSet]], span=0.1)
#   smoothedDscores[[dataSet]] = predict(lo)
#   plot(dscores[[dataSet]], main=dataSet)
#   lines(smoothedDscores[[dataSet]], col="red", lwd=2)
#   
#   meanDscores[[dataSet]] = mean(dscores[[dataSet]], na.rm=TRUE)
#   
#   meanBinnedDscores[[dataSet]] = tapply(dscores[[dataSet]], theIntervals[[dataSet]], mean, na.rm=TRUE)
#   lo2 = loess( meanBinnedDscores[[dataSet]]~seq(1,(180 * (1000000/binSize)) + 1,1), span=0.1 )
#   nanBoolIdx = is.na(meanBinnedDscores[[dataSet]])
#   smoothedBinnedDscores[[dataSet]] = rep(NaN, length(meanBinnedDscores[[dataSet]]) )
#   smoothedBinnedDscoresMse[[dataSet]] = rep(NaN, length(meanBinnedDscores[[dataSet]]) )
#   smoothedBinnedDscoresPse[[dataSet]] = rep(NaN, length(meanBinnedDscores[[dataSet]]) )
#   lo2predict = predict(lo2, se=TRUE)
#   smoothedBinnedDscores[[dataSet]][!nanBoolIdx] = lo2predict$fit
#   smoothedBinnedDscoresPse[[dataSet]][!nanBoolIdx] = lo2predict$fit + lo2predict$se.fit
#   smoothedBinnedDscoresMse[[dataSet]][!nanBoolIdx] = lo2predict$fit - lo2predict$se.fit
#   plot(seq(1,(180 * (1000000/binSize)) + 1,1), meanBinnedDscores[[dataSet]], type="l", main=dataSet)
#   lines(seq(1,(180 * (1000000/binSize)) + 1,1), smoothedBinnedDscores[[dataSet]], col="red", lwd=4)
#   #lines(smoothedBinnedDscoresPse[[dataSet]], col="red", lwd=2, lty="dashed")
#   #lines(smoothedBinnedDscoresMse[[dataSet]] , col="red", lwd=2, lty="dashed")
#   x = c(seq(1,(180 * (1000000/binSize)) + 1,1), rev(seq(1,(180 * (1000000/binSize)) + 1,1)))
#   y = c(smoothedBinnedDscoresPse[[dataSet]], rev(smoothedBinnedDscoresMse[[dataSet]]))
#   polygon(x[!is.na(y)], y[!is.na(y)],
#           col="lightblue", border = NA)
#   lines(seq(1,(180 * (1000000/binSize)) + 1,1), meanBinnedDscores[[dataSet]])
#   lines(seq(1,(180 * (1000000/binSize)) + 1,1), smoothedBinnedDscores[[dataSet]], col="red", lwd=4)
#   
# }



















######################################################################################################################
######################################################################################################################
# Peak assignment  tables

peakAssignmentDir = file.path(outputDir, "peakAssignment")
if ( ! file.exists( peakAssignmentDir ) ) {
  if ( ! dir.create( peakAssignmentDir, showWarnings=FALSE, recursive=TRUE )  ) peakAssignmentDir = outputDir
}


###########################################################
# peak assignment table

minCov = 5
adjPropThresh = 0.7

peakCategoryTable = c()

for (dataSet in dataSets) {
  noCovIndices = which(peakLevelAScovAndRatios[[dataSet]]$totCov < minCov) 
  theseMinCovIndices = which(peakLevelAScovAndRatios[[dataSet]]$totCov >= minCov) 
  altPeakIndices = which(peakLevelAScovAndRatios[[dataSet]]$adjAltProp >= adjPropThresh & peakLevelAScovAndRatios[[dataSet]]$totCov >= minCov)  
  refPeakIndices = which(peakLevelAScovAndRatios[[dataSet]]$adjRefProp >= adjPropThresh & peakLevelAScovAndRatios[[dataSet]]$totCov >= minCov) 
  bothPeakIndices = which(peakLevelAScovAndRatios[[dataSet]]$adjAltProp < adjPropThresh &
                            peakLevelAScovAndRatios[[dataSet]]$adjRefProp < adjPropThresh & peakLevelAScovAndRatios[[dataSet]]$totCov >= minCov) 

  peakCategoryTable = rbind( peakCategoryTable,
                             c(dataSet, 
                               # Genomes
                               length(noCovIndices),
                               length(theseMinCovIndices),
                               length(theseMinCovIndices) / (length(theseMinCovIndices) + length(noCovIndices)),
                               length(altPeakIndices),
                               length(refPeakIndices),
                               length(bothPeakIndices),
                               length(altPeakIndices) / length(theseMinCovIndices), 
                               length(refPeakIndices) / length(theseMinCovIndices),
                               length(bothPeakIndices) / length(theseMinCovIndices),
                               # Autosomes
                               length(intersect(chrIndices[[dataSet]][["autosomes"]], noCovIndices)), 
                               length(intersect(chrIndices[[dataSet]][["autosomes"]], theseMinCovIndices)),
                               length(intersect(chrIndices[[dataSet]][["autosomes"]], theseMinCovIndices)) / 
                                 (length(intersect(chrIndices[[dataSet]][["autosomes"]], theseMinCovIndices)) + length(intersect(chrIndices[[dataSet]][["autosomes"]], noCovIndices))),
                               length(intersect(chrIndices[[dataSet]][["autosomes"]], altPeakIndices)),
                               length(intersect(chrIndices[[dataSet]][["autosomes"]], refPeakIndices)),
                               length(intersect(chrIndices[[dataSet]][["autosomes"]], bothPeakIndices)),
                               length(intersect(chrIndices[[dataSet]][["autosomes"]], altPeakIndices)) / length(intersect(chrIndices[[dataSet]][["autosomes"]], theseMinCovIndices)),
                               length(intersect(chrIndices[[dataSet]][["autosomes"]], refPeakIndices)) / length(intersect(chrIndices[[dataSet]][["autosomes"]], theseMinCovIndices)),
                               length(intersect(chrIndices[[dataSet]][["autosomes"]], bothPeakIndices)) / length(intersect(chrIndices[[dataSet]][["autosomes"]], theseMinCovIndices)),
                               # chrX
                               length(intersect(chrIndices[[dataSet]][["chrX"]], noCovIndices)), 
                               length(intersect(chrIndices[[dataSet]][["chrX"]], theseMinCovIndices)),
                               length(intersect(chrIndices[[dataSet]][["chrX"]], theseMinCovIndices)) / 
                                 (length(intersect(chrIndices[[dataSet]][["chrX"]], theseMinCovIndices)) + length(intersect(chrIndices[[dataSet]][["chrX"]], noCovIndices))),
                               length(intersect(chrIndices[[dataSet]][["chrX"]], altPeakIndices)),
                               length(intersect(chrIndices[[dataSet]][["chrX"]], refPeakIndices)),
                               length(intersect(chrIndices[[dataSet]][["chrX"]], bothPeakIndices)),
                               length(intersect(chrIndices[[dataSet]][["chrX"]], altPeakIndices)) / length(intersect(chrIndices[[dataSet]][["chrX"]], theseMinCovIndices)),
                               length(intersect(chrIndices[[dataSet]][["chrX"]], refPeakIndices)) / length(intersect(chrIndices[[dataSet]][["chrX"]], theseMinCovIndices)),
                               length(intersect(chrIndices[[dataSet]][["chrX"]], bothPeakIndices)) / length(intersect(chrIndices[[dataSet]][["chrX"]], theseMinCovIndices))                                    
                             ))
}  

colnames(peakCategoryTable) = c("dataset",
                                #
                                "genomeNotCovPeaks",
                                "genomeCovPeaks",
                                "genomeCovPeakProportion",
                                "genomeAltPeaks",
                                "genomeRefPeaks",
                                "genomeBothPeaks",
                                "genomeAltPeakProportion",
                                "genomeRefPeakProportion",
                                "genomeBothPeakProportion",
                                #
                                "autosomesNotCovPeaks",
                                "autosomesCovPeaks",
                                "autosomesCovPeakProportion",
                                "autosomesAltPeaks",
                                "autosomesRefPeaks",
                                "autosomesBothPeaks",
                                "autosomesAltPeakProportion",
                                "autosomesRefPeakProportion",
                                "autosomesBothPeakProportion", 
                                #
                                "chrXnotCovPeaks",
                                "chrXcovPeaks",
                                "chrXcovPeakProportion",
                                "chrXaltPeaks",
                                "chrXrefPeaks",
                                "chrXbothPeaks",
                                "chrXaltPeakProportion",
                                "chrXrefPeakProportion",
                                "chrXbothPeakProportion"                               
)

fileName = file.path(peakAssignmentDir, paste(projectName, ".peakCategories", ".minCov=", minCov, ".adjPropThresh=", adjPropThresh, ".tsv", sep=""))
write.table(peakCategoryTable, fileName, row.names=FALSE, col.names=TRUE, sep="\t")  


############
# count data

countdat = as.data.frame(peakCategoryTable[,c("dataset",
                                         "autosomesNotCovPeaks", "autosomesBothPeaks", "autosomesAltPeaks", "autosomesRefPeaks",   
                                         "chrXnotCovPeaks", "chrXbothPeaks", "chrXaltPeaks", "chrXrefPeaks" )])
# sapply(countdat, class)
# sapply(countdat, typeof)
# sapply(countdat, mode)

countdat[, -1] <- sapply(countdat[, -1], function(x) as.numeric(as.character(x)))
countdat[, 1] <- sapply(countdat[, 1], function(x) sub("_CTCF", "", sub("ChIPseq_Patski_", "", x)))
# Must factor datasets here so correctly ordered in plot
countdat[, 1] <- sapply(countdat[, 1], function(x) as.factor(x))
countdat <- melt(countdat)


############
# prop data

propdat = as.data.frame(peakCategoryTable[,c("dataset",
                                         "autosomesCovPeakProportion", "autosomesBothPeakProportion", "autosomesAltPeakProportion", "autosomesRefPeakProportion",   
                                         "chrXcovPeakProportion", "chrXbothPeakProportion", "chrXaltPeakProportion", "chrXrefPeakProportion" )])
# sapply(propdat, class)
# sapply(propdat, typeof)
# sapply(propdat, mode)

# Get proportion not covered
propdat[,"autosomesCovPeakProportion"] = 1 - as.numeric(as.character(propdat[,"autosomesCovPeakProportion"]))
propdat[,"chrXcovPeakProportion"] = 1 - as.numeric(as.character(propdat[,"chrXcovPeakProportion"]))

propdat[, -1] <- sapply(propdat[, -1], function(x) as.numeric(as.character(x)))
propdat[, 1] <- sapply(propdat[, 1], function(x) sub("_CTCF", "", sub("ChIPseq_Patski_", "", x)))
# Must factor datasets here so correctly ordered in plot
propdat[, 1] <- sapply(propdat[, 1], function(x) as.factor(x))
propdat <- melt(propdat)


##############
# Bar plots: counts of allelic peaks @ threshold
#
# 170424
# With percentages as labels:
# http://stackoverflow.com/questions/29831535/avoid-applying-alpha-aesthetic-to-geom-text-in-ggplot2

# Sort the way you want it
datA = countdat[ grepl("autosomes", countdat$variable), ]
datA$variable = factor(datA$variable, levels=c("autosomesNotCovPeaks", "autosomesBothPeaks",
                                               "autosomesAltPeaks", "autosomesRefPeaks"))
datA$percent = round(propdat[ grepl("autosomes", propdat$variable), "value"] * 100, 1)

# Just rename chrX columns b/c only using chrX columns for legend
datX = countdat[ grepl("chrX", countdat$variable), ]
datX$variable = sub("chrX", "", datX$variable)
datX$variable = sub("notCovPeaks", "No Coverage", datX$variable)
datX$variable = sub("bothPeaks", "Ambiguous", datX$variable)
datX$variable = sub("altPeaks", "Sp", datX$variable)
datX$variable = sub("refPeaks", "B6", datX$variable)
datX$variable = factor(datX$variable, levels=c("No Coverage", "Ambiguous", "Sp", "B6"))
datX$percent= round(propdat[ grepl("chrX", propdat$variable), "value"] * 100, 1)

theCols2use = sampleCols
names(theCols2use) = unique(countdat[,1])

ggp = list()
ggp[["autosomes"]] = ggplot(datA, aes(x=dataset, y=value)) + 
  geom_bar(aes(fill=dataset, alpha=factor(variable)), stat="identity", position="dodge", color="black") +
  geom_text(aes(fill=dataset, group=factor(variable), y=value, label=percent), position=position_dodge(width=0.9), vjust=-0.5) +
  scale_x_discrete(limits=c("WT", "2-4clone", "Del1", "InvDxz4")) +
  scale_fill_manual(values=theCols2use) +
  ggtitle(paste("autosomes", sep="")) +
  labs(x="") + labs(y="# peaks") +
  theme_bw() + theme( plot.title=element_text(size=16),
                      strip.background  = element_blank(),
                      panel.grid.major = element_line(colour = "grey80"),
                      panel.border = element_blank(),
                      axis.ticks = element_blank(),
                      panel.grid.minor.x=element_blank(),
                      panel.grid.major.x=element_blank(),
                      legend.position='none')
# ggp[["autosomes"]]$labels$fill="Data set"
# ggp[["autosomes"]]$labels$alpha="Minimum coverage"
# ggp[["autosomes"]]

ggp[["chrX"]] = ggplot(datX, aes(x=dataset, y=value)) + 
  geom_bar(aes(fill=dataset, alpha=factor(variable)), stat="identity", position="dodge", color="black") +
  geom_text(aes(fill=dataset, group=factor(variable), y=value, label=percent), position=position_dodge(width=0.9), vjust=-0.5) +
  scale_x_discrete(limits=c("WT", "2-4clone", "Del1", "InvDxz4")) +
  scale_fill_manual(values=theCols2use) +
  ggtitle(paste("chrX", sep="")) +
  labs(x="") + labs(y="# peaks") +
  theme_bw() + theme( plot.title=element_text(size=16),
                      strip.background  = element_blank(),
                      panel.grid.major = element_line(colour = "grey80"),
                      panel.border = element_blank(),
                      axis.ticks = element_blank(),
                      panel.grid.minor.x=element_blank(),
                      panel.grid.major.x=element_blank()                      
  ) 
ggp[["chrX"]]$labels$fill="Data set"
ggp[["chrX"]]$labels$alpha="Peak type"
# ggp[["chrX"]]

# grid.arrange(ggp[["autosomes"]], ggp[["chrX"]], nrow=1, ncol=2, 
#              top=textGrob(paste("Counts of CTCF peak types", sep=""),
#                            gp=gpar(fontsize=16)))
# 
# grid_arrange_shared_legend(ggp[["autosomes"]], ggp[["chrX"]], ncol=2, nrow=1, position="right",
#                            title=textGrob(paste("Counts of CTCF peak types", sep="")))

mylegend<-g_legend(ggp[["chrX"]])

PDFfilename = file.path(peakAssignmentDir, 
                        paste(projectName, ".barplots.peakCategories.counts.minCov=", minCov, ".adjPropThresh=", adjPropThresh, ".pdf", sep=""))
pdf(PDFfilename, height=6, width=12)

grid.arrange(ggp[["autosomes"]] + theme(legend.position="none"), 
             ggp[["chrX"]] + theme(legend.position="none"),
             mylegend, nrow=1, ncol=3, widths=c(10, 10, 5),
             top=textGrob(paste("Counts of CTCF peak types", sep=""), gp=gpar(fontsize=16)))

dev.off()






##############
# Bar plots: proportions of allelic peaks proportions @ threshold

# Sort the way you want it
datA = propdat[ grepl("autosomes", propdat$variable), ]
datA$variable = factor(datA$variable, levels=c("autosomesCovPeakProportion", "autosomesBothPeakProportion",
                                               "autosomesAltPeakProportion", "autosomesRefPeakProportion"))
datA$value = datA$value * 100

# Just rename chrX columns b/c only using chrX columns for legend
datX = propdat[ grepl("chrX", propdat$variable), ]
datX$variable = sub("chrX", "", datX$variable)
datX$variable = sub("covPeakProportion", "No Coverage", datX$variable)
datX$variable = sub("bothPeakProportion", "Ambiguous", datX$variable)
datX$variable = sub("altPeakProportion", "Sp", datX$variable)
datX$variable = sub("refPeakProportion", "B6", datX$variable)
datX$variable = factor(datX$variable, levels=c("No Coverage", "Ambiguous", "Sp", "B6"))
datX$value = datX$value * 100

theCols2use = sampleCols
names(theCols2use) = unique(propdat[,1])

maxProp = max(propdat$value, na.rm=TRUE) * 100

ggp = list()
ggp[["autosomes"]] = ggplot(datA, aes(x=dataset, y=value)) + 
  geom_bar(aes(fill=dataset, alpha=factor(variable)), stat="identity", position="dodge", color="black") +
  scale_x_discrete(limits=c("WT", "2-4clone", "Del1", "InvDxz4")) +
  scale_fill_manual(values=theCols2use) +
  ggtitle(paste("autosomes", sep="")) +
  labs(x="") + labs(y="Percentage of peaks") +
  coord_cartesian(ylim = c(0, maxProp)) + 
  theme_bw() + theme( plot.title=element_text(size=16),
                      strip.background  = element_blank(),
                      panel.grid.major = element_line(colour = "grey80"),
                      panel.border = element_blank(),
                      axis.ticks = element_blank(),
                      panel.grid.minor.x=element_blank(),
                      panel.grid.major.x=element_blank(),
                      legend.position='none')
# ggp[["autosomes"]]$labels$fill="Data set"
# ggp[["autosomes"]]$labels$alpha="Minimum coverage"
# ggp[["autosomes"]]

ggp[["chrX"]] = ggplot(datX, aes(x=dataset, y=value)) + 
  geom_bar(aes(fill=dataset, alpha=factor(variable)), stat="identity", position="dodge", color="black") +
  scale_x_discrete(limits=c("WT", "2-4clone", "Del1", "InvDxz4")) +
  scale_fill_manual(values=theCols2use) +
  ggtitle(paste("chrX", sep="")) +
  labs(x="") + labs(y="Percentage of peaks") +
  coord_cartesian(ylim = c(0, maxProp)) + 
  theme_bw() + theme( plot.title=element_text(size=16),
                      strip.background  = element_blank(),
                      panel.grid.major = element_line(colour = "grey80"),
                      panel.border = element_blank(),
                      axis.ticks = element_blank(),
                      panel.grid.minor.x=element_blank(),
                      panel.grid.major.x=element_blank()
  ) 
ggp[["chrX"]]$labels$fill="Data set"
ggp[["chrX"]]$labels$alpha="Peak type"
# ggp[["chrX"]]

# grid.arrange(ggp[["autosomes"]], ggp[["chrX"]], nrow=1, ncol=2, 
#              top=textGrob(paste("Counts of CTCF peak types", sep=""),
#                            gp=gpar(fontsize=16)))
# 
# grid_arrange_shared_legend(ggp[["autosomes"]], ggp[["chrX"]], ncol=2, nrow=1, position="right",
#                            title=textGrob(paste("Counts of CTCF peak types", sep="")))

mylegend<-g_legend(ggp[["chrX"]])
#dev.off()

PDFfilename = file.path(peakAssignmentDir, 
                        paste(projectName, ".barplots.peakCategories.proportions.minCov=", minCov, ".adjPropThresh=", adjPropThresh, ".pdf", sep=""))
pdf(PDFfilename, height=6, width=12)

grid.arrange(ggp[["autosomes"]] + theme(legend.position="none"), 
             ggp[["chrX"]] + theme(legend.position="none"),
             mylegend, nrow=1, ncol=3, widths=c(10, 10, 5),
             top=textGrob(paste("Proportion of CTCF peak types", sep=""), gp=gpar(fontsize=16)))

dev.off()


###########################################################
# peak assignment beds

peakAssignmentBedDir = file.path(peakAssignmentDir, "BEDfiles")
if ( ! file.exists( peakAssignmentBedDir ) ) {
  if ( ! dir.create( peakAssignmentBedDir, showWarnings=FALSE, recursive=TRUE )  ) peakAssignmentBedDir = peakAssignmentDir
}

minCov = 5
adjPropThresh = 0.7

for (dataSet in dataSets) {
  
  cellLine = sub("ChIPseq_Patski_", "", dataSet)
  cellLine = sub("_CTCF", "", cellLine)
  cellLine = sub("1", "", cellLine)
  
  noCovIndices = which(peakLevelAScovAndRatios[[dataSet]]$totCov < minCov) 
  theseMinCovIndices = which(peakLevelAScovAndRatios[[dataSet]]$totCov >= minCov) 
  altPeakIndices = which(peakLevelAScovAndRatios[[dataSet]]$adjAltProp >= adjPropThresh & peakLevelAScovAndRatios[[dataSet]]$totCov >= minCov)  
  refPeakIndices = which(peakLevelAScovAndRatios[[dataSet]]$adjRefProp >= adjPropThresh & peakLevelAScovAndRatios[[dataSet]]$totCov >= minCov) 
  bothPeakIndices = which(peakLevelAScovAndRatios[[dataSet]]$adjAltProp < adjPropThresh &
                            peakLevelAScovAndRatios[[dataSet]]$adjRefProp < adjPropThresh & peakLevelAScovAndRatios[[dataSet]]$totCov >= minCov) 

  peakCatList = list(noCovIndices, theseMinCovIndices, altPeakIndices, refPeakIndices, bothPeakIndices)
  names(peakCatList) = c("noCoverage", "Covered", "Sp-specific", "B6-specific", "ambiguous")
  
  for (peakType in names(peakCatList)) {
    print(peakType)
    bedTable = cbind( as.character(peakLevelAScovAndRatios[[dataSet]][peakCatList[[peakType]],1]),
                      format(peakLevelAScovAndRatios[[dataSet]][peakCatList[[peakType]],2], scientific=FALSE, trim=TRUE, justify="none"),
                      format(peakLevelAScovAndRatios[[dataSet]][peakCatList[[peakType]],3], scientific=FALSE, trim=TRUE, justify="none") )
    
    bedFilename = file.path(peakAssignmentBedDir, paste(projectName, ".", cellLine, ".", peakType, ".minCov=", minCov, ".adjPropThresh=", adjPropThresh, ".bed", sep = ""))
  
    cat("track type=bed ",
        "name=\"", paste(projectName, ".", cellLine, ".", peakType, ".minCov=", minCov, ".adjPropThresh=", adjPropThresh, sep = ""), "\"", 
        "description=\"", paste(projectName, ".", cellLine, ".", peakType, ".minCov=", minCov, ".adjPropThresh=", adjPropThresh, sep = ""), "\"", 
        "visibility=full color=200,100,0 altColor=0,100,200 priority=20 viewLimits=0:1000 autoScale=off\n", 
        file=bedFilename, sep="")

    write.table(bedTable, bedFilename, 
                append = TRUE, sep = "\t", 
                col.names = FALSE, row.names = FALSE, quote = FALSE)
    try(system(paste("gzip -f \"", sub("~/", "$HOME/", bedFilename), "\"", sep = "")))
  }

  # chrX
  chrXnoCovIndices = intersect(noCovIndices, chrIndices[[dataSet]][["chrX"]])
  chrXminCovIndices = intersect(theseMinCovIndices, chrIndices[[dataSet]][["chrX"]])
  chrXaltPeakIndices = intersect(altPeakIndices, chrIndices[[dataSet]][["chrX"]]) 
  chrXrefPeakIndices= intersect(refPeakIndices, chrIndices[[dataSet]][["chrX"]])
  chrXbothPeakIndices = intersect(bothPeakIndices, chrIndices[[dataSet]][["chrX"]])
  
  chrXpeakCatList = list(chrXnoCovIndices, chrXminCovIndices, chrXaltPeakIndices, chrXrefPeakIndices, chrXbothPeakIndices)
  names(chrXpeakCatList) = c("noCoverage", "Covered", "Sp-specific", "B6-specific", "ambiguous")
  
  for (peakType in names(chrXpeakCatList)) {
    print(peakType)
    bedTable = cbind( as.character(peakLevelAScovAndRatios[[dataSet]][chrXpeakCatList[[peakType]],1]),
                      format(peakLevelAScovAndRatios[[dataSet]][chrXpeakCatList[[peakType]],2], scientific=FALSE, trim=TRUE, justify="none"),
                      format(peakLevelAScovAndRatios[[dataSet]][chrXpeakCatList[[peakType]],3], scientific=FALSE, trim=TRUE, justify="none") )
    
    bedFilename = file.path(peakAssignmentBedDir, paste(projectName, ".", cellLine, ".", peakType, ".minCov=", minCov, ".adjPropThresh=", adjPropThresh, ".chrX.bed", sep = ""))
    
    cat("track type=bed ",
        "name=\"", paste(projectName, ".", cellLine, ".", peakType, ".minCov=", minCov, ".adjPropThresh=", adjPropThresh, ".chrX", sep = ""), "\"", 
        "description=\"", paste(projectName, ".", cellLine, ".", peakType, ".minCov=", minCov, ".adjPropThresh=", adjPropThresh, ".chrX", sep = ""), "\"", 
        "visibility=full color=200,100,0 altColor=0,100,200 priority=20 viewLimits=0:1000 autoScale=off\n", 
        file=bedFilename, sep="")
    
    write.table(bedTable, bedFilename, 
                append = TRUE, sep = "\t", 
                col.names = FALSE, row.names = FALSE, quote = FALSE)
    try(system(paste("gzip -f \"", sub("~/", "$HOME/", bedFilename), "\"", sep = "")))
  }
  
}
