##################################################
##################################################
##################################################
#
# 170511
# Traditionallly binned size contact distances:
#
# 171019
# Investigating scaling coefficients
#
# 171109
# Regions
#
##################################################
##################################################
##################################################

sessionInfo()

require(plotrix)


######################################################################################################################
######################################################################################################################
# Variables and Functions

projectDir="<projectDirectory>"
workDir = "decayCurves"
inputDir=file.path(projectDir, workDir, "data")

outputDir=file.path(projectDir, workDir, "decayCurves")
if ( ! file.exists( outputDir ) ) {
    if ( ! dir.create( outputDir, showWarnings=FALSE, recursive=TRUE )  ) outputDir = dataDir
}
setwd(outputDir)
getwd()

#save.image() # To .RData
#load(file = ".RData")

rezs = c("500000")
#rezs = c("500000", "100000", "40000")

theThreeSamples = c("WT", "Del", "invDxz4")
altSamplesDesc = c("WT*", "Del-Hinge", "invDxz4")
names(altSamplesDesc) = theThreeSamples
altaltSamplesDesc = c("pooledWT", "pooledDel", "InvDxz4")
names(altaltSamplesDesc) = theThreeSamples

sampleCols = c("blue", "red", "black")
names(sampleCols) = theThreeSamples

lineStyles = c("solid", "dotted")
names(lineStyles) = c("ref", "alt")

lineWidths = c("2", "2")
names(lineWidths) = c("ref", "alt")
# theLineWidth = 3

genomeDesc = list()
genomeDesc[["autos"]] = c("B6", "Sp")
names(genomeDesc[["autos"]]) = c("ref", "alt")
genomeDesc[["chrX"]] = c("Xi", "Xa")
names(genomeDesc[["chrX"]]) = c("ref", "alt")

fitRanges = list()
fitRanges[["0.1-10"]] = c(0.1, 10) * 10^6
fitRanges[["0.5-7"]] = c(0.5, 7) * 10^6
fitRanges[["1-7"]] = c(1, 7) * 10^6
fitRanges[["1-4"]] = c(1, 4) * 10^6

maxdist = 70000000

# https://stackoverflow.com/questions/8047668/transparent-equivalent-of-given-color
makeTransparent = function(..., alpha=0.5) {

        if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
    
    alpha = floor(255*alpha)  
    newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
    
    .makeTransparent = function(col, alpha) {
        rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
    }
    
    newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
    
    return(newColor)
    
}
# polygon(log10(entireRange * rezInt)[which(fitIdx)],
#         log10(diagData[[rez]][["WT"]][["alt"]][[chrom]]$meanProp[range])[which(fitIdx)], 
#         col=makeTransparent("yellow", alpha=0.25) )

fitRangeCols = list()
fitRangeCols[["0.1-10"]] = makeTransparent(gray(0.85), alpha=1)
fitRangeCols[["0.5-7"]] = makeTransparent(gray(0.70), alpha=1)
fitRangeCols[["1-7"]] = makeTransparent("lightblue", alpha=1)
fitRangeCols[["1-4"]]  = makeTransparent("yellow", alpha=1)

fitRangeNames = c("Naumova", "Lieberman", "alt1", "alt2")
names(fitRangeNames) = names(fitRangeCols)

plot2Distance = list()
plot2Distance[["500000"]] = 100 * 10^6
plot2Distance[["100000"]] = 100 * 10^6
plot2Distance[["40000"]] = 100 * 10^6

theAnnotationSize = 1.0

whichChroms = list()
whichChroms[["chrX"]] = c("chrX")
whichChroms[["individual"]] = paste( "chr", c(1:19, "X"), sep="")
whichChroms[["autoVsChrX"]] = c( "autosomes", "chrX" )
numAutos=19

reset <- function() {
    par(mfrow=c(1, 1), oma=c(1,1,1,1), new=TRUE)
    plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)
}

rangeLabel = ""

######################################################################################################################
# ChrX regions

chrXlandmarks=list()
chrXlandmarks[["Start"]] = 1
chrXlandmarks[["Firre"]] = (50563120+(abs(50563120-50635321)/2))
chrXlandmarks[["Dxz4"]] = (75764754+(abs(75637518-75764754)/2)) 
chrXlandmarks[["Xist"]] = (103460373+(abs(103460373-103483233)/2))
chrXlandmarks[["End"]] = 171031299

chrXregions = list()
chrXregions[["Superdomain1"]] = c("Start", "Dxz4")
chrXregions[["Superdomain2"]] = c("Dxz4", "End")
chrXregions[["Start-Firre"]] = c("Start", "Firre")
chrXregions[["Start-Xist"]] = c("Start", "Xist")
chrXregions[["Firre-Dxz4"]] = c("Firre", "Dxz4")
chrXregions[["Firre-Xist"]] = c("Firre", "Xist")
chrXregions[["Dxz4-Xist"]] = c("Dxz4", "Xist")
chrXregions[["Xist-End"]] = c("Xist", "End")

chrXregionBins = list()
for (rez in rezs) {
    rezInt = as.numeric(rez)
    chrXregionBins[[rez]] = list()
    for (regionLabel in names(chrXregions)) {
        print(regionLabel)
        aLOI = chrXregions[[regionLabel]][1]
        bLOI = chrXregions[[regionLabel]][2]
        chrXregionBins[[rez]][[regionLabel]] = seq(ceiling(chrXlandmarks[[aLOI]]/rezInt),
                                             ceiling(chrXlandmarks[[bLOI]]/rezInt))
    }
}

##################
# offDiagonalRegions

offDiagonalRegions = list()
offDiagonalRegions[["Dxz4-V"]] = c("Dxz4")

chrXoffDiagonalRegionBins = list()
for (rez in rezs) {
    rezInt = as.numeric(rez)
    chrXoffDiagonalRegionBins[[rez]] = list()
    for (regionLabel in names(offDiagonalRegions)) {
        print(regionLabel)
        aLOI = offDiagonalRegions[[regionLabel]][1]
        chrXoffDiagonalRegionBins[[rez]][[regionLabel]] = list()
        chrXoffDiagonalRegionBins[[rez]][[regionLabel]][["rows"]] = seq(ceiling(chrXlandmarks[["Start"]]/rezInt), ceiling(chrXlandmarks[[aLOI]]/rezInt)) 
        chrXoffDiagonalRegionBins[[rez]][[regionLabel]][["cols"]] = seq(ceiling(chrXlandmarks[[aLOI]]/rezInt), ceiling(chrXlandmarks[["End"]]/rezInt))
    }
}

whichChroms[["chrXregions"]] = c("chrX", 
                                 names(chrXregionBins[[1]]), 
                                 names(chrXoffDiagonalRegionBins[[1]]))


######################################################################################################################
# Load data

countData = list()

for (rez in rezs) {
    rezInt = as.numeric(rez)
    countData[[rez]] = list()
    
    for (sample in theThreeSamples) {
        
        altaltSample = altaltSamplesDesc[sample]
        countData[[rez]][[sample]] = list()
        
        for ( genome in c( "alt", "ref" )) {

            countData[[rez]][[sample]][[genome]] = list()
            
            for (chrom in paste( "chr", c(1:19, "X"), sep="")) {
            # chrom = "chrX"

                cat( rez, "; ", sample, "; ", genome, "; ", chrom, "\n", sep="" )

                if (chrom == "chrX") {
                    altChrom = "chr20"
                } else {
                    altChrom = chrom
                } 
                
                # WT.500000.ref.IC.chr9.txt.gz
            	binnedContactCountFile = file.path(inputDir, 
            	                                   paste(altaltSample, ".",
            	                                         format(rez, scientific=FALSE), ".", 
            	                                         genome, ".IC.", 
            	                                         altChrom, ".txt.gz", sep=""))
            	countData[[rez]][[sample]][[genome]][[chrom]] = read.table(gzfile(binnedContactCountFile), 
            	                                                           sep="\t", header=FALSE, comment.char = "#")
            }
        }
    }
}
        

######################################################################################################################
# Contacts vs Distance

diagData = list()

for (rez in rezs) {
    rezInt = as.numeric(rez)
    diagData[[rez]] = list()
    
    for (sample in theThreeSamples) {
        
        altaltSample = altaltSamplesDesc[sample]
        diagData[[rez]][[sample]] = list()
        
            
        for ( genome in c( "alt", "ref" )) {
                
            diagData[[rez]][[sample]][[genome]] = list()

            for (chrom in paste( "chr", c(1:19, "X"), sep="")) {
                # chrom = "chrX"
                
                cat( rez, "; ", sample, "; ", genome, "; ", chrom, "\n", sep="" )
                    
                if (chrom == "chrX") {
                    altChrom = "chr20"
                } else {
                    altChrom = chrom
                } 
                
                countMatrix = data.matrix(countData[[rez]][[sample]][[genome]][[chrom]])
            	nbins = nrow(countMatrix)
            	# print(nbins)
            	entireRange = seq(0, nbins-1)
            	
            	diagTallies = c()
            	diagBins = c()

            	# vectorized method 4: *** LIGHTENING ***
            	# https://stackoverflow.com/questions/27935555/get-all-diagonal-vectors-from-matrix
            	diagIds = -1 * (row(countMatrix) - col(countMatrix))
            	# diagIds[1:5,1:5]
            	diagSplits = split(countMatrix[upper.tri(diagIds, diag=TRUE)], diagIds[upper.tri(diagIds, diag=TRUE)])
            	# names(diagSplits)
            	for (i in seq(1, nbins)) {
            	    diagTallies[i] = sum(diagSplits[[as.character(i-1)]], na.rm=TRUE)
            	    diagBins[i] = nbins - i + 1
            	    # cat(i, " ", diagTallies[i], " ", diagBins[i], "\n", sep="")
            	}
            	
            	diagData[[rez]][[sample]][[genome]][[chrom]] = cbind(diagTallies, diagBins)
            	colnames(diagData[[rez]][[sample]][[genome]][[chrom]]) = c("counts", "numbins")
            	diagData[[rez]][[sample]][[genome]][[chrom]] = as.data.frame(diagData[[rez]][[sample]][[genome]][[chrom]])
            
            	# Get chrX subregions
            	if (chrom == "chrX") {
            	    for(regionName in names(chrXregionBins[[rez]])) {
            	        subCountMatrix = countMatrix[chrXregionBins[[rez]][[regionName]], 
            	                                  chrXregionBins[[rez]][[regionName]]]
            	        nbins = nrow(subCountMatrix)
            	        # cat(regionName, ": ", nbins, "\n", sep="")
            	        entireRange = seq(0, nbins-1)
            	        
            	        diagTallies = c()
            	        diagBins = c()
            	        
            	        # vectorized method 4: *** LIGHTENING ***
            	        # https://stackoverflow.com/questions/27935555/get-all-diagonal-vectors-from-matrix
            	        diagIds = -1 * (row(subCountMatrix) - col(subCountMatrix))
            	        # diagIds[1:5,1:5]
            	        diagSplits = split(subCountMatrix[upper.tri(diagIds, diag=TRUE)], diagIds[upper.tri(diagIds, diag=TRUE)])
            	        # names(diagSplits)
            	        for (i in seq(1, nbins)) {
            	            diagTallies[i] = sum(diagSplits[[as.character(i-1)]], na.rm=TRUE)
            	            diagBins[i] = length(diagSplits[[as.character(i-1)]])
            	            # cat(i, " ", diagTallies[i], " ", diagBins[i], "\n", sep="")
            	        }
            	        
            	        diagData[[rez]][[sample]][[genome]][[regionName]] = cbind(diagTallies, diagBins)
            	        colnames(diagData[[rez]][[sample]][[genome]][[regionName]]) = c("counts", "numbins")
            	        diagData[[rez]][[sample]][[genome]][[regionName]] = as.data.frame(diagData[[rez]][[sample]][[genome]][[regionName]])
            	    }
            	    
            	    for(regionName in names(chrXoffDiagonalRegionBins[[rez]])) {
            	        subCountMatrix = countMatrix[chrXoffDiagonalRegionBins[[rez]][[regionName]][["rows"]],
            	                                     chrXoffDiagonalRegionBins[[rez]][[regionName]][["cols"]]]
            	        # subCountMatrix[,1:3]
            	        cat( rez, "; ", sample, "; ", genome, "; ", chrom, ": ", sep="" )
            	        cat( round(sum(subCountMatrix, na.rm=TRUE), 1), " / ", round(sum(countMatrix, na.rm=TRUE), 1), " = ",
            	             round(sum(subCountMatrix, na.rm=TRUE) / sum(countMatrix, na.rm=TRUE),2), "\n"
            	             )
            	        nbins = max(dim(subCountMatrix))
            	        # vectorized method 4: *** LIGHTENING ***
            	        # https://stackoverflow.com/questions/27935555/get-all-diagonal-vectors-from-matrix
            	        diagIds = -1 * (row(subCountMatrix) - col(subCountMatrix))
            	        uniqDiagIds = sort(unique(as.vector(diagIds)))
            	        nbins = length(uniqDiagIds)
            	        cat(regionName, ": ", nbins, "\n", sep="")
            	        entireRange = seq(0, nbins-1)
            	        
            	        diagTallies = c()
            	        diagBins = c()

            	        # diagIds[1:5,1:5]
            	        diagSplits = split(subCountMatrix, diagIds)
            	        # names(diagSplits)
            	        i = 0 
            	        for (uniqDiagId in uniqDiagIds) {
            	            i = i + 1
            	            diagTallies[i] = sum(diagSplits[[as.character(uniqDiagId)]], na.rm=TRUE)
            	            diagBins[i] = length(diagSplits[[as.character(uniqDiagId)]])
            	            # cat(i, " ", uniqDiagId, " ", diagTallies[i], " ", diagBins[i], "\n", sep="")
            	        }
            	        
            	        diagData[[rez]][[sample]][[genome]][[regionName]] = cbind(diagTallies, diagBins)
            	        colnames(diagData[[rez]][[sample]][[genome]][[regionName]]) = c("counts", "numbins")
            	        diagData[[rez]][[sample]][[genome]][[regionName]] = as.data.frame(diagData[[rez]][[sample]][[genome]][[regionName]])
            	    }
            	    
            	}
            }
        }
    }
}


######################################################################################################################
# Decay coefficients

for (rez in rezs) {
   rezInt = as.numeric(rez)
    for (sample in theThreeSamples) {
        altaltSample = altaltSamplesDesc[sample]
        for ( genome in c( "alt", "ref" )) {
            for ( chrom in names(diagData[[rez]][[sample]][[genome]]) ) {
                # chrom = "chrX"
            
                cat( rez, "; ", sample, "; ", genome, "; ", chrom, "\n", sep="" )

                if (chrom == "chrX") {
                    altChrom = "chr20"
                } else {
                    altChrom = chrom
                } 
                    
                # Based on entire matrix
                diagCounts = diagData[[rez]][[sample]][[genome]][[chrom]]$counts
                diagBins = diagData[[rez]][[sample]][[genome]][[chrom]]$numbins
                totalCounts = sum(diagCounts, na.rm=TRUE)
                diagProp = diagCounts/totalCounts
                diagMeanCounts = diagCounts / diagBins
                totalMeanCounts = sum(diagMeanCounts, na.rm=TRUE)
                diagMeanProp = diagMeanCounts/totalMeanCounts
                diagData[[rez]][[sample]][[genome]][[chrom]]$meanProp = diagMeanProp
            }
        }
    }
}
    
########
# Check:

for ( chrom in names(diagData[[rez]][[sample]][[genome]]) ) {
    cat("\n", chrom, ":\n", sep="")
    print(diagData[[rez]][["WT"]][["alt"]][[chrom]]$meanProp[1:5])
}

    
######################################################################################################################
# Aggregate all autosomes

Achrom = "autosomes"

for (rez in rezs) {
    rezInt = as.numeric(rez)
    for (sample in theThreeSamples) {
        cat( rez, "; ", sample, "; ", genome, "; ", Achrom, "\n", sep="" )
        for (genome in c("alt", "ref")) {
            diagData[[rez]][[sample]][[genome]][[Achrom]] = diagData[[rez]][[sample]][[genome]][["chr1"]]
        }       
        for (genome in c("alt", "ref")) {
            for (chrom in paste( "chr", c(2:19), sep="")) {
                chromLen = length(diagData[[rez]][[sample]][[genome]][[chrom]][,"counts"])
                diagData[[rez]][[sample]][[genome]][[Achrom]]$numbins[1:chromLen] = diagData[[rez]][[sample]][[genome]][[Achrom]]$numbins[1:chromLen] + 
                    diagData[[rez]][[sample]][[genome]][[chrom]]$numbins
                diagData[[rez]][[sample]][[genome]][[Achrom]][,"counts"][1:chromLen] = diagData[[rez]][[sample]][[genome]][[Achrom]][,"counts"][1:chromLen] + 
                    diagData[[rez]][[sample]][[genome]][[chrom]][,"counts"]
            }
        }
        for (genome in c("alt", "ref")) {
            nbins = nrow(diagData[[rez]][[sample]][[genome]][[Achrom]])
            entireRange = seq(0, nbins-1)
            
            # Based on entire matrix
            diagData[[rez]][[sample]][[genome]][[Achrom]]$chr = Achrom
            totalAutoCounts = sum(diagData[[rez]][[sample]][[genome]][[Achrom]][,"counts"], na.rm=TRUE)
            diagData[[rez]][[sample]][[genome]][[Achrom]]$prop = diagData[[rez]][[sample]][[genome]][[Achrom]][,"counts"]/totalAutoCounts
            diagData[[rez]][[sample]][[genome]][[Achrom]]$meanCounts = diagData[[rez]][[sample]][[genome]][[Achrom]][,"counts"] / diagData[[rez]][[sample]][[genome]][[Achrom]]$numbins
            totalMeanCounts = sum(diagData[[rez]][[sample]][[genome]][[Achrom]]$meanCounts, na.rm=TRUE)
            diagData[[rez]][[sample]][[genome]][[Achrom]]$meanProp = diagData[[rez]][[sample]][[genome]][[Achrom]]$meanCounts/totalMeanCounts
        }
    }       
}




######################################################################################################################
# Correct decay rates by scale normalizing across diagonals based on autosomal ratios
#
#  Samples/Diag:  A   B   C   ...  N
#             0  x0A x0B x0C  ... x0N
#             1  x1A x1B x1C  ... x1N
#             ...
#             n  xnA xnB xnC  ... xnN 

Achrom = "autosomes"
meanPropPerBandPerSample = list()
correctionFactor = list()
for (rez in rezs) {
    rezInt = as.numeric(rez)
    meanPropPerBandPerSample[[rez]] =c()
    sampleNames = c()
    for (sample in theThreeSamples) {
        correctionFactor[[rez]][[sample]] = list()
        for (genome in c("alt", "ref")) {
            cat( rez, "; ", sample, "; ", genome, "\n", sep="" )
            meanPropPerBandPerSample[[rez]] = cbind(meanPropPerBandPerSample[[rez]], 
                                                    diagData[[rez]][[sample]][[genome]][[Achrom]]$meanProp)
            sampleNames = c(sampleNames, paste(sample, genome, sep="_"))
        }
    }
    meanPropPerBandPerSample[[rez]] = cbind(meanPropPerBandPerSample[[rez]], 
                                            rowMeans(meanPropPerBandPerSample[[rez]], na.rm=TRUE))
    colnames(meanPropPerBandPerSample[[rez]]) = c(sampleNames, "mean")
    rownames(meanPropPerBandPerSample[[rez]]) = as.character(0:(nrow(meanPropPerBandPerSample[[rez]])-1))
    
    correctionFactor[[rez]] = list()
    for (sample in theThreeSamples) {
        correctionFactor[[rez]][[sample]] = list()
        for (genome in c("alt", "ref")) {
            cat( rez, "; ", sample, "; ", genome, "\n", sep="" )
            correctionFactor[[rez]][[sample]][[genome]] = meanPropPerBandPerSample[[rez]][,"mean"] / diagData[[rez]][[sample]][[genome]][[Achrom]]$meanProp 
            correctionFactor[[rez]][[sample]][[genome]][is.infinite(correctionFactor[[rez]][[sample]][[genome]])] = NaN
        }
    }
    
}

########
# Check:

head(meanPropPerBandPerSample[[rez]])

for (rez in rezs) {
    for (sample in theThreeSamples) {
        for (genome in c("alt", "ref")) {
            cat( rez, "; ", sample, "; ", genome, "\n", sep="" )
            print(correctionFactor[[rez]][[sample]][[genome]][3:5])
        }
    }
}


######################################################################################################################
# Fits

theFits = list()
for (rez in rezs) {
    theFits[[rez]] = list()
    rezInt = as.numeric(rez)
    for (sample in theThreeSamples) {
        theFits[[rez]][[sample]] = list()
        for ( genome in c( "alt", "ref" )) {
            theFits[[rez]][[sample]][[genome]] = list()
            for ( chrom in names(diagData[[rez]][[sample]][[genome]]) ) {
                # chrom = "chrX"
                cat( rez, "; ", sample, "; ", chrom, "\n", sep="" )
                if (chrom == "chrX") {
                    altChrom = "chr20"
                } else {
                    altChrom = chrom
                }
                nbins = nrow(diagData[[rez]][[sample]][[genome]][[chrom]])
                entireRange = seq(0, nbins-1)
                # Fits
                theFits[[rez]][[sample]][[genome]][[chrom]] = list()
                for( fitRangeName in names(fitRanges)) {
                    fitRange = fitRanges[[fitRangeName]]
                    fitIdx = (entireRange * rezInt) >= fitRange[1] & (entireRange * rezInt) <= fitRange[2] 
                    y = log10(diagData[[rez]][[sample]][[genome]][[chrom]]$meanProp)
                    y[is.infinite(y)] = NaN
                    x = log10(entireRange * rezInt)
                    fitRes = lm(y[fitIdx] ~ x[fitIdx])
                    theFits[[rez]][[sample]][[genome]][[chrom]][[fitRangeName]] = fitRes
                }
            }
        }
    }
}


######################################################################################################################
# Fits corrected

theFitsCorrected = list()
for (rez in rezs) {
    theFitsCorrected[[rez]] = list()
    rezInt = as.numeric(rez)
    for (sample in theThreeSamples) {
        theFitsCorrected[[rez]][[sample]] = list()
        for ( genome in c( "alt", "ref" )) {
            theFitsCorrected[[rez]][[sample]][[genome]] = list()
            for ( chrom in names(diagData[[rez]][[sample]][[genome]]) ) {
                # chrom = "chrX"
                cat( rez, "; ", sample, "; ", chrom, "\n", sep="" )
                
                if (chrom == "chrX") {
                    altChrom = "chr20"
                } else {
                    altChrom = chrom
                }
                nbins = nrow(diagData[[rez]][[sample]][[genome]][[chrom]])
                entireRange = seq(0, nbins-1)
                # Fits
                theFitsCorrected[[rez]][[sample]][[genome]][[chrom]] = list()
                for( fitRangeName in names(fitRanges)) {
                    fitRange = fitRanges[[fitRangeName]]
                    fitIdx = (entireRange * rezInt) >= fitRange[1] & (entireRange * rezInt) <= fitRange[2] 
                    y = log10(diagData[[rez]][[sample]][[genome]][[chrom]]$meanProp * 
                                  correctionFactor[[rez]][[sample]][[genome]][entireRange+1])
                    y[is.infinite(y)] = NaN
                    x = log10(entireRange * rezInt)
                    fitRes = lm(y[fitIdx] ~ x[fitIdx])
                    theFitsCorrected[[rez]][[sample]][[genome]][[chrom]][[fitRangeName]] = fitRes
                }
            }
        }
    }
}



######################################################################################################################
# Contact Decay Log-Log Profiles: proportions (probability)

plotDir = file.path(outputDir, "loglog")
if ( ! file.exists( plotDir ) ) {
    if ( ! dir.create( plotDir, showWarnings=FALSE, recursive=TRUE )  ) plotDir = outputDir
}

for (rez in rezs) {

    rezInt = as.numeric(rez)
    
    for (whichChrom in names(whichChroms)) {
        # whichChrom = "autoVsChrX"
        
        pdfName = file.path(plotDir, paste(workDir, ".", format(rez, scientific=FALSE), 
                                             ".loglog.", whichChrom, rangeLabel, ".pdf", sep=""))  
        if (whichChrom=="autoVsChrX") {
            pdf(file=pdfName, height=8+2, width=8*2)
            par(mfcol=c(1, 2), mar=par()$mar + c(2,2,0,0), oma=c(1,1,2,1))
        } else if (whichChrom=="individual") {
            pdf(file=pdfName, height=(4*5)+2, width=4*4)
            par(mfrow=c(5, 4), mar=par()$mar + c(2,2,0,0), oma=c(1,1,2,1))
        } else if (whichChrom=="chrXregions") {
            nRows = ceiling(length(whichChroms[[whichChrom]])/3)
            pdf(file=pdfName, height=(4*nRows)+2, width=4*3)
            par(mfrow=c(nRows, 3), mar=par()$mar + c(2,2,0,0), oma=c(1,1,2,1))
        } else {
            pdf(file=pdfName, height=8+2, width=8*1)
            par(mfcol=c(1, 1), mar=par()$mar + c(2,2,0,0), oma=c(1,1,2,1))
        }
        
        for (chrom in whichChroms[[whichChrom]]) {
            # chrom="chrX"
            cat( rez, "; ", chrom, "\n", sep="" )
            
            # maxlen = length(diagData[[rez]][["WT"]][["alt"]][[chrom]]$meanProp)
            # range = c(3:maxlen)
            maxPlotBin = ceiling(plot2Distance[[rez]] / rezInt)
            entireRange = c(1:maxPlotBin)
        
            maxdistAsBin = ceiling(maxdist/rezInt)
            allMeanProps= c(diagData[[rez]][["WT"]][["alt"]][[chrom]]$meanProp[0:maxdistAsBin],
                            diagData[[rez]][["WT"]][["ref"]][[chrom]]$meanProp[0:maxdistAsBin],
                            diagData[[rez]][["Del"]][["alt"]][[chrom]]$meanProp[0:maxdistAsBin],
                            diagData[[rez]][["Del"]][["ref"]][[chrom]]$meanProp[0:maxdistAsBin],
                            diagData[[rez]][["invDxz4"]][["alt"]][[chrom]]$meanProp[0:maxdistAsBin],
                            diagData[[rez]][["invDxz4"]][["ref"]][[chrom]]$meanProp[0:maxdistAsBin])
            maxprop = max(allMeanProps, na.rm=TRUE)
            minprop = min(allMeanProps[allMeanProps>0], na.rm=TRUE)
            #maxprop = 0.05
            #minprop = 0.0005

            # Figure frame
            y = log10(diagData[[rez]][["WT"]][["alt"]][[chrom]]$meanProp[entireRange])
            y[is.infinite(y)] = NaN
            x = log10(entireRange * rezInt)
            plot(x, y, type="l", col="white", lwd=0, main=chrom, ylab="", xlab="", xaxt="n", yaxt="n", 
                 ylim=log10(c(minprop, maxprop)), xlim=c(log10(rezInt * 2), log10(maxdist)))

            # # Fit regions & 
            # for( fitRangeName in names(fitRanges)) {
            #     fitRange = fitRanges[[fitRangeName]]
            #     fitIdx = (entireRange * rezInt) >= fitRange[1] & (entireRange * rezInt) <= fitRange[2] 
            #     
            #     # x0 = log10(entireRange * rezInt)[which(fitIdx)[1]]
            #     # x1 = log10(entireRange * rezInt)[tail(which(fitIdx),1)]
            #     # y0 = -100
            #     # y1 = +100
            #     #rect(x0, y0, x1, y1, col=fitRangeCols[[fitRangeName]], border=NA)  
            #     
            #     x = log10(entireRange * rezInt)
            #     # polygon(c(x[which(fitIdx)[1]],x[which(fitIdx)[1]], x[tail(which(fitIdx),1)] ,x[tail(which(fitIdx),1)]),  
            #     #         c(-100, y[which(fitIdx)[1]], y[tail(which(fitIdx),1)], -100),
            #     #         col=fitRangeCols[[fitRangeName]], border=NA)
            #     #
            #     # y2 = m(x2-x1) + y1
            #     polygon(c(x[which(fitIdx)[1]], x[which(fitIdx)[1]], x[tail(which(fitIdx),1)], x[tail(which(fitIdx),1)]),  
            #             c(-100, -1*(x[which(fitIdx)[1]] - x[3]) + log10(maxprop), -1*(x[tail(which(fitIdx),1)] - x[3]) + log10(maxprop), -100),
            #             col=fitRangeCols[[fitRangeName]], border=NA)
            #     
            #     # text(x[which(fitIdx)[1]], -1*(x[which(fitIdx)[1]] - log10(100000)) + log10(0.2), 
            #     #      fitRangeNames[fitRangeName], pos=4)
            # }

            # plot lines
            for (sample in theThreeSamples) {
                for (genome in c("alt", "ref")) {
                    y = log10(diagData[[rez]][[sample]][[genome]][[chrom]]$meanProp[entireRange])
                    y[is.infinite(y)] = NaN
                    x = log10(entireRange * rezInt)
                    lines(x, y, col=sampleCols[sample], lty=lineStyles[genome], lwd=lineWidths[genome]) 
                }
            }
            
            # # meana = mean(theFitAs, na.rm=TRUE)
            # meana=-1
            # mitoa=-0.5
            # text(x0, maxY+0.05, paste("Observed P(s) vs. s", sep=""), 
            #      srt=atan(meana)*180/pi, pos=4, cex=theAnnotationSize, xpd=TRUE)
            # text(x0, maxY-1.45, paste("Fitted slopes (", rezInt*2/1000, "kb - ", round(fitRange[2]/10^6), "Mb)", sep=""), 
            #      srt=atan(meana)*180/pi, pos=4, cex=theAnnotationSize)
            # text(x0, minY-0.7, paste("Fractal globule decay rate (m = -1)", sep=""), 
            #      srt=atan(meana)*180/pi, pos=4, cex=theAnnotationSize)
            # segments(x0, minY-0.75, x1, meana*(x1-x0)+minY-0.75, col="gray", lty="solid", lwd=1)
            # text(x0, minY-2.45, paste("Mitotic decay rate (m = -0.5)", sep=""), 
            #      srt=atan(mitoa)*180/pi, pos=4, cex=theAnnotationSize)
            # segments(x0, minY-2.5, x1, mitoa*(x1-x0)+minY-2.5, col="gray", lty="solid", lwd=1)
            
            prettyValues = pretty(c(log10(minprop), log10(maxprop)), n=10)
            prettyValuesForLabels = signif(10^(prettyValues), 1)                                                     
            axis(2, prettyValues, labels=prettyValuesForLabels, las=2)
            title(ylab="P(s)", line=5)
            
            prettyValues = pretty(log10(entireRange * rezInt), n=10)
            prettyValuesForLabels = 10^(prettyValues)                                                      
            kbIdx = which(prettyValuesForLabels>=1000)
            MbIdx = which(prettyValuesForLabels>=1000000)
            prettyLabels = prettyValuesForLabels
            prettyLabels[kbIdx] = paste(round(prettyValuesForLabels[kbIdx]/1000), "kb", sep="")                                                      
            prettyLabels[MbIdx] = paste(round(prettyValuesForLabels[MbIdx]/1000000), "Mb", sep="")                                                      
            axis(1, prettyValues, labels=prettyLabels, las=1)
            title(xlab="distance (s)", line=5)
            
            legend("bottomleft", horiz=FALSE,  bty='n',
                   legend=c(
                       paste("WT* ", ifelse(chrom=="chrX", "Xa", "Sp"), sep=""),
                       paste("Del-hinge ", ifelse(chrom=="chrX", "Xa", "Sp"), sep=""),
                       paste("Inv-Dxz4 ", ifelse(chrom=="chrX", "Xa", "Sp"), sep=""),
                       paste("WT* ", ifelse(chrom=="chrX", "Xi", "B6"), sep=""),
                       paste("Del-hinge ", ifelse(chrom=="chrX", "Xi", "B6"), sep=""),
                       paste("Inv-Dxz4 ", ifelse(chrom=="chrX", "Xi", "B6"), sep="")
                       ),
                   col=sampleCols[theThreeSamples],
                   lty=c(rep(lineStyles["alt"], 3), rep(lineStyles["ref"], 3)),
                   lwd=c(rep(lineWidths["alt"], 3), rep(lineWidths["ref"], 3))
                   )
            
            # Fit result table
            fitResTable = c()
            fitResTableRowNames = c()
            DFcols = c()
            for (genome in c("alt", "ref")) {
                for (sample in theThreeSamples) {
                    fitResTableRowNames = c(fitResTableRowNames, 
                                            paste(altSamplesDesc[sample], 
                                                  ifelse(chrom=="chrX", genomeDesc[["chrX"]][genome], genomeDesc[["autos"]][genome]),
                                                  sep=" ") )
                    fitResTableColNames = c()
                    fitResRow = c()
                    DFcolsRow = c()
                    for( fitRangeName in names(fitRanges)) {
                        fitRes = theFits[[rez]][[sample]][[genome]][[chrom]][[fitRangeName]]
                        b = fitRes$coefficients[1]
                        a = fitRes$coefficients[2]
                        fitResRow = c(fitResRow, round(a,2))
                        fitResTableColNames = c(fitResTableColNames, fitRangeName)
                        DFcolsRow = c(DFcolsRow, fitRangeCols[[fitRangeName]])
                    }
                    fitResTable = rbind(fitResTable, fitResRow)
                    DFcols = rbind(DFcols, DFcolsRow)
                }
            }
            rownames(fitResTable) = fitResTableRowNames
            colnames(fitResTable) = fitResTableColNames
            fitResDF = as.data.frame(fitResTable)
            # fitResDFList[[chrom]] = fitResDF
            
            # add table
            # https://stackoverflow.com/questions/18663159/conditional-coloring-of-cells-in-table
            addtable2plot("topright", table=fitResDF, bg=DFcols, title="decay rates over range (Mb)",
                          bty="o", display.rownames=TRUE, display.colnames=TRUE,
                          hlines=TRUE, vlines=TRUE,
                          xjust=2, yjust=1, cex=1)
        }
        
        reset()
        title(main=paste("WT* vs Del-hinge vs Inv-Dxz4 contact decay profiles", 
                         gsub("\\.", " ", rangeLabel),
                         " (", format(rez, scientific=FALSE), ")", sep=""), line=3)
        dev.off()
            
    }
}


######################################################################################################################
# Corrected Contact Decay Log-Log Profiles: proportions (probability)

plotDir = file.path(outputDir, "loglog.corrected")
if ( ! file.exists( plotDir ) ) {
    if ( ! dir.create( plotDir, showWarnings=FALSE, recursive=TRUE )  ) plotDir = outputDir
}

for (rez in rezs) {
    
    rezInt = as.numeric(rez)
    
    for (whichChrom in names(whichChroms)) {
        # whichChrom = "autoVsChrX"
        
        pdfName = file.path(plotDir, paste(workDir, ".", format(rez, scientific=FALSE), 
                                           ".loglog.corrected.", whichChrom, rangeLabel, ".pdf", sep=""))  
        if (whichChrom=="autoVsChrX") {
            pdf(file=pdfName, height=8+2, width=8*2)
            par(mfcol=c(1, 2), mar=par()$mar + c(2,2,0,0), oma=c(1,1,2,1))
        } else if (whichChrom=="individual") {
            pdf(file=pdfName, height=(4*5)+2, width=4*4)
            par(mfrow=c(5, 4), mar=par()$mar + c(2,2,0,0), oma=c(1,1,2,1))
        } else if (whichChrom=="chrXregions") {
            nRows = ceiling(length(whichChroms[[whichChrom]])/3)
            pdf(file=pdfName, height=(4*nRows)+2, width=4*3)
            par(mfrow=c(nRows, 3), mar=par()$mar + c(2,2,0,0), oma=c(1,1,2,1))
        } else {
            pdf(file=pdfName, height=8+2, width=8*1)
            par(mfcol=c(1, 1), mar=par()$mar + c(2,2,0,0), oma=c(1,1,2,1))
        }
        
        for (chrom in whichChroms[[whichChrom]]) {
            # chrom="chrX"
            cat( rez, "; ", chrom, "\n", sep="" )
            
            # maxlen = length(diagData[[rez]][["WT"]][["alt"]][[chrom]]$meanProp)
            # range = c(3:maxlen)
            maxPlotBin = ceiling(plot2Distance[[rez]] / rezInt)
            entireRange = c(1:maxPlotBin)
            
            maxdistAsBin = ceiling(maxdist/rezInt)
            allMeanProps= c(diagData[[rez]][["WT"]][["alt"]][[chrom]]$meanProp[0:maxdistAsBin],
                            diagData[[rez]][["WT"]][["ref"]][[chrom]]$meanProp[0:maxdistAsBin],
                            diagData[[rez]][["Del"]][["alt"]][[chrom]]$meanProp[0:maxdistAsBin],
                            diagData[[rez]][["Del"]][["ref"]][[chrom]]$meanProp[0:maxdistAsBin],
                            diagData[[rez]][["invDxz4"]][["alt"]][[chrom]]$meanProp[0:maxdistAsBin],
                            diagData[[rez]][["invDxz4"]][["ref"]][[chrom]]$meanProp[0:maxdistAsBin])
            maxprop = max(allMeanProps, na.rm=TRUE)
            minprop = min(allMeanProps[allMeanProps>0], na.rm=TRUE)
            #maxprop = 0.05
            #minprop = 0.0005
            
            # Figure frame
            y = log10(diagData[[rez]][["WT"]][["alt"]][[chrom]]$meanProp[entireRange] * 
                          correctionFactor[[rez]][["WT"]][["alt"]][entireRange])
            y[is.infinite(y)] = NaN
            x = log10(entireRange * rezInt)
            plot(x, y, type="l", col="white", lwd=0, main=chrom, ylab="", xlab="", xaxt="n", yaxt="n", 
                 ylim=log10(c(minprop, maxprop)), xlim=c(log10(rezInt * 2), log10(maxdist)))
            
            # # Fit regions & 
            # for( fitRangeName in names(fitRanges)) {
            #     fitRange = fitRanges[[fitRangeName]]
            #     fitIdx = (entireRange * rezInt) >= fitRange[1] & (entireRange * rezInt) <= fitRange[2] 
            #     
            #     # x0 = log10(entireRange * rezInt)[which(fitIdx)[1]]
            #     # x1 = log10(entireRange * rezInt)[tail(which(fitIdx),1)]
            #     # y0 = -100
            #     # y1 = +100
            #     #rect(x0, y0, x1, y1, col=fitRangeCols[[fitRangeName]], border=NA)  
            #     
            #     x = log10(entireRange * rezInt)
            #     # polygon(c(x[which(fitIdx)[1]],x[which(fitIdx)[1]], x[tail(which(fitIdx),1)] ,x[tail(which(fitIdx),1)]),  
            #     #         c(-100, y[which(fitIdx)[1]], y[tail(which(fitIdx),1)], -100),
            #     #         col=fitRangeCols[[fitRangeName]], border=NA)
            #     #
            #     # y2 = m(x2-x1) + y1
            #     polygon(c(x[which(fitIdx)[1]], x[which(fitIdx)[1]], x[tail(which(fitIdx),1)], x[tail(which(fitIdx),1)]),  
            #             c(-100, -1*(x[which(fitIdx)[1]] - x[3]) + log10(maxprop), -1*(x[tail(which(fitIdx),1)] - x[3]) + log10(maxprop), -100),
            #             col=fitRangeCols[[fitRangeName]], border=NA)
            #     
            #     # text(x[which(fitIdx)[1]], -1*(x[which(fitIdx)[1]] - log10(100000)) + log10(0.2), 
            #     #      fitRangeNames[fitRangeName], pos=4)
            # }
            
            # plot lines
            for (sample in theThreeSamples) {
                for (genome in c("alt", "ref")) {
                    cat( "\n", rez, "; ", sample, "; ", genome, "\n", sep="" )
                    print(diagData[[rez]][[sample]][[genome]][[chrom]]$meanProp[3:5])
                    print(correctionFactor[[rez]][[sample]][[genome]][3:5])
                    print((diagData[[rez]][[sample]][[genome]][[chrom]]$meanProp[entireRange] * 
                               correctionFactor[[rez]][[sample]][[genome]][entireRange])[3:5])
                    print(log10(diagData[[rez]][[sample]][[genome]][[chrom]]$meanProp[entireRange] * 
                                    correctionFactor[[rez]][[sample]][[genome]][entireRange])[3:5])
                    
                    y = log10(diagData[[rez]][[sample]][[genome]][[chrom]]$meanProp[entireRange] * 
                                  correctionFactor[[rez]][[sample]][[genome]][entireRange])
                    y[is.infinite(y)] = NaN
                    x = log10(entireRange * rezInt)
                    lines(x, y, col=sampleCols[sample], lty=lineStyles[genome], lwd=lineWidths[genome]) 
                }
            }
            
            # # meana = mean(theFitAs, na.rm=TRUE)
            # meana=-1
            # mitoa=-0.5
            # text(x0, maxY+0.05, paste("Observed P(s) vs. s", sep=""), 
            #      srt=atan(meana)*180/pi, pos=4, cex=theAnnotationSize, xpd=TRUE)
            # text(x0, maxY-1.45, paste("Fitted slopes (", rezInt*2/1000, "kb - ", round(fitRange[2]/10^6), "Mb)", sep=""), 
            #      srt=atan(meana)*180/pi, pos=4, cex=theAnnotationSize)
            # text(x0, minY-0.7, paste("Fractal globule decay rate (m = -1)", sep=""), 
            #      srt=atan(meana)*180/pi, pos=4, cex=theAnnotationSize)
            # segments(x0, minY-0.75, x1, meana*(x1-x0)+minY-0.75, col="gray", lty="solid", lwd=1)
            # text(x0, minY-2.45, paste("Mitotic decay rate (m = -0.5)", sep=""), 
            #      srt=atan(mitoa)*180/pi, pos=4, cex=theAnnotationSize)
            # segments(x0, minY-2.5, x1, mitoa*(x1-x0)+minY-2.5, col="gray", lty="solid", lwd=1)
            
            prettyValues = pretty(c(log10(minprop), log10(maxprop)), n=10)
            prettyValuesForLabels = signif(10^(prettyValues), 1)                                                     
            axis(2, prettyValues, labels=prettyValuesForLabels, las=2)
            title(ylab="P(s)", line=5)
            
            prettyValues = pretty(log10(entireRange * rezInt), n=10)
            prettyValuesForLabels = 10^(prettyValues)                                                      
            kbIdx = which(prettyValuesForLabels>=1000)
            MbIdx = which(prettyValuesForLabels>=1000000)
            prettyLabels = prettyValuesForLabels
            prettyLabels[kbIdx] = paste(round(prettyValuesForLabels[kbIdx]/1000), "kb", sep="")                                                      
            prettyLabels[MbIdx] = paste(round(prettyValuesForLabels[MbIdx]/1000000), "Mb", sep="")                                                      
            axis(1, prettyValues, labels=prettyLabels, las=1)
            title(xlab="distance (s)", line=5)
            
            legend("bottomleft", horiz=FALSE,  bty='n',
                   legend=c(
                       paste("WT* ", ifelse(chrom=="chrX", "Xa", "Sp"), sep=""),
                       paste("Del-hinge ", ifelse(chrom=="chrX", "Xa", "Sp"), sep=""),
                       paste("Inv-Dxz4 ", ifelse(chrom=="chrX", "Xa", "Sp"), sep=""),
                       paste("WT* ", ifelse(chrom=="chrX", "Xi", "B6"), sep=""),
                       paste("Del-hinge ", ifelse(chrom=="chrX", "Xi", "B6"), sep=""),
                       paste("Inv-Dxz4 ", ifelse(chrom=="chrX", "Xi", "B6"), sep="")
                   ),
                   col=sampleCols[theThreeSamples],
                   lty=c(rep(lineStyles["alt"], 3), rep(lineStyles["ref"], 3)),
                   lwd=c(rep(lineWidths["alt"], 3), rep(lineWidths["ref"], 3))
            )
            
            # Fit result table
            fitResTable = c()
            fitResTableRowNames = c()
            DFcols = c()
            for (genome in c("alt", "ref")) {
                for (sample in theThreeSamples) {
                    fitResTableRowNames = c(fitResTableRowNames, 
                                            paste(altSamplesDesc[sample], 
                                                  ifelse(chrom=="chrX", genomeDesc[["chrX"]][genome], genomeDesc[["autos"]][genome]),
                                                  sep=" ") )
                    fitResTableColNames = c()
                    fitResRow = c()
                    DFcolsRow = c()
                    for( fitRangeName in names(fitRanges)) {
                        fitRes = theFits[[rez]][[sample]][[genome]][[chrom]][[fitRangeName]]
                        b = fitRes$coefficients[1]
                        a = fitRes$coefficients[2]
                        fitResRow = c(fitResRow, round(a,2))
                        fitResTableColNames = c(fitResTableColNames, fitRangeName)
                        DFcolsRow = c(DFcolsRow, fitRangeCols[[fitRangeName]])
                    }
                    fitResTable = rbind(fitResTable, fitResRow)
                    DFcols = rbind(DFcols, DFcolsRow)
                }
            }
            rownames(fitResTable) = fitResTableRowNames
            colnames(fitResTable) = fitResTableColNames
            fitResDF = as.data.frame(fitResTable)
            # fitResDFList[[chrom]] = fitResDF
            
            # add table
            # https://stackoverflow.com/questions/18663159/conditional-coloring-of-cells-in-table
            addtable2plot("topright", table=fitResDF, bg=DFcols, title="decay rates over range (Mb)",
                          bty="o", display.rownames=TRUE, display.colnames=TRUE,
                          hlines=TRUE, vlines=TRUE,
                          xjust=2, yjust=1, cex=1)
        }
        
        reset()
        title(main=paste("WT* vs Del-hinge vs Inv-Dxz4 corrected contact decay profiles", 
                         gsub("\\.", " ", rangeLabel),
                         " (", format(rez, scientific=FALSE), ")", sep=""), line=3)
        dev.off()
        
    }
}


######################################################################################################################
# Contact Decay Profiles: unlogged distance: proportions (probability)

plotDir = file.path(outputDir, "logProb")
if ( ! file.exists( plotDir ) ) {
    if ( ! dir.create( plotDir, showWarnings=FALSE, recursive=TRUE )  ) plotDir = outputDir
}

for (rez in rezs) {
    
    rezInt = as.numeric(rez)
    
    for (whichChrom in names(whichChroms)) {
        # whichChrom = "autoVsChrX"
        
        pdfName = file.path(plotDir, paste(workDir, ".", format(rez, scientific=FALSE), 
                                           ".loggedProb.", whichChrom, rangeLabel, ".pdf", sep=""))  
        if (whichChrom=="autoVsChrX") {
            pdf(file=pdfName, height=8+2, width=8*2)
            par(mfcol=c(1, 2), mar=par()$mar + c(2,2,0,0), oma=c(1,1,2,1))
        } else if (whichChrom=="individual") {
            pdf(file=pdfName, height=(4*5)+2, width=4*4)
            par(mfrow=c(5, 4), mar=par()$mar + c(2,2,0,0), oma=c(1,1,2,1))
        } else if (whichChrom=="chrXregions") {
            nRows = ceiling(length(whichChroms[[whichChrom]])/3)
            pdf(file=pdfName, height=(4*nRows)+2, width=4*3)
            par(mfrow=c(nRows, 3), mar=par()$mar + c(2,2,0,0), oma=c(1,1,2,1))
        } else {
            pdf(file=pdfName, height=8+2, width=8*1)
            par(mfcol=c(1, 1), mar=par()$mar + c(2,2,0,0), oma=c(1,1,2,1))
        }
        
        for (chrom in whichChroms[[whichChrom]]) {
            # chrom="chrX"
            cat( rez, "; ", chrom, "\n", sep="" )
            
            # maxlen = length(countData[[rez]][["WT"]][["alt"]][[chrom]]$meanProp)
            # range = c(3:maxlen)
            maxPlotBin = ceiling(plot2Distance[[rez]] / rezInt)
            entireRange = c(1:maxPlotBin)
            
            maxdistAsBin = ceiling(maxdist/rezInt)
            allMeanProps= c(diagData[[rez]][["WT"]][["alt"]][[chrom]]$meanProp[0:maxdistAsBin],
                            diagData[[rez]][["WT"]][["ref"]][[chrom]]$meanProp[0:maxdistAsBin],
                            diagData[[rez]][["Del"]][["alt"]][[chrom]]$meanProp[0:maxdistAsBin],
                            diagData[[rez]][["Del"]][["ref"]][[chrom]]$meanProp[0:maxdistAsBin],
                            diagData[[rez]][["invDxz4"]][["alt"]][[chrom]]$meanProp[0:maxdistAsBin],
                            diagData[[rez]][["invDxz4"]][["ref"]][[chrom]]$meanProp[0:maxdistAsBin])
            maxprop = max(allMeanProps, na.rm=TRUE)
            minprop = min(allMeanProps[allMeanProps>0], na.rm=TRUE)
            #maxprop = 0.05
            #minprop = 0.0005
            
            # Figure frame
            y = log10(diagData[[rez]][["WT"]][["alt"]][[chrom]]$meanProp[entireRange])
            y[is.infinite(y)] = NaN
            x = entireRange * rezInt
            plot(x, y, type="l", col="white", lwd=0, main=chrom, ylab="", xlab="", xaxt="n", yaxt="n", 
                 ylim=log10(c(minprop, maxprop)), xlim=c(0, maxdist))
            
            # plot lines
            for (sample in theThreeSamples) {
                for (genome in c("alt", "ref")) {
                    y = log10(diagData[[rez]][[sample]][[genome]][[chrom]]$meanProp[entireRange])
                    y[is.infinite(y)] = NaN
                    x = entireRange * rezInt
                    lines(x, y, col=sampleCols[sample], lty=lineStyles[genome], lwd=lineWidths[genome]) 
                }
            }   
            
            prettyValues = pretty(c(log10(minprop), log10(maxprop)), n=10)
            prettyValuesForLabels = signif(10^(prettyValues), 1)                                                     
            axis(2, prettyValues, labels=prettyValuesForLabels, las=2)
            title(ylab="P(s)", line=5)
            
            prettyValues = pretty(c(0, maxdist), n=10)
            prettyValuesForLabels = prettyValues                                                    
            kbIdx = which(prettyValuesForLabels>=1000)
            MbIdx = which(prettyValuesForLabels>=1000000)
            prettyLabels = prettyValuesForLabels
            prettyLabels[kbIdx] = paste(round(prettyValuesForLabels[kbIdx]/1000, 2), "kb", sep="")                                                      
            prettyLabels[MbIdx] = paste(round(prettyValuesForLabels[MbIdx]/1000000, 2), sep="")                                                      
            axis(1, prettyValues, labels=prettyLabels, las=1)
            title(xlab="distance (s)", line=5)
            
            legend("bottomleft", horiz=FALSE,  bty='n',
                   legend=c(
                       paste("WT* ", ifelse(chrom=="chrX", "Xa", "Sp"), sep=""),
                       paste("Del-hinge ", ifelse(chrom=="chrX", "Xa", "Sp"), sep=""),
                       paste("Inv-Dxz4 ", ifelse(chrom=="chrX", "Xa", "Sp"), sep=""),
                       paste("WT* ", ifelse(chrom=="chrX", "Xi", "B6"), sep=""),
                       paste("Del-hinge ", ifelse(chrom=="chrX", "Xi", "B6"), sep=""),
                       paste("Inv-Dxz4 ", ifelse(chrom=="chrX", "Xi", "B6"), sep="")
                   ),
                   col=sampleCols[theThreeSamples],
                   lty=c(rep(lineStyles["alt"], 3), rep(lineStyles["ref"], 3)),
                   lwd=c(rep(lineWidths["alt"], 3), rep(lineWidths["ref"], 3))
            )
            
            # Fit result table
            fitResTable = c()
            fitResTableRowNames = c()
            DFcols = c()
            for (genome in c("alt", "ref")) {
                for (sample in theThreeSamples) {
                    fitResTableRowNames = c(fitResTableRowNames, 
                                            paste(altSamplesDesc[sample], 
                                                  ifelse(chrom=="chrX", genomeDesc[["chrX"]][genome], genomeDesc[["autos"]][genome]),
                                                  sep=" ") )
                    fitResTableColNames = c()
                    fitResRow = c()
                    DFcolsRow = c()
                    for( fitRangeName in names(fitRanges)) {
                        fitRes = theFits[[rez]][[sample]][[genome]][[chrom]][[fitRangeName]]
                        b = fitRes$coefficients[1]
                        a = fitRes$coefficients[2]
                        fitResRow = c(fitResRow, round(a,2))
                        fitResTableColNames = c(fitResTableColNames, fitRangeName)
                        DFcolsRow = c(DFcolsRow, fitRangeCols[[fitRangeName]])
                    }
                    fitResTable = rbind(fitResTable, fitResRow)
                    DFcols = rbind(DFcols, DFcolsRow)
                }
            }
            rownames(fitResTable) = fitResTableRowNames
            colnames(fitResTable) = fitResTableColNames
            fitResDF = as.data.frame(fitResTable)
            # fitResDFList[[chrom]] = fitResDF
            
            # add table
            # https://stackoverflow.com/questions/18663159/conditional-coloring-of-cells-in-table
            addtable2plot("topright", table=fitResDF, bg=DFcols, title="decay rates over range (Mb)",
                          bty="o", display.rownames=TRUE, display.colnames=TRUE,
                          hlines=TRUE, vlines=TRUE,
                          xjust=2, yjust=1, cex=1)
        }
        
        reset()
        title(main=paste("WT* vs Del-hinge vs Inv-Dxz4 contact decay profiles", 
                         gsub("\\.", " ", rangeLabel),
                         " (", format(rez, scientific=FALSE), ")", sep=""), line=3)
        dev.off()
        
    }
}





######################################################################################################################
# Corrected Unlogged Contact Decay Profiles: proportions (probability)

plotDir = file.path(outputDir, "logProb_corrected")
if ( ! file.exists( plotDir ) ) {
    if ( ! dir.create( plotDir, showWarnings=FALSE, recursive=TRUE )  ) plotDir = outputDir
}

for (rez in rezs) {
    
    rezInt = as.numeric(rez)
    
    for (whichChrom in names(whichChroms)) {
        # whichChrom = "autoVsChrX"
        
        pdfName = file.path(plotDir, paste(workDir, ".", format(rez, scientific=FALSE), 
                                           ".loggedProb.corrected.", whichChrom, rangeLabel, ".pdf", sep=""))  
        if (whichChrom=="autoVsChrX") {
            pdf(file=pdfName, height=8+2, width=8*2)
            par(mfcol=c(1, 2), mar=par()$mar + c(2,2,0,0), oma=c(1,1,2,1))
        } else if (whichChrom=="individual") {
            pdf(file=pdfName, height=(4*5)+2, width=4*4)
            par(mfrow=c(5, 4), mar=par()$mar + c(2,2,0,0), oma=c(1,1,2,1))
        } else if (whichChrom=="chrXregions") {
            nRows = ceiling(length(whichChroms[[whichChrom]])/3)
            pdf(file=pdfName, height=(4*nRows)+2, width=4*3)
            par(mfrow=c(nRows, 3), mar=par()$mar + c(2,2,0,0), oma=c(1,1,2,1))
        } else {
            pdf(file=pdfName, height=8+2, width=8*1)
            par(mfcol=c(1, 1), mar=par()$mar + c(2,2,0,0), oma=c(1,1,2,1))
        }
        
        for (chrom in whichChroms[[whichChrom]]) {
            # chrom="chrX"
            cat( rez, "; ", chrom, "\n", sep="" )
            
            # maxlen = length(diagData[[rez]][["WT"]][["alt"]][[chrom]]$meanProp)
            # range = c(3:maxlen)
            maxPlotBin = ceiling(plot2Distance[[rez]] / rezInt)
            entireRange = c(1:maxPlotBin)
            
            maxdistAsBin = ceiling(maxdist/rezInt)
            allMeanProps= c(diagData[[rez]][["WT"]][["alt"]][[chrom]]$meanProp[0:maxdistAsBin],
                            diagData[[rez]][["WT"]][["ref"]][[chrom]]$meanProp[0:maxdistAsBin],
                            diagData[[rez]][["Del"]][["alt"]][[chrom]]$meanProp[0:maxdistAsBin],
                            diagData[[rez]][["Del"]][["ref"]][[chrom]]$meanProp[0:maxdistAsBin],
                            diagData[[rez]][["invDxz4"]][["alt"]][[chrom]]$meanProp[0:maxdistAsBin],
                            diagData[[rez]][["invDxz4"]][["ref"]][[chrom]]$meanProp[0:maxdistAsBin])
            maxprop = max(allMeanProps, na.rm=TRUE)
            minprop = min(allMeanProps[allMeanProps>0], na.rm=TRUE)
            #maxprop = 0.05
            #minprop = 0.0005
            
            # Figure frame
            y = log10(diagData[[rez]][["WT"]][["alt"]][[chrom]]$meanProp[entireRange]
                      * correctionFactor[[rez]][[sample]][[genome]][entireRange])
            y[is.infinite(y)] = NaN
            x = entireRange * rezInt
            plot(x, y, type="l", col="white", lwd=0, main=chrom, ylab="", xlab="", xaxt="n", yaxt="n", 
                 ylim=log10(c(minprop, maxprop)), xlim=c(0, maxdist))
            
            # plot lines
            for (sample in theThreeSamples) {
                for (genome in c("alt", "ref")) {
                    y = log10(diagData[[rez]][[sample]][[genome]][[chrom]]$meanProp[entireRange] * 
                                  correctionFactor[[rez]][[sample]][[genome]][entireRange])
                    y[is.infinite(y)] = NaN
                    x = entireRange * rezInt
                    lines(x, y, col=sampleCols[sample], lty=lineStyles[genome], lwd=lineWidths[genome]) 
                }
            }   
            
            prettyValues = pretty(c(log10(minprop), log10(maxprop)), n=10)
            prettyValuesForLabels = signif(10^(prettyValues), 1)                                                     
            axis(2, prettyValues, labels=prettyValuesForLabels, las=2)
            title(ylab="P(s)", line=5)
            
            prettyValues = pretty(c(0, maxdist), n=10)
            prettyValuesForLabels = prettyValues                                                    
            kbIdx = which(prettyValuesForLabels>=1000)
            MbIdx = which(prettyValuesForLabels>=1000000)
            prettyLabels = prettyValuesForLabels
            prettyLabels[kbIdx] = paste(round(prettyValuesForLabels[kbIdx]/1000, 2), "kb", sep="")                                                      
            prettyLabels[MbIdx] = paste(round(prettyValuesForLabels[MbIdx]/1000000, 2), sep="")                                                      
            axis(1, prettyValues, labels=prettyLabels, las=1)
            title(xlab="distance (s)", line=5)
            
            legend("bottomleft", horiz=FALSE,  bty='n',
                   legend=c(
                       paste("WT* ", ifelse(chrom=="chrX", "Xa", "Sp"), sep=""),
                       paste("Del-hinge ", ifelse(chrom=="chrX", "Xa", "Sp"), sep=""),
                       paste("Inv-Dxz4 ", ifelse(chrom=="chrX", "Xa", "Sp"), sep=""),
                       paste("WT* ", ifelse(chrom=="chrX", "Xi", "B6"), sep=""),
                       paste("Del-hinge ", ifelse(chrom=="chrX", "Xi", "B6"), sep=""),
                       paste("Inv-Dxz4 ", ifelse(chrom=="chrX", "Xi", "B6"), sep="")
                   ),
                   col=sampleCols[theThreeSamples],
                   lty=c(rep(lineStyles["alt"], 3), rep(lineStyles["ref"], 3)),
                   lwd=c(rep(lineWidths["alt"], 3), rep(lineWidths["ref"], 3))
            )
            
            # Fit result table
            fitResTable = c()
            fitResTableRowNames = c()
            DFcols = c()
            for (genome in c("alt", "ref")) {
                for (sample in theThreeSamples) {
                    fitResTableRowNames = c(fitResTableRowNames, 
                                            paste(altSamplesDesc[sample], 
                                                  ifelse(chrom=="chrX", genomeDesc[["chrX"]][genome], genomeDesc[["autos"]][genome]),
                                                  sep=" ") )
                    fitResTableColNames = c()
                    fitResRow = c()
                    DFcolsRow = c()
                    for( fitRangeName in names(fitRanges)) {
                        fitRes = theFitsCorrected[[rez]][[sample]][[genome]][[chrom]][[fitRangeName]]
                        b = fitRes$coefficients[1]
                        a = fitRes$coefficients[2]
                        fitResRow = c(fitResRow, round(a,2))
                        fitResTableColNames = c(fitResTableColNames, fitRangeName)
                        DFcolsRow = c(DFcolsRow, fitRangeCols[[fitRangeName]])
                    }
                    fitResTable = rbind(fitResTable, fitResRow)
                    DFcols = rbind(DFcols, DFcolsRow)
                }
            }
            rownames(fitResTable) = fitResTableRowNames
            colnames(fitResTable) = fitResTableColNames
            fitResDF = as.data.frame(fitResTable)
            # fitResDFList[[chrom]] = fitResDF
            
            # add table
            # https://stackoverflow.com/questions/18663159/conditional-coloring-of-cells-in-table
            addtable2plot("topright", table=fitResDF, bg=DFcols, title="decay rates over range (Mb)",
                          bty="o", display.rownames=TRUE, display.colnames=TRUE,
                          hlines=TRUE, vlines=TRUE,
                          xjust=2, yjust=1, cex=1)
        }
        
        reset()
        title(main=paste("WT* vs Del-hinge vs Inv-Dxz4 corrected contact decay profiles", 
                         gsub("\\.", " ", rangeLabel),
                         " (", format(rez, scientific=FALSE), ")", sep=""), line=3)
        dev.off()
        
    }
}


######################################################################################################################
# Alt plots: Corrected Unlogged Contact Decay Profiles: proportions (probability)


maxdist = 30000000

plotDir = file.path(outputDir, "logProb_corrected.alt")
if ( ! file.exists( plotDir ) ) {
    if ( ! dir.create( plotDir, showWarnings=FALSE, recursive=TRUE )  ) plotDir = outputDir
}

for (rez in rezs) {
    
    rezInt = as.numeric(rez)
    
    for (whichChrom in names(whichChroms)) {
        # whichChrom = "autoVsChrX"
        
        for (chrom in whichChroms[[whichChrom]]) {
            # chrom="chrX"
            cat( rez, "; ", chrom, "\n", sep="" )
            
            # maxlen = length(diagData[[rez]][["WT"]][["alt"]][[chrom]]$meanProp)
            # range = c(3:maxlen)
            maxPlotBin = ceiling(plot2Distance[[rez]] / rezInt)
            entireRange = c(1:maxPlotBin)
            
            maxdistAsBin = ceiling(maxdist/rezInt)
            # allMeanProps= c(diagData[[rez]][["WT"]][["alt"]][[chrom]]$meanProp[0:maxdistAsBin],
            #                 diagData[[rez]][["WT"]][["ref"]][[chrom]]$meanProp[0:maxdistAsBin],
            #                 diagData[[rez]][["Del"]][["alt"]][[chrom]]$meanProp[0:maxdistAsBin],
            #                 diagData[[rez]][["Del"]][["ref"]][[chrom]]$meanProp[0:maxdistAsBin],
            #                 diagData[[rez]][["invDxz4"]][["alt"]][[chrom]]$meanProp[0:maxdistAsBin],
            #                 diagData[[rez]][["invDxz4"]][["ref"]][[chrom]]$meanProp[0:maxdistAsBin])
            # maxprop = max(allMeanProps, na.rm=TRUE)
            # minprop = min(allMeanProps[allMeanProps>0], na.rm=TRUE)
            maxprop = 0.1
            minprop = 0.003
            
            # Figure frame
            y = log10(diagData[[rez]][["WT"]][["alt"]][[chrom]]$meanProp[entireRange]
                      * correctionFactor[[rez]][[sample]][[genome]][entireRange])
            y[is.infinite(y)] = NaN
            x = entireRange * rezInt
            
            
            pdfName = file.path(plotDir, paste(workDir, ".", format(rez, scientific=FALSE), 
                                               ".loggedProb.corrected.", chrom, rangeLabel, ".pdf", sep=""))  
            pdf(file=pdfName, height=4+2, width=5*1)
            par(mfcol=c(1, 1), mar=par()$mar + c(2,2,0,0))
            
            plot(x, y, type="l", col="white", lwd=0, main=chrom, ylab="", xlab="", xaxt="n", yaxt="n", 
                 ylim=log10(c(minprop, maxprop)), xlim=c(0, maxdist))
            
            # plot lines
            for (sample in theThreeSamples) {
                for (genome in c("alt", "ref")) {
                    y = log10(diagData[[rez]][[sample]][[genome]][[chrom]]$meanProp[entireRange] * 
                                  correctionFactor[[rez]][[sample]][[genome]][entireRange])
                    y[is.infinite(y)] = NaN
                    x = entireRange * rezInt
                    lines(x, y, col=sampleCols[sample], lty=lineStyles[genome], lwd=lineWidths[genome]) 
                }
            }   
            
            prettyValues = pretty(c(log10(minprop), log10(maxprop)), n=10)
            prettyValuesForLabels = signif(10^(prettyValues), 1)                                                     
            axis(2, prettyValues, labels=prettyValuesForLabels, las=2)
            title(ylab="P(s)", line=5)
            
            prettyValues = pretty(c(0, maxdist), n=10)
            prettyValuesForLabels = prettyValues                                                    
            kbIdx = which(prettyValuesForLabels>=1000)
            MbIdx = which(prettyValuesForLabels>=1000000)
            prettyLabels = prettyValuesForLabels
            prettyLabels[kbIdx] = paste(round(prettyValuesForLabels[kbIdx]/1000, 2), "kb", sep="")                                                      
            prettyLabels[MbIdx] = paste(round(prettyValuesForLabels[MbIdx]/1000000, 2), sep="")                                                      
            axis(1, prettyValues, labels=prettyLabels, las=1)
            title(xlab="distance (s)", line=5)
            
            # legend("bottomleft", horiz=FALSE,  bty='n',
            #        legend=c(
            #            paste("WT* ", ifelse(chrom=="chrX", "Xa", "Sp"), sep=""),
            #            paste("Del-hinge ", ifelse(chrom=="chrX", "Xa", "Sp"), sep=""),
            #            paste("Inv-Dxz4 ", ifelse(chrom=="chrX", "Xa", "Sp"), sep=""),
            #            paste("WT* ", ifelse(chrom=="chrX", "Xi", "B6"), sep=""),
            #            paste("Del-hinge ", ifelse(chrom=="chrX", "Xi", "B6"), sep=""),
            #            paste("Inv-Dxz4 ", ifelse(chrom=="chrX", "Xi", "B6"), sep="")
            #        ),
            #        col=sampleCols[theThreeSamples],
            #        lty=c(rep(lineStyles["alt"], 3), rep(lineStyles["ref"], 3)),
            #        lwd=c(rep(lineWidths["alt"], 3), rep(lineWidths["ref"], 3))
            # )
            
            #     # Fit result table
            #     fitResTable = c()
            #     fitResTableRowNames = c()
            #     DFcols = c()
            #     for (genome in c("alt", "ref")) {
            #         for (sample in theThreeSamples) {
            #             fitResTableRowNames = c(fitResTableRowNames, 
            #                                     paste(altSamplesDesc[sample], 
            #                                           ifelse(chrom=="chrX", genomeDesc[["chrX"]][genome], genomeDesc[["autos"]][genome]),
            #                                           sep=" ") )
            #             fitResTableColNames = c()
            #             fitResRow = c()
            #             DFcolsRow = c()
            #             for( fitRangeName in names(fitRanges)) {
            #                 fitRes = theFitsCorrected[[rez]][[sample]][[genome]][[chrom]][[fitRangeName]]
            #                 b = fitRes$coefficients[1]
            #                 a = fitRes$coefficients[2]
            #                 fitResRow = c(fitResRow, round(a,2))
            #                 fitResTableColNames = c(fitResTableColNames, fitRangeName)
            #                 DFcolsRow = c(DFcolsRow, fitRangeCols[[fitRangeName]])
            #             }
            #             fitResTable = rbind(fitResTable, fitResRow)
            #             DFcols = rbind(DFcols, DFcolsRow)
            #         }
            #     }
            #     rownames(fitResTable) = fitResTableRowNames
            #     colnames(fitResTable) = fitResTableColNames
            #     fitResDF = as.data.frame(fitResTable)
            #     # fitResDFList[[chrom]] = fitResDF
            #     
            #     # add table
            #     # https://stackoverflow.com/questions/18663159/conditional-coloring-of-cells-in-table
            #     addtable2plot("topright", table=fitResDF, bg=DFcols, title="decay rates over range (Mb)",
            #                   bty="o", display.rownames=TRUE, display.colnames=TRUE,
            #                   hlines=TRUE, vlines=TRUE,
            #                   xjust=2, yjust=1, cex=1)
            
            # reset()
            # title(main=paste("WT* vs Del-hinge vs Inv-Dxz4 corrected contact decay profiles", 
            #                  gsub("\\.", " ", rangeLabel),
            #                  " (", format(rez, scientific=FALSE), ")", sep=""), line=3)
            dev.off()
            
        }
    }
}
