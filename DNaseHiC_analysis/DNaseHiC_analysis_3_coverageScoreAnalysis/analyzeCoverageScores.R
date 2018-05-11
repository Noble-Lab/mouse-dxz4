##################################################
##################################################
##################################################
#
# 170711
# Analyze coverage scores
#
##################################################
##################################################
##################################################

sessionInfo()


######################################################################################################################
######################################################################################################################
# Preliminaries

projectDir="<projectDirectory>"
subDir="20170610_pooledWTDelDxz4Inv_plotAndCompareContacts_vSNPs"
workDir = "coverageScoreAnalysisAroundDxz4"

outputDir = file.path(projectDir, subDir,
                      paste(workDir, sep=""))
if ( ! file.exists( outputDir ) ) {
  if ( ! dir.create( outputDir, showWarnings=FALSE, recursive=TRUE )  ) outputDir = projectDir
}
setwd(outputDir)
getwd()

#save.image() # To .RData
#load(file = ".RData")


######################################################################################################################
# Variables

whichChroms = list()

rezs = c("500000", "100000", "40000")
# rezs = c("40000")

dataDir = list()
dataDir[["500000"]] = file.path(projectDir, subDir, "output/PatskiWTvsDelvsDxz4InvFixed.500000/coverageScores.chrX/bedgraphs")
dataDir[["100000"]] = file.path(projectDir, subDir, "output/PatskiWTvsDelvsDxz4InvFixed.100000/coverageScores.chrX/bedgraphs")
dataDir[["40000"]] = file.path(projectDir, subDir, "output/PatskiWTvsDelvsDxz4InvFixed.40000/coverageScores.chrX/bedgraphs")

theSamples = c("WT", "Del", "Dxz4Inv")

theSampleCols = list()
theSampleCols[["WTXi"]] = "blue"
theSampleCols[["WTXa"]] = "blue"
theSampleCols[["DelXi"]] = "red"    
theSampleCols[["DelXa"]] = "red"    
theSampleCols[["Dxz4InvXi"]] = "black"
theSampleCols[["Dxz4InvXa"]] = "black"

theSampleLtype = list()
theSampleLtype[["WTXi"]] = "solid"
theSampleLtype[["WTXa"]] = "dashed"
theSampleLtype[["DelXi"]] = "solid"    
theSampleLtype[["DelXa"]] = "dashed"    
theSampleLtype[["Dxz4InvXi"]] = "solid"
theSampleLtype[["Dxz4InvXa"]] = "dashed"

theSampleAltCols = list()
theSampleAltCols[["WTXi"]] = rgb(t(col2rgb("blue")/255), alph=0.25)
theSampleAltCols[["WTXa"]] = rgb(t(col2rgb("blue")/255), alph=0.25)
theSampleAltCols[["DelXi"]] = rgb(t(col2rgb("red")/255), alph=0.25)    
theSampleAltCols[["DelXa"]] = rgb(t(col2rgb("red")/255), alph=0.25)   
theSampleAltCols[["Dxz4InvXi"]] = rgb(t(col2rgb("black")/255), alph=0.25)
theSampleAltCols[["Dxz4InvXa"]] = rgb(t(col2rgb("black")/255), alph=0.25)

######################################################################################################################
# Load coverage data

covScores = list()
for (rez in rezs) {
    rezInt = as.numeric(rez)
    covScores[[rez]] = list()
    for (covScoreType in c("coverageScores", "coverageZscores", "rawCoverageScores", "rescaledCoverageScores" )) {
    #for (covScoreType in c("coverageScores", "coverageZscores", "rawCoverageScores" )) {
    covScores[[rez]][[covScoreType]] = list()
        for (sample in theSamples) {
        # sample = "invDxz4"
            for (allele in c("Xi", "Xa")) {
                cat(rez, ", ", covScoreType, ", ", sample, ", ", allele,  "\n", sep="")
                covScoreFile = file.path(dataDir[[rez]], 
                                         paste(covScoreType, ".", sample, "_", allele, ".chrX.bedGraph.gz", sep=""))
            	covScores[[rez]][[covScoreType]][[paste(sample, allele, sep="")]] <- read.table(gzfile(covScoreFile), sep="\t", header=FALSE, comment.char = "#")
            	# Fill with NaNs
            	starts = covScores[[rez]][[covScoreType]][[paste(sample, allele, sep="")]][,2]
            	deltas = diff(starts) 
            	missingIdx = which(deltas > rezInt)
                nanRows=c()
            	for (i in missingIdx) {
            	    # i =256
            	    for(j in 1:((deltas[i]/rezInt)-1)){
            	        nanRows = rbind(nanRows, c("chrX", starts[i] + j*rezInt, starts[i] + (j+1)*rezInt, NaN))
            	    }
            	}
                nanDF = data.frame(nanRows[,1],  as.numeric(as.character(nanRows[,2])),  as.numeric(as.character(nanRows[,3])),  as.numeric(as.character(nanRows[,4])))
                colnames(nanDF) = c("chrs", "start", "end", "score")
                colnames(covScores[[rez]][[covScoreType]][[paste(sample, allele, sep="")]]) = c("chrs", "start", "end", "score")
                covScores[[rez]][[covScoreType]][[paste(sample, allele, sep="")]] = rbind(covScores[[rez]][[covScoreType]][[paste(sample, allele, sep="")]], nanDF)
                starts = as.numeric(as.character(covScores[[rez]][[covScoreType]][[paste(sample, allele, sep="")]]$start))
                covScores[[rez]][[covScoreType]][[paste(sample, allele, sep="")]] = covScores[[rez]][[covScoreType]][[paste(sample, allele, sep="")]][order(starts),]
            	# midpoints
            	covScores[[rez]][[covScoreType]][[paste(sample, allele, sep="")]]$midpoint = covScores[[rez]][[covScoreType]][[paste(sample, allele, sep="")]]$start + 
            	    (covScores[[rez]][[covScoreType]][[paste(sample, allele, sep="")]]$end- covScores[[rez]][[covScoreType]][[paste(sample, allele, sep="")]]$start)/2
            }
        }
    }
}


########    
# Plot scores

scoreDir = file.path(outputDir, "scoreLinePlots")
if ( ! file.exists( scoreDir ) ) {
    if ( ! dir.create( scoreDir, showWarnings=FALSE, recursive=TRUE )  ) scoreDir = outputDir
}

# Del1chrX:75637518-75764754
# Dxz4=chrX:75721083-75764760
dxz4Start = 75721083
dxz4End = 75764760
window = 10000000

dxz4Mid = round(dxz4Start + (dxz4End-dxz4Start)/2 + 1)
dxz4WindowStart = round(dxz4Mid - (window/2))
dxz4WindowEnd = round(dxz4Mid + (window/2))

# # https://stackoverflow.com/questions/34205515/finding-local-maxima-and-minima-in-r
# # https://www.rdocumentation.org/packages/data.table/versions/1.10.4/topics/shift
# # detach("package:data.table", unload=TRUE)
# require(data.table)
# minimums <- function(x, minDiff) which(x - data.table::shift(x, 1, type='lag') < -minDiff & 
#                                            x - data.table::shift(x, 1, type='lead') < -minDiff)
# minPropDiff = 0.01

for (rez in rezs) {
    # rez="500000"
    rezInt = as.numeric(rez)

    for (covScoreType in names(covScores[[rez]])) {
        
        cat(rez, ", ", covScoreType, "\n", sep="")
        
        windowBoolIdx = covScores[[rez]][[covScoreType]][["WTXi"]]$midpoint >= dxz4WindowStart & 
            covScores[[rez]][[covScoreType]][["WTXi"]]$midpoint  <= dxz4WindowEnd  
        xcoords = covScores[[rez]][[covScoreType]][["WTXi"]]$midpoint[windowBoolIdx]
        theScores = covScores[[rez]][[covScoreType]][["WTXi"]]$score[windowBoolIdx]
        
        minMax = c(Inf, -Inf)
        for (sampsamp in names(covScores[[rez]][[covScoreType]])) {
            theRange = range(covScores[[rez]][[covScoreType]][[sampsamp]]$score, na.rm=TRUE)
            if (theRange[1] < minMax[1]) { minMax[1] = theRange[1] }
            if (theRange[2] > minMax[2]) { minMax[2] = theRange[2] }        
        }
        
        pdfName = file.path(scoreDir, paste(workDir, ".", format(rez, scientific=FALSE), 
                                             ".", covScoreType, ".pdf", sep=""))  
        pdf(file=pdfName, height=10, width=5)
        plot(xcoords, covScores[[rez]][[covScoreType]][[1]]$score[windowBoolIdx], type="n", 
             xlim=c(dxz4WindowStart, dxz4WindowEnd), ylim=minMax,
             main="Coverage scores", ylab=covScoreType)
        for (sampsamp in names(covScores[[rez]][[covScoreType]])) {
            yVals = covScores[[rez]][[covScoreType]][[sampsamp]]$score[windowBoolIdx]
            lines(xcoords, yVals, col=theSampleCols[[sampsamp]], lty=theSampleLtype[[sampsamp]], lwd=2)
        }
        rect(dxz4Start, minMax[1], dxz4End, minMax[2], col='black', border=NA)
        dev.off()
        
        if (covScoreType == "rescaledCoverageScores") {
            
            # coverage threshold
            # dxz4EndBin = which(xcoords+(rezInt/2) > dxz4End & !is.na(theScores))[1]
            # lowThresh = theScores[dxz4EndBin]
            # dxz4End; xcoords[dxz4EndBin]
            lowThresh = 0.015
            
            pdfName = file.path(scoreDir, paste(workDir, ".", format(rez, scientific=FALSE), 
                                                ".", covScoreType, ".alt.pdf", sep=""))  
            pdf(file=pdfName, height=7, width=10)
            par(mar = par()$mar + c(0,0,6,0))
            plot(xcoords, covScores[[rez]][[covScoreType]][[1]]$score[windowBoolIdx], type="n", 
                 xlim=c(dxz4Mid-3000000, dxz4Mid+5000000), ylim=c(0, 0.2),
                 main="", ylab=covScoreType)
            legtext = c()
            legcols = c()
            for (sampsamp in grep("Xi", names(covScores[[rez]][[covScoreType]]), value=TRUE)) {
                yVals = covScores[[rez]][[covScoreType]][[sampsamp]]$score[windowBoolIdx]
                lines(xcoords, yVals, col=theSampleCols[[sampsamp]], lty=theSampleLtype[[sampsamp]], lwd=2)
                minIdx = which(yVals <= lowThresh)
                if (length(minIdx) > 0) {
                    rect(xcoords[head(minIdx, 1)], minMax[1],
                         xcoords[tail(minIdx, 1)], minMax[2],
                         col=theSampleAltCols[[sampsamp]], border=NA)
                    legtext = c(legtext, 
                                paste(sampsamp, ": ", xcoords[head(minIdx, 1)], "-", xcoords[tail(minIdx, 1)],
                                      " = ", (xcoords[tail(minIdx, 1)] - xcoords[head(minIdx, 1)] + rezInt)/1000, "kb", sep=""))
                    legcols = c(legcols, theSampleAltCols[[sampsamp]])
                }
            }
            if (length(legtext) > 0){
                legend(x=dxz4Mid-3000000, y=0.25, 
                       legend=legtext, fill=legcols, xpd=TRUE, bty="n")
            }
            rect(dxz4Start, minMax[1], dxz4End, minMax[2], col='black', border=NA)
            abline(h=lowThresh, lty="dashed", col="red")
            text(x=dxz4Mid-3000000, y=lowThresh, round(lowThresh,3), pos=1, col="red")
            legend(x=dxz4Mid+1000000, y=0.25, 
                   legend=grep("Xi", names(covScores[[rez]][[covScoreType]]), value=TRUE), 
                   lty=c("solid"),
                   col=unlist(theSampleCols[grep("Xi", names(covScores[[rez]][[covScoreType]]), value=TRUE)]), 
                   xpd=TRUE, bty="n")
            title(main="Xi coverage plots", line=+8)
            dev.off()
       }
    }
}

