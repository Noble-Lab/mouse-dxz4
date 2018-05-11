##################################################
##################################################
##################################################
#
# 170928
# Analyze insualtion scores
#
##################################################
##################################################
##################################################

sessionInfo()


######################################################################################################################
######################################################################################################################
# Preliminaries

projectDir = "<projectDirectory>"
subDir = "20170610_pooledWTDelDxz4Inv_plotAndCompareContacts_vSNPs"
workDir = "dekkerInsulationScoreAnalysis"

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

dataDir = file.path(projectDir, subDir, "PatskiWTvsDelvsDxz4InvFixed/dekkerIS/output")

whichChroms = list()

rezs = c("500000", "40000")
# rezs = c("40000")

theSamples = c("pooledWT", "pooledDel", "Dxz4Inv")

insScoreParms = list()
insScoreParms[["500000"]] = "--is3500001--nt0.01--ids2000001--ss1000001--immean"
insScoreParms[["40000"]] = "--is520001--nt0.01--ids320001--ss160001--immean"

theSampleCols = list()
theSampleCols[["pooledWTXi"]] = "blue"
theSampleCols[["pooledWTXa"]] = "blue"
theSampleCols[["pooledDelXi"]] = "red"    
theSampleCols[["pooledDelXa"]] = "red"    
theSampleCols[["Dxz4InvXi"]] = "black"
theSampleCols[["Dxz4InvXa"]] = "black"

theSampleLtype = list()
theSampleLtype[["pooledWTXi"]] = "solid"
theSampleLtype[["pooledWTXa"]] = "dashed"
theSampleLtype[["pooledDelXi"]] = "solid"    
theSampleLtype[["pooledDelXa"]] = "dashed"    
theSampleLtype[["Dxz4InvXi"]] = "solid"
theSampleLtype[["Dxz4InvXa"]] = "dashed"

theSampleAltCols = list()
theSampleAltCols[["pooledWTXi"]] = rgb(t(col2rgb("blue")/255), alph=0.25)
theSampleAltCols[["pooledWTXa"]] = rgb(t(col2rgb("blue")/255), alph=0.25)
theSampleAltCols[["pooledDelXi"]] = rgb(t(col2rgb("red")/255), alph=0.25)    
theSampleAltCols[["pooledDelXa"]] = rgb(t(col2rgb("red")/255), alph=0.25)   
theSampleAltCols[["Dxz4InvXi"]] = rgb(t(col2rgb("black")/255), alph=0.25)
theSampleAltCols[["Dxz4InvXa"]] = rgb(t(col2rgb("black")/255), alph=0.25)

######################################################################################################################
# Load insualtion score data

insulationScores = list()
for (rez in rezs) {
    rezInt = as.numeric(rez)
    insulationScores[[rez]] = list()
    for (sample in theSamples) {
    # sample = "Dxz4Inv"
        for (allele in c("Xi", "Xa")) {
            cat(rez, ", ", insScoreParms[[rez]], ", ", sample, ", ", allele,  "\n", sep="")
            covScoreFile = file.path(dataDir,  
                                     # Dxz4Inv_Xa.40000.chrX.IC
                                     paste(sample, "_", allele, ".", rez, ".chrX.IC", sep=""),
                                     # Dxz4Inv_Xa.40000.chrX.IC--is520001--nt0.01--ids320001--ss160001--immean.insulation.bedGraph.gz
                                     paste(sample, "_", allele, ".", rez, ".chrX.IC", insScoreParms[[rez]], ".insulation.bedGraph.gz", sep=""))
        	insulationScores[[rez]][[paste(sample, allele, sep="")]] <- read.table(gzfile(covScoreFile), sep="\t", header=FALSE, skip=1, comment.char = "#")
            colnames(insulationScores[[rez]][[paste(sample, allele, sep="")]]) = c("chrs", "start", "end", "score")
        }
    }
}

names(insulationScores[[rez]])

########    
# Plot scores

scorePlotDir = file.path(outputDir, "scatterPlots")
if ( ! file.exists( scorePlotDir ) ) {
    if ( ! dir.create( scorePlotDir, showWarnings=FALSE, recursive=TRUE )  ) scorePlotDir = outputDir
}

for (rez in rezs) {
    # rez="500000"
    rezInt = as.numeric(rez)
    
    for( sampleA in names(insulationScores[[rez]]) ) {
        for(sampleB in names(insulationScores[[rez]]) ) {
            if (sampleA != sampleB) {
                print(paste(sampleA, "-vs-", sampleB, sep=""))
                
                x = insulationScores[[rez]][[sampleA]]$score 
                y = insulationScores[[rez]][[sampleB]]$score
                
                pdfName = file.path(scorePlotDir, 
                                    paste(workDir, ".", 
                                          format(rez, scientific=FALSE), 
                                                    ".", insScoreParms[[rez]], ".", 
                                          sampleA, "-vs-", sampleB, ".pdf", sep=""))  
                pdf(file=pdfName, height=6, width=6)
                plot(x, y, ylim=c(-1.5,1.5), xlim=c(-1.5,1.5), 
                     col=rgb(0,0,0,0.1,maxColorValue=1), pch=16,
                     main=paste(sampleA, "-vs-", sampleB, sep=""),
                     xlab=sampleA, ylab=sampleB)
                t = length(x)
                n = sum(complete.cases(x, y))
                r = cor(x, y, method="pearson", use="pairwise.complete.obs")
                p = cor(x, y, method="spearman", use="pairwise.complete.obs")
                reg<-lm(y~x)
                abline(a=0, b=1)
                abline(reg, col="red", lwd=2, lty="dashed")
                legend("bottomright", bty='n',  
                       legend = c( paste("total bins = ", t, sep=""),
                                   paste("paired bins = ", n, sep=""),
                                   paste("m = ", round(reg$coefficients[["x"]],2), sep=""), 
                                   paste("pearson = ", round(r,2), sep=""),
                                   paste("spearman = ", round(p,2), sep="") ),
                       text.col= c("black", "black", "red", "black", "black"))
                dev.off()
            }
        }
    }
}

