##################################################
##################################################
##################################################
#
# 170512
# Calculated compartments
#
# 171027
# New data
#
##################################################
##################################################
##################################################

sessionInfo()


######################################################################################################################
######################################################################################################################
# Preliminaries

projectDir="<projectDirectory>"
workDir = "compartments"

outputDir = file.path(projectDir,
                      paste(workDir, ifelse(spearmanCorMatrix, ".spearman", ""), sep=""))
if ( ! file.exists( outputDir ) ) {
  if ( ! dir.create( outputDir, showWarnings=FALSE, recursive=TRUE )  ) outputDir = projectDir
}
setwd(outputDir)
getwd()

#save.image() # To .RData
#load(file = ".RData")


######################################################################################################################
# Variables

# spearmanCorMatrix = FALSE
spearmanCorMatrix = TRUE

dataDir = list()
dataDir[["WT"]] = file.path(projectDir, "heatmap/insituDNaseHiC-WG-patski-pooled")
dataDir[["Del"]] = file.path(projectDir, "heatmap/insituDNaseHiC-WG-DelDxz4Patski-pooled")
dataDir[["invDxz4"]] = file.path(projectDir, "heatmap/insituDNaseHiC-WG-InvDxz4Patski")

whichChroms = list()
whichChroms[["individual"]] = paste( "chr", c(1:19, "X"), sep="")

rezs = c("500000")
# rezs = c("500000", "100000", "40000")

theSamples = c("WT", "Del", "invDxz4")
nSamples = length(theSamples)
sampleCompDesc = "PooledWTvsDelVsInvDxz4" # For file names
sampleCompDesc2 = "WT* vs Del-hinge vs Inv-Dxz4" # For figure text

######################################################################################################################
# Load *** ICed *** data

contactMaps = list()
for (rez in rezs) {
    contactMaps[[rez]] = list()
    for (sample in theSamples) {
        # sample = "invDxz4"
        contactMaps[[rez]][[sample]] = list()
        for (chrom in paste( "chr", c(1:19, "X"), sep="")) {
            # chrom = "chrX"
            contactMaps[[rez]][[sample]][[chrom]] = list()
    
            if (chrom == "chrX") {
                altChrom = "chr20"
            } else {
                altChrom = chrom
            }
            
            cat(rez, ", ", sample, ", ", chrom, "\n", sep="")
            
            # 500000.ref.chr3.txt.gz
            genome = "alt"
        	contactMapFile = file.path(paste(dataDir[[sample]], ".", format(rez, scientific=FALSE), sep=""), 
        	                                 paste(format(rez, scientific=FALSE), ".", genome, ".IC.", altChrom, ".txt.gz", sep=""))
        	contactMaps[[rez]][[sample]][[chrom]][[genome]] <- read.table(gzfile(contactMapFile), sep="\t", header=FALSE, comment.char = "#")
        	
        	genome = "ref"
        	contactMapFile = file.path(paste(dataDir[[sample]], ".", format(rez, scientific=FALSE), sep=""), 
        	                                 paste(format(rez, scientific=FALSE), ".", genome, ".IC.", altChrom, ".txt.gz", sep=""))
        	contactMaps[[rez]][[sample]][[chrom]][[genome]] <- read.table(gzfile(contactMapFile), sep="\t", header=FALSE, comment.char = "#")
        }
    }
}


######################################################################################################################
# PCA

whichChrom="individual"

princompResultsList = list()
rmIdxList = list()

for (rez in rezs) {
    rezInt = as.numeric(rez)
    
    princompResultsList[[rez]] = list()
    rmIdxList[[rez]] = list()
    
    for (chrom in whichChroms[[whichChrom]]) {
        # rez = "500000"; rezInt = 500000; chrom="chrX"

        # A common set of indices
        rmIdxAll = rep(TRUE, length(contactMaps[[rez]][[1]][[chrom]][[1]]))
        for (sample in theSamples) {
            for (genome in c("alt", "ref")) {
                dat = contactMaps[[rez]][[sample]][[chrom]][[genome]]
                rmIdx = sapply(dat,function(x) length(unique(x)))==1
                rmIdxAll = rmIdxAll & rmIdx
            }
        }
        
        # PCA
        princompResults = list()
        for (sample in theSamples) {
            for (genome in c("alt", "ref")) {
                
                cat(chrom, "; ", sample, "; ", genome, "\n", sep="")

                dat = contactMaps[[rez]][[sample]][[chrom]][[genome]]

                # Distance
                gridsize = ncol(dat)
                m = matrix(0, gridsize, gridsize)
                for(i in 1:(gridsize-1)) {
                    # i = 1
                    band = (row(m)==col(m)+i)
                    m[band] = mean(dat[band], na.rm=TRUE)
                }
                m = m + t(m)
                diag(m) = mean(diag(data.matrix(dat)), na.rm=TRUE)

                # Distance normalize
                dat.m = (dat+0.000001)/(m+0.000001)

                # Remove invariant rows/cols
                # rmIdx = sapply(dat,function(x) length(unique(x)))==1
                #length(rmIdx)
                dat.rm <- dat[!rmIdxAll,!rmIdxAll] 
                dat.m.rm <- dat.m[!rmIdxAll,!rmIdxAll] 
            
                princompResults[[paste(sample, "_", genome, sep="")]] = princomp(cor(dat.m.rm, 
                                                                                     method=ifelse(spearmanCorMatrix, "spearman", "pearson")))
            }
        }
        
        princompResultsList[[rez]][[chrom]] = princompResults
        rmIdxList[[rez]][[chrom]] = rmIdxAll
    }
}


######################################################################################################################
# PCA analysis (yes, yes, analysis of the analysis)
# Plots etc...

########    
# Variance

varDir = file.path(outputDir, "varianceExplained")
if ( ! file.exists( varDir ) ) {
    if ( ! dir.create( varDir, showWarnings=FALSE, recursive=TRUE )  ) varDir = outputDir
}

nAutosomes = length(grep("chrX", whichChroms[[whichChrom]], invert=TRUE))
nPCs = 5

for (rez in rezs) {
    rezInt = as.numeric(rez)

    princompResults = princompResultsList[[rez]][[1]] 
    princompVarTotals = list()
    for (sampsamp in names(princompResults)) {
        princompVarTotals[[sampsamp]] = rep(0, nPCs)
    }
    
    for (chrom in whichChroms[[whichChrom]]) {
        # rez = "500000"; rezInt = 500000; chrom="chr1"
        
        princompResults = princompResultsList[[rez]][[chrom]] 
        rmIdxAll = rmIdxList[[rez]][[chrom]]       
        
        # Autosomal totals
        if (chrom!="chrX") {
            for (sampsamp in names(princompResults)) {
                princompVarTotals[[sampsamp]] = princompVarTotals[[sampsamp]] + 
                    ((princompResults[[sampsamp]]$sdev^2) / sum(princompResults[[sampsamp]]$sdev^2))[1:nPCs]
            }
        }
        
        maxVar=-Inf
        for (sampsamp in names(princompResults)) {
            PC1var = princompResults[[sampsamp]]$sdev[1]^2 / sum(princompResults[[sampsamp]]$sdev^2)
            if (PC1var > maxVar) { maxVar = PC1var }        
        }
        
        pdfName = file.path(varDir, paste(workDir, ".", format(rez, scientific=FALSE), 
                                             ".PCA.varianceExplained.", sampleCompDesc, ".", chrom, ".pdf", sep=""))  
        pdf(file=pdfName, height=(2*3)+2, width=(nSamples*3)+2)
        par(mfcol=c(2, nSamples), mar=par()$mar + c(1,1,0,0), oma=c(0,0,2,0))

        for (sampsamp in names(princompResults)) {
            # If Xi, then PC1 is black
            if (grepl("_ref$", sampsamp) && chrom=="chrX") {
                barplot(princompResults[[sampsamp]]$sdev[1:nPCs]^2 / sum(princompResults[[sampsamp]]$sdev^2), 
                        main=sampsamp, ylab="Proportion of variance", las=2,
                        col=c("black", "red", "blue", "green", rep("gray", 6)),
                        ylim=c(0, maxVar))
            } else {
                barplot(princompResults[[sampsamp]]$sdev[1:nPCs]^2 / sum(princompResults[[sampsamp]]$sdev^2),  
                        main=sampsamp, ylab="Proportion of variance", las=2,
                        col=c("red", "blue", "green", rep("gray", 7)),
                        ylim=c(0, maxVar))
            }
        }
        
        reset <- function() {
            par(mfrow=c(1, 1), oma=c(1,1,1,1), new=TRUE)
            plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)
        }
        reset()
        title(main=paste(chrom, " PCA variance explained (", sampleCompDesc2, "; ", 
                         format(rez, scientific=FALSE), ")", sep=""), line=3 )
        dev.off()
    }

    princompVarMean = list()
    for (sampsamp in names(princompResults)) {
        princompVarMean[[sampsamp]] = princompVarTotals[[sampsamp]] / nAutosomes
    }
    
    maxVar=-Inf
    for (sampsamp in names(princompVarMean)) {
        PC1var = princompVarMean[[sampsamp]][1]
        if (PC1var > maxVar) { maxVar = PC1var }        
    }
    
    pdfName = file.path(varDir, paste(workDir, ".", format(rez, scientific=FALSE), 
                                      ".PCA.varianceExplained.", sampleCompDesc, ".autosomalAverages.pdf", sep=""))  
    pdf(file=pdfName, height=(2*3)+2, width=(nSamples*4)+2)
    par(mfcol=c(2, nSamples), mar=par()$mar + c(1,1,0,0), oma=c(0,0,2,0))
    for (sampsamp in names(princompResults)) {
        barplot(princompVarMean[[sampsamp]], 
                main=sampsamp, ylab="Proportion of variance", las=2,
                col=c("red", "blue", "green", rep("gray", 7)),
                ylim=c(0, maxVar))
    }
    reset <- function() {
        par(mfrow=c(1, 1), oma=c(1,1,1,1), new=TRUE)
        plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)
    }
    reset()
    title(main=paste("average autosomal PCA variance explained (", sampleCompDesc2, "; ", 
                     format(rez, scientific=FALSE), ")", sep=""), line=3 )
    dev.off()
}


########    
# Scores

scoreDir = file.path(outputDir, "scoreLinePlots")
if ( ! file.exists( scoreDir ) ) {
    if ( ! dir.create( scoreDir, showWarnings=FALSE, recursive=TRUE )  ) scoreDir = outputDir
}

signChanger = list()
    
for (rez in rezs) {
    rezInt = as.numeric(rez)
    
    signChanger[[rez]] = list()
    
    for (chrom in whichChroms[[whichChrom]]) {
        # rez = "500000"; rezInt = 500000; chrom="chr1"

        princompResults = princompResultsList[[rez]][[chrom]] 
        rmIdxAll = rmIdxList[[rez]][[chrom]]
        AllNans = rep(NaN, length(rmIdxAll))
        
        # Get signs of first component and used for all
        # This is pretty much guaranted to work for PC1, and higher components perhaps
        # But for lower components the sign og correlation may be not reflect similar or opposite natures
        signChanger[[rez]][[chrom]] = list()
        signChanger[[rez]][[chrom]][[names(princompResults)[1]]] = rep(1, nPCs)
        for (sampsamp in names(princompResults)[-1]) {
            signChanger[[rez]][[chrom]][[sampsamp]] = c()
            for (whichPC in 1:nPCs){
                signChanger[[rez]][[chrom]][[sampsamp]][whichPC] = sign(cor(princompResults[[1]]$scores[,whichPC], princompResults[[sampsamp]]$scores[,whichPC], method="pearson"))
            }
        }
        
        pdfName = file.path(scoreDir, paste(workDir, ".", format(rez, scientific=FALSE), 
                                             ".PCA.scores.", sampleCompDesc, ".", chrom, ".pdf", sep=""))  
        pdf(file=pdfName, height=(2*3)+2, width=(nSamples*4)+2)
        par(mfcol=c(2, nSamples), mar=par()$mar + c(1,1,0,0), oma=c(0,0,4,0))

        minMax = c(Inf, -Inf)
        for (sampsamp in names(princompResults)) {
            theRange = range(princompResults[[sampsamp]]$scores[,1])
            if (theRange[1] < minMax[1]) { minMax[1] = theRange[1] }
            if (theRange[2] > minMax[2]) { minMax[2] = theRange[2] }        
        }
        
        for (sampsamp in names(princompResults)) {
            # If Xi, then PC1 is black
            if (grepl("_ref$", sampsamp) && chrom=="chrX") {
                scores2Plot = AllNans
                scores2Plot[!rmIdxAll] = princompResults[[sampsamp]]$scores[,1] * signChanger[[rez]][[chrom]][[sampsamp]][1]
                plot(scores2Plot, 
                     type="l", col="black", lwd=2, 
                     main="", ylab="PC scores", ylim=minMax)
                scores2Plot = AllNans
                scores2Plot[!rmIdxAll] = princompResults[[sampsamp]]$scores[,2] * signChanger[[rez]][[chrom]][[sampsamp]][2]
                lines(scores2Plot, 
                      type="l", col="red", lwd=2)
                scores2Plot = AllNans
                scores2Plot[!rmIdxAll] = princompResults[[sampsamp]]$scores[,3] * signChanger[[rez]][[chrom]][[sampsamp]][3]
                lines(scores2Plot, 
                      type="l", col="blue", lwd=2)
                scores2Plot = AllNans
                scores2Plot[!rmIdxAll] = princompResults[[sampsamp]]$scores[,4] * signChanger[[rez]][[chrom]][[sampsamp]][4]
                lines(scores2Plot, 
                      type="l", col="green", lwd=2)
                legend("top", horiz=TRUE,  bty='n', xpd=TRUE, inset=c(0,-0.15),
                       legend=c("PC1", "PC2", "PC3", "PC4"),
                       #fill=c("black", "red", "blue", "green"),
                       col=c("black", "red", "blue", "green"),
                       lty=c("solid", "solid", "solid", "solid"),
                       lwd=c(2, 2, 2, 2))
            } else {
                scores2Plot = AllNans
                scores2Plot[!rmIdxAll] = princompResults[[sampsamp]]$scores[,1] * signChanger[[rez]][[chrom]][[sampsamp]][1]
                plot(scores2Plot, 
                     type="l", col="red", lwd=2, 
                     main="", ylab="PC scores", ylim=minMax)
                scores2Plot = AllNans
                scores2Plot[!rmIdxAll] = princompResults[[sampsamp]]$scores[,2] * signChanger[[rez]][[chrom]][[sampsamp]][2]
                lines(scores2Plot, 
                      type="l", col="blue", lwd=2)
                scores2Plot = AllNans
                scores2Plot[!rmIdxAll] = princompResults[[sampsamp]]$scores[,3] * signChanger[[rez]][[chrom]][[sampsamp]][3]
                lines(scores2Plot, 
                      type="l", col="green", lwd=2)
                legend("top", horiz=TRUE,  bty='n', xpd=TRUE, inset=c(0,-0.15),
                       legend=c("PC1", "PC2", "PC3"),
                       #fill=c("red", "blue", "green"),
                       col=c("red", "blue", "green"),
                       lty=c("solid", "solid", "solid"),
                       lwd=c(2, 2, 2))
            }
            title(main=sampsamp, line=3)
        }
        
        reset <- function() {
            par(mfrow=c(1, 1), oma=c(1,1,1,1), new=TRUE)
            plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)
        }
        reset()
        title(main=paste(chrom, " PCA scores (", sampleCompDesc2, "; ", 
                         format(rez, scientific=FALSE), ")", sep=""), line=3)
        dev.off()
    }
}


# ########    
# # Correlation matrices
# 
# for (rez in rezs) {
#     rezInt = as.numeric(rez)
#     
#     for (chrom in whichChroms[[whichChrom]]) {
#         # rez = "500000"; rezInt = 500000; chrom="chr1"
#         
#         princompResults = princompResultsList[[rez]][[chrom]] 
#         rmIdxAll = rmIdxList[[rez]][[chrom]]       
#         
#         corMatrixList=list()
#         for (PC in 1:3) {
#             corMatrixList[[as.character(PC)]]=c()
#             for (sampsampA in names(princompResults)) {
#                 for (sampsampB in names(princompResults)) {
#                     corMatrixList[[as.character(PC)]] = c(corMatrixList[[as.character(PC)]],
#                                                           abs(cor(princompResults[[sampsampA]]$scores[,PC], princompResults[[sampsampB]]$scores[,PC],
#                                                                   method="spearman"))
#                     )
#                 }
#             }
#             corMatrixList[[as.character(PC)]] = matrix(corMatrixList[[as.character(PC)]], ncol=4, nrow=4)
#             colnames(corMatrixList[[as.character(PC)]]) = names(princompResults)
#             rownames(corMatrixList[[as.character(PC)]]) = names(princompResults)
#         }
#         
#     }
# }


########    
# Pair plots - PCs per chrom

whichPCs = c(1:3)

# *** Must have run score line plot code above to populate signChanger[[rez]][[chrom]][[sampsamp]][]

pairplotDir = file.path(outputDir, "pairPlots")
if ( ! file.exists( pairplotDir ) ) {
    if ( ! dir.create( pairplotDir, showWarnings=FALSE, recursive=TRUE )  ) pairplotDir = outputDir
}

# http://r.789695.n4.nabble.com/color-in-pairs-td894955.html
# http://stackoverflow.com/questions/15625510/different-data-in-upper-and-lower-panel-of-scatterplot-matrix
panel.spearman <- function(x, y, cexbase, cexext, ...) {
    horizontal <- (par("usr")[1] + par("usr")[2]) / 2;
    vertical <- (par("usr")[3] + par("usr")[4]) / 2;
    corrValue = abs(cor(x,y, method="spearman")) 
    bgcolor <- colorRampPalette(c("white", "red"))(100)[round(corrValue*100)]
    u <- par("usr") # The coordinates of the plot area
    rect(u[1], u[3], u[2], u[4], col=bgcolor, border=NA)
    text(horizontal, vertical, 
         format(corrValue, digits=2),
         cex=cexbase + cexext*corrValue)
}

for (rez in rezs) {
    rezInt = as.numeric(rez)
    
    princompResults = princompResultsList[[rez]][[1]] 
    autosomalData4PairsList = list()
    for (PC in whichPCs) {
        autosomalData4PairsList[[as.character(PC)]]=c()
    }

    for (chrom in whichChroms[[whichChrom]]) {
        # rez = "500000"; rezInt = 500000; chrom="chr1"
  
        print(chrom)
        
        princompResults = princompResultsList[[rez]][[chrom]] 
        rmIdxAll = rmIdxList[[rez]][[chrom]]       

        data4PairsList = list()
        for (PC in whichPCs) {
            data4PairsList[[as.character(PC)]]=c()
            for (sampsamp in names(princompResults)) {
                data4PairsList[[as.character(PC)]] = cbind(data4PairsList[[as.character(PC)]], 
                                                           princompResults[[sampsamp]]$scores[,PC] * signChanger[[rez]][[chrom]][[sampsamp]][PC])
            }
            data4PairsList[[as.character(PC)]]  = as.data.frame(data4PairsList[[as.character(PC)]])
            colnames(data4PairsList[[as.character(PC)]]) = names(princompResults)

            if (chrom != "chrX") {
                autosomalData4PairsList[[as.character(PC)]] = rbind(autosomalData4PairsList[[as.character(PC)]], 
                                                                    data4PairsList[[as.character(PC)]])
            }        
        }
        
        for (PC in whichPCs) {
            # PC = 1
            pngName = file.path(pairplotDir, paste(workDir, ".", format(rez, scientific=FALSE), 
                                                   ".PCA.pairplots.Spearman.", sampleCompDesc, ".", chrom, ".PC", PC, ".png", sep=""))  
            png(file=pngName, height=(1*6)+2, width=(1*8)+2, units="in", res=300)
            pairs(data4PairsList[[as.character(PC)]], pch=20, col=gray(0, alpha=0.25),
              labels = sub("_ref", " B6", sub("_alt", " Sp", names(data4PairsList[[as.character(PC)]]))),
              cex.labels=2, font.labels=2,
              upper.panel=panel.spearman,
              cexbase=0.5, cexext=2,
              main=paste(chrom, " PC", as.character(PC), " pairwise Spearman correlations (", rez, "bp)", sep=""), cex.main=1)
            dev.off()
            
            pdfName = file.path(pairplotDir, paste(workDir, ".", format(rez, scientific=FALSE), 
                                                   ".PCA.pairplots.Spearman.", sampleCompDesc, ".", chrom, ".PC", PC, ".pdf", sep=""))  
            pdf(file=pdfName, height=(1*6)+2, width=(1*8)+2)
            pairs(data4PairsList[[as.character(PC)]], pch=20, col=gray(0, alpha=0.25),
                  labels = sub("_ref", " B6", sub("_alt", " Sp", names(data4PairsList[[as.character(PC)]]))),
                  cex.labels=2, font.labels=2,
                  upper.panel=panel.spearman,
                  cexbase=0.5, cexext=2,
                  main=paste(chrom, " PC", as.character(PC), " pairwise Spearman correlations (", rez, "bp)", sep=""), cex.main=1)
            dev.off()
            
        }
    }
    
    for (PC in whichPCs) {
        # PC = 1
        # Autosomes
        pngName = file.path(pairplotDir, paste(workDir, ".", format(rez, scientific=FALSE), 
                                               ".PCA.pairplots.Spearman.", sampleCompDesc, ".concatenatedAutosomes.PC", PC, ".png", sep=""))  
        png(file=pngName, height=(1*6)+2, width=(1*8)+2, units="in", res=300)
        pairs(autosomalData4PairsList[[as.character(PC)]], pch=20, col=gray(0, alpha=0.05),
              labels = sub("_ref", " B6", sub("_alt", " Sp", names(autosomalData4PairsList[[as.character(PC)]]))),
              cex.labels=2, font.labels=2,
              upper.panel=panel.spearman,
              cexbase=0.5, cexext=2,
              main=paste("concatenated autosomes PC", as.character(PC), " pairwise Spearman correlations (", rez, "bp)", sep=""), cex.main=1)
        dev.off()
        
        pdfName = file.path(pairplotDir, paste(workDir, ".", format(rez, scientific=FALSE), 
                                               ".PCA.pairplots.Spearman.", sampleCompDesc, ".concatenatedAutosomes.PC", PC, ".pdf", sep=""))  
        pdf(file=pdfName, height=(1*6)+2, width=(1*8)+2)
        pairs(autosomalData4PairsList[[as.character(PC)]], pch=20, col=gray(0, alpha=0.05),
              labels = sub("_ref", " B6", sub("_alt", " Sp", names(autosomalData4PairsList[[as.character(PC)]]))),
              cex.labels=2, font.labels=2,
              upper.panel=panel.spearman,
              cexbase=0.5, cexext=2,
              main=paste("concatenated autosomes PC", as.character(PC), " pairwise Spearman correlations (", rez, "bp)", sep=""), cex.main=1)
        dev.off()
    }
}


########    
# Pair plots - all PC1s in one plot 
# http://stackoverflow.com/questions/13081310/combining-multiple-complex-plots-as-panels-in-a-single-figure

PC = 1
# PC = 3

# unloadNamespace("gridGraphics")
# unloadNamespace("grid")
# unloadNamespace("gplots")
# unloadNamespace("gridExtra")
# unloadNamespace("ggplot2")
# unloadNamespace("lattice")

require(gridGraphics) # grid.echo, grid.grab
require(gridExtra)
require(grid)
require(ggplot2)
require(lattice)

# http://r.789695.n4.nabble.com/color-in-pairs-td894955.html
# http://stackoverflow.com/questions/15625510/different-data-in-upper-and-lower-panel-of-scatterplot-matrix
panel.spearman <- function(x, y, cexbase, cexext, ...) {
    horizontal <- (par("usr")[1] + par("usr")[2]) / 2;
    vertical <- (par("usr")[3] + par("usr")[4]) / 2;
    corrValue = abs(cor(x,y, method="spearman")) 
    bgcolor <- colorRampPalette(c("white", "red"))(100)[round(corrValue*100)]
    u <- par("usr") # The coordinates of the plot area
    rect(u[1], u[3], u[2], u[4], col=bgcolor, border=NA)
    text(horizontal, vertical, 
         format(corrValue, digits=2),
         cex=cexbase + cexext*corrValue)
}

for (rez in rezs) {
    rezInt = as.numeric(rez)

    gl = list()
    gltally=0

    for (chrom in whichChroms[[whichChrom]]) {
        # rez = "500000"; rezInt = 500000; chrom="chr1"
        
        print(chrom)
        
        gltally = gltally + 1

        princompResults = princompResultsList[[rez]][[chrom]] 
        rmIdxAll = rmIdxList[[rez]][[chrom]]       
        
        data4PairsList = list()
        data4PairsList[[as.character(PC)]]=c()
        for (sampsamp in names(princompResults)) {
            data4PairsList[[as.character(PC)]] = cbind(data4PairsList[[as.character(PC)]], 
                                                       princompResults[[sampsamp]]$scores[,PC] * signChanger[[rez]][[chrom]][[sampsamp]][1])
        }
        data4PairsList[[as.character(PC)]]  = as.data.frame(data4PairsList[[as.character(PC)]])
        colnames(data4PairsList[[as.character(PC)]]) = names(princompResults)

        pairs(data4PairsList[[as.character(PC)]], pch=20, col=gray(0, alpha=0.2),
              labels = sub("_ref", " B6", sub("_alt", " Sp", names(data4PairsList[[as.character(PC)]]))),
              cex.labels=2, font.labels=2,
              upper.panel=panel.spearman,
              cexbase=0.5, cexext=2,
              main=paste(chrom, sep=""), cex.main=1)
        grid.echo()
        gl[[gltally]] = grid.grab()     
    }
    
    # Test/Check:
    # grid.draw(gl[[1]])
    # grid.draw(gl[[2]])
    # grid.arrange(gl[[1]], gl[[2]], ncol=2, as.table=TRUE, clip=TRUE, main="test")
    # grid.arrange(arrangeGrob(gl[[1]], gl[[2]], ncol=2), gl[[3]], ncol=2, as.table=TRUE, clip=TRUE, main="test")
    
    pdfName = file.path(pairplotDir, paste(workDir, ".", format(rez, scientific=FALSE), 
                                         ".PC", as.character(PC), ".pairplots.Spearman.", sampleCompDesc, ".pdf", sep=""))  
    pdf(file=pdfName, height=5*4+2, width=4*5+2)
    #grid.newpage()
    grid.arrange(do.call("arrangeGrob", c(gl, ncol=4)), 
                 main=paste("PC", as.character(PC), " pairwise Spearman correlations (", rez, "bp)", sep=""))
    dev.off()
}





########    
# Pair plots - ChrX PCs
# http://stackoverflow.com/questions/13081310/combining-multiple-complex-plots-as-panels-in-a-single-figure

# http://r.789695.n4.nabble.com/color-in-pairs-td894955.html
# http://stackoverflow.com/questions/15625510/different-data-in-upper-and-lower-panel-of-scatterplot-matrix
panel.spearman <- function(x, y, cexbase, cexext, ...) {
    horizontal <- (par("usr")[1] + par("usr")[2]) / 2;
    vertical <- (par("usr")[3] + par("usr")[4]) / 2;
    corrValue = abs(cor(x,y, method="spearman")) 
    bgcolor <- colorRampPalette(c("white", "red"))(100)[round(corrValue*100)]
    u <- par("usr") # The coordinates of the plot area
    rect(u[1], u[3], u[2], u[4], col=bgcolor, border=NA)
    text(horizontal, vertical, 
         format(corrValue, digits=2),
         cex=cexbase + cexext*corrValue)
}

# panel.diagcolor <- function(x, ...){
#     u <- par("usr") # The coordinates of the plot area
#     rect(u[1], u[3], u[2], u[4], col="#FDFF65") 
# }

for (rez in rezs) {
    rezInt = as.numeric(rez)
    
    for (chrom in whichChroms[[whichChrom]]) {
        # rez = "500000"; rezInt = 500000; chrom="chr1"
        
        princompResults = princompResultsList[[rez]][[chrom]] 
        rmIdxAll = rmIdxList[[rez]][[chrom]]     
        
        data4Pairs=c()
        for (PC in whichPCs) {
            for (sampsamp in names(princompResults)) {
                data4Pairs= cbind(data4Pairs, 
                                  princompResults[[sampsamp]]$scores[,PC])
            }
        }
        data4Pairs = as.data.frame(data4Pairs)
        colnames(data4Pairs) = c( paste( "PC1", names(princompResults), sep="_" ),
                                  paste( "PC2", names(princompResults), sep="_" ),
                                  paste( "PC3", names(princompResults), sep="_" ))
        
        
        pdfName = file.path(outputDir, paste(workDir, ".", format(rez, scientific=FALSE), 
                                             ".", chrom, ".PCA.pairplots.Spearman.", sampleCompDesc, ".pdf", sep=""))  
        pdf(file=pdfName, height=3*4+2, width=3*4+2)
        #grid.newpage()
        
        if (chrom=="chrX") {
            diagLabels = sub("_", " ", sub("_ref", " Xi", sub("_alt", " Xa", names(data4Pairs))))
        } else {
            diagLabels = sub("_", " ", sub("_ref", " B6", sub("_alt", " Sp", names(data4Pairs))))
        }
        
        pairs(data4Pairs, pch=20, col=gray(0, alpha=0.2),
              labels = diagLabels,
              cex.labels=1, font.labels=2,
              upper.panel=panel.spearman,
              cexbase=0.5, cexext=2,
              # diag.panel=panel.diagcolor,
              main=paste(chrom, " PCA pairwise Spearman correlations (", rez, "bp)", sep=""), cex.main=1)
        
        dev.off()
    }
}
    



########    
# Save PC1 - PC5 scores

scoreTableDir = file.path(outputDir, "scoreTables")
if ( ! file.exists( scoreTableDir ) ) {
    if ( ! dir.create( scoreTableDir, showWarnings=FALSE, recursive=TRUE )  ) scoreTableDir = outputDir
}

signChanger = list()

for (rez in rezs) {
    rezInt = as.numeric(rez)
    
    signChanger[[rez]] = list()
    
    for (chrom in whichChroms[[whichChrom]]) {
        # rez = "500000"; rezInt = 500000; chrom="chr1"
        
        princompResults = princompResultsList[[rez]][[chrom]] 
        rmIdxAll = rmIdxList[[rez]][[chrom]]
        AllNans = rep(NaN, length(rmIdxAll))
        
        # Get signs of first component and used for all
        # This is pretty much guaranted to work for PC1, and higher components perhaps
        # But for lower components the sign og correlation may be not reflect similar or opposite natures
        signChanger[[rez]][[chrom]] = list()
        signChanger[[rez]][[chrom]][[names(princompResults)[1]]] = rep(1, nPCs)
        for (sampsamp in names(princompResults)[-1]) {
            signChanger[[rez]][[chrom]][[sampsamp]] = c()
            for (whichPC in 1:nPCs){
                signChanger[[rez]][[chrom]][[sampsamp]][whichPC] = sign(cor(princompResults[[1]]$scores[,whichPC], princompResults[[sampsamp]]$scores[,whichPC], method="pearson"))
            }
        }
 
        for (sampsamp in names(princompResults)) {
            scoresTable = c()
            for (PCnum in 1:5) {
                scores2Save = AllNans
                scores2Save[!rmIdxAll] = princompResults[[sampsamp]]$scores[,PCnum] * signChanger[[rez]][[chrom]][[sampsamp]][PCnum]
                scoresTable = cbind(scoresTable, scores2Save)
            }
            scoresTable = cbind( cbind(chrom,  
                                       format( seq(0, length(scores2Save)-1) * rezInt, scientific=FALSE ), 
                                       format( seq(1, length(scores2Save)) * rezInt, scientific=FALSE )
                                       ),
                                 scoresTable ) 
            colnames(scoresTable) = c("chr", "start", "end", paste("PC", c(1:5), sep=""))
            fileName = file.path(scoreTableDir, paste(workDir, ".", format(rez, scientific=FALSE), 
                                                    ".PCA.scores.", sampsamp, ".", chrom, ".tsv", sep=""))  
            write.table(scoresTable, fileName, row.names=FALSE, col.names=TRUE, sep="\t")
        }
    }
}



