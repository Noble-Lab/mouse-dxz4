######################################################################################################################
#
# Code to generate RPKMs from HTSeq counts
# Environment: Server
# Author: GAB
# 01/14/2015
#
######################################################################################################################


###########################################################
###########################################################
# Load master counts file

##########################################
#
# If uploading individual .out files from HTseq-count, use the code below:
#
# Create a merged dataframe from HTseq-count output files
# file_list        <- list.files()
# file_list
# countData       <- do.call("cbind", lapply(file_list, FUN=function(files){read.table(files,header=F, row.names=1, sep="\t")}))
# countData       <- head(countData, -5) 
# names(countData) <- list.files()
#
##########################################

inputDir="~/proj/2017_RNAseq_analysis/data/rawReads-20170620_RNAseq_InvDxz4abAndWTab/HTSeq"

masterFilePrefix = "Patski.combo.AllData"

# Read in raw counts table and remove the last 5 lines
countFile = file.path(inputDir, paste(masterFilePrefix, "ReadCounts.tsv.gz", sep="."))

countData <- read.table(gzfile(countFile), sep="\t", header=TRUE)

rownames(countData) = countData[,1]
countData = countData[,-1]

tail(countData)
# http://www-huber.embl.de/HTSeq/doc/count.html#count
# __no_feature: reads (or read pairs) which could not be assigned to any feature (set S as described above was empty).
# __ambiguous: reads (or read pairs) which could have been assigned to more than one feature and hence were not counted for any of these (set S had mroe than one element).
# __too_low_aQual: reads (or read pairs) which were skipped due to the -a option, see below
# __not_aligned: reads (or read pairs) in the SAM file without alignment
# __alignment_not_unique: reads (or read pairs) with more than one reported alignment. These reads are recognized from the NH optional SAM field tag. 
# (If the aligner does not set this field, multiply aligned reads will be counted multiple times, unless they getv filtered out by due to the -a option.)


# Total counts
colSums(countData)
# #               PatskiWT   PatskiWT_0uMaza_rep1   PatskiWT_0uMaza_rep2 
# 124286860              107645201              110258008 
# PatskiWT_4uMaza_rep1   PatskiWT_4uMaza_rep2             PatskiDel1 
# 114842003              109060673              129285405 
# PatskiDel1b PatskiDel1_0uMaza_rep1 PatskiDel1_0uMaza_rep2 
# 119134211              106336456              101304811 
# PatskiDel1_4uMaza_rep1 PatskiDel1_4uMaza_rep2             PatskiDel2 
# 104636061               97202505              128875899 
# PatskiDel5      PatskiWT_DelFirre    PatskiDel1_DelFirre 
# 120257943              119684577              122696310 
# PatskiDel5_DelFirre        patski.Fan.rep1        patski.Fan.rep2 
# 125185209               16652817               18293361 
# patski.FirreKDcntl         patski.FirreKD     patski.ENCODE.rep1 
# 21267438               22508856               16758515 
# patski.ENCODE.rep2         PatskiInvDxz4a         PatskiInvDxz4b 
# 21517619               56059600               51133985 
# PatskiWTa              PatskiWTb 
# 57104888               62945774 


# % ambiguous
countData["__ambiguous",] / colSums(countData) * 100
# PatskiWT PatskiWT_0uMaza_rep1 PatskiWT_0uMaza_rep2
# __ambiguous 1.319609             1.411233              1.44122
# PatskiWT_4uMaza_rep1 PatskiWT_4uMaza_rep2 PatskiDel1 PatskiDel1b
# __ambiguous             1.388972             1.382409   1.300359    1.375451
# PatskiDel1_0uMaza_rep1 PatskiDel1_0uMaza_rep2
# __ambiguous               1.361762               1.343605
# PatskiDel1_4uMaza_rep1 PatskiDel1_4uMaza_rep2 PatskiDel2 PatskiDel5
# __ambiguous               1.309667               1.306154   1.335713   1.599063
# PatskiWT_DelFirre PatskiDel1_DelFirre PatskiDel5_DelFirre
# __ambiguous          1.530178            1.442504            1.417194
# patski.Fan.rep1 patski.Fan.rep2 patski.FirreKDcntl patski.FirreKD
# __ambiguous        1.208354        1.225117           1.100433       1.136988
# patski.ENCODE.rep1 patski.ENCODE.rep2 PatskiInvDxz4a PatskiInvDxz4b
# __ambiguous           1.161165           1.169335       1.289392       1.331547
# PatskiWTa PatskiWTb
# __ambiguous  1.274551  1.282267



# % __alignment_not_unique
countData["__alignment_not_unique",] / colSums(countData) * 100
# PatskiWT PatskiWT_0uMaza_rep1 PatskiWT_0uMaza_rep2
# __alignment_not_unique 3.476783             3.416634             3.514535
# PatskiWT_4uMaza_rep1 PatskiWT_4uMaza_rep2 PatskiDel1
# __alignment_not_unique             5.645245             5.465469    3.55517
# PatskiDel1b PatskiDel1_0uMaza_rep1
# __alignment_not_unique    3.825733               3.297212
# PatskiDel1_0uMaza_rep2 PatskiDel1_4uMaza_rep1
# __alignment_not_unique               3.270219               3.595536
# PatskiDel1_4uMaza_rep2 PatskiDel2 PatskiDel5
# __alignment_not_unique               3.574519   3.512087   4.603424
# PatskiWT_DelFirre PatskiDel1_DelFirre
# __alignment_not_unique          4.594343            3.712883
# PatskiDel5_DelFirre patski.Fan.rep1 patski.Fan.rep2
# __alignment_not_unique             3.78554        10.68291         10.4653
# patski.FirreKDcntl patski.FirreKD patski.ENCODE.rep1
# __alignment_not_unique           8.332781       8.286934           5.690653
# patski.ENCODE.rep2 PatskiInvDxz4a PatskiInvDxz4b
# __alignment_not_unique           4.236152        3.48655       3.684857
# PatskiWTa PatskiWTb
# __alignment_not_unique  3.431065  3.508587


countData <- head(countData, -5) 
countData = data.matrix(countData)
head(countData)
nrow(countData)
# 24420

# ######################################################################################################################
# ######################################################################################################################
# # sample merging & RPKMs
# 
# expTableDir = file.path(inputDir, "expTables")
# if ( ! file.exists( expTableDir ) ) {
#   if ( ! dir.create( expTableDir, showWarnings=FALSE, recursive=TRUE )  ) expTableDir = inputDir
# }
# 
# # Link the HTSeq-generated counts table
# system( paste("ln -sf ", countFile, " ", expTableDir, sep="") )


###########################################################
###########################################################
# Get Gene lengths

# This is generated using the 'geneLengths.R' script
# geneLengthsFile = "~/refdata/hg19/hg19_UCSC_geneLengths.tsv"
geneLengthsFile = "~/refdata/iGenomes/Mus_musculus/UCSC/mm10/Annotation/Genes/GeneLengths.tsv"

# Read in gene lengths
# Ensure same annotation was used
geneLengths = read.table(file = geneLengthsFile, sep="\t", header=TRUE)
geneLengths = geneLengths[order(as.character(geneLengths$geneId)),]
head(geneLengths);  tail(geneLengths)
nrow(geneLengths)

# For some reason geneId for first entry of GTF is blank.
geneLengths = geneLengths[-1,]
nrow(geneLengths)
nrow(countData)

# Ensure the gene ids match up
matchRes = match(rownames(countData), as.character(geneLengths$geneId))
#all(1:nrow(countData) ==  matchRes)

# matchDiffs = which(1:nrow(countData) != matchRes)
# rownames(countData)[matchDiffs]
# as.character(geneLengths$geneId)[matchRes[matchDiffs]]
#all(rownames(countData)[matchDiffs] == as.character(geneLengths$geneId)[matchRes[matchDiffs]])

# Re-sort
geneLengths = geneLengths[matchRes, ]

# Check
matchResAfter = match(rownames(countData), as.character(geneLengths$geneId))
all(1:nrow(countData) ==  matchResAfter)
# TRUE


##########################################
##########################################
# RPKMS 

RPKMtable = c()
for (sample in colnames(countData)) {
  # sample = "gm0227_rep1_WAP021L8"
  cat(sample, "\n", sep="")
  totReads = sum(countData[, sample], na.rm=TRUE)
  RPKMs = (countData[, sample] / geneLengths[,2] ) *  (10^9 / totReads) 
  RPKMtable = cbind(RPKMtable, RPKMs)
}
colnames(RPKMtable) = colnames(countData)
head(RPKMtable)
nrow(RPKMtable)

colSums(RPKMtable)
# PatskiWT   PatskiWT_0uMaza_rep1   PatskiWT_0uMaza_rep2 
# 409243.8               421274.0               429559.1 
# PatskiWT_4uMaza_rep1   PatskiWT_4uMaza_rep2             PatskiDel1 
# 421080.8               418231.4               397743.6 
# PatskiDel1b PatskiDel1_0uMaza_rep1 PatskiDel1_0uMaza_rep2 
# 415124.4               410073.7               402560.5 
# PatskiDel1_4uMaza_rep1 PatskiDel1_4uMaza_rep2             PatskiDel2 
# 396275.0               392943.4               406357.4 
# PatskiDel5      PatskiWT_DelFirre    PatskiDel1_DelFirre 
# 433902.3               450953.7               430742.4 
# PatskiDel5_DelFirre        patski.Fan.rep1        patski.Fan.rep2 
# 427524.4               421118.2               416794.7 
# patski.FirreKDcntl         patski.FirreKD     patski.ENCODE.rep1 
# 389243.7               392041.3               380362.3 
# patski.ENCODE.rep2         PatskiInvDxz4a         PatskiInvDxz4b 
# 374564.8               378506.2               394378.1 
# PatskiWTa              PatskiWTb 
# 376705.4               379883.7 

outputFile = file.path(inputDir, paste(masterFilePrefix, "RPKMs.tsv", sep="."))
write.table(RPKMtable, file=outputFile, sep="\t", row.names=TRUE, col.names=NA)
system(paste('gzip -f ', outputFile, '\n', sep=""), intern=TRUE)



##########################################
##########################################
# TPMs
# Proportional RPKMs scaled to 1 million

TPMtable = RPKMtable 
for (colIdx in 1:ncol(TPMtable) ) {
  TPMtable[,colIdx] = TPMtable[,colIdx]/sum(TPMtable[,colIdx], na.rm=TRUE) * 1000000
}
colnames(TPMtable) = colnames(countData)
head(TPMtable)
nrow(TPMtable)

colSums(TPMtable)
# 1e+06                  1e+06                  1e+06 

outputFile = file.path(inputDir, paste(masterFilePrefix, "TPMs.tsv", sep="."))
write.table(TPMtable, file=outputFile, sep="\t", row.names=TRUE, col.names=NA)
system(paste('gzip -f ', outputFile, '\n', sep=""), intern=TRUE)



# ##########################################
# ##########################################
# # Merged RPKMS 
# 
# ##################
# # First count data
# 
# mergedCountTable = c()
# uniqueCellTypes = unique( sub( "_.*$", "", colnames(countData)) )
# 
# for (cellType in uniqueCellTypes) {
#   cat(cellType, "\t", sep="")
#   cellTypeIdx = grep( paste( "^", cellType, "_", sep=""), colnames(countData))
#   cat(length(cellTypeIdx), "\n", sep="")
#   if (length(cellTypeIdx)==1) {
#     mergedCounts = countData[, cellTypeIdx]
#   } else {
#     mergedCounts = rowSums(countData[, cellTypeIdx])
#   }
#   mergedCountTable = cbind(mergedCountTable, mergedCounts)
# }
# colnames(mergedCountTable) = uniqueCellTypes
# head(mergedCountTable)
# nrow(mergedCountTable)
# 
# outputFile = file.path(expTableDir, paste(projectName, ".MergedReadCounts.tsv", sep=""))
# write.table(mergedCountTable, file=outputFile, sep="\t", row.names=TRUE, col.names=NA)
# 
# ##################
# # Merge data RPKMs
# 
# mergedRPKMtable = c()
# for (sample in colnames(mergedCountTable)) {
#   cat(sample, "\n", sep="")
#   totReads = sum(mergedCountTable[, sample], na.rm=TRUE)
#   mergedRPKMs = (mergedCountTable[, sample] / geneLengths[,2]) * (10^9 / totReads) 
#   mergedRPKMtable = cbind(mergedRPKMtable, mergedRPKMs)
# }
# colnames(mergedRPKMtable) = colnames(mergedCountTable)
# head(mergedRPKMtable)
# nrow(mergedRPKMtable)
# 
# outputFile = file.path(expTableDir, paste(projectName, ".MergedRPKMs.tsv", sep=""))
# write.table(mergedRPKMtable, file=outputFile, sep="\t", row.names=TRUE, col.names=NA)
