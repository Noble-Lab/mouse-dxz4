##################################################
##################################################
##################################################
#
# 170707
# Overlap of CTCF and ATAC peaks
#
##################################################

source ~/.bashrc
source ~/.bash_profile

CTCFPROJDIR="<CTCFprojDirName>"
ATACPROJDIR="<ATACprojDirName>"

##################################################
# Working folder

WORKDIR="${CTCFPROJDIR}"/diploidPeakSegregation/CTCFvsATACpeaks
mkdir -p "${WORKDIR}"


##################################################
# Variables

sampleGenomes=( WTAmbig WTXi WTXa Patski2-4Ambig Patski2-4Xi Patski2-4Xa DelAmbig DelXi DelXa InvDxz4Ambig InvDxz4Xi InvDxz4Xa )

CTCFBEDDIR="${CTCFPROJDIR}/diploidPeakSegregation/peakAssignment/BEDfiles"
declare -A CTCF
CTCF=(
["WTAmbig"]="patskiWTvsDelvsInv.allelicCTCF.WT.ambiguous.minCov=5.adjPropThresh=0.7.chrX.bed.gz"
["WTXi"]="patskiWTvsDelvsInv.allelicCTCF.WT.B6-specific.minCov=5.adjPropThresh=0.7.chrX.bed.gz"
["WTXa"]="patskiWTvsDelvsInv.allelicCTCF.WT.Sp-specific.minCov=5.adjPropThresh=0.7.chrX.bed.gz"
["Patski2-4Ambig"]="patskiWTvsDelvsInv.allelicCTCF.2-4clone.ambiguous.minCov=5.adjPropThresh=0.7.chrX.bed.gz"
["Patski2-4Xi"]="patskiWTvsDelvsInv.allelicCTCF.2-4clone.B6-specific.minCov=5.adjPropThresh=0.7.chrX.bed.gz"
["Patski2-4Xa"]="patskiWTvsDelvsInv.allelicCTCF.2-4clone.Sp-specific.minCov=5.adjPropThresh=0.7.chrX.bed.gz"
["DelAmbig"]="patskiWTvsDelvsInv.allelicCTCF.Del.ambiguous.minCov=5.adjPropThresh=0.7.chrX.bed.gz"
["DelXi"]="patskiWTvsDelvsInv.allelicCTCF.Del.B6-specific.minCov=5.adjPropThresh=0.7.chrX.bed.gz"
["DelXa"]="patskiWTvsDelvsInv.allelicCTCF.Del.Sp-specific.minCov=5.adjPropThresh=0.7.chrX.bed.gz"
["InvDxz4Ambig"]="patskiWTvsDelvsInv.allelicCTCF.InvDxz4.ambiguous.minCov=5.adjPropThresh=0.7.chrX.bed.gz"
["InvDxz4Xi"]="patskiWTvsDelvsInv.allelicCTCF.InvDxz4.B6-specific.minCov=5.adjPropThresh=0.7.chrX.bed.gz"
["InvDxz4Xa"]="patskiWTvsDelvsInv.allelicCTCF.InvDxz4.Sp-specific.minCov=5.adjPropThresh=0.7.chrX.bed.gz"
)

ATACBEDDIR="${ATACPROJDIR}/diploidPeakSegregation/peakAssignment/BEDfiles"
declare -A ATAC
ATAC=(
["WTAmbig"]="patskiWTvsDelvsInvDxz4.allelicATAC.WT.ambiguous.minCov=5.adjPropThresh=0.7.chrX.bed.gz"
["WTXi"]="patskiWTvsDelvsInvDxz4.allelicATAC.WT.B6-specific.minCov=5.adjPropThresh=0.7.chrX.bed.gz"
["WTXa"]="patskiWTvsDelvsInvDxz4.allelicATAC.WT.Sp-specific.minCov=5.adjPropThresh=0.7.chrX.bed.gz"
["Patski2-4Ambig"]="patskiWTvsDelvsInvDxz4.allelicATAC.2-4clone.ambiguous.minCov=5.adjPropThresh=0.7.chrX.bed.gz"
["Patski2-4Xi"]="patskiWTvsDelvsInvDxz4.allelicATAC.2-4clone.B6-specific.minCov=5.adjPropThresh=0.7.chrX.bed.gz"
["Patski2-4Xa"]="patskiWTvsDelvsInvDxz4.allelicATAC.2-4clone.Sp-specific.minCov=5.adjPropThresh=0.7.chrX.bed.gz"
["DelAmbig"]="patskiWTvsDelvsInvDxz4.allelicATAC.Del.ambiguous.minCov=5.adjPropThresh=0.7.chrX.bed.gz"
["DelXi"]="patskiWTvsDelvsInvDxz4.allelicATAC.Del.B6-specific.minCov=5.adjPropThresh=0.7.chrX.bed.gz"
["DelXa"]="patskiWTvsDelvsInvDxz4.allelicATAC.Del.Sp-specific.minCov=5.adjPropThresh=0.7.chrX.bed.gz"
["InvDxz4Ambig"]="patskiWTvsDelvsInvDxz4.allelicATAC.InvDxz4.ambiguous.minCov=5.adjPropThresh=0.7.chrX.bed.gz"
["InvDxz4Xi"]="patskiWTvsDelvsInvDxz4.allelicATAC.InvDxz4.B6-specific.minCov=5.adjPropThresh=0.7.chrX.bed.gz"
["InvDxz4Xa"]="patskiWTvsDelvsInvDxz4.allelicATAC.InvDxz4.Sp-specific.minCov=5.adjPropThresh=0.7.chrX.bed.gz"
)


samples=( WT  Patski2-4 Del InvDxz4 )

CTCFDSCOREDIR="${CTCFPROJDIR}/diploidPeakSegregation/allelicRatios/adjAllelicPropBedGraphs"
declare -A CTCFDSCORE
CTCFDSCORE=(
["WT"]="patskiWTvsDelvsInv.allelicCTCF.WT.dScore.minCov=5.chrX.bedgraph.gz"
["Patski2-4"]="patskiWTvsDelvsInv.allelicCTCF.2-4clone.dScore.minCov=5.chrX.bedgraph.gz"
["Del"]="patskiWTvsDelvsInv.allelicCTCF.Del.dScore.minCov=5.chrX.bedgraph.gz"
["InvDxz4"]="patskiWTvsDelvsInv.allelicCTCF.InvDxz4.dScore.minCov=5.chrX.bedgraph.gz"
)

ATACDSCOREDIR="${ATACPROJDIR}/diploidPeakSegregation/allelicRatios/adjAllelicPropBedGraphs"
declare -A ATACDSCORE
ATACDSCORE=(
["WT"]="patskiWTvsDelvsInvDxz4.allelicATAC.WT.dScore.minCov=5.chrX.bedgraph.gz"
["Patski2-4"]="patskiWTvsDelvsInvDxz4.allelicATAC.2-4clone.dScore.minCov=5.chrX.bedgraph.gz"
["Del"]="patskiWTvsDelvsInvDxz4.allelicATAC.Del.dScore.minCov=5.chrX.bedgraph.gz"
["InvDxz4"]="patskiWTvsDelvsInvDxz4.allelicATAC.InvDxz4.dScore.minCov=5.chrX.bedgraph.gz"
)


##################################################
# chrX overlaps between samples - unsegregated

mkdir "${WORKDIR}"/chrXoverlapsSamplesUnseg
cd "${WORKDIR}"/chrXoverlapsSamplesUnseg

# module load bedtools/2.26.0

OUTFILE="CTCFoverlaps.tsv"

printf "GenomeA\tGenomeB\tAmeanPeakWidth\tAtally\tAoverlapTally\tAoverlap%%\tBmeanPeakWidth\tBtally\tBoverlapTally\tBoverlap%%\n" > ${OUTFILE}

for genomeA in ${samples[@]}; do
for genomeB in ${samples[@]}; do

Atally=$( zcat "${CTCFDSCOREDIR}"/"${CTCFDSCORE[${genomeA}]}" | wc -l )
Atally=$(( Atally - 1 ))
Ampw=$( zcat "${CTCFDSCOREDIR}"/"${CTCFDSCORE[${genomeA}]}" | awk '{totSize+=$3-$2+1} END{ print totSize/NR}' )
Btally=$( zcat "${CTCFDSCOREDIR}"/"${CTCFDSCORE[${genomeB}]}" | wc -l )
Btally=$(( Btally - 1 ))
Bmpw=$( zcat "${CTCFDSCOREDIR}"/"${CTCFDSCORE[${genomeB}]}" | awk '{totSize+=$3-$2+1} END{ print totSize/NR}' )
AoverlapTally=$( bedtools window -u -a "${CTCFDSCOREDIR}"/"${CTCFDSCORE[${genomeA}]}" -b "${CTCFDSCOREDIR}"/"${CTCFDSCORE[${genomeB}]}"  | wc -l )
BoverlapTally=$( bedtools window -u -a "${CTCFDSCOREDIR}"/"${CTCFDSCORE[${genomeB}]}" -b "${CTCFDSCOREDIR}"/"${CTCFDSCORE[${genomeA}]}" | wc -l )
AoverlapPerc=$( bc <<< "scale=2; ${AoverlapTally}/${Atally} * 100" )
BoverlapPerc=$( bc <<< "scale=2; ${BoverlapTally}/${Btally} * 100" )
printf "%s\t%s\t%0.2f\t%d\t%d\t%0.2f\t%0.2f\t%d\t%d\t%0.2f\n" ${genomeA} ${genomeB} ${Ampw} ${Atally} ${AoverlapTally} ${AoverlapPerc} ${Bmpw} ${Btally} ${BoverlapTally} ${BoverlapPerc}
done; done >> ${OUTFILE}

column -t ${OUTFILE} > ${OUTFILE}.txt
less  ${OUTFILE}.txt
grep  "^WT" ${OUTFILE}.txt



OUTFILE="ATACoverlaps.tsv"

printf "GenomeA\tGenomeB\tAmeanPeakWidth\tAtally\tAoverlapTally\tAoverlap%%\tBmeanPeakWidth\tBtally\tBoverlapTally\tBoverlap%%\n" > ${OUTFILE}

for genomeA in ${samples[@]}; do
for genomeB in ${samples[@]}; do

Atally=$( zcat "${ATACDSCOREDIR}"/"${ATACDSCORE[${genomeA}]}" | wc -l )
Atally=$(( Atally - 1 ))
Ampw=$( zcat "${ATACDSCOREDIR}"/"${ATACDSCORE[${genomeA}]}" | awk '{totSize+=$3-$2+1} END{ print totSize/NR}' )
Btally=$( zcat "${ATACDSCOREDIR}"/"${ATACDSCORE[${genomeB}]}" | wc -l )
Btally=$(( Btally - 1 ))
Bmpw=$( zcat "${ATACDSCOREDIR}"/"${ATACDSCORE[${genomeB}]}" | awk '{totSize+=$3-$2+1} END{ print totSize/NR}' )
AoverlapTally=$( bedtools window -u -a "${ATACDSCOREDIR}"/"${ATACDSCORE[${genomeA}]}" -b "${ATACDSCOREDIR}"/"${ATACDSCORE[${genomeB}]}"  | wc -l )
BoverlapTally=$( bedtools window -u -a "${ATACDSCOREDIR}"/"${ATACDSCORE[${genomeB}]}" -b "${ATACDSCOREDIR}"/"${ATACDSCORE[${genomeA}]}" | wc -l )
AoverlapPerc=$( bc <<< "scale=2; ${AoverlapTally}/${Atally} * 100" )
BoverlapPerc=$( bc <<< "scale=2; ${BoverlapTally}/${Btally} * 100" )
printf "%s\t%s\t%0.2f\t%d\t%d\t%0.2f\t%0.2f\t%d\t%d\t%0.2f\n" ${genomeA} ${genomeB} ${Ampw} ${Atally} ${AoverlapTally} ${AoverlapPerc} ${Bmpw} ${Btally} ${BoverlapTally} ${BoverlapPerc}
done; done >> ${OUTFILE}

column -t ${OUTFILE} > ${OUTFILE}.txt
less  ${OUTFILE}.txt
grep  "^WT" ${OUTFILE}.txt




##################################################
# chrX overlaps between samples - segregated

mkdir "${WORKDIR}"/chrXoverlapsSamples
cd "${WORKDIR}"/chrXoverlapsSamples

# module load bedtools/2.26.0

OUTFILE="CTCFoverlaps.tsv"

printf "GenomeA\tGenomeB\tAmeanPeakWidth\tAtally\tAoverlapTally\tAoverlap%%\tBmeanPeakWidth\tBtally\tBoverlapTally\tBoverlap%%\n" > ${OUTFILE}

for genomeA in ${sampleGenomes[@]}; do
for genomeB in ${sampleGenomes[@]}; do

Atally=$( zcat "${CTCFBEDDIR}"/"${CTCF[${genomeA}]}" | wc -l )
Atally=$(( Atally - 1 ))
Ampw=$( zcat "${CTCFBEDDIR}"/"${CTCF[${genomeA}]}" | awk '{totSize+=$3-$2+1} END{ print totSize/NR}' )
Btally=$( zcat "${CTCFBEDDIR}"/"${CTCF[${genomeB}]}" | wc -l )
Btally=$(( Btally - 1 ))
Bmpw=$( zcat "${CTCFBEDDIR}"/"${CTCF[${genomeB}]}" | awk '{totSize+=$3-$2+1} END{ print totSize/NR}' )
AoverlapTally=$( bedtools intersect -u -a "${CTCFBEDDIR}"/"${CTCF[${genomeA}]}" -b "${CTCFBEDDIR}"/"${CTCF[${genomeB}]}"  | wc -l )
BoverlapTally=$( bedtools intersect -u -a "${CTCFBEDDIR}"/"${CTCF[${genomeB}]}" -b "${CTCFBEDDIR}"/"${CTCF[${genomeA}]}" | wc -l )
AoverlapPerc=$( bc <<< "scale=2; ${AoverlapTally}/${Atally} * 100" )
BoverlapPerc=$( bc <<< "scale=2; ${BoverlapTally}/${Btally} * 100" )
printf "%s\t%s\t%0.2f\t%d\t%d\t%0.2f\t%0.2f\t%d\t%d\t%0.2f\n" ${genomeA} ${genomeB} ${Ampw} ${Atally} ${AoverlapTally} ${AoverlapPerc} ${Bmpw} ${Btally} ${BoverlapTally} ${BoverlapPerc}
done; done >> ${OUTFILE}

column -t ${OUTFILE} > ${OUTFILE}.txt
# less  ${OUTFILE}.txt
# grep  "^WT" ${OUTFILE}.txt



OUTFILE="ATACoverlaps.tsv"

printf "GenomeA\tGenomeB\tAmeanPeakWidth\tAtally\tAoverlapTally\tAoverlap%%\tBmeanPeakWidth\tBtally\tBoverlapTally\tBoverlap%%\n" > ${OUTFILE}

for genomeA in ${sampleGenomes[@]}; do
for genomeB in ${sampleGenomes[@]}; do

Atally=$( zcat "${ATACBEDDIR}"/"${ATAC[${genomeA}]}" | wc -l )
Atally=$(( Atally - 1 ))
Ampw=$( zcat "${ATACBEDDIR}"/"${ATAC[${genomeA}]}" | awk '{totSize+=$3-$2+1} END{ print totSize/NR}' )
Btally=$( zcat "${ATACBEDDIR}"/"${ATAC[${genomeB}]}" | wc -l )
Btally=$(( Btally - 1 ))
Bmpw=$( zcat "${ATACBEDDIR}"/"${ATAC[${genomeB}]}" | awk '{totSize+=$3-$2+1} END{ print totSize/NR}' )
AoverlapTally=$( bedtools intersect -u -a "${ATACBEDDIR}"/"${ATAC[${genomeA}]}" -b "${ATACBEDDIR}"/"${ATAC[${genomeB}]}"  | wc -l )
BoverlapTally=$( bedtools intersect -u -a "${ATACBEDDIR}"/"${ATAC[${genomeB}]}" -b "${ATACBEDDIR}"/"${ATAC[${genomeA}]}" | wc -l )
AoverlapPerc=$( bc <<< "scale=2; ${AoverlapTally}/${Atally} * 100" )
BoverlapPerc=$( bc <<< "scale=2; ${BoverlapTally}/${Btally} * 100" )
printf "%s\t%s\t%0.2f\t%d\t%d\t%0.2f\t%0.2f\t%d\t%d\t%0.2f\n" ${genomeA} ${genomeB} ${Ampw} ${Atally} ${AoverlapTally} ${AoverlapPerc} ${Bmpw} ${Btally} ${BoverlapTally} ${BoverlapPerc}
done; done >> ${OUTFILE}

column -t ${OUTFILE} > ${OUTFILE}.txt
# less  ${OUTFILE}.txt
# grep  "^WT" ${OUTFILE}.txt



##################################################
# chrX overlaps

mkdir "${WORKDIR}"/chrXoverlaps 
cd "${WORKDIR}"/chrXoverlaps

# module load bedtools/2.26.0

OUTFILE="CTCFvsATACoverlaps.tsv"

printf "Genome\tCTCFmeanPeakWidth\tCTCFtally\tCTCFoverlapTally\tCTCFoverlap%%\tATACmeanPeakWidth\tATACtally\tATACoverlapTally\tATACoverlap%%\n" > ${OUTFILE}

for genome in ${sampleGenomes[@]}; do
# genome=WTXi

#echo 
#echo "${CTCF[${genome}]}"
#ls -ltrh "${CTCFBEDDIR}"/"${CTCF[${genome}]}"
#echo "${ATACBEDDIR}"/"${ATAC[${genome}]}"
#ls -ltrh "${ATACBEDDIR}"/"${ATAC[${genome}]}"
#done 2>&1 | more

CTCFtally=$( zcat "${CTCFBEDDIR}"/"${CTCF[${genome}]}" | wc -l )
CTCFtally=$(( CTCFtally - 1 ))
CTCFmpw=$( zcat "${CTCFBEDDIR}"/"${CTCF[${genome}]}" | awk '{totSize+=$3-$2+1} END{ print totSize/NR}' )
ATACtally=$( zcat "${ATACBEDDIR}"/"${ATAC[${genome}]}" | wc -l )
ATACtally=$(( ATACtally - 1 ))
ATACmpw=$( zcat "${ATACBEDDIR}"/"${ATAC[${genome}]}" | awk '{totSize+=$3-$2+1} END{ print totSize/NR}' )
CTCFoverlapTally=$( bedtools window -w 1000 -u -a "${CTCFBEDDIR}"/"${CTCF[${genome}]}" -b "${ATACBEDDIR}"/"${ATAC[${genome}]}"  | wc -l )
ATACoverlapTally=$( bedtools window -w 1000 -u -a "${ATACBEDDIR}"/"${ATAC[${genome}]}" -b "${CTCFBEDDIR}"/"${CTCF[${genome}]}" | wc -l )
CTCFoverlapPerc=$( bc <<< "scale=2; ${CTCFoverlapTally}/${CTCFtally} * 100" )
ATACoverlapPerc=$( bc <<< "scale=2; ${ATACoverlapTally}/${ATACtally} * 100" )
printf "%s\t%0.2f\t%d\t%d\t%0.2f\t%0.2f\t%d\t%d\t%0.2f\n" ${genome} ${CTCFmpw} ${CTCFtally} ${CTCFoverlapTally} ${CTCFoverlapPerc} ${ATACmpw} ${ATACtally} ${ATACoverlapTally} ${ATACoverlapPerc}
done >> ${OUTFILE}

column -t ${OUTFILE} > ${OUTFILE}.txt
less  ${OUTFILE}.txt





##################################################
# Relative distance
# http://bedtools.readthedocs.io/en/latest/content/tools/reldist.html

mkdir "${WORKDIR}"/relDist 
cd "${WORKDIR}"/relDist

# Generate relative distances
for genome in ${sampleGenomes[@]}; do
# genome=WTXi
OUTFILE="CTCFvsATACreldist_${genome}.tsv"
bedtools reldist -a "${CTCFBEDDIR}"/"${CTCF[${genome}]}" -b "${ATACBEDDIR}"/"${ATAC[${genome}]}" > ${OUTFILE}
done 


# expected
bc <<< "scale=2; 1/50"
# 0.2


# Generate R script to plot figures
scriptName="plotRelDist.R"
cat << EOF > "${scriptName}"
options(echo=TRUE) # if you want see commands in output file

args <- commandArgs(trailingOnly = TRUE)
print(args)
# trailingOnly=TRUE means that only your arguments are returned, check:
# print(commandArgs(trailingOnly=FALSE))

fileName <- as.character(args[1])
figFolder <- as.character(args[2])
rm(args)

# Load data
relDistData = as.data.frame(read.table(fileName, header=TRUE, sep="\t")) #, fill = TRUE)

# Plots:
pdf(file.path(figFolder, paste(fileName,".pdf",sep="")), width=6, height=10)
par(mfrow=c(3,1), mar=c(2,4,0,1), oma=c(3,0,2,0))

plot(relDistData\$reldist, relDistData\$count, type="b",
main="", xlab="", ylab="counts", axes=TRUE, labels=FALSE)
axis(2, labels=TRUE)
abline(h=0.02*relDistData\$total[1], lty="dashed", col="red")
legend("topright", legend=paste("Total CTCF sites: ", relDistData\$total[1]))

plot(relDistData\$reldist, relDistData\$fraction, type="b",
main="", xlab="", ylab="proportion", axes=TRUE, labels=FALSE)
axis(2, labels=TRUE)
abline(h=0.02, lty="dashed", col="red")

plot(relDistData\$reldist, relDistData\$fraction/0.02, type="h",
main="", xlab="", ylab="O/E")
abline(h=1, lty="dashed", col="red")

title(main=fileName, outer=TRUE)
title(xlab="relative distance", outer=TRUE, line=1)

dev.off()

EOF
chmod 755 "${scriptName}"

# Execute scripts
mkdir pdfs
for genome in ${sampleGenomes[@]}; do
INFILE="CTCFvsATACreldist_${genome}.tsv"
Rscript "${scriptName}" "${INFILE}" "pdfs"
done




##################################################
# Absolute distance from Dxz4
#
# Dxz4=chrX:75637518-75764754
# chrX	171031299

mkdir "${WORKDIR}"/absDistFromDxz4
cd "${WORKDIR}"/absDistFromDxz4

# R script to calc distances
scriptName="plotDxz4Dist.R"
cat << EOF > "${scriptName}"
options(echo=TRUE) # if you want see commands in output file

args <- commandArgs(trailingOnly = TRUE)
print(args)
# trailingOnly=TRUE means that only your arguments are returned, check:
# print(commandArgs(trailingOnly=FALSE))

fileName <- as.character(args[1])
figFolder <- as.character(args[2])
dataDesc <- as.character(args[3])
binSize <- as.numeric(as.character(args[4]))
rm(args)

# Load data
peakData = as.data.frame(read.table(gzfile(fileName), header=FALSE, skip=1, sep="\t")) #, fill = TRUE)

# head(peakData)

Dxz4dists = cbind(
abs(as.numeric(as.character(peakData[,2])) - 75637518), 
abs(as.numeric(as.character(peakData[,2])) - 75764754),
abs(as.numeric(as.character(peakData[,3])) - 75637518), 
abs(as.numeric(as.character(peakData[,3])) - 75764754))

# head(Dxz4dists)

Dxz4dists = apply(Dxz4dists, 1, min, na.rm=TRUE)
Dxz4dists = Dxz4dists/binSize

# head(Dxz4dists)
# max(Dxz4dists)

# 171-75.5 = 95.5
expectedCounts = length(Dxz4dists)/(96 * (1000000/binSize))

# Plots:
pdf(file.path(figFolder, paste(dataDesc,".pdf",sep="")), width=12, height=8)
par(mfrow=c(2,1), mar=c(2,4,0,1), oma=c(3,0,2,0))

# Binned Counts
histData = hist(Dxz4dists, breaks=seq(0,(100 * (1000000/binSize)) + 1,1), main="", xlab="")
abline(h=expectedCounts, lty="dashed", col="red")

# Enrichment
mp = barplot(as.numeric(as.character(histData\$counts/expectedCounts)), ylab="O/E")
# length(mp)
# https://stackoverflow.com/questions/15595452/x-axis-in-barplot-in-r
axis(1, at=mp[seq(0, 100 * (1000000/binSize), 20 * (1000000/binSize)) + 1], labels=seq(0, 100, 20))
abline(h=1, lty="dashed", col="red") 

title(main=dataDesc, outer=TRUE)
title(xlab="distance to Dxz4 (Mb)", outer=TRUE, line=1)

dev.off()

EOF
chmod 755 "${scriptName}"


# plots
mkdir pdfs_1Mb
for genome in ${sampleGenomes[@]}; do
# genome=WTXi
INFILE="${CTCFBEDDIR}"/"${CTCF[${genome}]}"
Rscript "${scriptName}" "${INFILE}" "pdfs_1Mb" "Dxz4dists_CTCF_${genome}" "1000000"

INFILE="${ATACBEDDIR}"/"${ATAC[${genome}]}"
Rscript "${scriptName}" "${INFILE}" "pdfs_1Mb" "Dxz4dists_ATAC_${genome}" "1000000"
done

# plots
mkdir pdfs_100kb
for genome in ${sampleGenomes[@]}; do
# genome=WTXi
INFILE="${CTCFBEDDIR}"/"${CTCF[${genome}]}"
Rscript "${scriptName}" "${INFILE}" "pdfs_100kb" "Dxz4dists_CTCF_${genome}" "100000"

INFILE="${ATACBEDDIR}"/"${ATAC[${genome}]}"
Rscript "${scriptName}" "${INFILE}" "pdfs_100kb" "Dxz4dists_ATAC_${genome}" "100000"
done

# plots
mkdir pdfs_10kb
for genome in ${sampleGenomes[@]}; do
# genome=WTXi
INFILE="${CTCFBEDDIR}"/"${CTCF[${genome}]}"
Rscript "${scriptName}" "${INFILE}" "pdfs_10kb" "Dxz4dists_CTCF_${genome}" "10000"

INFILE="${ATACBEDDIR}"/"${ATAC[${genome}]}"
Rscript "${scriptName}" "${INFILE}" "pdfs_10kb" "Dxz4dists_ATAC_${genome}" "10000"
done






##################################################
# Positions along the X chromosome
#
# Dxz4=chrX:75637518-75764754
# Firre=chrX:50563120-50635321
# Xist=chrX:103460373-103483233 
# chrX	171031299

mkdir "${WORKDIR}"/binnedCounts
cd "${WORKDIR}"/binnedCounts

# R script to calc distances
scriptName="plotBinnedCountsAlongChrX.R"
cat << EOF > "${scriptName}"
options(echo=TRUE) # if you want see commands in output file

args <- commandArgs(trailingOnly = TRUE)
print(args)
# trailingOnly=TRUE means that only your arguments are returned, check:
# print(commandArgs(trailingOnly=FALSE))

fileName <- as.character(args[1])
figFolder <- as.character(args[2])
dataDesc <- as.character(args[3])
binSize <- as.numeric(as.character(args[4]))
rm(args)

# Load data
peakData = as.data.frame(read.table(gzfile(fileName), header=FALSE, skip=1, sep="\t")) #, fill = TRUE)

# head(peakData)

midPoints = as.numeric(as.character(peakData[,2])) + (as.numeric(as.character(peakData[,3])) - as.numeric(as.character(peakData[,2])))/2
midPoints = midPoints/binSize

firrebin = (50563120+(abs(50563120-50635321)/2)) / binSize
dxz4bin = (75764754+(abs(75637518-75764754)/2)) / binSize
xistbin = (103460373+(abs(103460373-103483233)/2)) / binSize

# head(midPoints)

# 171-75.5 = 95.5
expectedCounts = length(midPoints)/(171 * (1000000/binSize))

# Plots:
pdf(file.path(figFolder, paste(dataDesc,".pdf",sep="")), width=12, height=8)
par(mfrow=c(2,1), mar=c(2,4,3,1), oma=c(3,0,2,0))

# Binned counts
histData = hist(midPoints, breaks=seq(0,(180 * (1000000/binSize)) + 1,1), main="", xlab="")
abline(h=expectedCounts, lty="dashed", col="red")
# LOI
abline(v=firrebin, lty="dashed", lw=0.5, col="blue")
abline(v=dxz4bin, lty="dashed", lw=0.5, col="blue")
abline(v=xistbin, lty="dashed", lw=0.5, col="blue")
axis(3, at=c(firrebin, dxz4bin, xistbin), labels=c("Firre", "Dxz4", "Xist"), tick = FALSE)

# Enrichment
mp = barplot(as.numeric(as.character(histData\$counts/expectedCounts)), ylab="O/E")
# length(mp)
# https://stackoverflow.com/questions/15595452/x-axis-in-barplot-in-r
axis(1, at=mp[seq(0, 180 * (1000000/binSize), 20 * (1000000/binSize)) + 1], labels=seq(0, 180, 20))
abline(h=1, lty="dashed", col="red") 
# LOI
abline(v=mp[firrebin + 1], lty="dashed", lw=0.5, col="blue")
abline(v=mp[dxz4bin + 1], lty="dashed", lw=0.5, col="blue")
abline(v=mp[xistbin + 1], lty="dashed", lw=0.5, col="blue")
axis(3, at=mp[c(firrebin, dxz4bin, xistbin) + 1], labels=c("Firre", "Dxz4", "Xist"), tick = FALSE)

title(main=paste(dataDesc, " (", format(binSize, scientific=FALSE), "bp bins)", sep=""), outer=TRUE)
title(xlab="chrX (Mb)", outer=TRUE, line=1)

dev.off()

EOF
chmod 755 "${scriptName}"


# plots
mkdir pdfs_1Mb
for genome in ${sampleGenomes[@]}; do
# genome=WTXi
INFILE="${CTCFBEDDIR}"/"${CTCF[${genome}]}"
Rscript "${scriptName}" "${INFILE}" "pdfs_1Mb" "binnedCounts_CTCF_${genome}" "1000000"

INFILE="${ATACBEDDIR}"/"${ATAC[${genome}]}"
Rscript "${scriptName}" "${INFILE}" "pdfs_1Mb" "binnedCounts_ATAC_${genome}" "1000000"
done

# plots
mkdir pdfs_100kb
for genome in ${sampleGenomes[@]}; do
# genome=WTXi
INFILE="${CTCFBEDDIR}"/"${CTCF[${genome}]}"
Rscript "${scriptName}" "${INFILE}" "pdfs_100kb" "binnedCounts_CTCF_${genome}" "100000"

INFILE="${ATACBEDDIR}"/"${ATAC[${genome}]}"
Rscript "${scriptName}" "${INFILE}" "pdfs_100kb" "binnedCounts_ATAC_${genome}" "100000"
done

# plots
mkdir pdfs_10kb
for genome in ${sampleGenomes[@]}; do
# genome=WTXi
INFILE="${CTCFBEDDIR}"/"${CTCF[${genome}]}"
Rscript "${scriptName}" "${INFILE}" "pdfs_10kb" "binnedCounts_CTCF_${genome}" "10000"

INFILE="${ATACBEDDIR}"/"${ATAC[${genome}]}"
Rscript "${scriptName}" "${INFILE}" "pdfs_10kb" "binnedCounts_ATAC_${genome}" "10000"
done









##################################################
# Overlaid positions along the X chromosome
#
# Dxz4=chrX:75637518-75764754
# Firre=chrX:50563120-50635321
# Xist=chrX:103460373-103483233 
# chrX	171031299

mkdir "${WORKDIR}"/overlaidBinnedCounts
cd "${WORKDIR}"/overlaidBinnedCounts

# R script
scriptName="plotBinnedCountsAlongChrX.R"
cat << EOF > "${scriptName}"
options(echo=TRUE) # if you want see commands in output file

args <- commandArgs(trailingOnly = TRUE)
print(args)
# trailingOnly=TRUE means that only your arguments are returned, check:
# print(commandArgs(trailingOnly=FALSE))

ctcfFileName <- as.character(args[1])
atacFileName <- as.character(args[2])
figFolder <- as.character(args[3])
dataDesc <- as.character(args[4])
binSize <- as.numeric(as.character(args[5]))
rm(args)

# Load data
peakData = list()
peakData[["CTCF"]] = as.data.frame(read.table(gzfile(ctcfFileName), header=FALSE, skip=1, sep="\t")) #, fill = TRUE)
peakData[["ATAC"]] = as.data.frame(read.table(gzfile(atacFileName), header=FALSE, skip=1, sep="\t")) #, fill = TRUE)

peakCols =  list()
peakCols[["CTCF"]] = "blue"
peakCols[["ATAC"]] = "red"

peakAxis = list()
peakAxis[["CTCF"]] = 2
peakAxis[["ATAC"]] = 4

# head(peakData)

midPoints = list()
for (peakType in names(peakData)) {
midPoints[[peakType]] = as.numeric(as.character(peakData[[peakType]][,2])) + (as.numeric(as.character(peakData[[peakType]][,3])) - as.numeric(as.character(peakData[[peakType]][,2])))/2
midPoints[[peakType]] = midPoints[[peakType]]/binSize
}

firrebin = (50563120+(abs(50563120-50635321)/2)) / binSize
dxz4bin = (75764754+(abs(75637518-75764754)/2)) / binSize
xistbin = (103460373+(abs(103460373-103483233)/2)) / binSize

# head(midPoints)

# 171-75.5 = 95.5
expectedCounts = list()
for (peakType in names(midPoints)) {
expectedCounts[[peakType]] = length(midPoints[[peakType]])/(171 * (1000000/binSize))
}

# Plots:
pdf(file.path(figFolder, paste(dataDesc,".pdf",sep="")), width=15, height=10)
par(mfrow=c(2,1), mar=c(2,4,6,4), oma=c(3,0,2,0))

# Bin counts
histData = list()
for (peakType in names(midPoints)) {
histData[[peakType]] = hist(midPoints[[peakType]], breaks=seq(0,(180 * (1000000/binSize)) + 1,1), plot=FALSE)
}

# Plot counts
for (peakType in names(histData)) {
if (peakType=="ATAC") { par(new=TRUE) }
plot(as.numeric(as.character(histData[[peakType]]\$counts)), type="l", lty="solid", lwd=2, col=peakCols[[peakType]], axes=F, main=NA, xlab=NA, ylab=NA)
axis(side=peakAxis[[peakType]])
mtext(side=peakAxis[[peakType]], line=2.5, text=peakType)
abline(h=expectedCounts[[peakType]], lty="dashed", col=peakCols[[peakType]])
}
axis(1, at=seq(0, 180 * (1000000/binSize), 20 * (1000000/binSize)) + 1, labels=seq(0, 180, 20))
title(main="Counts", line=3)

Pcorr = cor(histData[["CTCF"]]\$counts, histData[["ATAC"]]\$counts, method="pearson")
Scorr = cor(histData[["CTCF"]]\$counts, histData[["ATAC"]]\$counts, method="spearman")

legend("topright", legend=c("CTCF", "ATAC", paste("Pearson = ", round(Pcorr, 2), sep=""), paste("Spearman = ", round(Scorr, 2), sep="")),
lty=c(1,1,NA,NA), col=c(peakCols[["CTCF"]], peakCols[["ATAC"]], NA, NA), bty='n', cex=0.75)

# LOI
abline(v=firrebin, lty="dashed", lw=0.5, col="black")
abline(v=dxz4bin, lty="dashed", lw=0.5, col="black")
abline(v=xistbin, lty="dashed", lw=0.5, col="black")
axis(3, at=c(firrebin, dxz4bin, xistbin), labels=c("Firre", "Dxz4", "Xist"), tick = FALSE)


# Plot enrichment
for (peakType in names(histData)) {
if (peakType=="ATAC") { par(new=TRUE) }
plot(as.numeric(as.character(histData[[peakType]]\$counts/expectedCounts[[peakType]])), type="l", lty="solid", lwd=2, col=peakCols[[peakType]], axes=F, main=NA, xlab=NA, ylab=NA)
axis(side=peakAxis[[peakType]])
mtext(side=peakAxis[[peakType]], line=2.5, text=peakType)
abline(h=1, lty="dashed", col=peakCols[[peakType]])
}
axis(1, at=seq(0, 180 * (1000000/binSize), 20 * (1000000/binSize)) + 1, labels=seq(0, 180, 20))
title(main="O/E", line=3)

# LOI
abline(v=firrebin, lty="dashed", lw=0.5, col="black")
abline(v=dxz4bin, lty="dashed", lw=0.5, col="black")
abline(v=xistbin, lty="dashed", lw=0.5, col="black")
axis(3, at=c(firrebin, dxz4bin, xistbin), labels=c("Firre", "Dxz4", "Xist"), tick = FALSE)

title(main=paste(dataDesc, " (", format(binSize, scientific=FALSE), "bp bins)", sep=""), outer=TRUE)
title(xlab="chrX (Mb)", outer=TRUE, line=1)

dev.off()

EOF
chmod 755 "${scriptName}"


# plots
mkdir pdfs_1Mb
for genome in ${sampleGenomes[@]}; do
# genome=WTXi
CTCFINFILE="${CTCFBEDDIR}"/"${CTCF[${genome}]}"
ATACINFILE="${ATACBEDDIR}"/"${ATAC[${genome}]}"
Rscript "${scriptName}" "${CTCFINFILE}" "${ATACINFILE}" "pdfs_1Mb" "binnedCounts_CTCFandATAC_${genome}" "1000000"
done

# plots
mkdir pdfs_100kb
for genome in ${sampleGenomes[@]}; do
# genome=WTXi
CTCFINFILE="${CTCFBEDDIR}"/"${CTCF[${genome}]}"
ATACINFILE="${ATACBEDDIR}"/"${ATAC[${genome}]}"
Rscript "${scriptName}" "${CTCFINFILE}" "${ATACINFILE}" "pdfs_100kb" "binnedCounts_CTCFandATAC_${genome}" "100000"
done

# plots
mkdir pdfs_10kb
for genome in ${sampleGenomes[@]}; do
# genome=WTXi
CTCFINFILE="${CTCFBEDDIR}"/"${CTCF[${genome}]}"
ATACINFILE="${ATACBEDDIR}"/"${ATAC[${genome}]}"
Rscript "${scriptName}" "${CTCFINFILE}" "${ATACINFILE}" "pdfs_10kb" "binnedCounts_CTCFandATAC_${genome}" "10000"
done





















##################################################
# Overlaid positions along the X chromosome V2
#
# Dxz4=chrX:75637518-75764754
# Firre=chrX:50563120-50635321
# Xist=chrX:103460373-103483233 
# chrX	171031299

mkdir "${WORKDIR}"/overlaidBinnedCountsV2
cd "${WORKDIR}"/overlaidBinnedCountsV2

# R script
scriptName="plotBinnedCountsAlongChrX.R"
cat << EOF > "${scriptName}"
options(echo=TRUE) # if you want see commands in output file

args <- commandArgs(trailingOnly = TRUE)
print(args)
# trailingOnly=TRUE means that only your arguments are returned, check:
# print(commandArgs(trailingOnly=FALSE))

WTfileName <- as.character(args[1])
Patski24FileName <- as.character(args[2])
DelFileName <- as.character(args[3])
InvDxz4FileName <- as.character(args[4])
figFolder <- as.character(args[5])
dataDesc <- as.character(args[6])
binSize <- as.numeric(as.character(args[7]))
rm(args)

# Load data
peakData = list()
peakData[["WT"]] = as.data.frame(read.table(gzfile(WTfileName), header=FALSE, skip=1, sep="\t")) #, fill = TRUE)
peakData[["Patski24"]] = as.data.frame(read.table(gzfile(Patski24FileName), header=FALSE, skip=1, sep="\t")) #, fill = TRUE)
peakData[["Del"]] = as.data.frame(read.table(gzfile(DelFileName), header=FALSE, skip=1, sep="\t")) #, fill = TRUE)
peakData[["InvDxz4"]] = as.data.frame(read.table(gzfile(InvDxz4FileName), header=FALSE, skip=1, sep="\t")) #, fill = TRUE)

peakCols =  list()
peakCols[["WT"]] = "blue"
peakCols[["Patski24"]] = "green"
peakCols[["Del"]] = "red"
peakCols[["InvDxz4"]] = "black"

# head(peakData)

midPoints = list()
for (peakType in names(peakData)) {
midPoints[[peakType]] = as.numeric(as.character(peakData[[peakType]][,2])) + (as.numeric(as.character(peakData[[peakType]][,3])) - as.numeric(as.character(peakData[[peakType]][,2])))/2
midPoints[[peakType]] = midPoints[[peakType]]/binSize
}

firrebin = (50563120+(abs(50563120-50635321)/2)) / binSize
dxz4bin = (75764754+(abs(75637518-75764754)/2)) / binSize
xistbin = (103460373+(abs(103460373-103483233)/2)) / binSize

# head(midPoints)

# Plots:
pdf(file.path(figFolder, paste(dataDesc,".pdf",sep="")), width=15, height=10)
par(mfrow=c(2,1), mar=c(2,4,6,4), oma=c(3,0,2,0))
theXlims = c(-5, 185 * (1000000/binSize))

# Bin counts
histData = list()
for (peakType in names(midPoints)) {
histData[[peakType]] = hist(midPoints[[peakType]], breaks=seq(0,(180 * (1000000/binSize)) + 1,1), plot=FALSE)
}

# Get ylims, total counts and expected counts
totalCounts = list()
expectedCounts = list()
ymax = -Inf 
Eymax = -Inf
for (peakType in names(histData)) {
totalCounts[[peakType]] = sum(as.numeric(as.character(histData[[peakType]]\$counts)))
expectedCounts[[peakType]] = length(midPoints[[peakType]])/(171 * (1000000/binSize)) # 171-75.5 = 95.5
if ( ymax < max(as.numeric(as.character(histData[[peakType]]\$counts)), na.rm=TRUE) ) {
ymax = max(as.numeric(as.character(histData[[peakType]]\$counts)), na.rm=TRUE)
}
if (Eymax < max(as.numeric(as.character(histData[[peakType]]\$counts/expectedCounts[[peakType]])), na.rm=TRUE)) {
Eymax = max(as.numeric(as.character(histData[[peakType]]\$counts/expectedCounts[[peakType]])), na.rm=TRUE)
}
}

# Plot WT counts
# Plot an empty graph and legend to get the size of the legend
# https://stackoverflow.com/questions/8929663/r-legend-placement-in-a-plot
plot(as.numeric(as.character(histData[["WT"]]\$counts)), type="n", xaxt="n", yaxt="n",  xlab=NA, ylab=NA, ylim=c(0,ymax), xlim=theXlims)
my.legend.size = legend("topleft", 
legend=c(
paste("WT (", totalCounts[["WT"]], " peaks)", sep=""),
paste("Patski2-4 (", totalCounts[["Patski24"]], " peaks)", sep=""),
paste("Del (", totalCounts[["Del"]], " peaks)", sep=""),
paste("InvDxz4 (", totalCounts[["InvDxz4"]], " peaks)", sep="")
),
lty=c(1,1,1,1), col=c(peakCols[["WT"]], peakCols[["Patski24"]], peakCols[["Del"]], peakCols[["InvDxz4"]]), bty='n', cex=1,
plot = FALSE)
# Custom ylim. Add the height of legend to upper bound of the range
my.range <- c(0, ymax)
my.range[2] <- 1.04*(my.range[2]+my.legend.size\$rect\$h)

# Draw the WT plot with custom ylim
par(new=TRUE)
plot(as.numeric(as.character(histData[["WT"]]\$counts)), type="l", lty="solid", lwd=2, col=peakCols[["WT"]], axes=F, main=NA, xlab=NA, ylab=NA, ylim=my.range, xlim=theXlims)
axis(side=2)
mtext(side=2, line=2.5, text="Counts")
axis(1, at=seq(0, 180 * (1000000/binSize), 20 * (1000000/binSize)) + 1, labels=seq(0, 180, 20))
title(main="Counts", line=3)

# Plot other samples
for (peakType in names(histData)) {
if (peakType != "WT") {
lines(as.numeric(as.character(histData[[peakType]]\$counts)), type="l", lty="solid", lwd=2, col=peakCols[[peakType]])
}
abline(h=expectedCounts[[peakType]], lty="dashed", col=peakCols[[peakType]])
}

# Legend
legend("topleft", 
legend=c(
paste("WT (", totalCounts[["WT"]], " peaks)", sep=""),
paste("Patski2-4 (", totalCounts[["Patski24"]], " peaks)", sep=""),
paste("Del (", totalCounts[["Del"]], " peaks)", sep=""),
paste("InvDxz4 (", totalCounts[["InvDxz4"]], " peaks)", sep="")
),
lty=c(1,1,1,1), col=c(peakCols[["WT"]], peakCols[["Patski24"]], peakCols[["Del"]], peakCols[["InvDxz4"]]), bty='n', cex=1)

Pcorr = c()
Scorr = c()
for (peakTypeA in names(peakData)) {
PcorrVect = c()
ScorrVect = c()
for (peakTypeB in names(peakData)) {
PcorrVect = c(PcorrVect, round(cor(histData[[peakTypeA]]\$counts, histData[[peakTypeB]]\$counts, method="pearson"),2))
ScorrVect = c(ScorrVect, round(cor(histData[[peakTypeA]]\$counts, histData[[peakTypeB]]\$counts, method="spearman"),2))
}
Pcorr =rbind(Pcorr, PcorrVect)
Scorr =rbind(Scorr, ScorrVect)
}
rownames(Pcorr) = names(peakData)
colnames(Pcorr) = names(peakData)
rownames(Scorr) = names(peakData)
colnames(Scorr) = names(peakData)

# add table
require(plotrix)
# https://stackoverflow.com/questions/18663159/conditional-coloring-of-cells-in-table
addtable2plot("topright", table=Scorr, bg="white", title="Spearman",
bty="o", display.rownames=TRUE, display.colnames=TRUE,
hlines=TRUE, vlines=TRUE,
xjust=2, yjust=1, cex=1)

# LOI
abline(v=firrebin, lty="dashed", lw=0.5, col="black")
abline(v=dxz4bin, lty="dashed", lw=0.5, col="black")
abline(v=xistbin, lty="dashed", lw=0.5, col="black")
axis(3, at=c(firrebin, dxz4bin, xistbin), labels=c("Firre", "Dxz4", "Xist"), tick = FALSE)


# Plot WT enrichment
# Plot an empty graph and legend to get the size of the legend
# https://stackoverflow.com/questions/8929663/r-legend-placement-in-a-plot
plot(as.numeric(as.character(histData[["WT"]]\$counts/expectedCounts[["WT"]])), type="n", xaxt="n", yaxt="n",  xlab=NA, ylab=NA, ylim=c(0,Eymax), xlim=theXlims)
my.legend.size = legend("topleft", 
legend=c(
paste("WT (", round(expectedCounts[["WT"]], 2), " peaks/bin)", sep=""),
paste("Patski2-4 (", round(expectedCounts[["Patski24"]], 2), " peaks/bin)", sep=""),
paste("Del (", round(expectedCounts[["Del"]], 2), " peaks/bin)", sep=""),
paste("InvDxz4 (", round(expectedCounts[["InvDxz4"]], 2), " peaks/bin)", sep="")
),
lty=c(1,1,1,1), col=c(peakCols[["WT"]], peakCols[["Patski24"]], peakCols[["Del"]], peakCols[["InvDxz4"]]), bty='n', cex=1,
plot = FALSE)
# Custom ylim. Add the height of legend to upper bound of the range
my.range <- c(0, Eymax)
my.range[2] <- 1.04*(my.range[2]+my.legend.size\$rect\$h)

# Draw the WT plot with custom ylim
par(new=TRUE)
plot(as.numeric(as.character(histData[["WT"]]\$counts/expectedCounts[["WT"]])), type="l", lty="solid", lwd=2, col=peakCols[["WT"]], axes=F, main=NA, xlab=NA, ylab=NA, ylim=my.range, xlim=theXlims)
axis(side=2)
mtext(side=2, line=2.5, text="O/E")
axis(1, at=seq(0, 180 * (1000000/binSize), 20 * (1000000/binSize)) + 1, labels=seq(0, 180, 20))
title(main="O/E", line=3)

# Legend
legend("topleft", 
legend=c(
paste("WT (", round(expectedCounts[["WT"]], 2),  " peaks/bin)", sep=""),
paste("Patski2-4 (", round(expectedCounts[["Patski24"]], 2), " peaks/bin)", sep=""),
paste("Del (", round(expectedCounts[["Del"]], 2),  " peaks/bin)", sep=""),
paste("InvDxz4 (", round(expectedCounts[["InvDxz4"]], 2),  " peaks/bin)", sep="")
),
lty=c(1,1,1,1), col=c(peakCols[["WT"]], peakCols[["Patski24"]], peakCols[["Del"]], peakCols[["InvDxz4"]]), bty='n', cex=1)

# Plot other enrichment
for (peakType in names(histData)) {
if (peakType != "WT") {
lines(as.numeric(as.character(histData[[peakType]]\$counts/expectedCounts[[peakType]])), type="l", lty="solid", lwd=2, col=peakCols[[peakType]])
}
abline(h=1, lty="dashed", col=peakCols[[peakType]])
}

# LOI
abline(v=firrebin, lty="dashed", lw=0.5, col="black")
abline(v=dxz4bin, lty="dashed", lw=0.5, col="black")
abline(v=xistbin, lty="dashed", lw=0.5, col="black")
axis(3, at=c(firrebin, dxz4bin, xistbin), labels=c("Firre", "Dxz4", "Xist"), tick = FALSE)

title(main=paste(dataDesc, " (", format(binSize, scientific=FALSE), "bp bins)", sep=""), outer=TRUE)
title(xlab="chrX (Mb)", outer=TRUE, line=1)

dev.off()

EOF
chmod 755 "${scriptName}"


# plots
mkdir pdfs_1Mb
# CTCF
Rscript "${scriptName}" "${CTCFBEDDIR}"/"${CTCF['WTXa']}" "${CTCFBEDDIR}"/"${CTCF['Patski2-4Xa']}" "${CTCFBEDDIR}"/"${CTCF['DelXa']}" "${CTCFBEDDIR}"/"${CTCF['InvDxz4Xa']}" "pdfs_1Mb" "binnedCounts_CTCF_Xa" "1000000"
Rscript "${scriptName}" "${CTCFBEDDIR}"/"${CTCF['WTXi']}" "${CTCFBEDDIR}"/"${CTCF['Patski2-4Xi']}" "${CTCFBEDDIR}"/"${CTCF['DelXi']}" "${CTCFBEDDIR}"/"${CTCF['InvDxz4Xi']}" "pdfs_1Mb" "binnedCounts_CTCF_Xi" "1000000"
Rscript "${scriptName}" "${CTCFBEDDIR}"/"${CTCF['WTAmbig']}" "${CTCFBEDDIR}"/"${CTCF['Patski2-4Ambig']}" "${CTCFBEDDIR}"/"${CTCF['DelAmbig']}" "${CTCFBEDDIR}"/"${CTCF['InvDxz4Ambig']}" "pdfs_1Mb" "binnedCounts_CTCF_Ambig" "1000000"
# ATAC
Rscript "${scriptName}" "${ATACBEDDIR}"/"${ATAC['WTXa']}" "${ATACBEDDIR}"/"${ATAC['Patski2-4Xa']}" "${ATACBEDDIR}"/"${ATAC['DelXa']}" "${ATACBEDDIR}"/"${ATAC['InvDxz4Xa']}" "pdfs_1Mb" "binnedCounts_ATAC_Xa" "1000000"
Rscript "${scriptName}" "${ATACBEDDIR}"/"${ATAC['WTXi']}" "${ATACBEDDIR}"/"${ATAC['Patski2-4Xi']}" "${ATACBEDDIR}"/"${ATAC['DelXi']}" "${ATACBEDDIR}"/"${ATAC['InvDxz4Xi']}" "pdfs_1Mb" "binnedCounts_ATAC_Xi" "1000000"
Rscript "${scriptName}" "${ATACBEDDIR}"/"${ATAC['WTAmbig']}" "${ATACBEDDIR}"/"${ATAC['Patski2-4Ambig']}" "${ATACBEDDIR}"/"${ATAC['DelAmbig']}" "${ATACBEDDIR}"/"${ATAC['InvDxz4Ambig']}" "pdfs_1Mb" "binnedCounts_ATAC_Ambig" "1000000"


# plots
mkdir pdfs_100kb
# CTCF
Rscript "${scriptName}" "${CTCFBEDDIR}"/"${CTCF['WTXa']}" "${CTCFBEDDIR}"/"${CTCF['Patski2-4Xa']}" "${CTCFBEDDIR}"/"${CTCF['DelXa']}" "${CTCFBEDDIR}"/"${CTCF['InvDxz4Xa']}" "pdfs_100kb" "binnedCounts_CTCF_Xa" "100000"
Rscript "${scriptName}" "${CTCFBEDDIR}"/"${CTCF['WTXi']}" "${CTCFBEDDIR}"/"${CTCF['Patski2-4Xi']}" "${CTCFBEDDIR}"/"${CTCF['DelXi']}" "${CTCFBEDDIR}"/"${CTCF['InvDxz4Xi']}" "pdfs_100kb" "binnedCounts_CTCF_Xi" "100000"
Rscript "${scriptName}" "${CTCFBEDDIR}"/"${CTCF['WTAmbig']}" "${CTCFBEDDIR}"/"${CTCF['Patski2-4Ambig']}" "${CTCFBEDDIR}"/"${CTCF['DelAmbig']}" "${CTCFBEDDIR}"/"${CTCF['InvDxz4Ambig']}" "pdfs_100kb" "cuminnedCounts_CTCF_Ambig" "100000"
# ATAC
Rscript "${scriptName}" "${ATACBEDDIR}"/"${ATAC['WTXa']}" "${ATACBEDDIR}"/"${ATAC['Patski2-4Xa']}" "${ATACBEDDIR}"/"${ATAC['DelXa']}" "${ATACBEDDIR}"/"${ATAC['InvDxz4Xa']}" "pdfs_100kb" "binnedCounts_ATAC_Xa" "100000"
Rscript "${scriptName}" "${ATACBEDDIR}"/"${ATAC['WTXi']}" "${ATACBEDDIR}"/"${ATAC['Patski2-4Xi']}" "${ATACBEDDIR}"/"${ATAC['DelXi']}" "${ATACBEDDIR}"/"${ATAC['InvDxz4Xi']}" "pdfs_100kb" "binnedCounts_ATAC_Xi" "100000"
Rscript "${scriptName}" "${ATACBEDDIR}"/"${ATAC['WTAmbig']}" "${ATACBEDDIR}"/"${ATAC['Patski2-4Ambig']}" "${ATACBEDDIR}"/"${ATAC['DelAmbig']}" "${ATACBEDDIR}"/"${ATAC['InvDxz4Ambig']}" "pdfs_100kb" "binnedCounts_ATAC_Ambig" "100000"












##################################################
# 'CDFs' wrt distance along the X chromosome
# Dxz4=chrX:75637518-75764754
# Firre=chrX:50563120-50635321
# Xist=chrX:103460373-103483233 
# chrX	171031299

mkdir "${WORKDIR}"/overlaidCumBinnedCounts
cd "${WORKDIR}"/overlaidCumBinnedCounts

# R script
# https://stackoverflow.com/questions/2151212/how-can-i-read-command-line-parameters-from-an-r-script#2154190
scriptName="plotCumulativeBinnedCountsAlongChrX.R"
cat << EOF > "${scriptName}"
options(echo=TRUE) # if you want see commands in output file

args <- commandArgs(trailingOnly = TRUE)
print(args)
# trailingOnly=TRUE means that only your arguments are returned, check:
# print(commandArgs(trailingOnly=FALSE))

WTfileName <- as.character(args[1])
Patski24FileName <- as.character(args[2])
DelFileName <- as.character(args[3])
InvDxz4FileName <- as.character(args[4])
figFolder <- as.character(args[5])
dataDesc <- as.character(args[6])
binSize <- as.numeric(as.character(args[7]))
rm(args)

# Load data
peakData = list()
peakData[["WT"]] = as.data.frame(read.table(gzfile(WTfileName), header=FALSE, skip=1, sep="\t")) #, fill = TRUE)
peakData[["Patski24"]] = as.data.frame(read.table(gzfile(Patski24FileName), header=FALSE, skip=1, sep="\t")) #, fill = TRUE)
peakData[["Del"]] = as.data.frame(read.table(gzfile(DelFileName), header=FALSE, skip=1, sep="\t")) #, fill = TRUE)
peakData[["InvDxz4"]] = as.data.frame(read.table(gzfile(InvDxz4FileName), header=FALSE, skip=1, sep="\t")) #, fill = TRUE)

peakCols =  list()
peakCols[["WT"]] = "blue"
peakCols[["Patski24"]] = "green"
peakCols[["Del"]] = "red"
peakCols[["InvDxz4"]] = "black"

# head(peakData)

midPoints = list()
for (peakType in names(peakData)) {
midPoints[[peakType]] = as.numeric(as.character(peakData[[peakType]][,2])) + (as.numeric(as.character(peakData[[peakType]][,3])) - as.numeric(as.character(peakData[[peakType]][,2])))/2
midPoints[[peakType]] = midPoints[[peakType]]/binSize
}

firrebin = (50563120+(abs(50563120-50635321)/2)) / binSize
dxz4bin = (75764754+(abs(75637518-75764754)/2)) / binSize
xistbin = (103460373+(abs(103460373-103483233)/2)) / binSize

# head(midPoints)

# Plots:
pdf(file.path(figFolder, paste(dataDesc,".pdf",sep="")), width=15, height=10)
par(mfrow=c(2,1), mar=c(2,4,6,4), oma=c(3,0,2,0))
theXlims = c(-5, 185 * (1000000/binSize))

# Bin counts
histData = list()
for (peakType in names(midPoints)) {
histData[[peakType]] = hist(midPoints[[peakType]], breaks=seq(0,(180 * (1000000/binSize)) + 1,1), plot=FALSE)
}

# Get ylims, total counts and expected counts
totalCounts = list()
expectedCounts = list()
cumCounts = list()
ymax = -Inf 
Eymax = -Inf
for (peakType in names(histData)) {
totalCounts[[peakType]] = sum(as.numeric(as.character(histData[[peakType]]\$counts)))
expectedCounts[[peakType]] = length(midPoints[[peakType]])/(171 * (1000000/binSize)) # 171-75.5 = 95.5
cumCounts[[peakType]] = cumsum(as.numeric(as.character(histData[[peakType]]\$counts)))
if ( ymax < max(as.numeric(as.character(histData[[peakType]]\$counts)), na.rm=TRUE) ) {
ymax = max(as.numeric(as.character(histData[[peakType]]\$counts)), na.rm=TRUE)
}
if (Eymax < max(as.numeric(as.character(histData[[peakType]]\$counts/expectedCounts[[peakType]])), na.rm=TRUE)) {
Eymax = max(as.numeric(as.character(histData[[peakType]]\$counts/expectedCounts[[peakType]])), na.rm=TRUE)
}
}

cumProp = list()
totMax = -Inf
for (peakType in names(histData)) {
cumProp[[peakType]] = cumCounts[[peakType]] / totalCounts[[peakType]]
if ( totMax < totalCounts[[peakType]] ) {
totMax = totalCounts[[peakType]]
}
}

# Plot WT counts
# Plot an empty graph and legend to get the size of the legend
# https://stackoverflow.com/questions/8929663/r-legend-placement-in-a-plot
plot(cumCounts[["WT"]], type="n", xaxt="n", yaxt="n",  xlab=NA, ylab=NA, ylim=c(0,ymax), xlim=theXlims)
my.legend.size = legend("topleft", 
legend=c(
paste("WT (", totalCounts[["WT"]], " peaks)", sep=""),
paste("Del (", totalCounts[["Del"]], " peaks)", sep=""),
paste("InvDxz4 (", totalCounts[["InvDxz4"]], " peaks)", sep="")
),
lty=c(1,1,1,1), col=c(peakCols[["WT"]], peakCols[["Patski24"]], peakCols[["Del"]], peakCols[["InvDxz4"]]), bty='n', cex=1,
plot = FALSE)
# Custom ylim. Add the height of legend to upper bound of the range
my.range <- c(0, totMax)
my.range[2] <- 1.04*(my.range[2]+my.legend.size\$rect\$h)
my.range <- c(0, totMax)

# Draw the WT plot with custom ylim
par(new=TRUE)
plot(cumCounts[["WT"]], type="l", lty="solid", lwd=2, col=peakCols[["WT"]], axes=F, main=NA, xlab=NA, ylab=NA, ylim=my.range, xlim=theXlims)
axis(side=2)
mtext(side=2, line=2.5, text="Counts")
axis(1, at=seq(0, 180 * (1000000/binSize), 20 * (1000000/binSize)) + 1, labels=seq(0, 180, 20))
title(main="Counts", line=3)

# Plot other samples
for (peakType in names(histData)) {
if (peakType != "WT") {
lines(cumCounts[[peakType]], type="l", lty="solid", lwd=2, col=peakCols[[peakType]])
}
abline(h=expectedCounts[[peakType]], lty="dashed", col=peakCols[[peakType]])
}

# Legends
legend("topleft", 
legend=c(
paste("WT (", totalCounts[["WT"]], " peaks)", sep=""),
paste("Patski2-4 (", totalCounts[["Patski24"]], " peaks)", sep=""),
paste("Del (", totalCounts[["Del"]], " peaks)", sep=""),
paste("InvDxz4 (", totalCounts[["InvDxz4"]], " peaks)", sep="")
),
lty=c(1,1,1,1), col=c(peakCols[["WT"]], peakCols[["Patski24"]], peakCols[["Del"]], peakCols[["InvDxz4"]]), bty='n', cex=1)

Pcorr = c()
Scorr = c()
for (peakTypeA in names(peakData)) {
PcorrVect = c()
ScorrVect = c()
for (peakTypeB in names(peakData)) {
PcorrVect = c(PcorrVect, round(cor(histData[[peakTypeA]]\$counts, histData[[peakTypeB]]\$counts, method="pearson"),2))
ScorrVect = c(ScorrVect, round(cor(histData[[peakTypeA]]\$counts, histData[[peakTypeB]]\$counts, method="spearman"),2))
}
Pcorr =rbind(Pcorr, PcorrVect)
Scorr =rbind(Scorr, ScorrVect)
}
rownames(Pcorr) = names(peakData)
colnames(Pcorr) = names(peakData)
rownames(Scorr) = names(peakData)
colnames(Scorr) = names(peakData)

# add table
require(plotrix)
# https://stackoverflow.com/questions/18663159/conditional-coloring-of-cells-in-table
addtable2plot("left", table=Scorr, bg="white", title="Spearman",
bty="o", display.rownames=TRUE, display.colnames=TRUE,
hlines=TRUE, vlines=TRUE,
xjust=2, yjust=1, cex=1)

# LOI
abline(v=firrebin, lty="dashed", lw=0.5, col="black")
abline(v=dxz4bin, lty="dashed", lw=0.5, col="black")
abline(v=xistbin, lty="dashed", lw=0.5, col="black")
axis(3, at=c(firrebin, dxz4bin, xistbin), labels=c("Firre", "Dxz4", "Xist"), tick = FALSE)


# Plot WT Proportions
# Plot an empty graph and legend to get the size of the legend
# https://stackoverflow.com/questions/8929663/r-legend-placement-in-a-plot
plot(cumProp[["WT"]], type="n", xaxt="n", yaxt="n",  xlab=NA, ylab=NA, ylim=c(0,1), xlim=theXlims)
my.legend.size = legend("topleft", 
legend=c(
paste("WT (", round(expectedCounts[["WT"]], 2), " peaks/bin)", sep=""),
paste("Patski2-4 (", round(expectedCounts[["Patski24"]], 2), " peaks/bin)", sep=""),
paste("Del (", round(expectedCounts[["Del"]], 2), " peaks/bin)", sep=""),
paste("InvDxz4 (", round(expectedCounts[["InvDxz4"]], 2), " peaks/bin)", sep="")
),
lty=c(1,1,1,1), col=c(peakCols[["WT"]], peakCols[["Patski24"]], peakCols[["Del"]], peakCols[["InvDxz4"]]), bty='n', cex=1,
plot = FALSE)
# Custom ylim. Add the height of legend to upper bound of the range
my.range <- c(0, 1)
my.range[2] <- 1.04*(my.range[2]+my.legend.size\$rect\$h)
my.range <- c(0, 1)

# Draw the WT plot with custom ylim
par(new=TRUE)
plot(cumProp[["WT"]], type="l", lty="solid", lwd=2, col=peakCols[["WT"]], axes=F, main=NA, xlab=NA, ylab=NA, ylim=my.range, xlim=theXlims)
axis(side=2)
mtext(side=2, line=2.5, text="Proportion")
axis(1, at=seq(0, 180 * (1000000/binSize), 20 * (1000000/binSize)) + 1, labels=seq(0, 180, 20))
title(main="Proportion", line=3)

# Legend
legend("topleft", 
legend=c(
paste("WT (", round(expectedCounts[["WT"]], 2),  " peaks/bin)", sep=""),
paste("Patski2-4 (", round(expectedCounts[["Patski24"]], 2), " peaks/bin)", sep=""),
paste("Del (", round(expectedCounts[["Del"]], 2),  " peaks/bin)", sep=""),
paste("InvDxz4 (", round(expectedCounts[["InvDxz4"]], 2),  " peaks/bin)", sep="")
),
lty=c(1,1,1,1), col=c(peakCols[["WT"]], peakCols[["Patski24"]], peakCols[["Del"]], peakCols[["InvDxz4"]]), bty='n', cex=1)

# Plot other enrichment
for (peakType in names(histData)) {
if (peakType != "WT") {
lines(cumProp[[peakType]], type="l", lty="solid", lwd=2, col=peakCols[[peakType]])
}
abline(h=1, lty="dashed", col=peakCols[[peakType]])
}

# LOI
abline(v=firrebin, lty="dashed", lw=0.5, col="black")
abline(v=dxz4bin, lty="dashed", lw=0.5, col="black")
abline(v=xistbin, lty="dashed", lw=0.5, col="black")
axis(3, at=c(firrebin, dxz4bin, xistbin), labels=c("Firre", "Dxz4", "Xist"), tick = FALSE)

title(main=paste(dataDesc, " (", format(binSize, scientific=FALSE), "bp bins)", sep=""), outer=TRUE)
title(xlab="chrX (Mb)", outer=TRUE, line=1)

dev.off()

EOF
chmod 755 "${scriptName}"


# plots
mkdir pdfs_1Mb
# CTCF
Rscript "${scriptName}" "${CTCFBEDDIR}"/"${CTCF['WTXa']}" "${CTCFBEDDIR}"/"${CTCF['Patski2-4Xa']}" "${CTCFBEDDIR}"/"${CTCF['DelXa']}" "${CTCFBEDDIR}"/"${CTCF['InvDxz4Xa']}" "pdfs_1Mb" "cumbinnedCounts_CTCF_Xa" "1000000"
Rscript "${scriptName}" "${CTCFBEDDIR}"/"${CTCF['WTXi']}" "${CTCFBEDDIR}"/"${CTCF['Patski2-4Xi']}" "${CTCFBEDDIR}"/"${CTCF['DelXi']}" "${CTCFBEDDIR}"/"${CTCF['InvDxz4Xi']}" "pdfs_1Mb" "cumbinnedCounts_CTCF_Xi" "1000000"
Rscript "${scriptName}" "${CTCFBEDDIR}"/"${CTCF['WTAmbig']}" "${CTCFBEDDIR}"/"${CTCF['Patski2-4Ambig']}" "${CTCFBEDDIR}"/"${CTCF['DelAmbig']}" "${CTCFBEDDIR}"/"${CTCF['InvDxz4Ambig']}" "pdfs_1Mb" "cumbinnedCounts_CTCF_Ambig" "1000000"
# ATAC
Rscript "${scriptName}" "${ATACBEDDIR}"/"${ATAC['WTXa']}" "${ATACBEDDIR}"/"${ATAC['Patski2-4Xa']}" "${ATACBEDDIR}"/"${ATAC['DelXa']}" "${ATACBEDDIR}"/"${ATAC['InvDxz4Xa']}" "pdfs_1Mb" "cumbinnedCounts_ATAC_Xa" "1000000"
Rscript "${scriptName}" "${ATACBEDDIR}"/"${ATAC['WTXi']}" "${ATACBEDDIR}"/"${ATAC['Patski2-4Xi']}" "${ATACBEDDIR}"/"${ATAC['DelXi']}" "${ATACBEDDIR}"/"${ATAC['InvDxz4Xi']}" "pdfs_1Mb" "cumbinnedCounts_ATAC_Xi" "1000000"
Rscript "${scriptName}" "${ATACBEDDIR}"/"${ATAC['WTAmbig']}" "${ATACBEDDIR}"/"${ATAC['Patski2-4Ambig']}" "${ATACBEDDIR}"/"${ATAC['DelAmbig']}" "${ATACBEDDIR}"/"${ATAC['InvDxz4Ambig']}" "pdfs_1Mb" "cumbinnedCounts_ATAC_Ambig" "1000000"


# plots
mkdir pdfs_100kb
# CTCF
Rscript "${scriptName}" "${CTCFBEDDIR}"/"${CTCF['WTXa']}" "${CTCFBEDDIR}"/"${CTCF['Patski2-4Xa']}" "${CTCFBEDDIR}"/"${CTCF['DelXa']}" "${CTCFBEDDIR}"/"${CTCF['InvDxz4Xa']}" "pdfs_100kb" "cumbinnedCounts_CTCF_Xa" "100000"
Rscript "${scriptName}" "${CTCFBEDDIR}"/"${CTCF['WTXi']}" "${CTCFBEDDIR}"/"${CTCF['Patski2-4Xi']}" "${CTCFBEDDIR}"/"${CTCF['DelXi']}" "${CTCFBEDDIR}"/"${CTCF['InvDxz4Xi']}" "pdfs_100kb" "cumbinnedCounts_CTCF_Xi" "100000"
Rscript "${scriptName}" "${CTCFBEDDIR}"/"${CTCF['WTAmbig']}" "${CTCFBEDDIR}"/"${CTCF['Patski2-4Ambig']}" "${CTCFBEDDIR}"/"${CTCF['DelAmbig']}" "${CTCFBEDDIR}"/"${CTCF['InvDxz4Ambig']}" "pdfs_100kb" "cuminnedCounts_CTCF_Ambig" "100000"
# ATAC
Rscript "${scriptName}" "${ATACBEDDIR}"/"${ATAC['WTXa']}" "${ATACBEDDIR}"/"${ATAC['Patski2-4Xa']}" "${ATACBEDDIR}"/"${ATAC['DelXa']}" "${ATACBEDDIR}"/"${ATAC['InvDxz4Xa']}" "pdfs_100kb" "cumbinnedCounts_ATAC_Xa" "100000"
Rscript "${scriptName}" "${ATACBEDDIR}"/"${ATAC['WTXi']}" "${ATACBEDDIR}"/"${ATAC['Patski2-4Xi']}" "${ATACBEDDIR}"/"${ATAC['DelXi']}" "${ATACBEDDIR}"/"${ATAC['InvDxz4Xi']}" "pdfs_100kb" "cumbinnedCounts_ATAC_Xi" "100000"
Rscript "${scriptName}" "${ATACBEDDIR}"/"${ATAC['WTAmbig']}" "${ATACBEDDIR}"/"${ATAC['Patski2-4Ambig']}" "${ATACBEDDIR}"/"${ATAC['DelAmbig']}" "${ATACBEDDIR}"/"${ATAC['InvDxz4Ambig']}" "pdfs_100kb" "cumbinnedCounts_ATAC_Ambig" "100000"


# plots
mkdir pdfs_10kb
# CTCF
Rscript "${scriptName}" "${CTCFBEDDIR}"/"${CTCF['WTXa']}" "${CTCFBEDDIR}"/"${CTCF['Patski2-4Xa']}" "${CTCFBEDDIR}"/"${CTCF['DelXa']}" "${CTCFBEDDIR}"/"${CTCF['InvDxz4Xa']}" "pdfs_10kb" "cumbinnedCounts_CTCF_Xa" "10000"
Rscript "${scriptName}" "${CTCFBEDDIR}"/"${CTCF['WTXi']}" "${CTCFBEDDIR}"/"${CTCF['Patski2-4Xi']}" "${CTCFBEDDIR}"/"${CTCF['DelXi']}" "${CTCFBEDDIR}"/"${CTCF['InvDxz4Xi']}" "pdfs_10kb" "cumbinnedCounts_CTCF_Xi" "10000"
Rscript "${scriptName}" "${CTCFBEDDIR}"/"${CTCF['WTAmbig']}" "${CTCFBEDDIR}"/"${CTCF['Patski2-4Ambig']}" "${CTCFBEDDIR}"/"${CTCF['DelAmbig']}" "${CTCFBEDDIR}"/"${CTCF['InvDxz4Ambig']}" "pdfs_10kb" "cumbinnedCounts_CTCF_Ambig" "10000"
# ATAC
Rscript "${scriptName}" "${ATACBEDDIR}"/"${ATAC['WTXa']}" "${ATACBEDDIR}"/"${ATAC['Patski2-4Xa']}" "${ATACBEDDIR}"/"${ATAC['DelXa']}" "${ATACBEDDIR}"/"${ATAC['InvDxz4Xa']}" "pdfs_10kb" "cumbinnedCounts_ATAC_Xa" "10000"
Rscript "${scriptName}" "${ATACBEDDIR}"/"${ATAC['WTXi']}" "${ATACBEDDIR}"/"${ATAC['Patski2-4Xi']}" "${ATACBEDDIR}"/"${ATAC['DelXi']}" "${ATACBEDDIR}"/"${ATAC['InvDxz4Xi']}" "pdfs_10kb" "cumbinnedCounts_ATAC_Xi" "10000"
Rscript "${scriptName}" "${ATACBEDDIR}"/"${ATAC['WTAmbig']}" "${ATACBEDDIR}"/"${ATAC['Patski2-4Ambig']}" "${ATACBEDDIR}"/"${ATAC['DelAmbig']}" "${ATACBEDDIR}"/"${ATAC['InvDxz4Ambig']}" "pdfs_10kb" "cumbinnedCounts_ATAC_Ambig" "10000"



























##################################################
# Overlaid binned d-scores along the X chromosome
#
# Take average d-score within a bin.
# Allows one to correlate.
#
# Dxz4=chrX:75637518-75764754
# Firre=chrX:50563120-50635321
# Xist=chrX:103460373-103483233 
# chrX	171031299

mkdir "${WORKDIR}"/binnedDscoresAlongX
cd "${WORKDIR}"/binnedDscoresAlongX

# R script
# https://stackoverflow.com/questions/2151212/how-can-i-read-command-line-parameters-from-an-r-script#2154190
scriptName="plotBinnedDscoresAlongChrX.R"
cat << EOF > "${scriptName}"
options(echo=TRUE) # if you want see commands in output file

args <- commandArgs(trailingOnly = TRUE)
print(args)
# trailingOnly=TRUE means that only your arguments are returned, check:
# print(commandArgs(trailingOnly=FALSE))

WTfileName <- as.character(args[1])
Patski24FileName <- as.character(args[2])
DelFileName <- as.character(args[3])
InvDxz4FileName <- as.character(args[4])
figFolder <- as.character(args[5])
dataDesc <- as.character(args[6])
binSize <- as.numeric(as.character(args[7]))
loessspan <- as.numeric(as.character(args[8]))
rm(args)

# Load data
peakData = list()
peakData[["WT"]] = as.data.frame(read.table(gzfile(WTfileName), header=FALSE, skip=1, sep="\t")) #, fill = TRUE)
peakData[["Patski24"]] = as.data.frame(read.table(gzfile(Patski24FileName), header=FALSE, skip=1, sep="\t")) #, fill = TRUE)
peakData[["Del"]] = as.data.frame(read.table(gzfile(DelFileName), header=FALSE, skip=1, sep="\t")) #, fill = TRUE)
peakData[["InvDxz4"]] = as.data.frame(read.table(gzfile(InvDxz4FileName), header=FALSE, skip=1, sep="\t")) #, fill = TRUE)

peakCols =  list()
peakCols[["WT"]] = "blue"
peakCols[["Patski24"]] = "green"
peakCols[["Del"]] = "red"
peakCols[["InvDxz4"]] = "black"

firrebin = (50563120+(abs(50563120-50635321)/2)) / binSize
dxz4bin = (75764754+(abs(75637518-75764754)/2)) / binSize
xistbin = (103460373+(abs(103460373-103483233)/2)) / binSize

# Bin D-scores
midPoints = list()
totalCounts = list()
theIntervals = list()
dscores = list()
meanDscores = list()
meanBinnedDscores = list()
nanBoolIdx = list()
meanBinnedDscoresLoess = list()
meanBinnedDscoresLoessPredict = list()
meanBinnedDscoresSmoothed = list()
meanBinnedDscoresSmoothedPse = list()
meanBinnedDscoresSmoothedMse = list()

for (peakType in names(peakData)) {
#print(head(peakData[[peakType]]))

midPoints[[peakType]] = as.numeric(as.character(peakData[[peakType]][,2])) + (as.numeric(as.character(peakData[[peakType]][,3])) - as.numeric(as.character(peakData[[peakType]][,2])))/2
midPoints[[peakType]] = midPoints[[peakType]]/binSize

#print(head(midPoints[[peakType]]))

totalCounts[[peakType]] = length(midPoints[[peakType]])

theIntervals[[peakType]] = cut(as.numeric(as.character(midPoints[[peakType]])), breaks=seq(0,(180 * (1000000/binSize)) + 1,1))

dscores[[peakType]] = as.numeric(as.character(peakData[[peakType]][, 4])) 
meanDscores[[peakType]] = mean(dscores[[peakType]], na.rm=TRUE)
meanBinnedDscores[[peakType]] = tapply(dscores[[peakType]], theIntervals[[peakType]], mean, na.rm=TRUE)

meanBinnedDscoresLoess[[peakType]] = loess( meanBinnedDscores[[peakType]]~seq(1,(180 * (1000000/binSize)) + 1,1), span=loessspan )
meanBinnedDscoresLoessPredict[[peakType]] = predict(meanBinnedDscoresLoess[[peakType]], se=TRUE)

nanBoolIdx[[peakType]] = is.na(meanBinnedDscores[[peakType]])

meanBinnedDscoresSmoothed[[peakType]] = rep(NaN, length(meanBinnedDscores[[peakType]]))
meanBinnedDscoresSmoothed[[peakType]][!nanBoolIdx[[peakType]]] = meanBinnedDscoresLoessPredict[[peakType]]\$fit

meanBinnedDscoresSmoothedPse[[peakType]] = rep(NaN, length(meanBinnedDscores[[peakType]]))
meanBinnedDscoresSmoothedPse[[peakType]][!nanBoolIdx[[peakType]]] = meanBinnedDscoresLoessPredict[[peakType]]\$fit + meanBinnedDscoresLoessPredict[[peakType]]\$se
meanBinnedDscoresSmoothedMse[[peakType]] = rep(NaN, length(meanBinnedDscores[[peakType]]))
meanBinnedDscoresSmoothedMse[[peakType]][!nanBoolIdx[[peakType]]] = meanBinnedDscoresLoessPredict[[peakType]]\$fit - meanBinnedDscoresLoessPredict[[peakType]]\$se
}

# stop()

# Plots:
pdf(file.path(figFolder, paste(dataDesc,".pdf",sep="")), width=15, height=10)
par(mfrow=c(1,1), mar=c(2,4,6,4), oma=c(3,0,2,0))
theXlims = c(-5, 185 * (1000000/binSize))

# Plot WT d-scores
# Plot an empty graph and legend to get the size of the legend
# https://stackoverflow.com/questions/8929663/r-legend-placement-in-a-plot
plot(seq(1,(180 * (1000000/binSize)) + 1,1), meanBinnedDscores[["WT"]], type="n", xaxt="n", yaxt="n",  xlab=NA, ylab=NA, ylim=c(-0.5, 0.5), xlim=theXlims)
my.legend.size = legend("topleft", 
legend=c(
paste("WT (", totalCounts[["WT"]], " peaks; ", round(meanDscores[["WT"]], 2), ")", sep=""),
paste("Patski2-4 (", totalCounts[["Patski24"]], " peaks; ", round(meanDscores[["Patski24"]], 2), ")", sep=""),
paste("Del (", totalCounts[["Del"]], " peaks; ", round(meanDscores[["Del"]], 2), ")", sep=""),
paste("InvDxz4 (", totalCounts[["InvDxz4"]], " peaks; ", round(meanDscores[["InvDxz4"]], 2), ")", sep="")
),
lty=c(1,1,1,1), col=c(peakCols[["WT"]], peakCols[["Patski24"]], peakCols[["Del"]], peakCols[["InvDxz4"]]), bty='n', cex=1,
plot = FALSE)
# Custom ylim. Add the height of legend to upper bound of the range
my.range <- c(-0.5, 0.5)
my.range[2] <- 1.04*(my.range[2]+my.legend.size\$rect\$h)
my.range <- c(-0.5, 0.5)

# Draw the WT plot with custom ylim
par(new=TRUE)
plot(seq(1,(180 * (1000000/binSize)) + 1,1), meanBinnedDscores[["WT"]], type="p", pch=20, cex=0.75, col=rgb(t(col2rgb(peakCols[["WT"]]))/255, alpha=0.5), axes=F, main=NA, xlab=NA, ylab=NA, ylim=my.range, xlim=theXlims)
px = c(seq(1,(180 * (1000000/binSize)) + 1,1), rev(seq(1,(180 * (1000000/binSize)) + 1,1)))
py = c(meanBinnedDscoresSmoothedPse[["WT"]], rev(meanBinnedDscoresSmoothedMse[["WT"]]))
print(length(px))
print(length(py))
polygon(px[!is.na(py)], py[!is.na(py)], col=rgb(t(col2rgb(peakCols[["WT"]]))/255, alpha=0.2), border=NA)
points(seq(1,(180 * (1000000/binSize)) + 1,1), meanBinnedDscores[["WT"]], pch=20, cex=0.75, col=rgb(t(col2rgb(peakCols[["WT"]]))/255, alpha=0.5))
lines(seq(1,(180 * (1000000/binSize)) + 1,1), meanBinnedDscoresSmoothed[["WT"]], col=peakCols[["WT"]], lty="solid", lwd=4)

lines(meanBinnedDscoresSmoothed[["WT"]], col=peakCols[["WT"]], lty="solid", lwd=5)
axis(side=2)
mtext(side=2, line=2.5, text="mean d-score")
axis(1, at=seq(0, 180 * (1000000/binSize), 20 * (1000000/binSize)) + 1, labels=seq(0, 180, 20))
title(main="Mean d-score", line=3)

# Plot other samples
for (peakType in names(meanDscores)) {
if (peakType != "WT") {
px = c(seq(1,(180 * (1000000/binSize)) + 1,1), rev(seq(1,(180 * (1000000/binSize)) + 1,1)))
py = c(meanBinnedDscoresSmoothedPse[[peakType]], rev(meanBinnedDscoresSmoothedMse[[peakType]]))
polygon(px[!is.na(py)], py[!is.na(py)], col=rgb(t(col2rgb(peakCols[[peakType]]))/255, alpha=0.2), border=NA)
points(as.numeric(as.character(meanBinnedDscores[[peakType]])), type="p", pch=20, cex=0.75, col=rgb(t(col2rgb(peakCols[[peakType]]))/255, alpha=0.5))
lines(meanBinnedDscoresSmoothed[[peakType]], col=peakCols[[peakType]], lty="solid", lwd=5)
}
abline(h=meanDscores[[peakType]], lty="dashed", col=peakCols[[peakType]])
}

# Legend
legend("topleft", 
legend=c(
paste("WT (", totalCounts[["WT"]], " peaks; ", round(meanDscores[["WT"]], 2),  " avg. d-score)",  sep=""),
paste("Patski2-4 (", totalCounts[["Patski24"]], " peaks; ", round(meanDscores[["Patski24"]], 2),  " avg. d-score)",  sep=""),
paste("Del (", totalCounts[["Del"]], " peaks; ", round(meanDscores[["Del"]], 2),  " avg. d-score)", sep=""),
paste("InvDxz4 (", totalCounts[["InvDxz4"]], " peaks; ", round(meanDscores[["InvDxz4"]], 2),  " avg. d-score)", sep="")
),
lty=c(1,1,1,1), col=c(peakCols[["WT"]], peakCols[["Patski24"]], peakCols[["Del"]], peakCols[["InvDxz4"]]), bty='n', cex=1)

# Corr
Pcorr = c()
Scorr = c()
for (peakTypeA in names(peakData)) {
PcorrVect = c()
ScorrVect = c()
for (peakTypeB in names(peakData)) {
PcorrVect = c(PcorrVect, round(cor(as.numeric(as.character(meanBinnedDscores[[peakTypeA]])), as.numeric(as.character(meanBinnedDscores[[peakTypeB]])), method="pearson", use="pairwise.complete.obs"),2))
ScorrVect = c(ScorrVect, round(cor(as.numeric(as.character(meanBinnedDscores[[peakTypeA]])), as.numeric(as.character(meanBinnedDscores[[peakTypeB]])), method="spearman", use="pairwise.complete.obs"),2))
}
Pcorr =rbind(Pcorr, PcorrVect)
Scorr =rbind(Scorr, ScorrVect)
}
rownames(Pcorr) = names(peakData)
colnames(Pcorr) = names(peakData)
rownames(Scorr) = names(peakData)
colnames(Scorr) = names(peakData)

print(Scorr)

# add table
require(plotrix)
# https://stackoverflow.com/questions/18663159/conditional-coloring-of-cells-in-table
addtable2plot("bottomleft", table=Scorr, bg="white", title="Spearman (unsmoothed)",
bty="o", display.rownames=TRUE, display.colnames=TRUE,
hlines=TRUE, vlines=TRUE,
xjust=2, yjust=1, cex=1)

# Smoothed Corr
sPcorr = c()
sScorr = c()
for (peakTypeA in names(peakData)) {
sPcorrVect = c()
sScorrVect = c()
for (peakTypeB in names(peakData)) {
sPcorrVect = c(sPcorrVect, round(cor(as.numeric(as.character(meanBinnedDscoresSmoothed[[peakTypeA]])), as.numeric(as.character(meanBinnedDscoresSmoothed[[peakTypeB]])), method="pearson", use="pairwise.complete.obs"),2))
sScorrVect = c(sScorrVect, round(cor(as.numeric(as.character(meanBinnedDscoresSmoothed[[peakTypeA]])), as.numeric(as.character(meanBinnedDscoresSmoothed[[peakTypeB]])), method="spearman", use="pairwise.complete.obs"),2))
}
sPcorr =rbind(sPcorr, sPcorrVect)
sScorr =rbind(sScorr, sScorrVect)
}
rownames(sPcorr) = names(peakData)
colnames(sPcorr) = names(peakData)
rownames(sScorr) = names(peakData)
colnames(sScorr) = names(peakData)

print(Scorr)

# add table
addtable2plot("bottomright", table=sScorr, bg="white", title="Spearman (smoothed)",
bty="o", display.rownames=TRUE, display.colnames=TRUE,
hlines=TRUE, vlines=TRUE,
xjust=2, yjust=1, cex=1)


# LOI
abline(v=firrebin, lty="dashed", lw=0.5, col="black")
abline(v=dxz4bin, lty="dashed", lw=0.5, col="black")
abline(v=xistbin, lty="dashed", lw=0.5, col="black")
axis(3, at=c(firrebin, dxz4bin, xistbin), labels=c("Firre", "Dxz4", "Xist"), tick = FALSE)

title(main=paste(dataDesc, " (", format(binSize, scientific=FALSE), "bp bins)", sep=""), outer=TRUE)
title(xlab="chrX (Mb)", outer=TRUE, line=1)

dev.off()

EOF
chmod 755 "${scriptName}"


# plots
mkdir pdfs_1Mb
# CTCF
Rscript "${scriptName}" "${CTCFDSCOREDIR}"/"${CTCFDSCORE['WT']}" "${CTCFDSCOREDIR}"/"${CTCFDSCORE['Patski2-4']}" "${CTCFDSCOREDIR}"/"${CTCFDSCORE['Del']}" "${CTCFDSCOREDIR}"/"${CTCFDSCORE['InvDxz4']}" "pdfs_1Mb" "binnedDscores_CTCF" "1000000" "0.15"
# ATAC
Rscript "${scriptName}" "${ATACDSCOREDIR}"/"${ATACDSCORE['WT']}" "${ATACDSCOREDIR}"/"${ATACDSCORE['Patski2-4']}" "${ATACDSCOREDIR}"/"${ATACDSCORE['Del']}" "${ATACDSCOREDIR}"/"${ATACDSCORE['InvDxz4']}" "pdfs_1Mb" "binnedDscores_ATAC" "1000000" "0.15"


# plots
mkdir pdfs_100kb
# CTCF
Rscript "${scriptName}" "${CTCFDSCOREDIR}"/"${CTCFDSCORE['WT']}" "${CTCFDSCOREDIR}"/"${CTCFDSCORE['Patski2-4']}"  "${CTCFDSCOREDIR}"/"${CTCFDSCORE['Del']}" "${CTCFDSCOREDIR}"/"${CTCFDSCORE['InvDxz4']}" "pdfs_100kb" "binnedDscores_CTCF" "100000" "0.15"
# ATAC
Rscript "${scriptName}" "${ATACDSCOREDIR}"/"${ATACDSCORE['WT']}" "${ATACDSCOREDIR}"/"${ATACDSCORE['Patski2-4']}" "${ATACDSCOREDIR}"/"${ATACDSCORE['Del']}" "${ATACDSCOREDIR}"/"${ATACDSCORE['InvDxz4']}" "pdfs_100kb" "binnedDscores_ATAC" "100000" "0.15"




















##################################################
# Overlaid unbinned d-scores along the X chromosome
#
# This is directly from the d-scores.
# But b/c the positions of the peak are all unique,
# one cannot correlate the vectors. 
#
# Dxz4=chrX:75637518-75764754
# Firre=chrX:50563120-50635321
# Xist=chrX:103460373-103483233 
# chrX	171031299

mkdir "${WORKDIR}"/dScoresAlongX
cd "${WORKDIR}"/dScoresAlongX

# R script
# https://stackoverflow.com/questions/2151212/how-can-i-read-command-line-parameters-from-an-r-script#2154190
scriptName="plotBinnedDscoresAlongChrX.R"
cat << EOF > "${scriptName}"
options(echo=TRUE) # if you want see commands in output file

args <- commandArgs(trailingOnly = TRUE)
print(args)
# trailingOnly=TRUE means that only your arguments are returned, check:
# print(commandArgs(trailingOnly=FALSE))

WTfileName <- as.character(args[1])
Patski24FileName <- as.character(args[2])
DelFileName <- as.character(args[3])
InvDxz4FileName <- as.character(args[4])
figFolder <- as.character(args[5])
dataDesc <- as.character(args[6])
loessspan <- as.numeric(as.character(args[7]))
rm(args)

# Load data
peakData = list()
peakData[["WT"]] = as.data.frame(read.table(gzfile(WTfileName), header=FALSE, skip=1, sep="\t")) #, fill = TRUE)
peakData[["Patski24"]] = as.data.frame(read.table(gzfile(Patski24FileName), header=FALSE, skip=1, sep="\t")) #, fill = TRUE)
peakData[["Del"]] = as.data.frame(read.table(gzfile(DelFileName), header=FALSE, skip=1, sep="\t")) #, fill = TRUE)
peakData[["InvDxz4"]] = as.data.frame(read.table(gzfile(InvDxz4FileName), header=FALSE, skip=1, sep="\t")) #, fill = TRUE)

peakCols =  list()
peakCols[["WT"]] = "blue"
peakCols[["Patski24"]] = "green"
peakCols[["Del"]] = "red"
peakCols[["InvDxz4"]] = "black"

firrebin = (50563120+(abs(50563120-50635321)/2)) 
dxz4bin = (75764754+(abs(75637518-75764754)/2))
xistbin = (103460373+(abs(103460373-103483233)/2))

# Bin D-scores
midPoints = list()
totalCounts = list()
theIntervals = list()
dscores = list()
meanDscores = list()
nanBoolIdx = list()
dscoresLoess = list()
dscoresLoessPredict = list()
dscoresSmoothed = list()
dscoresSmoothedPse = list()
dscoresSmoothedMse = list()
for (peakType in names(peakData)) {
#print(head(peakData[[peakType]]))

midPoints[[peakType]] = as.numeric(as.character(peakData[[peakType]][,2])) + (as.numeric(as.character(peakData[[peakType]][,3])) - as.numeric(as.character(peakData[[peakType]][,2])))/2

#print(head(midPoints[[peakType]]))

totalCounts[[peakType]] = length(midPoints[[peakType]])

dscores[[peakType]] = as.numeric(as.character(peakData[[peakType]][, 4])) 
meanDscores[[peakType]] = mean(dscores[[peakType]], na.rm=TRUE)

dscoresLoess[[peakType]] = loess( dscores[[peakType]]~midPoints[[peakType]], span=loessspan )
dscoresLoessPredict[[peakType]] = predict(dscoresLoess[[peakType]], se=TRUE)

nanBoolIdx[[peakType]] = is.na(dscoresLoess[[peakType]])

dscoresSmoothed[[peakType]] = rep(NaN, length(dscores[[peakType]]) )
dscoresSmoothed[[peakType]][!nanBoolIdx[[peakType]]] = dscoresLoessPredict[[peakType]]\$fit
dscoresSmoothedPse[[peakType]] = rep(NaN, length(dscores[[peakType]]) )
dscoresSmoothedPse[[peakType]][!nanBoolIdx[[peakType]]] = dscoresLoessPredict[[peakType]]\$fit + dscoresLoessPredict[[peakType]]\$se
dscoresSmoothedMse[[peakType]] = rep(NaN, length(dscores[[peakType]]) )
dscoresSmoothedMse[[peakType]][!nanBoolIdx[[peakType]]] = dscoresLoessPredict[[peakType]]\$fit - dscoresLoessPredict[[peakType]]\$se
}

# stop()

# Plots:
pdf(file.path(figFolder, paste(dataDesc,".pdf",sep="")), width=20, height=10)
par(mfrow=c(1,1), mar=c(2,4,6,4), oma=c(3,0,2,0))
theXlims = c(0, 180 * 1000000)

# Plot WT d-scores
# Plot an empty graph and legend to get the size of the legend
# https://stackoverflow.com/questions/8929663/r-legend-placement-in-a-plot
plot(midPoints[["WT"]], dscores[["WT"]], type="n", xaxt="n", yaxt="n",  xlab=NA, ylab=NA, ylim=c(-0.5, 0.5), xlim=theXlims)
my.legend.size = legend("topleft", 
legend=c(
paste("WT (", totalCounts[["WT"]], " peaks; ", round(meanDscores[["WT"]], 2), ")", sep=""),
paste("Patski2-4 (", totalCounts[["Patski24"]], " peaks; ", round(meanDscores[["Patski24"]], 2), ")", sep=""),
paste("Del (", totalCounts[["Del"]], " peaks; ", round(meanDscores[["Del"]], 2), ")", sep=""),
paste("InvDxz4 (", totalCounts[["InvDxz4"]], " peaks; ", round(meanDscores[["InvDxz4"]], 2), ")", sep="")
),
lty=c(1,1,1,1), col=c(peakCols[["WT"]], peakCols[["Patski24"]], peakCols[["Del"]], peakCols[["InvDxz4"]]), bty='n', cex=1, 
plot = FALSE)
# Custom ylim. Add the height of legend to upper bound of the range
my.range <- c(-0.5, 0.5)
my.range[2] <- 1.04*(my.range[2]+my.legend.size\$rect\$h)
my.range <- c(-0.5, 0.5)

# Draw the WT plot with custom ylim
par(new=TRUE)
plot(midPoints[["WT"]], dscores[["WT"]], type="p", pch=20, cex=0.75, col=rgb(t(col2rgb(peakCols[["WT"]]))/255, alpha=0.5), axes=F, main=NA, xlab=NA, ylab=NA, ylim=my.range, xlim=theXlims)
px = c(midPoints[["WT"]], rev(midPoints[["WT"]]))
py = c(dscoresSmoothedPse[["WT"]], rev(dscoresSmoothedMse[["WT"]]))
polygon(px[!is.na(py)], py[!is.na(py)], col=rgb(t(col2rgb(peakCols[["WT"]]))/255, alpha=0.2), border=NA)
points(midPoints[["WT"]], dscores[["WT"]], pch=20, cex=0.75, col=rgb(t(col2rgb(peakCols[["WT"]]))/255, alpha=0.5))
lines(midPoints[["WT"]], dscoresSmoothed[["WT"]], col=peakCols[["WT"]], lty="solid", lwd=4)
axis(side=2)
mtext(side=2, line=2.5, text="d-scores")
axis(1, at=seq(0, 180, 20)*1000000, labels=seq(0, 180, 20))
title(main="d-scores", line=3)

# Plot other samples
for (peakType in names(dscores)) {
if (peakType != "WT") {
px = c(midPoints[[peakType]], rev(midPoints[[peakType]]))
py = c(dscoresSmoothedPse[[peakType]], rev(dscoresSmoothedMse[[peakType]]))
polygon(px[!is.na(py)], py[!is.na(py)], col=rgb(t(col2rgb(peakCols[[peakType]]))/255, alpha=0.2), border=NA)
points(midPoints[[peakType]], dscores[[peakType]], type="p", pch=20, cex=0.75, col=rgb(t(col2rgb(peakCols[[peakType]]))/255, alpha=0.35))
lines(midPoints[[peakType]], dscoresSmoothed[[peakType]], col=peakCols[[peakType]], lty="solid", lwd=4)
}
abline(h=meanDscores[[peakType]], lty="dashed", col=peakCols[[peakType]])
}

# Legend
legend("bottomright", 
legend=c(
paste("WT (", totalCounts[["WT"]], " peaks; ", round(meanDscores[["WT"]], 2),  " avg. d-score)",  sep=""),
paste("Patski2-4 (", totalCounts[["Patski24"]], " peaks; ", round(meanDscores[["Patski24"]], 2),  " avg. d-score)", sep=""),
paste("Del (", totalCounts[["Del"]], " peaks; ", round(meanDscores[["Del"]], 2),  " avg. d-score)", sep=""),
paste("InvDxz4 (", totalCounts[["InvDxz4"]], " peaks; ", round(meanDscores[["InvDxz4"]], 2),  " avg. d-score)", sep="")
),
lty=c(1,1,1,1), col=c(peakCols[["WT"]], peakCols[["Patski24"]], peakCols[["Del"]], peakCols[["InvDxz4"]]), bty='n', cex=1)

# LOI
abline(v=firrebin, lty="dashed", lw=0.5, col="black")
abline(v=dxz4bin, lty="dashed", lw=0.5, col="black")
abline(v=xistbin, lty="dashed", lw=0.5, col="black")
axis(3, at=c(firrebin, dxz4bin, xistbin), labels=c("Firre", "Dxz4", "Xist"), tick = FALSE)

title(main=paste(dataDesc, sep=""), outer=TRUE)
title(xlab="chrX (Mb)", outer=TRUE, line=1)

dev.off()

EOF
chmod 755 "${scriptName}"


# plots
mkdir pdfs
# CTCF
Rscript "${scriptName}" "${CTCFDSCOREDIR}"/"${CTCFDSCORE['WT']}" "${CTCFDSCOREDIR}"/"${CTCFDSCORE['Patski2-4']}" "${CTCFDSCOREDIR}"/"${CTCFDSCORE['Del']}" "${CTCFDSCOREDIR}"/"${CTCFDSCORE['InvDxz4']}" "pdfs" "dScores_CTCF" "0.2"
# ATAC
Rscript "${scriptName}" "${ATACDSCOREDIR}"/"${ATACDSCORE['WT']}" "${ATACDSCOREDIR}"/"${ATACDSCORE['Patski2-4']}" "${ATACDSCOREDIR}"/"${ATACDSCORE['Del']}" "${ATACDSCOREDIR}"/"${ATACDSCORE['InvDxz4']}" "pdfs" "dScores_ATAC" "0.2"








