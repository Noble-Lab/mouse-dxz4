##################################################
##################################################
##################################################
#
# 160609
# HTseq to get RPKMs and analysis
#
# Must first align data:
# 'RNAseq_dataQCandMapping_workflow.sh'
# And run allele-specific analysis:
# 'RNAseq_allelicExpressionAnalysis_workflow.sh'
#
##################################################

source ~/.bashrc
source ~/.bash_profile

PROJDIR="<projDirName>"
WORKDIR="<workDirName>"
DATADIR="<dataDirName>"
REFDIR="<refdataDirName>"
BINDIR="<binariesDirName"


##################################################
# Working folder
mkdir -p "${PROJDIR}"/"${WORKDIR}"
cd  "${PROJDIR}"/"${WORKDIR}"





#################################################
#################################################
#################################################
# HTSeq tallying
# http://www-huber.embl.de/HTSeq/doc/overview.html

WDDIR="${PROJDIR}"/"${WORKDIR}/HTSeqCounts"
mkdir "${WDDIR}"
cd "${WDDIR}"

DATADIR=$( cat "${PROJDIR}"/"${WORKDIR}"/DATADIR )
INPUTS=$( cat "${PROJDIR}"/"${WORKDIR}"/libIDs ) 
INPUTDIR="${DATADIR}/sortedBAMs"
HTDIR="${DATADIR}/HTSeq"
mkdir -p "${HTDIR}"

GTF="${REFDIR}/iGenomes/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf"

PARAMS="--format=bam --order=pos --stranded=no --minaqual=0 --type=exon --mode=union --idattr=gene_id"


#################################################
# Generate jobs

for INPUT in ${INPUTS[@]}; do
BASEINPUT=$( basename "${INPUT}" ".fq.gz" )
echo "${INPUT}" "${BASEINPUT}"

jobFile="${BASEINPUT}.HTSeq.job"
cat << EOF > "${jobFile}"
#!/bin/bash -x
#\$ -l h_rt=7:59:59
#\$ -l mfree=4G
#\$ -cwd
#\$ -j y

source /etc/profile.d/modules.sh
module load modules modules-init modules-gs modules-noble
module load python/2.7.3

hostname

printf "\n\nstart: %s\n\n" "\$( date )"

${BINDIR}/htseq-count ${PARAMS} "${INPUTDIR}/${BASEINPUT}".sorted.bam "${GTF}" > "${HTDIR}"/"${BASEINPUT}".HTseq.txt 2> "${BASEINPUT}".HTseq.err

printf "\n\nfinish: %s\n" "\$( date )"

EOF
chmod 755 "${jobFile}"
done


#################################################
# Submit jobs
#sshgrid
for INPUT in ${INPUTS[@]}; do
echo "${INPUT}"
BASEINPUT=$( basename "${INPUT/.fq.gz/}" )

jobFile="${BASEINPUT}.HTSeq.job"
qsub "${jobFile}"
done 


#################################################                                                                                                        
# Check results

ls -ltrh "${HTDIR}"



#################################################                                                                                                        
# Tabulate counts

MASTERFILE="${HTDIR}"/"Patski.AllData.ReadCounts.tsv"

TMPDIR=$(mktemp -dt "GBFILE.XXXXXXX")

###
FILECOUNT=0
PASTECOMMAND="paste"
HEADER="geneID"
for INPUT in ${INPUTS[@]}; do
echo "${INPUT}"
BASEINPUT=$( basename "${INPUT/.fq.gz/}" )

FILECOUNT=$(( FILECOUNT + 1 ))

# Get geneIDs from one of the source files
if [[ "${FILECOUNT}" == 1 ]]; then
# Make sure to zero base data to be consistent with previous 
cat "${HTDIR}"/"${BASEINPUT}".HTseq.txt > "${TMPDIR}"/"${BASEINPUT}"_temptable.tsv
else
# Otherwise just paste the relevant column
cat "${HTDIR}"/"${BASEINPUT}".HTseq.txt | cut -f2 > "${TMPDIR}"/"${BASEINPUT}"_temptable.tsv 
fi
#
PASTECOMMAND="${PASTECOMMAND} ${TMPDIR}/${BASEINPUT}_temptable.tsv"
HEADER="${HEADER}\t${BASEINPUT}"

done

echo "${FILECOUNT} data sets included"
echo "${PASTECOMMAND}"
printf "${HEADER}\n"
###

eval "${PASTECOMMAND}" > "${TMPDIR}"/temptable.tsv
printf "${HEADER}\n" | cat - "${TMPDIR}"/temptable.tsv | awk 'BEGIN{IFS="\t"; OFS="\t"} NR!=2{print $0}' > "${MASTERFILE}"

gzip  -f "${MASTERFILE}"

rm -rf "${TMPDIR}"



zcat -f "${MASTERFILE}" | head 














#################################################
#################################################
#################################################
# Aggregate with previous RNA-seq data sets

PREVRNASEQDIR=${DATADIR}/proj/2017_RNAseq_analysis/data/rawReads-20170220_RNAseq_reprocessAllLibs
ls -ltrh "${PREVRNASEQDIR}"

CURRNASEQDIR=${DATADIR}/data/rawReads-20170620_RNAseq_InvDxz4abAndWTab
ls -ltrh "${CURRNASEQDIR}"

FILE1="${PREVRNASEQDIR}"/HTSeq/Patski.AllData.ReadCounts.tsv.gz
FILE2="${CURRNASEQDIR}"/HTSeq/Patski.AllData.ReadCounts.tsv.gz

zcat "${FILE1}" | wc -l
# 24426
zcat "${FILE2}" | wc -l
# 24426

zcat "${FILE1}" | head -n 2 | awk '{print NF}'
# 23
zcat "${FILE2}" | head -n 2 | awk '{print NF}'
# 5

# Gene info columns
zcat "${FILE1}" | head | cut -f 1
zcat "${FILE2}" | head | cut -f 1

# Genes are the same, so just paste together.
paste <( zcat "${FILE1}" ) <( zcat "${FILE2}" | cut -f 2- ) | gzip > "${CURRNASEQDIR}"/HTSeq/Patski.combo.AllData.ReadCounts.tsv.gz

zcat "${CURRNASEQDIR}"/HTSeq/Patski.combo.AllData.ReadCounts.tsv.gz | wc -l
# 24426

zcat "${CURRNASEQDIR}"/HTSeq/Patski.combo.AllData.ReadCounts.tsv.gz | head -n 2 | awk '{print NF}'
# 27
# = 23 + 5 - 1
















#################################################
#################################################
#################################################
# RPKMs and gene expression analysis


##################################################
# Working on cluster node
sshgrid
# setup using above script

# Interactive session with lots of memory
qlogin -l mfree=64G 
# OR
qlogin -l mfree=8G -pe serial 8 -q noble-long.q 
# setup using above script

WDDIR="${PROJDIR}"/"${WORKDIR}/expressionAnalysis"
mkdir "${WDDIR}"
cd "${WDDIR}"

module load R/3.2.1

R


#################################################                                                                                                        
# Gene lengths (for RPKMs)

# R code in 'geneLengths.R'


#################################################                                                                                                        
# RPKMs

# R code in 'HTSeqCounts2RPKMs.R'


#################################################                                                                                                        
# Expression analysis

# R code in 'allelicExpressionAnalysis.R'


#################################################
# chrX-specific expression

cd "${PROJDIR}"/"${WORKDIR}"/expressionAnalysis/expressionTables
cd "${PROJDIR}"/"${WORKDIR}"/expressionAnalysis/expressionTables.withStrand

# while read FILE; do echo "${FILE/.tsv.gz/.chrX.tsv.gz}"; done < <( ls -1 *.tsv.gz )
while read FILE; do 
# FILE="counts.tsv.gz"
zcat "${FILE}" | head -n 1 > "${FILE/.tsv.gz/.chrX.tsv}"; 
zcat "${FILE}" | grep "chrX" >> "${FILE/.tsv.gz/.chrX.tsv}"; 
gzip -f "${FILE/.tsv.gz/.chrX.tsv}"; 
done < <( ls -1 *.tsv.gz | grep -v "chrX" )


#################################################
# Postprocessing of spreadsheets

cd 	/expressionAnalysis/escapeGenes

cat Patski.AllData_geneLists_overlapMatrix_annotated_minXiProp0.01_minXiTPM1_minXiCov1.tsv | \
awk 'BEGIN{FS="\t"} NR==1{for(i=1; i<=NF; i++){if($i~/Xi\+Xa/){foi[i]=1} else {foi[i]=0}}} \
{for(j=1; j<=NF; j++){if(foi[j]==1){printf("%s", $j); if(j<NF){printf("\t")}}}; printf("\n")}' > \
Patski.AllData_geneLists_overlapMatrix_XiProp_minXiProp0.01_minXiTPM1_minXiCov1.tsv

cat Patski.AllData_geneLists_overlapMatrix_annotated_minXiProp0.01_minXiTPM1_minXiCov1.tsv | \
awk 'BEGIN{FS="\t"} NR==1{for(i=1; i<=NF; i++){if($i~/Xi SNPcov/){foi[i]=1} else {foi[i]=0}}} \
{for(j=1; j<=NF; j++){if(foi[j]==1){printf("%s", $j); if(j<NF){printf("\t")}}}; printf("\n")}' > \
Patski.AllData_geneLists_overlapMatrix_XiSNPcov_minXiProp0.01_minXiTPM1_minXiCov1.tsv

cat Patski.AllData_geneLists_overlapMatrix_annotated_minXiProp0.01_minXiTPM1_minXiCov1.tsv | \
awk 'BEGIN{FS="\t"} NR==1{for(i=1; i<=NF; i++){if($i~/XiRPKM/){foi[i]=1} else {foi[i]=0}}} \
{for(j=1; j<=NF; j++){if(foi[j]==1){printf("%s", $j); if(j<NF){printf("\t")}}}; printf("\n")}' > \
Patski.AllData_geneLists_overlapMatrix_XiRPKM_minXiProp0.01_minXiTPM1_minXiCov1.tsv
