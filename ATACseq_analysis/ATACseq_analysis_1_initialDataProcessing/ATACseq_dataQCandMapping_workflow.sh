##################################################
##################################################
##################################################
# 
# 20171005
# Giancarlo
#
##################################################

source ~/.bashrc
source ~/.bash_profile

PROJDIR="<projDirName>"
WORKDIR="<workDirName>"
DATADIR="<dataDirName>"
REFDIR="<refdataDirName>"
BINDIR="<binariesDirName>"


##################################################
# Working folder
mkdir -p "${PROJDIR}"/"${WORKDIR}"
cd  "${PROJDIR}"/"${WORKDIR}"





#################################################
#################################################
#################################################
# Paired end alignments

mkdir "${PROJDIR}"/"${WORKDIR}"/alignmentsPE
cd "${PROJDIR}"/"${WORKDIR}"/alignmentsPE


##################################################
# Run BWA
# NOTE: PAIRED-END READS 
# ***       2 STEPS:
# ***       1) STEP 1: INITIAL BWA-MEM ALIGNMENT AND 
# ***       2) STEP 2: READ FILTERING FOR PRIMARY READS
./runall.sh


#####################################
# For MACs peak calling, need to filter data 
# Take primary, mapped, paired Reads, MAPQ>=10 
# samtools view -F 256 -F 4 -f 2 -q 10 
#
# https://www.biostars.org/p/5500/
# samtools view -F 0x04 -f 0x02 -b in.bam > out.aligned.bam

############
# genJobs.sh

SAMPLES=( $( ls -1 "${DATADIR}"/fullReads ) )
ASSEMBLY="black6"

for SAMPLE in ${SAMPLES[@]}; do 

# SAMPLE="ChIPseq_Patski_WT_CTCF"
printf "\n----------\n${SAMPLE}\t${ASSEMBLY}\n"

INPUTFILE="${DATADIR}"/mappedReadsPE/"${SAMPLE}"."${ASSEMBLY}"/L1.bam
OUTPUTFILE="${DATADIR}"/mappedReadsPE/"${SAMPLE}"."${ASSEMBLY}"/L1.filtered.bam

jobFile="${SAMPLE}.${ASSEMBLY}.filter.job"
cat << EOF > "${jobFile}"
#!/bin/bash -x
#\$ -l mfree=4G
#\$ -cwd
#\$ -j y

source /etc/profile.d/modules.sh
module load modules modules-init modules-gs modules-noble
module load samtools/1.3

hostname

printf "\n\nstart samtools: %s\n\n" "\$( date )"

samtools view -h -F 256 -F 4 -f 2 -q 10 "${INPUTFILE}" | samtools sort -m 3G -n -T "${SAMPLE}"."${ASSEMBLY}".tmp - > "${OUTPUTFILE}"

printf "\n\nend: %s\n\n" "\$( date )"

EOF

done


############
# Submit jobs
#sshgrid
while read JOBFILE; do 
echo  "${JOBFILE}"
qsub "${JOBFILE}"
done  < <(ls -1 *.filter.job)


############
# Check
ls -ltrh "${DATADIR}"/mappedReadsPE/*

module load samtools/1.3

ASSEMBLY="black6"
for SAMPLE in ${SAMPLES[@]}; do 
printf "\n----------\n${SAMPLE}\t${ASSEMBLY}\n"
OUTPUTFILE="${DATADIR}"/mappedReadsPE/"${SAMPLE}"."${ASSEMBLY}"/L1.filtered.bam
samtools view -h "${OUTPUTFILE}" | head -n 25
done 

###############
# Flag stats before and after filtering
module load samtools/1.3

########
# Before
while read BAM; do
#echo "${BAM}"
samtools flagstat "${BAM}" > "${BAM}.flagstat.txt" 
done < <( find "${DATADIR}"/mappedReadsPE -name "L1.bam" ) &

while read SAMF; do
printf "\n\n${SAMF}\n"
cat "${SAMF}"  
done < <( find "${DATADIR}"/mappedReadsPE -name "*.bam.flagstat.txt" | sort ) | grep -v "^0 + 0" > "${WORKDIR}"_flagstats.txt

cat "${WORKDIR}"_flagstats.txt

########
# After

while read BAM; do
#echo "${BAM}"
samtools flagstat "${BAM}" > "${BAM}.flagstat.filtered.txt" 
done < <( find "${DATADIR}"/mappedReadsPE -name "L1.filtered.bam" ) &

while read SAMF; do
printf "\n\n${SAMF}\n"
cat "${SAMF}"  
done < <( find "${DATADIR}"/mappedReadsPE -name "*.bam.flagstat.filtered.txt" | sort ) | grep -v "^0 + 0" > "${WORKDIR}"_flagstats.filtered.txt

cat "${WORKDIR}"_flagstats.filtered.txt





##################################################
##################################################
##################################################
# Peak calls on black (i.e. unsegregated) PE data
# Using MACS2 https://pypi.python.org/pypi/MACS2

mkdir "${PROJDIR}"/"${WORKDIR}"/callPeaksPE
cd "${PROJDIR}"/"${WORKDIR}"/callPeaksPE


############
# genJobs.sh

INPUTDIR="mappedReadsPE"
OUTPUTDIR="peaksPE"

SAMPLES=( $( ls -1 "${DATADIR}"/mappedReadsPE | grep black6 ) )

for SAMPLE in ${SAMPLES[@]}; do 
# SAMPLE="ChIPseq_Patski_WT_CTCF"
printf "\n----------\n${SAMPLE}\n"

# filenames will be prefix with sample names
OUTDIR="${DATADIR}"/"${OUTPUTDIR}"/"${SAMPLE}"
mkdir -p "${OUTDIR}"

IPFILE="${DATADIR}"/"${INPUTDIR}"/"${SAMPLE}"/L1.filtered.bam
#CNTFILE="${DATADIR}"/"${INPUTDIR}"/"${SAMPLE/CTCF/input}"/L1.filtered.bam

#CHROMINFO="${HOME}/refdata/mm9/chromInfo.txt"
#CHROMINFO="${HOME}/refdata/mm10/chromInfo.txt"
GSIZE=1.87e9 # mm10 genome size

jobFile="${SAMPLE}.MACS2.job"
cat << EOF > "${jobFile}"
#!/bin/bash -x
#\$ -l mfree=4G
#\$ -l h_rt=7:59:59
#\$ -cwd
#\$ -j y

source /etc/profile.d/modules.sh
module load modules modules-init modules-gs modules-noble
module load python/2.7.3
module load cython/0.18.0
module load numpy/1.11.0

hostname

printf "\n\nstart peak calls: %s\n\n" "\$( date )"

export PYTHONPATH=/net/noble/vol2/home/gbonora/.local/lib/python2.7/site-packages:\$PYTHONPATH

# macs2 callpeak -t "${IPFILE}" -c "${CNTFILE}" -f BAMPE -g "${GSIZE}" -B --call-summits --outdir "${OUTDIR}" -n "${SAMPLE}" --verbose 3
macs2 callpeak -t "${IPFILE}" -f BAMPE -g "${GSIZE}" -B --call-summits --outdir "${OUTDIR}" -n "${SAMPLE/.black6/}" --verbose 3

printf "\n\nend: %s\n\n" "\$( date )"
EOF

done


############
# Submit jobs
#sshgrid
while read JOBFILE; do 
echo "${JOBFILE}"
qsub "${JOBFILE}"
done  < <(ls -1 *.MACS2.job)


############
# Check
OUTPUTDIR="peaksPE"
ls -ltrh "${DATADIR}"/"${OUTPUTDIR}"/*


############
# Extract chrX data

while read FILE; do
echo ${FILE}
cat "${FILE}" | grep chrX | gzip -f > "${FILE/.bdg/.chrX.bedgraph.gz}"
done < <( find "${DATADIR}"/"${OUTPUTDIR}"/ -name "*.bdg" )

while read FILE; do
echo ${FILE}
cat "${FILE}" | cut -f 1-3 | grep chrX | gzip -f > "${FILE/.narrowPeak/.chrX.bed.gz}"
done < <( find "${DATADIR}"/"${OUTPUTDIR}"/ -name "*.narrowPeak" )

while read FILE; do
echo ${FILE}
cat "${FILE}" | cut -f 1-3 | grep chrX | gzip -f > "${FILE/summits.bed/summits.chrX.bed.gz}"
done < <( find "${DATADIR}"/"${OUTPUTDIR}"/ -name "*summits.bed" )


############
# Generate BED files for peaks

while read FILE; do
cat "${FILE}" | cut -f 1-3 | gzip -f > "${FILE/.narrowPeak/.bed.gz}"
done < <( find "${DATADIR}"/"${OUTPUTDIR}"/ -name "*.narrowPeak" )

############
# Zip summits BED

while read FILE; do
gzip -f "${FILE}"
done < <( find "${DATADIR}"/"${OUTPUTDIR}"/ -name "*summits.bed" )

############
# Zip narrowpeaks file

while read FILE; do
gzip -f "${FILE}"
done < <( find "${DATADIR}"/"${OUTPUTDIR}"/ -name "*.narrowPeak" )

############
# Rename bedgraph pileups

while read FILE; do
mv "${FILE}" "${FILE/.bdg/.bedgraph}"
gzip -f "${FILE/.bdg/.bedgraph}" &
done < <( find "${DATADIR}"/"${OUTPUTDIR}"/ -name "*.bdg" )

############
# Parse xls file

SAMPLES=( $( ls -1 "${DATADIR}"/mappedReadsPE  | grep black6 ) )

(printf "DataSet\tMappedReads\tDedupReads\tTotalPeaks\tchrXPeaks\n"
for SAMPLE in ${SAMPLES[@]}; do 
printf "${SAMPLE}\t"
XLSFILE="${DATADIR}"/"peaksPE/${SAMPLE}/${SAMPLE/.black6/}_peaks.xls"
#printf "${XLSFILE}\n"
cat "${XLSFILE}" | perl -ne 'print "$1\t" while m|# total fragments in treatment: (.*)|g'
cat "${XLSFILE}" | perl -ne 'print "$1\t" while m|# fragments after filtering in treatment: (.*)|g'
cat "${XLSFILE}" | grep -c "^chr[0-9X]" | tr '\n' '\t'
cat "${XLSFILE}" | grep -c "^chrX"
done ) > "${WORKDIR}".MACS2_summary.txt
$BINDIR/rdb2html -noformatline "${WORKDIR}".MACS2_summary.txt > "${WORKDIR}".MACS2_summary.html
