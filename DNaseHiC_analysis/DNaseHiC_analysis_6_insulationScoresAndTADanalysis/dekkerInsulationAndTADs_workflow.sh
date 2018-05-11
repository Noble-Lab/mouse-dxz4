##################################################
##################################################
#
# 170329
# Crane version of insulation score to call TADs
#
# From the giorgietti paper:
# 
# "Insulation and boundary calculation. TAD structure (insulation/boundaries) was defined via the insulation method as previously described with minor modifications32. 
# The code used to calculate the insulation score is publicly available on Github (matrix2insulation.pl): https://github.com/dekkerlab/giorgettinature-2016. 
# Insulation vectors were detected using the following options: (-is 480000 -ids 320000 -im iqrMean -nt 0 -ss 160000 -yb 1.5 -nt 0 -bmoe 0). 
# The output of the insulation script is a vector of insulation scores, and a list of minima along the insulation vector (inferred as TAD boundaries). 
# The TAD boundaries were not used in this study"
# 
# Here's an ouput of the script that lists all possible parameters:
# "Options:
#     -v         []         FLAG, verbose mode
#     -o         []         prefix for output file(s)
#     --is       [500000]   optional, insulation square size, size of the insulation square in bp
#     --ss       [0]        optional, insulation smooth size, size of insulation vector smoothing vector in bp
#     --ids      [250000]   optional, insulation delta span (window size of insulationd delta - blue line)
#     --im       [mean]     optional, insulation mode (how to aggregrate signal within insulation square), mean,sum,median
#     --nt       [0.1]      optional, noise threshold, minimum depth of valley
#     --bmoe     [3]        optional, boundary margin of error (specified in number of BINS), added to each side of the boundary
#     --yb       [auto]     optional, -yBound - +yBound of insulation plot
#     --bg       []         FLAG, use transparent background of insulation plot
#     --minDist  []         minimum allowed interaction distance, exclude < N distance (in BP)
#     --maxDist  []         maximum allowed interaction distance, exclude > N distance (in BP)
#     --ez       []         FLAG, ignore 0s in all calculations"
# 
##################################################
##################################################

source ~/.bashrc
source ~/.bash_profile

PROJDIR="<projDirName>"
WKDIR="20170610_pooledWTDelDxz4Inv_plotAndCompareContacts_vSNPs"
SUBDIR="PatskiWTvsDelvsDxz4Inv"

SRCDIR="${PROJDIR}/${WKDIR}/${SUBDIR}"
DKKRDIR="${PROJDIR}/${WKDIR}/${SUBDIR}/dekkerIS"
mkdir -p "${DKKRDIR}" 
cd "${DKKRDIR}"




##################################################
##################################################
##################################################
# 170329
# Insulation score code
# *** New version ***
# https://github.com/dekkerlab/cworld-dekker

cd ~/src
git clone https://github.com/dekkerlab/cworld-dekker.git

cd cworld-dekker

perl Build.PL
./Build
./Build install --install_base ~/perl

# Test
perl ~/src/cworld-dekker/scripts/perl/matrix2insulation.pl




##################################################
##################################################
##################################################
# Convert matrix files to Crane format

mkdir "${DKKRDIR}"/input

CHROM=chr20
ALTCHROM=chrX

SAMPLES=($( ls -1 "${SRCDIR}"/input ))
REZS=(500000 40000)

# For each bin size
for REZ in ${REZS[@]}; do 
for SAMPLE in ${SAMPLES[@]}; do

printf "${SAMPLE}.${REZ}.${ALTCHROM}\n"


if [[ ${SAMPLE/_*/} == "pooledBrain" ]]; then

if [[ ${SAMPLE/*_/} == "Xi" ]]; then
GENOME=alt
else
GENOME=ref 
fi

else

if [[ ${SAMPLE/*_/} == "Xi" ]]; then
GENOME=ref
else
GENOME=alt 
fi

fi


zcat "${SRCDIR}"/input/${SAMPLE}/${REZ}.${GENOME}.IC.${CHROM}.txt.gz | \
awk -v sample=${SAMPLE} -v rez=${REZ} -v chr=${ALTCHROM} \
'NR==1{for(i=1;i<=NF;i++){printf("\tbin%s|%s|%s:%d-%d", i, sample, chr, ((i-1)*rez)+1, (i*rez)+1)}; printf("\n")} \
{printf("bin%s|%s|%s:%d-%d\t", NR, sample, chr, ((NR-1)*rez)+1, (NR*rez)+1); print $0}' | \
gzip > "${DKKRDIR}"/input/${SAMPLE}.${REZ}.${ALTCHROM}.IC.txt.gz

done
done


# Check
ls -ltrh "${DKKRDIR}"/input/

zcat "${DKKRDIR}"/input/*.txt.gz | head



##################################################
##################################################
##################################################
# Create cluster jobs
#
# -is 480000 -ids 320000 -im iqrMean -nt 0 -ss 160000 -yb 1.5 -nt 0 -bmoe 0

SCRIPTDIR="${HOME}/src/cworld-dekker/scripts/perl/"
# perl "${SCRIPTDIR}"/matrix2insulation.pl 

CHROM=chr20
ALTCHROM=chrX

REZS=(500000 40000)
#REZS=(40000)

ISWINDOWS=(7 13) 
# 7 * 0.5 = 3.5
# # 7 * 0.1 = 0.7
# 13 * 0.04 = 0.52

# For each bin size
REZTALLY=-1
for REZ in ${REZS[@]}; do 
REZTALLY=$(( REZTALLY + 1))
ISWINDOW=${ISWINDOWS[$REZTALLY]}

ISPARM=$(( ${REZ} * ${ISWINDOW} ))
IDSPARM=$(( (${REZ} * ${ISWINDOW}) / 2 ))
#SSPARM=$(( ${REZ} * 3 ))
SSPARM=$(( (${REZ} * ${ISWINDOW}) / 4 ))
NOISETHRESH=0.01

SAMPLES=($( ls -1 "${DKKRDIR}"/input/*.${REZ}.${ALTCHROM}.IC.txt.gz ))
for SAMPLE in ${SAMPLES[@]}; do

BNSAMPLE=$(basename ${SAMPLE/.txt.gz/})
echo ${BNSAMPLE}

jobName="dekkerIS.${BNSAMPLE}.job"

cat << EOF > "${jobName}"
source ${HOME}/.bashrc
source ${HOME}/.bash_profile
source /etc/profile.d/modules.sh

module load perl/5.14.2
module load bedtools/2.26.0

mkdir -p output/${BNSAMPLE}

perl "${SCRIPTDIR}"/matrix2insulation.pl -i ${SAMPLE} --is ${ISPARM} --ids ${IDSPARM} --ss ${SSPARM} --nt ${NOISETHRESH} --im mean  -o output/${BNSAMPLE}/${BNSAMPLE}

EOF
chmod 755 "${jobName}"

done # SAMPLE
done # REZ


#############
# Submit jobs
while read jobName; do
qsub  -l h_vmem=16G -l h_rt=7:59:59 -cwd -j y "${jobName}"
done < <( ls -1 dekkerIS.*.job )



#############
# Run jons locally

while read jobName; do
./"${jobName}" > "${jobName}".oxo 2>&1 &
done < <( ls -1 insulationScores.*.job )



#############
# Check

ls -ltrh output/*
ls -ltrh output/*/*.log

for REZ in ${REZS[@]}; do 
SAMPLES=($( ls -1 "${DKKRDIR}"/input/*.${REZ}.${ALTCHROM}.IC.txt.gz ))
for SAMPLE in ${SAMPLES[@]}; do
BNSAMPLE=$(basename ${SAMPLE/.txt.gz/})
#echo ${BNSAMPLE}

find output -name "${BNSAMPLE}.log" | while read FILE; do
printf "\n${FILE}\n"
cat ${FILE} | grep -v "^#"
done

done
done




#############
# Fix and Zip bedgraphs

for REZ in ${REZS[@]}; do 
SAMPLES=($( ls -1 "${DKKRDIR}"/input/*.${REZ}.${ALTCHROM}.IC.txt.gz ))
for SAMPLE in ${SAMPLES[@]}; do
BNSAMPLE=$(basename ${SAMPLE/.txt.gz/})

while read FILE; do
echo "${FILE}"  

mv "${FILE}" "${FILE}".txt
gzip -f "${FILE}".txt

mv "${FILE}".boundaries "${FILE}".boundaries.txt
gzip -f "${FILE}".boundaries.txt

gzip -f "${FILE}".boundaries.bed

cat "${FILE}".bedGraph | sed 's/NA/NaN/g' | gzip > "${FILE}".bedGraph.gz 
rm "${FILE}".bedGraph

done < <( ls -1 output/${BNSAMPLE}/${BNSAMPLE}*.insulation )

done
done

ls -ltrh output/*




#############
# Generate TADs

for REZ in ${REZS[@]}; do 
SAMPLES=($( ls -1 "${DKKRDIR}"/input/*.${REZ}.${ALTCHROM}.IC.txt.gz ))
for SAMPLE in ${SAMPLES[@]}; do
BNSAMPLE=$(basename ${SAMPLE/.txt.gz/})

while read FILE; do
echo "${FILE}"  
gunzip "${FILE}"
gunzip "${FILE/.txt.gz/.boundaries.txt.gz}"
perl "${SCRIPTDIR}"/insulation2tads.pl -i "${FILE/.gz/}" -b "${FILE/.txt.gz/.boundaries.txt}" -o "${FILE/.txt.gz}" 
gzip -f "${FILE/.gz/}"
gzip -f "${FILE/.txt.gz/.boundaries.txt}"
gzip -f "${FILE/.txt.gz/}"*.tads.bed
gzip -f "${FILE/.txt.gz/}"*.tads.bedGraph
gzip -f "${FILE/.txt.gz/}"*.tadHeaders.map
done < <( ls -1 output/${BNSAMPLE}/${BNSAMPLE}*.insulation.txt.gz )

done
done

# Check:
for REZ in ${REZS[@]}; do 
SAMPLES=($( ls -1 "${DKKRDIR}"/input/*.${REZ}.${ALTCHROM}.IC.txt.gz ))
for SAMPLE in ${SAMPLES[@]}; do
BNSAMPLE=$(basename ${SAMPLE/.txt.gz/})

find output -name "${BNSAMPLE}*.tads.bed.gz" | while read FILE; do
printf "${FILE}\t"
zcat ${FILE} | wc -l
done

done
done





