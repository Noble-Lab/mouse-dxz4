##################################################
##################################################
#
# 161101
# Begun
#
# 161129
# Refactored
# 
# 170112
# For remapped HiC data
#
# 170119
# Reworked
#
# 170124
# Alternative plots and additions
#
##################################################
##################################################

source ~/.bashrc
source ~/.bash_profile

PROJDIR="2017_HiC_analysis"
WKDIR="20170610_pooledWTDelDxz4Inv_plotAndCompareContacts_vSNPs"

SUBDIR="PatskiWTvsDelvsDxz4InvFixed"


##################################################
##################################################
##################################################
# Create config files

OUTDIR="${HOME}/proj/${PROJDIR}/results/gbonora/${WKDIR}/${SUBDIR}"
mkdir -p "${OUTDIR}"/input
mkdir -p "${OUTDIR}"/configs
mkdir -p "${OUTDIR}"/output 
cd "${OUTDIR}"

TORUN=(
quantnorm=True
plotrawmaps=False
colormapoptions=False
vanillamaps=False
rot45maps=False
diffmaps=True
v4cplots=False
v4cplots4allLoci=False
v4cplots4alllocilines=False
v4cplots4alllocilinesxaxi=False
v4ccormatrixplots4alllocilinesxaxi=False
)

# #DATASETS=(pooledBrain pooledWT pooledDel)
DATASETS=(pooledWT_Xa pooledWT_Xi pooledDel_Xa pooledDel_Xi Dxz4Inv_Xa Dxz4Inv_Xi )
DATACOLS=('blue' 'blue' 'red' 'red' 'black' 'black')
DATALWS=(1 4 1 4 1 4)
DATALSS=('dashed' 'solid' 'dashed' 'solid' 'dashed' 'solid')

GENES=(
Dxz4=chrX:75637518-75764754
Firre=chrX:50563120-50635321
# #x75_5530601H04Rik=chrX:105048700-105083985
# #X56=chrX:56781370-56833369
ICCE_2210013O21Rik=chrX:153723554-153741296
Xist=chrX:103460373-103483233 
# #WTF=chrX:63221825-67171825 # chrX:60475000-64425000
# TUDR=chrX:65000000-65000001 # chrX:60475000-64425000
# TDDR=chrX:85000000-85000001
Ammecr1=chrX:142853473-142966728
Clcn5=chrX:7158411-7319358
Ddx3x=chrX:13281021-13293983
Eda2r=chrX:97333840-97377209
Eif2s3x=chrX:94188708-94212651
Fhl1=chrX:56731760-56793346
Flna=chrX:74223460-74246534
Mid1=chrX:169685198-169990797
Msn=chrX:96096044-96168553
Ngfrap1=chrX:136270252-136271978
Pcdh19=chrX:133582860-133688993
Rhox5=chrX:37754607-37808878
Shroom4=chrX:6400144-6633717
Tspan7=chrX:10485115-10596604
Wbp5=chrX:136245079-136247139
Zfx=chrX:94074630-94123407
Avpr2=chrX:73892102-73894428
Car5b=chrX:163976822-164028010 
)

#CHROMS=(
#$(seq 19)
#X
#)
#CHROMS=( $(seq 21) )
CHROMS=(20)

REZS=(500000 100000 40000)

# ISWINDOWS=(6 6 12) 
# # 7 * 0.5 = 3.5
# # 7 * 0.1 = 0.7
# # 12 * 0.04 = 0.52

(
# For each bin size
REZTALLY=-1
for REZ in ${REZS[@]}; do 
REZTALLY=$(( REZTALLY + 1))
ISWINDOW=${ISWINDOWS[$REZTALLY]}

for CHROM in ${CHROMS[@]}; do 

CFGfile="${OUTDIR}/configs/heatmaps.${REZ}.chr${CHROM}.cfg"

echo "[general]" > "${CFGfile}"
echo "desc = ${SUBDIR}.${REZ}" >> "${CFGfile}"
echo "assembly = mm10" >> "${CFGfile}"
echo "chrom = chr${CHROM}" >> "${CFGfile}"
echo "rez = ${REZ}" >> "${CFGfile}"

# what to run
echo "" >> "${CFGfile}"
echo "[toRun]" >> "${CFGfile}"
for TRLINE in "${TORUN[@]}"; do
echo "${TRLINE}" >> ${CFGfile}
done

# LOI
echo "" >> "${CFGfile}"
echo "[LOI]" >> "${CFGfile}"
for GENE in "${GENES[@]}"; do
echo "${GENE}" >> ${CFGfile}
done

# # CTCF peaks:
# echo "" >> "${CFGfile}"
# echo "[CTCF]" >> "${CFGfile}"
# echo "ctcf_WT_Xa = ${HOME}/proj/2017_ChIPseq_analysis/results/gbonora/20170101_ChIPseq_CTCF_WTdel1/peaksSegregatedPE/peaksSegregatedPE/ChIPseq_Patski_WT_CTCF.alt/ChIPseq_Patski_WT_CTCF.alt_peaks.chrX.bed.gz" >> "${CFGfile}"
# echo "ctcf_WT_Xi = ${HOME}/proj/2017_ChIPseq_analysis/results/gbonora/20170101_ChIPseq_CTCF_WTdel1/peaksSegregatedPE/peaksSegregatedPE/ChIPseq_Patski_WT_CTCF.ref/ChIPseq_Patski_WT_CTCF.ref_peaks.chrX.bed.gz" >> "${CFGfile}"
# echo "ctcf_Del_Xa = ${HOME}/proj/2017_ChIPseq_analysis/results/gbonora/20170101_ChIPseq_CTCF_WTdel1/peaksSegregatedPE/peaksSegregatedPE/ChIPseq_Patski_Del1_CTCF.alt/ChIPseq_Patski_Del1_CTCF.alt_peaks.chrX.bed.gz" >> "${CFGfile}"
# echo "ctcf_Del_Xi = ${HOME}/proj/2017_ChIPseq_analysis/results/gbonora/20170101_ChIPseq_CTCF_WTdel1/peaksSegregatedPE/peaksSegregatedPE/ChIPseq_Patski_Del1_CTCF.ref/ChIPseq_Patski_Del1_CTCF.ref_peaks.chrX.bed.gz" >> "${CFGfile}"
# 
# # RNAseq
# echo "" >> "${CFGfile}"
# echo "[RNAseq]" >> "${CFGfile}"
# echo "TPM_Xa = ${HOME}/proj/2016_CRISPR_RNAseq/results/gbonora/2016-09-20_CRISPR_5azaRNAseq/dataAndAnalysis/expressionAnalysis/expressionTables/altTPM.chrX.tsv.gz" >> "${CFGfile}"
# echo "TPM_Xi = ${HOME}/proj/2016_CRISPR_RNAseq/results/gbonora/2016-09-20_CRISPR_5azaRNAseq/dataAndAnalysis/expressionAnalysis/expressionTables/refTPM.chrX.tsv.gz" >> "${CFGfile}"

#######################
# For each data set

TALLY=-1
for DATASET in ${DATASETS[@]}; do
TALLY=$(( TALLY + 1))
DATACOL=${DATACOLS[${TALLY}]}
DATALW=${DATALWS[${TALLY}]}
DATALS=${DATALSS[${TALLY}]}

echo "${REZ}.chr${CHROM} ${DATASET} ${DATACOL} ${DATALW}"

if [[ ${DATASET} == "pooledWT_Xa" ]]; then 
dsDir=${HOME}/proj/${PROJDIR}/results/gbonora/20170106_PatskiHiC_poolData_mm10pseudoSp/heatmap.vSNPs
dsFile=insituDNaseHiC-WG-patski-pooled
genome=alt
sampleLabel="WT_Xa"

elif [[ ${DATASET} == "pooledWT_Xi" ]]; then 
dsDir=${HOME}/proj/${PROJDIR}/results/gbonora/20170106_PatskiHiC_poolData_mm10pseudoSp/heatmap.vSNPs
dsFile=insituDNaseHiC-WG-patski-pooled
genome=ref
sampleLabel="WT_Xi"

elif [[ ${DATASET} == "pooledDel_Xa" ]]; then
dsDir=${HOME}/proj/${PROJDIR}/results/gbonora/20170106_DelPatskiHiC_poolData_mm10pseudoSp/heatmap.vSNPs
dsFile=insituDNaseHiC-WG-DelDxz4Patski-pooled
genome=alt
sampleLabel="Del_Xa"

elif [[ ${DATASET} == "pooledDel_Xi" ]]; then
dsDir=${HOME}/proj/${PROJDIR}/results/gbonora/20170106_DelPatskiHiC_poolData_mm10pseudoSp/heatmap.vSNPs
dsFile=insituDNaseHiC-WG-DelDxz4Patski-pooled
genome=ref
sampleLabel="Del_Xi"

elif [[ ${DATASET} == "Dxz4Inv_Xa" ]]; then
dsDir=${HOME}/proj/${PROJDIR}/results/gbonora/20170607_Dxz4InvPatskiHiC_remap_bwamem_mm10pseudoSp/heatmap.vSNPs.CIGARsed # <--- *** NOTE ***
dsFile=insituDNaseHiC-WG-InvDxz4Patski
genome=alt
sampleLabel="Dxz4Inv_Xa"

elif [[ ${DATASET} == "Dxz4Inv_Xi" ]]; then
dsDir=${HOME}/proj/${PROJDIR}/results/gbonora/20170607_Dxz4InvPatskiHiC_remap_bwamem_mm10pseudoSp/heatmap.vSNPs.CIGARsed # <--- *** NOTE ***
dsFile=insituDNaseHiC-WG-InvDxz4Patski
genome=ref
sampleLabel="Dxz4Inv_Xi"

fi
echo ${dsDir} ${dsFile} ${genome}

# Data files:
# E.g.
# ~/NobleLab/proj/2016_CRISPR_HiC/results/gbonora/2016-08-05_CRISPR_HiC_del2PlusWT_pooled/heatmap.vSNPs/insituDNaseHiC-WG-patski-Del2PlusWT-pooled.50000
DATADIR="${dsDir}/${dsFile}.${REZ}"
DATAFILE="${REZ}.${genome}.IC.chr${CHROM}.txt.gz"
echo "--- ${DATADIR}/${DATAFILE} ---"
if [[ ! -f "${DATADIR}/${DATAFILE}"  ]]; then
echo "*** FILE MISSING! ***"
continue
fi
# Link data files
LINKDATADIR="${OUTDIR}/input/${DATASET}"
mkdir -p "${LINKDATADIR}"
ln -sf "${DATADIR}/${DATAFILE}" "${LINKDATADIR}"

# Make config entries
# E.g.
# [WTrefIC]
# inputdir = input/WTpooled.IC
# inputfile = 1000000.ref.IC.chr20.txt.gz
# outputdir = WTpooled.IC
# outputlabel = WT_1Mb_Xi
echo "" >> "${CFGfile}"
echo "[${sampleLabel}]" >> "${CFGfile}"
echo "inputdir = ${LINKDATADIR}" >> "${CFGfile}"
echo "inputfile = ${DATAFILE}" >> "${CFGfile}"
echo "outputdir = heatmaps.${REZ}.chr${CHROM}" >> "${CFGfile}" # Local subfolder under  "${SUBDIR}"/output 
echo "outputlabel = ${DATASET}.${REZ}.chr${CHROM}" >> "${CFGfile}"
echo "datacolor = ${DATACOL}" >> "${CFGfile}"
echo "datalinewidth = ${DATALW}" >> "${CFGfile}"
echo "datalinestyle = ${DATALS}" >> "${CFGfile}"
# echo "iswindows = ${ISWINDOW}" >> "${CFGfile}"

done # DATASET



#######################
# For differential plots
DIFFSETA=(pooledWT_Xa pooledDel_Xi Dxz4Inv_Xi  pooledWT_Xi pooledDel_Xi Dxz4Inv_Xi pooledDel_Xa Dxz4Inv_Xa pooledDel_Xi Dxz4Inv_Xi)
DIFFSETB=(pooledWT_Xi pooledWT_Xi  pooledWT_Xi pooledWT_Xa pooledDel_Xa Dxz4Inv_Xa pooledWT_Xa pooledWT_Xa Dxz4Inv_Xi   pooledDel_Xi)
# DIFFSETA=(pooledWT_Xa)
# DIFFSETB=(pooledWT_Xi)

declare -A DS2LABEL
DS2LABEL[pooledWT_Xa]="WT_Xa"
DS2LABEL[pooledWT_Xi]="WT_Xi"
DS2LABEL[pooledDel_Xa]="Del_Xa"
DS2LABEL[pooledDel_Xi]="Del_Xi"
DS2LABEL[Dxz4Inv_Xa]="Dxz4Inv_Xa"
DS2LABEL[Dxz4Inv_Xi]="Dxz4Inv_Xi"

echo "" >> "${CFGfile}"
echo "[diffMaps]" >> "${CFGfile}"

TALLY=-1
for DATASETA in ${DIFFSETA[@]}; do
TALLY=$(( TALLY + 1))
DATASETB=${DIFFSETB[${TALLY}]}
echo "${TALLY} = ${DS2LABEL[${DATASETA}]}<<<--->>>${DS2LABEL[${DATASETB}]}"
echo "${TALLY} = ${DS2LABEL[${DATASETA}]}<<<--->>>${DS2LABEL[${DATASETB}]}" >> "${CFGfile}"
done # DIFFSETS


done # CHROM
done # REZ
)



OUTDIR="${HOME}/proj/${PROJDIR}/results/gbonora/${WKDIR}/${SUBDIR}"
cd "${OUTDIR}"



##################################################
##################################################
##################################################
# Chromsizes

export PATH=$PATH:"${BINDIR}"/UCSC.utilities
assembly=mm10
if [ ! -e $assembly.sizes ]; then
	fetchChromSizes $assembly > $assembly.sizes
fi



##################################################
##################################################
##################################################
# Plot figures

##################################################
# Create cluster jobs

# 160829
SCRIPTDIR="${HOME}/proj/${PROJDIR}/results/gbonora/${WKDIR}"

#CHROMS=(
#$(seq 19)
#X
#)
CHROMS=( 20 )

REZS=(500000 100000 40000)
#REZS=(40000)


# For each bin size
for REZ in ${REZS[@]}; do 
for CHROM in ${CHROMS[@]}; do 

while read configFile; do
printf "\n$configFile\n"

jobName="heatmaps.${REZ}.chr${CHROM}.job"

cat << EOF > "${jobName}"
source ${HOME}/.bashrc
source ${HOME}/.bash_profile
source /etc/profile.d/modules.sh

# module load python/2.7.11
# module load numpy/1.8.1
# module load scipy/0.17.0
# module load scikit-learn/0.16.1

# 161117 Renamed /net/noble/vol2/home/gbonora/.local/lib/python2.7/site-packages for compatibility reasons.
# For things like Pandas, matplotlib
PYPACKAGEPATH="/net/noble/vol2/home/gbonora/.local/lib/python2.7/site-packages.bak"
export PYTHONPATH=\${PYPACKAGEPATH}:$PYTHONPATH

# Plot matrices
python -u "${SCRIPTDIR}"/plotContactMapsV2.py ${configFile} 

EOF
chmod 755 "${jobName}"

done < <( ls -1 "configs/heatmaps.${REZ}.chr${CHROM}.cfg" )
done # CHROM
done # REZ



#############
# Submit jobs
while read jobName; do
HIRES=$( echo ${jobName} | grep "\.40000\." | wc -l )
if [[ "${HIRES}" -gt 0 ]]; then
qsub  -l h_vmem=32G -l h_rt=47:59:59 -cwd -j y "${jobName}"
else
qsub  -l h_vmem=16G -l h_rt=7:59:59 -cwd -j y "${jobName}"
fi
#./"${jobName}"
done < <( ls -1 heatmaps.*.job | grep -v 500000 )
#done < <( ls -1 heatmaps.500000.*.job )
#done < <( ls -1 heatmaps.*.job | grep -v 40000 )


#############
# Run jons locally

while read jobName; do
./"${jobName}" > "${jobName}".oxo 2>&1 &
done < <( ls -1 heatmaps.*.job )
