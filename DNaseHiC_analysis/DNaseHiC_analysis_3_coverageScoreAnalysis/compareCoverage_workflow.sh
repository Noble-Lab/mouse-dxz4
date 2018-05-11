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
# 170214
# Added genes Eda2r, Ddx3x, Zfx
#
##################################################
##################################################

source ~/.bashrc
source ~/.bash_profile

PROJDIR="2017_HiC_analysis"
WKDIR="20170610_pooledWTDelDxz4Inv_plotAndCompareContacts_vSNPs"

SUBDIR="PatskiWTvsDelvsDxz4InvFixed"
usecormatrices=False 

# SUBDIR="PatskiWTvsDelvsDxz4InvFixed_Corr"
# usecormatrices=True


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
usecormatrices="${usecormatrices}"
plotindividuallineplots=False
plotlineplots=False
plotgroupedlineplots=False
plotlineplotzooms=False
genbedgraphs=True
plotdistributions=False
plotpairplots=False
plotsimilaritymaps=False
plotsimilaritymapssansbraindata=False
differentialanalysis=False
differentiallineplots=False
plotdiffplotzooms=False
gendiffbedgraphs=False
plotdiffpairplots=False
delwtdifferentialenrichmentanalysis=False
)

# #DATASETS=(pooledBrain pooledWT pooledDel)
DATASETS=(pooledWT_Xa pooledWT_Xi pooledDel_Xa pooledDel_Xi Dxz4Inv_Xa Dxz4Inv_Xi)
DATACOLS=('blue' 'blue' 'red' 'red'  'black' 'black')
DATALWS=(1 4 1 4 1 4)
DATALSS=('dashed' 'solid' 'dashed' 'solid' 'dashed' 'solid')

GENES=(
Dxz4=chrX:75637518-75764754
Firre=chrX:50563120-50635321
##x75_5530601H04Rik=chrX:105048700-105083985
##X56=chrX:56781370-56833369
#ICCE_2210013O21Rik=chrX:153723554-153741296
Xist=chrX:103460373-103483233 
##WTF=chrX:63221825-67171825 # chrX:60475000-64425000
#TUDR=chrX:65000000-65000001 # chrX:60475000-64425000
#TDDR=chrX:85000000-85000001
#Ammecr1=chrX:142853473-142966728
#Clcn5=chrX:7158411-7319358
#Ddx3x=chrX:13281021-13293983
Eda2r=chrX:97333840-97377209
#Eif2s3x=chrX:94188708-94212651
#Fhl1=chrX:56731760-56793346
#Flna=chrX:74223460-74246534
#Mid1=chrX:169685198-169990797
#Msn=chrX:96096044-96168553
#Ngfrap1=chrX:136270252-136271978
#Pcdh19=chrX:133582860-133688993
#Rhox5=chrX:37754607-37808878
#Shroom4=chrX:6400144-6633717
#Tspan7=chrX:10485115-10596604
#Wbp5=chrX:136245079-136247139
Zfx=chrX:94074630-94123407
)

#CHROMS=(
#$(seq 19)
#X
#)
#CHROMS=( $(seq 21) )
CHROMS=(20)

REZS=(500000 100000 40000)

# ISWINDOWS=(7 7 13) 
# # 7 * 0.5 = 3.5
# # 7 * 0.1 = 0.7
# # 12 * 0.04 = 0.52

# For each bin size
REZTALLY=-1
for REZ in ${REZS[@]}; do 
REZTALLY=$(( REZTALLY + 1))
# ISWINDOW=${ISWINDOWS[$REZTALLY]}

for CHROM in ${CHROMS[@]}; do 

CFGfile="${OUTDIR}/configs/coverageScores.${REZ}.chr${CHROM}.cfg"

echo "[general]" > "${CFGfile}"
echo "desc = ${SUBDIR}.${REZ}" >> "${CFGfile}"
echo "assembly = mm10" >> "${CFGfile}"
echo "chrom = chr${CHROM}" >> "${CFGfile}"
echo "rez = ${REZ}" >> "${CFGfile}"
echo "covscoreoffset = 10" >> "${CFGfile}" # In Mb

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

# CTCF peaks:
echo "" >> "${CFGfile}"
echo "[CTCF]" >> "${CFGfile}"
echo "ctcf_wt_xa = ${HOME}/proj/2017_ChIPseq_analysis/results/gbonora/20170101_ChIPseq_CTCF_WTdel1/callPeaksSegregatedPE/peaksSegregatedPE/ChIPseq_Patski_WT_CTCF.alt/ChIPseq_Patski_WT_CTCF.alt_peaks.chrX.bed.gz" >> "${CFGfile}"
echo "ctcf_wt_xi = ${HOME}/proj/2017_ChIPseq_analysis/results/gbonora/20170101_ChIPseq_CTCF_WTdel1/callPeaksSegregatedPE/peaksSegregatedPE/ChIPseq_Patski_WT_CTCF.ref/ChIPseq_Patski_WT_CTCF.ref_peaks.chrX.bed.gz" >> "${CFGfile}"
echo "ctcf_del_xa = ${HOME}/proj/2017_ChIPseq_analysis/results/gbonora/20170101_ChIPseq_CTCF_WTdel1/callPeaksSegregatedPE/peaksSegregatedPE/ChIPseq_Patski_Del1_CTCF.alt/ChIPseq_Patski_Del1_CTCF.alt_peaks.chrX.bed.gz" >> "${CFGfile}"
echo "ctcf_del_xi = ${HOME}/proj/2017_ChIPseq_analysis/results/gbonora/20170101_ChIPseq_CTCF_WTdel1/callPeaksSegregatedPE/peaksSegregatedPE/ChIPseq_Patski_Del1_CTCF.ref/ChIPseq_Patski_Del1_CTCF.ref_peaks.chrX.bed.gz" >> "${CFGfile}"

# ATAC peaks:
echo "" >> "${CFGfile}"
echo "[ATAC]" >> "${CFGfile}"
echo "atac_wt_xa = ${HOME}/proj/2017_ATACseq_analysis/results/gbonora/20170213_ATACseq_poolData/callPeaksSegregatedPE/peaksSegregatedPE/WT_rep1.alt/WT_rep1.alt_peaks.chrX.bed.gz" >> "${CFGfile}"
echo "atac_wt_xi = ${HOME}/proj/2017_ATACseq_analysis/results/gbonora/20170213_ATACseq_poolData/callPeaksSegregatedPE/peaksSegregatedPE/WT_rep1.ref/WT_rep1.ref_peaks.chrX.bed.gz" >> "${CFGfile}"
echo "atac_del_xa = ${HOME}/proj/2017_ATACseq_analysis/results/gbonora/20170213_ATACseq_poolData/callPeaksSegregatedPE/peaksSegregatedPE/Del1_rep1.alt/Del1_rep1.alt_peaks.chrX.bed.gz" >> "${CFGfile}"
echo "atac_del_xi = ${HOME}/proj/2017_ATACseq_analysis/results/gbonora/20170213_ATACseq_poolData/callPeaksSegregatedPE/peaksSegregatedPE/Del1_rep1.ref/Del1_rep1.ref_peaks.chrX.bed.gz" >> "${CFGfile}"

# RNAseq
echo "" >> "${CFGfile}"
echo "[RNAseq]" >> "${CFGfile}"
echo "tpm_xa = ${HOME}/proj/2016_CRISPR_RNAseq/results/gbonora/2016-09-20_CRISPR_5azaRNAseq/dataAndAnalysis/expressionAnalysis/expressionTables/altTPM.chrX.tsv.gz" >> "${CFGfile}"
echo "tpm_xi = ${HOME}/proj/2016_CRISPR_RNAseq/results/gbonora/2016-09-20_CRISPR_5azaRNAseq/dataAndAnalysis/expressionAnalysis/expressionTables/refTPM.chrX.tsv.gz" >> "${CFGfile}"

# Published TADs:
echo "" >> "${CFGfile}"
echo "[publishedTADs]" >> "${CFGfile}"
echo "dixon_mesc = ${HOME}/refdata/DixonDomains/mESC.combined.domains.mm10.bed.gz" >> "${CFGfile}"
echo "dixon_cortex = ${HOME}/refdata/DixonDomains/cortex.combined.domains.mm10.bed.gz" >> "${CFGfile}"
echo "marks_tsixstop_mesc = ${HOME}/refdata/MarksDomains/Tsix-stop.female.mESC.domains.mm10.bed.gz" >> "${CFGfile}"

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
#echo "outputdir = insulationScores.${REZ}.chr${CHROM}" >> "${CFGfile}" # Local subfolder under  "${SUBDIR}"/output 
#echo "outputlabel = ${DATASET}.${REZ}.chr${CHROM}" >> "${CFGfile}"
echo "datacolor = ${DATACOL}" >> "${CFGfile}"
echo "datalinewidth = ${DATALW}" >> "${CFGfile}"
echo "datalinestyle = ${DATALS}" >> "${CFGfile}"
#echo "iswindows = ${ISWINDOW}" >> "${CFGfile}"

done # DATASET
done # CHROM
done # REZ



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
# Coverage scores

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

jobName="coverageScores.${REZ}.chr${CHROM}.job"

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

## Unlogged count comparison w/o smoothing
#python "${SCRIPTDIR}"/calcCoverageScoresWTvsBrain.py \
#${configFile} \
#"none"

# Unlogged count comparison w/ smoothing
python -u "${SCRIPTDIR}"/calcCoverageScoresV2.py ${configFile} 

EOF
chmod 755 "${jobName}"

done < <( ls -1 "configs/coverageScores.${REZ}.chr${CHROM}.cfg" )
done # CHROM
done # REZ



#############
# Submit jobs
while read jobName; do
HIRES=$( echo ${jobName} | grep "\.40000\." | wc -l )
if [[ "${HIRES}" -gt 0 ]]; then
qsub  -l h_vmem=16G -l h_rt=7:59:59 -cwd -j y "${jobName}"
else
qsub  -l h_vmem=8G -l h_rt=7:59:59 -cwd -j y "${jobName}"
fi
#./"${jobName}"
done < <( ls -1 coverageScores.*.job )



#############
# Run jons locally

while read jobName; do
./"${jobName}" > "${jobName}".oxo 2>&1 &
done < <( ls -1 coverageScores.*.job )









##################################################
##################################################
##################################################
# Download
# Some images

PROJDIR="2017_HiC_analysis"
WKDIR="20170610_pooledWTDelDxz4Inv_plotAndCompareContacts_vSNPs"
SUBDIR="PatskiWTvsDelvsDxz4InvFixed"
#SUBDIR="PatskiWTvsDelvsDxz4InvFixed_Corr"
WDIR="${PROJDIR}/results/gbonora/${WKDIR}/${SUBDIR}"

EXCLUDIR="pooled_ambiguous_40xx000"
EXCLUDIR2="pooled_ambiguous_1000xx00"

# ***NOTE: ORDER MATTERS ***
#rsync -auvPn --include='*nsulationScores*100kb.png' --include='*/' --prune-empty-dirs  --exclude='*'  gbonora@nexus.gs.washington.edu:proj/"${WDIR}"/ ~/NobleLab/proj/"${WDIR}"/
#rsync -auvPn --include='*pm4	0bins*' --include='*/' --prune-empty-dirs  --exclude='*'  gbonora@nexus.gs.washington.edu:proj/"${WDIR}"/ ~/NobleLab/proj/"${WDIR}"/
rsync -auvPn --exclude "${EXCLUDIR}" --exclude "${EXCLUDIR2}" --include='*.bed*' --include='*.png' --include='*.pdf' --include='*/' --prune-empty-dirs  --exclude='*'  gbonora@nexus.gs.washington.edu:proj/"${WDIR}"/ ~/NobleLab/proj/"${WDIR}"/

# Specific
rsync -auvPn  --include='*/corrMatrices/*.pdf' --include='*/' --prune-empty-dirs  --exclude='*'  gbonora@nexus.gs.washington.edu:proj/"${WDIR}"/ ~/NobleLab/proj/"${WDIR}"/
rsync -auvPn  --include='*/distributions/*.pdf' --include='*/' --prune-empty-dirs  --exclude='*'  gbonora@nexus.gs.washington.edu:proj/"${WDIR}"/ ~/NobleLab/proj/"${WDIR}"/
rsync -auvPn  --include='*/linePlotsGrouped/*.pdf' --include='*/' --prune-empty-dirs  --exclude='*'  gbonora@nexus.gs.washington.edu:proj/"${WDIR}"/ ~/NobleLab/proj/"${WDIR}"/
rsync -auvPn  --include='*/coverageScores.chrX/bedgraphs/*.chrX.bedGraph.gz' --include='*/' --prune-empty-dirs  --exclude='*'  gbonora@nexus.gs.washington.edu:proj/"${WDIR}"/ ~/NobleLab/proj/"${WDIR}"/








##################################################
##################################################
##################################################
# HTML

PROJDIR="2017_HiC_analysis"
WKDIR="20170610_pooledWTDelDxz4Inv_plotAndCompareContacts_vSNPs"
SUBDIR="PatskiWTvsDelvsDxz4InvFixed"
#SUBDIR="PatskiWTvsDelvsDxz4InvFixed_Corr"

OUTDIR="${HOME}/proj/${PROJDIR}/results/gbonora/${WKDIR}/${SUBDIR}"
cd "${OUTDIR}"

htmlDir="${HOME}/proj/${PROJDIR}/results/gbonora/${WKDIR}/${SUBDIR}/output/html4CoverageScoreAnalysis"
mkdir "${htmlDir}"

REZS=(500000 100000 40000)

GENES=(
Dxz4=chrX:75637518-75764754
Firre=chrX:50563120-50635321
##x75_5530601H04Rik=chrX:105048700-105083985
##X56=chrX:56781370-56833369
#ICCE_2210013O21Rik=chrX:153723554-153741296
Xist=chrX:103460373-103483233 
##WTF=chrX:63221825-67171825 # chrX:60475000-64425000
#TUDR=chrX:65000000-65000001 # chrX:60475000-64425000
#TDDR=chrX:85000000-85000001
#Ammecr1=chrX:142853473-142966728
#Clcn5=chrX:7158411-7319358
#Ddx3x=chrX:13281021-13293983
Eda2r=chrX:97333840-97377209
#Eif2s3x=chrX:94188708-94212651
#Fhl1=chrX:56731760-56793346
#Flna=chrX:74223460-74246534
#Mid1=chrX:169685198-169990797
#Msn=chrX:96096044-96168553
#Ngfrap1=chrX:136270252-136271978
#Pcdh19=chrX:133582860-133688993
#Rhox5=chrX:37754607-37808878
#Shroom4=chrX:6400144-6633717
#Tspan7=chrX:10485115-10596604
#Wbp5=chrX:136245079-136247139
Zfx=chrX:94074630-94123407
)


################
# Grouped line Plots

plotWhat="linePlotsGrouped"

HTMLfile="${htmlDir}"/coverageScoreAnalysis."${plotWhat}.chrX-wide.html"

echo "<html>" > "${HTMLfile}"

echo "<table>" >> "${HTMLfile}"
echo "<tr><td><h3>chrX-wide ${plotWhat}</h3></td></tr>" >> "${HTMLfile}"
echo "</table>"  >> "${HTMLfile}"

for rez in "${REZS[@]}"; do 
# rez=1000000

echo "<table>" >> "${HTMLfile}"
echo "<tr><td align=\"center\" colspan=2><h1>----------------------------------</h1></td></tr>" >> "${HTMLfile}"
echo "<tr><td align=\"center\" colspan=2><h2>${rez} bp bins</h2></td></tr>" >> "${HTMLfile}"
echo "<tr></tr>" >> "${HTMLfile}"

# Header
echo "<tr>" >> "${HTMLfile}"
for Xtype in Xa Xi; do
echo "<td align=\"center\"><h4>${Xtype}</h4></td>" >> "${HTMLfile}"
done
echo "</tr>" >> "${HTMLfile}"

# PNG
echo "<tr>" >> "${HTMLfile}"
for Xtype in Xa Xi; do
subDir="output/${SUBDIR}.${rez}/coverageScores.chrX"
while read FILE; do
ls -1 ${FILE}
echo "<td align=\"center\"><img width=400 src=\"../../${FILE}\"></img></td>" >> "${HTMLfile}"
done < <( find "${subDir}" -name "coverageZscores.${Xtype}.${SUBDIR}.${rez}.chrX.png" )
done 
echo "</tr>" >> "${HTMLfile}"

# PDF
echo "<tr>" >> "${HTMLfile}"
for Xtype in Xa Xi; do
subDir="output/${SUBDIR}.${rez}/coverageScores.chrX"
while read FILE; do
ls -1 ${FILE}
#echo "<td align=\"center\"><embed src=\"../../${FILE}#page=1&zoom=50\" width=1200 height=550 type='application/pdf'></td>" >> "${HTMLfile}"
echo "<td align=\"center\"><a href=\"../../${FILE}\"><font size="2">PDF</font></a></td>" >> "${HTMLfile}"
done < <( find "${subDir}" -name "coverageZscores.${Xtype}.${SUBDIR}.${rez}.chrX.pdf" )
done 
echo "</tr>" >> "${HTMLfile}"

echo "</table>" >> "${HTMLfile}"

done # rez

echo "</html>" >> "${HTMLfile}"


#--# ################
#--# # line Plots around loci
#--# 
#--# plotWhat="linePlotsGrouped"
#--# 
#--# for geneName in "${GENES[@]}"; do
#--# geneName=$( echo "${geneName}" | cut -d '=' -f 1 | awk '{print tolower($1)}' )
#--# 
#--# 
#--# HTMLfile="${htmlDir}"/coverageScoreAnalysis."${plotWhat}.AroundLocus.${geneName}.html"
#--# 
#--# echo "<html>" > "${HTMLfile}"
#--# 
#--# echo "<table>" >> "${HTMLfile}"
#--# echo "<tr><td><h3>${plotWhat} around ${geneName}</h3></td></tr>" >> "${HTMLfile}"
#--# echo "</table>"  >> "${HTMLfile}"
#--# 
#--# for rez in "${REZS[@]}"; do 
#--# # rez=1000000
#--# 
#--# echo "<table>" >> "${HTMLfile}"
#--# echo "<tr><td align=\"center\" colspan=3><h1>----------------------------------</h1></td></tr>" >> "${HTMLfile}"
#--# echo "<tr><td align=\"center\" colspan=3><h2>${rez} bp bins</h2></td></tr>" >> "${HTMLfile}"
#--# echo "<tr></tr>" >> "${HTMLfile}"
#--# 
#--# # Header
#--# echo "<tr><td></td>" >> "${HTMLfile}"
#--# for Xtype in Xa Xi; do
#--# echo "<td align=\"center\"><h4>${Xtype}</h4></td>" >> "${HTMLfile}"
#--# done
#--# echo "</tr>" >> "${HTMLfile}"
#--# 
#--# # PNG
#--# echo "<tr>" >> "${HTMLfile}"
#--# echo "<td align=\"center\"><h4>${geneName}</h4></td>" >> "${HTMLfile}"
#--# for Xtype in Xa Xi; do
#--# subDir="output/${SUBDIR}.${rez}/coverageScores.chrX"
#--# while read FILE; do
#--# ls -1 ${FILE}
#--# echo "<td align=\"center\"><img width=400 src=\"../../${FILE}\"></img></td>" >> "${HTMLfile}"
#--# done < <( find "${subDir}" -name "coverageZscores.${Xtype}.${geneName}.*.${SUBDIR}.${rez}.chrX.png" )
#--# done 
#--# echo "</tr>" >> "${HTMLfile}"
#--# 
#--# # PDF
#--# echo "<tr>" >> "${HTMLfile}"
#--# echo "<td align=\"center\"><h4>${geneName}</h4></td>" >> "${HTMLfile}"
#--# for Xtype in Xa Xi; do
#--# subDir="output/${SUBDIR}.${rez}/coverageScores.chrX"
#--# while read FILE; do
#--# ls -1 ${FILE}
#--# #echo "<td align=\"center\"><embed src=\"../../${FILE}#page=1&zoom=50\" width=1000 height=550 type='application/pdf'></td>" >> "${HTMLfile}"
#--# echo "<td align=\"center\"><a href=\"../../${FILE}\"><font size="2">PDF</font></a></td>" >> "${HTMLfile}"
#--# done < <( find "${subDir}" -name "coverageZscores.${Xtype}.${geneName}.*.${SUBDIR}.${rez}.chrX.pdf" )
#--# done 
#--# echo "</tr>" >> "${HTMLfile}"
#--# 
#--# echo "</table>" >> "${HTMLfile}"
#--# 
#--# done # rez
#--# 
#--# echo "</html>" >> "${HTMLfile}"
#--# 
#--# done # geneName










################ -------------------------------------------------------------------------------------------------------------- ################
# Often not run.




################
# Differential line Plots

plotWhat="diffCoverageScores"

HTMLfile="${htmlDir}"/coverageScoreAnalysis."${plotWhat}.chrX-wide.html"

echo "<html>" > "${HTMLfile}"

echo "<table>" >> "${HTMLfile}"
echo "<tr><td><h3>chrX-wide ${plotWhat}</h3></td></tr>" >> "${HTMLfile}"
echo "</table>"  >> "${HTMLfile}"

for rez in "${REZS[@]}"; do 
# rez=1000000

echo "<table>" >> "${HTMLfile}"
echo "<tr><td align=\"center\" colspan=2><h1>----------------------------------</h1></td></tr>" >> "${HTMLfile}"
echo "<tr><td align=\"center\" colspan=2><h2>${rez} bp bins</h2></td></tr>" >> "${HTMLfile}"
echo "<tr></tr>" >> "${HTMLfile}"

# Header
echo "<tr>" >> "${HTMLfile}"
for Xtype in Xa Xi; do
echo "<td align=\"center\"><h4>${Xtype}</h4></td>" >> "${HTMLfile}"
done
echo "</tr>" >> "${HTMLfile}"

# PNG
echo "<tr>" >> "${HTMLfile}"
for Xtype in Xa Xi; do
subDir="output/${SUBDIR}.${rez}/coverageScores.chrX"
while read FILE; do
ls -1 ${FILE}
echo "<td align=\"center\"><img width=400 src=\"../../${FILE}\"></img></td>" >> "${HTMLfile}"
done < <( find "${subDir}" -name "${plotWhat}.${Xtype}.${SUBDIR}.${rez}.chrX.png" )
# diffCoverageScores.Xi.PatskiWTvsDelvsDxz4InvFixed.500000.chrX.png
done 
echo "</tr>" >> "${HTMLfile}"

# PDF
echo "<tr>" >> "${HTMLfile}"
for Xtype in Xa Xi; do
subDir="output/${SUBDIR}.${rez}/coverageScores.chrX"
while read FILE; do
ls -1 ${FILE}
#echo "<td align=\"center\"><embed src=\"../../${FILE}#page=1&zoom=50\" width=1200 height=550 type='application/pdf'></td>" >> "${HTMLfile}"
echo "<td align=\"center\"><a href=\"../../${FILE}\"><font size="2">PDF</font></a></td>" >> "${HTMLfile}"
done < <( find "${subDir}" -name "${plotWhat}.${Xtype}.${SUBDIR}.${rez}.chrX.pdf" )
done 
echo "</tr>" >> "${HTMLfile}"

echo "</table>" >> "${HTMLfile}"

done # rez

echo "</html>" >> "${HTMLfile}"





#--# ################
#--# # Differential line Plots around loci
#--# 
#--# plotWhat="diffCoverageScores"
#--# 
#--# for geneName in "${GENES[@]}"; do
#--# geneName=$( echo "${geneName}" | cut -d '=' -f 1 | awk '{print tolower($1)}' )
#--# 
#--# 
#--# HTMLfile="${htmlDir}"/coverageScoreAnalysis."${plotWhat}.AroundLocus.${geneName}.html"
#--# 
#--# echo "<html>" > "${HTMLfile}"
#--# 
#--# echo "<table>" >> "${HTMLfile}"
#--# echo "<tr><td><h3>${plotWhat} around ${geneName}</h3></td></tr>" >> "${HTMLfile}"
#--# echo "</table>"  >> "${HTMLfile}"
#--# 
#--# for rez in "${REZS[@]}"; do 
#--# # rez=1000000
#--# 
#--# echo "<table>" >> "${HTMLfile}"
#--# echo "<tr><td align=\"center\" colspan=3><h1>----------------------------------</h1></td></tr>" >> "${HTMLfile}"
#--# echo "<tr><td align=\"center\" colspan=3><h2>${rez} bp bins</h2></td></tr>" >> "${HTMLfile}"
#--# echo "<tr></tr>" >> "${HTMLfile}"
#--# 
#--# # Header
#--# echo "<tr><td></td>" >> "${HTMLfile}"
#--# for Xtype in Xa Xi; do
#--# echo "<td align=\"center\"><h4>${Xtype}</h4></td>" >> "${HTMLfile}"
#--# done
#--# echo "</tr>" >> "${HTMLfile}"
#--# 
#--# # PNG
#--# echo "<tr>" >> "${HTMLfile}"
#--# echo "<td align=\"center\"><h4>${geneName}</h4></td>" >> "${HTMLfile}"
#--# for Xtype in Xa Xi; do
#--# subDir="output/${SUBDIR}.${rez}/coverageScores.chrX"
#--# while read FILE; do
#--# ls -1 ${FILE}
#--# echo "<td align=\"center\"><img width=400 src=\"../../${FILE}\"></img></td>" >> "${HTMLfile}"
#--# done < <( find "${subDir}" -name "${plotWhat}.${Xtype}.${geneName}.*.${SUBDIR}.${rez}.chrX.png" )
#--# done 
#--# echo "</tr>" >> "${HTMLfile}"
#--# 
#--# # PDF
#--# echo "<tr>" >> "${HTMLfile}"
#--# echo "<td align=\"center\"><h4>${geneName}</h4></td>" >> "${HTMLfile}"
#--# for Xtype in Xa Xi; do
#--# subDir="output/${SUBDIR}.${rez}/coverageScores.chrX"
#--# while read FILE; do
#--# ls -1 ${FILE}
#--# #echo "<td align=\"center\"><embed src=\"../../${FILE}#page=1&zoom=50\" width=1000 height=550 type='application/pdf'></td>" >> "${HTMLfile}"
#--# echo "<td align=\"center\"><a href=\"../../${FILE}\"><font size="2">PDF</font></a></td>" >> "${HTMLfile}"
#--# done < <( find "${subDir}" -name "${plotWhat}.${Xtype}.${geneName}.*.${SUBDIR}.${rez}.chrX.pdf" )
#--# done 
#--# echo "</tr>" >> "${HTMLfile}"
#--# 
#--# echo "</table>" >> "${HTMLfile}"
#--# 
#--# done # rez
#--# 
#--# echo "</html>" >> "${HTMLfile}"
#--# 
#--# done # geneName


################ -------------------------------------------------------------------------------------------------------------- ################









################
# difference Matrices

plotWhat="diffMatrices"

HTMLfile="${htmlDir}"/coverageScoreAnalysis."${plotWhat}.html"

echo "<html>" > "${HTMLfile}"

echo "<table>" >> "${HTMLfile}"
echo "<tr><td><h3>${plotWhat}</h3></td></tr>" >> "${HTMLfile}"
echo "</table>"  >> "${HTMLfile}"

for rez in "${REZS[@]}"; do 
# rez=1000000

echo "<table>" >> "${HTMLfile}"
echo "<tr><td align=\"center\" colspan=5><h1>----------------------------------</h1></td></tr>" >> "${HTMLfile}"
echo "<tr><td align=\"center\" colspan=5><h2>${rez} bp bins</h2></td></tr>" >> "${HTMLfile}"
echo "<tr></tr>" >> "${HTMLfile}"

# Header
echo "<tr>" >> "${HTMLfile}"
for diffType in Euclidean.coverageZscores Pearson.coverageZscores TAD.withNans.adjRI TAD.adjRI TAD.plusPublishedTADs.adjRI; do
echo "<td align=\"center\"><h4>${diffType}</h4></td>" >> "${HTMLfile}"
done
echo "</tr>" >> "${HTMLfile}"

# PNG
echo "<tr>" >> "${HTMLfile}"
for diffType in Euclidean.coverageZscores Pearson.coverageZscores TAD.withNans.adjRI TAD.adjRI TAD.plusPublishedTADs.adjRI; do
subDir="output/${SUBDIR}.${rez}/coverageScores.chrX"
while read FILE; do
ls -1 ${FILE}
echo "<td align=\"center\"><img width=400 src=\"../../${FILE}\"></img></td>" >> "${HTMLfile}"
done < <( find "${subDir}" -name "dendroMatrix.${diffType}.${SUBDIR}.${rez}.chrX.png" )
done 
echo "</tr>" >> "${HTMLfile}"

# PDF
echo "<tr>" >> "${HTMLfile}"
for diffType in Euclidean.coverageZscores Pearson.coverageZscores TAD.withNans.adjRI TAD.adjRI TAD.plusPublishedTADs.adjRI; do
subDir="output/${SUBDIR}.${rez}/coverageScores.chrX"
while read FILE; do
ls -1 ${FILE}
#if [[ "${diffType}" == "TAD.plusPublishedTADs.adjRI" ]]; then
#echo "<td align=\"center\"><embed src=\"../../${FILE}#page=1&zoom=45\" width=550 height=550 type='application/pdf'></td>" >> "${HTMLfile}"
#else
#echo "<td align=\"center\"><embed src=\"../../${FILE}#page=1&zoom=50\" width=550 height=550 type='application/pdf'></td>" >> "${HTMLfile}"
#fi
echo "<td align=\"center\"><a href=\"../../${FILE}\"><font size="2">PDF</font></a></td>" >> "${HTMLfile}"
done < <( find "${subDir}" -name "dendroMatrix.${diffType}.${SUBDIR}.${rez}.chrX.pdf" )
done 
echo "</tr>" >> "${HTMLfile}"

echo "</table>" >> "${HTMLfile}"

done # rez

echo "</html>" >> "${HTMLfile}"




################
# Coverage score violin plots

plotWhat="violinPlot"

for plotType in coverageZscores split.coverageZscores; do

if [[ "${plotType}" == "coverageZscores" ]]; then
WIDTH=950
HEIGHT=500
else
WIDTH=750
HEIGHT=750
fi

HTMLfile="${htmlDir}"/coverageScoreAnalysis."${plotWhat}.${plotType}.html"

echo "<html>" > "${HTMLfile}"

echo "<table>" >> "${HTMLfile}"
echo "<tr><td><h3>${plotWhat}.${plotType}</h3></td></tr>" >> "${HTMLfile}"
echo "</table>"  >> "${HTMLfile}"

for rez in "${REZS[@]}"; do 
# rez=1000000

echo "<table>" >> "${HTMLfile}"
echo "<tr><td align=\"center\" colspan=1><h1>----------------------------------</h1></td></tr>" >> "${HTMLfile}"
echo "<tr><td align=\"center\" colspan=1><h2>${rez} bp bins</h2></td></tr>" >> "${HTMLfile}"
echo "<tr></tr>" >> "${HTMLfile}"

# PNG
echo "<tr>" >> "${HTMLfile}"
subDir="output/${SUBDIR}.${rez}/coverageScores.chrX"
while read FILE; do
ls -1 ${FILE}
echo "<td align=\"center\"><img width=400 src=\"../../${FILE}\"></img></td>" >> "${HTMLfile}"
done < <( find "${subDir}" -name "${plotWhat}.${plotType}.${SUBDIR}.${rez}.chrX.png" )
echo "</tr>" >> "${HTMLfile}"

# PDF
echo "<tr>" >> "${HTMLfile}"
subDir="output/${SUBDIR}.${rez}/coverageScores.chrX"
while read FILE; do
ls -1 ${FILE}
#echo "<td align=\"center\"><embed src=\"../../${FILE}#page=1&zoom=75\" width=${WIDTH} height=${HEIGHT} type='application/pdf'></td>" >> "${HTMLfile}"
echo "<td align=\"center\"><a href=\"../../${FILE}\"><font size="2">PDF</font></a></td>" >> "${HTMLfile}"
done < <( find "${subDir}" -name "${plotWhat}.${plotType}.${SUBDIR}.${rez}.chrX.pdf" )
echo "</tr>" >> "${HTMLfile}"

echo "</table>" >> "${HTMLfile}"

done # rez

echo "</html>" >> "${HTMLfile}"

done # plotType




################
# Coverage pairGrids

plotWhat="pairGrid"

#for plotType in "ByDist" "ByDiffCat"; do
for plotType in "ByDist"; do

HTMLfile="${htmlDir}"/coverageScoreAnalysis."${plotWhat}s.${plotType}.html"

echo "<html>" > "${HTMLfile}"

echo "<table>" >> "${HTMLfile}"
echo "<tr><td><h3>${plotWhat}.${plotType}</h3></td></tr>" >> "${HTMLfile}"
echo "</table>"  >> "${HTMLfile}"

for rez in "${REZS[@]}"; do 
# rez=1000000

echo "<table>" >> "${HTMLfile}"
echo "<tr><td align=\"center\" colspan=2><h1>----------------------------------</h1></td></tr>" >> "${HTMLfile}"
echo "<tr><td align=\"center\" colspan=2><h2>${rez} bp bins</h2></td></tr>" >> "${HTMLfile}"
echo "<tr></tr>" >> "${HTMLfile}"

# Header
echo "<tr>" >> "${HTMLfile}"
for Xtype in "Xa" "Xi"; do
echo "<td align=\"center\" colspan=1><h4>${Xtype}</h4></td>" >> "${HTMLfile}"
done
echo "</tr>" >> "${HTMLfile}"

# PNG
echo "<tr>" >> "${HTMLfile}"
for Xtype in "Xa" "Xi"; do
subDir="output/${SUBDIR}.${rez}/coverageScores.chrX"
while read FILE; do
ls -1 ${FILE}
echo "<td align=\"center\"><img width=400 src=\"../../${FILE}\"></img></td>" >> "${HTMLfile}"
done < <( find "${subDir}" -name "${plotWhat}${plotType}.coverageZscores.${Xtype}.${SUBDIR}.${rez}.chrX.png" )
done # Xtype
echo "</tr>" >> "${HTMLfile}"

# PDF
echo "<tr>" >> "${HTMLfile}"
for Xtype in "Xa" "Xi"; do
subDir="output/${SUBDIR}.${rez}/coverageScores.chrX"
while read FILE; do
ls -1 ${FILE}
#echo "<td align=\"center\"><embed src=\"../../${FILE}#page=1&zoom=75\" width=650 height=650 type='application/pdf'></td>" >> "${HTMLfile}"
echo "<td align=\"center\"><a href=\"../../${FILE}\"><font size="2">PDF</font></a></td>" >> "${HTMLfile}"
done < <( find "${subDir}" -name "${plotWhat}${plotType}.coverageZscores.${Xtype}.${SUBDIR}.${rez}.chrX.pdf" )
done # Xtype
echo "</tr>" >> "${HTMLfile}"

echo "</table>" >> "${HTMLfile}"

done # rez

echo "</html>" >> "${HTMLfile}"

done # plotType




#--# ################
#--# # Barplots of diff regions
#--# 
#--# plotWhat="barplot"
#--# 
#--# for plotType in CTCFpeaks genes; do
#--# 
#--# if [[ "${plotType}" == "CTCFpeaks" ]]; then
#--# WIDTH=600
#--# HEIGHT=650
#--# else
#--# WIDTH=500
#--# HEIGHT=500
#--# fi
#--# 
#--# HTMLfile="${htmlDir}"/coverageScoreAnalysis.differentialXiregions."${plotWhat}.${plotType}.html"
#--# 
#--# echo "<html>" > "${HTMLfile}"
#--# 
#--# echo "<table>" >> "${HTMLfile}"
#--# echo "<tr><td><h3>differentialXiregions.${plotWhat}.${plotType}</h3></td></tr>" >> "${HTMLfile}"
#--# echo "</table>"  >> "${HTMLfile}"
#--# 
#--# for rez in "${REZS[@]}"; do 
#--# # rez=1000000
#--# 
#--# echo "<table>" >> "${HTMLfile}"
#--# echo "<tr><td align=\"center\" colspan=2><h1>----------------------------------</h1></td></tr>" >> "${HTMLfile}"
#--# echo "<tr><td align=\"center\" colspan=2><h2>${rez} bp bins</h2></td></tr>" >> "${HTMLfile}"
#--# echo "<tr></tr>" >> "${HTMLfile}"
#--# 
#--# echo "<tr>" >> "${HTMLfile}"
#--# for plotVersion in "XXX" "OE"; do
#--# echo "<td align=\"center\" colspan=1><h4></h4>${plotVersion/XXX/counts}</td>" >> "${HTMLfile}"
#--# done 
#--# echo "</tr>" >> "${HTMLfile}"
#--# 
#--# # PNG
#--# echo "<tr>" >> "${HTMLfile}"
#--# for plotVersion in "XXX" "OE."; do 
#--# subDir="output/${SUBDIR}.${rez}/coverageScores.chrX"
#--# while read FILE; do
#--# ls -1 ${FILE}
#--# echo "<td align=\"center\"><img width=400 src=\"../../${FILE}\"></img></td>" >> "${HTMLfile}"
#--# done < <( find "${subDir}" -name "${plotWhat}.${plotType}.${plotVersion/XXX/}Del-WTXiDifferentialRegions.${SUBDIR}.${rez}.chrX.png" )
#--# done 
#--# echo "</tr>" >> "${HTMLfile}"
#--# 
#--# # PDF
#--# echo "<tr>" >> "${HTMLfile}"
#--# for plotVersion in "XXX" "OE."; do 
#--# subDir="output/${SUBDIR}.${rez}/coverageScores.chrX"
#--# while read FILE; do
#--# ls -1 ${FILE}
#--# #echo "<td align=\"center\"><embed src=\"../../${FILE}#page=1&zoom=75\" width=${WIDTH} height=${HEIGHT} type='application/pdf'></td>" >> "${HTMLfile}"
#--# echo "<td align=\"center\"><a href=\"../../${FILE}\"><font size="2">PDF</font></a></td>" >> "${HTMLfile}"
#--# done < <( find "${subDir}" -name "${plotWhat}.${plotType}.${plotVersion/XXX/}Del-WTXiDifferentialRegions.${SUBDIR}.${rez}.chrX.pdf" )
#--# done 
#--# echo "</tr>" >> "${HTMLfile}"
#--# 
#--# echo "</table>" >> "${HTMLfile}"
#--# 
#--# done # rez
#--# 
#--# echo "</html>" >> "${HTMLfile}"
#--# 
#--# done # plotType
#--# 
#--# 
#--# 
#--# 
#--# ################
#--# # Expression of diff regions
#--# 
#--# WIDTH=750
#--# HEIGHT=500
#--# 
#--# plotWhat="violinplots.expression"
#--# 
#--# for plotType in "NegDel" "PosDel"; do
#--# 
#--# HTMLfile="${htmlDir}"/coverageScoreAnalysis.differentialXiregions."${plotWhat}.${plotType}.html"
#--# 
#--# echo "<html>" > "${HTMLfile}"
#--# 
#--# echo "<table>" >> "${HTMLfile}"
#--# echo "<tr><td><h3>differentialXiregions.${plotWhat}.${plotType}</h3></td></tr>" >> "${HTMLfile}"
#--# echo "</table>"  >> "${HTMLfile}"
#--# 
#--# for rez in "${REZS[@]}"; do 
#--# # rez=1000000
#--# 
#--# echo "<table>" >> "${HTMLfile}"
#--# echo "<tr><td align=\"center\" colspan=2><h1>----------------------------------</h1></td></tr>" >> "${HTMLfile}"
#--# echo "<tr><td align=\"center\" colspan=2><h2>${rez} bp bins</h2></td></tr>" >> "${HTMLfile}"
#--# echo "<tr></tr>" >> "${HTMLfile}"
#--# 
#--# echo "<tr>" >> "${HTMLfile}"
#--# for Xtype in "xa" "xi"; do
#--# echo "<td align=\"center\" colspan=1><h4></h4>${Xtype}</td>" >> "${HTMLfile}"
#--# done 
#--# echo "</tr>" >> "${HTMLfile}"
#--# 
#--# # PNG
#--# echo "<tr>" >> "${HTMLfile}"
#--# for Xtype in "xa" "xi"; do
#--# subDir="output/${SUBDIR}.${rez}/coverageScores.chrX"
#--# while read FILE; do
#--# ls -1 ${FILE}
#--# echo "<td align=\"center\"><img width=400 src=\"../../${FILE}\"></img></td>" >> "${HTMLfile}"
#--# done < <( find "${subDir}" -name "${plotWhat}.tpm_${Xtype}.${plotType}-WTXiDifferentialRegions.${SUBDIR}.${rez}.chrX.png" )
#--# done 
#--# echo "</tr>" >> "${HTMLfile}"
#--# 
#--# # PDF
#--# echo "<tr>" >> "${HTMLfile}"
#--# for Xtype in "xa" "xi"; do
#--# subDir="output/${SUBDIR}.${rez}/coverageScores.chrX"
#--# while read FILE; do
#--# ls -1 ${FILE}
#--# #echo "<td align=\"center\"><embed src=\"../../${FILE}#page=1&zoom=75\" width=${WIDTH} height=${HEIGHT} type='application/pdf'></td>" >> "${HTMLfile}"
#--# echo "<td align=\"center\"><a href=\"../../${FILE}\"><font size="2">PDF</font></a></td>" >> "${HTMLfile}"
#--# done < <( find "${subDir}" -name "${plotWhat}.tpm_${Xtype}.${plotType}-WTXiDifferentialRegions.${SUBDIR}.${rez}.chrX.pdf" )
#--# done 
#--# echo "</tr>" >> "${HTMLfile}"
#--# 
#--# echo "</table>" >> "${HTMLfile}"
#--# 
#--# done # rez
#--# 
#--# echo "</html>" >> "${HTMLfile}"
#--# 
#--# done # plotType

