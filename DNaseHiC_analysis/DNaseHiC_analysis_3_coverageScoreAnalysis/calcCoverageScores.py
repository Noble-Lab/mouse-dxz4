#############################################################################
#
# Giancarlo Bonora
#
# 170215
# Started by on 'plotContactMaps.py'
#
# 170803
# http://pandas.pydata.org/pandas-docs/version/0.19.2/generated/pandas.DataFrame.icol.html
# DataFrame.icol(i)[source] DEPRECATED. Use .iloc[:, i] instead
# 
# 170804
# XtypeString = '_' + Xtype
# whichData = [XtypeString in i for i in dataSets]
#
#############################################################################

import sys
import os
# import gzip
# import csv
import ConfigParser
import gzip

import numpy as np
import scipy.stats as sps
import scipy.stats.mstats as mstats

# NOT IDEAL
sys.path.insert(0,"$HOME/.local/lib/python2.7/site-packages.bak")
import pandas as pd
print( pd.__version__ )

import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
# matplotlib.style.use('ggplot')
# from matplotlib.colors import LogNorm
import matplotlib.cm

import seaborn as sns

from sklearn.utils.graph import graph_laplacian
from numpy.linalg import eigh
from numpy.linalg import inv

from sklearn.decomposition import PCA

from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist, squareform

# 160829
# For image rotation
from scipy import ndimage
# To resize color bars
from mpl_toolkits.axes_grid1 import make_axes_locatable

# 161130
from matplotlib.font_manager import FontProperties

# 170119
from scipy import signal

# 170124
from itertools import compress
import itertools as it
# For heatmaps
import pylab
import scipy.cluster.hierarchy as sch

# 170201
from sklearn.metrics.cluster import adjusted_rand_score


USAGE = """USAGE: calccoverageScores.py <config file name>
    Calculates coverage scores for contact maps pointed to by config file.

"""


#############################################################################
#############################################################################
##############################################################################
#  FUNCTIONS

##############################################################################
##############################################################################
# Misc functions

##############################################################################
# http://pythoncentral.io/pythons-range-function-explained/
# Note: All arguments are required.
# We're not fancy enough to implement all.
def frange(start, stop, step):
    i = start
    while i < stop:
        yield i
        i += step


##############################################################################
# http://stackoverflow.com/questions/715417/converting-from-a-string-to-boolean-in-python
def str2bool(v):
    return v.lower() in ("yes", "true", "t", "1")


##############################################################################
# http://stackoverflow.com/questions/1094841/reusable-library-to-get-human-readable-version-of-file-size
def sizeof_fmt(num, suffix='b'):
    for unit in ['', 'k', 'M', 'G', 'T', 'P', 'E', 'Z']:
        if abs(num) < 1000.0:
            # return "%3.1f%s%s" % (num, unit, suffix)
            return "%d%s%s" % (num, unit, suffix)
        num /= 1000.0
    # return "%.1f%s%s" % (num, 'Y', suffix)
    return "%d%s%s" % (num, 'Y', suffix)


##############################################################################
##############################################################################
#  Quantile normalization functions

##############################################################################
# 160711
# http://stackoverflow.com/questions/37935920/quantile-normalization-on-pandas-dataframe
# https://github.com/ShawnLYU/Quantile_Normalize/blob/master/quantile_norm.py
def quantileNormalize(df_input):
    df = df_input.copy()
    # compute rank
    dic = {}
    for col in df:
        dic.update({col: sorted(df[col])})
    sorted_df = pd.DataFrame(dic)
    rank = sorted_df.mean(axis=1).tolist()
    # sort
    for col in df:
        t = np.searchsorted(np.sort(df[col]), df[col])
        df[col] = [rank[i] for i in t]
    return df


##############################################################################
def quantileNormalize_2D(arrA, arrB):
    # arrA = dataAlt; arrB = dataRef

    arrAshape = np.shape(arrA)
    arrBshape = np.shape(arrB)

    vecA = arrA.flatten()
    vecB = arrB.flatten()

    dfAB = pd.DataFrame(data=np.column_stack((vecA, vecB)))
    dfAB_QN = quantileNormalize(dfAB)

    vecA_QN = dfAB_QN.iloc[:, 0]
    vecB_QN = dfAB_QN.iloc[:, 1]

    arrA_QN = vecA_QN.reshape(arrAshape)
    arrB_QN = vecB_QN.reshape(arrBshape)

    return [arrA_QN, arrB_QN]


##############################################################################
# 160820
# For 2+ arrays
def quantileNormalize_2Dplus(dataDict):
    # arrA = dataAlt; arrB = dataRef

    tally = -1
    arrShape = {}
    dfArr = []
    for whatData in dataDict:
        # whatData = WTrefIC
        tally += 1

        arr = dataDict[whatData]
        #print arr.shape
        #print np.mean(arr)

        arrShape[whatData] = np.shape(arr)
        vecArr = arr.flatten()

        if tally == 0:
            dfArr = vecArr
        else:
            dfArr = pd.DataFrame(data=np.column_stack((dfArr, vecArr)))

        #print dfArr.shape

    dfArr_QN = quantileNormalize(dfArr)

    dataDict_QN = {}
    tally = -1
    for whatData in dataDict:
        tally += 1
        vecArr_QN = dfArr_QN.iloc[:, tally]
        arr_QN = vecArr_QN.reshape(arrShape[whatData])
        dataDict_QN[whatData] = arr_QN

    return dataDict_QN


##############################################################################
#############################################################################
#############################################################################
# MAIN

##############################################################################
#  Parse the command line.
if len(sys.argv) != 2:
    #sys.stderr.write(str(len(sys.argv)))
    sys.stderr.write(USAGE)
    sys.exit(1)
configFileName = sys.argv[1]


##############################################################################
#  Parse the config file.

print "\n\n\nParse config\n"

# For testing:
# os.getcwd()
# os.chdir('PatskiWTvsDelvsDxz4InvFixed'); configFileName = 'configs/coverageScores.500000.chr20.cfg'

if not os.path.isfile(configFileName):
    sys.stderr.write(configFileName + " does not exist.")
    sys.exit(2)

config = ConfigParser.RawConfigParser()
config.read(configFileName)
config.sections()

print("\n" + config.get('general', 'desc'))
# # List all contents
# print("List all contents")
# for section in config.sections():
#     print("Section: %s" % section)
#     for options in config.options(section):
#         print("x %s:::%s:::%s" % (options,
#                                   config.get(section, options),
#                                   str(type(options))))

sectionTally=0

# -- general  --
projectDesc = config.get('general', 'desc')
assembly = config.get('general', 'assembly')
chrom = config.get('general', 'chrom')
theRez = config.get('general', 'rez')
OffsetMb = int(config.get('general', 'covscoreoffset'))
sectionTally+=1

# -- What to run --
toRun = {}
for option in config.options('toRun'):
    toRun[option] = str2bool(config.get('toRun', option))
sectionTally+=1
# toRun

# -- LOI --
LOI = {}
for option in config.options('LOI'):
    LOI[option] = config.get('LOI', option)
sectionTally+=1
# LOI

# -- CTCF --
CTCF = {}
for option in config.options('CTCF'):
    CTCF[option] = config.get('CTCF', option)
sectionTally+=1
# CTCF

# -- ATAC --
ATAC = {}
for option in config.options('ATAC'):
    ATAC[option] = config.get('ATAC', option)
sectionTally+=1
# ATAC

# -- RNAseq --
RNAseq = {}
for option in config.options('RNAseq'):
    RNAseq[option] = config.get('RNAseq', option)
sectionTally+=1
# RNAseq

# -- publishedTADs --
publishedTADs = {}
for option in config.options('publishedTADs'):
    publishedTADs[option] = config.get('publishedTADs', option)
sectionTally+=1
# publishedTADs

# -- Contact Frequency Data --
dataSets = config.sections()[sectionTally:]
contactMapIOmetadata = {}
for whatData in dataSets:
    contactMapIOmetadata[whatData] = {}
    contactMapIOmetadata[whatData]['inputDir'] = config.get(whatData, 'inputDir')
    contactMapIOmetadata[whatData]['inputFile'] = config.get(whatData, 'inputFile')
    #contactMapIOmetadata[whatData]['outputDir'] = config.get(whatData, 'outputDir')
    #contactMapIOmetadata[whatData]['outputLabel'] = config.get(whatData, 'outputLabel')
    contactMapIOmetadata[whatData]['dataColor'] = config.get(whatData, 'dataColor')
    contactMapIOmetadata[whatData]['dataLineWidth'] = config.get(whatData, 'dataLineWidth')
    contactMapIOmetadata[whatData]['dataLineStyle'] = config.get(whatData, 'dataLineStyle')
    #contactMapIOmetadata[whatData]['isWindows'] = config.get(whatData, 'isWindows')
    sectionTally+=1
# contactMapIOmetadata


################################################
# Other variables

intTheRez = int(theRez)
theRezTxt = sizeof_fmt(intTheRez)
zoomWin = intTheRez * 50

# Contact map smoothing window
smoothingBins = 5

# txformMethod = "none"
# txformMethod = "logged"

# Smooth the contact matrix
# Hard-coded for now:
smoothMethod = "uniform5"
# smoothMethod = "none"
# smoothMethod = "uniform5"
# smoothMethod = "median5"
# smoothMethod = "maximum5"

#theLoi = ['dxz4', 'firre', 'icce_2210013o21rik', 'x75_5530601h04rik', 'xist', 'mecp2', 'kdm5c', 'col4a5']
theLoi = LOI.keys()

# minima filter parameters
yDiffPropThresh = 0.01 # 1% difference
xDiffThresh = 0 # Ignored

# top and bottom 5% for bed files
outThresh = 5

# For plotting & chrom sizes & loi
thechrom = chrom
if assembly[0:2] == 'hg' and chrom == 'chr23':
    thechrom = 'chrX'
elif assembly[0:2] == 'mm' and chrom == 'chr20':
    thechrom = 'chrX'

# output folder
subFolderName =  "coverageScores." + thechrom
outDir = "output/" + projectDesc + "/" + subFolderName
if not os.path.exists(outDir):
    os.makedirs(outDir)

# figure width
# for theRez in (1000000, 500000, 100000, 40000):
#    print 20 + (10 * (500000/intTheRez))
figureWidth = 20 + (10 * (500000/intTheRez))

# # 161130
# # Can make this a parameter in future
# # fileFormat = 'png'
# fileFormat = 'pdf'

# Use TrueType fonts that work in illustrator:
#    http://stackoverflow.com/questions/5956182/cannot-edit-text-in-chart-exported-by-matplotlib-and-opened-in-illustrator
matplotlib.rcParams['pdf.fonttype'] = 42

# Tidy up and save
figTitleSuffix = projectDesc + " " + thechrom
figNameSuffix = projectDesc + "." + thechrom
bedNameSuffix = thechrom


################################################
################################################
################################################
# load chrom sizes
chromSizes = {}
with open(assembly + '.sizes') as f:
  for line in f:
    tok = line.split()
    chromSizes[tok[0]] = tok[1]
#print(chromSizes)


################################################
################################################
################################################
# load CTCF peaks
CTCFpeaks = {}
for CTCFset in CTCF:
    print CTCFset + "; " + CTCF[CTCFset]
    CTCFpeaks[CTCFset] = pd.read_csv(CTCF[CTCFset], compression='gzip', header=None, names=('chr', 'start', 'end'), sep='\t', quotechar='"')


################################################
################################################
################################################
# load ATAC peaks
ATACpeaks = {}
for ATACset in ATAC:
    print ATACset + "; " + ATAC[ATACset]
    ATACpeaks[ATACset] = pd.read_csv(ATAC[ATACset], compression='gzip', header=None, names=('chr', 'start', 'end'), sep='\t', quotechar='"')


################################################
################################################
################################################
# load RNAseq
TPMdata = {}
for RNAset in RNAseq:
    print RNAset + "; " + RNAseq[RNAset]
    TPMdata[RNAset] = pd.read_csv(RNAseq[RNAset], compression='gzip', header=0, sep='\t', quotechar='"')


################################################
################################################
################################################
# Load published TADs

publishedDomains = {}
for dTADset in publishedTADs:
    print dTADset + "; " + publishedTADs[dTADset]
    publishedDomains[dTADset] = pd.read_csv(publishedTADs[dTADset], compression='gzip', header=None, names=('chr', 'start', 'end'), sep='\t', quotechar='"')

publishedDomainsX = {}
for dTADset in publishedDomains:
    publishedDomainsX[dTADset] = publishedDomains[dTADset][publishedDomains[dTADset]['chr'] == "chrX"]



################################################
################################################
################################################
# Bin genomic features

######
# CTCF

CTCFdatasets = sorted(CTCFpeaks.keys())

CTCFpeakDensity = {}
CTCFpeakEnrichment = {}

for ctcfSet in CTCFdatasets:
    # ctcfSet = 'ctcf_wt_xa'

    CTCFpeakDensity[ctcfSet] = [0] * (1 + int(chromSizes[thechrom]) / intTheRez)
    total_peaks = 0
    for index, row in CTCFpeaks[ctcfSet].iterrows():
        # print row
        peakmidpoint = row['start'] + (row['end'] - row['start'])/2
        bin_idx  = peakmidpoint/intTheRez
        CTCFpeakDensity[ctcfSet][bin_idx] += 1
        total_peaks += 1
    CTCFpeakEnrichment[ctcfSet] = np.array(CTCFpeakDensity[ctcfSet]) / np.nanmean(np.array(CTCFpeakDensity[ctcfSet]))


######
# ATAC

ATACdatasets = sorted(ATACpeaks.keys())

ATACpeakDensity = {}
ATACpeakEnrichment = {}

for atacSet in ATACdatasets:
    # atacSet = 'atac_wt_xa'

    ATACpeakDensity[atacSet] = [0] * (1 + int(chromSizes[thechrom]) / intTheRez)
    total_peaks = 0
    for index, row in ATACpeaks[atacSet].iterrows():
        # print row
        peakmidpoint = row['start'] + (row['end'] - row['start'])/2
        bin_idx  = peakmidpoint/intTheRez
        ATACpeakDensity[atacSet][bin_idx] += 1
        total_peaks += 1
    ATACpeakEnrichment[atacSet] = np.array(ATACpeakDensity[atacSet]) / np.nanmean(np.array(ATACpeakDensity[atacSet]))



################################################
################################################
################################################
# load contact maps

print "\n\n\nLoad data\n"

dataDict = {}  # Dictionary of contact data

for whatData in dataSets:
    inputDir = contactMapIOmetadata[whatData]["inputDir"]
    inputFile = contactMapIOmetadata[whatData]["inputFile"]
    #outputDir = contactMapIOmetadata[whatData]["outputDir"]
    #outputLabel = contactMapIOmetadata[whatData]["outputLabel"]

    #print whatData + ": " + inputDir + ", " + inputFile + ", " + outputDir + ", " + outputLabel
    print whatData + ": " + inputDir + ", " + inputFile

    ################################################
    # load data

    # if fileType(inputFile) == "gz":
    #     contactData = gzip.open(inputFile, "r")
    # else:
    #     contactData = open(inputFile, "r")
    # #data = np.loadtxt(inputFile)
    # data = np.array( contactData.read())
    data = np.loadtxt(inputDir + "/" + inputFile)
    #print data.shape
    #print np.mean(data)
    # print np.nanmean(data)

    # Millionify the data
    sum = np.nansum(data)
    mdata = np.array([[z / sum * 1000000 for z in y] for y in data])
    print "millionify:", sum, " --> ", np.nansum(mdata)
    data = mdata
    #print data.shape
    #print np.mean(data)

    # save data
    dataDict[whatData] = data

    # pos_i = np.where(data>0)
    # neg_i = np.where(data<0)

    # # load other stuff
    # chrom_lengths = get_chrom_lengths(chrfile)
    # size = chrom_lengths[chrom]

    # ################################################
    # # Raw data with diags filled with max values
    #
    # # Fill in diag values
    # filldata = np.empty_like(data)
    # filldata[:] = data
    # maxdata = data.max()
    # for i in range(len(data)):
    #     filldata[i, i] = maxdata
    # for i in range(len(data)-1):
    #     filldata[i, i+1] = maxdata
    #     filldata[i+1, i] = maxdata
    #
    # # outFilenameBase =  outDir + "/" + outputLabel + "." + chrom + ".filled"
    # # title = chrom + " " + outputLabel + " filled"
    # # plot_corrMatrices_selectionA(filldata, title, outFilenameBase)



################################################
################################################
################################################
# Quantile normalized across samples

if 'quantnorm' in toRun and toRun['quantnorm']:
    print "\n\n\nQuantile normalization"
    dataDict_QN = quantileNormalize_2Dplus(dataDict)
else:
    dataDict_QN = dataDict



################################################
################################################
################################################
print "\n\n\nClean and smooth contact data\n"
# https://rorasa.wordpress.com/2012/05/13/l0-norm-l1-norm-l2-norm-l-infinity-norm/


################################################
# NaNorator AND DATA SMOOTHERORATOR
# Must remove 0 rows/cols and off/main diags
# Do so across samples, not per sample

dataDict_QN_NaN = {}

# http://stackoverflow.com/questions/15114843/accessing-dictionary-value-by-index-in-python
contactMapSz = dataDict_QN.values()[1].shape[0]
# http://stackoverflow.com/questions/13382774/initialize-list-with-same-bool-value
unMappleRegionBool = [False] * contactMapSz

# Need two passes to get unmappable regions across samples:
# Pass 1: get unmappable regions per samples and aggregate in unMappleRegionBool
for whatData in dataSets:
    #inputDir = contactMapIOmetadata[whatData]["inputDir"]
    #inputFile = contactMapIOmetadata[whatData]["inputFile"]
    #outputDir = contactMapIOmetadata[whatData]["outputDir"]
    #outputLabel = contactMapIOmetadata[whatData]["outputLabel"]

    #print whatData + ": " + inputDir + ", " + inputFile + ", " + outputDir + ", " + outputLabel

    # The QN Data --> FOR SMOOTHING
    data = dataDict_QN[whatData].copy()  # VERY IMPORTANT TO COPY OTHERWISE dataDict gets changed!
    #print data.shape
    print np.mean(data)  # Note: data should have no NaNs so can just use np.mean (a.o.t. np.nanmean)

    # The RAW Data!
    # Need  this to define unmappable regions.
    data_raw = dataDict[whatData].copy()  # VERY IMPORTANT TO COPY OTHERWISE dataDict gets changed!
    #print data_raw.shape
    #print np.mean(data_raw)  # Note: data should have no NaNs so can just use np.mean (a.o.t. np.nanmean)

    # Unmappable regions: rows/cols with all zero counts
    dataDF = pd.DataFrame(data_raw)

    # Per sample, not across samples
    # unMappableRegions = np.where((dataDF.apply(pd.Series.nunique, axis=1) == 1) & (dataDF[1] == 0))
    # data[unMappableRegions, :] = np.nan
    # data[:, unMappableRegions] = np.nan
    # # Diagonal & off diags
    # # data = [[0 for x in range(5)] for x in range(5)] # For testing
    # for i in range(min(len(data[0]), len(data))):
    #     data[i][i] = np.nan
    #     if i > 0:
    #         data[i][i - 1] = np.nan
    #         data[i - 1][i] = np.nan
    # # dataDF = pd.DataFrame(dataDF)
    #
    # dataDict_QN_NaN[whatData] = data

    # For across sample unmappability, keep boolean
    unMappleRegionBool = unMappleRegionBool | ((dataDF.apply(pd.Series.nunique, axis=1) == 1) & (dataDF[1] == 0))

# Indices
unMappableRegions = np.where((unMappleRegionBool))
#print unMappableRegions

# Pass 2: nan unmappable regions given by unMappleRegionBool
#         AND SMOOTH DATA

# Gaussian smoothing affected by nans
# Also, difficult to interpret sigma
# smoothSz = 0.5
# data = ndimage.gaussian_filter(data, sigma=smoothSz, mode='nearest')
if smoothMethod == "uniform5":
    smoothSz = 5
    data = ndimage.uniform_filter(data, size=smoothSz, mode='nearest')
elif smoothMethod == "median5":
    smoothSz = 5
    data = ndimage.median_filter(data, size=smoothSz, mode='nearest')
elif smoothMethod == "maximum5":
    smoothSz = 5
    data = ndimage.maximum_filter(data, size=smoothSz, mode='nearest')

for whatData in dataSets:
    # The QN Data!
    data = dataDict_QN[whatData].copy()  # VERY IMPORTANT TO COPY OTHERWISE dataDict gets changed!

    data[unMappableRegions, :] = np.nan
    data[:, unMappableRegions] = np.nan
    # Diagonal & off diags
    # data = [[0 for x in range(5)] for x in range(5)] # For testing
    for i in range(min(len(data[0]), len(data))):
        data[i][i] = np.nan
        if i > 0:
            data[i][i - 1] = np.nan
            data[i - 1][i] = np.nan
    # dataDF = pd.DataFrame(dataDF)

    # if txformMethod == "logged":
    #     data = np.log([[z + 1 for z in y] for y in data])
    #     #np.nansum(data)
    #     #np.nansum(np.log([[z + 1 for z in y] for y in data]))

    dataDict_QN_NaN[whatData] = data



################################################
################################################
################################################
# 170623
#  Generate correlation matrices


if 'usecormatrices' in toRun and toRun['usecormatrices']:

    print "\n\n\nGenerate correlation matrices\n"

    for whatData in dataSets:
        # The QN Data!
        data = dataDict_QN_NaN[whatData].copy()  # VERY IMPORTANT TO COPY OTHERWISE dataDict gets changed!
        dataDF = pd.DataFrame(data)

        # corData =  np.corrcoef(data)
        # corData =  np.ma.corrcoef(data)
        # corData = pd.DataFrame(data).corr
        corData = pd.DataFrame.corr(dataDF)
        # Rediagonal & off diags
        # data = [[0 for x in range(5)] for x in range(5)] # For testing
        for i in range(min(len(corData[0]), len(corData))):
            corData[i][i] = np.nan
            if i > 0:
                corData[i][i - 1] = np.nan
                corData[i - 1][i] = np.nan

        dataDict_QN_NaN[whatData] = corData.as_matrix()



################################################
################################################
################################################
# Coverage scores

print "\n\n\nCoverage Scores\n"

binOffset = (OffsetMb * 10**6) / intTheRez

covScoresRaw = {}
covScores = {}
covScoresRescaled = {}
covScoresSmoothed = {}
covScoresZscored = {}

for whatDataA in dataDict_QN_NaN.keys():
    #print whatDataA

    # get current contact map
    A = dataDict_QN_NaN[whatDataA].copy()

    Alims = [0, A.shape[0]]

    # From start to finish - minus offset
    covAlims = [0+binOffset, A.shape[0]-binOffset]
    covScore = []

    # pad vector
    for i in range(binOffset):
        covScore.append(np.nan)

    # Bin sliding along the diagonal
    for midBin in range(covAlims[0], covAlims[1]):
        #print midBin
        #print len(covScore)

        # Get coverage array Ids
        # *** NOTE: Ids based on a flattened matrix ***
        # midBin=1
        regionIds = []
        for i in range(0, midBin+1):
            #i = 0
            #i = midBin+1
            # print i
            regionIds.extend( (i*Alims[1]) + np.array(range(midBin, Alims[1])) )
                              # offset            # extent
        # print regionIds
        # print len(regionIds)
        # 6783 ?!
        #  |||
        #  VVV
        # 170621
        # At 500kb, for covAlims == [20, 323], should be 20rows * 323cols = 6460 indices
        # ' for i in range(0, midBin+1):' should be ' for i in range(0, midBin):
        # but former used in all analysis: Unlikey to make much of difference,
        # as each midBin will have the same issue.
        # But strictly off by one.

        # Flatten A
        AF = np.ndarray.flatten(A)

        AFS = AF[regionIds]
        #print(np.nansum(AFS))

        # Remove  unmappable regions
        AnotNanBoolIdx = np.logical_not(np.isnan(AFS))
        AFM = AFS[AnotNanBoolIdx]
        #print(np.nansum(AFM))

        #print(np.sum(AFM)/len(AFM))

        # Coverage score: average count normalized for window size
        # 161129
        # Check % of unmappable regions
        #if len(AFM) > 0:
        # if usuable regions >= 50% then ignore
        if (len(AFM) / float(len(AFS))) >= 0.5:
            covScore.append(np.sum(AFM)/len(AFM))
        else:
            covScore.append(np.nan)

        #print len(covScore)
    # print len(covScore)

    # pad vector
    for i in range(binOffset):
        covScore.append(np.nan)

    # Nan coluns
    AcolSums = np.nansum(A, axis=1)
    ANanBoolIdx = AcolSums==0
    covScore = np.array(covScore)
    covScore[ANanBoolIdx] = np.nan

    # 170621
    # Set zeros to small value to avoid log2(0) == -Inf
    covScore[covScore==0] = 1
    
    # 170711
    covScoresRaw[whatDataA] = covScore
    
    covScores[whatDataA] = np.log2( covScore / np.nanmean(covScore) )

    # # Rescale between 0 and 1... as done in coverage score paper
    # OldMax = np.nanmax(covScore)
    # OldMin = np.nanmin(covScore)
    # NewMax = 1
    # NewMin = 0
    # OldRange = (OldMax - OldMin)
    # NewRange = (NewMax - NewMin)
    # covScoreR = []
    # for OldValue in covScore:
    #     NewValue = (((OldValue - OldMin) * NewRange) / OldRange) + NewMin
    #     covScoreR.append(NewValue)
    # covScoreR = np.array(covScoreR)
    # covScoresRescaled[whatDataA] = covScoreR

    # http://stackoverflow.com/questions/20618804/how-to-smooth-a-curve-in-the-right-way
    # http://scipy.github.io/devdocs/generated/scipy.signal.savgol_filter.html#scipy.signal.savgol_filter
    # window size, polynomial order
    covScoresSmoothed[whatDataA] = signal.savgol_filter(covScores[whatDataA], smoothingBins, 2)

    # Standardize score
    # Doesn't play well with nans
    #  covScoresZscored[whatDataA] = mstats.zscore(covScoresSmoothed[whatDataA])
    covScoresZscored[whatDataA] = (covScoresSmoothed[whatDataA] - np.nanmean(covScoresSmoothed[whatDataA])) \
                                                              / np.nanstd(covScoresSmoothed[whatDataA])

    # Rescale between 0 and 1 after smoothing and z-scoring
    OldMax = np.nanmax(covScoresZscored[whatDataA])
    OldMin = np.nanmin(covScoresZscored[whatDataA])
    NewMax = 1
    NewMin = 0
    OldRange = (OldMax - OldMin)
    NewRange = (NewMax - NewMin)
    covScoreR = []
    for OldValue in covScoresZscored[whatDataA]:
        NewValue = (((OldValue - OldMin) * NewRange) / OldRange) + NewMin
        covScoreR.append(NewValue)
    covScoreR = np.array(covScoreR)
    covScoresRescaled[whatDataA] = covScoreR

# Coverage score loci
covScoreLoci = np.r_[0:len(covScores[whatDataA])]
covScoreLoci = (covScoreLoci * intTheRez) / (1.0*10**6)

# Get union of NaN bins
nandices = []
nonnandices = []
for whatData in covScores:
    nandices.extend(np.argwhere(np.isnan(covScores[whatData])).flatten().tolist())
    nonnandices.extend(np.argwhere(~np.isnan(covScores[whatData])).flatten().tolist())
nandices = np.array(sorted(set(nandices)))
nonnandices = np.array(sorted(set(nonnandices)))

# To denanify coverage scores
allTheNaNsBoolIdx = [False] * covScoresZscored[covScoresZscored.keys()[0]].shape[0]
for dataSet in dataSets:
    allTheNaNsBoolIdx = allTheNaNsBoolIdx | np.isnan(covScoresZscored[dataSet])

# Coverage score loci relative to Dxz4
aloi='dxz4'
chrCoords = LOI[aloi].split(":")
theChr = chrCoords[0]
if theChr == thechrom:
    #print "\tPlot it!"
    startStop = chrCoords[1].split("-")
    startStop = map(int, startStop)
    #startStop = [round(x / float(theRez)) for x in startStop]
    geneMidPoint = (startStop[0] + (startStop[1] - startStop[0]))/10.0**6
    covScoreLociRTD = abs(covScoreLoci - geneMidPoint)
covScoreLociRTDclasses = ['>20Mb'] * len(covScoreLociRTD)
# covScoreLociRTDclasses = ['>20Mb' for x in range(len(covScoreLociRTD))]
for i in np.nditer(np.where(np.array(covScoreLociRTD < 1))):
    #print i
    covScoreLociRTDclasses[i] = '<1Mb'
for i in np.nditer(np.where(np.bitwise_and(np.array(covScoreLociRTD >= 1), np.array(covScoreLociRTD < 5)))):
    covScoreLociRTDclasses[i] = '1-5Mb'
for i in np.nditer(np.where(np.bitwise_and(np.array(covScoreLociRTD >= 5), np.array(covScoreLociRTD < 10)))):
    covScoreLociRTDclasses[i] = '5-10Mb'
for i in np.nditer(np.where(np.bitwise_and(np.array(covScoreLociRTD >= 10), np.array(covScoreLociRTD < 20)))):
    covScoreLociRTDclasses[i] = '10-20Mb'


################################################
################################################
################################################
# TADs based on local minima of Coverage z-scores
# Based on local filter minima
# http://stackoverflow.com/questions/4624970/finding-local-maxima-minima-with-numpy-in-a-1d-numpy-array#9667121

print "\n\n\ncoverage score TADs\n"

TADdir = outDir + "/TADs"
if not os.path.exists(TADdir):
    os.makedirs(TADdir)

# Y-axis limits
# minY = -2
# maxY = 2
minY = 10**6
maxY = 0
for dataSet in dataSets:
    # dataSet = 'pooledWT_alt_IC'
    #print dataSet
    if minY > np.nanmin(covScoresZscored[dataSet]):
        minY = np.nanmin(covScoresZscored[dataSet])
    if maxY < np.nanmax(covScoresZscored[dataSet]):
        maxY = np.nanmax(covScoresZscored[dataSet])
# Hard-coded for now
yDiffThresh = (maxY - minY) * yDiffPropThresh

filteredMin = {}

idxTally = -1
for dataSet in dataSets:
    # dataSet = 'pooledWT_alt_IC'
    #print dataSet

    idxTally += 1

    xVals = covScoreLoci
    yValsS = covScoresZscored[dataSet]

    # that's the line, you need:
    # a = np.diff(np.sign(np.diff(yValsS))).nonzero()[0] + 1 # local min+max -- NANs an issue here
    b = (np.diff(np.sign(np.diff(yValsS))) > 0).nonzero()[0] + 1 # local min
    c = (np.diff(np.sign(np.diff(yValsS))) < 0).nonzero()[0] + 1 # local max

    if len(b) > 0:
        TADfileName = TADdir + "/TADs.ByCoverageScore." \
                      + dataSet + "." + bedNameSuffix + ".bed.gz"
        f = gzip.open(TADfileName, "wb")

        bSmoothed = []
        # old_i = 0
        old_i = np.nanmin(nonnandices)
        for i in np.nditer(b):
            #print "\n---"
            # i = 14
            i = int(i)

            # Check for intervening NaNs
            if len(nandices[np.bitwise_and(nandices<i, nandices>old_i)]) >= 1:
                 minNonNanIdx = nandices[np.bitwise_and(nandices<i, nandices>old_i)][0] - 1
                 maxNonNanIdx = nandices[np.bitwise_and(nandices<i, nandices>old_i)][-1] + 1
                 # print "NANs: old_i = " + str(old_i)  + "; i = " + str(i) + "; minNonNanIdx = " + str(minNonNanIdx) + "; maxNonNanIdx = " + str(maxNonNanIdx)
                 f.write(
                    thechrom + "\t" \
                    + str(int(round(xVals[old_i] * 1.0*10**6))) + '\t' \
                    + str(int(round(xVals[minNonNanIdx] * 1.0*10**6))) + "\n")
                 old_i = maxNonNanIdx

            # Previous max
            # j = c[(c < i) & (c > old_i)]
            # print j
            # if len(j) > 1:
            #     j = int(j[-1])
            #old_i = i
            if len(c[c<i]) >= 1:
                j = c[c<i][-1]
            else:
                j = np.nan
            # Next max
            if len(c[c>i]) >= 1:
                k = c[(c>i)][0]
            else:
                k = np.nan
            # print "min = " + str(i)  + "; prevMax = " + str(j) + "; nexMax = " + str(k)

            # Just compare to previous max.
            # print str((xVals[i] - xVals[old_i]) * 1.0*10**6) + " - " + str(xThresh)
            # print str(yVals[i] - yVals[old_i]) + " - " + str(yThresh)
            # if ( (((xVals[i] - xVals[old_i]) * 1.0*10**6) >= xThresh) \
            #         and (abs(yVals[i] - yVals[old_i]) >= yThresh) )\
            #    or \

            # Compare to previous and next max.
            if not np.isnan(j) and not np.isnan(k):
                xMinMaxDiff = np.nanmin( np.array((xVals[i] - xVals[j]) * 1.0*10**6,
                                                  (xVals[k] - xVals[i]) * 1.0*10**6 ))
                yMinMaxDiff = np.nanmin( np.array(abs(yValsS[i] - yValsS[j]),
                                                  abs(yValsS[k] - yValsS[i]) ))
            elif not np.isnan(j):
                xMinMaxDiff = (xVals[i] - xVals[j]) * 1.0*10**6
                yMinMaxDiff = abs(yValsS[i] - yValsS[j])
            elif not np.isnan(k):
                xMinMaxDiff = (xVals[k] - xVals[i]) * 1.0*10**6
                yMinMaxDiff = abs(yValsS[k] - yValsS[i])
            else:
                xMinMaxDiff = np.inf
                yMinMaxDiff = np.inf

            if ( ( xMinMaxDiff >= xDiffThresh) and (yMinMaxDiff >= yDiffThresh) ):
                bSmoothed.append(i)
                f.write(
                    thechrom + "\t" \
                    + str(int(round(xVals[old_i] * 1.0*10**6))) + '\t' \
                    + str(int(round(xVals[i] * 1.0*10**6))) + "\n")
                old_i = i

        maxNonnanIdx = np.nanmax(nonnandices)
        f.write(
            thechrom + "\t" \
            + str(int(round(xVals[old_i] * 1.0*10**6))) + '\t' \
            #+ str(int(round(xVals[-1] * 1.0*10**6))) + "\n")
            + str(int(round(xVals[maxNonnanIdx] * 1.0*10**6))) + "\n")
        f.close()

        filteredMin[dataSet] = bSmoothed
        #print("**")
        #print(bSmoothed)

#print("***")
#print(filteredMin)



################################################
################################################
################################################
# Analysis and plots

print "\n\n\nPlots and analysis\n"


if 'plotindividuallineplots' in toRun and toRun['plotindividuallineplots']:

    ################################################
    # Line plots of coverage scores with local optima
    # and with filtered local minima
    # and LOI
    # http://stackoverflow.com/questions/4624970/finding-local-maxima-minima-with-numpy-in-a-1d-numpy-array#9667121

    print "\nLine plots of coverage scores\n"

    lpDir = outDir + "/linePlots"
    if not os.path.exists(lpDir):
        os.makedirs(lpDir)

    minY = 10**6
    maxY = 0
    for dataSet in dataSets:
        # dataSet = 'pooledWT_alt_IC'
        #print dataSet
        if minY > np.nanmin(covScores[dataSet]):
            minY = np.nanmin(covScores[dataSet])
        if maxY < np.nanmax(covScores[dataSet]):
            maxY = np.nanmax(covScores[dataSet])
    #minY, maxY = conditionDistanceLimits(minY, maxY)
    # Make it symmetrical, if opposite signs
    if maxY > 0 and minY < 0 :
        maxMinMaxY = np.nanmax([abs(minY), abs(maxY)])
        maxY = maxMinMaxY
        minY = -maxMinMaxY
    minY *= 1.02
    maxY *= 1.02
    # Bounds
    if minY < -2:
        minY = -2
    if maxY > 2:
        maxY = 2

    # X-axis limits
    minX = 0
    maxX = int(np.ceil(float(chromSizes[thechrom]) / 1000000.0))

    # Plot
    idxTally = -1
    for dataSet in dataSets:
        # dataSet = 'pooledWT_alt_IC'
        #print dataSet

        idxTally += 1

        fig, ax = plt.subplots(1, 1, figsize=(figureWidth, 10))
        plots = {}
        plots2 = {}

        xVals = covScoreLoci
        yVals = covScores[dataSet]
        yValsS = covScoresSmoothed[dataSet]

        plots[dataSet], = ax.plot(xVals, yVals,
                                  '-', label=dataSet,
                                  color=contactMapIOmetadata[dataSet]['dataColor'],
                                  linewidth=2,
                                  linestyle='dashed')
        plots2[dataSet], = ax.plot(xVals, yValsS,
                                  '-', label=dataSet + ".smoothed",
                                  color=contactMapIOmetadata[dataSet]['dataColor'],
                                  linewidth=3,
                                  linestyle='solid')
        ax.grid(True)
        ax.set_ylabel('log2(coverage score / mean score)', fontsize=14)
        ax.set_xlabel('Mb', fontsize=14)
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(12)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(12)

        # a = np.diff(np.sign(np.diff(yValsS))).nonzero()[0] + 1 # local min+max -- NANs an issue here
        b = (np.diff(np.sign(np.diff(yValsS))) > 0).nonzero()[0] + 1 # local min
        c = (np.diff(np.sign(np.diff(yValsS))) < 0).nonzero()[0] + 1 # local max
        ax.plot(xVals[b], yValsS[b], 'o', markersize=8, markeredgewidth=2,
                markeredgecolor='black', markerfacecolor=contactMapIOmetadata[dataSet]['dataColor'] )
        ax.plot(xVals[c], yValsS[c], 'o', markersize=8, markeredgewidth=2,
                markeredgecolor='black', markerfacecolor=contactMapIOmetadata[dataSet]['dataColor'] )

        # local filtered minima
        if dataSet in filteredMin.keys():
            ax.plot(xVals[filteredMin[dataSet]], yValsS[filteredMin[dataSet]], 'v', markersize=10, markeredgewidth=2,
                    markeredgecolor='black', markerfacecolor=contactMapIOmetadata[dataSet]['dataColor'] )
            ax.plot(xVals[filteredMin[dataSet]], yValsS[filteredMin[dataSet]], 'v', markersize=10, markeredgewidth=2,
                    markeredgecolor='black', markerfacecolor=contactMapIOmetadata[dataSet]['dataColor'] )

        for aloi in theLoi:
            #print aloi
            #print LOI[aloi]

            chrCoords = LOI[aloi].split(":")
            theChr = chrCoords[0]
            if theChr == thechrom:
                #print "\tPlot it!"
                startStop = chrCoords[1].split("-")
                startStop = map(int, startStop)
                #startStop = [round(x / float(theRez)) for x in startStop]
                geneMidPoint = (startStop[0] + (startStop[1] - startStop[0]))/10.0**6

                ax.axvline(geneMidPoint, color='k', linestyle='--', linewidth=0.5)
                # http://matplotlib.org/examples/pylab_examples/fonts_demo.html
                font0 = FontProperties()
                font1 = font0.copy()
                font1.set_size('12')
                ax.text(geneMidPoint, maxY, aloi, rotation=90, verticalalignment="top", fontproperties=font1)

        #startStop = [round(x / float(theRez)) for x in startStop]
        #ax.set_title(aloi + "\n(" + chrCoords[0] + ":" + sizeof_fmt(startStop[0]) + ")", fontsize=14)

        ax.set_xlim(minX, maxX)
        ax.set_ylim(minY, maxY)
        #ax.legend(loc='best', fontsize=14)
        ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=2, mode="expand", borderaxespad=0., fontsize=14)

        # Tidy up and save
        fig.suptitle(dataSet + " coverage Scores +localOptima\n" + figTitleSuffix, fontsize=20, y=1.01)
        #fig.subplots_adjust(top=1)
        #fig.tight_layout()
        fig.savefig(lpDir + "/coverageScores." + dataSet +".localOptima." + figNameSuffix + ".pdf", dpi=150, bbox_inches='tight')
        fig.savefig(lpDir + "/coverageScores." + dataSet +".localOptima." + figNameSuffix + ".png", dpi=150, bbox_inches='tight')

        if 'plotlineplotzooms' in toRun and toRun['plotlineplotzooms']:
            # Zoom around LOI

            lpDirZ = lpDir + "/Zoomed"
            if not os.path.exists(lpDirZ):
                os.makedirs(lpDirZ)

            for aloi in theLoi:
                #print aloi
                #print LOI[aloi]

                chrCoords = LOI[aloi].split(":")
                theChr = chrCoords[0]
                if theChr == thechrom:
                    #print "\tPlot it!"

                    startStop = chrCoords[1].split("-")
                    startStop = map(int, startStop)
                    #startStop = [round(x / float(theRez)) for x in startStop]
                    # theLocus = startStop[0]
                    # #print str(theLocus)
                    geneMidPoint = (startStop[0] + (startStop[1] - startStop[0]))/10.0**6
                    plusMinusDist = round((zoomWin / 2.0) / 1000000.0)

                    # Zoom in
                    ax.set_xlim(max(minX, geneMidPoint - plusMinusDist),
                                 min(maxX, geneMidPoint + plusMinusDist))

                    ax.legend(loc='best', fontsize=14)
                    # ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                    #        ncol=2, mode="expand", borderaxespad=0., fontsize=14)

                    fig.set_size_inches(16, 8)

                    # Tidy up and save
                    fig.suptitle(dataSet + " coverage scores +localOptima\n" \
                            + aloi + " @ " + str(round(geneMidPoint, 1)) + "Mb; " \
                            + " zoom" + str(round(zoomWin/1000000.0, 1)) + "Mb; \n" \
                            + figTitleSuffix, fontsize=20, y=1.1)
                    #fig.subplots_adjust(top=1)
                    fig.tight_layout()
                    fig.savefig(lpDirZ + "/coverageScores." + dataSet + ".localOptima."
                                + aloi + "." + "zoom" + str(round(zoomWin/1000000.0,1)) + "Mb." \
                                + figNameSuffix + ".pdf" , dpi=150,
                                #bbox_inches='tight', )
                                bbox_inches=0, orientation='landscape', pad_inches=0)
                    fig.savefig(lpDirZ + "/coverageScores." + dataSet + ".localOptima."
                                + aloi + "." + "zoom" + str(round(zoomWin/1000000.0,1)) + "Mb." \
                                + figNameSuffix + ".png" , dpi=150,
                                #bbox_inches='tight', )
                                bbox_inches=0, orientation='landscape', pad_inches=0)

        plt.close()


    ################################################
    # Line plots of rescaled coverage scores with local optima
    # and with filtered local minima
    # and LOI
    # http://stackoverflow.com/questions/4624970/finding-local-maxima-minima-with-numpy-in-a-1d-numpy-array#9667121

    print "\nLine plots of rescaled coverage scores\n"

    lpDir = outDir + "/linePlots"
    if not os.path.exists(lpDir):
        os.makedirs(lpDir)

    minY = 0
    maxY = 1

    # X-axis limits
    minX = 0
    maxX = int(np.ceil(float(chromSizes[thechrom]) / 1000000.0))

    # Plot
    idxTally = -1
    for dataSet in dataSets:
        # dataSet = 'pooledWT_alt_IC'
        #print dataSet

        idxTally += 1

        fig, ax = plt.subplots(1, 1, figsize=(figureWidth, 10))
        plots = {}
        plots2 = {}

        xVals = covScoreLoci
        yVals = covScoresRescaled[dataSet]

        plots[dataSet], = ax.plot(xVals, yVals,
                                  '-', label=dataSet,
                                  color=contactMapIOmetadata[dataSet]['dataColor'],
                                  linewidth=2,
                                  linestyle='solid')
        ax.grid(True)
        ax.set_ylabel('coverage score [rescaled to 0 and 1]', fontsize=14)
        ax.set_xlabel('Mb', fontsize=14)
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(12)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(12)

        # a = np.diff(np.sign(np.diff(yValsS))).nonzero()[0] + 1 # local min+max -- NANs an issue here
        b = (np.diff(np.sign(np.diff(yValsS))) > 0).nonzero()[0] + 1 # local min
        c = (np.diff(np.sign(np.diff(yValsS))) < 0).nonzero()[0] + 1 # local max
        ax.plot(xVals[b], yVals[b], 'o', markersize=8, markeredgewidth=2,
                markeredgecolor='black', markerfacecolor=contactMapIOmetadata[dataSet]['dataColor'] )
        ax.plot(xVals[c], yVals[c], 'o', markersize=8, markeredgewidth=2,
                markeredgecolor='black', markerfacecolor=contactMapIOmetadata[dataSet]['dataColor'] )

        # local filtered minima
        if dataSet in filteredMin.keys():
            ax.plot(xVals[filteredMin[dataSet]], yVals[filteredMin[dataSet]], 'v', markersize=10, markeredgewidth=2,
                    markeredgecolor='black', markerfacecolor=contactMapIOmetadata[dataSet]['dataColor'] )
            ax.plot(xVals[filteredMin[dataSet]], yVals[filteredMin[dataSet]], 'v', markersize=10, markeredgewidth=2,
                    markeredgecolor='black', markerfacecolor=contactMapIOmetadata[dataSet]['dataColor'] )

        for aloi in theLoi:
            #print aloi
            #print LOI[aloi]

            chrCoords = LOI[aloi].split(":")
            theChr = chrCoords[0]
            if theChr == thechrom:
                #print "\tPlot it!"
                startStop = chrCoords[1].split("-")
                startStop = map(int, startStop)
                #startStop = [round(x / float(theRez)) for x in startStop]
                geneMidPoint = (startStop[0] + (startStop[1] - startStop[0]))/10.0**6

                ax.axvline(geneMidPoint, color='k', linestyle='--', linewidth=0.5)
                # http://matplotlib.org/examples/pylab_examples/fonts_demo.html
                font0 = FontProperties()
                font1 = font0.copy()
                font1.set_size('12')
                ax.text(geneMidPoint, maxY, aloi, rotation=90, verticalalignment="top", fontproperties=font1)

        #startStop = [round(x / float(theRez)) for x in startStop]
        #ax.set_title(aloi + "\n(" + chrCoords[0] + ":" + sizeof_fmt(startStop[0]) + ")", fontsize=14)

        ax.set_xlim(minX, maxX)
        ax.set_ylim(minY, maxY)
        #ax.legend(loc='best', fontsize=14)
        ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=2, mode="expand", borderaxespad=0., fontsize=14)

        # Tidy up and save
        fig.suptitle(dataSet + " coverage Scores rescaled +localOptima\n" + figTitleSuffix, fontsize=20, y=1.01)
        #fig.subplots_adjust(top=1)
        #fig.tight_layout()
        fig.savefig(lpDir + "/coverageScoresRescaled." + dataSet +".localOptima." + figNameSuffix + ".pdf", dpi=150, bbox_inches='tight')
        fig.savefig(lpDir + "/coverageScoresRescaled." + dataSet +".localOptima." + figNameSuffix + ".png", dpi=150, bbox_inches='tight')

        if 'plotlineplotzooms' in toRun and toRun['plotlineplotzooms']:

            lpDirZ = lpDir + "/Zoomed"
            if not os.path.exists(lpDirZ):
                os.makedirs(lpDirZ)

            # Zoom around LOI
            for aloi in theLoi:
                #print aloi
                #print LOI[aloi]

                chrCoords = LOI[aloi].split(":")
                theChr = chrCoords[0]
                if theChr == thechrom:
                    #print "\tPlot it!"

                    startStop = chrCoords[1].split("-")
                    startStop = map(int, startStop)
                    #startStop = [round(x / float(theRez)) for x in startStop]
                    # theLocus = startStop[0]
                    # #print str(theLocus)
                    geneMidPoint = (startStop[0] + (startStop[1] - startStop[0]))/10.0**6
                    plusMinusDist = round((zoomWin / 2.0) / 1000000.0)

                    # Zoom in
                    ax.set_xlim(max(minX, geneMidPoint - plusMinusDist),
                                 min(maxX, geneMidPoint + plusMinusDist))

                    ax.legend(loc='best', fontsize=14)
                    # ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                    #        ncol=2, mode="expand", borderaxespad=0., fontsize=14)

                    fig.set_size_inches(16, 8)

                    # Tidy up and save
                    fig.suptitle(dataSet + " coverage scores rescaled +localOptima\n" \
                            + aloi + " @ " + str(round(geneMidPoint, 1)) + "Mb; " \
                            + " zoom" + str(round(zoomWin/1000000.0, 1)) + "Mb; \n" \
                            + figTitleSuffix, fontsize=20, y=1.1)
                    #fig.subplots_adjust(top=1)
                    fig.tight_layout()
                    fig.savefig(lpDirZ + "/coverageScoresRescaled." + dataSet + ".localOptima."
                                + aloi + "." + "zoom" + str(round(zoomWin/1000000.0,1)) + "Mb." \
                                + figNameSuffix + ".pdf" , dpi=150,
                                #bbox_inches='tight', )
                                bbox_inches=0, orientation='landscape', pad_inches=0)
                    fig.savefig(lpDirZ + "/coverageScoresRescaled." + dataSet + ".localOptima."
                                + aloi + "." + "zoom" + str(round(zoomWin/1000000.0,1)) + "Mb." \
                                + figNameSuffix + ".png" , dpi=150,
                                #bbox_inches='tight', )
                                bbox_inches=0, orientation='landscape', pad_inches=0)

        plt.close()


if 'plotgroupedlineplots' in toRun and toRun['plotgroupedlineplots']:

    ################################################
    # Line plots of coverage zscores with filtered local minima
    # and LOI
    # http://stackoverflow.com/questions/4624970/finding-local-maxima-minima-with-numpy-in-a-1d-numpy-array#9667121

    print "\nLine plots of coverage z-scores with minima\n"

    lpDir = outDir + "/linePlotsGrouped"
    if not os.path.exists(lpDir):
        os.makedirs(lpDir)

    # Y-axis limits
    # minY = -2
    # maxY = 2
    minY = np.inf
    maxY = -np.inf
    for dataSet in dataSets:
        # dataSet = 'pooledWT_alt_IC'
        #print dataSet
        if minY > np.nanmin(covScoresZscored[dataSet]):
            minY = np.nanmin(covScoresZscored[dataSet])
        if maxY < np.nanmax(covScoresZscored[dataSet]):
            maxY = np.nanmax(covScoresZscored[dataSet])
    #minY, maxY = conditionDistanceLimits(minY, maxY)
    # Make it symmetrical, if opposite signs
    if maxY > 0 and minY < 0 :
        maxMinMaxY = np.nanmax([abs(minY), abs(maxY)])
        maxY = maxMinMaxY
        minY = -maxMinMaxY
    minY *= 1.02
    maxY *= 1.02
    # Bounds
    #if minY < -4:
    #    minY = -4
    #if maxY > 4:
    #    maxY = 4

    # X-axis limits
    minX = 0
    maxX = int(np.ceil(float(chromSizes[thechrom]) / 1000000.0))

    for Xtype in ['Xa', 'Xi']:
        #print(Xtype)
        XtypeString = '_' + Xtype
        whichData = [XtypeString in i for i in dataSets]
        np.where(whichData)
        dataSubSets = list(compress(dataSets, whichData))

        fig, ax = plt.subplots(1, 1, figsize=(figureWidth, 10))
        plots = {}

        for dataSet in dataSubSets:
            # dataSet = 'pooledWT_alt_IC'
            #print dataSet

            xVals = covScoreLoci
            yVals = covScoresZscored[dataSet]

            plots[dataSet], = ax.plot(xVals, yVals,
                                      '-', label=dataSet,
                                      color=contactMapIOmetadata[dataSet]['dataColor'],
                                      linewidth=4)
            ax.grid(True)
            ax.set_ylabel('log2(coverage score / mean score) [z-scored]', fontsize=14)
            ax.set_xlabel('Mb', fontsize=14)
            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(12)
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(12)

            # local filtered minima
            if dataSet in filteredMin.keys():
                ax.plot(xVals[filteredMin[dataSet]], yVals[filteredMin[dataSet]], 'v', markersize=10, markeredgewidth=1,
                        markeredgecolor='black', markerfacecolor=contactMapIOmetadata[dataSet]['dataColor'] )

        for aloi in theLoi:
            #print aloi
            #print LOI[aloi]

            chrCoords = LOI[aloi].split(":")
            theChr = chrCoords[0]
            if theChr == thechrom:
                #print "\tPlot it!"
                startStop = chrCoords[1].split("-")
                startStop = map(int, startStop)
                #startStop = [round(x / float(theRez)) for x in startStop]
                geneMidPoint = (startStop[0] + (startStop[1] - startStop[0]))/10.0**6

                ax.axvline(geneMidPoint, color='k', linestyle='--', linewidth=0.5)
                # http://matplotlib.org/examples/pylab_examples/fonts_demo.html
                font0 = FontProperties()
                font1 = font0.copy()
                font1.set_size('12')
                ax.text(geneMidPoint, maxY, aloi, rotation=90, verticalalignment="top", fontproperties=font1)

        #startStop = [round(x / float(theRez)) for x in startStop]
        #ax.set_title(aloi + "\n(" + chrCoords[0] + ":" + sizeof_fmt(startStop[0]) + ")", fontsize=14)

        ax.set_xlim(minX, maxX)
        ax.set_ylim(minY, maxY)
        #ax.legend(loc='best', fontsize=14)
        ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=2, mode="expand", borderaxespad=0., fontsize=14)

        # Tidy up and save
        fig.suptitle(Xtype + " coverage Z-scores +filteredMinima\n" + figTitleSuffix, fontsize=20, y=1.01)
        #fig.subplots_adjust(top=1)
        #fig.tight_layout()
        fig.savefig(lpDir + "/coverageZscores." + Xtype + ".filteredMinima." + figNameSuffix + ".pdf", dpi=150, bbox_inches='tight')
        fig.savefig(lpDir + "/coverageZscores." + Xtype + ".filteredMinima." + figNameSuffix + ".png", dpi=150, bbox_inches='tight')

        if 'plotlineplotzooms' in toRun and toRun['plotlineplotzooms']:
            # Zoom around LOI

            lpDirZ = lpDir + "/Zoomed"
            if not os.path.exists(lpDirZ):
                os.makedirs(lpDirZ)

            for aloi in theLoi:
                #print aloi
                #print LOI[aloi]

                chrCoords = LOI[aloi].split(":")
                theChr = chrCoords[0]
                if theChr == thechrom:
                    #print "\tPlot it!"

                    startStop = chrCoords[1].split("-")
                    startStop = map(int, startStop)
                    #startStop = [round(x / float(theRez)) for x in startStop]
                    # theLocus = startStop[0]
                    # #print str(theLocus)
                    geneMidPoint = (startStop[0] + (startStop[1] - startStop[0]))/10.0**6
                    plusMinusDist = round((zoomWin / 2.0) / 1000000.0)

                    # Zoom in
                    ax.set_xlim(max(minX, geneMidPoint - plusMinusDist),
                                 min(maxX, geneMidPoint + plusMinusDist))

                    ax.legend(loc='best', fontsize=14)
                    # ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                    #        ncol=2, mode="expand", borderaxespad=0., fontsize=14)

                    fig.set_size_inches(20, 10)

                    # Tidy up and save
                    fig.suptitle(Xtype + " coverage Z-scores +filteredMinima\n" \
                            + aloi + " @ " + str(round(geneMidPoint, 1)) + "Mb; " \
                            + " zoom" + str(round(zoomWin/1000000.0, 1)) + "Mb; \n" \
                            + figTitleSuffix, fontsize=20, y=1.1)
                    #fig.subplots_adjust(top=1)
                    fig.tight_layout()
                    fig.savefig(lpDirZ + "/coverageZscores." + Xtype + ".filteredMinima."
                                + aloi + "." + "zoom" + str(round(zoomWin/1000000.0,1)) + "Mb." \
                                + figNameSuffix + ".pdf" , dpi=150,
                                #bbox_inches='tight', )
                                bbox_inches=0, orientation='landscape', pad_inches=0)
                    fig.savefig(lpDirZ + "/coverageZscores." + Xtype + ".filteredMinima."
                                + aloi + "." + "zoom" + str(round(zoomWin/1000000.0,1)) + "Mb." \
                                + figNameSuffix + ".png" , dpi=150,
                                #bbox_inches='tight', )
                                bbox_inches=0, orientation='landscape', pad_inches=0)

            plt.close()


    ################################################
    # Line plots of coverage zscores w/o filtered local minima
    # and LOI

    print "\nLine plots of coverage z-scores\n"

    lpDir = outDir + "/linePlotsGrouped"
    if not os.path.exists(lpDir):
        os.makedirs(lpDir)

    # Y-axis limits
    # minY = -2
    # maxY = 2
    minY = np.inf
    maxY = -np.inf
    for dataSet in dataSets:
        # dataSet = 'pooledWT_alt_IC'
        #print dataSet
        if minY > np.nanmin(covScoresZscored[dataSet]):
            minY = np.nanmin(covScoresZscored[dataSet])
        if maxY < np.nanmax(covScoresZscored[dataSet]):
            maxY = np.nanmax(covScoresZscored[dataSet])
    #minY, maxY = conditionDistanceLimits(minY, maxY)
    # Make it symmetrical, if opposite signs
    if maxY > 0 and minY < 0 :
        maxMinMaxY = np.nanmax([abs(minY), abs(maxY)])
        maxY = maxMinMaxY
        minY = -maxMinMaxY
    minY *= 1.02
    maxY *= 1.02
    # Bounds
    #if minY < -4:
    #    minY = -4
    #if maxY > 4:
    #    maxY = 4

    # X-axis limits
    minX = 0
    maxX = int(np.ceil(float(chromSizes[thechrom]) / 1000000.0))

    for Xtype in ['Xa', 'Xi']:
        #print(Xtype)
        XtypeString = '_' + Xtype
        whichData = [XtypeString in i for i in dataSets]
        np.where(whichData)
        dataSubSets = list(compress(dataSets, whichData))

        fig, ax = plt.subplots(1, 1, figsize=(figureWidth, 10))
        plots = {}

        for dataSet in dataSubSets:
            # dataSet = 'pooledWT_alt_IC'
            #print dataSet

            xVals = covScoreLoci
            yVals = covScoresZscored[dataSet]

            plots[dataSet], = ax.plot(xVals, yVals,
                                      '-', label=dataSet,
                                      color=contactMapIOmetadata[dataSet]['dataColor'],
                                      linewidth=4)
            ax.grid(True)
            ax.set_ylabel('log2(coverage score / mean score) [z-scored]', fontsize=14)
            ax.set_xlabel('Mb', fontsize=14)
            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(12)
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(12)

            ## local filtered minima
            #if dataSet in filteredMin.keys():
            #    ax.plot(xVals[filteredMin[dataSet]], yVals[filteredMin[dataSet]], 'v', markersize=10, markeredgewidth=1,
            #            markeredgecolor='black', markerfacecolor=contactMapIOmetadata[dataSet]['dataColor'] )

        for aloi in theLoi:
            #print aloi
            #print LOI[aloi]

            chrCoords = LOI[aloi].split(":")
            theChr = chrCoords[0]
            if theChr == thechrom:
                #print "\tPlot it!"
                startStop = chrCoords[1].split("-")
                startStop = map(int, startStop)
                #startStop = [round(x / float(theRez)) for x in startStop]
                geneMidPoint = (startStop[0] + (startStop[1] - startStop[0]))/10.0**6

                ax.axvline(geneMidPoint, color='k', linestyle='--', linewidth=0.5)
                # http://matplotlib.org/examples/pylab_examples/fonts_demo.html
                font0 = FontProperties()
                font1 = font0.copy()
                font1.set_size('12')
                ax.text(geneMidPoint, maxY, aloi, rotation=90, verticalalignment="top", fontproperties=font1)

        #startStop = [round(x / float(theRez)) for x in startStop]
        #ax.set_title(aloi + "\n(" + chrCoords[0] + ":" + sizeof_fmt(startStop[0]) + ")", fontsize=14)

        ax.set_xlim(minX, maxX)
        ax.set_ylim(minY, maxY)
        #ax.legend(loc='best', fontsize=14)
        ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=2, mode="expand", borderaxespad=0., fontsize=14)

        # Tidy up and save
        fig.suptitle(Xtype + " coverage Z-scores\n" + figTitleSuffix, fontsize=20, y=1.01)
        #fig.subplots_adjust(top=1)
        #fig.tight_layout()
        fig.savefig(lpDir + "/coverageZscores." + Xtype + "." + figNameSuffix + ".pdf", dpi=150, bbox_inches='tight')
        fig.savefig(lpDir + "/coverageZscores." + Xtype + "." + figNameSuffix + ".png", dpi=150, bbox_inches='tight')

        if 'plotlineplotzooms' in toRun and toRun['plotlineplotzooms']:
            # Zoom around LOI

            lpDirZ = lpDir + "/Zoomed"
            if not os.path.exists(lpDirZ):
                os.makedirs(lpDirZ)

            for aloi in theLoi:
                #print aloi
                #print LOI[aloi]

                chrCoords = LOI[aloi].split(":")
                theChr = chrCoords[0]
                if theChr == thechrom:
                    #print "\tPlot it!"

                    startStop = chrCoords[1].split("-")
                    startStop = map(int, startStop)
                    #startStop = [round(x / float(theRez)) for x in startStop]
                    # theLocus = startStop[0]
                    # #print str(theLocus)
                    geneMidPoint = (startStop[0] + (startStop[1] - startStop[0]))/10.0**6
                    plusMinusDist = round((zoomWin / 2.0) / 1000000.0)

                    # Zoom in
                    ax.set_xlim(max(minX, geneMidPoint - plusMinusDist),
                                 min(maxX, geneMidPoint + plusMinusDist))

                    ax.legend(loc='best', fontsize=14)
                    # ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                    #        ncol=2, mode="expand", borderaxespad=0., fontsize=14)

                    fig.set_size_inches(20, 10)

                    # Tidy up and save
                    fig.suptitle(Xtype + " coverage Z-scores +filteredMinima\n" \
                            + aloi + " @ " + str(round(geneMidPoint, 1)) + "Mb; " \
                            + " zoom" + str(round(zoomWin/1000000.0, 1)) + "Mb; \n" \
                            + figTitleSuffix, fontsize=20, y=1.1)
                    #fig.subplots_adjust(top=1)
                    fig.tight_layout()
                    fig.savefig(lpDirZ + "/coverageZscores." + Xtype + ".filteredMinima."
                                + aloi + "." + "zoom" + str(round(zoomWin/1000000.0,1)) + "Mb." \
                                + figNameSuffix + ".pdf" , dpi=150,
                                #bbox_inches='tight', )
                                bbox_inches=0, orientation='landscape', pad_inches=0)
                    fig.savefig(lpDirZ + "/coverageZscores." + Xtype + ".filteredMinima."
                                + aloi + "." + "zoom" + str(round(zoomWin/1000000.0,1)) + "Mb." \
                                + figNameSuffix + ".png" , dpi=150,
                                #bbox_inches='tight', )
                                bbox_inches=0, orientation='landscape', pad_inches=0)

            plt.close()


if 'genbedgraphs' in toRun and toRun['genbedgraphs']:

    ################################################
    # Bedgraphs: regular and z-scored
    # Just bedgraphs for now
    # #  and beds of topBottom 5%

    print "\nBedgraphs of Coverage scores\n"

    bedgraphDir = outDir + "/bedgraphs"
    if not os.path.exists(bedgraphDir):
        os.makedirs(bedgraphDir)

    # bedDir = outDir + "/beds"
    # if not os.path.exists(bedDir):
    #     os.makedirs(bedDir)

    idxTally = -1
    for dataSet in dataSets:
        # dataSet = 'pooledWT_alt_IC'
        #print dataSet
    
        idxTally += 1

        xVals = covScoreLoci

        ##########
        # bedgraphs of raw coverages scores

        yVals = covScoresRaw[dataSet]

        bedgraphFileName = bedgraphDir + "/rawCoverageScores." + dataSet + "." + bedNameSuffix + ".bedGraph.gz"
        f = gzip.open(bedgraphFileName, "wb")
        for i in range(len(xVals)):
            if ~np.isnan(yVals[i]):
                f.write(
                    thechrom + "\t" + \
                    str(int(round(xVals[i] * 1.0*10**6) - (intTheRez/2))) + "\t" + \
                    str(int(round(xVals[i] * 1.0*10**6) + (intTheRez/2))) + "\t" + \
                    str(yVals[i]) + "\n")
        f.close()

        ##########
        # bedgraphs of rescaled coverages scores

        yVals = covScoresRescaled[dataSet]

        bedgraphFileName = bedgraphDir + "/rescaledCoverageScores." + dataSet + "." + bedNameSuffix + ".bedGraph.gz"
        f = gzip.open(bedgraphFileName, "wb")
        for i in range(len(xVals)):
            if ~np.isnan(yVals[i]):
                f.write(
                    thechrom + "\t" + \
                    str(int(round(xVals[i] * 1.0*10**6) - (intTheRez/2))) + "\t" + \
                    str(int(round(xVals[i] * 1.0*10**6) + (intTheRez/2))) + "\t" + \
                    str(yVals[i]) + "\n")
        f.close()

        ##########
        # bedgraphs of Coverages scores

        yVals = covScoresSmoothed[dataSet]

        bedgraphFileName = bedgraphDir + "/coverageScores." + dataSet + "." + bedNameSuffix + ".bedGraph.gz"
        f = gzip.open(bedgraphFileName, "wb")
        for i in range(len(xVals)):
            if ~np.isnan(yVals[i]):
                f.write(
                    thechrom + "\t" + \
                    str(int(round(xVals[i] * 1.0*10**6) - (intTheRez/2))) + "\t" + \
                    str(int(round(xVals[i] * 1.0*10**6) + (intTheRez/2))) + "\t" + \
                    str(yVals[i]) + "\n")
        f.close()

        ##########
        # bedgraphs of Coverage zscores

        yVals = covScoresZscored[dataSet]

        bedgraphFileName = bedgraphDir + "/coverageZscores." + dataSet + "." + bedNameSuffix + ".bedGraph.gz"
        f = gzip.open(bedgraphFileName, "wb")
        for i in range(len(xVals)):
            if ~np.isnan(yVals[i]):
                f.write(
                    thechrom + "\t" + \
                    str(int(round(xVals[i] * 1.0*10**6) - (intTheRez/2))) + "\t" + \
                    str(int(round(xVals[i] * 1.0*10**6) + (intTheRez/2))) + "\t" + \
                    str(yVals[i]) + "\n")
        f.close()


        # ##########
        # # Co-ordinates of outlier values and save as BED files
        #
        # bedFileName = bedDir + "/" + "topBot" + str(outThresh) + "pct.coverageScores." \
        #               + dataSet + "." + bedNameSuffix + ".bedGraph.gz"
        # f = gzip.open(bedFileName, "wb")
        #
        # # yVals = np.linspace(1, 100, num=100)
        # # xVals = yVals
        #
        # lowYval = np.nanpercentile(yVals, outThresh)
        # highYval = np.nanpercentile(yVals, 100 - outThresh)
        #
        # lowYvalsBoolIdx = np.less_equal(yVals, lowYval)
        # highYvalsBoolIdx = np.greater_equal(yVals, highYval)
        # hiloYvalsIdx = np.where(lowYvalsBoolIdx | highYvalsBoolIdx)
        #
        # old_i = np.nanmin(hiloYvalsIdx) - 1
        # start = str(int((xVals[np.nanmin(hiloYvalsIdx)] * 1.0*10**6) - (intTheRez/2)))
        # for i in np.nditer(hiloYvalsIdx):
        #     if i != (old_i + 1) or i == np.nanmax(hiloYvalsIdx):
        #         f.write(
        #             thechrom + "\t" + start + "\t" + \
        #             str(int(round(xVals[old_i] * 1.0*10**6) + (intTheRez/2))) + "\n")
        #         start = str(int(round(xVals[i] * 1.0*10**6) - (intTheRez/2)))
        #     old_i = i
        # f.close()


if 'plotdistributions' in toRun and toRun['plotdistributions']:

    ################################################
    # Distributions as violin plots

    print "\nDistributions of coverage scores\n"

    distDir = outDir + "/distributions"
    if not os.path.exists(distDir):
        os.makedirs(distDir)

    # Y-axis limits
    # minY = -2
    # maxY = 2
    minY = np.inf
    maxY = -np.inf
    for dataSet in dataSets:
        # dataSet = 'pooledWT_alt_IC'
        #print dataSet
        if minY > np.nanmin(covScoresZscored[dataSet]):
            minY = np.nanmin(covScoresZscored[dataSet])
        if maxY < np.nanmax(covScoresZscored[dataSet]):
            maxY = np.nanmax(covScoresZscored[dataSet])
    #minY, maxY = conditionDistanceLimits(minY, maxY)
    # Make it symmetrical, if opposite signs
    # if maxY > 0 and minY < 0 :
    #     maxMinMaxY = np.nanmax([abs(minY), abs(maxY)])
    #     maxY = maxMinMaxY
    #     minY = -maxMinMaxY
    minY *= 1.02
    maxY *= 1.02
    # Bounds
    #if minY < -4:
    #    minY = -4
    #if maxY > 4:
    #    maxY = 4

    sns.set(style="darkgrid")
    # Set up the matplotlib figure
    fig, axes = plt.subplots(1, len(covScoresZscored), sharey=True, figsize=(12, 6))

    axIdx = -1
    for dataSet in dataSets:
        axIdx += 1

        # Draw a violinplot with a narrower bandwidth than the default
        sns.violinplot(data=covScoresZscored[dataSet], color=contactMapIOmetadata[dataSet]['dataColor'],
                       cut=0, linewidth=1, ax=axes[axIdx])
                       # bw=0.2

        # Finalize the figure
        axes[axIdx].set(ylim=(minY, maxY))
        sns.despine(left=True, bottom=True)
        axes[axIdx].set_title(dataSet)
        # axes[axIdx].set_xlabel('')
        axes[axIdx].set_xticklabels('')
        #fig.tight_layout()

    axes[0].set_ylabel("Coverage Z-scores\n" + " (n=" + str(covScoresZscored[dataSet].shape[0]) + ")" )

    # Tidy up and save
    fig.suptitle("Violin plot of coverage Z-score distributions\n" + figTitleSuffix, fontsize=12, y=1.05)
    #fig.subplots_adjust(top=1)
    fig.tight_layout()
    fig.savefig(distDir + "/violinPlot.coverageZscores." + figNameSuffix + ".pdf", dpi=150, bbox_inches='tight')
    fig.savefig(distDir + "/violinPlot.coverageZscores." + figNameSuffix + ".png", dpi=150, bbox_inches='tight')

    plt.close()


    ################################################
    # Distributions as violin plots II
    # http://seaborn.pydata.org/tutorial/categorical.html#categorical-scatterplots
    #
    # tips = sns.load_dataset("tips")

    # Denanify
    # Done above now.
    #  allTheNaNsBoolIdx = [False] * covScoresZscored[covScoresZscored.keys()[0]].shape[0]
    # for dataSet in dataSets:
    #     allTheNaNsBoolIdx = allTheNaNsBoolIdx | np.isnan(covScoresZscored[dataSet])
    covScoresZscoredDeNaned = {}
    for dataSet in dataSets:
        covScoresZscoredDeNaned[dataSet] = covScoresZscored[dataSet][~allTheNaNsBoolIdx]
    covScoreLociDeNaned = covScoreLoci[~allTheNaNsBoolIdx]
    covScoreLociRTDdeNaned = covScoreLociRTD[~allTheNaNsBoolIdx]
    # covScoreLociRTDclassesDeNaned = covScoreLociRTDclasses[~allTheNaNsBoolIdx]
    #       TypeError: only integer arrays with one element can be converted to an index
    covScoreLociRTDclassesDeNaned = list(compress(covScoreLociRTDclasses, ~allTheNaNsBoolIdx))

    # Get cell types
    cellTypes = [w.replace('_Xi', '') for w in dataSets]
    cellTypes = [w.replace('_Xa', '') for w in cellTypes]
    cellTypes = list(set(cellTypes))

    # Construct a dataframe
    # alltheData = pd.DataFrame()
    theData = []
    theType = []
    theXstate = []
    for cellType in cellTypes:
        for Xtype in ['Xa', 'Xi']:
            thisData = covScoresZscoredDeNaned[cellType + "_" + Xtype]
            thisType = [cellType] * len(thisData)
            thisXstate = [Xtype] * len(thisData)
            theData.extend(thisData)
            theType.extend(thisType)
            theXstate.extend(thisXstate)
    alltheData = pd.DataFrame()
    alltheData  = alltheData.append(theData)
    alltheData.column = 'IS'
    alltheData['celltype'] = theType
    alltheData['Xstate'] = theXstate
    alltheData.columns = ['IS', 'celltype', 'Xstate']

    # Draw a violinplot with a narrower bandwidth than the default
    sns.set(style="darkgrid")
    # Set up the matplotlib figure
    fig = plt.figure(figsize=(10,10))

    #sns.violinplot(y='IS', x='celltype', hue='Xstate', data=alltheData, split=True, palette=sns.color_palette("pastel", 2))
    sns.violinplot(y='IS', x='celltype', hue='Xstate', data=alltheData, palette=sns.color_palette("pastel", 2),
                   cut=0, linewidth=1, )
    sns.swarmplot(y='IS', x='celltype', hue='Xstate', data=alltheData, split=True)

    # Finalize the figure
    sns.despine(left=True, bottom=True)

    # Tidy up and save
    fig.suptitle("Violin plot of coverage Z-score distributions\n" + figTitleSuffix, fontsize=12, y=1.01)
    #fig.subplots_adjust(top=1)
    #fig.tight_layout()
    fig.savefig(distDir + "/violinPlot.split.coverageZscores." + figNameSuffix + ".pdf", dpi=150, bbox_inches='tight')
    fig.savefig(distDir + "/violinPlot.split.coverageZscores." + figNameSuffix + ".png", dpi=150, bbox_inches='tight')

    plt.close()


if 'plotpairplots' in toRun and toRun['plotpairplots']:

    ################################################
    # Distributions as pair plots
    # https://stanford.edu/~mwaskom/software/seaborn/generated/seaborn.pairplot.html#seaborn.pairplot
    # https://stanford.edu/~mwaskom/software/seaborn/generated/seaborn.PairGrid.html#seaborn.PairGrid
    # https://stanford.edu/~mwaskom/software/seaborn/generated/seaborn.FacetGrid.html#seaborn.FacetGrid
    # http://web.stanford.edu/~mwaskom/software/seaborn/generated/seaborn.kdeplot.html

    # CAN't PRINT THESE SIDE-BY-SIDE
    # http://stackoverflow.com/questions/23969619/plotting-with-seaborn-using-the-matplotlib-object-oriented-interface

    print "\nDistributions of as pair plots\n"

    # Y-axis limits
    # minY = -2
    # maxY = 2
    minY = np.inf
    maxY = -np.inf
    for dataSet in dataSets:
        # dataSet = 'pooledWT_alt_IC'
        #print dataSet
        if minY > np.nanmin(covScoresZscored[dataSet]):
            minY = np.nanmin(covScoresZscored[dataSet])
        if maxY < np.nanmax(covScoresZscored[dataSet]):
            maxY = np.nanmax(covScoresZscored[dataSet])
    #minY, maxY = conditionDistanceLimits(minY, maxY)
    # Make it symmetrical, if opposite signs
    if maxY > 0 and minY < 0 :
         maxMinMaxY = np.nanmax([abs(minY), abs(maxY)])
         maxY = maxMinMaxY
         minY = -maxMinMaxY
    minY *= 1.02
    maxY *= 1.02
    # Bounds
    #if minY < -4:
    #    minY = -4
    #if maxY > 4:
    #    maxY = 4

    # Denanify
    # Done above now.
    # allTheNaNsBoolIdx = [False] * covScoresZscored[covScoresZscored.keys()[0]].shape[0]
    # for dataSet in dataSets:
    #     allTheNaNsBoolIdx = allTheNaNsBoolIdx | np.isnan(covScoresZscored[dataSet])
    covScoresZscoredDeNaned = {}
    for dataSet in dataSets:
        covScoresZscoredDeNaned[dataSet] = covScoresZscored[dataSet][~allTheNaNsBoolIdx]
    covScoreLociDeNaned = covScoreLoci[~allTheNaNsBoolIdx]
    covScoreLociRTDdeNaned = covScoreLociRTD[~allTheNaNsBoolIdx]
    # covScoreLociRTDclassesDeNaned = covScoreLociRTDclasses[~allTheNaNsBoolIdx]
    #       TypeError: only integer arrays with one element can be converted to an index
    covScoreLociRTDclassesDeNaned = list(compress(covScoreLociRTDclasses, ~allTheNaNsBoolIdx))

    for Xtype in ['Xa', 'Xi']:
        #print(Xtype)
        XtypeString = '_' + Xtype
        whichData = [XtypeString in i for i in dataSets]
        np.where(whichData)
        dataSubSets = list(compress(dataSets, whichData))

        diagcolors = []
        for dataSet in dataSubSets:
            diagcolors.append(contactMapIOmetadata[dataSet]['dataColor'])

        theData = pd.DataFrame(covScoresZscoredDeNaned)[dataSubSets]


        #################
        # All data points
        sns.set(style="darkgrid")
        plt.figure(figsize=(10,10))

        # Draw a violinplot with a narrower bandwidth than the default
        # sns.pairplot(data=altRatioData, kind="reg", markers="+", size=3)
        # sns.pairplot(data=altRatioData, size=3, diag_kind="kde")
        g = sns.PairGrid(data=theData)
        g = g.map_upper(plt.scatter, edgecolor="w", s=35)
        g = g.map_lower(sns.kdeplot, gridsize=50, n_levels=200,
                        #shade=True, shade_lowest=False, cmap="Blues")
                        lw=1, cmap="Blues_r")
        #g = g.map_lower(plt.hexbin, gridsize=50, cmap='inferno')
        g = g.map_diag(plt.hist, edgecolor="w",
                       bins=np.arange(minY, maxY, 0.5))
        #g = g.map_diag(sns.kdeplot, lw=2, legend=False)
        g = g.set(xlim=(minY, maxY),
                  ylim=(minY, maxY))

        # import pprint
        # pp = pprint.PrettyPrinter(indent=4)
        # for ax in g.diag_axes:
        #     #pp.pprint(vars(ax))
        #     diag_axes_axes = getattr(ax, '_axes')
        #     #pp.pprint(vars(diag_axes_axes))
        #     daGraph = ax.figure
        #     daGraph.set_facecolor('red')

        # Tidy up and save
        g.fig.suptitle(Xtype + " pair grid of coverage Z-score distributions\n" + figTitleSuffix +
        #               "\n(n=" + str(covScoresZscored[dataSet].shape[0]) + ")", fontsize=12, y=1.01)
                       "; n=" + str((~allTheNaNsBoolIdx).sum()), fontsize=12, y=1.05)
        #fig.subplots_adjust(top=1)
        g.fig.tight_layout()
        g.fig.savefig(distDir + "/pairGrid.coverageZscores." + Xtype + "." + figNameSuffix + ".pdf", dpi=150, bbox_inches='tight')
        g.fig.savefig(distDir + "/pairGrid.coverageZscores." + Xtype + "." + figNameSuffix + ".png", dpi=150, bbox_inches='tight')

        plt.close()


        ###############################
        # Grouped by distance from dxz4
        theData['FromDxz4'] = covScoreLociRTDclassesDeNaned

        sns.set(style="darkgrid")
        plt.figure(figsize=(10,10))

        # Draw a violinplot with a narrower bandwidth than the default
        # sns.pairplot(data=altRatioData, kind="reg", markers="+", size=3)
        # sns.pairplot(data=altRatioData, size=3, diag_kind="kde")
        g = sns.PairGrid(data=theData, hue="FromDxz4", palette=sns.diverging_palette(255, 133, l=60, n=5, center="dark"))
        g = g.map_offdiag(plt.scatter, edgecolor="w", s=35)
        g = g.map_diag(plt.hist, edgecolor="w",
                       bins=np.arange(minY, maxY, 0.5))
        #g = g.map_diag(sns.kdeplot, lw=2, legend=False)
        g = g.set(xlim=(minY, maxY),
                  ylim=(minY, maxY))
        g.add_legend();

        # import pprint
        # pp = pprint.PrettyPrinter(indent=4)
        # for ax in g.diag_axes:
        #     #pp.pprint(vars(ax))
        #     diag_axes_axes = getattr(ax, '_axes')
        #     #pp.pprint(vars(diag_axes_axes))
        #     daGraph = ax.figure
        #     daGraph.set_facecolor('red')

        # Tidy up and save
        g.fig.suptitle(Xtype + " pair grid of coverage Z-score distributions by distance from Dxz4\n" + figTitleSuffix +
        #               "\n(n=" + str(covScoresZscored[dataSet].shape[0]) + ")", fontsize=12, y=1.01)
                       "; n=" + str((~allTheNaNsBoolIdx).sum()), fontsize=12, y=1.05)
        #fig.subplots_adjust(top=1)
        #g.fig.tight_layout()
        g.fig.savefig(distDir + "/pairGridByDist.coverageZscores." + Xtype + "." + figNameSuffix + ".pdf", dpi=150, bbox_inches='tight')
        g.fig.savefig(distDir + "/pairGridByDist.coverageZscores." + Xtype + "." + figNameSuffix + ".png", dpi=150, bbox_inches='tight')

        plt.close()


if 'plotsimilaritymaps' in toRun and toRun['plotsimilaritymaps']:

    ################################################
    # heat maps with dendrograms
    # http://stackoverflow.com/questions/2982929/plotting-results-of-hierarchical-clustering-ontop-of-a-matrix-of-data-in-python
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html
    # http://stackoverflow.com/questions/23451028/matplotlib-pyplot-vs-matplotlib-pylab#23451106

    print "\ncoverage scores clustering\n"

    corrMatrixDirr = outDir + "/corrMatrices"
    if not os.path.exists(corrMatrixDirr):
        os.makedirs(corrMatrixDirr)

    # *** THIS DOES NOT INTITIALIZE TWO INSTANCES OF AN ARRAY,
    # *** BUT ONE INSTANCE WITH TWO ALIASES!
    # corrMatrix = distMatrix = []
    covScorComps = {}

    covScorComps['Pearson'] = []
    covScorComps['Euclidean'] = []
    #covScorComps['meanAbsErr'] = []
    #covScorComps['meanSquErr'] = []

    for dataSetA in dataSets:
        for dataSetB in dataSets:

            #print str(idxTally) + ": " + dataSetA + " " + dataSetB

            scorzA = covScoresZscored[dataSetA].copy()
            scorzB = covScoresZscored[dataSetB].copy()

            # Remove mutually unscorable regions
            AnanBoolIdx = np.logical_not(np.isnan(scorzA))
            BnanBoolIdx = np.logical_not(np.isnan(scorzB))
            scorzAnn = scorzA[AnanBoolIdx & BnanBoolIdx]
            scorzBnn = scorzB[AnanBoolIdx & BnanBoolIdx]

            covScorComps['Pearson'].append(sps.pearsonr(scorzAnn, scorzBnn)[0])

            diff = scorzAnn-scorzBnn
            # distMatrix.append(np.sum(diff**2)**0.5) # euclidean
            covScorComps['Euclidean'].append(np.linalg.norm(diff)) # euclidean

            #covScorComps['meanAbsErr'].append(np.sum(abs(diff))/len(scorzAnn))  # Mean absolute error
            #covScorComps['meanSquErr'].append(np.sum(diff**2)/len(scorzAnn)) # Mean Squared error

    # Reshape into matrices
    for metricType in covScorComps:
        #print metricType

        covScorComps[metricType] = np.reshape(covScorComps[metricType], (-1, len(dataSets)))
        covScorComps[metricType] = pd.DataFrame(covScorComps[metricType], index=dataSets, columns=dataSets)

    # Heat map plotting
    for metricType in covScorComps:
        #print metricType

        C = np.asarray(covScorComps[metricType])

        # # Try to order in the way you want the leaves sorted
        # # This does not appear to work
        # C = C[[0,3,1,2],:]
        # C = C[:,[0,3,1,2]]
        if metricType == 'Pearson':
            D = 1 - C
        else:
            D = C

        # sampleNames = np.asarray(corrMatrixList['chromosome'].columns.values.tolist())[[0,3,1,2]]
        # sampleNames = np.asarray(['WT Xa', 'Brain Xi', 'WT Xi', 'Brain Xa' ])[[0,3,1,2]]
        sampleNamesReordered = np.asarray(dataSets)[0:]

        # Compute and plot first dendrogram.
        fig = pylab.figure(figsize=(10,10))
        ax1 = fig.add_axes([0.09,0.1,0.2,0.6])
        #Y = sch.linkage(D, method='single')
        Y = sch.linkage(D, method='average')
        Z1 = sch.dendrogram(Y, count_sort='ascending', orientation='left')
        ax1.set_xticks([])
        ax1.set_yticks([])
        # http://stackoverflow.com/questions/30009062/get-rid-of-grey-background-in-python-matplotlib-bar-chart
        ax1.set_axis_bgcolor('white')

        # Compute and plot second dendrogram.
        ax2 = fig.add_axes([0.3,0.71,0.6,0.2])
        Y = sch.linkage(D, method='average')
        Z2 = sch.dendrogram(Y, count_sort='ascending')
        ax2.set_xticks([])
        ax2.set_yticks([])
        # http://stackoverflow.com/questions/30009062/get-rid-of-grey-background-in-python-matplotlib-bar-chart
        ax2.set_axis_bgcolor('white')

        # Plot distance matrix.
        axmatrix = fig.add_axes([0.3,0.1,0.6,0.6])
        idx1 = Z1['leaves']
        idx2 = Z2['leaves']
        C = C[idx1,:]
        C = C[:,idx2]
        im = axmatrix.matshow(C, aspect='auto', origin='lower', cmap=pylab.cm.YlGnBu)
        axmatrix.set_xticks([])
        axmatrix.set_yticks([])

        # Sample names
        axmatrix.set_xticks(range(D.shape[0]))
        axmatrix.set_xticklabels(np.asarray(sampleNamesReordered)[idx1], minor=False, fontsize=14, rotation=90)
        axmatrix.xaxis.set_label_position('bottom')
        axmatrix.xaxis.tick_bottom()
        # for tic in axmatrix.xaxis.get_major_ticks():
        #     tic.tick1On = tic.tick2On = False
        for t in axmatrix.xaxis.get_ticklines():
            t.set_visible(False)
        # http://stackoverflow.com/questions/16074392/getting-vertical-gridlines-to-appear-in-line-plot-in-matplotlib
        axmatrix.grid(False)

        # pylab.xticks(rotation=-90, fontsize=14)
        #
        # axmatrix.set_yticks(range(D.shape[0]))
        # axmatrix.set_yticklabels(np.asarray(sampleNamesReordered)[idx2], minor=False)
        # axmatrix.yaxis.set_label_position('right')
        # axmatrix.yaxis.tick_right()
        #
        # axcolor = fig.add_axes([0.94,0.1,0.02,0.6])

        # Plot colorbar.
        axcolor = fig.add_axes([0.91,0.1,0.02,0.6])
        pylab.colorbar(im, cax=axcolor)

        fig.suptitle(metricType + " matrix with hierarchical clustering using coverage Z-scores\n" + figTitleSuffix, fontsize=20, y=1.01)
        #fig.subplots_adjust(top=1)
        #fig.tight_layout()
        fig.savefig(corrMatrixDirr + "/dendroMatrix." + metricType + ".coverageZscores." + figNameSuffix + ".pdf", dpi=150, bbox_inches='tight')
        fig.savefig(corrMatrixDirr + "/dendroMatrix." + metricType + ".coverageZscores." + figNameSuffix + ".png", dpi=150, bbox_inches='tight')
        plt.close()


    ################################################
    # Rand Index of TADs

    print "\nTAD rand indices\n"

    # Convert boundaries to membership vectors
    TADvectors = {}
    for dataSet in filteredMin:
        #print "\n" + dataSet
        old_i = 0
        TADid = 0
        TADvectors[dataSet] = []
        for i in filteredMin[dataSet]:
            #print str(i)
            TADvectors[dataSet].extend( [TADid] * (i - old_i) )
            TADid += 1
            old_i = i
        TADvectors[dataSet].extend( [TADid] * (len(xVals) - old_i) )
        #print "\n" + str(len(TADvectors[dataSet]))
        # print "\n"+ str(len(xVals)) + "\n"

    # TADARI = {}
    # for (a, b) in it.combinations(dataSets,2):
    #     combo = a + "-vs-" + b
    #     #print combo
    #     X =  TADvectors[a]
    #     Y =  TADvectors[b]
    #     TADARI[combo] = adjusted_rand_score(X,Y)

    TADARI = []
    for dataSetA in dataSets:
        for dataSetB in dataSets:

            #print str(idxTally) + ": " + dataSetA + " " + dataSetB

            X = TADvectors[dataSetA]
            Y = TADvectors[dataSetB]

            TADARI.append(adjusted_rand_score(X,Y))

    # Reshape into matriz
    TADARI= np.reshape(TADARI, (-1, len(dataSets)))

    # Heat map plotting
    C = np.asarray(TADARI)
    D = 1 - C

    sampleNamesReordered = np.asarray(dataSets)[0:]

    # Compute and plot first dendrogram.
    fig = pylab.figure(figsize=(10,10))
    ax1 = fig.add_axes([0.09,0.1,0.2,0.6])
    #Y = sch.linkage(D, method='single')
    Y = sch.linkage(D, method='average')
    Z1 = sch.dendrogram(Y, count_sort='ascending', orientation='left')
    ax1.set_xticks([])
    ax1.set_yticks([])
    # http://stackoverflow.com/questions/30009062/get-rid-of-grey-background-in-python-matplotlib-bar-chart
    ax1.set_axis_bgcolor('white')

    # Compute and plot second dendrogram.
    ax2 = fig.add_axes([0.3,0.71,0.6,0.2])
    Y = sch.linkage(D, method='average')
    Z2 = sch.dendrogram(Y, count_sort='ascending')
    ax2.set_xticks([])
    ax2.set_yticks([])
    # http://stackoverflow.com/questions/30009062/get-rid-of-grey-background-in-python-matplotlib-bar-chart
    ax2.set_axis_bgcolor('white')

    # Plot distance matrix.
    axmatrix = fig.add_axes([0.3,0.1,0.6,0.6])
    idx1 = Z1['leaves']
    idx2 = Z2['leaves']
    C = C[idx1,:]
    C = C[:,idx2]
    im = axmatrix.matshow(C, aspect='auto', origin='lower', cmap=pylab.cm.YlGnBu)
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])

    # Sample names
    axmatrix.set_xticks(range(D.shape[0]))
    axmatrix.set_xticklabels(np.asarray(sampleNamesReordered)[idx1], minor=False, fontsize=14, rotation=90)
    axmatrix.xaxis.set_label_position('bottom')
    axmatrix.xaxis.tick_bottom()
    # for tic in axmatrix.xaxis.get_major_ticks():
    #     tic.tick1On = tic.tick2On = False
    for t in axmatrix.xaxis.get_ticklines():
        t.set_visible(False)
    # http://stackoverflow.com/questions/16074392/getting-vertical-gridlines-to-appear-in-line-plot-in-matplotlib
    axmatrix.grid(False)

    # pylab.xticks(rotation=-90, fontsize=14)
    #
    # axmatrix.set_yticks(range(D.shape[0]))
    # axmatrix.set_yticklabels(np.asarray(sampleNamesReordered)[idx2], minor=False)
    # axmatrix.yaxis.set_label_position('right')
    # axmatrix.yaxis.tick_right()
    #
    # axcolor = fig.add_axes([0.94,0.1,0.02,0.6])

    # Plot colorbar.
    axcolor = fig.add_axes([0.91,0.1,0.02,0.6])
    pylab.colorbar(im, cax=axcolor)

    fig.suptitle("Matrix of adjusted Rand Index of TADs with hierarchical clustering\n" + figTitleSuffix, fontsize=20, y=1.01)
    #fig.subplots_adjust(top=1)
    #fig.tight_layout()
    fig.savefig(corrMatrixDirr + "/dendroMatrix.TAD.withNans.adjRI." + figNameSuffix + ".pdf", dpi=150, bbox_inches='tight')
    fig.savefig(corrMatrixDirr + "/dendroMatrix.TAD.withNans.adjRI." + figNameSuffix + ".png", dpi=150, bbox_inches='tight')
    plt.close()


    ################################################
    # Rand Index of TADs after removing NaNs

    # Remove NaNs
    for dataSet in TADvectors:
        #print len(TADvectors[dataSet])
        TADvectors[dataSet] = np.array(TADvectors[dataSet], dtype=float)
        TADvectors[dataSet][nandices] = np.nan

    # TADARI = {}
    # for (a, b) in it.combinations(dataSets,2):
    #     combo = a + "-vs-" + b
    #     #print combo
    #     X =  TADvectors[a]
    #     Y =  TADvectors[b]
    #     TADARI[combo] = adjusted_rand_score(X,Y)

    TADARI = []
    for dataSetA in dataSets:
        for dataSetB in dataSets:

            #print str(idxTally) + ": " + dataSetA + " " + dataSetB

            X = TADvectors[dataSetA]
            Y = TADvectors[dataSetB]

            TADARI.append(adjusted_rand_score(X,Y))

    # Reshape into matriz
    TADARI= np.reshape(TADARI, (-1, len(dataSets)))

    # Heat map plotting
    C = np.asarray(TADARI)
    D = 1 - C

    sampleNamesReordered = np.asarray(dataSets)[0:]

    # Compute and plot first dendrogram.
    fig = pylab.figure(figsize=(10,10))
    ax1 = fig.add_axes([0.09,0.1,0.2,0.6])
    #Y = sch.linkage(D, method='single')
    Y = sch.linkage(D, method='average')
    Z1 = sch.dendrogram(Y, count_sort='ascending', orientation='left')
    ax1.set_xticks([])
    ax1.set_yticks([])
    # http://stackoverflow.com/questions/30009062/get-rid-of-grey-background-in-python-matplotlib-bar-chart
    ax1.set_axis_bgcolor('white')

    # Compute and plot second dendrogram.
    ax2 = fig.add_axes([0.3,0.71,0.6,0.2])
    Y = sch.linkage(D, method='average')
    Z2 = sch.dendrogram(Y, count_sort='ascending')
    ax2.set_xticks([])
    ax2.set_yticks([])
    # http://stackoverflow.com/questions/30009062/get-rid-of-grey-background-in-python-matplotlib-bar-chart
    ax2.set_axis_bgcolor('white')

    # Plot distance matrix.
    axmatrix = fig.add_axes([0.3,0.1,0.6,0.6])
    idx1 = Z1['leaves']
    idx2 = Z2['leaves']
    C = C[idx1,:]
    C = C[:,idx2]
    im = axmatrix.matshow(C, aspect='auto', origin='lower', cmap=pylab.cm.YlGnBu)
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])

    # Sample names
    axmatrix.set_xticks(range(D.shape[0]))
    axmatrix.set_xticklabels(np.asarray(sampleNamesReordered)[idx1], minor=False, fontsize=14, rotation=90)
    axmatrix.xaxis.set_label_position('bottom')
    axmatrix.xaxis.tick_bottom()
    # for tic in axmatrix.xaxis.get_major_ticks():
    #     tic.tick1On = tic.tick2On = False
    for t in axmatrix.xaxis.get_ticklines():
        t.set_visible(False)
    # http://stackoverflow.com/questions/16074392/getting-vertical-gridlines-to-appear-in-line-plot-in-matplotlib
    axmatrix.grid(False)

    # pylab.xticks(rotation=-90, fontsize=14)
    #
    # axmatrix.set_yticks(range(D.shape[0]))
    # axmatrix.set_yticklabels(np.asarray(sampleNamesReordered)[idx2], minor=False)
    # axmatrix.yaxis.set_label_position('right')
    # axmatrix.yaxis.tick_right()
    #
    # axcolor = fig.add_axes([0.94,0.1,0.02,0.6])

    # Plot colorbar.
    axcolor = fig.add_axes([0.91,0.1,0.02,0.6])
    pylab.colorbar(im, cax=axcolor)

    fig.suptitle("Matrix of adjusted Rand Index of TADs with hierarchical clustering\n" + figTitleSuffix, fontsize=20, y=1.01)
    #fig.subplots_adjust(top=1)
    #fig.tight_layout()
    fig.savefig(corrMatrixDirr + "/dendroMatrix.TAD.adjRI." + figNameSuffix + ".pdf", dpi=150, bbox_inches='tight')
    fig.savefig(corrMatrixDirr + "/dendroMatrix.TAD.adjRI." + figNameSuffix + ".png", dpi=150, bbox_inches='tight')
    plt.close()


    ################################################
    # Rand Index of TADs w/ published TADs

    # Add published Domains
    for dataSet in publishedDomainsX.keys():
        # dataSet = 'marks_tsixstop_mesc'
        #print "\n" + dataSet
        #old_i = 0

        # Plot heat maps
        A = dataDict_QN_NaN[whatDataA].copy()
        Alims = [0, A.shape[0]]

        old_i = Alims[0]
        maxLen = Alims[1] - Alims[0] + 1
        TADid = 0
        TADvectors[dataSet] = []
        for dixidx in range(publishedDomainsX[dataSet].shape[0]):
            start = publishedDomainsX[dataSet]['start'].iloc[dixidx]
            startBin = int(start / intTheRez)
            end = publishedDomainsX[dataSet]['end'].iloc[dixidx]
            endBin = int(end / intTheRez)
            if endBin > Alims[1]:
                endBin = Alims[1]
            #print str(start) + "-" + str(end) + " = " + str(startBin) + "-" + str(endBin)
            TADvectors[dataSet].extend( [TADid] * (startBin - old_i) )
            #print(len(TADvectors[dataSet]))
            old_i = np.nanmax((old_i, startBin))
            if (old_i > Alims[1]):
                old_i = Alims[1]
            TADid += 1
            TADvectors[dataSet].extend( [TADid] * (endBin - old_i + 1) )
            #print(len(TADvectors[dataSet]))
            TADid += 1
            old_i = endBin + 1
            if (old_i > Alims[1]):
                old_i = Alims[1]
        #TADvectors[dataSet].extend( [TADid] * (len(xVals) - old_i) )
        TADvectors[dataSet].extend( [TADid] * (len(xVals) - len(TADvectors[dataSet])) )
        #print(len(TADvectors[dataSet]))

    # Remove NaNs
    for dataSet in TADvectors:
        print dataSet + ": " + str(len(TADvectors[dataSet]))
        TADvectors[dataSet] = np.array(TADvectors[dataSet], dtype=float)
        TADvectors[dataSet][nandices] = np.nan

    # TADARI = {}
    # for (a, b) in it.combinations(dataSets,2):
    #     combo = a + "-vs-" + b
    #     #print combo
    #     X =  TADvectors[a]
    #     Y =  TADvectors[b]
    #     TADARI[combo] = adjusted_rand_score(X,Y)

    dataPubSets = sorted(set(TADvectors.keys()))

    TADARI = []
    for dataSetA in dataPubSets:
        for dataSetB in dataPubSets:

            #print str(idxTally) + ": " + dataSetA + " " + dataSetB

            X = TADvectors[dataSetA]
            Y = TADvectors[dataSetB]

            TADARI.append(adjusted_rand_score(X,Y))

    # Reshape into matriz
    TADARI= np.reshape(TADARI, (-1, len(dataPubSets)))

    # Heat map plotting
    C = np.asarray(TADARI)
    D = 1 - C

    sampleNamesReordered = np.asarray(dataPubSets)[0:]

    # Compute and plot first dendrogram.
    fig = pylab.figure(figsize=(10,10))
    ax1 = fig.add_axes([0.09,0.1,0.2,0.6])
    #Y = sch.linkage(D, method='single')
    Y = sch.linkage(D, method='average')
    Z1 = sch.dendrogram(Y, count_sort='ascending', orientation='left')
    ax1.set_xticks([])
    ax1.set_yticks([])
    # http://stackoverflow.com/questions/30009062/get-rid-of-grey-background-in-python-matplotlib-bar-chart
    ax1.set_axis_bgcolor('white')

    # Compute and plot second dendrogram.
    ax2 = fig.add_axes([0.3,0.71,0.6,0.2])
    Y = sch.linkage(D, method='average')
    Z2 = sch.dendrogram(Y, count_sort='ascending')
    ax2.set_xticks([])
    ax2.set_yticks([])
    # http://stackoverflow.com/questions/30009062/get-rid-of-grey-background-in-python-matplotlib-bar-chart
    ax2.set_axis_bgcolor('white')

    # Plot distance matrix.
    axmatrix = fig.add_axes([0.3,0.1,0.6,0.6])
    idx1 = Z1['leaves']
    idx2 = Z2['leaves']
    C = C[idx1,:]
    C = C[:,idx2]
    im = axmatrix.matshow(C, aspect='auto', origin='lower', cmap=pylab.cm.YlGnBu)
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])

    # Sample names
    axmatrix.set_xticks(range(D.shape[0]))
    axmatrix.set_xticklabels(np.asarray(sampleNamesReordered)[idx1], minor=False, fontsize=14, rotation=90)
    axmatrix.xaxis.set_label_position('bottom')
    axmatrix.xaxis.tick_bottom()
    # for tic in axmatrix.xaxis.get_major_ticks():
    #     tic.tick1On = tic.tick2On = False
    for t in axmatrix.xaxis.get_ticklines():
        t.set_visible(False)
    # http://stackoverflow.com/questions/16074392/getting-vertical-gridlines-to-appear-in-line-plot-in-matplotlib
    axmatrix.grid(False)

    # pylab.xticks(rotation=-90, fontsize=14)
    #
    # axmatrix.set_yticks(range(D.shape[0]))
    # axmatrix.set_yticklabels(np.asarray(sampleNamesReordered)[idx2], minor=False)
    # axmatrix.yaxis.set_label_position('right')
    # axmatrix.yaxis.tick_right()
    #
    # axcolor = fig.add_axes([0.94,0.1,0.02,0.6])

    # Plot colorbar.
    axcolor = fig.add_axes([0.91,0.1,0.02,0.6])
    pylab.colorbar(im, cax=axcolor)

    fig.suptitle("Matrix of adjusted Rand Index of TADs with hierarchical clustering including published domains\n" + figTitleSuffix, fontsize=20, y=1.01)
    #fig.subplots_adjust(top=1)
    #fig.tight_layout()
    fig.savefig(corrMatrixDirr + "/dendroMatrix.TAD.plusPublishedTADs.adjRI." + figNameSuffix + ".pdf", dpi=150, bbox_inches='tight')
    fig.savefig(corrMatrixDirr + "/dendroMatrix.TAD.plusPublishedTADs.adjRI." + figNameSuffix + ".png", dpi=150, bbox_inches='tight')
    plt.close()


if 'plotsimilaritymapssansbraindata' in toRun and toRun['plotsimilaritymapssansbraindata']:

    ################################################
    # heat maps with dendrograms
    # *** excluding brain data ***

    noBrainsDataSets = filter(lambda x:'Brain' not in x, dataSets)

    covScorComps = {}

    covScorComps['Pearson'] = []
    covScorComps['Euclidean'] = []
    #covScorComps['meanAbsErr'] = []
    #covScorComps['meanSquErr'] = []

    for dataSetA in noBrainsDataSets:
        for dataSetB in noBrainsDataSets:

            #print str(idxTally) + ": " + dataSetA + " " + dataSetB

            scorzA = covScoresZscored[dataSetA].copy()
            scorzB = covScoresZscored[dataSetB].copy()

            # Remove mutually unscorable regions
            AnanBoolIdx = np.logical_not(np.isnan(scorzA))
            BnanBoolIdx = np.logical_not(np.isnan(scorzB))
            scorzAnn = scorzA[AnanBoolIdx & BnanBoolIdx]
            scorzBnn = scorzB[AnanBoolIdx & BnanBoolIdx]

            covScorComps['Pearson'].append(sps.pearsonr(scorzAnn, scorzBnn)[0])

            diff = scorzAnn-scorzBnn
            # distMatrix.append(np.sum(diff**2)**0.5) # euclidean
            covScorComps['Euclidean'].append(np.linalg.norm(diff)) # euclidean

            #covScorComps['meanAbsErr'].append(np.sum(abs(diff))/len(scorzAnn))  # Mean absolute error
            #covScorComps['meanSquErr'].append(np.sum(diff**2)/len(scorzAnn)) # Mean Squared error

    # Reshape into matrices
    for metricType in covScorComps:
        #print metricType

        covScorComps[metricType] = np.reshape(covScorComps[metricType], (-1, len(noBrainsDataSets)))
        covScorComps[metricType] = pd.DataFrame(covScorComps[metricType], index=noBrainsDataSets, columns=noBrainsDataSets)

    # Heat map plotting
    for metricType in covScorComps:
        #print metricType

        C = np.asarray(covScorComps[metricType])

        # # Try to order in the way you want the leaves sorted
        # # This does not appear to work
        # C = C[[0,3,1,2],:]
        # C = C[:,[0,3,1,2]]
        if metricType == 'Pearson':
            D = 1 - C
        else:
            D = C

        # sampleNames = np.asarray(corrMatrixList['chromosome'].columns.values.tolist())[[0,3,1,2]]
        # sampleNames = np.asarray(['WT Xa', 'Brain Xi', 'WT Xi', 'Brain Xa' ])[[0,3,1,2]]
        sampleNamesReordered = np.asarray(noBrainsDataSets)[0:]

        # Compute and plot first dendrogram.
        fig = pylab.figure(figsize=(10,10))
        ax1 = fig.add_axes([0.09,0.1,0.2,0.6])
        #Y = sch.linkage(D, method='single')
        Y = sch.linkage(D, method='average')
        Z1 = sch.dendrogram(Y, count_sort='ascending', orientation='left')
        ax1.set_xticks([])
        ax1.set_yticks([])
        # http://stackoverflow.com/questions/30009062/get-rid-of-grey-background-in-python-matplotlib-bar-chart
        ax1.set_axis_bgcolor('white')

        # Compute and plot second dendrogram.
        ax2 = fig.add_axes([0.3,0.71,0.6,0.2])
        Y = sch.linkage(D, method='average')
        Z2 = sch.dendrogram(Y, count_sort='ascending')
        ax2.set_xticks([])
        ax2.set_yticks([])
        # http://stackoverflow.com/questions/30009062/get-rid-of-grey-background-in-python-matplotlib-bar-chart
        ax2.set_axis_bgcolor('white')

        # Plot distance matrix.
        axmatrix = fig.add_axes([0.3,0.1,0.6,0.6])
        idx1 = Z1['leaves']
        idx2 = Z2['leaves']
        C = C[idx1,:]
        C = C[:,idx2]
        im = axmatrix.matshow(C, aspect='auto', origin='lower', cmap=pylab.cm.YlGnBu)
        axmatrix.set_xticks([])
        axmatrix.set_yticks([])

        # Sample names
        axmatrix.set_xticks(range(D.shape[0]))
        axmatrix.set_xticklabels(np.asarray(sampleNamesReordered)[idx1], minor=False, fontsize=14, rotation=90)
        axmatrix.xaxis.set_label_position('bottom')
        axmatrix.xaxis.tick_bottom()
        # for tic in axmatrix.xaxis.get_major_ticks():
        #     tic.tick1On = tic.tick2On = False
        for t in axmatrix.xaxis.get_ticklines():
            t.set_visible(False)
        # http://stackoverflow.com/questions/16074392/getting-vertical-gridlines-to-appear-in-line-plot-in-matplotlib
        axmatrix.grid(False)

        # pylab.xticks(rotation=-90, fontsize=14)
        #
        # axmatrix.set_yticks(range(D.shape[0]))
        # axmatrix.set_yticklabels(np.asarray(sampleNamesReordered)[idx2], minor=False)
        # axmatrix.yaxis.set_label_position('right')
        # axmatrix.yaxis.tick_right()
        #
        # axcolor = fig.add_axes([0.94,0.1,0.02,0.6])

        # Plot colorbar.
        axcolor = fig.add_axes([0.91,0.1,0.02,0.6])
        pylab.colorbar(im, cax=axcolor)

        fig.suptitle(metricType + " matrix with hierarchical clustering using coverage Z-scores\n" + figNameSuffix.replace('BrainVs', ''), fontsize=20, y=1.01)
        #fig.subplots_adjust(top=1)
        #fig.tight_layout()
        fig.savefig(corrMatrixDirr + "/dendroMatrix." + metricType + ".coverageZscores." + figNameSuffix.replace('BrainVs', '') + ".pdf", dpi=150, bbox_inches='tight')
        fig.savefig(corrMatrixDirr + "/dendroMatrix." + metricType + ".coverageZscores." + figNameSuffix.replace('BrainVs', '') + ".png", dpi=150, bbox_inches='tight')
        plt.close()


    ################################################
    # Rand Index of TADs
    # *** excluding brain data ***

    TADARI = []
    for dataSetA in noBrainsDataSets:
        for dataSetB in noBrainsDataSets:

            #print str(idxTally) + ": " + dataSetA + " " + dataSetB

            X = TADvectors[dataSetA]
            Y = TADvectors[dataSetB]

            TADARI.append(adjusted_rand_score(X,Y))

    # Reshape into matriz
    TADARI= np.reshape(TADARI, (-1, len(noBrainsDataSets)))

    # Heat map plotting
    C = np.asarray(TADARI)
    D = 1 - C

    sampleNamesReordered = np.asarray(noBrainsDataSets)[0:]

    # Compute and plot first dendrogram.
    fig = pylab.figure(figsize=(10,10))
    ax1 = fig.add_axes([0.09,0.1,0.2,0.6])
    #Y = sch.linkage(D, method='single')
    Y = sch.linkage(D, method='average')
    Z1 = sch.dendrogram(Y, count_sort='ascending', orientation='left')
    ax1.set_xticks([])
    ax1.set_yticks([])
    # http://stackoverflow.com/questions/30009062/get-rid-of-grey-background-in-python-matplotlib-bar-chart
    ax1.set_axis_bgcolor('white')

    # Compute and plot second dendrogram.
    ax2 = fig.add_axes([0.3,0.71,0.6,0.2])
    Y = sch.linkage(D, method='average')
    Z2 = sch.dendrogram(Y, count_sort='ascending')
    ax2.set_xticks([])
    ax2.set_yticks([])
    # http://stackoverflow.com/questions/30009062/get-rid-of-grey-background-in-python-matplotlib-bar-chart
    ax2.set_axis_bgcolor('white')

    # Plot distance matrix.
    axmatrix = fig.add_axes([0.3,0.1,0.6,0.6])
    idx1 = Z1['leaves']
    idx2 = Z2['leaves']
    C = C[idx1,:]
    C = C[:,idx2]
    im = axmatrix.matshow(C, aspect='auto', origin='lower', cmap=pylab.cm.YlGnBu)
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])

    # Sample names
    axmatrix.set_xticks(range(D.shape[0]))
    axmatrix.set_xticklabels(np.asarray(sampleNamesReordered)[idx1], minor=False, fontsize=14, rotation=90)
    axmatrix.xaxis.set_label_position('bottom')
    axmatrix.xaxis.tick_bottom()
    # for tic in axmatrix.xaxis.get_major_ticks():
    #     tic.tick1On = tic.tick2On = False
    for t in axmatrix.xaxis.get_ticklines():
        t.set_visible(False)
    # http://stackoverflow.com/questions/16074392/getting-vertical-gridlines-to-appear-in-line-plot-in-matplotlib
    axmatrix.grid(False)

    # pylab.xticks(rotation=-90, fontsize=14)
    #
    # axmatrix.set_yticks(range(D.shape[0]))
    # axmatrix.set_yticklabels(np.asarray(sampleNamesReordered)[idx2], minor=False)
    # axmatrix.yaxis.set_label_position('right')
    # axmatrix.yaxis.tick_right()
    #
    # axcolor = fig.add_axes([0.94,0.1,0.02,0.6])

    # Plot colorbar.
    axcolor = fig.add_axes([0.91,0.1,0.02,0.6])
    pylab.colorbar(im, cax=axcolor)

    fig.suptitle("Matrix of adjusted Rand Index of TADs with hierarchical clustering\n" + figNameSuffix.replace('BrainVs', ''), fontsize=20, y=1.01)
    #fig.subplots_adjust(top=1)
    #fig.tight_layout()
    fig.savefig(corrMatrixDirr + "/dendroMatrix.TAD.withNans.adjRI." + figNameSuffix.replace('BrainVs', '') + ".pdf", dpi=150, bbox_inches='tight')
    fig.savefig(corrMatrixDirr + "/dendroMatrix.TAD.withNans.adjRI." + figNameSuffix.replace('BrainVs', '') + ".png", dpi=150, bbox_inches='tight')
    plt.close()


    ################################################
    # Rand Index of TADs after removing NaNs
    # *** excluding brain data ***

    TADARI = []
    for dataSetA in noBrainsDataSets:
        for dataSetB in noBrainsDataSets:

            #print str(idxTally) + ": " + dataSetA + " " + dataSetB

            X = TADvectors[dataSetA]
            Y = TADvectors[dataSetB]

            TADARI.append(adjusted_rand_score(X,Y))

    # Reshape into matriz
    TADARI= np.reshape(TADARI, (-1, len(noBrainsDataSets)))

    # Heat map plotting
    C = np.asarray(TADARI)
    D = 1 - C

    sampleNamesReordered = np.asarray(noBrainsDataSets)[0:]

    # Compute and plot first dendrogram.
    fig = pylab.figure(figsize=(10,10))
    ax1 = fig.add_axes([0.09,0.1,0.2,0.6])
    #Y = sch.linkage(D, method='single')
    Y = sch.linkage(D, method='average')
    Z1 = sch.dendrogram(Y, count_sort='ascending', orientation='left')
    ax1.set_xticks([])
    ax1.set_yticks([])
    # http://stackoverflow.com/questions/30009062/get-rid-of-grey-background-in-python-matplotlib-bar-chart
    ax1.set_axis_bgcolor('white')

    # Compute and plot second dendrogram.
    ax2 = fig.add_axes([0.3,0.71,0.6,0.2])
    Y = sch.linkage(D, method='average')
    Z2 = sch.dendrogram(Y, count_sort='ascending')
    ax2.set_xticks([])
    ax2.set_yticks([])
    # http://stackoverflow.com/questions/30009062/get-rid-of-grey-background-in-python-matplotlib-bar-chart
    ax2.set_axis_bgcolor('white')

    # Plot distance matrix.
    axmatrix = fig.add_axes([0.3,0.1,0.6,0.6])
    idx1 = Z1['leaves']
    idx2 = Z2['leaves']
    C = C[idx1,:]
    C = C[:,idx2]
    im = axmatrix.matshow(C, aspect='auto', origin='lower', cmap=pylab.cm.YlGnBu)
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])

    # Sample names
    axmatrix.set_xticks(range(D.shape[0]))
    axmatrix.set_xticklabels(np.asarray(sampleNamesReordered)[idx1], minor=False, fontsize=14, rotation=90)
    axmatrix.xaxis.set_label_position('bottom')
    axmatrix.xaxis.tick_bottom()
    # for tic in axmatrix.xaxis.get_major_ticks():
    #     tic.tick1On = tic.tick2On = False
    for t in axmatrix.xaxis.get_ticklines():
        t.set_visible(False)
    # http://stackoverflow.com/questions/16074392/getting-vertical-gridlines-to-appear-in-line-plot-in-matplotlib
    axmatrix.grid(False)

    # pylab.xticks(rotation=-90, fontsize=14)
    #
    # axmatrix.set_yticks(range(D.shape[0]))
    # axmatrix.set_yticklabels(np.asarray(sampleNamesReordered)[idx2], minor=False)
    # axmatrix.yaxis.set_label_position('right')
    # axmatrix.yaxis.tick_right()
    #
    # axcolor = fig.add_axes([0.94,0.1,0.02,0.6])

    # Plot colorbar.
    axcolor = fig.add_axes([0.91,0.1,0.02,0.6])
    pylab.colorbar(im, cax=axcolor)

    fig.suptitle("Matrix of adjusted Rand Index of TADs with hierarchical clustering\n" + figNameSuffix.replace('BrainVs', ''), fontsize=20, y=1.01)
    #fig.subplots_adjust(top=1)
    #fig.tight_layout()
    fig.savefig(corrMatrixDirr + "/dendroMatrix.TAD.adjRI." + figNameSuffix.replace('BrainVs', '') + ".pdf", dpi=150, bbox_inches='tight')
    fig.savefig(corrMatrixDirr + "/dendroMatrix.TAD.adjRI." + figNameSuffix.replace('BrainVs', '') + ".png", dpi=150, bbox_inches='tight')
    plt.close()


if 'differentialanalysis' in toRun and toRun['differentialanalysis']:

    ################################################
    # Differential coverage scores: Xa - Xi

    print "\nDifferential coverage scores\n"

    # Get cell types
    cellTypes = [w.replace('_Xi', '') for w in dataSets]
    cellTypes = [w.replace('_Xa', '') for w in cellTypes]
    cellTypes = list(set(cellTypes))

    deltaIS = {}
    for (a, b) in it.combinations(cellTypes,2):
        combo = b + "-vs-" + a
        #print combo
        for Xstate in ("Xa", "Xi"):
            deltaIS[combo + "_" + Xstate] = covScoresZscored[b + "_" + Xstate] - covScoresZscored[a + "_" + Xstate]

    flatDeltas = []
    for combo in deltaIS.keys():
        # print combo
        flatDeltas.extend(deltaIS[combo])
    # Get 'outlier' values: top bottom
    lowYval = np.nanpercentile(flatDeltas, outThresh)
    highYval = np.nanpercentile(flatDeltas, 100 - outThresh)
    #tZval = np.nanpercentile(abs(np.array(flatDeltas)), 100-outThresh)


    ################################################
    # Top/bot most differential regions

    sigDiffBinsLow = {}
    sigDiffBinsHigh = {}
    sigDiffBins = {}

    for Xstate in ("Xa", "Xi"):

        idxTally = -1
        for combo in deltaIS.keys():
            if Xstate in combo:
                #print combo

                idxTally += 1

                xVals = covScoreLoci
                yVals = deltaIS[combo]

                # Get 'outlier' positions
                lowYvalsBoolIdx = np.less_equal(yVals, lowYval)
                highYvalsBoolIdx = np.greater_equal(yVals, highYval)
                #lowYvalsBoolIdx = np.less_equal(yVals, -tZval)
                #highYvalsBoolIdx = np.greater_equal(yVals, tZval)
                hiloYvalsIdx = lowYvalsBoolIdx | highYvalsBoolIdx
                #print(hiloYvalsIdx)
                #print len(hiloYvalsIdx)
                #print np.size(hiloYvalsIdx)

                sigDiffBinsLow[combo] = lowYvalsBoolIdx
                sigDiffBinsHigh[combo] = highYvalsBoolIdx
                sigDiffBins[combo] = hiloYvalsIdx


    ################################################
    # Plot differential Z-score distribution

    diffDir = outDir + "/diffPlots"
    if not os.path.exists(diffDir):
        os.makedirs(diffDir)
    sns.set(style="darkgrid")

    # Set up the matplotlib figure
    fig = plt.figure(figsize=(10,10))

    #sns.violinplot(y='IS', x='celltype', hue='Xstate', data=alltheData, split=True, palette=sns.color_palette("pastel", 2))
    ax = sns.distplot(list(compress(flatDeltas, ~np.isnan(flatDeltas))), hist=True, rug=False)

    # http://matplotlib.org/examples/pylab_examples/fonts_demo.html
    font0 = FontProperties()
    font1 = font0.copy()
    font1.set_size('12')
    # ax.axvline(-tZval, color='k', linestyle='--', linewidth=1)
    # ax.text(-tZval, 0.5, str(-round(tZval,2)), rotation=90, verticalalignment="top", fontproperties=font1)
    # ax.axvline(tZval, color='k', linestyle='--', linewidth=1)
    # ax.text(tZval, 0.5, str(round(tZval,2)), rotation=90, verticalalignment="top", fontproperties=font1)
    ax.axvline(lowYval, color='k', linestyle='--', linewidth=1)
    ax.text(lowYval, 0.5, str(-round(lowYval,2)), rotation=90, verticalalignment="top", fontproperties=font1)
    ax.axvline(highYval, color='k', linestyle='--', linewidth=1)
    ax.text(highYval, 0.5, str(round(highYval,2)), rotation=90, verticalalignment="top", fontproperties=font1)

    ax.set(xlabel='Differential coverage z-score', ylabel='Density')

    # Finalize the figure
    # sns.despine(left=True, bottom=True)

    # Tidy up and save
    fig.suptitle("Distribution of differential coverage Z-scores\n" + figTitleSuffix, fontsize=12, y=1.01)
    #fig.subplots_adjust(top=1)
    #fig.tight_layout()
    fig.savefig(diffDir + "/differentialcoverageZscores." + figNameSuffix + ".pdf", dpi=150, bbox_inches='tight')
    fig.savefig(diffDir + "/differentialcoverageZscores." + figNameSuffix + ".png", dpi=150, bbox_inches='tight')

    plt.close()


    if 'differentiallineplots' in toRun and toRun['differentiallineplots']:

        ################################################
        # Plot differential Z-score line plots

        diffDir = outDir + "/diffPlots"
        if not os.path.exists(diffDir):
            os.makedirs(diffDir)

        # Y-axis limits
        # minY = -2
        # maxY = 2
        minY = 10**6
        maxY = 0
        if minY > np.nanmin(flatDeltas):
            minY = np.nanmin(flatDeltas)
        if maxY < np.nanmax(flatDeltas):
            maxY = np.nanmax(flatDeltas)
        #minY, maxY = conditionDistanceLimits(minY, maxY)
        # Make it symmetrical, if opposite signs
        if maxY > 0 and minY < 0 :
            maxMinMaxY = np.nanmax([abs(minY), abs(maxY)])
            maxY = maxMinMaxY
            minY = -maxMinMaxY
        minY *= 1.02
        maxY *= 1.02
        # Bounds
        # if minY < -2:
        #     minY = -2
        # if maxY > 2:
        #     maxY = 2

        # X-axis limits
        minX = 0
        maxX = int(np.ceil(float(chromSizes[thechrom]) / 1000000.0))

        for Xstate in ("Xa", "Xi"):

            fig, ax = plt.subplots(1, 1, figsize=(figureWidth, 10))
            plots = {}

            idxTally = -1
            for combo in deltaIS.keys():
                if Xstate in combo:
                    #print combo

                    idxTally += 1

                    xVals = covScoreLoci
                    yVals = deltaIS[combo]

                    plots[combo], = ax.plot(xVals, yVals,
                                              '-', label=combo,
                                              color=sns.color_palette('muted', len(deltaIS.keys()))[idxTally],
                                              linewidth=4, linestyle='solid')

                    # Get 'outlier' positions
                    lowYvalsBoolIdx = np.less_equal(yVals, lowYval)
                    highYvalsBoolIdx = np.greater_equal(yVals, highYval)
                    #lowYvalsBoolIdx = np.less_equal(yVals, -tZval)
                    #highYvalsBoolIdx = np.greater_equal(yVals, tZval)
                    hiloYvalsIdx = lowYvalsBoolIdx | highYvalsBoolIdx
                    # print np.size(hiloYvalsIdx)  / np.size(yVals) * 100.0

                    hiloYvalsIdx = np.where(hiloYvalsIdx)
                    if np.size(hiloYvalsIdx) > 0:
                        ax.plot(xVals[hiloYvalsIdx], yVals[hiloYvalsIdx], 'o', markersize=12, markeredgewidth=0,
                                markeredgecolor='black', markerfacecolor=sns.color_palette('bright', len(deltaIS.keys()))[idxTally],
                                label=combo + " topbot" + str(outThresh) + "%")

                    ax.grid(True)
                    ax.set_ylabel('Delta log2(coverage score / mean score)', fontsize=14)
                    ax.set_xlabel('Mb', fontsize=14)
                    for tick in ax.xaxis.get_major_ticks():
                        tick.label.set_fontsize(12)
                    for tick in ax.yaxis.get_major_ticks():
                        tick.label.set_fontsize(12)

            for aloi in theLoi:
                #print aloi
                #print LOI[aloi]

                chrCoords = LOI[aloi].split(":")
                theChr = chrCoords[0]
                if theChr == thechrom:
                    #print "\tPlot it!"
                    startStop = chrCoords[1].split("-")
                    startStop = map(int, startStop)
                    #startStop = [round(x / float(theRez)) for x in startStop]
                    geneMidPoint = (startStop[0] + (startStop[1] - startStop[0]))/10.0**6

                    ax.axvline(geneMidPoint, color='k', linestyle='--', linewidth=0.5)
                    # http://matplotlib.org/examples/pylab_examples/fonts_demo.html
                    font0 = FontProperties()
                    font1 = font0.copy()
                    font1.set_size('12')
                    ax.text(geneMidPoint, maxY, aloi, rotation=90, verticalalignment="top", fontproperties=font1)

            #startStop = [round(x / float(theRez)) for x in startStop]
            #ax.set_title(aloi + "\n(" + chrCoords[0] + ":" + sizeof_fmt(startStop[0]) + ")", fontsize=14)

            ax.set_xlim(minX, maxX)
            ax.set_ylim(minY, maxY)
            #ax.legend(loc='best', fontsize=14)
            ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                   ncol=2, mode="expand", borderaxespad=0., fontsize=14)

            # Tidy up and save
            fig.suptitle(Xstate + "Differential coverageScores\n" + figTitleSuffix, fontsize=20, y=1.01)
            #fig.subplots_adjust(top=1)
            #fig.tight_layout()
            fig.savefig(diffDir + "/diffCoverageScores." + Xstate + "." + figNameSuffix + ".pdf", dpi=150, bbox_inches='tight')
            fig.savefig(diffDir + "/diffCoverageScores." + Xstate + "." + figNameSuffix + ".png", dpi=150, bbox_inches='tight')

            if 'plotdiffplotzooms' in toRun and toRun['plotdiffplotzooms']:
                # Zoom around LOI

                diffDirZ = diffDir + "/Zoomed"
                if not os.path.exists(diffDirZ):
                    os.makedirs(diffDirZ)

                for aloi in theLoi:
                    #print aloi
                    #print LOI[aloi]

                    chrCoords = LOI[aloi].split(":")
                    theChr = chrCoords[0]
                    if theChr == thechrom:
                        #print "\tPlot it!"

                        startStop = chrCoords[1].split("-")
                        startStop = map(int, startStop)
                        #startStop = [round(x / float(theRez)) for x in startStop]
                        # theLocus = startStop[0]
                        # #print str(theLocus)
                        geneMidPoint = (startStop[0] + (startStop[1] - startStop[0]))/10.0**6
                        plusMinusDist = round((zoomWin / 2.0) / 1000000.0)

                        # Zoom in
                        ax.set_xlim(max(minX, geneMidPoint - plusMinusDist),
                                     min(maxX, geneMidPoint + plusMinusDist))

                        ax.legend(loc='best', fontsize=14)
                        # ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                        #        ncol=2, mode="expand", borderaxespad=0., fontsize=14)

                        fig.set_size_inches(20, 10)

                        # Tidy up and save
                        fig.suptitle(Xstate + " differential coverageScores\n" \
                                + aloi + " @ " + str(round(geneMidPoint, 1)) + "Mb; " \
                                + " zoom" + str(round(zoomWin/1000000.0, 1)) + "Mb; \n" \
                                + figTitleSuffix, fontsize=20, y=1.1)
                        #fig.subplots_adjust(top=1)
                        fig.tight_layout()
                        fig.savefig(diffDirZ + "/diffCoverageScores." + Xstate + "."
                                    + aloi + "." + "zoom" + str(round(zoomWin/1000000.0,1)) + "Mb." \
                                    + figNameSuffix + ".pdf" , dpi=150,
                                    #bbox_inches='tight', )
                                    bbox_inches=0, orientation='landscape', pad_inches=0)
                        fig.savefig(diffDirZ + "/diffCoverageScores." + Xstate + "."
                                    + aloi + "." + "zoom" + str(round(zoomWin/1000000.0,1)) + "Mb." \
                                    + figNameSuffix + ".png" , dpi=150,
                                    #bbox_inches='tight', )
                                    bbox_inches=0, orientation='landscape', pad_inches=0)

            plt.close()


    if 'gendiffbedgraphs' in toRun and toRun['gendiffbedgraphs']:

        ################################################
        #  Differential bedgraphs and beds of topBottom 5%

        print "\nDifferential bedgraphs and beds\n"

        bedgraphDir = outDir + "/diffBedgraphs"
        if not os.path.exists(bedgraphDir):
            os.makedirs(bedgraphDir)

        bedDir = outDir + "/diffBeds"
        if not os.path.exists(bedDir):
            os.makedirs(bedDir)

        for Xstate in ("Xa", "Xi"):

            idxTally = -1
            for combo in deltaIS.keys():
                if Xstate in combo:
                    #print combo

                    idxTally += 1

                    xVals = covScoreLoci
                    yVals = deltaIS[combo]

                    ##########
                    # plot bedgraphs

                    bedgraphFileName = bedgraphDir + "/diffCoverageScores." + combo + "." + bedNameSuffix + ".bedGraph.gz"
                    f = gzip.open(bedgraphFileName, "wb")
                    for i in range(len(xVals)):
                        if ~np.isnan(yVals[i]):
                            f.write(
                                thechrom + "\t" + \
                                str(int(round(xVals[i] * 1.0*10**6) - (intTheRez/2))) + "\t" + \
                                str(int(round(xVals[i] * 1.0*10**6) + (intTheRez/2))) + "\t" + \
                                str(yVals[i]) + "\n")
                    f.close()

                    ##########
                    # plot co-ordinates of outlier values and save as BED files
                    hiloYvalsIdx = np.where(sigDiffBins[combo])
                    if np.size(hiloYvalsIdx) > 0:
                        bedFileName = bedDir + "/" + "topBot" + str(outThresh) + "pct." \
                                      + "diffCoverageScores." \
                                      + combo + "." + bedNameSuffix + ".bed.gz"
                        f = gzip.open(bedFileName, "wb")

                        # yVals = np.linspace(1, 100, num=100)
                        # xVals = yVals

                        old_i = np.nanmin(hiloYvalsIdx) - 1
                        start = str(int((xVals[np.nanmin(hiloYvalsIdx)] * 1.0*10**6) - (intTheRez/2)))
                        for i in np.nditer(hiloYvalsIdx):
                            if i != (old_i + 1) or i == np.nanmax(hiloYvalsIdx):
                                f.write(
                                    thechrom + "\t" + start + "\t" + \
                                    str(int(round(xVals[old_i] * 1.0*10**6) + (intTheRez/2))) + "\n")
                                start = str(int(round(xVals[i] * 1.0*10**6) - (intTheRez/2)))
                            old_i = i
                        f.close()


    if 'plotdiffpairplots' in toRun and toRun['plotdiffpairplots']:

        ################################################
        # Distributions as pair plots
        # with most differential Xi Del - WT

        print "\nDistributions of as pair plots with differential regions\n"

        distDir = outDir + "/distributions"
        if not os.path.exists(distDir):
            os.makedirs(distDir)

        # Y-axis limits
        # minY = -2
        # maxY = 2
        minY = np.inf
        maxY = -np.inf
        for dataSet in dataSets:
            # dataSet = 'pooledWT_alt_IC'
            #print dataSet
            if minY > np.nanmin(covScoresZscored[dataSet]):
                minY = np.nanmin(covScoresZscored[dataSet])
            if maxY < np.nanmax(covScoresZscored[dataSet]):
                maxY = np.nanmax(covScoresZscored[dataSet])
        #minY, maxY = conditionDistanceLimits(minY, maxY)
        # Make it symmetrical, if opposite signs
        if maxY > 0 and minY < 0 :
             maxMinMaxY = np.nanmax([abs(minY), abs(maxY)])
             maxY = maxMinMaxY
             minY = -maxMinMaxY
        minY *= 1.02
        maxY *= 1.02
        # Bounds
        #if minY < -4:
        #    minY = -4
        #if maxY > 4:
        #    maxY = 4

        # Denanify
        # Done above now.
        # allTheNaNsBoolIdx = [False] * covScoresZscored[covScoresZscored.keys()[0]].shape[0]
        # for dataSet in dataSets:
        #     allTheNaNsBoolIdx = allTheNaNsBoolIdx | np.isnan(covScoresZscored[dataSet])
        covScoresZscoredDeNaned = {}
        for dataSet in dataSets:
            covScoresZscoredDeNaned[dataSet] = covScoresZscored[dataSet][~allTheNaNsBoolIdx]
        covScoreLociDeNaned = covScoreLoci[~allTheNaNsBoolIdx]
        covScoreLociRTDdeNaned = covScoreLociRTD[~allTheNaNsBoolIdx]
        # covScoreLociRTDclassesDeNaned = covScoreLociRTDclasses[~allTheNaNsBoolIdx]
        #       TypeError: only integer arrays with one element can be converted to an index
        covScoreLociRTDclassesDeNaned = list(compress(covScoreLociRTDclasses, ~allTheNaNsBoolIdx))

        #for combo in sigDiffBins.keys():
        combo = "Del-vs-WT_Xi"
        HiBins = sigDiffBinsHigh[combo][~allTheNaNsBoolIdx]
        LoBins = sigDiffBinsLow[combo][~allTheNaNsBoolIdx]
        DiffBinType = ["Intermediate"] * len(HiBins)
        if len(np.where(HiBins)[0]) > 0:
            for i in np.nditer(np.where(HiBins)):
                #print i
                DiffBinType[i] = 'High'
        if len(np.where(LoBins)[0]) > 0:
            for i in np.nditer(np.where(LoBins)):
                #print i
                DiffBinType[i] = 'Low'

        for Xtype in ['Xa', 'Xi']:
            #print(Xtype)
            XtypeString = '_' + Xtype
            whichData = [XtypeString in i for i in dataSets]
            np.where(whichData)
            dataSubSets = list(compress(dataSets, whichData))

            diagcolors = []
            for dataSet in dataSubSets:
                diagcolors.append(contactMapIOmetadata[dataSet]['dataColor'])

            theData = pd.DataFrame(covScoresZscoredDeNaned)[dataSubSets]

            ###############################
            # Grouped by differential category
            theData[combo] = DiffBinType

            sns.set(style="darkgrid")
            plt.figure(figsize=(10,10))

            # Draw a violinplot with a narrower bandwidth than the default
            # sns.pairplot(data=altRatioData, kind="reg", markers="+", size=3)
            # sns.pairplot(data=altRatioData, size=3, diag_kind="kde")
            g = sns.PairGrid(data=theData, hue=combo, palette=sns.diverging_palette(255, 133, l=60, n=5, center="dark"))
            g = g.map_offdiag(plt.scatter, edgecolor="w", s=35)
            g = g.map_diag(plt.hist, edgecolor="w",
                           bins=np.arange(minY, maxY, 0.5))
            #g = g.map_diag(sns.kdeplot, lw=2, legend=False)
            g = g.set(xlim=(minY, maxY),
                      ylim=(minY, maxY))
            g.add_legend()

            # import pprint
            # pp = pprint.PrettyPrinter(indent=4)
            # for ax in g.diag_axes:
            #     #pp.pprint(vars(ax))
            #     diag_axes_axes = getattr(ax, '_axes')
            #     #pp.pprint(vars(diag_axes_axes))
            #     daGraph = ax.figure
            #     daGraph.set_facecolor('red')

            # Tidy up and save
            g.fig.suptitle(Xtype + " pair grid of coverage Z-score distributions by differential category\n" + figTitleSuffix +
            #               "\n(n=" + str(covScoresZscored[dataSet].shape[0]) + ")", fontsize=12, y=1.01)
                           "; n=" + str((~allTheNaNsBoolIdx).sum()), fontsize=12, y=1.05)
            #fig.subplots_adjust(top=1)
            #g.fig.tight_layout()
            g.fig.savefig(distDir + "/pairGridByDiffCat.coverageZscores." + Xtype + "." + figNameSuffix + ".pdf", dpi=150, bbox_inches='tight')
            g.fig.savefig(distDir + "/pairGridByDiffCat.coverageZscores." + Xtype + "." + figNameSuffix + ".png", dpi=150, bbox_inches='tight')

            plt.close()


    if 'delwtdifferentialenrichmentanalysis' in toRun and toRun['delwtdifferentialenrichmentanalysis']:

        ################################################
        # CTCF in Differential IS regions

        print "\nCTCF analysis\n"

        #for combo in sigDiffBins.keys():
        combo = "Del-vs-WT_Xi"

        CTCFsigDiffBinsLow = {}
        CTCFsigDiffBinsHigh = {}

        resultTable = []
        for whichCTCF in CTCFpeaks.keys():
            # whichCTCF = 'ctcf_del_xi'
            # print "\n" + whichCTCF

            # High: Del - WT greatest
            HighIdx = np.where(sigDiffBinsHigh[combo])
            diffLoci = covScoreLoci[HighIdx]
            CTCFsigDiffBinsHigh[whichCTCF] = [False] * np.size(CTCFpeaks[whichCTCF]['chr'])
            for diffLocus in diffLoci:
                startPos = (diffLocus * 10**6) - intTheRez/2
                endPos = (diffLocus * 10**6) + intTheRez/2
                overlapBool = np.bitwise_or(
                    np.bitwise_and(CTCFpeaks[whichCTCF]['start'] > startPos,  CTCFpeaks[whichCTCF]['start'] < endPos),
                    np.bitwise_and(CTCFpeaks[whichCTCF]['end'] < endPos,  CTCFpeaks[whichCTCF]['end'] > startPos))
                CTCFsigDiffBinsHigh[whichCTCF] = np.bitwise_or(CTCFsigDiffBinsHigh[whichCTCF], overlapBool)
            resultVector = [whichCTCF, "BigPosDiff",
                                          np.sum(CTCFsigDiffBinsHigh[whichCTCF]),
                                          np.sum(~CTCFsigDiffBinsHigh[whichCTCF]),
                                          np.sum(CTCFsigDiffBinsHigh[whichCTCF])/(np.size(CTCFsigDiffBinsHigh[whichCTCF]) *1.0),
                                          np.size(HighIdx)/(np.size(covScoreLoci) *1.0),
                                          (np.sum(CTCFsigDiffBinsHigh[whichCTCF])/(np.size(CTCFsigDiffBinsHigh[whichCTCF]) *1.0)) / \
                                          (np.size(HighIdx)/(np.size(covScoreLoci) *1.0)) ]
            resultTable.append(resultVector)

            # Low
            LowIdx = np.where(sigDiffBinsLow[combo])
            diffLoci = covScoreLoci[LowIdx]
            CTCFsigDiffBinsLow[whichCTCF] = [False] * np.size(CTCFpeaks[whichCTCF]['chr'])
            for diffLocus in diffLoci:
                startPos = (diffLocus * 10**6) - intTheRez/2
                endPos = (diffLocus * 10**6) + intTheRez/2
                overlapBool = np.bitwise_or(
                    np.bitwise_and(CTCFpeaks[whichCTCF]['start'] > startPos,  CTCFpeaks[whichCTCF]['start'] < endPos),
                    np.bitwise_and(CTCFpeaks[whichCTCF]['end'] < endPos,  CTCFpeaks[whichCTCF]['end'] > startPos))
                CTCFsigDiffBinsLow[whichCTCF] = np.bitwise_or(CTCFsigDiffBinsLow[whichCTCF], overlapBool)
            resultVector = [ whichCTCF, "BigNegDiff",
                             np.sum(CTCFsigDiffBinsLow[whichCTCF]),
                             np.sum(~CTCFsigDiffBinsLow[whichCTCF]),
                             np.sum(CTCFsigDiffBinsLow[whichCTCF])/(np.size(CTCFsigDiffBinsLow[whichCTCF]) *1.0),
                             np.size(LowIdx)/(np.size(covScoreLoci) *1.0),
                             (np.sum(CTCFsigDiffBinsLow[whichCTCF])/(np.size(CTCFsigDiffBinsLow[whichCTCF]) *1.0)) / \
                             (np.size(LowIdx)/(np.size(covScoreLoci) *1.0)) ]
            resultTable.append(resultVector)

        resultTable = pd.DataFrame(resultTable, columns=('CTCFdataset', 'diffCategory',
                                        'CTCFpeaksIn', 'CTCFpeaksOut',
                                        'observedProp', 'expectedProp',
                                        'O/E'))
        print(resultTable)

        ##########
        # Barplots

        #####
        # Counts

        diffDir = outDir + "/diffPlots"
        if not os.path.exists(diffDir):
            os.makedirs(diffDir)

        sns.set(style="darkgrid")
        fig = plt.figure(figsize=(8,8))

        #ax = sns.factorplot(y="CTCFpeaksIn", hue="diffCategory", x="CTCFdataset", data=resultTable, kind='bar')
        ax = sns.barplot(y="CTCFpeaksIn", hue="diffCategory", x="CTCFdataset", data=resultTable)

        ax.set_ylabel('CTCF peaks within differential regions', fontsize=10)
        ax.set_xlabel('CTCF ChIP-seq data set', fontsize=10)

        # Tidy up and save
        fig.suptitle("CTCF peaks within Del-WT Xi differential regions \n" + figTitleSuffix +
        #               "\n(n=" + str(covScoresZscored[dataSet].shape[0]) + ")", fontsize=12, y=1.01)
                       "; n=" + str((~allTheNaNsBoolIdx).sum()), fontsize=12, y=1.01)
        #fig.subplots_adjust(top=1)
        #g.fig.tight_layout()
        fig.savefig(diffDir + "/barplot.CTCFpeaks.Del-WTXiDifferentialRegions." + figNameSuffix + ".pdf", dpi=150, bbox_inches='tight')
        fig.savefig(diffDir + "/barplot.CTCFpeaks.Del-WTXiDifferentialRegions." + figNameSuffix + ".png", dpi=150, bbox_inches='tight')

        plt.close()

        #####
        # O/E

        sns.set(style="darkgrid")
        fig = plt.figure(figsize=(8,8))

        #ax = sns.factorplot(y="CTCFpeaksIn", hue="diffCategory", x="CTCFdataset", data=resultTable, kind='bar')
        ax = sns.barplot(y="O/E", hue="diffCategory", x="CTCFdataset", data=resultTable)

        ax.set_ylabel('O/E CTCF peaks within differential regions', fontsize=10)
        ax.set_xlabel('CTCF ChIP-seq data set', fontsize=10)

        # Tidy up and save
        fig.suptitle("O/E CTCF peaks within Del-WT Xi differential regions \n" + figTitleSuffix +
        #               "\n(n=" + str(covScoresZscored[dataSet].shape[0]) + ")", fontsize=12, y=1.01)
                       "; n=" + str((~allTheNaNsBoolIdx).sum()), fontsize=12, y=1.01)
        #fig.subplots_adjust(top=1)
        #g.fig.tight_layout()
        fig.savefig(diffDir + "/barplot.CTCFpeaks.OE.Del-WTXiDifferentialRegions." + figNameSuffix + ".pdf", dpi=150, bbox_inches='tight')
        fig.savefig(diffDir + "/barplot.CTCFpeaks.OE.Del-WTXiDifferentialRegions." + figNameSuffix + ".png", dpi=150, bbox_inches='tight')

        plt.close()


        ################################################
        # Expression in differential IS regions

        print "\nExpression analysis\n"

        #for combo in sigDiffBins.keys():
        combo = "Del-vs-WT_Xi"

        TPMsigDiffBinsHigh = {}
        TPMsigDiffBinsLow = {}

        resultTable = []
        for whichRNA in TPMdata.keys():
            # whichRNA = 'tpm_xi'
            print "\n" + whichRNA

            # High: Del - WT greatest
            HighIdx = np.where(sigDiffBinsHigh[combo])
            diffLoci = covScoreLoci[HighIdx]
            TPMsigDiffBinsHigh[whichRNA] = [False] * np.size(TPMdata[whichRNA]['chr'])
            for diffLocus in diffLoci:
                startPos = (diffLocus * 10**6) - intTheRez/2
                endPos = (diffLocus * 10**6) + intTheRez/2
                overlapBool = np.bitwise_or(
                    np.bitwise_and(TPMdata[whichRNA]['start'] > startPos,  TPMdata[whichRNA]['start'] < endPos),
                    np.bitwise_and(TPMdata[whichRNA]['stop'] < endPos,  TPMdata[whichRNA]['stop'] > startPos))
                TPMsigDiffBinsHigh[whichRNA] = np.bitwise_or(TPMsigDiffBinsHigh[whichRNA], overlapBool)
            resultVector = [whichRNA, "BigPosDiff",
                                          np.sum(TPMsigDiffBinsHigh[whichRNA]),
                                          np.sum(~TPMsigDiffBinsHigh[whichRNA]),
                                          np.sum(TPMsigDiffBinsHigh[whichRNA])/(np.size(TPMsigDiffBinsHigh[whichRNA]) *1.0),
                                          np.size(HighIdx)/(np.size(covScoreLoci) *1.0),
                                          (np.sum(TPMsigDiffBinsHigh[whichRNA])/(np.size(TPMsigDiffBinsHigh[whichRNA]) *1.0)) / \
                                          (np.size(HighIdx)/(np.size(covScoreLoci) *1.0)) ]
            resultTable.append(resultVector)

            # Low
            LowIdx = np.where(sigDiffBinsLow[combo])
            diffLoci = covScoreLoci[LowIdx]
            TPMsigDiffBinsLow[whichRNA] = [False] * np.size(TPMdata[whichRNA]['chr'])
            for diffLocus in diffLoci:
                startPos = (diffLocus * 10**6) - intTheRez/2
                endPos = (diffLocus * 10**6) + intTheRez/2
                overlapBool = np.bitwise_or(
                    np.bitwise_and(TPMdata[whichRNA]['start'] > startPos,  TPMdata[whichRNA]['start'] < endPos),
                    np.bitwise_and(TPMdata[whichRNA]['stop'] < endPos,  TPMdata[whichRNA]['stop'] > startPos))
                TPMsigDiffBinsLow[whichRNA] = np.bitwise_or(TPMsigDiffBinsLow[whichRNA], overlapBool)
            resultVector = [whichRNA, "BigNegDiff",
                                          np.sum(TPMsigDiffBinsLow[whichRNA]),
                                          np.sum(~TPMsigDiffBinsLow[whichRNA]),
                                          np.sum(TPMsigDiffBinsLow[whichRNA])/(np.size(TPMsigDiffBinsLow[whichRNA]) *1.0),
                                          np.size(LowIdx)/(np.size(covScoreLoci) *1.0),
                                          (np.sum(TPMsigDiffBinsLow[whichRNA])/(np.size(TPMsigDiffBinsLow[whichRNA]) *1.0)) / \
                                          (np.size(LowIdx)/(np.size(covScoreLoci) *1.0)) ]
            resultTable.append(resultVector)

        resultTable = pd.DataFrame(resultTable, columns=('TPMdataset', 'diffCategory',
                                        'genesIn', 'genesOut',
                                        'observedProp', 'expectedProp',
                                        'O/E'))
        print(resultTable)

        ##########
        # Barplots

        ####
        # Counts

        diffDir = outDir + "/diffPlots"
        if not os.path.exists(diffDir):
            os.makedirs(diffDir)

        sns.set(style="darkgrid")
        fig = plt.figure(figsize=(6,6))

        ax = sns.barplot(y="genesIn", x="diffCategory", data=resultTable[resultTable['TPMdataset'] == 'tpm_xi'])

        ax.set_ylabel('Genes within differential regions', fontsize=10)
        ax.set_xlabel('RNA-seq data set', fontsize=10)

        # Tidy up and save
        fig.suptitle("Genes within Del-WT Xi differential regions \n" + figTitleSuffix +
        #               "\n(n=" + str(covScoresZscored[dataSet].shape[0]) + ")", fontsize=12, y=1.01)
                       "; n=" + str((~allTheNaNsBoolIdx).sum()), fontsize=12, y=1.01)
        #fig.subplots_adjust(top=1)
        #g.fig.tight_layout()
        fig.savefig(diffDir + "/barplot.genes.Del-WTXiDifferentialRegions." + figNameSuffix + ".pdf", dpi=150, bbox_inches='tight')
        fig.savefig(diffDir + "/barplot.genes.Del-WTXiDifferentialRegions." + figNameSuffix + ".png", dpi=150, bbox_inches='tight')

        plt.close()

        #####
        # O/E

        sns.set(style="darkgrid")
        fig = plt.figure(figsize=(6,6))

        ax = sns.barplot(y="O/E", x="diffCategory", data=resultTable[resultTable['TPMdataset'] == 'tpm_xi'])

        ax.set_ylabel('O/E genes within differential regions', fontsize=10)
        ax.set_xlabel('RNA-seq data set', fontsize=10)

        # Tidy up and save
        fig.suptitle("O/E genes within Del-WT Xi differential regions \n" + figTitleSuffix +
        #               "\n(n=" + str(covScoresZscored[dataSet].shape[0]) + ")", fontsize=12, y=1.01)
                       "; n=" + str((~allTheNaNsBoolIdx).sum()), fontsize=12, y=1.01)
        #fig.subplots_adjust(top=1)
        #g.fig.tight_layout()
        fig.savefig(diffDir + "/barplot.genes.OE.Del-WTXiDifferentialRegions." + figNameSuffix + ".pdf", dpi=150, bbox_inches='tight')
        fig.savefig(diffDir + "/barplot.genes.OE.Del-WTXiDifferentialRegions." + figNameSuffix + ".png", dpi=150, bbox_inches='tight')

        plt.close()


        ##########
        # Expression

        # TPMdata[whichRNA][TPMsigDiffBinsHigh[whichRNA]].head()
        # TPMdata[whichRNA][TPMsigDiffBinsHigh[whichRNA]].ix[:,5:]

        #####
        # High

        for whichRNA in TPMdata.keys():
            sns.set(style="darkgrid")
            fig = plt.figure(figsize=(8,6))

            ax = sns.violinplot(data=TPMdata[whichRNA][TPMsigDiffBinsHigh[whichRNA]].ix[:,4:], orient="h")

            ax.set_xlabel('Genes expression within differential regions (TPM)', fontsize=10)
            ax.set_ylabel('RNA-seq data set', fontsize=10)

            # Tidy up and save
            fig.suptitle(whichRNA + " expression within positive Del-WT Xi differential regions \n" + figTitleSuffix +
            #               "\n(n=" + str(covScoresZscored[dataSet].shape[0]) + ")", fontsize=12, y=1.01)
                           "; n=" + str(np.size(np.where(TPMsigDiffBinsHigh[whichRNA]))), fontsize=12, y=1.01)
            #fig.subplots_adjust(top=1)
            #g.fig.tight_layout()
            fig.savefig(diffDir + "/violinplots.expression." + whichRNA + ".PosDel-WTXiDifferentialRegions." + figNameSuffix + ".pdf", dpi=150, bbox_inches='tight')
            fig.savefig(diffDir + "/violinplots.expression." + whichRNA + ".PosDel-WTXiDifferentialRegions." + figNameSuffix + ".png", dpi=150, bbox_inches='tight')

            plt.close()

        #####
        # In Low

        for whichRNA in TPMdata.keys():
            sns.set(style="darkgrid")
            fig = plt.figure(figsize=(8,6))

            ax = sns.violinplot(data=TPMdata[whichRNA][TPMsigDiffBinsLow[whichRNA]].ix[:,4:], orient="h")

            ax.set_xlabel('Genes expression within differential regions (TPM)', fontsize=10)
            ax.set_ylabel('RNA-seq data set', fontsize=10)

            # Tidy up and save
            fig.suptitle(whichRNA + " expression within negative Del-Wt Xi differential regions \n" + figTitleSuffix +
            #               "\n(n=" + str(covScoresZscored[dataSet].shape[0]) + ")", fontsize=12, y=1.01)
                           "; n=" + str(np.size(np.where(TPMsigDiffBinsLow[whichRNA]))), fontsize=12, y=1.01)
            #fig.subplots_adjust(top=1)
            #g.fig.tight_layout()
            fig.savefig(diffDir + "/violinplots.expression." + whichRNA + ".NegDel-WTXiDifferentialRegions." + figNameSuffix + ".pdf", dpi=150, bbox_inches='tight')
            fig.savefig(diffDir + "/violinplots.expression." + whichRNA + ".NegDel-WTXiDifferentialRegions." + figNameSuffix + ".png", dpi=150, bbox_inches='tight')

            plt.close()


# #############################################################################
# #############################################################################
# # MAIN CALL
# if __name__=="__main__":
#     main()
