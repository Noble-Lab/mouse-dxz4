#############################################################################
#
# Giancarlo Bonora
#
# 160427
# Custom contact matrix plotting:
# Various contact map visualizations
# Absolute and differential maps
# PCA
#
# 160502
# Improved and streamlined.
# Not closely related to Kate's code anymore.
#
# 160711
# Added quantile normalization prior to pairwise comparisons
#
# 160821
# Quant norm across multiple samples
# Inferno color scheme
# Higher resolutions
#
# 160823
# Changed input and output folder structure
# Implemented a config file parser
# Removed grid lines from contact maps
# Unmappable regions as NaNs: rows/cols with zero counts
# Diagonal as NaNs, as well as off-diagonals
# NaN-aware code
# log(+1) data
# Millionify total counts
# Spearman as Pearson of ranks (faster)
# All map types called as a function
#
# 160829
# Annotated heat maps with loci of interest (LOCI)
# Rotated heat maps by 45'
# Rotated heat maps with zoom
#
# 161102
# Plot contact matrix ZOOMED IN *BUT* WITHOUT ROTATION
# Discovered how to prevent blurred image after zooming:
#   http://stackoverflow.com/questions/20010882/how-to-turn-off-blurry-effect-of-imshow-in-matplotlib
# Plotting calls to different plot types now more consistent.
#
# 170202
# Differential plots
#
# 170210
# Different color map options
# virtual 4C plots
#
# 170214
# Made defaultAbsCmap and 'seismic_r' the default sequential and divergent color schemes respectively
# Rotation only plots (without annotations)
#
# 170307
# Made defaultAbsCmap, sns.color_palette("RdBu", n_colors=51)[::-1] and 'seismic'
#   the default absolute, correlation divergent and differential divergent color schemes respectively
# Added vmin=vmin, vmax=vmax to all imshow() calls
# Symmetrical correlation color maps
#
# 170623
# Correlation matrices up front to save time
# And to produce delta correlation matrices
#
# 170803
# http://pandas.pydata.org/pandas-docs/version/0.19.2/generated/pandas.DataFrame.icol.html
# DataFrame.icol(i)[source] DEPRECATED. Use .iloc[:, i] instead
#
# 170804
# Differential correlation matrices w/ wider white band
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

import numpy as np
# import scipy.stats as sps
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
from matplotlib.colors import LogNorm
import matplotlib.cm

import seaborn as sns

# from sklearn.utils.graph import graph_laplacian
# from numpy.linalg import eigh
# from numpy.linalg import inv
#
# from sklearn.decomposition import PCA
#
# from scipy.cluster.hierarchy import dendrogram, linkage
# from scipy.spatial.distance import pdist, squareform

# 160829
# For image rotation
from scipy import ndimage
# To resize color bars
from mpl_toolkits.axes_grid1 import make_axes_locatable

# 161130
from matplotlib.font_manager import FontProperties

# 170227
from scipy import signal

# 170322
from itertools import compress

USAGE = """USAGE: plotContactMapsV2.py <config file name>
    Plots a variety of heat maps given frequency matrices pointed to by config file.
"""


#############################################################################
#############################################################################
##############################################################################
#  FUNCTIONS

##############################################################################
##############################################################################
# Misc functions

##############################################################################
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
        print arr.shape
        print np.mean(arr)

        arrShape[whatData] = np.shape(arr)
        vecArr = arr.flatten()

        if tally == 0:
            dfArr = vecArr
        else:
            dfArr = pd.DataFrame(data=np.column_stack((dfArr, vecArr)))

        print dfArr.shape

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
##############################################################################
# Contact map plotting functions

defaultAbsCmap = matplotlib.cm.gist_heat_r
defaultCorrDivCmap = matplotlib.colors.ListedColormap(sns.color_palette("RdBu", n_colors=256)[::-1])
defaultDiffDivCmap = matplotlib.cm.seismic

#############################################################################
# http://stackoverflow.com/questions/7404116/defining-the-midpoint-of-a-colormap-in-matplotlib
def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower ofset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax/(vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highets point in the colormap's range.
          Defaults to 1.0 (no upper ofset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False),
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap


#############################################################################
# Adapted:
# http://stackoverflow.com/questions/7404116/defining-the-midpoint-of-a-colormap-in-matplotlib
def shiftedColorMapWithMidBand(cmap, midBandStart=0.5, midBandStop=0.5, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower ofset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax/(vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highets point in the colormap's range.
          Defaults to 1.0 (no upper ofset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    # reg_index = np.linspace(0, 1, 257)
    reg_index = np.hstack([
        np.linspace(0, 0.5, 85, endpoint=False),
        np.linspace(0.5, 0.5, 86, endpoint=False),
        np.linspace(0.5, 1, 86, endpoint=True)
    ])
    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0, midBandStart, 85, endpoint=False),
        np.linspace(midBandStart, midBandStop, 86, endpoint=False),
        np.linspace(midBandStop, 1, 86, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap


#############################################################################
# matplot lib

# 160916 - This was not necessary before, but for the LOI, I mess w/ the color bar
#          And that seems to affect subsequent plots.
#          So just borrowed code from plot_matplotlib_heatmap_with_LOI().
#//# def plot_matplotlib_heatmap(data, title, fileName, cmap=defaultAbsCmap):
#//#     # data = data[::-1]
#//#     plt.imshow(data, interpolation='nearest', cmap=cmap)
#//#     plt.grid(b=False)
#//#     plt.colorbar()
#//#     plt.title(title)
#//#     plt.savefig(fileName + ".png", dpi=300, bbox_inches='tight')
#//#     plt.clf()
def plot_matplotlib_heatmap(data, rez, title, fileName, cmap=defaultAbsCmap, vmin=None, vmax=None):
    # data = data[::-1]

    # 161102
    # Blurred figures due to the default interpolation which is set to 'bilinear'.
    # I think 'none' would be a more intuitive default.
    # http://stackoverflow.com/questions/20010882/how-to-turn-off-blurry-effect-of-imshow-in-matplotlib
    img = plt.imshow(data, cmap=cmap, interpolation='none', vmin=vmin, vmax=vmax)
    plt.grid(b=False)

    # Reset tick marks
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='major',     # both major and minor ticks are affected
        bottom='on',       # ticks along the bottom edge are on
        top='off',         # ticks along the top edge are off
        labelbottom='on',  # labels along the bottom edge are off
        color='black',
        length=5,
        width=0.5)
    plt.tick_params(
        axis='y',
        which='major',
        left='on',
        right='off',
        labelleft='on',
        direction='out',
        color='black',
        length=5,
        width=0.5)

    # Reset tick labels to be Mb
    ax = plt.gca()
    # http://stackoverflow.com/questions/11244514/modify-tick-label-text
    # See CT zhu response
    #labels = [item.get_text() for item in ax.get_yticklabels()]
    labels = ax.get_xticks().tolist()
    deciFlag = False
    for labelIdx in range(len(labels)-2):
        labelIdx += 1
        oldLabel = labels[labelIdx]
        newLabel = round(float(oldLabel),1) / (1000000.0/ rez)
        if (newLabel < 1 and newLabel > 0) or deciFlag:
            deciFlag = True
        else:
            newLabel = int(newLabel)
        labels[labelIdx] = str(newLabel)
        #print oldLabel, newLabel
    ax.set_yticklabels(labels)
    ax.set_ylabel('Mb')
    ax.set_xticklabels(labels)
    ax.set_xlabel('Mb')

    # More Matplotlib subplot weirdness!
    # http://stackoverflow.com/questions/18195758/set-matplotlib-colorbar-size-to-match-graph
    # create an axes on the right side of ax. The width of cax will be 5%
    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
    # from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    #fig.colorbar(img, cax=cax)
    plt.colorbar(img, cax=cax)

    # 160916 - Because of colorbar hack, must hack this too.
	# plt.title(title)
    # 161104 - Added to match other plot functions
    #  plt.suptitle(title)
    ax.set_title(title)

    # 161104 - Added to match other plot functions
    fig = plt.gcf()
    fig.tight_layout()
    fig.set_size_inches(8, 8)

    plt.savefig(fileName + ".png", dpi=300, bbox_inches='tight')
    plt.savefig(fileName + ".pdf", dpi=300, bbox_inches='tight')
    plt.close()


#############################################################################
# matplot lib: Entire contact matrix with LOI
def plot_matplotlib_heatmap_with_LOI(data, loi, rez, title, fileName, cmap=defaultAbsCmap, vmin=None, vmax=None):
    rez = int(rez)
    # data = data[::-1]

    # 161102
    # Blurred figures due to the default interpolation which is set to 'bilinear'.
    # I think 'none' would be a more intuitive default.
    # http://stackoverflow.com/questions/20010882/how-to-turn-off-blurry-effect-of-imshow-in-matplotlib
    img = plt.imshow(data, cmap=cmap, interpolation='none', vmin=vmin, vmax=vmax)
    # ax.set_axis_bgcolor('gray')
    plt.grid(b=False)

    # Reset tick marks
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='on',          # ticks along the top edge are on
        labelbottom='off', # labels along the bottom edge are off
        color='black',
        length=5,
        width=0.5)
    plt.tick_params(
        axis='y',
        which='major',
        left='on',
        right='on',
        labelleft='on',
        direction='out',
        color='black',
        length=5,
        width=0.5)

    # Reset tick labels to be Mb
    ax = plt.gca()
    # http://stackoverflow.com/questions/11244514/modify-tick-label-text
    # See CT zhu response
    #labels = [item.get_text() for item in ax.get_yticklabels()]
    labels = ax.get_xticks().tolist()
    deciFlag = False
    for labelIdx in range(len(labels)-2):
        labelIdx += 1
        oldLabel = labels[labelIdx]
        newLabel = round(float(oldLabel),1) / (1000000.0/ rez)
        if (newLabel < 1 and newLabel > 0) or deciFlag:
            deciFlag = True
        else:
            newLabel = int(newLabel)
        labels[labelIdx] = str(newLabel)
        print oldLabel, newLabel
    ax.set_yticklabels(labels)
    ax.set_ylabel('Mb')

    # Set up default x and y limits
    xylims = [0, data.shape[0]]

    # Plot loi lines
    tally = 0
    for aloi in loi:
        tally += 1

        # aloi = "dzx4"
        print aloi
        chrCoords = loi[aloi].split(":")
        startStop = chrCoords[1].split("-")
        startStop = map(int, startStop)
        startStop = [round(x / rez) for x in startStop]
        print chrCoords, ";  ", startStop
        plt.axhline(int(startStop[0]), xylims[0], xylims[1],
                    color='white', linewidth=0.5, linestyle='dashdot', alpha=1)
        plt.axvline(int(startStop[0]), xylims[0], xylims[1],
                    color='white', linewidth=0.5, linestyle='dashdot', alpha=1)
        plt.annotate(aloi, xy=(startStop[0], xylims[1]-1), # minus one so that within plot bounds
                     xycoords='data', xytext=(-50, -15 + (-10 * tally)), textcoords='offset points',
                     arrowprops=dict(arrowstyle="->", connectionstyle="angle,angleA=0,angleB=-90,rad=10"))

    # Color bar
    # More Matplotlib subplot weirdness!
    # http://stackoverflow.com/questions/18195758/set-matplotlib-colorbar-size-to-match-graph
    # create an axes on the right side of ax. The width of cax will be 5%
    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
    # from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    #fig.colorbar(img, cax=cax)
    plt.colorbar(img, cax=cax)

    # 160916 - Because of colorbar hack, must hack this too.
	# plt.title(title)
    # 161104 - Added to match other plot functions
    #  plt.suptitle(title)
    ax.set_title(title)

    # 161104 - Added to match other plot functions
    fig = plt.gcf()
    fig.tight_layout()
    fig.set_size_inches(8, 8)

    plt.savefig(fileName + ".LOI.png", dpi=300, bbox_inches='tight')
    plt.savefig(fileName + ".LOI.pdf", dpi=300, bbox_inches='tight')
    plt.close()


#############################################################################
# matplot lib: Rotated contact matrix
def plot_matplotlib_heatmap_r45(data, rez, title, fileName, cmap=defaultAbsCmap, vmin=None, vmax=None):
    rez = int(rez)
    # data = data[::-1]

    # Set up default x and y limits
    xylims = [0, data.shape[0]]

    # data = data[::-1]
    # from scipy import ndimage
    data_r45 = ndimage.rotate(data, 45, reshape=True, order=0, mode='constant', cval=float('nan'))
    # np.where(~np.isnan(data)); np.where(~np.isnan(data_r45))

    # Plot
    #fig, ax = plt.subplots()
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    # 161102
    # Blurred figures due to the default interpolation which is set to 'bilinear'.
    # I think 'none' would be a more intuitive default.
    # http://stackoverflow.com/questions/20010882/how-to-turn-off-blurry-effect-of-imshow-in-matplotlib
    img = ax.imshow(data_r45, cmap=cmap, interpolation='none', vmin=vmin, vmax=vmax)

    # Set up new x and y limits, which are now different:
    # *** X AND Y CO-ORDINATES HAVe BEEN RESCALED BY THE ROTATION ***
    # *** SO ALL LOCI CO-ORDS NEED TO BE RESCALED TOO! NOT IDEAL! ***
    xlims_r45 = [0, data_r45.shape[0]]
    ylims_r45 = [int(data_r45.shape[1] / 2), data_r45.shape[1]]
    ax.set_ylim(ylims_r45[0], ylims_r45[1])
    # ax.set_axis_off()
    # ax.set_axis_bgcolor('gray')
    # fig.patch.set_facecolor('gray')
    ax.axes.xaxis.set_ticklabels([])
    ax.axes.yaxis.set_ticklabels([])
    plt.grid(b=False)

    # Plot lines around matrix
    xmaxx = int((xlims_r45[1] * np.cos(45 * np.pi / 180.))
                * np.cos(45 * np.pi / 180.)) + 1
    x = np.array(range(0, xmaxx, 1))
    ax.plot(x, int(data_r45.shape[1] / 2) + x,
            color='black', linewidth=1.5, alpha=1)
    ax.plot(int(xlims_r45[1]) - x, int(ylims_r45[0]) + x,
            color='black', linewidth=1.5, alpha=1)

    # Color bar
    # More Matplotlib subplot weirdness!
    # http://stackoverflow.com/questions/18195758/set-matplotlib-colorbar-size-to-match-graph
    # create an axes on the right side of ax. The width of cax will be 5%
    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
    # from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    fig.colorbar(img, cax=cax)

    # Tidy up and save
    # ax.set_xlim(xlims_r45[0], xlims_r45[1])
    # ax.set_ylim(ylims_r45[0], ylims_r45[1])
    ax.set_title(title)
    fig.tight_layout()
    fig.set_size_inches(8, 8)
    # http://stackoverflow.com/questions/7906365/matplotlib-savefig-plots-different-from-show
    fig.savefig(fileName + ".ROT45.png", dpi=300, bbox_inches='tight',
                facecolor=fig.get_facecolor(), edgecolor='none')
    fig.savefig(fileName + ".ROT45.pdf", dpi=300, bbox_inches='tight',
                facecolor=fig.get_facecolor(), edgecolor='none')
    plt.close()


#############################################################################
# matplot lib: Rotated contact matrix with LOI
def plot_matplotlib_heatmap_with_LOI_r45(data, loi, rez, title, fileName, cmap=defaultAbsCmap, vmin=None, vmax=None):
    rez = int(rez)
    # data = data[::-1]

    # Set up default x and y limits
    xylims = [0, data.shape[0]]

    # data = data[::-1]
    # from scipy import ndimage
    data_r45 = ndimage.rotate(data, 45, reshape=True, order=0, mode='constant', cval=float('nan'))
    # np.where(~np.isnan(data)); np.where(~np.isnan(data_r45))

    # Plot
    #fig, ax = plt.subplots()
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    # 161102
    # Blurred figures due to the default interpolation which is set to 'bilinear'.
    # I think 'none' would be a more intuitive default.
    # http://stackoverflow.com/questions/20010882/how-to-turn-off-blurry-effect-of-imshow-in-matplotlib
    img = ax.imshow(data_r45, cmap=cmap, interpolation='none', vmin=vmin, vmax=vmax)

    # Set up new x and y limits, which are now different:
    # *** X AND Y CO-ORDINATES HAVe BEEN RESCALED BY THE ROTATION ***
    # *** SO ALL LOCI CO-ORDS NEED TO BE RESCALED TOO! NOT IDEAL! ***
    xlims_r45 = [0, data_r45.shape[0]]
    ylims_r45 = [int(data_r45.shape[1] / 2), data_r45.shape[1]]
    ax.set_ylim(ylims_r45[0], ylims_r45[1])
    # ax.set_axis_off()
    # ax.set_axis_bgcolor('gray')
    # fig.patch.set_facecolor('gray')
    ax.axes.xaxis.set_ticklabels([])
    ax.axes.yaxis.set_ticklabels([])
    plt.grid(b=False)

    # Plot loi lines
    tally = 0
    for aloi in loi:
        tally += 1

        # aloi = "dzx4"
        # print aloi
        chrCoords = loi[aloi].split(":")
        startStop = chrCoords[1].split("-")
        startStop = map(int, startStop)
        startStop = [round(x / rez) for x in startStop]
        # *** X AND Y CO-ORDINATES HAVe BEEN RESCALED BY THE ROTATION ***
        # *** SO ALL LOCI CO-ORDS NEED TO BE RESCALED TOO! NOT IDEAL! ***
        # startStop = np.array(startStop) * (xlims_r45[1] / float(xylims[1]))
        startStop = [round(x * (xlims_r45[1] / float(xylims[1]))) for x in startStop]

        # print chrCoords, ";  ", startStop
        # Matplotlib subplot weirdness!
        # http://stackoverflow.com/questions/11864975/a-bug-in-axvline-of-matplotlib
        # ax.axvline(x=int(startStop[0]),ymin=ylims_r45[0],ymax=ylims_r45[1],
        #           color='white',linewidth=0.5,linestyle='dashdot',alpha=1)
        # ax.axvline(x=int(startStop[0]), ymin=0, ymax=1,
        #             color='white', linewidth=1, linestyle='dashdot', alpha=1)
        # But didn't use axvline afterall:
        xmaxx = int(((xlims_r45[1] - startStop[0]) * np.cos(45 * np.pi / 180.))
                    * np.cos(45 * np.pi / 180.)) + 1
        x = np.array(range(0, xmaxx, 1))
        ax.plot(startStop[0] + x, int(data_r45.shape[1] / 2) + x,
                color='white', linewidth=0.5, alpha=1)
        # ax.plot(startStop[0] + x, int(data_r45.shape[1] / 2) + x,
        #         color='black', linewidth=0.25, alpha=1)
        xmaxx = int(startStop[0] - (((startStop[0] * np.cos(45 * np.pi / 180.))
                                     * np.cos(45 * np.pi / 180.)))) + 1
        x = np.array(range(0, xmaxx, 1))
        ax.plot(startStop[0] - x, int(ylims_r45[0]) + x,
                color='white', linewidth=0.5, alpha=1)
        # ax.plot(startStop[0] - x, int(ylims_r45[0]) + x,
        #         color='black', linewidth=0.25, alpha=1)
        plt.annotate(aloi, xy=(startStop[0], ylims_r45[0]),
                     xycoords='data', xytext=(-50, -15 + (-10 * tally)), textcoords='offset points',
                     arrowprops=dict(arrowstyle="->", connectionstyle="angle,angleA=0,angleB=-90,rad=10"))

    # Plot lines around matrix
    xmaxx = int((xlims_r45[1] * np.cos(45 * np.pi / 180.))
                * np.cos(45 * np.pi / 180.)) + 1
    x = np.array(range(0, xmaxx, 1))
    ax.plot(x, int(data_r45.shape[1] / 2) + x,
            color='black', linewidth=1.5, alpha=1)
    ax.plot(int(xlims_r45[1]) - x, int(ylims_r45[0]) + x,
            color='black', linewidth=1.5, alpha=1)

    # Color bar
    # More Matplotlib subplot weirdness!
    # http://stackoverflow.com/questions/18195758/set-matplotlib-colorbar-size-to-match-graph
    # create an axes on the right side of ax. The width of cax will be 5%
    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
    # from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    fig.colorbar(img, cax=cax)

    # Tidy up and save
    # ax.set_xlim(xlims_r45[0], xlims_r45[1])
    # ax.set_ylim(ylims_r45[0], ylims_r45[1])
    ax.set_title(title)
    fig.tight_layout()
    fig.set_size_inches(8, 8)
    # http://stackoverflow.com/questions/7906365/matplotlib-savefig-plots-different-from-show
    fig.savefig(fileName + ".LOI.ROT45.png", dpi=300, bbox_inches='tight',
                facecolor=fig.get_facecolor(), edgecolor='none')
    fig.savefig(fileName + ".LOI.ROT45.pdf", dpi=300, bbox_inches='tight',
                facecolor=fig.get_facecolor(), edgecolor='none')
    plt.close()


#############################################################################
# matplot lib: Rotated contact matrix zoomed in on LOI
def plot_matplotlib_heatmap_with_LOI_r45_zoom(data, loi, rez, zoomWin, title, fileName, cmap=defaultAbsCmap, vmin=None, vmax=None):
    rez = int(rez)
    zoomWin = int(zoomWin)
    if zoomWin >= rez * 3:

        # Set up default x and y limits
        xylims = [0, data.shape[0]]

        # data = data[::-1]
        # from scipy import ndimage
        data_r45 = ndimage.rotate(data, 45, reshape=True, order=0, mode='constant', cval=float('nan'))
        # np.where(~np.isnan(data)); np.where(~np.isnan(data_r45))

        # Plot loi regions
        for aloi in loi:
            #fig, ax = plt.subplots()
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)

            # Plot
            # 161102
            # Blurred figures due to the default interpolation which is set to 'bilinear'.
            # I think 'none' would be a more intuitive default.
            # http://stackoverflow.com/questions/20010882/how-to-turn-off-blurry-effect-of-imshow-in-matplotlib
            img = ax.imshow(data_r45, cmap=cmap, interpolation='none', vmin=vmin, vmax=vmax)

            # Set up new x and y limits, which are now different:
            # *** X AND Y CO-ORDINATES HAVe BEEN RESCALED BY THE ROTATION ***
            # *** SO ALL LOCI CO-ORDS NEED TO BE RESCALED TOO! NOT IDEAL! ***
            xlims_r45 = [0, data_r45.shape[0]]
            ylims_r45 = [int(data_r45.shape[1] / 2), data_r45.shape[1]]
            ax.set_ylim(ylims_r45[0], ylims_r45[1])
            # ax.set_axis_off()
            # ax.set_axis_bgcolor('gray')
            # fig.patch.set_facecolor('gray')
            ax.axes.xaxis.set_ticklabels([])
            ax.axes.yaxis.set_ticklabels([])
            plt.grid(b=False)

            # aloi = "dzx4"
            # print aloi
            chrCoords = loi[aloi].split(":")
            startStop = chrCoords[1].split("-")
            startStop = map(int, startStop)
            startStop = [round(x / rez) for x in startStop]
            # *** X AND Y CO-ORDINATES HAVe BEEN RESCALED BY THE ROTATION ***
            # *** SO ALL LOCI CO-ORDS NEED TO BE RESCALED TOO! NOT IDEAL! ***
            # startStop = np.array(startStop) * (xlims_r45[1] / float(xylims[1]))
            startStop = [round(x * (xlims_r45[1] / float(xylims[1]))) for x in startStop]

            # print chrCoords, ";  ", startStop
            # Matplotlib subplot weirdness!
            # http://stackoverflow.com/questions/11864975/a-bug-in-axvline-of-matplotlib
            # ax.axvline(x=int(startStop[0]),ymin=ylims_r45[0],ymax=ylims_r45[1],
            #           color='white',linewidth=0.5,linestyle='dashdot',alpha=1)
            # ax.axvline(x=int(startStop[0]), ymin=0, ymax=1,
            #             color='white', linewidth=1, linestyle='dashdot', alpha=1)
            # But didn't use axvline afterall:
            xmaxx = int(((xlims_r45[1] - startStop[0]) * np.cos(45 * np.pi / 180.))
                        * np.cos(45 * np.pi / 180.)) + 1
            x = np.array(range(0, xmaxx, 1))
            ax.plot(startStop[0] + x, int(data_r45.shape[1] / 2) + x,
                    color='white', linewidth=1, alpha=1)
            ax.plot(startStop[0] + x, int(data_r45.shape[1] / 2) + x,
                    color='black', linewidth=0.25, alpha=1)
            xmaxx = int(startStop[0] - (((startStop[0] * np.cos(45 * np.pi / 180.))
                                         * np.cos(45 * np.pi / 180.)))) + 1
            x = np.array(range(0, xmaxx, 1))
            ax.plot(startStop[0] - x, int(ylims_r45[0]) + x,
                    color='white', linewidth=1, alpha=1)
            ax.plot(startStop[0] - x, int(ylims_r45[0]) + x,
                    color='black', linewidth=0.25, alpha=1)

            plt.annotate(aloi, xy=(startStop[0], ylims_r45[0]),
                         xycoords='data', xytext=(-50, -15), textcoords='offset points',
                         arrowprops=dict(arrowstyle="->", connectionstyle="angle,angleA=0,angleB=-90,rad=10"))

            # Plot lines around matrix
            xmaxx = int((xlims_r45[1] * np.cos(45 * np.pi / 180.))
                        * np.cos(45 * np.pi / 180.)) + 1
            x = np.array(range(0, xmaxx, 1))
            ax.plot(x, int(data_r45.shape[1] / 2) + x,
                    color='black', linewidth=1.5, alpha=1)
            ax.plot(int(xlims_r45[1]) - x, int(ylims_r45[0]) + x,
                    color='black', linewidth=1.5, alpha=1)

            # Zoom in
            # In terms of windows, and scaled for rot45
            plusMinusDist = ((zoomWin / 2) / rez) * (xlims_r45[1] / float(xylims[1]))
            ax.set_xlim(max(0, startStop[0] - plusMinusDist),
                        min(xlims_r45[1], startStop[0] + plusMinusDist))
            ax.set_ylim(int(data_r45.shape[1] / 2), int(data_r45.shape[1] / 2) + plusMinusDist)

            # Color bar
            # More Matplotlib subplot weirdness!
            # http://stackoverflow.com/questions/18195758/set-matplotlib-colorbar-size-to-match-graph
            # create an axes on the right side of ax. The width of cax will be 5%
            # of ax and the padding between cax and ax will be fixed at 0.05 inch.
            # from mpl_toolkits.axes_grid1 import make_axes_locatable
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            fig.colorbar(img, cax=cax)

            # Tidy up and save
            aloiTitle = title + " " + aloi + " " + sizeof_fmt(zoomWin)
            ax.set_title(aloiTitle)
            fig.tight_layout()
            fig.set_size_inches(8, 8)
            # http://stackoverflow.com/questions/7906365/matplotlib-savefig-plots-different-from-show
            fig.savefig(fileName + ".LOI.ROT45." +
                        sizeof_fmt(zoomWin) + "." + aloi + ".png",
                        dpi=300, bbox_inches='tight',
                        facecolor=fig.get_facecolor(), edgecolor='none')
            fig.savefig(fileName + ".LOI.ROT45." +
                        sizeof_fmt(zoomWin) + "." + aloi + ".pdf",
                        dpi=300, bbox_inches='tight',
                        facecolor=fig.get_facecolor(), edgecolor='none')
            plt.close()


#############################################################################
# 161102
# matplot lib: Contact matrix zoomed in on LOI

# THIS NEEDED TO RESET AXES AFTER ZOOMING
# http://stackoverflow.com/questions/11086724/matplotlib-linked-x-axes-with-autoscaled-y-axes-on-zoom
def on_xlim_changed(ax):
    xlim = ax.get_xlim()
    for a in ax.figure.axes:
        # shortcuts: last avoids n**2 behavior when each axis fires event
        if a is ax or len(a.lines) == 0 or getattr(a, 'xlim', None) == xlim:
            continue

        ylim = np.inf, -np.inf
        for l in a.lines:
            x, y = l.get_data()
            # faster, but assumes that x is sorted
            start, stop = np.searchsorted(x, xlim)
            yc = y[max(start-1,0):(stop+1)]
            ylim = min(ylim[0], np.nanmin(yc)), max(ylim[1], np.nanmax(yc))

        # TODO: update limits from Patches, Texts, Collections, ...

        # x axis: emit=False avoids infinite loop
        a.set_xlim(xlim, emit=False)

        # y axis: set dataLim, make sure that autoscale in 'y' is on
        corners = (xlim[0], ylim[0]), (xlim[1], ylim[1])
        a.dataLim.update_from_data_xy(corners, ignore=True, updatex=False)
        a.autoscale(enable=True, axis='y')
        # cache xlim to mark 'a' as treated
        a.xlim = xlim


def plot_matplotlib_heatmap_with_LOI_zoom(data, loi, rez, zoomWin, title, fileName, cmap=defaultAbsCmap, vmin=None, vmax=None):
    rez = int(rez)
    zoomWin = int(zoomWin)

    # Set up default x and y limits
    xylims = [0, data.shape[0]]

    if zoomWin >= rez * 3:

        # Plot loi regions
        for aloi in loi:
            # aloi = 'dxz4'
            # print aloi

            aloiTitle = title + " " + aloi + " " + sizeof_fmt(zoomWin)
            aloiFileName = fileName + "." + aloi + "." + sizeof_fmt(zoomWin)

            # *** THIS METHOD ZOOMS INTO A PLOT BASED ON ENTIRE MATRIX, ***
            # *** BUT AXES AND COLOR BAR STILL PERTAIN TO ENTIRE MATRIX ***
            # *** SO IS IT BETTER TO JUST PLOT THE RELEVANT BINS??????? ***
            # Plot
            fig, ax = plt.subplots()
            # 161102
            # Blurred figures due to the default interpolation which is set to 'bilinear'.
            # I think 'none' would be a more intuitive default.
            # http://stackoverflow.com/questions/20010882/how-to-turn-off-blurry-effect-of-imshow-in-matplotlib
            img = ax.imshow(data, cmap=cmap, interpolation='none', vmin=vmin, vmax=vmax)
            # ax.set_axis_bgcolor('gray')
            plt.grid(b=False)

            # Reset tick marks
            plt.tick_params(
                axis='x',          # changes apply to the x-axis
                which='both',      # both major and minor ticks are affected
                bottom='on',      # ticks along the bottom edge are on
                top='off',          # ticks along the top edge are off
                labeltop='off', # labels along the bottom edge are off
                color='black',
                length=5,
                width=0.5)
            plt.tick_params(
                axis='y',
                which='major',
                left='on',
                right='off',
                labelleft='on',
                direction='out',
                color='black',
                length=5,
                width=0.5)

            chrCoords = loi[aloi].split(":")
            startStop = chrCoords[1].split("-")
            startStop = map(int, startStop)
            startStop = [round(x / rez) for x in startStop]

            # Zoom in
            plusMinusDist = (zoomWin / 2) / rez
            ax.set_xlim(max(xylims[0], startStop[0] - plusMinusDist),
                        min(xylims[1], int(startStop[0] + plusMinusDist)))
            ax.set_ylim(min(xylims[1], int(startStop[0] + plusMinusDist)),
                        max(xylims[0], startStop[0] - plusMinusDist))

            # http://stackoverflow.com/questions/11086724/matplotlib-linked-x-axes-with-autoscaled-y-axes-on-zoom
            for ax in fig.axes:
                ax.callbacks.connect('xlim_changed', on_xlim_changed)

            # Reset tick labels to be Mb
            # ax = plt.gca()
            # http://stackoverflow.com/questions/11244514/modify-tick-label-text
            # See CT zhu response
            #labels = [item.get_text() for item in ax.get_yticklabels()]
            labels = ax.get_xticks().tolist()
            deciFlag = False
            for labelIdx in range(len(labels)-2):
                labelIdx += 1
                oldLabel = labels[labelIdx]
                newLabel = round(float(oldLabel),1) / (1000000.0/ rez)
                if (newLabel < 1 and newLabel > 0) or deciFlag:
                    deciFlag = True
                else:
                    newLabel = int(newLabel)
                labels[labelIdx] = str(newLabel)
                print oldLabel, newLabel
            ax.set_yticklabels(labels)
            ax.set_ylabel('Mb')
            ax.set_xticklabels(labels)
            ax.set_xlabel('Mb')

            # Color bar
            # More Matplotlib subplot weirdness!
            # http://stackoverflow.com/questions/18195758/set-matplotlib-colorbar-size-to-match-graph
            # create an axes on the right side of ax. The width of cax will be 5%
            # of ax and the padding between cax and ax will be fixed at 0.05 inch.
            # from mpl_toolkits.axes_grid1 import make_axes_locatable
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            fig.colorbar(img, cax=cax)
            # plt.colorbar(img, cax=cax)

            # Tidy up and save
            ax.set_title(aloiTitle)
            fig.tight_layout()
            fig.set_size_inches(8, 8)
            # # http://stackoverflow.com/questions/7906365/matplotlib-savefig-plots-different-from-show
            fig.savefig(aloiFileName + ".png", dpi=300, bbox_inches='tight',
                         facecolor=fig.get_facecolor(), edgecolor='none')
            fig.savefig(aloiFileName + ".pdf", dpi=300, bbox_inches='tight',
                         facecolor=fig.get_facecolor(), edgecolor='none')
            plt.close()


#############################################################################
# Different color maps

def plot_contact_map(data, rez, title, outFilenameBase, cmap=defaultAbsCmap):
    # outFilename = outFilenameBase + ".jet"
    # plot_matplotlib_heatmap(data, title, outFilename)
    # outFilename = outFilenameBase + ".CMRmap"
    # plot_matplotlib_heatmap(data, title, outFilename,
    #                         cmap="CMRmap")
    outFilename = outFilenameBase
    plot_matplotlib_heatmap(data, rez, title, outFilename, cmap=cmap)


def plot_divergent_maps(data, rez, title, outFilenameBase, cmap=defaultDiffDivCmap):
    outFilename = outFilenameBase
    plot_matplotlib_heatmap(data, rez, title, outFilename, cmap=cmap)


##############################################################################
# Print 'em all!
# Do correlation in function
def printAllMapsAndPerformCorrelation(data, chrom, outputLabel, outDir, theRez,
                 cmapSeq=defaultAbsCmap, cmapDiv=defaultCorrDivCmap):
    print "\nprintAllMaps" + ", " + chrom + ", " + outputLabel + ", " + outDir + ", " + theRez
    theRez = int(theRez)

    dataDF = pd.DataFrame(data)

    ################################################
    # Raw data

    print "Raw"
    title = chrom + " " + outputLabel + " linear"
    outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".linear"
    plot_contact_map(data, theRez, title, outFilenameBase, cmap=cmapSeq)


    ################################################
    # Clipped data

    print "Clipped"
    title = chrom + " " + outputLabel + " linear Q5-Q95"
    outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".linear.Q5-Q95"
    plot_contact_map(np.clip(data,
                             a_min=np.nanpercentile(data, 5),
                             a_max=np.nanpercentile(data, 95)),
                     theRez, title, outFilenameBase, cmap=cmapSeq)

    ################################################
    # Clipped variations

    if theRez < 500000:
        print "More clipped plots"
        theMin = 0.05
        for theMax in frange(0.96, 1, 0.01):
            rangeStr = "Q" + str(int(round(theMin * 100))) + "-Q" + str(int(round(theMax * 100)))
            print rangeStr
            title = chrom + " " + outputLabel + " linear " + rangeStr
            outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".linear." + rangeStr
            plot_contact_map(np.clip(data,
                                     a_min=np.nanpercentile(data, round(theMin * 100)),
                                     a_max=np.nanpercentile(data, round(theMax * 100))),
                             theRez, title, outFilenameBase, cmap=cmapSeq)


    ################################################
    # Logged data

    print "Logged"
    title = chrom + " " + outputLabel + " logged"
    outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".logged"
    # data = [[100 for x in range(5)] for x in range(5)] # For testing
    # for i in range(min(len(data[0]),len(data))):
    #     data[i][i] = np.nan
    #     if i > 0:
    #         data[i][i-1] = np.nan
    #         data[i-1][i] = np.nan
    # # np.log(data+1) # Python is not as clever as R in this respect
    # [[z+1 for z in y] for y in data]
    # np.log([[z+1 for z in y] for y in data])
    plot_contact_map(np.log([[z + 1 for z in y] for y in data]), theRez, title, outFilenameBase, cmap=cmapSeq)


    ################################################
    # Ranked data (quantiles)

    if theRez >= 500000:
        print "Ranked"
        # dataRanks = sps.rankdata(data).reshape(len(data), len(data))
        # http://stackoverflow.com/questions/12519866/rank-array-in-python-while-ignoring-missing-values
        dataRanks = mstats.rankdata(np.ma.masked_invalid(data)).reshape(len(data), len(data))
        dataRanks[dataRanks == 0] = np.nan
        dataRanks -= 1
        title = chrom + " " + outputLabel + " ranks"
        outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".ranks"
        plot_contact_map(dataRanks, theRez, title, outFilenameBase, cmap=cmapSeq)


    ################################################
    # Inverse hyperbolicc sine (arcsinh)

    print "IHS"
    title = chrom + " " + outputLabel + " IHS"
    outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".IHS"
    plot_contact_map(np.arcsinh(data), theRez, title, outFilenameBase, cmap=cmapSeq)


    ################################################
    # observed/expected

    # Expected counts
    expectedDiagCounts = []
    for i in range(len(data)):
        thisDiag = []
        for j in np.array(range(len(data) - i)) + i:
            thisDiag.append(data[j - i, j])
        expectedDiagCounts.append(np.nanmean(thisDiag))

    # OE
    Edata = np.empty_like(data)
    OEdata = np.empty_like(data)
    for i in range(len(data)):
        for j in range(len(data)):
            Edata[i, j] = expectedDiagCounts[abs(i - j)]
            # print str(data[i, j]) + " " + str(Edata[i, j]) + " " + str(expectedDiagCounts[abs(i-j)])
            OEdata[i, j] = np.log2(float(data[i, j]) / Edata[i, j])
            # print str(np.log2(float(data[i, j]+0.01)/float(expectedDiagCounts[abs(i-j)]+0.01)))
            # print str(OEdata[i, j])

    # title = chrom + outputLabel + " expected clipped"
    # outFilenameBase =  outDir + "/" + outputLabel + "." + chrom + ".expected.clipped.Q5-Q95"
    # plot_contact_map(np.clip(Edata, a_min=np.percentile(data, 5), a_max=np.percentile(data, 95)), theRez, title, outFilenameBase)
    #
    # title = chrom + outputLabel + " expected logged"
    # outFilenameBase =  outDir + "/" + outputLabel + "." + chrom + ".expected.logged"
    # plot_contact_map(np.log(Edata), theRez, title, outFilenameBase)

    title = chrom + outputLabel + " OE"
    outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".OE"
    plot_divergent_maps(OEdata, theRez, title, outFilenameBase, cmap=cmapDiv)


    ################################################
    # Correlation matrices

    print "Pearson correlation"
    # dataDF = pd.DataFrame(data)
    # mask = datadf.isnull()

    #  corData =  np.corrcoef(data)
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
    outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".corrPearson"
    title = chrom + outputLabel + " Pearson"
    # plot_divergent_maps(corData, title, outFilenameBase)
    # Symmetrical color map using vmin and vmax
    vmax = np.nanmax(corData)
    vmin = np.nanmin(corData)
    vabs = np.nanmax([abs(vmax), abs(vmin)])
    plot_matplotlib_heatmap(corData, theRez, title, outFilenameBase, cmap=cmapDiv, vmax=vabs, vmin=-vabs)

    # if theRez >= 500000:
    #     print "Spearman correlation"
    #     dataRanksDF = pd.DataFrame(dataRanks)
    #     # spearcorData =  pd.DataFrame.corr(df, method="spearman")
    #     # spearcorData = pd.DataFrame(dataRanks).corr
    #     spearcorData = pd.DataFrame.corr(dataRanksDF)  # Pearson of ranks == Spearman
    #     # Rediagonal & off diags
    #     # data = [[0 for x in range(5)] for x in range(5)] # For testing
    #     for i in range(min(len(spearcorData[0]), len(spearcorData))):
    #         spearcorData[i][i] = np.nan
    #         if i > 0:
    #             spearcorData[i][i - 1] = np.nan
    #             spearcorData[i - 1][i] = np.nan
    #     outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".corrSpearman"
    #     title = chrom + outputLabel + " Spearman"
    #     plot_divergent_maps(spearcorData, theRez, title, outFilenameBase, cmap=cmapDiv)


    ################################################
    # Clipped correlations

    if theRez < 500000:
        print "More Pearson correlation plots"
        for theMin in frange(-0.2, -0.1, 0.05):
            for theMax in frange(0.3, 0.5, 0.05):
                rangeStr = str(int(round(theMin * 100))) + "to" + str(int(round(theMax * 100)))
                print rangeStr

                title = chrom + " " + outputLabel + " Pearson " + rangeStr
                outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".corrPearson." + rangeStr
                # plot_divergent_maps(np.clip(corData, a_min=theMin, a_max=theMax), theRez, title, outFilenameBase, cmap=cmapDiv)
                # Symmetrical color map using vmin and vmax
                vmax = np.nanmax(corData)
                vmin = np.nanmin(corData)
                vabs = np.nanmax([abs(vmax), abs(vmin)])
                plot_matplotlib_heatmap(corData, theRez, title, outFilenameBase, cmap=cmapDiv, vmax=theMin, vmin=-theMax)

                # title = chrom + " " + outputLabel + " Spearman " + rangeStr
                # outFilenameBase =  outDir + "/" + outputLabel + "." + chrom + ".corrSpearman." + rangeStr
                # plot_divergent_maps(np.clip(spearcorData, a_min=theMin, a_max=theMax), theRez, title, outFilenameBase, cmap=cmapDiv)

    # ################################################
    # # Diffusion kernel
    #
    # #  calculation the eigendecomposition of the laplacian
    # laplacian = -graph_laplacian(data, normed=True, return_diag=False)
    # eigvals, V = eigh(laplacian)
    # inverse_V = inv(V)
    # # do the graph diffusion
    # betas = [2,5,10] # controls the amount of diffusion
    # for beta in betas:
    #     # beta = 2
    #     #logging.debug("running beta %.1f" % beta)
    #     Wb = np.diag( np.exp( np.multiply( eigvals, beta ) ) )
    #     kernel = np.dot( np.dot( V, Wb), inverse_V )
    #     #np.savetxt("kernel_"+str(beta)+".txt",kernel,delimiter='\t')
    #     ax = sns.heatmap(kernel,
    #                      linewidths=0, square=True, xticklabels=False, yticklabels=False,
    #                      cmap=cmap)
    #

    # outFilenameBase =  outDir + "/" + outputLabel + "." + chrom
    # title = chrom + " " + outputLabel
    # plot_dkheatmap(data, title, outFilenameBase)

    # outFilenameBase =  outDir + "/" + outputLabel + "." + chrom + ".OE"
    # title = chrom + " " + outputLabel + " OE"
    # plot_dkheatmap(OEdata, title, outFilenameBase)


##############################################################################
# Print 'em all!
def printAllMaps(data, corData, chrom, outputLabel, outDir, theRez,
                 cmapSeq=defaultAbsCmap, cmapDiv=defaultCorrDivCmap):
    print "\nprintAllMaps" + ", " + chrom + ", " + outputLabel + ", " + outDir + ", " + theRez
    theRez = int(theRez)

    dataDF = pd.DataFrame(data)

    ################################################
    # Raw data

    print "Raw"
    title = chrom + " " + outputLabel + " linear"
    outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".linear"
    plot_contact_map(data, theRez, title, outFilenameBase, cmap=cmapSeq)


    ################################################
    # Clipped data

    print "Clipped"
    title = chrom + " " + outputLabel + " linear Q5-Q95"
    outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".linear.Q5-Q95"
    plot_contact_map(np.clip(data,
                             a_min=np.nanpercentile(data, 5),
                             a_max=np.nanpercentile(data, 95)),
                     theRez, title, outFilenameBase, cmap=cmapSeq)

    ################################################
    # Clipped variations

    if theRez < 500000:
        print "More clipped plots"
        theMin = 0.05
        for theMax in frange(0.96, 1, 0.01):
            rangeStr = "Q" + str(int(round(theMin * 100))) + "-Q" + str(int(round(theMax * 100)))
            print rangeStr
            title = chrom + " " + outputLabel + " linear " + rangeStr
            outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".linear." + rangeStr
            plot_contact_map(np.clip(data,
                                     a_min=np.nanpercentile(data, round(theMin * 100)),
                                     a_max=np.nanpercentile(data, round(theMax * 100))),
                             theRez, title, outFilenameBase, cmap=cmapSeq)


    ################################################
    # Logged data

    print "Logged"
    title = chrom + " " + outputLabel + " logged"
    outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".logged"
    # data = [[100 for x in range(5)] for x in range(5)] # For testing
    # for i in range(min(len(data[0]),len(data))):
    #     data[i][i] = np.nan
    #     if i > 0:
    #         data[i][i-1] = np.nan
    #         data[i-1][i] = np.nan
    # # np.log(data+1) # Python is not as clever as R in this respect
    # [[z+1 for z in y] for y in data]
    # np.log([[z+1 for z in y] for y in data])
    plot_contact_map(np.log([[z + 1 for z in y] for y in data]), theRez, title, outFilenameBase, cmap=cmapSeq)


    ################################################
    # Ranked data (quantiles)

    if theRez >= 500000:
        print "Ranked"
        # dataRanks = sps.rankdata(data).reshape(len(data), len(data))
        # http://stackoverflow.com/questions/12519866/rank-array-in-python-while-ignoring-missing-values
        dataRanks = mstats.rankdata(np.ma.masked_invalid(data)).reshape(len(data), len(data))
        dataRanks[dataRanks == 0] = np.nan
        dataRanks -= 1
        title = chrom + " " + outputLabel + " ranks"
        outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".ranks"
        plot_contact_map(dataRanks, theRez, title, outFilenameBase, cmap=cmapSeq)


    ################################################
    # Inverse hyperbolicc sine (arcsinh)

    print "IHS"
    title = chrom + " " + outputLabel + " IHS"
    outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".IHS"
    plot_contact_map(np.arcsinh(data), theRez, title, outFilenameBase, cmap=cmapSeq)


    ################################################
    # observed/expected

    # Expected counts
    expectedDiagCounts = []
    for i in range(len(data)):
        thisDiag = []
        for j in np.array(range(len(data) - i)) + i:
            thisDiag.append(data[j - i, j])
        expectedDiagCounts.append(np.nanmean(thisDiag))

    # OE
    Edata = np.empty_like(data)
    OEdata = np.empty_like(data)
    for i in range(len(data)):
        for j in range(len(data)):
            Edata[i, j] = expectedDiagCounts[abs(i - j)]
            # print str(data[i, j]) + " " + str(Edata[i, j]) + " " + str(expectedDiagCounts[abs(i-j)])
            OEdata[i, j] = np.log2(float(data[i, j]) / Edata[i, j])
            # print str(np.log2(float(data[i, j]+0.01)/float(expectedDiagCounts[abs(i-j)]+0.01)))
            # print str(OEdata[i, j])

    # title = chrom + outputLabel + " expected clipped"
    # outFilenameBase =  outDir + "/" + outputLabel + "." + chrom + ".expected.clipped.Q5-Q95"
    # plot_contact_map(np.clip(Edata, a_min=np.percentile(data, 5), a_max=np.percentile(data, 95)), theRez, title, outFilenameBase)
    #
    # title = chrom + outputLabel + " expected logged"
    # outFilenameBase =  outDir + "/" + outputLabel + "." + chrom + ".expected.logged"
    # plot_contact_map(np.log(Edata), theRez, title, outFilenameBase)

    title = chrom + outputLabel + " OE"
    outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".OE"
    plot_divergent_maps(OEdata, theRez, title, outFilenameBase, cmap=cmapDiv)


    ################################################
    # Correlation matrices

    print "Pearson correlation"
    # dataDF = pd.DataFrame(data)
    # mask = datadf.isnull()

    # 170623
    # #  corData =  np.corrcoef(data)
    # # corData =  np.ma.corrcoef(data)
    # # corData = pd.DataFrame(data).corr
    # corData = pd.DataFrame.corr(dataDF)
    # # Rediagonal & off diags
    # # data = [[0 for x in range(5)] for x in range(5)] # For testing
    # for i in range(min(len(corData[0]), len(corData))):
    #     corData[i][i] = np.nan
    #     if i > 0:
    #         corData[i][i - 1] = np.nan
    #         corData[i - 1][i] = np.nan
    outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".corrPearson"
    title = chrom + outputLabel + " Pearson"
    # plot_divergent_maps(corData, title, outFilenameBase)
    # Symmetrical color map using vmin and vmax
    vmax = np.nanmax(corData)
    vmin = np.nanmin(corData)
    vabs = np.nanmax([abs(vmax), abs(vmin)])
    plot_matplotlib_heatmap(corData, theRez, title, outFilenameBase, cmap=cmapDiv, vmax=vabs, vmin=-vabs)

    # if theRez >= 500000:
    #     print "Spearman correlation"
    #     dataRanksDF = pd.DataFrame(dataRanks)
    #     # spearcorData =  pd.DataFrame.corr(df, method="spearman")
    #     # spearcorData = pd.DataFrame(dataRanks).corr
    #     spearcorData = pd.DataFrame.corr(dataRanksDF)  # Pearson of ranks == Spearman
    #     # Rediagonal & off diags
    #     # data = [[0 for x in range(5)] for x in range(5)] # For testing
    #     for i in range(min(len(spearcorData[0]), len(spearcorData))):
    #         spearcorData[i][i] = np.nan
    #         if i > 0:
    #             spearcorData[i][i - 1] = np.nan
    #             spearcorData[i - 1][i] = np.nan
    #     outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".corrSpearman"
    #     title = chrom + outputLabel + " Spearman"
    #     plot_divergent_maps(spearcorData, theRez, title, outFilenameBase, cmap=cmapDiv)


    ################################################
    # Clipped correlations

    if theRez < 500000:
        print "More Pearson correlation plots"
        for theMin in frange(-0.2, -0.1, 0.05):
            for theMax in frange(0.3, 0.5, 0.05):
                rangeStr = str(int(round(theMin * 100))) + "to" + str(int(round(theMax * 100)))
                print rangeStr

                title = chrom + " " + outputLabel + " Pearson " + rangeStr
                outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".corrPearson." + rangeStr
                # plot_divergent_maps(np.clip(corData, a_min=theMin, a_max=theMax), theRez, title, outFilenameBase, cmap=cmapDiv)
                # Symmetrical color map using vmin and vmax
                vmax = np.nanmax(corData)
                vmin = np.nanmin(corData)
                vabs = np.nanmax([abs(vmax), abs(vmin)])
                plot_matplotlib_heatmap(corData, theRez, title, outFilenameBase, cmap=cmapDiv, vmax=theMin, vmin=-theMax)

                # title = chrom + " " + outputLabel + " Spearman " + rangeStr
                # outFilenameBase =  outDir + "/" + outputLabel + "." + chrom + ".corrSpearman." + rangeStr
                # plot_divergent_maps(np.clip(spearcorData, a_min=theMin, a_max=theMax), theRez, title, outFilenameBase, cmap=cmapDiv)

    # ################################################
    # # Diffusion kernel
    #
    # #  calculation the eigendecomposition of the laplacian
    # laplacian = -graph_laplacian(data, normed=True, return_diag=False)
    # eigvals, V = eigh(laplacian)
    # inverse_V = inv(V)
    # # do the graph diffusion
    # betas = [2,5,10] # controls the amount of diffusion
    # for beta in betas:
    #     # beta = 2
    #     #logging.debug("running beta %.1f" % beta)
    #     Wb = np.diag( np.exp( np.multiply( eigvals, beta ) ) )
    #     kernel = np.dot( np.dot( V, Wb), inverse_V )
    #     #np.savetxt("kernel_"+str(beta)+".txt",kernel,delimiter='\t')
    #     ax = sns.heatmap(kernel,
    #                      linewidths=0, square=True, xticklabels=False, yticklabels=False,
    #                      cmap=cmap)
    #

    # outFilenameBase =  outDir + "/" + outputLabel + "." + chrom
    # title = chrom + " " + outputLabel
    # plot_dkheatmap(data, title, outFilenameBase)

    # outFilenameBase =  outDir + "/" + outputLabel + "." + chrom + ".OE"
    # title = chrom + " " + outputLabel + " OE"
    # plot_dkheatmap(OEdata, title, outFilenameBase)


##############################################################################
# Print some of differenr color map options
def printSomeColorMapOptions(data, chrom, outputLabel, outDir, theRez):
    print "\nprintSomeMaps" + ", " + chrom + ", " + outputLabel + ", " + outDir \
          + ", " + str(theRez)
    theRez = int(theRez)

    dataDF = pd.DataFrame(data)

    cmapOptions = ['viridis', 'inferno', 'plasma', 'magma', 'afmhot', 'copper', 'bone', 'gist_heat', 'Reds', 'Blues']

    ################################################
    # Clipped data

    print "Clipped"
    for cmapOption in cmapOptions:
        for cmapSuffix in ['', '_r']:
            cmap2use = cmapOption + cmapSuffix
            print "\t" + cmap2use
            title = chrom + " " + outputLabel + " linear Q5-Q95 " + cmap2use
            outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".linear.Q5-Q95." + cmap2use
            plot_matplotlib_heatmap(np.clip(data,
                                            a_min=np.nanpercentile(data, 5),
                                            a_max=np.nanpercentile(data, 95)),
                                    theRez, title, outFilenameBase, cmap=cmap2use)

    ################################################
    # Logged data

    print "Logged"
    for cmapOption in cmapOptions:
        for cmapSuffix in ['', '_r']:
            cmap2use = cmapOption + cmapSuffix
            print "\t" + cmap2use
            title = chrom + " " + outputLabel + " logged " + cmap2use
            outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".logged." + cmap2use
            plot_matplotlib_heatmap(np.log([[z + 1 for z in y] for y in data]),
                                    theRez, title, outFilenameBase, cmap=cmap2use)


##############################################################################
# Print some of 'em (no LOI annotations)
def printSomeMaps(data, corData, chrom, outputLabel, outDir, theRez,
                  cmapSeq=defaultAbsCmap, cmapDiv=defaultCorrDivCmap):
    print "\nprintSomeMaps" + ", " + chrom + ", " + outputLabel + ", " + outDir \
          + ", " + str(theRez)
    theRez = int(theRez)

    dataDF = pd.DataFrame(data)

    ################################################
    # Clipped data

    print "Clipped"
    title = chrom + " " + outputLabel + " linear Q5-Q95"
    outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".linear.Q5-Q95"
    # plot_contact_map(np.clip(data,
    #                          a_min=np.nanpercentile(data, 5),
    #                          a_max=np.nanpercentile(data, 95)),
    #                  title, outFilenameBase)
    plot_matplotlib_heatmap(np.clip(data,
                                    a_min=np.nanpercentile(data, 5),
                                    a_max=np.nanpercentile(data, 95)),
                            theRez, title, outFilenameBase, cmap=cmapSeq)

    ################################################
    # Logged data

    print "Logged"
    title = chrom + " " + outputLabel + " logged"
    outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".logged"
    # data = [[100 for x in range(5)] for x in range(5)] # For testing
    # for i in range(min(len(data[0]),len(data))):
    #     data[i][i] = np.nan
    #     if i > 0:
    #         data[i][i-1] = np.nan
    #         data[i-1][i] = np.nan
    # # np.log(data+1) # Python is not as clever as R in this respect
    # [[z+1 for z in y] for y in data]
    # np.log([[z+1 for z in y] for y in data])
    #plot_contact_map(np.log([[z + 1 for z in y] for y in data]), title, outFilenameBase)
    plot_matplotlib_heatmap(np.log([[z + 1 for z in y] for y in data]),
                                     theRez, title, outFilenameBase, cmap=cmapSeq)

    ################################################
    # Correlation matrices

    print "Pearson correlation"
    # dataDF = pd.DataFrame(data)
    # mask = datadf.isnull()

    # 170623
    # # corData =  np.corrcoef(data)
    # # corData =  np.ma.corrcoef(data)
    # # corData = pd.DataFrame(data).corr
    # corData = pd.DataFrame.corr(dataDF)
    # # Rediagonal & off diags
    # # data = [[0 for x in range(5)] for x in range(5)] # For testing
    # for i in range(min(len(corData[0]), len(corData))):
    #     corData[i][i] = np.nan
    #     if i > 0:
    #         corData[i][i - 1] = np.nan
    #         corData[i - 1][i] = np.nan
    outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".corrPearson"
    title = chrom + outputLabel + " Pearson"

    # # plot_divergent_maps(corData, title, outFilenameBase)
    # plot_matplotlib_heatmap(corData, theRez, title, outFilenameBase, cmap=cmapDiv)
    # Symmetrical color map using vmin and vmax
    vmax = np.nanmax(corData)
    vmin = np.nanmin(corData)
    vabs = np.nanmax([abs(vmax), abs(vmin)])
    plot_matplotlib_heatmap(corData, theRez, title, outFilenameBase, cmap=cmapDiv, vmax=vabs, vmin=-vabs)

    '''
    Testing color map options:

    # Shifted midpoint:
    midRank = 1 - vmax/(vmax + abs(vmin))

    shifted_cmap = shiftedColorMap(cmapDiv, midpoint=midRank, name='shiftedcmapDiv')
    plot_matplotlib_heatmap(corData, theRez, title, outFilenameBase, cmap=shifted_cmap)

    # shifted midband:
    shifted_cmap = shiftedColorMapWithMidBand(cmapDiv,
                                              midBandStart=midRank-0.20,
                                              midBandStop=midRank+0.20,
                                              name='shiftedcmapDiv')
    plot_matplotlib_heatmap(corData, theRez, title, outFilenameBase, cmap=shifted_cmap)


    # shifted midband with vmin and vmax:
    midRank = 0.5
    shifted_cmap = shiftedColorMapWithMidBand(cmapDiv,
                                              midBandStart=midRank-0.25,
                                              midBandStop=midRank+0.25,
                                              name='shiftedcmapDiv')
    plot_matplotlib_heatmap(corData, theRez, title, outFilenameBase, cmap=shifted_cmap, vmax=vabs, vmin=-vabs)
    plot_matplotlib_heatmap(corData, theRez, title, outFilenameBase, cmap=cmapDiv, vmax=0.5, vmin=-0.5)
    '''

    # if theRez >= 500000:
    #     print "Spearman correlation"
    #     # dataRanks = sps.rankdata(data).reshape(len(data), len(data))
    #     # http://stackoverflow.com/questions/12519866/rank-array-in-python-while-ignoring-missing-values
    #     dataRanks = mstats.rankdata(np.ma.masked_invalid(data)).reshape(len(data), len(data))
    #     dataRanks[dataRanks == 0] = np.nan
    #     dataRanks -= 1
    #     dataRanksDF = pd.DataFrame(dataRanks)
    #     # spearcorData =  pd.DataFrame.corr(df, method="spearman")
    #     # spearcorData = pd.DataFrame(dataRanks).corr
    #     spearcorData = pd.DataFrame.corr(dataRanksDF)  # Pearson of ranks == Spearman
    #     # Rediagonal & off diags
    #     # data = [[0 for x in range(5)] for x in range(5)] # For testing
    #     for i in range(min(len(spearcorData[0]), len(spearcorData))):
    #         spearcorData[i][i] = np.nan
    #         if i > 0:
    #             spearcorData[i][i - 1] = np.nan
    #             spearcorData[i - 1][i] = np.nan
    #     outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".corrSpearman"
    #     title = chrom + outputLabel + " Spearman"
    #     # plot_divergent_maps(spearcorData, title, outFilenameBase)
    #     plot_matplotlib_heatmap(spearcorData, theRez, title, outFilenameBase, cmap=cmapDiv)

    # ################################################
    # # Clipped correlations
    #
    # if theRez < 500000:
    #     print "More Pearson correlation plots"
    #     for theMin in [-0.1]:
    #         #for theMax in frange(0.3, 0.6, 0.1):
    #         for theMax in [0.3]:
    #             rangeStr = str(int(round(theMin * 100))) + "to" + str(int(round(theMax * 100)))
    #             print rangeStr
    #
    #             title = chrom + " " + outputLabel + " Pearson " + rangeStr
    #             outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".corrPearson." + rangeStr
    #             # plot_divergent_maps(np.clip(corData, a_min=theMin, a_max=theMax), title, outFilenameBase)
    #             plot_matplotlib_heatmap(np.clip(corData, a_min=theMin, a_max=theMax),
    #                                              theRez, title, outFilenameBase, cmap=cmapDiv)


##############################################################################
# Print some of 'em with LOI annotations
def printSomeMapsWithLOI(data, corData, chrom, outputLabel, outDir, theRez, loi,
                         cmapSeq=defaultAbsCmap, cmapDiv=defaultCorrDivCmap):
    print "\nprintSomeMapsWithLOI" + ", " + chrom + ", " + outputLabel + ", " + outDir \
          + ", " + str(theRez)
    theRez = int(theRez)

    dataDF = pd.DataFrame(data)

    ################################################
    # Clipped data

    print "Clipped"
    title = chrom + " " + outputLabel + " linear Q5-Q95"
    outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".linear.Q5-Q95"
    # plot_contact_map(np.clip(data,
    #                          a_min=np.nanpercentile(data, 5),
    #                          a_max=np.nanpercentile(data, 95)),
    #                  title, outFilenameBase)
    plot_matplotlib_heatmap_with_LOI(np.clip(data,
                                             a_min=np.nanpercentile(data, 5),
                                             a_max=np.nanpercentile(data, 95)),
                                     loi, theRez, title, outFilenameBase, cmap=cmapSeq)

    ################################################
    # Logged data

    print "Logged"
    title = chrom + " " + outputLabel + " logged"
    outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".logged"
    # data = [[100 for x in range(5)] for x in range(5)] # For testing
    # for i in range(min(len(data[0]),len(data))):
    #     data[i][i] = np.nan
    #     if i > 0:
    #         data[i][i-1] = np.nan
    #         data[i-1][i] = np.nan
    # # np.log(data+1) # Python is not as clever as R in this respect
    # [[z+1 for z in y] for y in data]
    # np.log([[z+1 for z in y] for y in data])
    #plot_contact_map(np.log([[z + 1 for z in y] for y in data]), title, outFilenameBase)
    plot_matplotlib_heatmap_with_LOI(np.log([[z + 1 for z in y] for y in data]),
                                     loi, theRez, title, outFilenameBase, cmap=cmapSeq)

    ################################################
    # Correlation matrices

    print "Pearson correlation"
    # dataDF = pd.DataFrame(data)
    # mask = datadf.isnull()

    # 170623
    # # corData =  np.corrcoef(data)
    # # corData =  np.ma.corrcoef(data)
    # # corData = pd.DataFrame(data).corr
    # corData = pd.DataFrame.corr(dataDF)
    # # Rediagonal & off diags
    # # data = [[0 for x in range(5)] for x in range(5)] # For testing
    # for i in range(min(len(corData[0]), len(corData))):
    #     corData[i][i] = np.nan
    #     if i > 0:
    #         corData[i][i - 1] = np.nan
    #         corData[i - 1][i] = np.nan
    outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".corrPearson"
    title = chrom + outputLabel + " Pearson"
    # # plot_divergent_maps(corData, title, outFilenameBase)
    # plot_matplotlib_heatmap_with_LOI(corData,
    #                                  loi, theRez, title, outFilenameBase, cmap=cmapDiv)
    # Symmetrical color map using vmin and vmax
    vmax = np.nanmax(corData)
    vmin = np.nanmin(corData)
    vabs = np.nanmax([abs(vmax), abs(vmin)])
    plot_matplotlib_heatmap_with_LOI(corData, loi, theRez, title, outFilenameBase, cmap=cmapDiv, vmax=vabs, vmin=-vabs)

    # if theRez >= 500000:
    #     print "Spearman correlation"
    #     # dataRanks = sps.rankdata(data).reshape(len(data), len(data))
    #     # http://stackoverflow.com/questions/12519866/rank-array-in-python-while-ignoring-missing-values
    #     dataRanks = mstats.rankdata(np.ma.masked_invalid(data)).reshape(len(data), len(data))
    #     dataRanks[dataRanks == 0] = np.nan
    #     dataRanks -= 1
    #     dataRanksDF = pd.DataFrame(dataRanks)
    #     # spearcorData =  pd.DataFrame.corr(df, method="spearman")
    #     # spearcorData = pd.DataFrame(dataRanks).corr
    #     spearcorData = pd.DataFrame.corr(dataRanksDF)  # Pearson of ranks == Spearman
    #     # Rediagonal & off diags
    #     # data = [[0 for x in range(5)] for x in range(5)] # For testing
    #     for i in range(min(len(spearcorData[0]), len(spearcorData))):
    #         spearcorData[i][i] = np.nan
    #         if i > 0:
    #             spearcorData[i][i - 1] = np.nan
    #             spearcorData[i - 1][i] = np.nan
    #     outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".corrSpearman"
    #     title = chrom + outputLabel + " Spearman"
    #     # plot_divergent_maps(spearcorData, title, outFilenameBase)
    #     plot_matplotlib_heatmap_with_LOI(spearcorData,
    #                                      loi, theRez, title, outFilenameBase, cmap=cmapDiv)

    # ################################################
    # # Clipped correlations
    #
    # if theRez < 500000:
    #     print "More Pearson correlation plots"
    #     for theMin in [-0.1]:
    #         #for theMax in frange(0.3, 0.6, 0.1):
    #         for theMax in [0.3]:
    #             rangeStr = str(int(round(theMin * 100))) + "to" + str(int(round(theMax * 100)))
    #             print rangeStr
    #
    #             title = chrom + " " + outputLabel + " Pearson " + rangeStr
    #             outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".corrPearson." + rangeStr
    #             # plot_divergent_maps(np.clip(corData, a_min=theMin, a_max=theMax), title, outFilenameBase)
    #             plot_matplotlib_heatmap_with_LOI(np.clip(corData, a_min=theMin, a_max=theMax),
    #                                              loi, theRez, title, outFilenameBase, cmap=cmapDiv)


##############################################################################
# Print some of 'em rotated.
def printSomeMapsRot45(data, corData, chrom, outputLabel, outDir, theRez,
                       cmapSeq=defaultAbsCmap, cmapDiv=defaultCorrDivCmap):
    print "\nprintSomeMapsRot45" + ", " + chrom + ", " + outputLabel + ", " + outDir \
          + ", " + str(theRez)
    theRez = int(theRez)

    dataDF = pd.DataFrame(data)

    ################################################
    # Clipped data

    print "Clipped"
    title = chrom + " " + outputLabel + " linear Q5-Q95"
    outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".linear.Q5-Q95"
    plot_matplotlib_heatmap_r45(np.clip(data, a_min=np.nanpercentile(data, 5), a_max=np.nanpercentile(data, 95)),
                                theRez, title, outFilenameBase, cmap=cmapSeq)

    ################################################
    # Logged data

    print "Logged"
    title = chrom + " " + outputLabel + " logged"
    outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".logged"
    # data = [[100 for x in range(5)] for x in range(5)] # For testing
    # for i in range(min(len(data[0]),len(data))):
    #     data[i][i] = np.nan
    #     if i > 0:
    #         data[i][i-1] = np.nan
    #         data[i-1][i] = np.nan
    # # np.log(data+1) # Python is not as clever as R in this respect
    # [[z+1 for z in y] for y in data]
    # np.log([[z+1 for z in y] for y in data])
    #plot_contact_map(np.log([[z + 1 for z in y] for y in data]), title, outFilenameBase)
    plot_matplotlib_heatmap_r45(np.log([[z + 1 for z in y] for y in data]),
                                theRez, title, outFilenameBase, cmap=cmapSeq)

    ################################################
    # Correlation matrices

    print "Pearson correlation"
    # dataDF = pd.DataFrame(data)
    # mask = datadf.isnull()

    # 170623
    # # corData =  np.corrcoef(data)
    # # corData =  np.ma.corrcoef(data)
    # # corData = pd.DataFrame(data).corr
    # corData = pd.DataFrame.corr(dataDF)
    # # Rediagonal & off diags
    # # data = [[0 for x in range(5)] for x in range(5)] # For testing
    # for i in range(min(len(corData[0]), len(corData))):
    #     corData[i][i] = np.nan
    #     if i > 0:
    #         corData[i][i - 1] = np.nan
    #         corData[i - 1][i] = np.nan
    outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".corrPearson"
    title = chrom + outputLabel + " Pearson"
    # # plot_divergent_maps(corData, title, outFilenameBase)
    # plot_matplotlib_heatmap_r45(corData,theRez, title, outFilenameBase, cmap=cmapDiv)
    # Symmetrical color map using vmin and vmax
    vmax = np.nanmax(corData)
    vmin = np.nanmin(corData)
    vabs = np.nanmax([abs(vmax), abs(vmin)])
    plot_matplotlib_heatmap_r45(corData, theRez, title, outFilenameBase, cmap=cmapDiv, vmax=vabs, vmin=-vabs)

    # if theRez >= 500000:
    #     print "Spearman correlation"
    #     # dataRanks = sps.rankdata(data).reshape(len(data), len(data))
    #     # http://stackoverflow.com/questions/12519866/rank-array-in-python-while-ignoring-missing-values
    #     dataRanks = mstats.rankdata(np.ma.masked_invalid(data)).reshape(len(data), len(data))
    #     dataRanks[dataRanks == 0] = np.nan
    #     dataRanks -= 1
    #     dataRanksDF = pd.DataFrame(dataRanks)
    #     # spearcorData =  pd.DataFrame.corr(df, method="spearman")
    #     # spearcorData = pd.DataFrame(dataRanks).corr
    #     spearcorData = pd.DataFrame.corr(dataRanksDF)  # Pearson of ranks == Spearman
    #     # Rediagonal & off diags
    #     # data = [[0 for x in range(5)] for x in range(5)] # For testing
    #     for i in range(min(len(spearcorData[0]), len(spearcorData))):
    #         spearcorData[i][i] = np.nan
    #         if i > 0:
    #             spearcorData[i][i - 1] = np.nan
    #             spearcorData[i - 1][i] = np.nan
    #     outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".corrSpearman"
    #     title = chrom + outputLabel + " Spearman"
    #     # plot_divergent_maps(spearcorData, title, outFilenameBase)
    #     plot_matplotlib_heatmap_with_LOI_r45(spearcorData,
    #                                      loi, theRez, title, outFilenameBase, cmap=cmapDiv)

    # ################################################
    # # Clipped correlations
    #
    # if theRez < 500000:
    #     print "More Pearson correlation plots"
    #     for theMin in [-0.1]:
    #         #for theMax in frange(0.3, 0.6, 0.1):
    #         for theMax in [0.3]:
    #             rangeStr = str(int(round(theMin * 100))) + "to" + str(int(round(theMax * 100)))
    #             print rangeStr
    #
    #             title = chrom + " " + outputLabel + " Pearson " + rangeStr
    #             outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".corrPearson." + rangeStr
    #             # plot_divergent_maps(np.clip(corData, a_min=theMin, a_max=theMax), title, outFilenameBase)
    #             plot_matplotlib_heatmap_with_LOI_r45(np.clip(corData, a_min=theMin, a_max=theMax),
    #                                              loi, theRez, title, outFilenameBase, cmap=cmapDiv)


##############################################################################
# Print some of 'em with loci of interest and rotated.
def printSomeMapsWithLOIrot45(data, corData, chrom, outputLabel, outDir, theRez, loi,
                              cmapSeq=defaultAbsCmap, cmapDiv=defaultCorrDivCmap):
    print "\nprintSomeMapsWithLOIrot45" + ", " + chrom + ", " + outputLabel + ", " + outDir \
          + ", " + str(theRez)
    theRez = int(theRez)

    dataDF = pd.DataFrame(data)

    ################################################
    # Clipped data

    print "Clipped"
    title = chrom + " " + outputLabel + " linear Q5-Q95"
    outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".linear.Q5-Q95"
    # plot_contact_map(np.clip(data,
    #                          a_min=np.nanpercentile(data, 5),
    #                          a_max=np.nanpercentile(data, 95)),
    #                  title, outFilenameBase)
    plot_matplotlib_heatmap_with_LOI_r45(np.clip(data,
                                             a_min=np.nanpercentile(data, 5),
                                             a_max=np.nanpercentile(data, 95)),
                                     loi, theRez, title, outFilenameBase, cmap=cmapSeq)

    ################################################
    # Logged data

    print "Logged"
    title = chrom + " " + outputLabel + " logged"
    outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".logged"
    # data = [[100 for x in range(5)] for x in range(5)] # For testing
    # for i in range(min(len(data[0]),len(data))):
    #     data[i][i] = np.nan
    #     if i > 0:
    #         data[i][i-1] = np.nan
    #         data[i-1][i] = np.nan
    # # np.log(data+1) # Python is not as clever as R in this respect
    # [[z+1 for z in y] for y in data]
    # np.log([[z+1 for z in y] for y in data])
    #plot_contact_map(np.log([[z + 1 for z in y] for y in data]), title, outFilenameBase)
    plot_matplotlib_heatmap_with_LOI_r45(np.log([[z + 1 for z in y] for y in data]),
                                     loi, theRez, title, outFilenameBase, cmap=cmapSeq)

    ################################################
    # Correlation matrices

    print "Pearson correlation"
    # dataDF = pd.DataFrame(data)
    # mask = datadf.isnull()

    # 170623
    # #  corData =  np.corrcoef(data)
    # # corData =  np.ma.corrcoef(data)
    # # corData = pd.DataFrame(data).corr
    # corData = pd.DataFrame.corr(dataDF)
    # # Rediagonal & off diags
    # # data = [[0 for x in range(5)] for x in range(5)] # For testing
    # for i in range(min(len(corData[0]), len(corData))):
    #     corData[i][i] = np.nan
    #     if i > 0:
    #         corData[i][i - 1] = np.nan
    #         corData[i - 1][i] = np.nan
    outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".corrPearson"
    title = chrom + outputLabel + " Pearson"
    # # plot_divergent_maps(corData, title, outFilenameBase)
    # plot_matplotlib_heatmap_with_LOI_r45(corData,
    #                                  loi, theRez, title, outFilenameBase, cmap=cmapDiv)
    # Symmetrical color map using vmin and vmax
    vmax = np.nanmax(corData)
    vmin = np.nanmin(corData)
    vabs = np.nanmax([abs(vmax), abs(vmin)])
    plot_matplotlib_heatmap_with_LOI_r45(corData, loi, theRez, title, outFilenameBase, cmap=cmapDiv, vmax=vabs, vmin=-vabs)


    # if theRez >= 500000:
    #     print "Spearman correlation"
    #     # dataRanks = sps.rankdata(data).reshape(len(data), len(data))
    #     # http://stackoverflow.com/questions/12519866/rank-array-in-python-while-ignoring-missing-values
    #     dataRanks = mstats.rankdata(np.ma.masked_invalid(data)).reshape(len(data), len(data))
    #     dataRanks[dataRanks == 0] = np.nan
    #     dataRanks -= 1
    #     dataRanksDF = pd.DataFrame(dataRanks)
    #     # spearcorData =  pd.DataFrame.corr(df, method="spearman")
    #     # spearcorData = pd.DataFrame(dataRanks).corr
    #     spearcorData = pd.DataFrame.corr(dataRanksDF)  # Pearson of ranks == Spearman
    #     # Rediagonal & off diags
    #     # data = [[0 for x in range(5)] for x in range(5)] # For testing
    #     for i in range(min(len(spearcorData[0]), len(spearcorData))):
    #         spearcorData[i][i] = np.nan
    #         if i > 0:
    #             spearcorData[i][i - 1] = np.nan
    #             spearcorData[i - 1][i] = np.nan
    #     outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".corrSpearman"
    #     title = chrom + outputLabel + " Spearman"
    #     # plot_divergent_maps(spearcorData, title, outFilenameBase)
    #     plot_matplotlib_heatmap_with_LOI_r45(spearcorData,
    #                                      loi, theRez, title, outFilenameBase, cmap=cmapDiv)

    # ################################################
    # # Clipped correlations
    #
    # if theRez < 500000:
    #     print "More Pearson correlation plots"
    #     for theMin in [-0.1]:
    #         #for theMax in frange(0.3, 0.6, 0.1):
    #         for theMax in [0.3]:
    #             rangeStr = str(int(round(theMin * 100))) + "to" + str(int(round(theMax * 100)))
    #             print rangeStr
    #
    #             title = chrom + " " + outputLabel + " Pearson " + rangeStr
    #             outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".corrPearson." + rangeStr
    #             # plot_divergent_maps(np.clip(corData, a_min=theMin, a_max=theMax), title, outFilenameBase)
    #             plot_matplotlib_heatmap_with_LOI_r45(np.clip(corData, a_min=theMin, a_max=theMax),
    #                                              loi, theRez, title, outFilenameBase, cmap=cmapDiv)


##############################################################################
# Print some of 'em with loci of interest, rotated and zoomed!
def printSomeMapsWithLOIrot45andZoom(data, corData, chrom, outputLabel, outDir, theRez, loi, zoomWin,
                                     cmapSeq=defaultAbsCmap, cmapDiv=defaultCorrDivCmap):
    print "\nprintSomeMapsWithLOIandZoomAndRot45" + ", " + chrom + ", " + outputLabel + ", " + outDir + ", " \
          + str(theRez) + ", "  + str(zoomWin)
    theRez = int(theRez); zoomWin = int(zoomWin)
	
    dataDF = pd.DataFrame(data)

    ################################################
    # Clipped data

    print "Clipped"
    title = chrom + " " + outputLabel + " linear Q5-Q99"
    outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".linear.Q5-Q99"
    plot_matplotlib_heatmap_with_LOI_r45_zoom(np.clip(data,
                                             a_min=np.nanpercentile(data, 5),
                                             a_max=np.nanpercentile(data, 99)),
                                     loi, theRez, zoomWin, title, outFilenameBase, cmap=cmapSeq)

    ################################################
    # Logged data

    print "Logged"
    title = chrom + " " + outputLabel + " logged"
    outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".logged"
    plot_matplotlib_heatmap_with_LOI_r45_zoom(np.log([[z + 1 for z in y] for y in data]),
                                     loi, theRez, zoomWin, title, outFilenameBase, cmap=cmapSeq)

    ################################################
    # Correlation matrices

    print "Pearson correlation"
    # dataDF = pd.DataFrame(data)
    # mask = datadf.isnull()

    # 170623
    # #  corData =  np.corrcoef(data)
    # # corData =  np.ma.corrcoef(data)
    # # corData = pd.DataFrame(data).corr
    # corData = pd.DataFrame.corr(dataDF)
    # # Rediagonal & off diags
    # # data = [[0 for x in range(5)] for x in range(5)] # For testing
    # for i in range(min(len(corData[0]), len(corData))):
    #     corData[i][i] = np.nan
    #     if i > 0:
    #         corData[i][i - 1] = np.nan
    #         corData[i - 1][i] = np.nan
    outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".corrPearson"
    title = chrom + outputLabel + " Pearson"
    # # plot_divergent_maps(corData, title, outFilenameBase)
    # plot_matplotlib_heatmap_with_LOI_r45_zoom(corData,
    #                                  loi, theRez, zoomWin, title, outFilenameBase, cmap=cmapDiv)
    # Symmetrical color map using vmin and vmax
    vmax = np.nanmax(corData)
    vmin = np.nanmin(corData)
    vabs = np.nanmax([abs(vmax), abs(vmin)])
    plot_matplotlib_heatmap_with_LOI_r45_zoom(corData, loi, theRez, zoomWin, title, outFilenameBase, cmap=cmapDiv, vmax=vabs, vmin=-vabs)



    # if theRez >= 500000:
    #     print "Spearman correlation"
    #     # dataRanks = sps.rankdata(data).reshape(len(data), len(data))
    #     # http://stackoverflow.com/questions/12519866/rank-array-in-python-while-ignoring-missing-values
    #     dataRanks = mstats.rankdata(np.ma.masked_invalid(data)).reshape(len(data), len(data))
    #     dataRanks[dataRanks == 0] = np.nan
    #     dataRanks -= 1
    #     dataRanksDF = pd.DataFrame(dataRanks)
    #     # spearcorData =  pd.DataFrame.corr(df, method="spearman")
    #     # spearcorData = pd.DataFrame(dataRanks).corr
    #     spearcorData = pd.DataFrame.corr(dataRanksDF)  # Pearson of ranks == Spearman
    #     # Rediagonal & off diags
    #     # data = [[0 for x in range(5)] for x in range(5)] # For testing
    #     for i in range(min(len(spearcorData[0]), len(spearcorData))):
    #         spearcorData[i][i] = np.nan
    #         if i > 0:
    #             spearcorData[i][i - 1] = np.nan
    #             spearcorData[i - 1][i] = np.nan
    #     outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".corrSpearman"
    #     title = chrom + outputLabel + " Spearman"
    #     # plot_divergent_maps(spearcorData, title, outFilenameBase)
    #     plot_matplotlib_heatmap_with_LOI_r45_zoom(spearcorData,
    #                                      loi, theRez, zoomWin, title, outFilenameBase, cmap=cmapDiv)

    # ################################################
    # # Clipped correlations
    #
    # if theRez < 500000:
    #     print "More Pearson correlation plots"
    #     for theMin in [-0.1]:
    #         #for theMax in frange(0.3, 0.6, 0.1):
    #         for theMax in [0.3]:
    #             rangeStr = str(int(round(theMin * 100))) + "to" + str(int(round(theMax * 100)))
    #             print rangeStr
    #
    #             title = chrom + " " + outputLabel + " Pearson " + rangeStr
    #             outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".corrPearson." + rangeStr
    #             # plot_divergent_maps(np.clip(corData, a_min=theMin, a_max=theMax), title, outFilenameBase)
    #             plot_matplotlib_heatmap_with_LOI_zoom(np.clip(corData, a_min=theMin, a_max=theMax),
    #                                              loi, theRez, zoomWin, title, outFilenameBase, cmap=cmapDiv)


##############################################################################
# 161102
# Print some of 'em with loci of interest and Zoom!
def printSomeMapsWithLOIandZoom(data, corData, chrom, outputLabel, outDir, theRez, loi, zoomWin,
                                cmapSeq=defaultAbsCmap, cmapDiv=defaultCorrDivCmap):
    print "\nprintSomeMapsWithLOIandZoom" + ", " + chrom + ", " + outputLabel + ", " + outDir + ", " \
          + str(theRez) + ", "  + str(zoomWin)
    theRez = int(theRez); zoomWin = int(zoomWin)

    dataDF = pd.DataFrame(data)

    ################################################
    # Clipped data

    print "Clipped"
    title = chrom + " " + outputLabel + " linear Q5-Q99"
    outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".linear.Q5-Q99"
    plot_matplotlib_heatmap_with_LOI_zoom(np.clip(data,
                                             a_min=np.nanpercentile(data, 5),
                                             a_max=np.nanpercentile(data, 99)),
                                     loi, theRez, zoomWin, title, outFilenameBase, cmap=cmapSeq)

    ################################################
    # Logged data

    print "Logged"
    title = chrom + " " + outputLabel + " logged"
    outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".logged"
    plot_matplotlib_heatmap_with_LOI_zoom(np.log([[z + 1 for z in y] for y in data]),
                                     loi, theRez, zoomWin, title, outFilenameBase, cmap=cmapSeq)

    ################################################
    # Correlation matrices

    print "Pearson correlation"
    # dataDF = pd.DataFrame(data)
    # mask = datadf.isnull()

    # 170623
    # #  corData =  np.corrcoef(data)
    # # corData =  np.ma.corrcoef(data)
    # # corData = pd.DataFrame(data).corr
    # corData = pd.DataFrame.corr(dataDF)
    # # Rediagonal & off diags
    # # data = [[0 for x in range(5)] for x in range(5)] # For testing
    # for i in range(min(len(corData[0]), len(corData))):
    #     corData[i][i] = np.nan
    #     if i > 0:
    #         corData[i][i - 1] = np.nan
    #         corData[i - 1][i] = np.nan
    outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".corrPearson"
    title = chrom + outputLabel + " Pearson"
    # # plot_divergent_maps(corData, title, outFilenameBase)
    # plot_matplotlib_heatmap_with_LOI_zoom(corData,
    #                                  loi, theRez, zoomWin, title, outFilenameBase, cmap=cmapDiv)
    # Symmetrical color map using vmin and vmax
    vmax = np.nanmax(corData)
    vmin = np.nanmin(corData)
    vabs = np.nanmax([abs(vmax), abs(vmin)])
    plot_matplotlib_heatmap_with_LOI_zoom(corData, loi, theRez, zoomWin, title, outFilenameBase, cmap=cmapDiv, vmax=vabs, vmin=-vabs)

    # if theRez >= 500000:
    #     print "Spearman correlation"
    #     # dataRanks = sps.rankdata(data).reshape(len(data), len(data))
    #     # http://stackoverflow.com/questions/12519866/rank-array-in-python-while-ignoring-missing-values
    #     dataRanks = mstats.rankdata(np.ma.masked_invalid(data)).reshape(len(data), len(data))
    #     dataRanks[dataRanks == 0] = np.nan
    #     dataRanks -= 1
    #     dataRanksDF = pd.DataFrame(dataRanks)
    #     # spearcorData =  pd.DataFrame.corr(df, method="spearman")
    #     # spearcorData = pd.DataFrame(dataRanks).corr
    #     spearcorData = pd.DataFrame.corr(dataRanksDF)  # Pearson of ranks == Spearman
    #     # Rediagonal & off diags
    #     # data = [[0 for x in range(5)] for x in range(5)] # For testing
    #     for i in range(min(len(spearcorData[0]), len(spearcorData))):
    #         spearcorData[i][i] = np.nan
    #         if i > 0:
    #             spearcorData[i][i - 1] = np.nan
    #             spearcorData[i - 1][i] = np.nan
    #     outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".corrSpearman"
    #     title = chrom + outputLabel + " Spearman"
    #     # plot_divergent_maps(spearcorData, title, outFilenameBase)
    #     plot_matplotlib_heatmap_with_LOI_zoom(spearcorData,
    #                                      loi, theRez, zoomWin, title, outFilenameBase, cmap=cmapDiv)

    # ################################################
    # # Clipped correlations
    #
    # if theRez < 500000:
    #     print "More Pearson correlation plots"
    #     for theMin in [-0.1]:
    #         #for theMax in frange(0.3, 0.6, 0.1):
    #         for theMax in [0.3]:
    #             rangeStr = str(int(round(theMin * 100))) + "to" + str(int(round(theMax * 100)))
    #             print rangeStr
    #
    #             title = chrom + " " + outputLabel + " Pearson " + rangeStr
    #             outFilenameBase = outDir + "/" + outputLabel + "." + chrom + ".corrPearson." + rangeStr
    #             # plot_divergent_maps(np.clip(corData, a_min=theMin, a_max=theMax), title, outFilenameBase)
    #             plot_matplotlib_heatmap_with_LOI_zoom(np.clip(corData, a_min=theMin, a_max=theMax),
    #                                              loi, theRez, zoomWin, title, outFilenameBase, cmap=cmapDiv')


##############################################################################
#############################################################################
#############################################################################
# MAIN

##############################################################################
#  Parse the command line.
if len(sys.argv) != 2:
        #and len(sys.argv) != 3:
    sys.stderr.write(USAGE)
    sys.exit(1)
configFileName = sys.argv[1]
# if len(sys.argv) == 3:
#     plotRawMaps = str2bool(sys.argv[2])
# else:
#     plotRawMaps = False


##############################################################################
# For testing
# os.getcwd()
# os.chdir('PatskiWTvsDelvsDxz4InvFixed'); configFileName = 'configs/heatmaps.500000.chr20.cfg'
# plotRawMaps = False


##############################################################################
#  Parse the config file.

print "\n\n\nParse config"
config = ConfigParser.RawConfigParser()
config.read(configFileName)

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
contactMapOutputDesc = config.get('general', 'desc')
assembly = config.get('general', 'assembly')
chrom = config.get('general', 'chrom')
theRez = config.get('general', 'rez')
intTheRez = int(theRez)
sectionTally+=1

# -- What to run --
toRun = {}
for option in config.options('toRun'):
    toRun[option] = str2bool(config.get('toRun', option))
sectionTally+=1


# -- LOI --
LOI = {}
for option in config.options('LOI'):
    LOI[option] = config.get('LOI', option)
sectionTally+=1
# LOI

# -- Contact Frequency Data --
dataSets = config.sections()[sectionTally:-1]
contactMapIOmetadata = {}
for whatData in dataSets:
    contactMapIOmetadata[whatData] = {}
    contactMapIOmetadata[whatData]['inputDir'] = config.get(whatData, 'inputDir')
    contactMapIOmetadata[whatData]['inputFile'] = config.get(whatData, 'inputFile')
    contactMapIOmetadata[whatData]['outputDir'] = config.get(whatData, 'outputDir')
    contactMapIOmetadata[whatData]['outputLabel'] = config.get(whatData, 'outputLabel')
    contactMapIOmetadata[whatData]['dataColor'] = config.get(whatData, 'dataColor')
    contactMapIOmetadata[whatData]['dataLineWidth'] = config.get(whatData, 'dataLineWidth')
    contactMapIOmetadata[whatData]['dataLineStyle'] = config.get(whatData, 'dataLineStyle')
    sectionTally+=1

# contactMapIOmetadata
# -- Differential Maps --
diffMaps = {}
for option in config.options('diffMaps'):
    diffMaps[option] = config.get('diffMaps', option)
sectionTally+=1
# diffMaps


################################################
################################################
################################################
# Other variables and setup

# Use TrueType fonts that work in illustrator:
#    http://stackoverflow.com/questions/5956182/cannot-edit-text-in-chart-exported-by-matplotlib-and-opened-in-illustrator
matplotlib.rcParams['pdf.fonttype'] = 42


# For plotting & chrom sizes & loi
thechrom = chrom
if assembly[0:2] == 'hg' and chrom == 'chr23':
    thechrom = 'chrX'
elif assembly[0:2] == 'mm' and chrom == 'chr20':
    thechrom = 'chrX'


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
# load chrX contact maps

print "\n\n\nLoad data"

dataDict = {}  # Dictionary of contact data

for whatData in dataSets:
    inputDir = contactMapIOmetadata[whatData]["inputDir"]
    inputFile = contactMapIOmetadata[whatData]["inputFile"]
    outputDir = contactMapIOmetadata[whatData]["outputDir"]
    outputLabel = contactMapIOmetadata[whatData]["outputLabel"]

    print whatData + ": " + inputDir + ", " + inputFile + ", " + outputDir + ", " + outputLabel

    ################################################
    # load data

    # if fileType(inputFile) == "gz":
    #     contactData = gzip.open(inputFile, "r")
    # else:
    #     contactData = open(inputFile, "r")
    # #data = np.loadtxt(inputFile)
    # data = np.array( contactData.read())
    data = np.loadtxt(inputDir + "/" + inputFile)
    print data.shape
    print np.mean(data)
    # print np.nanmean(data)

    # Millionify the data
    sum = np.nansum(data)
    mdata = np.array([[z / sum * 1000000 for z in y] for y in data])
    print "millionify:", sum, " --> ", np.nansum(mdata)
    data = mdata
    print data.shape
    print np.mean(data)

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
    # # plot_heatmaps_selectionA(filldata, title, outFilenameBase)


################################################
################################################
################################################
# Plot heatmaps

if 'plotrawmaps' in toRun and toRun['plotrawmaps']:
# if plotRawMaps:

    print "\n\n\nRaw heatmaps"
    for whatData in dataSets:
        # whatData = 'WTrefIC'

        inputDir = contactMapIOmetadata[whatData]["inputDir"]
        inputFile = contactMapIOmetadata[whatData]["inputFile"]
        outputDir = contactMapIOmetadata[whatData]["outputDir"]
        outputLabel = contactMapIOmetadata[whatData]["outputLabel"]

        # print whatData + ": " + inputDir + ", " + inputFile + ", " + outputDir + ", " + outputLabel

        # The Data!
        data = dataDict[whatData].copy()  # VERY IMPORTANT TO COPY OTHERWISE dataDict gets changed!
        print data.shape
        print np.mean(data)  # Note: data should have no NaNs so can just use np.mean (a.o.t. np.nanmean)

        # Unmappable regions: rows/cols with all zero counts
        dataDF = pd.DataFrame(data)
        unMappableRegions = np.where((dataDF.apply(pd.Series.nunique, axis=1) == 1) & (dataDF[1] == 0))
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

        ################################################
        # Output folder

        outDir = "output/" + contactMapOutputDesc + "/" + outputDir + "/heatmaps"
        if not os.path.exists(outDir):
            os.makedirs(outDir)

        ################################################
        # Print 'em all!
        printAllMapsAndPerformCorrelation(data, chrom, outputLabel, outDir, theRez)



################################################
################################################
################################################
# Quantile normalization across samples

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
# NaNorator
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
    # data = dataDict_QN[whatData].copy()  # VERY IMPORTANT TO COPY OTHERWISE dataDict gets changed!
    #print data.shape
    #print np.mean(data)  # Note: data should have no NaNs so can just use np.mean (a.o.t. np.nanmean)

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
for whatData in dataSets:
    # The QN Data!
    data = dataDict_QN[whatData].copy()  # VERY IMPORTANT TO COPY OTHERWISE dataDict gets changed!

    print np.mean(data)

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

    print np.nanmean(data)

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

print "\n\n\nGenerate correlation matrices\n"

corDataDict_QN_NaN = {}

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

    corDataDict_QN_NaN[whatData] = corData.as_matrix()



################################################
################################################
################################################
# Print some different color map options

if 'colormapoptions' in toRun and toRun['colormapoptions']:

    print "\n\n\nQuantile normalized heatmaps"
    for whatData in dataSets:
        # whatData = 'pooledWT_ref_IC'
        # whatData = 'pooledWT_alt_IC'

        inputDir = contactMapIOmetadata[whatData]["inputDir"]
        inputFile = contactMapIOmetadata[whatData]["inputFile"]
        outputDir = contactMapIOmetadata[whatData]["outputDir"]
        outputLabel = contactMapIOmetadata[whatData]["outputLabel"]

        # print whatData + ": " + inputDir + ", " + inputFile + ", " + outputDir + ", " + outputLabel

        # The QN Data!
        data = dataDict_QN_NaN[whatData].copy()  # VERY IMPORTANT TO COPY OTHERWISE dataDict gets changed!
        print data.shape
        print np.nanmean(data)

        outDir = "output/" + contactMapOutputDesc + "/" + outputDir + "/heatmaps_QN/vanillaColorMapOptions"
        if not os.path.exists(outDir):
            os.makedirs(outDir)
        printSomeColorMapOptions(data, chrom, outputLabel, outDir, theRez)



################################################
################################################
################################################
# PLot NaNed QN heatmaps

if 'vanillamaps' in toRun and toRun['vanillamaps'] and intTheRez >= 100000:

    ################################################
    # 161102
    # *** NOTE ***
    # Calling each of the following functions is inefficient
    # b/c each data transformation has to be applied again each time.
    # It might be better to bin functions by transformation type, not plot type.


    ########################
    # Print some full maps!
    # Print all full maps first b/c rotation messed up full map plotting for some reason

    print "\n\n\nQuantile normalized heatmaps"
    for whatData in dataSets:
        # whatData = 'Del_Xi'

        inputDir = contactMapIOmetadata[whatData]["inputDir"]
        inputFile = contactMapIOmetadata[whatData]["inputFile"]
        outputDir = contactMapIOmetadata[whatData]["outputDir"]
        outputLabel = contactMapIOmetadata[whatData]["outputLabel"]

        # print whatData + ": " + inputDir + ", " + inputFile + ", " + outputDir + ", " + outputLabel

        # The QN Data!
        data = dataDict_QN_NaN[whatData].copy()  # VERY IMPORTANT TO COPY OTHERWISE dataDict gets changed!
        print data.shape
        print np.nanmean(data)
        # 170623
        corData = corDataDict_QN_NaN[whatData].copy()  # VERY IMPORTANT TO COPY OTHERWISE dataDict gets changed!

        # *** THIS IS NO LONGER NECESSARY B/C DONE ON ALL DATA SETS SIMULTANEOUSLY ABOVE ***
        # # The RAW Data!
        # # Need  this to define unmappable regions.
        # data_raw = dataDict[whatData].copy()  # VERY IMPORTANT TO COPY OTHERWISE dataDict gets changed!
        # print data_raw.shape
        # print np.mean(data_raw)  # Note: data should have no NaNs so can just use np.mean (a.o.t. np.nanmean)
        #
        # # Unmappable regions: rows/cols with all zero counts
        # dataDF = pd.DataFrame(data_raw)
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

        # ################################################
        # # Print 'em all!
        # outDir = "output/" + contactMapOutputDesc + "/" + outputDir + "/heatmaps_QN/unannotated"
        # if not os.path.exists(outDir):
        #     os.makedirs(outDir)
        # printAllMaps(data, corData, chrom, outputLabel, outDir, theRez)

        ################################################
        # Print some full maps!
        outDir = "output/" + contactMapOutputDesc + "/" + outputDir + "/heatmaps_QN/vanilla"
        if not os.path.exists(outDir):
            os.makedirs(outDir)
        printSomeMaps(data, corData, chrom, outputLabel, outDir, theRez)

        # ################################################
        # # Print some annotated maps!
        # if len(LOI) >= 0:
        #     outDir = "output/" + contactMapOutputDesc + "/" + outputDir + "/heatmaps_QN/annotated"
        #     if not os.path.exists(outDir):
        #         os.makedirs(outDir)
        #     printSomeMapsWithLOI(data, corData, chrom, outputLabel, outDir, theRez, LOI)

        # ################################################
        # # Print some (unrotated) zoomed maps per LOI!
        # # 161102
        # if intTheRez < 500000 and len(LOI) >= 0:
        #     outDir = "output/" + contactMapOutputDesc + "/" + outputDir + "/heatmaps_QN/zoom100"
        #     if not os.path.exists(outDir):
        #         os.makedirs(outDir)
        #     printSomeMapsWithLOIandZoom(data, corData, chrom, outputLabel, outDir, theRez, LOI, intTheRez*100)



################################################
################################################
################################################
# Print some rotated maps!

if 'rot45maps' in toRun and toRun['rot45maps']:

    print "\n\n\nQuantile normalized heatmaps - rotated"
    for whatData in dataSets:
        # whatData = 'pooledWT_ref_IC'
        # whatData = 'pooledWT_alt_IC'

        inputDir = contactMapIOmetadata[whatData]["inputDir"]
        inputFile = contactMapIOmetadata[whatData]["inputFile"]
        outputDir = contactMapIOmetadata[whatData]["outputDir"]
        outputLabel = contactMapIOmetadata[whatData]["outputLabel"]

        # print whatData + ": " + inputDir + ", " + inputFile + ", " + outputDir + ", " + outputLabel

        # The QN Data!
        data = dataDict_QN_NaN[whatData].copy()  # VERY IMPORTANT TO COPY OTHERWISE dataDict gets changed!
        print data.shape
        print np.nanmean(data)
        # 170623
        corData = corDataDict_QN_NaN[whatData].copy()  # VERY IMPORTANT TO COPY OTHERWISE dataDict gets changed!

        ################################################
        # Print some rotated unannotated maps!
        if intTheRez >= 100000:
            outDir = "output/" + contactMapOutputDesc + "/" + outputDir + "/heatmaps_QN/rot45"
            if not os.path.exists(outDir):
                os.makedirs(outDir)
            printSomeMapsRot45(data, corData, chrom, outputLabel, outDir, theRez)

        # ################################################
        # # Print some rotated annotated maps with LOI!
        # if len(LOI) >= 0:
        #     outDir = "output/" + contactMapOutputDesc + "/" + outputDir + "/heatmaps_QN/annotatedRot45"
        #     if not os.path.exists(outDir):
        #         os.makedirs(outDir)
        #     printSomeMapsWithLOIrot45(data, corData, chrom, outputLabel, outDir, theRez, LOI)

        ################################################
        # Print some zoomed and rot45 maps per LOI!
        if len(LOI) >= 0:
            outDir = "output/" + contactMapOutputDesc + "/" + outputDir + "/heatmaps_QN/zoom100Rot45"
            if not os.path.exists(outDir):
                os.makedirs(outDir)
            printSomeMapsWithLOIrot45andZoom(data, corData, chrom, outputLabel, outDir, theRez, LOI, intTheRez*100)



################################################
################################################
################################################
# PLot Differential heatmaps

if 'diffmaps' in toRun and toRun['diffmaps']:

    print "\n\n\nDifferential  heatmaps"

    diffDir = "output/" + contactMapOutputDesc + "/" + outputDir + "/heatmaps_QN/differentialPlots"
    if not os.path.exists(diffDir):
        os.makedirs(diffDir)

    allSums = []
    allDeltas = []
    allLog2ratios = []
    for diffMapID in range(len(diffMaps)):
        #print str(diffMapID) + ": " + diffMaps[str(diffMapID)]
        [ diffMapA, diffMapB ] = diffMaps[str(diffMapID)].split('<<<--->>>')

        # The QN Data!
        dataA = dataDict_QN_NaN[diffMapA].copy()  # VERY IMPORTANT TO COPY OTHERWISE dataDict gets changed!
        dataB = dataDict_QN_NaN[diffMapB].copy()  # VERY IMPORTANT TO COPY OTHERWISE dataDict gets changed!

        sumMap = dataA + dataB
        allSums.extend(sumMap.flatten())

        deltaDiffMap = dataA - dataB
        #print('\n')
        #print(stats.describe(deltaDiffMap.flatten(), nan_policy='omit'))
        #print('\n')
        allDeltas.extend(deltaDiffMap.flatten())

        # log the deltas
        log10deltaDiffMap = np.log10(np.abs(deltaDiffMap))
        # Remove fractional deltas
        log10deltaDiffMap[np.abs(deltaDiffMap)<1] = np.nan
        # Return sign to logged data
        log10deltaDiffMap[deltaDiffMap<0] = log10deltaDiffMap[deltaDiffMap<0] * -1

        log2ratioDiffMap = np.log2( dataA / dataB )
        allLog2ratios.extend(np.ndarray.flatten(log2ratioDiffMap))

    # totals stats
    allSums = np.array(allSums)
    log10variance = np.nanvar(np.log10(allSums[allSums>0]))
    lognormSigma = 10**(log10variance ** 0.5)
    # clip outliers and recalculate - try this b/c at higher rez, more outliers skew variance
    allSumsClipped = allSums
    allSumsClipped[allSums > lognormSigma * 4] = lognormSigma * 4
    log10varianceClipped = np.nanvar(np.log10(allSumsClipped[allSumsClipped>0]))
    lognormSigmaClipped = 10**(log10varianceClipped ** 0.5)
    minSum = lognormSigmaClipped # 1 sd based on 4-sigma-clipped lognormalized counts

    # Or filter by quantile levels
    percentiles = [25, 50, 75]

    #diffLim = 25
    #[nobs, minmax, mean, variance, skewness, kurtosis] = stats.describe(allDeltas, nan_policy='omit')

    #diffLim = np.nanmax([abs(minmax[0]), minmax[1]])
    #diffLim = np.nanmax([abs(minmax[0]), minmax[1]])*0.5
    #diffMidLim = diffLim * 0.5

    # Delta color scale limits
    # Log normalize deltas to get summary stats
    log10variance = np.nanvar(log10deltaDiffMap)
    lognormSigma = 10**(log10variance ** 0.5)
    # clip outliers and recalculate - try this b/c at higher rez, more outliers skew variance
    log10deltaDiffMapClipped = log10deltaDiffMap
    # NOTE: Working with logged delta here, so must use log10 of lognormSigma * 3
    log10deltaDiffMapClipped[log10deltaDiffMap > np.log10(lognormSigma * 4)] = np.log10(lognormSigma * 4)
    log10varianceClipped = np.nanvar(log10deltaDiffMapClipped)
    lognormSigmaClipped = 10**(log10varianceClipped ** 0.5)
    diffLim = lognormSigmaClipped * 4 # unlogged std dev * n
    diffMidLim = lognormSigmaClipped * 1 # unlogged std dev * 1
    # deltaDiffCmap WITH shifted midband:
    midBandStart = 1 - (diffLim+diffMidLim)/(2*diffLim)
    midBandStop = (diffLim+diffMidLim)/(2*diffLim)
    deltaDiffCmap = shiftedColorMapWithMidBand(defaultDiffDivCmap,
                                              midBandStart=midBandStart,
                                              midBandStop=midBandStop,
                                              name='deltaDiffCmap')

    # Log delta color scale limits
    log10diffLim = (log10varianceClipped ** 0.5) * 4 # logged std dev * n
    log10diffMidLim = (log10varianceClipped ** 0.5) * 1 # unlogged std dev * 1
    # deltaLog10DiffCmap WITH shifted midband:
    midBandStart = 1 - (log10diffLim+log10diffMidLim)/(2*log10diffLim)
    midBandStop = (log10diffLim+log10diffMidLim)/(2*log10diffLim)
    deltaLog10DiffCmap = shiftedColorMapWithMidBand(defaultDiffDivCmap,
                                              midBandStart=midBandStart,
                                              midBandStop=midBandStop,
                                              name='deltaLog10DiffCmap')

    # minDiff = np.nanpercentile(abs(np.array(allDeltas)), 50)
    minDiff = lognormSigmaClipped * 1 # unlogged std dev * 1
    log10minDiff = (log10varianceClipped ** 0.5) * 1 # unlogged std dev

    # Log2Ratio color scale limits
    #log2RatioLim = 3
    # [nobs, minmax, mean, variance, skewness, kurtosis] = stats.describe(allLog2ratios, nan_policy='omit')
    # log2RatioLim = (variance ** 0.5) * 3 # std dev * n
    # log2RatioMidLim = (variance ** 0.5) * 1.5  # std dev * 1.5
    log2RatioLim = 3.0
    log2RatioMidLim = 1.0
    # log2RatioDiffCMap WITH shifted midband:
    midBandStart = 1 - (log2RatioLim+log2RatioMidLim)/(2*log2RatioLim)
    midBandStop = (log2RatioLim+log2RatioMidLim)/(2*log2RatioLim)
    log2RatioCmap = shiftedColorMapWithMidBand(defaultDiffDivCmap,
                                              midBandStart=midBandStart,
                                              midBandStop=midBandStop,
                                              name='log2RatioCmap')
    # print str(minSum) + "; "  + str(minDiff) + "; " + str(diffLim) + "; " + str(log2RatioLim)
    
    # 170804
    # corDeltaDiffMap WITH shifted midband:
    corDeltaLim = 1.0
    corDeltaMidLim = 0.1
    midBandStart = 1.0 - (corDeltaLim+corDeltaMidLim)/(2*corDeltaLim)
    midBandStop = (corDeltaLim+corDeltaMidLim)/(2*corDeltaLim)
    deltaCorrDiffCmap = shiftedColorMapWithMidBand(defaultDiffDivCmap,
                                              midBandStart=midBandStart,
                                              midBandStop=midBandStop,
                                              name='deltaCorrDiffCmap')

    ########################
    # Print some full maps!
    # Print all full maps first b/c rotation messed up full map plotting for some reason
    print "Full maps"
    for diffMapID in range(len(diffMaps)):
        # diffMapID = 1
        print str(diffMapID) + ": " + diffMaps[str(diffMapID)]
        [ diffMapA, diffMapB ] = diffMaps[str(diffMapID)].split('<<<--->>>')

        # The QN Data!
        dataA = dataDict_QN_NaN[diffMapA].copy()  # VERY IMPORTANT TO COPY OTHERWISE dataDict gets changed!
        #print dataA.shape
        #print np.nanmean(dataA)

        dataB = dataDict_QN_NaN[diffMapB].copy()  # VERY IMPORTANT TO COPY OTHERWISE dataDict gets changed!
        #print dataB.shape
        #print np.nanmean(dataB)

        sumMap = dataA + dataB

        deltaDiffMap = dataA - dataB
        #print deltaDiffMap.shape
        #print np.nanmean(deltaDiffMap)

        # log delta
        log10deltaDiffMap = np.log10(np.abs(deltaDiffMap))
        # Remove fractional deltas
        log10deltaDiffMap[np.abs(deltaDiffMap)<1] = np.nan
        # Return sign to logged data
        log10deltaDiffMap[deltaDiffMap<0] = log10deltaDiffMap[deltaDiffMap<0] * -1

        # Log2 ratio
        log2ratioDiffMap = np.log2( dataA / dataB )
        #print log2ratioDiffMap.shape
        #print np.nanmean(log2ratioDiffMap)

        log2ratioDiffMapAdj = log2ratioDiffMap.copy()
        # smallChangeIdx = np.where(abs(deltaDiffMap) < minDiff)
        # for sci in range(len(smallChangeIdx[0])):
        #     log2ratioDiffMapAdj[smallChangeIdx[0][sci], smallChangeIdx[1][sci]] = np.nan
        smallCountsIdx = np.where(sumMap < minSum)
        for sci in range(len(smallCountsIdx[0])):
            log2ratioDiffMapAdj[smallCountsIdx[0][sci], smallCountsIdx[1][sci]] = np.nan

        log2ratioDiffMapAdjQuantiles = {}
        for percentile in percentiles:
            log2ratioDiffMapAdjQuantiles[str(percentile)] = log2ratioDiffMap.copy()
            smallCountsIdx = np.where(sumMap < np.nanpercentile(sumMap, percentile))
            for sci in range(len(smallCountsIdx[0])):
                log2ratioDiffMapAdjQuantiles[str(percentile)][smallCountsIdx[0][sci], smallCountsIdx[1][sci]] = np.nan

        ################################################
        # 170623
        # For differential correlation matrices
        # The cor QN Data!
        corDataA = corDataDict_QN_NaN[diffMapA].copy()  # VERY IMPORTANT TO COPY OTHERWISE dataDict gets changed!
        #print dataA.shape
        #print np.nanmean(dataA)
        corDataB = corDataDict_QN_NaN[diffMapB].copy()  # VERY IMPORTANT TO COPY OTHERWISE dataDict gets changed!
        #print dataB.shape
        #print np.nanmean(dataB)

        corDeltaDiffMap = corDataA - corDataB
        #print deltaDiffMap.shape
        #print np.nanmean(deltaDiffMap)

        corrDiffLim = 1


        ################################################
        # Deltas
        outDir = diffDir  + "/delta"
        if not os.path.exists(outDir):
            os.makedirs(outDir)

        outputLabel =  "delta." + "-V-".join(diffMaps[str(diffMapID)].split('<<<--->>>')) + "." + theRez + "." + chrom
        title = outputLabel.replace(".", " ")
        outFilenameBase = outDir + "/" + outputLabel
        # plot_divergent_maps(np.clip(deltaDiffMap, a_min=-diffLim, a_max=diffLim), intTheRez, title, outFilenameBase)
        # Default color scheme:
        # plot_matplotlib_heatmap(deltaDiffMap, intTheRez, title, outFilenameBase, cmap=defaultDiffDivCmap, vmax=diffLim, vmin=-diffLim)
        plot_matplotlib_heatmap(deltaDiffMap, intTheRez, title, outFilenameBase,
                                cmap=deltaDiffCmap, vmax=diffLim, vmin=-diffLim)

        ################################################
        # Log10 Deltas
        outDir = diffDir  + "/delta.log10"
        if not os.path.exists(outDir):
            os.makedirs(outDir)

        outputLabel =  "delta.log10." + "-V-".join(diffMaps[str(diffMapID)].split('<<<--->>>')) + "." + theRez + "." + chrom
        title = outputLabel.replace(".", " ")
        outFilenameBase = outDir + "/" + outputLabel
        plot_matplotlib_heatmap(log10deltaDiffMap, intTheRez, title, outFilenameBase,
                                cmap=deltaLog10DiffCmap, vmax=log10diffLim, vmin=-log10diffLim)

        ################################################
        # Log2Ratios
        outDir = diffDir + "/log2ratio"
        if not os.path.exists(outDir):
            os.makedirs(outDir)

        outputLabel =  "log2ratio." + "-V-".join(diffMaps[str(diffMapID)].split('<<<--->>>')) + "." + theRez + "." + chrom
        title = outputLabel.replace(".", " ")
        outFilenameBase = outDir + "/" + outputLabel
        # plot_divergent_maps(np.clip(log2ratioDiffMap, a_min=-log2RatioLim, a_max=log2RatioLim), intTheRez, title, outFilenameBase)
        plot_matplotlib_heatmap(log2ratioDiffMap, intTheRez, title, outFilenameBase,
                                cmap=log2RatioCmap, vmax=log2RatioLim, vmin=-log2RatioLim)

        ################################################
        # Filtered Log2Ratios
        outDir = diffDir  + "/filteredLog2ratio"
        if not os.path.exists(outDir):
            os.makedirs(outDir)

        outputLabel =  "log2ratio.filtered." + "-V-".join(diffMaps[str(diffMapID)].split('<<<--->>>')) + "." + theRez + "." + chrom
        title = outputLabel.replace(".", " ")
        outFilenameBase = outDir + "/" + outputLabel
        #plot_divergent_maps(np.clip(log2ratioDiffMapAdj, a_min=-log2RatioLim, a_max=log2RatioLim), intTheRez, title, outFilenameBase)
        plot_matplotlib_heatmap(log2ratioDiffMapAdj, intTheRez, title, outFilenameBase,
                                cmap=log2RatioCmap, vmax=log2RatioLim, vmin=-log2RatioLim)

        ################################################
        # Quantile filtered Log2Ratios
        for whichQuantile in log2ratioDiffMapAdjQuantiles.keys():
            outDir = diffDir  + "/filteredLog2ratio.q" + whichQuantile
            if not os.path.exists(outDir):
                os.makedirs(outDir)

            outputLabel =  "log2ratio.filtered.q" + whichQuantile + "." + "-V-".join(diffMaps[str(diffMapID)].split('<<<--->>>')) + "." + theRez + "." + chrom
            title = outputLabel.replace(".", " ")
            outFilenameBase = outDir + "/" + outputLabel
            #plot_divergent_maps(np.clip(log2ratioDiffMapAdj, a_min=-log2RatioLim, a_max=log2RatioLim), intTheRez, title, outFilenameBase)
            plot_matplotlib_heatmap(log2ratioDiffMapAdjQuantiles[whichQuantile], intTheRez, title, outFilenameBase,
                                    cmap=log2RatioCmap, vmax=log2RatioLim, vmin=-log2RatioLim)

        ################################################
        # 170623
        # Differential correlation matrices
        outDir = diffDir  + "/corDelta"
        if not os.path.exists(outDir):
            os.makedirs(outDir)

        outputLabel =  "corDelta." + "-V-".join(diffMaps[str(diffMapID)].split('<<<--->>>')) + "." + theRez + "." + chrom
        title = outputLabel.replace(".", " ")
        outFilenameBase = outDir + "/" + outputLabel
        plot_matplotlib_heatmap(corDeltaDiffMap, intTheRez, title, outFilenameBase,
                                cmap=defaultDiffDivCmap, vmax=corrDiffLim, vmin=-corrDiffLim)

        # 170804
        # Differential correlation matrices w/ wider white band
        outputLabel =  "corDelta.wideMidBand." + "-V-".join(diffMaps[str(diffMapID)].split('<<<--->>>')) + "." + theRez + "." + chrom
        title = outputLabel.replace(".", " ")
        outFilenameBase = outDir + "/" + outputLabel
        plot_matplotlib_heatmap(corDeltaDiffMap, intTheRez, title, outFilenameBase,
                                cmap=deltaCorrDiffCmap, vmax=corrDiffLim, vmin=-corrDiffLim)

    ########################
    #  Print some rotated and LOI differential maps!
    print "Print rot maps"
    for diffMapID in range(len(diffMaps)):
        print str(diffMapID) + ": " + diffMaps[str(diffMapID)]
        [ diffMapA, diffMapB ] = diffMaps[str(diffMapID)].split('<<<--->>>')

        # The QN Data!
        dataA = dataDict_QN_NaN[diffMapA].copy()  # VERY IMPORTANT TO COPY OTHERWISE dataDict gets changed!
        #print dataA.shape
        #print np.nanmean(dataA)

        dataB = dataDict_QN_NaN[diffMapB].copy()  # VERY IMPORTANT TO COPY OTHERWISE dataDict gets changed!
        #print dataB.shape
        #print np.nanmean(dataB)

        sumMap = dataA + dataB

        deltaDiffMap = dataA - dataB
        #print deltaDiffMap.shape
        #print np.nanmean(deltaDiffMap)

        # log delta
        log10deltaDiffMap = np.log10(np.abs(deltaDiffMap))
        # Remove fractional deltas
        log10deltaDiffMap[np.abs(deltaDiffMap)<1] = np.nan
        # Return sign to logged data
        log10deltaDiffMap[deltaDiffMap<0] = log10deltaDiffMap[deltaDiffMap<0] * -1

        # Log2 ratio
        log2ratioDiffMap = np.log2( dataA / dataB )
        #print log2ratioDiffMap.shape
        #print np.nanmean(log2ratioDiffMap)

        log2ratioDiffMapAdj = log2ratioDiffMap.copy()
        # smallChangeIdx = np.where(abs(deltaDiffMap) < minDiff)
        # for sci in range(len(smallChangeIdx[0])):
        #     log2ratioDiffMapAdj[smallChangeIdx[0][sci], smallChangeIdx[1][sci]] = np.nan
        smallCountsIdx = np.where(sumMap < minSum)
        for sci in range(len(smallCountsIdx[0])):
            log2ratioDiffMapAdj[smallCountsIdx[0][sci], smallCountsIdx[1][sci]] = np.nan

        ################################################
        # 170623
        # For differential correlation matrices
        # The cor QN Data!
        corDataA = corDataDict_QN_NaN[diffMapA].copy()  # VERY IMPORTANT TO COPY OTHERWISE dataDict gets changed!
        #print dataA.shape
        #print np.nanmean(dataA)
        corDataB = corDataDict_QN_NaN[diffMapB].copy()  # VERY IMPORTANT TO COPY OTHERWISE dataDict gets changed!
        #print dataB.shape
        #print np.nanmean(dataB)

        corDeltaDiffMap = corDataA - corDataB
        #print deltaDiffMap.shape
        #print np.nanmean(deltaDiffMap)

        corrDiffLim = 1

        if intTheRez >= 100000:

            ################################################
            # Deltas
            outDir = diffDir  + "/delta.rot45"
            if not os.path.exists(outDir):
                os.makedirs(outDir)

            outputLabel =  "delta." + "-V-".join(diffMaps[str(diffMapID)].split('<<<--->>>')) + "." + theRez + "." + chrom
            title = outputLabel.replace(".", " ")
            outFilenameBase = outDir + "/" + outputLabel
            # plot_matplotlib_heatmap_r45(np.clip(deltaDiffMap, a_min=-diffLim, a_max=diffLim),
            #                                     theRez, title, outFilenameBase, cmap=defaultDiffDivCmap)
            plot_matplotlib_heatmap_r45(deltaDiffMap, intTheRez, title, outFilenameBase,
                                        cmap=deltaDiffCmap, vmax=diffLim, vmin=-diffLim)

            ################################################
            # Log10 Deltas
            outDir = diffDir  + "/delta.log10.rot45"
            if not os.path.exists(outDir):
                os.makedirs(outDir)

            outputLabel =  "delta.log10." + "-V-".join(diffMaps[str(diffMapID)].split('<<<--->>>')) + "." + theRez + "." + chrom
            title = outputLabel.replace(".", " ")
            outFilenameBase = outDir + "/" + outputLabel
            plot_matplotlib_heatmap_r45(log10deltaDiffMap, intTheRez, title, outFilenameBase,
                                        cmap=deltaLog10DiffCmap, vmax=log10diffLim, vmin=-log10diffLim)

            ################################################
            # Log2Ratios
            outDir = diffDir  + "/log2ratio.rot45"
            if not os.path.exists(outDir):
                os.makedirs(outDir)

            outputLabel =  "log2ratio." + "-V-".join(diffMaps[str(diffMapID)].split('<<<--->>>')) + "." + theRez + "." + chrom
            title = outputLabel.replace(".", " ")
            outFilenameBase = outDir + "/" + outputLabel
            # plot_matplotlib_heatmap_r45(np.clip(log2ratioDiffMap, a_min=-log2RatioLim, a_max=log2RatioLim),
            #                                     theRez, title, outFilenameBase,  cmap=defaultDiffDivCmap)
            plot_matplotlib_heatmap_r45(log2ratioDiffMap, intTheRez, title, outFilenameBase,
                                        cmap=log2RatioCmap, vmax=log2RatioLim, vmin=-log2RatioLim)

            ################################################
            # Filtered Log2Ratios
            outDir = diffDir  + "/filteredLog2ratio.rot45"
            if not os.path.exists(outDir):
                os.makedirs(outDir)

            outputLabel =  "log2ratio.filtered." + "-V-".join(diffMaps[str(diffMapID)].split('<<<--->>>')) + "." + theRez + "." + chrom
            title = outputLabel.replace(".", " ")
            outFilenameBase = outDir + "/" + outputLabel
            # plot_matplotlib_heatmap_r45(np.clip(log2ratioDiffMapAdj, a_min=-log2RatioLim, a_max=log2RatioLim),
            #                                     theRez, title, outFilenameBase, cmap=defaultDiffDivCmap)
            plot_matplotlib_heatmap_r45(log2ratioDiffMapAdj, intTheRez, title, outFilenameBase,
                                        cmap=log2RatioCmap, vmax=log2RatioLim, vmin=-log2RatioLim)

            ################################################
            # 170623
            # Differential correlation matrices
            outDir = diffDir  + "/corDelta.rot45"
            if not os.path.exists(outDir):
                os.makedirs(outDir)

            outputLabel =  "corDelta." + "-V-".join(diffMaps[str(diffMapID)].split('<<<--->>>')) + "." + theRez + "." + chrom
            title = outputLabel.replace(".", " ")
            outFilenameBase = outDir + "/" + outputLabel
            plot_matplotlib_heatmap_r45(corDeltaDiffMap, intTheRez, title, outFilenameBase,
                                   cmap=defaultDiffDivCmap, vmax=corrDiffLim, vmin=-corrDiffLim)

            # 170804
            # Differential correlation matrices w/ wider white band
            outputLabel =  "corDelta.wideMidBand." + "-V-".join(diffMaps[str(diffMapID)].split('<<<--->>>')) + "." + theRez + "." + chrom
            title = outputLabel.replace(".", " ")
            outFilenameBase = outDir + "/" + outputLabel
            plot_matplotlib_heatmap_r45(corDeltaDiffMap, intTheRez, title, outFilenameBase,
                                   cmap=deltaCorrDiffCmap, vmax=corrDiffLim, vmin=-corrDiffLim)
    

        ################################################
        # LOI differential maps!

        if len(LOI) >= 0:
            ################################################
            # Deltas
            outDir = diffDir  + "/delta.LOIzoom100Rot45"
            if not os.path.exists(outDir):
                os.makedirs(outDir)

            outputLabel =  "delta." + "-V-".join(diffMaps[str(diffMapID)].split('<<<--->>>')) + "." + theRez + "." + chrom
            title = outputLabel.replace(".", " ")
            outFilenameBase = outDir + "/" + outputLabel
            # plot_matplotlib_heatmap_with_LOI_r45_zoom(np.clip(deltaDiffMap, a_min=-diffLim, a_max=diffLim),
            #                                     LOI, theRez, intTheRez*100, title, outFilenameBase, cmap=defaultDiffDivCmap)
            plot_matplotlib_heatmap_with_LOI_r45_zoom(deltaDiffMap, LOI, intTheRez, intTheRez*100, title, outFilenameBase,
                                                      cmap=deltaDiffCmap, vmax=diffLim, vmin=-diffLim)

            ################################################
            # Log10 Deltas
            outDir = diffDir  + "/delta.log10.LOIzoom100Rot45"
            if not os.path.exists(outDir):
                os.makedirs(outDir)

            outputLabel =  "delta.log10." + "-V-".join(diffMaps[str(diffMapID)].split('<<<--->>>')) + "." + theRez + "." + chrom
            title = outputLabel.replace(".", " ")
            outFilenameBase = outDir + "/" + outputLabel
            plot_matplotlib_heatmap_with_LOI_r45_zoom(log10deltaDiffMap, LOI, intTheRez, intTheRez*100, title, outFilenameBase,
                                        cmap=deltaLog10DiffCmap, vmax=log10diffLim, vmin=-log10diffLim)

            ################################################
            # Log2Ratios
            outDir = diffDir  + "/log2ratio.LOIzoom100Rot45"
            if not os.path.exists(outDir):
                os.makedirs(outDir)

            outputLabel =  "log2ratio." + "-V-".join(diffMaps[str(diffMapID)].split('<<<--->>>')) + "." + theRez + "." + chrom
            title = outputLabel.replace(".", " ")
            outFilenameBase = outDir + "/" + outputLabel
            # plot_matplotlib_heatmap_with_LOI_r45_zoom(np.clip(log2ratioDiffMap, a_min=-log2RatioLim, a_max=log2RatioLim),
            #                                     LOI, theRez, intTheRez*100, title, outFilenameBase, cmap=defaultDiffDivCmap)
            plot_matplotlib_heatmap_with_LOI_r45_zoom(log2ratioDiffMap, LOI, intTheRez, intTheRez*100, title, outFilenameBase,
                                                      cmap=log2RatioCmap, vmax=log2RatioLim, vmin=-log2RatioLim)

            ################################################
            # Filtered Log2Ratios
            outDir = diffDir  + "/filteredLog2ratio.LOIzoom100Rot45"
            if not os.path.exists(outDir):
                os.makedirs(outDir)

            outputLabel =  "log2ratio.filtered." + "-V-".join(diffMaps[str(diffMapID)].split('<<<--->>>')) + "." + theRez + "." + chrom
            title = outputLabel.replace(".", " ")
            outFilenameBase = outDir + "/" + outputLabel
            # plot_matplotlib_heatmap_with_LOI_r45_zoom(np.clip(log2ratioDiffMapAdj, a_min=-log2RatioLim, a_max=log2RatioLim),
            #                                     LOI, theRez, intTheRez*100, title, outFilenameBase, cmap=defaultDiffDivCmap)
            plot_matplotlib_heatmap_with_LOI_r45_zoom(log2ratioDiffMapAdj, LOI, intTheRez, intTheRez*100, title, outFilenameBase,
                                                      cmap=log2RatioCmap, vmax=log2RatioLim, vmin=-log2RatioLim)

            ################################################
            # 170623
            # Differential correlation matrices
            outDir = diffDir  + "/corDelta.LOIzoom100Rot45"
            if not os.path.exists(outDir):
                os.makedirs(outDir)

            outputLabel =  "corDelta." + "-V-".join(diffMaps[str(diffMapID)].split('<<<--->>>')) + "." + theRez + "." + chrom
            title = outputLabel.replace(".", " ")
            outFilenameBase = outDir + "/" + outputLabel
            plot_matplotlib_heatmap_with_LOI_r45_zoom(corDeltaDiffMap, LOI, intTheRez, intTheRez*100, title, outFilenameBase,
                                   cmap=defaultDiffDivCmap, vmax=corrDiffLim, vmin=-corrDiffLim)

            # 170804
            # Differential correlation matrices w/ wider white band
            outputLabel =  "corDelta.wideMidBand." + "-V-".join(diffMaps[str(diffMapID)].split('<<<--->>>')) + "." + theRez + "." + chrom
            title = outputLabel.replace(".", " ")
            outFilenameBase = outDir + "/" + outputLabel
            plot_matplotlib_heatmap_with_LOI_r45_zoom(corDeltaDiffMap, LOI, intTheRez, intTheRez*100, title, outFilenameBase,
                                    cmap=deltaCorrDiffCmap, vmax=corrDiffLim, vmin=-corrDiffLim)



################################################
################################################
################################################
# Virtual 4C plots for LOI

if 'v4cplots' in toRun and toRun['v4cplots']:

    if len(LOI) > 0:
        print "\n\n\nvirtual 4C line plots"

        # figure width
        # for theRez in (1000000, 500000, 100000, 40000):
        #    print 20 + (10 * (500000/intTheRez))
        figureWidth = 20 + (10 * (500000/intTheRez))

        # Use TrueType fonts that work in illustrator:
        #    http://stackoverflow.com/questions/5956182/cannot-edit-text-in-chart-exported-by-matplotlib-and-opened-in-illustrator
        matplotlib.rcParams['pdf.fonttype'] = 42

        v4CDir = "output/" + contactMapOutputDesc + "/" + outputDir + "/heatmaps_QN/virtual4C"
        if not os.path.exists(v4CDir):
            os.makedirs(v4CDir)

        vlineLoi = ['dxz4', 'firre']

        yClipPct = 0.50

        # Get max Y-vals
        maxYs = {}
        for aloi in LOI.keys():
            maxYs[aloi] = 0

        for whatData in dataSets:
            # The QN Data!
            data = dataDict_QN_NaN[whatData].copy()  # VERY IMPORTANT TO COPY OTHERWISE dataDict gets changed!

            for aloi in LOI.keys():
                # aloi = 'xist'
                #print aloi
                #print LOI[aloi]

                chrCoords = LOI[aloi].split(":")
                theChr = chrCoords[0]
                if theChr == thechrom:
                    #print "\tPlot it!"
                    startStop = chrCoords[1].split("-")
                    startStop = map(int, startStop)
                    #startStop = [round(x / float(theRez)) for x in startStop]
                    geneMidPoint = (startStop[0] + (startStop[1] - startStop[0]))
                    geneMidPointBin = geneMidPoint / intTheRez

                    # get data row
                    yVals = data[geneMidPointBin,]
                    maxY = np.nanmax(yVals)
                    if maxY > maxYs[aloi]:
                        maxYs[aloi] = maxY

        # Plot graphs
        for whatData in dataSets:
            # whatData = 'Brain_Xa'

            inputDir = contactMapIOmetadata[whatData]["inputDir"]
            inputFile = contactMapIOmetadata[whatData]["inputFile"]
            outputDir = contactMapIOmetadata[whatData]["outputDir"]
            outputLabel = contactMapIOmetadata[whatData]["outputLabel"]

            #print whatData + ": " + inputDir + ", " + inputFile + ", " + outputDir + ", " + outputLabel

            # The QN Data!
            data = dataDict_QN_NaN[whatData].copy()  # VERY IMPORTANT TO COPY OTHERWISE dataDict gets changed!
            #print data.shape
            #print np.nanmean(data)

            for aloi in LOI.keys():
                # aLOI = 'xist'
                #print aloi
                #print LOI[aloi]

                chrCoords = LOI[aloi].split(":")
                theChr = chrCoords[0]
                if theChr == thechrom:
                    #print "\tPlot it!"
                    startStop = chrCoords[1].split("-")
                    startStop = map(int, startStop)
                    #startStop = [round(x / float(theRez)) for x in startStop]
                    geneMidPoint = (startStop[0] + (startStop[1] - startStop[0]))
                    geneMidPointBin = geneMidPoint / intTheRez

                    #print(aloi + " midpoint: " + str(geneMidPoint) + "; Bin: " + str(geneMidPointBin))

                    # get data row
                    yVals = data[geneMidPointBin,]
                    # maxY = np.nanmax(yVals)

                    # Insulation score loci
                    xVals = np.r_[0:len(yVals)]
                    xVals = (xVals * intTheRez) / (1.0*10**6)

                    # X-axis limits
                    minX = 0
                    maxX = int(np.ceil(float(chromSizes[thechrom]) / 1000000.0))

                    fig, ax = plt.subplots(1, 1, figsize=(figureWidth, 10))
                    #plots = {}

                    # plots[whatData], = ax.plot(xVals, yVals,
                    #                           '-', label=whatData,
                    #                           color=contactMapIOmetadata[whatData]['dataColor'],
                    #                           linewidth=2,
                    #                           linestyle='solid')
                    ax.vlines(xVals, [0], yVals,
                                              label=whatData,
                                              colors=contactMapIOmetadata[whatData]['dataColor'] )

                    ax.grid(True)
                    ax.set_ylabel('Contact counts\n(clipped at ' + str(int(round(yClipPct * 100, 0))) + '% of max)', fontsize=14)
                    ax.set_xlabel('Mb\n(' + theRez + 'bp bins)', fontsize=14)
                    for tick in ax.xaxis.get_major_ticks():
                        tick.label.set_fontsize(12)
                    for tick in ax.yaxis.get_major_ticks():
                        tick.label.set_fontsize(12)

                    ax.set_xlim(minX, maxX)
                    ax.set_ylim(0, yClipPct*maxYs[aloi])

                    #ax.legend(loc='best', fontsize=14)
                    ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                           ncol=2, mode="expand", borderaxespad=0., fontsize=14)

                    ax.axvline(xVals[geneMidPointBin], color='k', linestyle='--', linewidth=0.5)
                    # http://matplotlib.org/examples/pylab_examples/fonts_demo.html
                    font0 = FontProperties()
                    font1 = font0.copy()
                    font1.set_size('16')
                    ax.text(xVals[geneMidPointBin], yClipPct*maxYs[aloi], aloi, rotation=90, verticalalignment="top", fontproperties=font1)

                    for vloi in vlineLoi:
                        if vloi in LOI.keys() and vloi != aloi:
                            v_chrCoords = LOI[vloi].split(":")
                            v_theChr = v_chrCoords[0]
                            if v_theChr == thechrom:
                                #print "\tPlot it!"
                                v_startStop = v_chrCoords[1].split("-")
                                v_startStop = map(int, v_startStop)
                                v_geneMidPoint = (v_startStop[0] + (v_startStop[1] - v_startStop[0]))
                                v_geneMidPointBin = v_geneMidPoint / intTheRez

                                ax.axvline(xVals[v_geneMidPointBin], color='k', linestyle='--', linewidth=0.5)
                                # http://matplotlib.org/examples/pylab_examples/fonts_demo.html
                                font0 = FontProperties()
                                font1 = font0.copy()
                                font1.set_size('12')
                                ax.text(xVals[v_geneMidPointBin], yClipPct*maxYs[aloi], vloi, rotation=90, verticalalignment="top", fontproperties=font1)

                    # Tidy up and save
                    #fig.subplots_adjust(top=1)
                    #fig.tight_layout()
                    outputLabel =  "virtual4Cplot." + whatData +"." + aloi + "." + theRez + "." + chrom
                    title = outputLabel.replace(".", " ")
                    fig.suptitle(title, fontsize=20, y=1.01)
                    fig.savefig(v4CDir + "/" + outputLabel + ".png", dpi=150, bbox_inches='tight')
                    fig.savefig(v4CDir + "/" + outputLabel + ".pdf", dpi=150, bbox_inches='tight')
                    plt.close()


################################################
################################################
################################################
# Virtual 4C plots for all loi

if 'v4cplots4allloci' in toRun and toRun['v4cplots4allloci']:

    print "\n\n\nvirtual 4C line plots for all loci"

    # figure width
    # for theRez in (1000000, 500000, 100000, 40000):
    #    print 20 + (10 * (500000/intTheRez))
    figureWidth = 20 + (10 * (500000/intTheRez))

    # Use TrueType fonts that work in illustrator:
    #    http://stackoverflow.com/questions/5956182/cannot-edit-text-in-chart-exported-by-matplotlib-and-opened-in-illustrator
    matplotlib.rcParams['pdf.fonttype'] = 42

    v4CDir = "output/" + contactMapOutputDesc + "/" + outputDir + "/heatmaps_QN/virtual4C.allLoci"
    if not os.path.exists(v4CDir):
        os.makedirs(v4CDir)

    vlineLoi = ['dxz4', 'firre']

    # Get max Y-vals
    yClipPct = 0.20
    maxMaxY = 0
    for whatData in dataSets:
        # The QN Data!
        data = dataDict_QN_NaN[whatData].copy()  # VERY IMPORTANT TO COPY OTHERWISE dataDict gets changed!
        maxY = np.nanmax(data)
        if maxY > maxMaxY:
            maxMaxY = maxY

    # Plot graphs
    for whatData in dataSets:
        # whatData = 'Brain_Xa'

        inputDir = contactMapIOmetadata[whatData]["inputDir"]
        inputFile = contactMapIOmetadata[whatData]["inputFile"]
        outputDir = contactMapIOmetadata[whatData]["outputDir"]
        outputLabel = contactMapIOmetadata[whatData]["outputLabel"]

        # print whatData + ": " + inputDir + ", " + inputFile + ", " + outputDir + ", " + outputLabel

        # The QN Data!
        data = dataDict_QN_NaN[whatData].copy()  # VERY IMPORTANT TO COPY OTHERWISE dataDict gets changed!
        #print data.shape
        #print np.nanmean(data)

        for currBin in range(data.shape[0]):

            # get data row
            yVals = data[currBin,]
            # maxY = np.nanmax(yVals)

            if all(np.isnan(yVals)):
                continue

            # Insulation score loci
            xVals = np.r_[0:len(yVals)]
            xVals = (xVals * intTheRez) / (1.0*10**6)

            # X-axis limits
            minX = 0
            maxX = int(np.ceil(float(chromSizes[thechrom]) / 1000000.0))

            fig, ax = plt.subplots(1, 1, figsize=(figureWidth, 10))
            #plots = {}

            # plots[whatData], = ax.plot(xVals, yVals,
            #                           '-', label=whatData,
            #                           color=contactMapIOmetadata[whatData]['dataColor'],
            #                           linewidth=2,
            #                           linestyle='solid')
            ax.vlines(xVals, [0], yVals,
                                      label=whatData,
                                      colors=contactMapIOmetadata[whatData]['dataColor'] )

            ax.grid(True)
            ax.set_ylabel('Contact counts\n(clipped at ' + str(int(round(yClipPct * 100, 0))) + '% of max)', fontsize=14)
            ax.set_xlabel('Mb\n(' + theRez + 'bp bins)', fontsize=14)
            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(12)
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(12)

            ax.set_xlim(minX, maxX)
            ax.set_ylim(0, yClipPct*maxMaxY)

            #ax.legend(loc='best', fontsize=14)
            ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                   ncol=2, mode="expand", borderaxespad=0., fontsize=14)

            ax.axvline(xVals[currBin], color='k', linestyle='--', linewidth=0.5)
            # http://matplotlib.org/examples/pylab_examples/fonts_demo.html
            font0 = FontProperties()
            font1 = font0.copy()
            font1.set_size('16')
            binMidpointLocus = (currBin * intTheRez) / (1.0*10**6)
            ax.text(xVals[currBin], yClipPct*maxMaxY, str(binMidpointLocus) + " Mb", rotation=90, verticalalignment="top", fontproperties=font1)

            for vloi in vlineLoi:
                if vloi in LOI.keys():
                    v_chrCoords = LOI[vloi].split(":")
                    v_theChr = v_chrCoords[0]
                    if v_theChr == thechrom:
                        #print "\tPlot it!"
                        v_startStop = v_chrCoords[1].split("-")
                        v_startStop = map(int, v_startStop)
                        v_geneMidPoint = (v_startStop[0] + (v_startStop[1] - v_startStop[0]))
                        v_currBin = v_geneMidPoint / intTheRez

                        ax.axvline(xVals[v_currBin], color='k', linestyle='--', linewidth=0.5)
                        # http://matplotlib.org/examples/pylab_examples/fonts_demo.html
                        font0 = FontProperties()
                        font1 = font0.copy()
                        font1.set_size('12')
                        ax.text(xVals[v_currBin], yClipPct*maxMaxY, vloi, rotation=90, verticalalignment="top", fontproperties=font1)

            # Tidy up and save
            #fig.subplots_adjust(top=1)
            #fig.tight_layout()
            outputLabel =  "virtual4Cplot." + whatData + ".Bin" + str(currBin).zfill(4) + "." + theRez + "." + chrom
            title = outputLabel.replace(".", " ")
            fig.suptitle(title, fontsize=20, y=1.01)
            fig.savefig(v4CDir + "/" + outputLabel + ".png", dpi=150, bbox_inches='tight')
            fig.savefig(v4CDir + "/" + outputLabel + ".pdf", dpi=150, bbox_inches='tight')
            plt.close()




################################################
################################################
################################################
# Virtual 4C plots for all loi - overlayed line plots

if 'v4cplots4alllocilines' in toRun and toRun['v4cplots4alllocilines']:

    print "\n\n\nvirtual 4C line plots for all loci - overlayed line plots"

    noBrainsDataSets = filter(lambda x:'Brain' not in x, dataSets)

    # figure width
    # for theRez in (1000000, 500000, 100000, 40000):
    #    print 20 + (10 * (500000/intTheRez))
    figureWidth = 20 + (10 * (500000/intTheRez))

    # Use TrueType fonts that work in illustrator:
    #    http://stackoverflow.com/questions/5956182/cannot-edit-text-in-chart-exported-by-matplotlib-and-opened-in-illustrator
    matplotlib.rcParams['pdf.fonttype'] = 42

    v4CDir = "output/" + contactMapOutputDesc + "/" + outputDir + "/heatmaps_QN/virtual4C.allLoci.lineplots"
    if not os.path.exists(v4CDir):
        os.makedirs(v4CDir)

    vlineLoi = ['dxz4', 'firre']

    # Get max Y-vals
    yClipPct = 0.20
    maxMaxY = 0
    binRange = 0
    for whatData in noBrainsDataSets:
        # The QN Data!
        data = dataDict_QN_NaN[whatData].copy()  # VERY IMPORTANT TO COPY OTHERWISE dataDict gets changed!
        binLen = data.shape[0]
        maxY = np.nanmax(data)
        if maxY > maxMaxY:
            maxMaxY = maxY

    xVals = np.r_[0:binLen]
    xVals = (xVals * intTheRez) / (1.0*10**6)

    for currBin in range(binLen):

        if unMappleRegionBool[currBin]:
            continue

        fig, ax = plt.subplots(1, 1, figsize=(figureWidth, 10))
        #plots = {}

        # Plot graphs
        for whatData in noBrainsDataSets:

            inputDir = contactMapIOmetadata[whatData]["inputDir"]
            inputFile = contactMapIOmetadata[whatData]["inputFile"]
            outputDir = contactMapIOmetadata[whatData]["outputDir"]
            outputLabel = contactMapIOmetadata[whatData]["outputLabel"]

            # print whatData + ": " + inputDir + ", " + inputFile + ", " + outputDir + ", " + outputLabel

            # The QN Data!
            data = dataDict_QN_NaN[whatData].copy()  # VERY IMPORTANT TO COPY OTHERWISE dataDict gets changed!
            #print data.shape
            #print np.nanmean(data)

            # get data row
            yVals = data[currBin,]
            # maxY = np.nanmax(yVals)

            # http://stackoverflow.com/questions/20618804/how-to-smooth-a-curve-in-the-right-way
            # http://scipy.github.io/devdocs/generated/scipy.signal.savgol_filter.html#scipy.signal.savgol_filter
            # window size, polynomial order
            yVals = signal.savgol_filter(yVals, 5, 3)

            # if all(np.isnan(yVals)):
            #     continue

            # # Insulation score loci
            # xVals = np.r_[0:len(yVals)]
            # xVals = (xVals * intTheRez) / (1.0*10**6)

            # plots[whatData], = ax.plot(xVals, yVals,
            #                           '-', label=whatData,
            #                           color=contactMapIOmetadata[whatData]['dataColor'],
            #                           linewidth=2,
            #                           linestyle='solid')
            ax.plot(xVals, yVals, '-', label=whatData,
                                    linewidth=contactMapIOmetadata[whatData]['dataLineWidth'],
                                    color=contactMapIOmetadata[whatData]['dataColor'])

        # X-axis limits
        minX = 0
        maxX = int(np.ceil(float(chromSizes[thechrom]) / 1000000.0))

        ax.grid(True)
        ax.set_ylabel('Contact counts\n(clipped at ' + str(int(round(yClipPct * 100, 0))) + '% of max)', fontsize=14)
        ax.set_xlabel('Mb\n(' + theRez + 'bp bins)', fontsize=14)
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(12)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(12)

        ax.set_xlim(minX, maxX)
        ax.set_ylim(0, yClipPct*maxMaxY)

        #ax.legend(loc='best', fontsize=14)
        ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=2, mode="expand", borderaxespad=0., fontsize=14)

        ax.axvline(xVals[currBin], color='k', linestyle='--', linewidth=0.5)
        # http://matplotlib.org/examples/pylab_examples/fonts_demo.html
        font0 = FontProperties()
        font1 = font0.copy()
        font1.set_size('16')
        binMidpointLocus = (currBin * intTheRez) / (1.0*10**6)
        ax.text(xVals[currBin], yClipPct*maxMaxY, str(binMidpointLocus) + " Mb", rotation=90, verticalalignment="top", fontproperties=font1)

        for vloi in vlineLoi:
            if vloi in LOI.keys():
                v_chrCoords = LOI[vloi].split(":")
                v_theChr = v_chrCoords[0]
                if v_theChr == thechrom:
                    #print "\tPlot it!"
                    v_startStop = v_chrCoords[1].split("-")
                    v_startStop = map(int, v_startStop)
                    v_geneMidPoint = (v_startStop[0] + (v_startStop[1] - v_startStop[0]))
                    v_currBin = v_geneMidPoint / intTheRez

                    ax.axvline(xVals[v_currBin], color='k', linestyle='--', linewidth=0.5)
                    # http://matplotlib.org/examples/pylab_examples/fonts_demo.html
                    font0 = FontProperties()
                    font1 = font0.copy()
                    font1.set_size('12')
                    ax.text(xVals[v_currBin], yClipPct*maxMaxY, vloi, rotation=90, verticalalignment="top", fontproperties=font1)

        # Tidy up and save
        #fig.subplots_adjust(top=1)
        #fig.tight_layout()
        outputLabel =  "virtual4Cplot.AllSamples.Bin" + str(currBin).zfill(4) + "." + theRez + "." + chrom
        title = outputLabel.replace(".", " ")
        fig.suptitle(title, fontsize=20, y=1.01)
        fig.savefig(v4CDir + "/" + outputLabel + ".png", dpi=150, bbox_inches='tight')
        fig.savefig(v4CDir + "/" + outputLabel + ".pdf", dpi=150, bbox_inches='tight')
        plt.close()


################################################
################################################
################################################
# Virtual 4C plots for all loi - overlayed line plots

if 'v4cplots4alllocilinesxaxi' in toRun and toRun['v4cplots4alllocilinesxaxi']:

    print "\n\n\nvirtual 4C line plots for all loci - overlayed line plots - seperate Xa and Xi"

    noBrainsDataSets = filter(lambda x:'Brain' not in x, dataSets)

    for Xtype in ['Xa', 'Xi']:
        #print(Xtype)
        XtypeString = '_' + Xtype
        whichData = [XtypeString in i for i in dataSets]
        # np.where(whichData)
        dataSubSets = list(compress(noBrainsDataSets, whichData))

        # figure width
        # for theRez in (1000000, 500000, 100000, 40000):
        #    print 20 + (10 * (500000/intTheRez))
        figureWidth = 20 + (10 * (500000/intTheRez))

        # Use TrueType fonts that work in illustrator:
        #    http://stackoverflow.com/questions/5956182/cannot-edit-text-in-chart-exported-by-matplotlib-and-opened-in-illustrator
        matplotlib.rcParams['pdf.fonttype'] = 42

        v4CDir = "output/" + contactMapOutputDesc + "/" + outputDir + "/heatmaps_QN/virtual4C.allLoci.lineplots.XaXi"
        if not os.path.exists(v4CDir):
            os.makedirs(v4CDir)

        vlineLoi = ['dxz4', 'firre']

        # Get max Y-vals
        yClipPct = 0.20
        maxMaxY = 0
        binRange = 0
        for whatData in dataSubSets:
            # The QN Data!
            data = dataDict_QN_NaN[whatData].copy()  # VERY IMPORTANT TO COPY OTHERWISE dataDict gets changed!
            binLen = data.shape[0]
            maxY = np.nanmax(data)
            if maxY > maxMaxY:
                maxMaxY = maxY

        xVals = np.r_[0:binLen]
        xVals = (xVals * intTheRez) / (1.0*10**6)

        for currBin in range(binLen):

            if unMappleRegionBool[currBin]:
                continue

            fig, ax = plt.subplots(1, 1, figsize=(figureWidth, 10))
            #plots = {}

            # Plot graphs
            for whatData in dataSubSets:

                inputDir = contactMapIOmetadata[whatData]["inputDir"]
                inputFile = contactMapIOmetadata[whatData]["inputFile"]
                outputDir = contactMapIOmetadata[whatData]["outputDir"]
                outputLabel = contactMapIOmetadata[whatData]["outputLabel"]

                # print whatData + ": " + inputDir + ", " + inputFile + ", " + outputDir + ", " + outputLabel

                # The QN Data!
                data = dataDict_QN_NaN[whatData].copy()  # VERY IMPORTANT TO COPY OTHERWISE dataDict gets changed!
                #print data.shape
                #print np.nanmean(data)

                # get data row
                yVals = data[currBin,]
                # maxY = np.nanmax(yVals)

                # http://stackoverflow.com/questions/20618804/how-to-smooth-a-curve-in-the-right-way
                # http://scipy.github.io/devdocs/generated/scipy.signal.savgol_filter.html#scipy.signal.savgol_filter
                # window size, polynomial order
                yVals = signal.savgol_filter(yVals, 5, 3)

                # if all(np.isnan(yVals)):
                #     continue

                # # Insulation score loci
                # xVals = np.r_[0:len(yVals)]
                # xVals = (xVals * intTheRez) / (1.0*10**6)

                # plots[whatData], = ax.plot(xVals, yVals,
                #                           '-', label=whatData,
                #                           color=contactMapIOmetadata[whatData]['dataColor'],
                #                           linewidth=2,
                #                           linestyle='solid')
                ax.plot(xVals, yVals, '-', label=whatData,
                                        linewidth=4,
                                        color=contactMapIOmetadata[whatData]['dataColor'])

            # X-axis limits
            minX = 0
            maxX = int(np.ceil(float(chromSizes[thechrom]) / 1000000.0))

            ax.grid(True)
            ax.set_ylabel('Contact counts\n(clipped at ' + str(int(round(yClipPct * 100, 0))) + '% of max)', fontsize=14)
            ax.set_xlabel('Mb\n(' + theRez + 'bp bins)', fontsize=14)
            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(12)
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(12)

            ax.set_xlim(minX, maxX)
            ax.set_ylim(0, yClipPct*maxMaxY)

            #ax.legend(loc='best', fontsize=14)
            ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                   ncol=2, mode="expand", borderaxespad=0., fontsize=14)

            ax.axvline(xVals[currBin], color='k', linestyle='--', linewidth=0.5)
            # http://matplotlib.org/examples/pylab_examples/fonts_demo.html
            font0 = FontProperties()
            font1 = font0.copy()
            font1.set_size('16')
            binMidpointLocus = (currBin * intTheRez) / (1.0*10**6)
            ax.text(xVals[currBin], yClipPct*maxMaxY, str(binMidpointLocus) + " Mb", rotation=90, verticalalignment="top", fontproperties=font1)

            for vloi in vlineLoi:
                if vloi in LOI.keys():
                    v_chrCoords = LOI[vloi].split(":")
                    v_theChr = v_chrCoords[0]
                    if v_theChr == thechrom:
                        #print "\tPlot it!"
                        v_startStop = v_chrCoords[1].split("-")
                        v_startStop = map(int, v_startStop)
                        v_geneMidPoint = (v_startStop[0] + (v_startStop[1] - v_startStop[0]))
                        v_currBin = v_geneMidPoint / intTheRez

                        ax.axvline(xVals[v_currBin], color='k', linestyle='--', linewidth=0.5)
                        # http://matplotlib.org/examples/pylab_examples/fonts_demo.html
                        font0 = FontProperties()
                        font1 = font0.copy()
                        font1.set_size('12')
                        ax.text(xVals[v_currBin], yClipPct*maxMaxY, vloi, rotation=90, verticalalignment="top", fontproperties=font1)

            # Tidy up and save
            #fig.subplots_adjust(top=1)
            #fig.tight_layout()
            outputLabel =  "virtual4Cplot.AllSamples." + Xtype + ".Bin" + str(currBin).zfill(4) + "." + theRez + "." + chrom
            title = outputLabel.replace(".", " ")
            fig.suptitle(title, fontsize=20, y=1.01)
            fig.savefig(v4CDir + "/" + outputLabel + ".png", dpi=150, bbox_inches='tight')
            fig.savefig(v4CDir + "/" + outputLabel + ".pdf", dpi=150, bbox_inches='tight')
            plt.close()


################################################
################################################
################################################
# 170623
# Virtual 4C plots based on correlation matrices for all loi - overlayed line plots

if 'v4ccormatrixplots4alllocilinesxaxi' in toRun and toRun['v4ccormatrixplots4alllocilinesxaxi']:

    print "\n\n\nvirtual 4C line plots based on correlation matrices for all loci - overlayed line plots - seperate Xa and Xi"

    noBrainsDataSets = filter(lambda x:'Brain' not in x, dataSets)

    for Xtype in ['Xa', 'Xi']:
        #print(Xtype)
        XtypeString = '_' + Xtype
        whichData = [XtypeString in i for i in dataSets]
        # np.where(whichData)
        dataSubSets = list(compress(noBrainsDataSets, whichData))

        # figure width
        # for theRez in (1000000, 500000, 100000, 40000):
        #    print 20 + (10 * (500000/intTheRez))
        figureWidth = 20 + (10 * (500000/intTheRez))

        # Use TrueType fonts that work in illustrator:
        #    http://stackoverflow.com/questions/5956182/cannot-edit-text-in-chart-exported-by-matplotlib-and-opened-in-illustrator
        matplotlib.rcParams['pdf.fonttype'] = 42

        v4CDir = "output/" + contactMapOutputDesc + "/" + outputDir + "/heatmaps_QN/virtual4C.corrMatrices.allLoci.lineplots.XaXi"
        if not os.path.exists(v4CDir):
            os.makedirs(v4CDir)

        vlineLoi = ['dxz4', 'firre']

        # # Get max Y-vals
        # yClipPct = 0.20
        # maxMaxY = 0
        # binRange = 0
        # *** NEED BIN LEN ***
        for whatData in dataSubSets:
            # The QN Data!
            data = corDataDict_QN_NaN[whatData].copy()  # VERY IMPORTANT TO COPY OTHERWISE dataDict gets changed!
            binLen = data.shape[0]
        #     maxY = np.nanmax(data)
        #     if maxY > maxMaxY:
        #         maxMaxY = maxY
        maxY = 0.5
        minY = -0.5

        xVals = np.r_[0:binLen]
        xVals = (xVals * intTheRez) / (1.0*10**6)

        for currBin in range(binLen):

            if unMappleRegionBool[currBin]:
                continue

            fig, ax = plt.subplots(1, 1, figsize=(figureWidth, 10))
            #plots = {}

            # Plot graphs
            for whatData in dataSubSets:

                inputDir = contactMapIOmetadata[whatData]["inputDir"]
                inputFile = contactMapIOmetadata[whatData]["inputFile"]
                outputDir = contactMapIOmetadata[whatData]["outputDir"]
                outputLabel = contactMapIOmetadata[whatData]["outputLabel"]

                # print whatData + ": " + inputDir + ", " + inputFile + ", " + outputDir + ", " + outputLabel

                # The QN Data!
                data = corDataDict_QN_NaN[whatData].copy()  # VERY IMPORTANT TO COPY OTHERWISE dataDict gets changed!
                #print data.shape
                #print np.nanmean(data)

                # get data row
                # yVals = data[currBin,]
                yVals = data[[currBin]]
                # maxY = np.nanmax(yVals)

                # http://stackoverflow.com/questions/20618804/how-to-smooth-a-curve-in-the-right-way
                # http://scipy.github.io/devdocs/generated/scipy.signal.savgol_filter.html#scipy.signal.savgol_filter
                # window size, polynomial order
                yVals = signal.savgol_filter(yVals, 5, 3)

                # if all(np.isnan(yVals)):
                #     continue

                # # Insulation score loci
                # xVals = np.r_[0:len(yVals)]
                # xVals = (xVals * intTheRez) / (1.0*10**6)

                # plots[whatData], = ax.plot(xVals, yVals,
                #                           '-', label=whatData,
                #                           color=contactMapIOmetadata[whatData]['dataColor'],
                #                           linewidth=2,
                #                           linestyle='solid')
                ax.plot(xVals, yVals, '-', label=whatData,
                                        linewidth=4,
                                        color=contactMapIOmetadata[whatData]['dataColor'])

            # X-axis limits
            minX = 0
            maxX = int(np.ceil(float(chromSizes[thechrom]) / 1000000.0))

            ax.grid(True)
            # ax.set_ylabel('Contact counts\n(clipped at ' + str(int(round(yClipPct * 100, 0))) + '% of max)', fontsize=14)
            ax.set_ylabel('Pearson correlation of contact matrix', fontsize=14)
            ax.set_xlabel('Mb\n(' + theRez + 'bp bins)', fontsize=14)
            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(12)
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(12)

            ax.set_xlim(minX, maxX)
            #ax.set_ylim(0, yClipPct*maxMaxY)
            ax.set_ylim(minY, maxY)

            #ax.legend(loc='best', fontsize=14)
            ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                   ncol=2, mode="expand", borderaxespad=0., fontsize=14)

            ax.axvline(xVals[currBin], color='k', linestyle='--', linewidth=0.5)
            # http://matplotlib.org/examples/pylab_examples/fonts_demo.html
            font0 = FontProperties()
            font1 = font0.copy()
            font1.set_size('16')
            binMidpointLocus = (currBin * intTheRez) / (1.0*10**6)
            #ax.text(xVals[currBin], yClipPct*maxMaxY, str(binMidpointLocus) + " Mb", rotation=90, verticalalignment="top", fontproperties=font1)
            ax.text(xVals[currBin], maxY,  str(binMidpointLocus) + " Mb", rotation=90, verticalalignment="top", fontproperties=font1)

            for vloi in vlineLoi:
                if vloi in LOI.keys():
                    v_chrCoords = LOI[vloi].split(":")
                    v_theChr = v_chrCoords[0]
                    if v_theChr == thechrom:
                        #print "\tPlot it!"
                        v_startStop = v_chrCoords[1].split("-")
                        v_startStop = map(int, v_startStop)
                        v_geneMidPoint = (v_startStop[0] + (v_startStop[1] - v_startStop[0]))
                        v_currBin = v_geneMidPoint / intTheRez

                        ax.axvline(xVals[v_currBin], color='k', linestyle='--', linewidth=0.5)
                        # http://matplotlib.org/examples/pylab_examples/fonts_demo.html
                        font0 = FontProperties()
                        font1 = font0.copy()
                        font1.set_size('12')
                        #ax.text(xVals[v_currBin], yClipPct*maxMaxY, vloi, rotation=90, verticalalignment="top", fontproperties=font1)
                        ax.text(xVals[v_currBin], yMax, vloi, rotation=90, verticalalignment="top", fontproperties=font1)

            # Tidy up and save
            #fig.subplots_adjust(top=1)
            #fig.tight_layout()
            outputLabel =  "virtual4Cplot.corrMatrices.AllSamples." + Xtype + ".Bin" + str(currBin).zfill(4) + "." + theRez + "." + chrom
            title = outputLabel.replace(".", " ")
            fig.suptitle(title, fontsize=20, y=1.01)
            fig.savefig(v4CDir + "/" + outputLabel + ".png", dpi=150, bbox_inches='tight')
            fig.savefig(v4CDir + "/" + outputLabel + ".pdf", dpi=150, bbox_inches='tight')
            plt.close()





# #############################################################################
# #############################################################################
# # MAIN CALL
# if __name__=="__main__":
#     main()
