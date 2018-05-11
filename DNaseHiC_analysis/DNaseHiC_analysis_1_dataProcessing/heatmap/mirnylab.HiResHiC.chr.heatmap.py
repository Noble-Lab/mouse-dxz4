#!/usr/bin/env python

'''
Created on Oct 30, 2012
 
@author: wenxiu
'''
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from mirnylib import genome
from mirnylib import h5dict
from mirnylib import plotting
from hiclib import binnedData
from hiclib import highResBinnedData

from genvis import mshow

USAGE = """USAGE: mirnylab.HiResHiC.chr.heatmap.py <bychr.hdf5file> <outputfilehead> <assembly> <resolution> <title>

  Use Mirnylab's code to plot intra-chr heatmap from <bychr.hdf5file>

"""

##############################################################################
# MAIN
##############################################################################
# Parse the command line.
if (len(sys.argv) != 6):
  sys.stderr.write(USAGE)
  sys.exit(1)
inputFileName = sys.argv[1] # hdf5file
outputFileHead = sys.argv[2]
assembly = sys.argv[3]
resolution = int(sys.argv[4])
title = sys.argv[5]

if assembly == "hg18":
  genome_db = genome.Genome('/net/noble/vol1/home/wenxiu/bin/mirnylab-hiclib-421eff863f12/fasta/hg18', readChrms=['#', 'X'])
elif assembly == "hg19":
  genome_db = genome.Genome('/net/noble/vol1/home/wenxiu/bin/mirnylab-hiclib-421eff863f12/fasta/hg19', readChrms=['#', 'X'])
elif assembly == "mm10":
  genome_db = genome.Genome('/net/noble/vol2/home/gbonora/bin/mirnylab-hiclib-9ed8d9e0ca7f/fasta/mm10', readChrms=['#', 'X'])
else:
  sys.stderr.write("Error! Do not have genome %s\n" % assembly)
  sys.exit(1)
genome_db.setResolution(resolution=resolution)

sys.stdout.write("self.chrmCount=%d\n" % genome_db.chrmCount)
for i in xrange(genome_db.chrmCount):
 sys.stdout.write("self.chrmLens[%d]=%d\n" % (i, genome_db.chrmLens[i]))
sys.stdout.write("self.numBins=%d\n" % genome_db.numBins)

# Read resolution from the dataset.
raw_heatmap = h5dict.h5dict(inputFileName, mode='r')
resolution = int(raw_heatmap['resolution'])

# Create a binnedData object, load the data.
#BD = binnedData.binnedData(resolution, genome_db)
#BD.simpleLoad(inputFileName,outputFileHead)
BD = highResBinnedData.HiResHiC(genome_db, resolution)
BD.loadData(inputFileName)
BD.removeDiagonal()
#plt.figure()
#plotting.plot_matrix(np.log(BD.dataDict[outputFileHead]), clip_min=-1, clip_max=6)
#plotting.plot_matrix(np.log(BD.dataDict[outputFileHead]))
#plt.savefig('%s.log.png' % outputFileHead)
#plt.clf()
#plotting.plot_matrix(BD.dataDict[outputFileHead], clip_min=-1, clip_max=64)
#plotting.plot_matrix(BD.dataDict[outputFileHead])
#plt.savefig('%s.png' % outputFileHead)
#plt.clf()

my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap', mshow.cdicts[0], 256)

#plot each chr
for i in xrange(genome_db.chrmCount):
  if (i,i) not in BD.data.keys():
    sys.stderr.write("Error! Can not find chr %d\n" %(d+1))
    sys.exit(1)
  chrdata = BD.data[(i,i)].getData() # get copy in dense format

  vmin = chrdata.min()
  vmax = chrdata.max()

#  plt.imshow(chrdata, interpolation="nearest", vmin=vmin, vmax=vmax)
#  plt.colorbar()
#  plt.xlabel("chromosome %d" % (i+1))
#  plt.ylabel("chromosome %d" % (i+1))
#  plt.savefig('%s.chr%d.imshow.png' % (outputFileHead, i+1))
#  plt.clf()

#   if resolution <= 500000:
#     plotting.plot_matrix(chrdata, cmap=my_cmap)
#   else:
#     plotting.plot_matrix(chrdata)
  plotting.plot_matrix(chrdata)
  plt.xlabel("chromosome %d" % (i+1))
  plt.ylabel("chromosome %d" % (i+1))
  plt.title("%s.chr%d" % (title, (i+1)))
  plt.savefig('%s.chr.chr%d.png' % (outputFileHead, i+1))
  plt.clf()
  
#   if resolution <= 500000:
#     plotting.plot_matrix(np.log(chrdata), cmap=my_cmap)
#   else:
#     plotting.plot_matrix(np.log(chrdata))
  plotting.plot_matrix(np.log(chrdata))
  plt.xlabel("chromosome %d" % (i+1))
  plt.ylabel("chromosome %d" % (i+1))
  plt.title("%s.chr%d.log" % (title, (i+1)))
  plt.savefig('%s.chr.chr%d.log.png' % (outputFileHead, i+1))
  plt.clf()

  chrFileName = "%s.chr%d.txt.gz" % (outputFileHead, i+1)
  np.savetxt(chrFileName, chrdata, delimiter="\t")
