#!/usr/bin/env python

'''
Created on Aug 19 2013
@author: wenxiu
'''
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.colors import LogNorm

from mirnylib import genome
from mirnylib import h5dict
from mirnylib import plotting
from hiclib import binnedData
from hiclib import highResBinnedData

USAGE = """USAGE: mirnylab.ICE.heatmap.py <matrix> <biases> <outputfilehead> <assembly> <resolution> <title>

  Modified the Mirnylab's ICE code, use the biases vector to correct the interaction matrix.
  Plot heatmap before and after normalization. 

"""
##############################################################################
def saveHeatmap(inputfilename, outputfilename, resolution=1000000,
                    countDiagonalReads="Once"):
        """
        Saves heatmap to filename at given resolution.

        Parameters
        ----------
        filename : str
            Filename of the output h5dict
        resolution : int
            Resolution of an all-by-all heatmap
        countDiagonalReads : "once" or "twice"
            How many times to count reads in the diagonal bin
        
        """

        try:
            os.remove(outputfilename)
        except:
            pass

        tosave = h5dict.h5dict(path=outputfilename, mode="w")
        heatmap = np.loadtxt(inputfilename, dtype=float)
        tosave["heatmap"] = heatmap
        del heatmap
        tosave["resolution"] = resolution
        print "----> Heatmap saved to '%s' at %d resolution" % (
            outputfilename, resolution)

##############################################################################
# MAIN
##############################################################################
# Parse the command line.
if (len(sys.argv) != 7):
  sys.stderr.write(USAGE)
  sys.exit(1)
inputFileName = sys.argv[1]
biasesFileName = sys.argv[2]
outputFileHead = sys.argv[3]
assembly = sys.argv[4]
resolution = int(sys.argv[5])
title = sys.argv[6]

hdf5FileName = "%s.chr.hdf5" % outputFileHead
#saveHeatmap(inputFileName, hdf5FileName, resolution=resolution)

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
#raw_heatmap = h5dict.h5dict(inputFileName, mode='r')
#resolution = int(raw_heatmap['resolution'])

# Create a binnedData object, load the data.
#BD = binnedData.binnedData(resolution, genome_db)
#BD.simpleLoad(inputFileName,outputFileHead)
BD = highResBinnedData.HiResHiC(genome_db, resolution)
BD.loadData(inputFileName)
BD.removeDiagonal()
BD.export(hdf5FileName)
#plt.figure()
if resolution > 100000:
  #plotting.plot_matrix(np.log(BD.dataDict[outputFileHead]), clip_min=-1, clip_max=6)
  #plotting.plot_matrix(np.log(BD.getCombinedMatrix()))
  plotting.plot_matrix(BD.getCombinedMatrix(), norm=LogNorm())
  plt.title("%s.raw.log" % title)
  plt.savefig('%s.log.png' % outputFileHead)
  plt.clf()
  #plotting.plot_matrix(BD.dataDict[outputFileHead], clip_min=-1, clip_max=64)
  plotting.plot_matrix(BD.getCombinedMatrix())
  plt.title("%s.raw" % title)
  plt.savefig('%s.png' % outputFileHead)
  plt.clf()

# Remove the contacts between loci located within the same bin.
BD.removeDiagonal()
# Remove bins with less than full of a bin sequenced.
#BD.removeBySequencedCount(0.5)
# Remove 2% of regions with low coverage.
BD.removePoorRegions(percent=2)
# Truncate top 0.05% of interchromosomal counts (possibly, PCR blowouts).
#BD.truncTrans(high=0.0005)

# Perform iterative correction.
#BD.iterativeCorrectWithoutSS(force=True)
#BD.iterativeCorrection()
biases = np.loadtxt(biasesFileName, dtype=float)
BD.divideByBiases(biases)
# Save the iteratively corrected heatmap.
hdf5CorrectedFileName = "%s.IC.byBiases.chr.hdf5" % outputFileHead
BD.export(hdf5CorrectedFileName)
#np.savetxt('%s.IC.full.matrix.txt' % outputFileHead, BD.dataDict[outputFileHead], fmt='%.2f', delimiter='\t')

# Plot the heatmap directly.
#plt.figure()
if resolution > 100000:
  #plotting.plot_matrix(np.log(BD.dataDict[outputFileHead]), clip_min=-1, clip_max=6)
  #plotting.plot_matrix(np.log(BD.getCombinedMatrix()))
  plotting.plot_matrix(BD.getCombinedMatrix(), norm=LogNorm())
  plt.title("%s.corrected.log" % title)
  #plt.show()
  plt.savefig('%s.IC.byBiases.log.png' % outputFileHead)
  #plt.figure()
  plt.clf()
  #plotting.plot_matrix(BD.dataDict[outputFileHead],clip_max=32)
  plotting.plot_matrix(BD.getCombinedMatrix())
  plt.title("%s.corrected" % title)
  plt.savefig('%s.IC.byBiases.png' % outputFileHead)
  plt.clf()

