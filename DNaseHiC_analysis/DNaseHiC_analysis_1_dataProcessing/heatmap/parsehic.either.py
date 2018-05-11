#!/usr/bin/env python

import sys
#import _normalization as ICE
import chromsizes
import numpy as np
import os
import gzip
from mirnylib import genome
from mirnylib import h5dict
from hiclib import binnedData

"""
This script is to be used as for compiling hi-c data into a matrix of
chromosome contacts. It can read in a hi-c data file, and parse the
locational-contact data into 'bins' of variable size, though the default
is to use 1mb bins to sub-divide each chromosome

Intended for 1-off use on hi-c data, so that the contact matrices may
simply be accessed by other scripts, and they may avoid the time-consuming
process of parsing the data sets.

NOTE: takes in an entire dir of hi-c files, and processes EACH
"""


###########
# helpers #
###########


# useful utility for bin-sizes; handles spill-over with integer division,
# returns a/b rounded toward infinity
def round_up_division(a, b): return (a+(-a%b))//b


def parse_contacts(hi_c_file, bin_starts, bin_size, chrom_sizes, genomePath, output_file, output_chr_file):
  # NOTE: makes assumption that we are using chromosome X, and that X is the
  #       last set of 'bins' in the matrix (i.e. no chrY, and after the numbers)
  matrix_size = bin_starts['X'] + round_up_division(chrom_sizes['X'], bin_size)
  # create a matrix_size * matrix_size ndarray filled with 0s, and use
  # unsigned 32-bit integers to save space (counts should never be big enough
  # to exceed 2^32 - 1
  counts = np.zeros((matrix_size,matrix_size), dtype=np.uint32)
  with gzip.open(hi_c_file, 'r') as file:
    i = 0
    for line in file:
      tokens = line.split()
      # extract chromosome name/letter of both contact points
      chrom1 = tokens[0].strip("chr")
      chrom2 = tokens[2].strip("chr")

      # we're ignoring the y-chromosome, and intra-chr interactions
      if chrom1 == 'Y' or chrom2 == 'Y':
        continue
#      if all(chromosome in bin_starts.keys() for chromosome in (chrom1, chrom2)):
      if chrom1 in bin_starts.keys() and chrom2 in bin_starts.keys():
#        sys.stderr.write(line)
        # figure out which bins we had contact in
        bin1 = bin_starts[chrom1] + (int(tokens[1])/bin_size)
        bin2 = bin_starts[chrom2] + (int(tokens[3])/bin_size)

        # mark the contact in our matrix
        counts[bin1, bin2] = counts[bin1, bin2] + int(tokens[4])
        counts[bin2, bin1] = counts[bin1, bin2]
#    np.savetxt(output_file, counts)
    toexport = {}
    toexport["heatmap"] = counts
    toexport["resolution"] = bin_size
    myh5dict = h5dict.h5dict(output_file, mode="w")
    myh5dict.update(toexport)
    genome_db = genome.Genome(genomePath, readChrms=['#', 'X'])
    genome_db.setResolution(resolution=bin_size)
    #Create a binnedData object, load the data and save by chr
    BD = binnedData.binnedData(bin_size, genome_db)
    BD.dataDict["test"] = counts
    BD.export("test", output_chr_file, byChromosome=True)

#############
# arguments #
#############
"""
ARGS:
  arg 1 - filepath to the file containing chromosome sizes
  arg 2 - either an "H", or an "M", specifying [H]uman or [M]ouse data
  arg 3 - filepath to the file containing hi-c contact data files
  arg 4 - filepth to a DIR where output files should be dumped
  arg 5 - genomePath for mirny lab genome object
  arg 6 - number representing size of bins in resulting contact
              matrix, measured in # of base pairs; default is 1 mb
"""

# check correcntess
if len(sys.argv) != 7:
  sys.exit("Error in " + sys.argv[0] + ":\n\t Usage: python " + sys.argv[0] + " <chrom sizes file> <[H]uman or [M]ouse> <dir containing hi-c files> <output dir> <bin size>)")

# first arg - filename of chromosome sizes
chrom_sizes_file = sys.argv[1]

# second arg - specifies human or mouse
H_or_M = sys.argv[2].upper()
if H_or_M not in ("H", "M"):
  sys.exit("Error in " + sys.argv[0] + ":\n\t 2nd arg must be an \"H\" or an \"M\"")

# third arg - directory containing hi-c data
interaction_count_file = sys.argv[3]

# fourth arg - directory for output
output_dir = sys.argv[4]

genomePath = sys.argv[5]

# optional fifth arg - bin size (default 1mb)
bin_size = int(sys.argv[6])


##########
# driver #
##########

# get the sizes of our chromosomes from the UCSC file
chrom_sizes = chromsizes.get_chrom_sizes(chrom_sizes_file, H_or_M)

# get the range of bins for each chromosome in our matrix
bin_starts = chromsizes.get_bin_starts(chrom_sizes, bin_size)

# calculate and output the contact matrices for every hi-c file
file = os.path.basename(interaction_count_file)
if file[-3:] == ".gz":
  file = file[:-3]
output_file = os.path.join(output_dir, file + ".hdf5")
output_chr_file = os.path.join(output_dir, file + ".chr.hdf5")
counts = parse_contacts(interaction_count_file, bin_starts, bin_size,
                        chrom_sizes, genomePath, output_file, output_chr_file)
