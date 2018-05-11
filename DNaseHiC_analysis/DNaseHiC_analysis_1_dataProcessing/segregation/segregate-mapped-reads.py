#!/usr/bin/env python
# AUTHOR: William Stafford Noble
# CREATE DATE: 18 June 2012

# # AUTHOR: Wenxiu Ma
# MODIFIED DATE: 10 August 2012
# MODIFIED DATE: 10 Feb 2015 # I am removing the ref/alt initially assignment since the read might not contain SNP at all.

# AUTHOR: Giancarlo Bonora
# MODIFIED DATE: 20 Apr 2016 # Parsed out BWA-MEM primary alignments and unmapped reads to ensure all reads in place.
#                            # Simplified mapping check to used only MAPQ.

# AUTHOR: Giancarlo Bonora
# MODIFIED DATE: 15 June 2017 # Consider CIGAR string

import sys
import string
import gzip
import re

USAGE = """USAGE: segregate-mapped-reads.py <SAM1.gz> <SAM2.gz> <SNPs.gz> <root>

  This program takes as input two SAM files corresponding to the same
  set of reads mapped to two different genomes.  The program
  segregates only the uniquely mapped and high quality reads into four 
  groups, based on whether they map to both genomes, neither genome, 
  or only one of the two.  This division
  is first accomplished by simpling checking the SAM file to see which
  genome(s) each read maps to.  For reads that map to both, the SNP
  file is used to distinguish between SNPs and mismatches.

  Specifically, given a read that initially maps to both genomes, the
  program considers each SNP that falls within the given read.  If the
  SAM file shows the reference allele at the given position, then the
  read is marked to map to the reference genome (the first SAM file on
  the command line).  Conversely, if the alternative allele appears in
  the read, then the read is marked to map to the alternative genome.
  If, after considering all SNPs within a read, it is marked to map to
  one genome but not the other, then the read is placed into that
  output file.

  The program produces the following output files:

    <root>.{either, both-ref,both-alt,neither,ref,alt}.sam

  The "both-ref" and "both-alt" files contain the same reads but with
  mapped to the two different genomes. The "either" is the union of 
  "both-ref", "both-alt", "ref" and "alt".

  The SNP file contains the following five columns:
    - chromosome number
    - position (indexed from 1)
    - strand (ignored)
    - reference base
    - SNP base

"""

#############################################################################
# Used to take the complement of a string.
complement = string.maketrans('ATCGN', 'TAGCN')


def reverseComplement(sequence):
    return sequence.lower().translate(complement)[::-1]


#############################################################################
# MAIN
#############################################################################

# Parse the command line.
if (len(sys.argv) != 5):
    sys.stderr.write(USAGE)
    sys.exit(1)
sam1FileName = sys.argv[1]
sam2FileName = sys.argv[2]
snpFileName = sys.argv[3]
root = sys.argv[4]

# Read the SNP file into a dictionary.
snps = {}  # key = (chrom, position), value = (ref allele, alt allele)
snpFile = gzip.open(snpFileName, "r")
for line in snpFile:
    if (line[0] == '#'):
        continue
    words = line.rstrip().split("\t")
    if (len(words) < 5):
        sys.stderr.write("Error reading from %s.\n%s" % (snpFileName, line))
        sys.exit(1)
    chrom = words[0]
    position = int(words[1])
    snps[(chrom, position)] = (words[3], words[4])
snpFile.close()
sys.stderr.write("Read %d SNPs from %s.\n" % (len(snps), snpFileName))

# Open the output files.
bothRefFile = open("%s.both-ref.sam" % root, "w")
# bothAltFile = open("%s.both-alt.sam" % root, "w")
# neitherFile = open("%s.neither.sam" % root, "w")
refFile = open("%s.ref.sam" % root, "w")
altFile = open("%s.alt.sam" % root, "w")

# Initialize counters.
numBoth = 0
numRef = 0
numAlt = 0
numSnpRef = 0
numSnpAlt = 0
numNeither = 0
numDiffLoc = 0

# Traverse the two SAM files in parallel.
sam1File = gzip.open(sam1FileName, "r")
sam2File = gzip.open(sam2FileName, "r")
lineNumber = 0
for refLine in sam1File:
    lineNumber += 1
    altLine = sam2File.readline()

    # Skip header lines.
    if ((refLine[0] == '@') or (altLine[0] == '@')):
        bothRefFile.write(refLine)
        #    bothAltFile.write(refLine)
        #    neitherFile.write(refLine)
        refFile.write(refLine)
        altFile.write(refLine)
        continue

    # Verify that we have the right number of fields.
    refWords = refLine.split("\t")
    if (len(refWords) < 11):
        sys.stderr.write("Error parsing line %d from %s.\n%s"
                         % (lineNumber, sam1FileName, refLine))
        sys.exit(1)
    altWords = altLine.split("\t")
    if (len(altWords) < 11):
        sys.stderr.write("Error parsing line %d from %s.\n%s"
                         % (lineNumber, sam2FileName, altLine))
        sys.exit(1)

    # Verify that the two lines are for the same read.
    if (refWords[0] != altWords[0]):
        sys.stderr.write("Read mismatch at line %d. (%s != %s).\n"
                         % (lineNumber, refWords[0], altWords[0]))
        sys.exit(1)

    # # Get the mapping status for each genome.
    # # (col#2: 4 = unmapped) (col#5: MAPQ, highest 37 or 60) (col#12 XT:A:U uniquelly mapped)
    # # (col#6 CIGAR) (col#17 XO gap openings) (col#18 XG gap extension)
    # Mapped or not?
    # Already parsed for primary or unmapped. Just filter on MAPQ.
    refUnmapped = bool(
        bool(int(refWords[1]) & 4) or
        bool(int(refWords[4]) < 30))
    altUnmapped = bool(
        bool(int(altWords[1]) & 4) or
        # bool(int(refWords[4]) < 30))
        # *** 170615 ***
        bool(int(altWords[4]) < 30))

    # Get the read and its coordinates.
    refStartPosition = int(refWords[3])
    refChrom = refWords[2][3:]
    altStartPosition = int(altWords[3])
    altChrom = altWords[2][3:]
    # no need check strand of the read since col#10 is already the complement
    refRead = refWords[9]
    altRead = altWords[9] # Should be same as refRead
    readLength = len(refRead)
    #  if (bool(int(refWords[1]) & 0x0010)):
    #    refRead = reverseComplement(refWords[9])
    #    sys.stderr.write("Reverse complement of is %s: %s" % (refRead, refLine))

    # Get CIGAR string
    refCIGAR = refWords[5]
    altCIGAR = altWords[5]

    # Parse ref CIGAR
    #refCIGAR = '42S108M'; refRead='AACCTTGTTGTCTATCAACTGATGAATGGAAAGCTGAGGGATAGTCTCTATAGGGTTAGGCATATCCTCTCCCACGAGGCCAAAAAAGGCAGTTCTCTGCTACATGTTCCTAGGGCCTTGGACCAGCCTATGTGTGCTTTTTGATTGGTG'
    #refCIGAR = '98M52S'; refRead='AAAACAACCAGAATGGCATAACGTGCATTTTTTCTCTCTGCTCAAGAGGACTGACTGAGTTGACCCTGTGTGATTTAAATTTGTGATGGAAGAATTTAAGCTGAGGGATAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAA'
    #refCIGAR = '1S1D50M'; refRead='AAAACAACCAGAATGGCATAACGTGCATTTTTTCTCTCTGCTCAAGAGGACTGACTGAGTTGACCCTGTGTGATTTAAATTTGTGATGGAAGAATTTAAGCTGAGGGATAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAA'
    #refCIGAR = '1S3I1D50M'; refRead='AAAACAACCAGAATGGCATAACGTGCATTTTTTCTCTCTGCTCAAGAGGACTGACTGAGTTGACCCTGTGTGATTTAAATTTGTGATGGAAGAATTTAAGCTGAGGGATAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAA'
    refReadCIGARsed = ''
    prevPos = 0
    if refCIGAR != '*':
        while len(refCIGAR)>0:
            m = re.search('^(\d+)(\D)', refCIGAR)
            mcount = int(m.group(1))
            mtype = m.group(2)
            if mtype == 'M':
                refReadCIGARsed = refReadCIGARsed + refRead[prevPos:(prevPos+mcount)]
                prevPos += mcount
            elif mtype in ['S', 'D']:
                prevPos += mcount
            elif mtype == 'I':
                refReadCIGARsed = refReadCIGARsed + 'N' * mcount
            refCIGAR = refCIGAR.replace(m.group(0), '', 1)
    else:
        refReadCIGARsed = refRead
    refReadCIGARsedLength = len(refReadCIGARsed)

    # Parse alt CIGAR
    altReadCIGARsed = ''
    prevPos = 0
    if altCIGAR != '*':
        while len(altCIGAR)>0:
            m = re.search('^(\d+)(\D)', altCIGAR)
            mcount = int(m.group(1))
            mtype = m.group(2)
            if mtype == 'M':
                altReadCIGARsed = altReadCIGARsed + altRead[prevPos:(prevPos+mcount)]
                prevPos += mcount
            elif mtype in ['S', 'D']:
                prevPos += mcount
            elif mtype == 'I':
                altReadCIGARsed = altReadCIGARsed + 'N' * mcount
            altCIGAR = altCIGAR.replace(m.group(0), '', 1)
    else:
        altReadCIGARsed = altRead
    altReadCIGARsedLength = len(altReadCIGARsed)

    if (refUnmapped) and (altUnmapped):
        #    neitherFile.write(refLine)
        numNeither += 1

    elif (refUnmapped) and (not altUnmapped):
        # sys.stderr.write("Alt Initially!\n%s%s\n" % (refLine, altLine))
        # check alt SNP
        refSNP = False
        for offset in range(0, refReadCIGARsedLength):
            key = (altChrom, altStartPosition + offset)
            if snps.has_key(key):
                if (snps[key][0] == refReadCIGARsed[offset]):
                    refSNP = True
                    # sys.stderr.write("Reference SNP %s at offset %d, position %d.\n"
                    #                  % (refRead[offset], offset, refStartPosition + offset))
        altSNP = False
        for offset in range(0, altReadCIGARsedLength):
            key = (altChrom, altStartPosition + offset)
            if snps.has_key(key):
                if (snps[key][1] == altReadCIGARsed[offset]):
                    altSNP = True
                    # sys.stderr.write("Alternative SNP %s at %d, position %d.\n"
                    #                  % (refRead[offset], offset, refStartPosition + offset))
        if (not refSNP) and altSNP:
            # sys.stderr.write("Alt SNP confirmed!\n%s%s\n" % (refLine, altLine))
            altFile.write(altLine)
            numAlt += 1

    elif (not refUnmapped) and (altUnmapped):
        # sys.stderr.write("Ref Initially!\n%s%s\n" % (refLine, altLine))
        # check ref SNP
        refSNP = False
        for offset in range(0, refReadCIGARsedLength):
            key = (refChrom, refStartPosition + offset)
            if snps.has_key(key):
                if (snps[key][0] == refReadCIGARsed[offset]):
                    refSNP = True
                    # sys.stderr.write("Reference SNP %s at offset %d, position %d.\n"
                    #                  % (refRead[offset], offset, refStartPosition + offset))
        altSNP = False
        for offset in range(0, altReadCIGARsedLength):
            key = (refChrom, refStartPosition + offset)
            if snps.has_key(key):
                if (snps[key][1] == altReadCIGARsed[offset]):
                    altSNP = True
                    # sys.stderr.write("Alternative SNP %s at %d, position %d.\n"
                    #                  % (refRead[offset], offset, refStartPosition + offset))
        if refSNP and (not altSNP):
            # sys.stderr.write("Ref SNP confirmed!\n%s%s\n" % (refLine, altLine))
            refFile.write(refLine)
            numRef += 1

    else:
        # If the read maps to both genomes, try to use SNPs to disambiguate....

        # check if the reads are mapped to the same chr/loc
        if (refChrom != altChrom or refStartPosition != altStartPosition):
            numDiffLoc += 1
            numNeither += 1
            continue

        # sys.stderr.write("Disambiguating %s ...\n" % refWords[0])
        #    sys.stderr.write("Read=%s.\n" % refRead)

        # Scan for SNPs in this read.
        refSNP = False
        for offset in range(0, refReadCIGARsedLength):
            key = (refChrom, refStartPosition + offset)
            if snps.has_key(key):
                if (snps[key][0] == refReadCIGARsed[offset]):
                    refSNP = True
                    # sys.stderr.write("Reference SNP %s at offset %d, position %d.\n"
                    #                  % (refRead[offset], offset, refStartPosition + offset))
        altSNP = False
        for offset in range(0, altReadCIGARsedLength):
            key = (refChrom, altStartPosition + offset)
            if snps.has_key(key):
                if (snps[key][1] == altReadCIGARsed[offset]):
                    altSNP = True
                    # sys.stderr.write("Alternative SNP %s at %d, position %d.\n"
                    #                  % (refRead[offset], offset, refStartPosition + offset))
                # else:
                #     sys.stderr.write("Neither Reference nor Alternative SNP %s at %d, position %d.\n"
                #                      % (refRead[offset], offset, refStartPosition + offset))

        # Decide where to print this read.
        if (refSNP and (not altSNP)):
            # sys.stderr.write("Reference.\n")
            refFile.write(refLine)
            numSnpRef += 1
        elif ((not refSNP) and altSNP):
            # sys.stderr.write("Alternative.\n")
            altFile.write(altLine)
            numSnpAlt += 1
        else:
            # sys.stderr.write("Both.\n")
            bothRefFile.write(refLine)
            # bothAltFile.write(altLine)
            numBoth += 1  # Close all files.

sam1File.close()
sam2File.close()
bothRefFile.close()
refFile.close()
altFile.close()

# Print a summary.
sys.stdout.write("%d mapped to both.\n" % numBoth)
sys.stdout.write("%d mapped to reference initially.\n" % numRef)
sys.stdout.write("%d mapped to reference based on SNPs.\n" % numSnpRef)
sys.stdout.write("%d mapped to alternative initially.\n" % numAlt)
sys.stdout.write("%d mapped to alternative based on SNPs.\n" % numSnpAlt)
sys.stdout.write("%d mapped to neither.\n" % numNeither)
sys.stdout.write("%d mapped to different chr/pos.\n" % numDiffLoc)
