#!/usr/bin/env python
'''
Created on Aug 08, 2012

@author: ferhat
'''
import sys

USAGE = """USAGE: sort_loci.py <interactionsFile> <numberOfLoci>

Reads the interactions file which contains rows that correspond to interactions between 2 or more loci and sorts these loci according to their chromosome names first and then according to their start indices and then according to their strands. 

<interactionsFile> is the main input file that contains interactions in this format (an example of 3 loci):
"FCD059NACXX:5:2107:18137:185323 +       chr10   100000018       ref  GAA..       +       chr5    66808096  alt ATG.. -       chr5    808096        both  ATC.." 

Output will be of the same format but just sorted. For teh above example:
"FCD059NACXX:5:2107:18137:185323 -       chr5    808096 both  ATC.. +       chr5    66808096    alt ATG.. +       chr10   100000018       ref GAA.."

"""


def main(infilename, numberOfLoci):
  if (infilename == "-"):
    infile = sys.stdin
  else:
    infile = open(infilename, "r")
  linecount=0
  for inline in infile:
    if len(inline)<2:
      sys.stderr.write(repr(linecount)+" lines were read from the file. ")
      return
    linecount+=1
    tokenizedStr=inline.rstrip().split()
    readID=tokenizedStr[0]
    strandSigns=[tokenizedStr[i] for i in range(1,5*numberOfLoci,5)]
    tempchrs=[tokenizedStr[i] for i in range(2,5*numberOfLoci,5)]
    readStarts=[int(tokenizedStr[i]) for i in range(3,5*numberOfLoci,5)]
    genomes=[tokenizedStr[i] for i in range(4,5*numberOfLoci+1,5)]
    seqs=[tokenizedStr[i] for i in range(5,5*numberOfLoci+1,5)]
    tempchrs,readStarts,strandSigns,genomes,seqs=zip(*sorted(zip(tempchrs,readStarts,strandSigns,genomes,seqs)))
    print(readID+"\t"),
    for i in range(numberOfLoci):
      print(strandSigns[i]+"\t"+tempchrs[i]+"\t"+repr(readStarts[i])+"\t"+genomes[i]+"\t"+seqs[i]),
    print
  sys.stderr.write(repr(linecount)+" lines were read from the file. \n")

if __name__ == "__main__":
  if (len(sys.argv) != 3):
    sys.stderr.write(USAGE)
    sys.exit(1)
  main(str(sys.argv[1]),int(sys.argv[2]))

