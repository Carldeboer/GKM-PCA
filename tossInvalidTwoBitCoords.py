#!/usr/bin/env python3
import warnings
import MYUTILS
import re
import sys
import argparse
parser = argparse.ArgumentParser(description='Goes through the provided set of genomic coordinates and throws out any that are for invalid chromosomes, or are outside the bounds of a chromosome.')
parser.add_argument('-i',dest='inFP', metavar='<inFile>',help='Input file of coordinates', required=True);
parser.add_argument('-s',dest='chromSizes', metavar='<chromSizes>',help='Input file of chrom.sizes', required=True);
parser.add_argument('-o',dest='outFP', metavar='<outFile>',help='Where to output results [default=stdout]', required=False);
parser.add_argument('-l',dest='logFP', metavar='<logFile>',help='Where to output errors/warnings [default=stderr]', required=False);
parser.add_argument('-v',dest='verbose', action='count',help='Verbose output?', required=False, default=0);

args = parser.parse_args();



chromSizes={}
chromSizesFile=MYUTILS.smartGZOpen(args.chromSizes,'r');
for line in chromSizesFile:
  if line is None or line == "" or line[0]=="#": continue
  data=line.rstrip().split("\t");
  chromSizes[data[0]]=int(data[1]);


if (args.logFP is not None):
  logFile=MYUTILS.smartGZOpen(args.logFP,'w');
  sys.stderr=logFile;

if (args.outFP is None):
  outFile= sys.stdout;
else:
  if args.verbose>0: warnings.warn("Outputting to file "+args.outFP);
  outFile = MYUTILS.smartGZOpen(args.outFP,'w');

inFile=MYUTILS.smartGZOpen(args.inFP,'r');
invalidChr= 0;
outOfBounds=0;
for line in inFile:
  if line is None or line == "" or line[0]=="#": continue
  m = re.search("([^:]*):([-0-9]*)-([-0-9])",line);
  if m:
    chr= m.group(1);
    st= int(m.group(2));
    en= int(m.group(3));
    if chr in chromSizes:
      if st>=0 and en<=chromSizes[chr]:
        outFile.write(line);
      else:
        outOfBounds+=1;
    else:
      invalidChr+=1;
  else:
    raise Exception("Failed to match regular expression at line '%s'" %(line));
inFile.close();
outFile.close();
sys.stderr.write("Tossed %i invalid chromosomes and %i out of bounds entries\n"%(invalidChr,outOfBounds))
if (args.logFP is not None):
  logFile.close();
