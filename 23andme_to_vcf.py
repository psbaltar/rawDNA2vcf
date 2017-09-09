#!/usr/bin/python3
#
# 23andme_to_vcf.py - Converts 23AndMe format to Variant Call Format
# Copyright (C) 2017 Philip Baltar (psbaltar@gmail.com)
#
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#
# 23andme_to_vcf.py [23andmeInputFile] [outputSampleName] [mapFile] [templateFile] [VCFoutputFile]


import csv
import re
import sys
from collections import defaultdict
import tabix
import gzip
from vcflib import *


infilename = sys.argv[1]
samplename = sys.argv[2]
mapfilename = sys.argv[3]
templatefilename = sys.argv[4]
vcffilename = sys.argv[5]


# open input file and load into memory
infile = open(infilename, 'r')
inreader = csv.reader(infile, delimiter='\t')

indata = {}
for line in inreader:
  if re.search("^#", line[0]):
    continue
  rsid = line[0]
  indata[rsid] = line

infile.close()



# open mapping file and read into memory
mapfile = open(mapfilename, 'r')
mapreader = csv.reader(mapfile, delimiter='\t')

markerMap = defaultdict(dict)
for line in mapreader:
  if re.search("^#", line[0]) or line[5] == '.':
    continue

  CHROM = line[0]
  POS = line[1]
  RSID = line[2]
  REFCHROM = line[3]
  REFPOS = line[4]
  REFRSID = line[5]

  markerMap[REFRSID] = RSID

mapfile.close()


# open template
if re.search("\.gz$", templatefilename):
  templatefile = gzip.open(templatefilename, 'rt')
else:
  templatefile = open(templatefilename, 'r')

templatereader = csv.reader(templatefile, delimiter='\t')

# open vcf and write headers to output
vcffile = open(vcffilename, 'w')

headers = ["##fileformat=VCFv4.1",
           "##FORMAT=<ID=GT,Number=1,Type=Integer,Description=\"Genotype\">",
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+samplename+'\n']
vcffile.write('\n'.join(headers))


# step through the template and output any entries that are also in the input file
for line in templatereader:
  if re.search("^#", line[0]):
    continue

  chrom = line[0]
  pos = line[1]
  rsid = line[2]
  ref = line[3]
  alt = line[4]

  # if template entry is in the map
  if rsid in markerMap:
    inrsid = markerMap[rsid]
    if inrsid in indata:
      genotype = indata[inrsid][3] # genotype from input file

      GT = calcGT(genotype, ref, alt)

      if 'E' not in GT:
        # blank out QUAL, FILTER, and INFO columns
        line[5] = '.'
        line[6] = '.'
        line[7] = '.'
        line.append('GT')
        line.append(GT)
        vcffile.write('\t'.join(line)+'\n')
      else:
        print("Error: "+' '.join([chrom,pos,rsid,ref,alt,inrsid,genotype]))



templatefile.close()
vcffile.close()

