#!/usr/bin/python3
#
# make_template.py - Creates templates for use with 23andme_to_vcf.py
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
# make_template.py [outputTemplateFile] [mapFile] [snpListVCF]



import csv
import re
import sys
from collections import defaultdict
import tabix
import gzip
from vcflib import *

templatefilename = sys.argv[1]
mapfilename = sys.argv[2]
reffilename = sys.argv[3]




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

  # replace chromosome id's with numbers for easy sorting
  if REFCHROM == "X":
    REFCHROM = 23
  elif REFCHROM == "Y":
    REFCHROM = 24
  elif REFCHROM == "MT":
    REFCHROM = 25
  else:
    REFCHROM = int(REFCHROM)

  REFPOS = int(REFPOS)  
  if REFRSID != '.':
    markerMap[REFCHROM][REFPOS] = REFRSID

mapfile.close()


# open output template
templatefile = open(templatefilename, 'w')

# open ref vcf and write headers to output
headerstring = getVCFheaders(reffilename)
templatefile.write(headerstring)
tb = tabix.open(reffilename)

for chrom in sorted(markerMap):
  for pos in sorted(markerMap[chrom]):
    rsid = markerMap[chrom][pos]

    # replace chrom numbers with canonical names
    qchrom = chrom
    qchrom = str(qchrom)
    if qchrom == '23':
      qchrom = 'X'
    elif qchrom == '24':
      qchrom = 'Y'
    elif qchrom == '25':
      qchrom = 'MT'

    query  = qchrom + ":" + str(pos) + "-" + str(pos)
    recorditer = tb.querys(query)

    # convert iterator to actual array of records
    records = []
    for record in recorditer:
      records.append(record)

    records = getIDmatches(rsid, records)
    if len(records) == 1:
      record = records[0]
      templatefile.write('\t'.join(record)+'\n')
    else:
      print("not found: " + ' '.join([chrom,pos,rsid]))

templatefile.close()

