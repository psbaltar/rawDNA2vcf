#!/usr/bin/python3
#
# make_map.py - Creates filters for use with make_template.py
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
# make_filter.py [outputMapFile] [filterFile] [snpListVCF]

import csv
import re
import sys
from collections import defaultdict
import tabix
import gzip
from vcflib import *


mapfilename = sys.argv[1]
infilename = sys.argv[2]
reffilename = sys.argv[3]



infile = open(infilename, 'r')
infilereader = csv.reader(infile, delimiter='\t')

tb = tabix.open(reffilename)

mapfile = open(mapfilename, 'w')
mapfile.write("#CHROM\tPOS\tRSID\tREFCHROM\tREFPOS\tREFRSID\n")




for line in infilereader:
  if re.search("^#", line[0]):
    continue

  chrom = line[0]
  pos = line[1]
  rsid = line[2]

  inputline = '\t'.join([chrom, pos, rsid])

  query = chrom + ":" + pos + "-" + pos
  recorditer = tb.querys(query)

  # convert iterator to actual array of records
  records = []
  for record in recorditer:
    records.append(record)

  refline = '\t'.join(['.','.','.'])

  if len(records) == 1:
    #do stuff for single match
    record = records[0]
    qchrom = record[0]
    qpos = record[1]
    qrsid = record[2]
    if (rsid == qrsid) or (getRSPOS(record) == pos):
      refline = '\t'.join([qchrom, qpos, qrsid])

  if '.' in refline:
    lowerpos = str(int(pos) - 100)
    upperpos = str(int(pos) + 100)
    query = chrom + ':' + lowerpos + '-' + upperpos
    recorditer = tb.querys(query)

    # convert iterator to actual array of records
    records = []
    for record in recorditer:
      records.append(record)

    # look for matching RSID and RSPOS
    filteredRecords = records
    if len(filteredRecords) > 1 and re.search("^rs", rsid):
      rsidMatched = getIDmatches(rsid, filteredRecords)
      if len(rsidMatched) > 0:
        filteredRecords = rsidMatched
    if len(filteredRecords) > 1:
      filteredRecords = getRSPOSmatches(pos, filteredRecords)
    if len(filteredRecords) > 1 and not re.search("^rs", rsid):
      # find lowest rsid in filteredRecords
      lrsidrecord = []
      for record in filteredRecords:
        qrsidnum = int(record[2][2:])
        if len(lrsidrecord) == 0:
          lrsidnum = sys.maxsize
        else:
          lrsidnum = int(lrsidrecord[2][2:])
        if qrsidnum < lrsidnum:
          lrsidrecord = record
      filteredRecords = [lrsidrecord]

    if len(filteredRecords) == 1:
      record = filteredRecords[0]
      qchrom = record[0]
      qpos = record[1]
      qrsid = record[2]
      refline = '\t'.join([qchrom, qpos, qrsid])
    else: # not a match, so skip
      mapfile.write(inputline + '\t' + refline + '\n')
      continue

  mapfile.write(inputline + '\t' + refline + '\n')


infile.close()
mapfile.close()
