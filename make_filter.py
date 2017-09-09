#!/usr/bin/python3
#
# make_filter.py - Creates filters for use with make_map.py
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
# make_filter.py [Sample1] [Sample2] [Sample3] ... [SampleN]

import csv
import re
import sys
from collections import defaultdict

files = sys.argv[1:]
varlist = defaultdict(dict)

for filename in files:
  infilename = filename
  infile = open(infilename, 'r')
  infilereader = csv.reader(infile, delimiter='\t')

  for line in infilereader:
    # skip comment lines
    if re.search("^#", line[0]):
      continue

    rsid = line[0]
    chrom = line[1]
    pos = int(line[2])
    genotype = line[3]

    # for AncestryDNA, convert numbers to canonical:
    # - 23==X, 24==Y, 25==Y, 26==MT
    if chrom == "23":
      chrom = "X"
    elif (chrom == "24") or (chrom == "25"):
      chrom = "Y"
    elif chrom == "26":
      chrom = "MT"

    # replace chromosome id's with numbers for easy sorting
    if chrom == "X":
      chrom = 23
    elif chrom == "Y":
      chrom = 24
    elif chrom == "MT":
      chrom = 25
    else:
      chrom = int(chrom)
  
    # if it's already on the list, skip
    if (chrom in varlist) and (pos in varlist[chrom]):
      continue
    else:
      varlist[chrom][pos] = [str(chrom), str(pos), rsid]

  infile.close()




# output
print("#CHROM\tPOS\tRSID")
for chrom in sorted(varlist):
  # chromosome 0 is showing up in output and don't know why
  # probably some sort of array initialization issue
  if chrom == 0:
    continue

  for pos in sorted(varlist[chrom]):
    variation = varlist[chrom][pos]

    # change the chromosome id's back to strings
    chromnum = variation[0]
    if chromnum == "23":
      chromnum = "X"
    elif chromnum == "24":
      chromnum = "Y"
    elif chromnum == "25":
      chromnum = "MT"
    variation[0] = chromnum

    print('\t'.join(variation))

