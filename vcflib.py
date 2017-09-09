# vcflib.py - Function library for rawDNA2vcf suite
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

import re
import gzip

def isDIV(record = [], *args):
  VC = getVC(record)
  if VC == 'DIV':
    return True
  else:
    return False

def isSNP(record = [], *args):
  VC = getVC(record)
  if VC == 'SNV':
    return True
  else:
    return False

def getVC(record = [], *args):
  INFO = record[7]
  VC = re.findall("(^.*;)(VC=\w*)(;.*$)", INFO)[0][1].split('=')[1]
  return VC

def getRSPOS(record = [], *args):
  INFO = record[7]
  RSPOS = re.findall("(^.*;)(RSPOS=\d*)(;.*$)", INFO)[0][1].split('=')[1]
  return RSPOS

def getRSPOSmatches(rspos, records = [], *args):
  matches = []
  for record in records:
    myRSPOS = getRSPOS(record)
    if myRSPOS == rspos:
      matches.append(record)
  return matches

def getVCmatches(vc, records = [], *args):
  matches = []
  for record in records:
    myVC = getVC(record)
    if myVC == vc:
      matches.append(record)
  return matches

def getIDmatches(ID, records = [], *args):
  matches = []
  for record in records:
    myID = record[2]
    if myID == ID:
      matches.append(record)
  return matches

def getPOSmatches(pos, records = [], *args):
  matches = []
  for record in records:
    mypos = record[1]
    if mypos == pos:
      matches.append(record)
  return matches

def getVCFheaders(vcffilename):
  headerstring = ''

  if re.search("\.gz$", vcffilename):
    vcffile = gzip.open(vcffilename, 'rt')
  else:
    vcffile = open(vcffilename, 'r')

  for line in vcffile:
    if re.search("^#", line):
      headerstring += line
    else:
      break

  vcffile.close()
  return headerstring


def calcGT(genotype, ref, alt):
  alts = alt.split(',')
  alleles = [ref]+alts

  # figure out what VC we have based on REF and ALT fields
  if len(ref) == 1 and len(alts[0]) == 1:
    vc = 'SNP'
  elif len(ref) > len(alts[0]):
    vc = 'DEL'
  else:
    vc = 'INS'

  gt = []
  for allele in genotype:
    if allele == '-': # no call
      gt.append('.')
    elif vc == 'SNP' and allele in ['A','T','C','G'] and allele in alleles:
      gt.append(str(alleles.index(allele)))
    elif vc == 'DEL' and allele in ['D','I']:
      if allele == 'D':
        gt.append('1')
      else:
        gt.append('0')
    elif vc == 'INS' and allele in ['D','I']:
      if allele == 'I':
        gt.append('1')
      else:
        gt.append('0')
    else:
      gt.append('E')

  if re.search("^\.+$", ''.join(gt)):
    GT = '.'
  else:
    GT = '/'.join(gt)

  return GT


def calcGTsnp(genotype, ref, alt):
  if genotype == '--':
    return '.'

  #double check that we're looking at a SNP
  if (len(ref.split(',')[0]) != 1) and (len(alt.split(',')[0]) != 1):
    return '.'

  alts = alt.split(',')
  if genotype[0] == ref:
    g0 = '0'
  elif genotype[0] in alts:
    g0 = str(alts.index(genotype[0])+1)
  else:
    g0 = '.'

  if genotype[1] == ref:
    g1 = '0'
  elif genotype[1] in alts:
    g1 = str(alts.index(genotype[1])+1)
  else:
    g1 = '.'

  if g0 != '.' and g1 != '.':
    return g0+'/'+g1
  else:
    return '.'

def calcGTdiv(genotype, ref, alt):
  if genotype == '--':
    return '.'

  if len(ref) < len(alt):
    vartype = 'I'
  else:
    vartype = 'D'

  if vartype == 'I':
    if genotype[0] == 'D':
      g0 = '0'
    else:
      g0 = '1'
    if genotype[1] == 'D':
      g1 = '0'
    else:
      g1 = '1'

  if vartype == 'D':
    if genotype[0] == 'I':
      g0 = '0'
    else:
      g0 = '1'
    if genotype[1] == 'I':
      g1 = '0'
    else:
      g1 = '1'

  return g0+'/'+g1

