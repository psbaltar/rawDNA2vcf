#!/usr/bin/awk -f
#
# AncestryDNA_to_23andme.awk - Converts AncestryDNA format to 23AndMe format
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
# AncestryDNA_to_23andme.awk [inputfile] > [outputfile]



BEGIN {print "# converted to 23andme format using AncestryDNA_to_23andme.awk\n#\n#"}
$1~"^#" {print; next}
$1~"^rsid" {print "# rsid\tchromosome\tposition\tgenotype"; next}
{print $1 "\t" $2 "\t" $3 "\t" $4 $5}
