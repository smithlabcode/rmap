This is the README file for the first release of RMAP version 2. RMAP
is a program for mapping reads from short-read sequencing technology
(such as Solexa/Illumina).


CONTACT INFORMATION:
========================================================================
Andrew D Smith
andrewds@usc.edu
http://www.cmb.usc.edu/people/andrewds


SYSTEM REQUIREMENTS:
========================================================================
The RMAP software will only run on UNIX-like operating systems, and
was developed on Linux systems. The RMAP software requires a fairly
recent C++ compiler (i.e. it must include tr1 headers). RMAP has been
compiled and tested on Linux and OS X operating systems using GCC v4.1
or greater. Also, RMAP will only run on 64-bit machines.


INSTALLATION:
========================================================================
This should be easy: unpack the archive and change into the archive
directory. Then type 'make install'. A 'bin' directory will be created
in the current directory, and it will contain the program
binaries. These can be moved around, and also do not depend on any
dynamic libraries, so they should simply work when executed.


USAGE EXAMPLES:
========================================================================
Each program included in this software package will print a list of
options if executed without any command line arguments. Many of the
programs use similar options (for example, output files are specified
with '-o'). For the most basic usage of rmap to map reads, use the
command:

     rmap -o output.bed -c chroms_dir input_reads.fa

The output will appear in output.bed, and the output is in BED format
(see the UCSC Genome Browser Help documentation for details of this
format). Each line of the file indicates the mapping location for a
read, and the 'score' field in each line indicates the number of
mismatches.


FEATURES:
========================================================================
The second version of RMAP includes several features that were lacking
in the original version:

* QUALITY SCORES: Full use of quality scores, meaning quality scores
  for each base at each position can be used directly in the mapping.

* PAIRED-END READS: Paired-end reads can be mapped, and the procedure
  considers both ends at the time of initial mapping, rather than
  trying to identify mapping positions for each end separately and
  then evaluating whether they have appropriate distance and
  orientation.

* BISULFITE SEQUENCING: RMAP can map reads from bisulfite sequencing
  to an ordinary reference genome. The algorithm can exploit
  unconverted cytosines at non-CpG bases to add mapping specificity.


HISTORY
========================================================================
RMAP was originally developed by Andrew D Smith and Zhenyu Xuan at
Cold Spring Harbor Laboratory (in the lab of Michael Q Zhang). The
current (second) version was written by Andrew D Smith (presently
Assistant Professor in the Molecular and Computational Biology
Section, Department of Biological Sciences at University of Southern
California).

PERFORMANCE
========================================================================
RMAP was able to map 30M reads, each of 32nt, allowing up to 2
mismatches in 11.5 hours on a single core, which is more than 2.6M
reads/hour. Reads were simulated and sampled uniformly from the human
genome. This particular run required roughly 11GB of memory. You might
see slightly different performance, depending on computing hardware
and sequencing application.

LICENSE
========================================================================
The RMAP software for mapping reads from short-read sequencers
Copyright (C) 2009 Andrew D Smith and University of Southern California

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
