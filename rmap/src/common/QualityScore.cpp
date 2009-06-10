/*
 *    Part of RMAP software
 *
 *    Copyright (C) 2009 University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Andrew D. Smith
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "QualityScore.hpp"

#include <fstream>
#include "rmap_utils.hpp"

using std::string;

static bool
check_formats(char c, bool &solexa, bool &phred) {
  solexa = solexa && valid_solexa_score(c);
  phred = phred && valid_phred_score(c);
  return (solexa && phred);
}

FASTQScoreType
fastq_score_type(const string filename) {
  static const size_t MAX_LINE_SIZE = 1000;
  std::ifstream f(filename.c_str());
  if (!f)
    throw RMAPException("cannot open input file " + string(filename));
  
  char line[MAX_LINE_SIZE];
  bool solexa = true, phred = true;
  size_t line_count = 0;
  while (f.getline(line, MAX_LINE_SIZE)) {
    if (line_count % 4 == 3) {
      char *c = line;
      while (*c != '\0' && check_formats(*c, solexa, phred)) ++c;
      if (!check_formats(*c, solexa, phred))
	return ((phred) ? FASTQ_Phred : FASTQ_Solexa);
    }
    ++line_count;
  }
  return (phred) ? FASTQ_Phred : FASTQ_Solexa;
}
