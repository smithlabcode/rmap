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

#include "load_color_space_reads.hpp"
#include "SeedMaker.hpp"
#include "rmap_os.hpp"

#include <cstring>
#include <fstream>

using std::vector;
using std::string;
using std::ptr_fun;
using std::not1;
using std::min;

using std::cerr;
using std::endl;

static const int INPUT_BUFFER_SIZE = 10000;

static char
to_base_symbol(char c) {return (isvalid(c)) ? toupper(c) : 'N';}

static size_t
get_read_word(const string &read) {
  const size_t trunc_to = min(read.length(), SeedMaker::max_seed_part);
  // Need to replace the Ns because otherwise they will destroy the
  // conversion from DNA to integers. Could replace with random
  // bases, but everyone hates non-deterministic programs.
  string s(read.begin(), read.begin() + trunc_to);
  replace(s.begin(), s.end(), 'N', 'A');
  return SeedMaker::make_read_word(s);
}

static void
check_and_add(string &read, const int max_diffs, size_t &read_width, 
	      size_t &prev_base, vector<FastReadCS> &fast_reads, 
	      vector<size_t> &read_words, vector<size_t> &read_index, 
	      size_t &read_count) {
  if (read_width == 0) read_width = read.length() - 1;
  else if (read.length() < read_width)
    throw RMAPException("Incorrect read width");
  else read.erase(read_width + 1);
  
  if (read_count == 0)
    FastReadCS::set_read_width(read_width);

  replace(read.begin(), read.end(), '0', 'A');
  replace(read.begin(), read.end(), '1', 'C');
  replace(read.begin(), read.end(), '2', 'G');
  replace(read.begin(), read.end(), '3', 'T');

  if (read_count == 0) {
    prev_base = base2int(read[0]);
  }
  else if (base2int(read[0]) != prev_base)
    throw RMAPException("inconsistent first base in color space reads");
  
  copy(read.begin() + 1, read.end(), read.begin());
  read.erase(read.begin() + read_width);
  // clean the read
  transform(read.begin(), read.end(), read.begin(), ptr_fun(&to_base_symbol));
  // check for quality
  // const bool good_read = (count(read.begin(), read.end(), 'N') <= max_diffs);
  // if (good_read) {
  fast_reads.push_back(FastReadCS(read));
  read_words.push_back(get_read_word(read));
  read_index.push_back(read_count);
  //}
  ++read_count;
}


void
load_color_space_reads_from_fasta_file(const string &filename, const size_t max_diffs,
				       size_t &read_width, size_t &prev_base,
				       vector<FastReadCS> &fast_reads,
				       vector<size_t> &read_words, vector<size_t> &read_index) {
  std::ifstream in(filename.c_str(), std::ios::binary);
  if (!in) throw RMAPException("cannot open input file " + filename);
  
  char buffer[INPUT_BUFFER_SIZE + 1];
  
  size_t read_count = 0;
  while (!in.eof()) {
    in.getline(buffer, INPUT_BUFFER_SIZE);
    if (in.gcount() > 1) {
      if (in.gcount() == INPUT_BUFFER_SIZE)
	throw RMAPException("Line in " + filename + "\nexceeds max length: " +
			    toa(INPUT_BUFFER_SIZE));
      // correct for dos carriage returns before newlines
      const size_t last_pos = in.gcount() - 2;//strlen(buffer) - 1;
      if (buffer[last_pos] == '\r') buffer[last_pos] = '\0';
      if (buffer[0] != '>') {
	string read(buffer);
	check_and_add(read, max_diffs, read_width, prev_base, fast_reads, 
		      read_words, read_index, read_count);
      }
    }
    in.peek();
  }
  if (fast_reads.empty())
    throw RMAPException("no high-quality reads in file:\"" + filename + "\"");
}

