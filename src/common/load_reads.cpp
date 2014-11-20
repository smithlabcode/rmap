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

#include "load_reads.hpp"
#include "SeedMaker.hpp"
#include "smithlab_os.hpp"
#include "QualityScore.hpp"
#include "clip_adaptor_from_reads.hpp"

#include <cstring>
#include <fstream>

using std::vector;
using std::string;
using std::ios_base;
using std::ptr_fun;
using std::not1;
using std::min;
using std::cerr;
using std::endl;


static const size_t MIN_NON_N_IN_READS = 17;

static char
to_base_symbol(char c) {return (isvalid(c)) ? toupper(c) : 'N';}

static size_t
get_read_word(const string &read) {
  const size_t trunc_to = min(read.length(), SeedMaker::max_seed_part);
  // Need to replace the Ns because otherwise they will destroy the
  // conversion from DNA to integers. Could replace with random
  // bases, but everyone hates non-deterministic programs.
  string s(read.begin(), read.begin() + trunc_to);
  for (size_t i = 0; i < s.length(); ++i)
    if (s[i] == 'N')
      s[i] = int2base(rand() % 4);
  return SeedMaker::make_read_word(s);
}

static void
check_and_add(const size_t read_count, string &read, const string &adaptor,
        size_t &read_width, vector<FastRead> &fast_reads,
        vector<size_t> &read_words, vector<unsigned int> &read_index) {

  if (read_width == 0)
    read_width = read.length();
  else if (read.length() < read_width)
    throw SMITHLABException("Incorrect read width:\n" + read + "\n");
  else read.erase(read_width);

  if (fast_reads.empty())
    FastRead::set_read_width(read_width);

  // clean the read
  transform(read.begin(), read.end(), read.begin(), ptr_fun(&to_base_symbol));

  if (!adaptor.empty())
    clip_adaptor_from_read(adaptor, read);

  // check for quality
  const bool good_read =
    (read_width - (count(read.begin(), read.end(), 'N')) >=
     MIN_NON_N_IN_READS);

  if (good_read) {
    fast_reads.push_back(FastRead(read));
    read_words.push_back(get_read_word(read));
    read_index.push_back(read_count);
  }
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

inline bool
is_fastq_name_line(size_t line_count) {
  return ((line_count & 3ul) == 0ul);
}

inline bool
is_fastq_sequence_line(size_t line_count) {
  return ((line_count & 3ul) == 1ul);
}

inline bool
is_fastq_score_name_line(size_t line_count) {
  return ((line_count & 3ul) == 2ul);
}

inline bool
is_fastq_score_line(size_t line_count) {
  return ((line_count & 3ul) == 3ul);
}

inline bool
meet_new_read(std::istream& in) {
  bool flag = false;
  string line;

  in.unget();
  if (in.get() == '\n') {
    getline(in, line);
    size_t size_to_go_back = strlen(line.c_str())+1;
    if (line[0] == '@') {
      size_to_go_back++;
      if (in.get() != '@') {
        flag = true;
      }
    }
    in.seekg(-size_to_go_back, std::ios::cur);
  }
  return flag;
}

static size_t
get_approx_read_size(const string &filename,
                      const size_t n_samples, size_t sample_size) {

  static const size_t megabyte = (1ul << 20);

  const size_t filesize = get_filesize(filename);

  if (sample_size == 0)
    sample_size = std::min(megabyte/10, filesize/n_samples);

  const size_t increment =
    std::floor((filesize - sample_size*n_samples)/
               (n_samples - 1.0)) + sample_size;

  assert(filesize > n_samples && filesize > sample_size &&
         filesize >= n_samples*sample_size);

  std::ifstream in(filename.c_str(), ios_base::binary);

  vector<char> buffer(sample_size);
  double total_lines = 0.0;
  for (size_t i = 0; i < filesize && in.good(); i += increment) {
    in.seekg(i, ios_base::beg);
    in.read(&buffer.front(), sample_size);
    if (in.good())
      total_lines += (0.5 + count(buffer.begin(), buffer.end(), '\n'));
  }
  size_t read_size = (4*n_samples*sample_size)/total_lines;
  return read_size;
}


size_t
divide_and_skip_reads(const string &filename, std::istream& in,
                      const size_t n_reads_to_process, const size_t read_start_idx) {

  // compute the size to read
  size_t read_size = get_approx_read_size(filename, 10, 0);
  size_t size_to_read = n_reads_to_process * read_size;

  // shift
  size_t offset = read_start_idx * read_size;
  in.seekg(offset, std::ios::beg);

  // locate the begining of the next read
  string line;
  size_t bytes_moved = 0;
  if (offset != 0) {
    while ( in.good() && !meet_new_read(in)) {
      getline(in, line);
      bytes_moved += strlen(line.c_str()) + 1;
    }
  }
  size_to_read -= bytes_moved;
  return size_to_read;
}

void
load_reads_from_fastq_file(const string &filename, const size_t read_start_idx,
         const size_t n_reads_to_process, const string &adaptor,
         size_t &read_width, vector<FastRead> &fast_reads,
         vector<size_t> &read_words, vector<unsigned int> &read_index) {
  std::fstream in(filename.c_str());
  if (!in)
    throw SMITHLABException("cannot open input file " + filename);

  in.seekg(0, std::ios::end);
  size_t size_to_read;

  size_to_read = divide_and_skip_reads(filename, in, n_reads_to_process, read_start_idx);

  string line;
  size_t size_readed = 0;
  size_t read_count = 0;
  size_t line_count = 0;
  // start to load reads
  while (size_readed < size_to_read && in.good() && getline(in, line)) {
    if (is_fastq_sequence_line(line_count)) {
      check_and_add(read_count, line, adaptor,read_width, fast_reads,
        read_words, read_index);
      ++read_count;
    }
    size_readed += strlen(line.c_str()) + 1;
    ++line_count;
  }

  // finish reading the last read
  while(in.good() && !meet_new_read(in)) {
    getline(in, line);
    if (is_fastq_sequence_line(line_count)) {
      check_and_add(read_count, line, adaptor,read_width, fast_reads,
        read_words, read_index);
      ++read_count;
    }
    ++line_count;
  }

  if (fast_reads.empty())
    throw SMITHLABException("no good reads between " +
                            toa(read_start_idx) + " and " +
                            toa(read_start_idx + n_reads_to_process) + " in " +
                            filename);
}
