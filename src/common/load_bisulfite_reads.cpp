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

#include "load_bisulfite_reads.hpp"
#include "SeedMaker.hpp"
#include "rmap_os.hpp"
#include "QualityScore.hpp"

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
get_read_word_bs(const string &read) {
  const size_t trunc_to = min(read.length(), SeedMaker::max_seed_part);
  // Need to replace the Ns because otherwise they will destroy the
  // conversion from DNA to integers. Could replace with random
  // bases, but everyone hates non-deterministic programs.
  string s(read.begin(), read.begin() + trunc_to);
  replace(s.begin(), s.end(), 'N', 'A');
  replace(s.begin(), s.end(), 'C', 'T');
  return SeedMaker::make_read_word(s);
}

static void
check_and_add(string &read, const int max_diffs,
	      size_t &read_width, vector<BisulfiteFastRead> &fast_reads, 
	      vector<size_t> &read_words, vector<size_t> &read_index, 
	      size_t &read_count) {
  if (read_width == 0) read_width = read.length();
  else if (read.length() < read_width)
    throw RMAPException("Incorrect read width:\n" + read + "\n");
  else read.erase(read_width);

  if (read_count == 0)
    BisulfiteFastRead::set_read_width(read_width);
  
  // clean the read
  transform(read.begin(), read.end(), read.begin(), ptr_fun(&to_base_symbol));
  // check for quality
  const bool good_read = (count(read.begin(), read.end(), 'N') <= max_diffs);
  if (good_read) {
    fast_reads.push_back(BisulfiteFastRead(read));
    read_words.push_back(get_read_word_bs(read));
    read_index.push_back(read_count);
  }
  ++read_count;
}


void
load_reads_from_fasta_file(const string &filename, const size_t max_diffs,
			   size_t &read_width, vector<BisulfiteFastRead> &fast_reads,
			   vector<size_t> &read_words, vector<size_t> &read_index) {
  std::ifstream in(filename.c_str(), std::ios::binary);
  if (!in) throw RMAPException("cannot open input file " + filename);
  
  char buffer[INPUT_BUFFER_SIZE + 1];
  
  size_t read_count = 0, line_count = 0;
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
	if ((line_count & 1ul) == 0)
	  throw RMAPException("empty/multi-line reads or bad FASTA header");
	string read(buffer);
	check_and_add(read, max_diffs, read_width, fast_reads, 
		      read_words, read_index, read_count);
      }
    }
    ++line_count;
    in.peek();
  }
  if (fast_reads.empty())
    throw RMAPException("no high-quality reads in file:\"" + filename + "\"");
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


inline bool
is_fastq_name_line(size_t line_count) {
  return ((line_count % 4) == 0);
}

inline bool
is_fastq_sequence_line(size_t line_count) {
  return ((line_count % 4) == 1);
}

inline bool
is_fastq_score_name_line(size_t line_count) {
  return ((line_count % 4) == 2);
}

inline bool
is_fastq_score_line(size_t line_count) {
  return ((line_count % 4) == 3);
}

void
load_reads_from_fastq_file(const string &filename, const size_t max_diffs,
			   size_t &read_width, vector<BisulfiteFastRead> &fast_reads,
			   vector<size_t> &read_words, vector<size_t> &read_index) {
  std::ifstream in(filename.c_str(), std::ios::binary);
  if (!in) throw RMAPException("cannot open input file " + filename);
  char buffer[INPUT_BUFFER_SIZE + 1];

  size_t read_count = 0, line_count = 0;
  while (!in.eof()) {
    in.getline(buffer, INPUT_BUFFER_SIZE);
    if (in.gcount() > 1) {
      if (in.gcount() == INPUT_BUFFER_SIZE)
	throw RMAPException("Line in " + filename + "\nexceeds max length: " +
			    toa(INPUT_BUFFER_SIZE));
      // correct for dos carriage returns before newlines
      const size_t last_pos = in.gcount() - 2;//strlen(buffer) - 1;
      if (buffer[last_pos] == '\r') buffer[last_pos] = '\0';
      
      // if (is_fastq_name_line(line_count))
      //   if (buffer[0] != '@')
      //     throw RMAPException("invalid FASTQ name line: " + string(buffer));
      if (is_fastq_sequence_line(line_count)) {
	string read(buffer);
	check_and_add(read, max_diffs, read_width, fast_reads, read_words, 
		      read_index, read_count);
      }
      //  if (is_fastq_score_name_line(line_count))
      //    if (buffer[0] != '+')
      //      throw RMAPException("invalid FASTQ score name line: " + string(buffer));
      //  if (is_fastq_score_line(line_count))
      //    ; //!!!!!!!!!!!!!
      ++line_count;
    }
    in.peek();
  }
  if (fast_reads.empty())
    throw RMAPException("no high-quality reads in file:\"" + filename + "\"");
}


static void
check_and_add(const FASTQScoreType score_format, const size_t max_diffs,
	      string &score_line, string &read, size_t &read_width, 
	      vector<BisulfiteFastReadWC> &fast_reads, vector<size_t> &read_words, 
	      vector<size_t> &read_index, size_t &read_count) {
  
  if (read_width == 0) read_width = read.length();
  else if (read.length() < read_width)
    throw RMAPException("Incorrect read width");
  else {
    read.erase(read_width);
    score_line.erase(read_width);
  }
  
  if (read_count == 0)
    BisulfiteFastReadWC::set_read_width(read_width);
  
  transform(read.begin(), read.end(), read.begin(), ptr_fun(&to_base_symbol));
  size_t bad_count = 0;
  vector<vector<double> > scores;
  for (size_t i = 0; i < read_width; ++i) {
    const double error_prob = 
      quality_char_to_error_probability(score_format, score_line[i]);
    const double other_probs = 1.0 - error_prob/(rmap::alphabet_size - 1);
    scores.push_back(vector<double>(rmap::alphabet_size, other_probs));
    scores[i][base2int(read[i])] = error_prob;
    bad_count += (error_prob > 0.5);
  }
  
  const bool good_read = (bad_count <= max_diffs);
  if (good_read) {
    fast_reads.push_back(BisulfiteFastReadWC(scores));
    read_words.push_back(get_read_word_bs(read));
    read_index.push_back(read_count);
  }
  ++read_count;
}

void
load_reads_from_fastq_file(const string &filename, const size_t max_diffs,
			   size_t &read_width, vector<BisulfiteFastReadWC> &fast_reads,
			   vector<size_t> &read_words, vector<size_t> &read_index) {

  FASTQScoreType score_format = fastq_score_type(filename);
  
  std::ifstream in(filename.c_str(), std::ios::binary);
  if (!in) throw RMAPException("cannot open input file " + filename);
  char buffer[INPUT_BUFFER_SIZE + 1];

  size_t read_count = 0, line_count = 0;
  string sequence;
  while (!in.eof()) {
    in.getline(buffer, INPUT_BUFFER_SIZE);
    if (in.gcount() > 1) {
      if (in.gcount() == INPUT_BUFFER_SIZE)
	throw RMAPException("Line in " + filename + "\nexceeds max length: " +
			    toa(INPUT_BUFFER_SIZE));
      // correct for dos carriage returns before newlines
      const size_t last_pos = in.gcount() - 2;//strlen(buffer) - 1;
      if (buffer[last_pos] == '\r') buffer[last_pos] = '\0';
      
      //       if (is_fastq_name_line(line_count))
      // 	if (buffer[0] != '@')
      // 	  throw RMAPException("invalid FASTQ name line: " + string(buffer));
      if (is_fastq_sequence_line(line_count)) {
	sequence = string(buffer);
      }
      //       if (is_fastq_score_name_line(line_count))
      // 	if (buffer[0] != '+')
      // 	  throw RMAPException("invalid FASTQ score name line: " + string(buffer));
      if (is_fastq_score_line(line_count)) {
	string score_line(buffer);
	check_and_add(score_format, max_diffs, score_line, sequence, read_width, 
		      fast_reads, read_words, read_index, read_count);
      }
      ++line_count;
    }
    in.peek();
  }
  if (fast_reads.empty())
    throw RMAPException("no high-quality reads in file:\"" + filename + "\"");
}


static void
check_and_add(const FASTQScoreType score_format, const size_t max_diffs,
	      string &score_line, string &read, size_t &read_width, 
	      vector<BisulfiteFastReadQuality> &fast_reads, vector<size_t> &read_words, 
	      vector<size_t> &read_index, size_t &read_count) {
  
  if (read_width == 0) read_width = read.length();
  else if (read.length() < read_width)
    throw RMAPException("Incorrect read width");
  else {
    read.erase(read_width);
    score_line.erase(read_width);
  }
  if (read_count == 0)
    BisulfiteFastReadQuality::set_read_width(read_width);
  
  // clean the read
  transform(read.begin(), read.end(), read.begin(), ptr_fun(&to_base_symbol));
  
  size_t bad_count = 0;
  vector<vector<double> > scores;
  for (size_t i = 0; i < read_width; ++i) {
    // convert to probability
    const double error_prob = 
      quality_char_to_error_probability(score_format, score_line[i]);
    const double other_probs = 1.0 - error_prob/(rmap::alphabet_size - 1);
    scores.push_back(vector<double>(rmap::alphabet_size, other_probs));
    scores[i][base2int(read[i])] = error_prob;
    bad_count += (error_prob > 0.5);
  }
  
  const bool good_read = (bad_count <= max_diffs);
  
  if (good_read) {
    fast_reads.push_back(BisulfiteFastReadQuality(scores));
    read_words.push_back(get_read_word_bs(read));
    read_index.push_back(read_count);
  }
  ++read_count;
}


void
load_reads_from_fastq_file(const string &filename, const size_t max_diffs,
			   size_t &read_width, vector<BisulfiteFastReadQuality> &fast_reads,
			   vector<size_t> &read_words, vector<size_t> &read_index) {
  FASTQScoreType score_format = fastq_score_type(filename);

  std::ifstream in(filename.c_str(), std::ios::binary);
  if (!in) throw RMAPException("cannot open input file " + filename);
  char buffer[INPUT_BUFFER_SIZE + 1];
  
  size_t read_count = 0, line_count = 0;
  string sequence;

  
  while (!in.eof()) {
    in.getline(buffer, INPUT_BUFFER_SIZE);
    if (in.gcount() > 1) {
      if (in.gcount() == INPUT_BUFFER_SIZE)
	throw RMAPException("Line in " + filename + "\nexceeds max length: " +
			    toa(INPUT_BUFFER_SIZE));
      // correct for dos/mac carriage returns before newlines
      const size_t last_pos = in.gcount() - 2; //strlen(buffer) - 1;
      if (buffer[last_pos] == '\r') buffer[last_pos] = '\0';
      
      //       if (is_fastq_name_line(line_count))
      // 	;
      if (is_fastq_sequence_line(line_count))
	sequence = string(buffer);
      //       if (is_fastq_score_name_line(line_count))
      // 	;
      if (is_fastq_score_line(line_count)) {
	string score_line(buffer);
	check_and_add(score_format, max_diffs, score_line, sequence, 
		      read_width, fast_reads, read_words, read_index, read_count);
      }
      ++line_count;
    }
    in.peek();
  }
  if (fast_reads.empty())
    throw RMAPException("no high-quality reads in file:\"" + filename + "\"");
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

static void
check_and_add(const FASTQScoreType score_format, const size_t max_diffs,
	      const string &score_line, size_t &read_width, 
	      vector<BisulfiteFastReadWC> &fast_reads, vector<size_t> &read_words, 
	      vector<size_t> &read_index, size_t &read_count) {
  
  // parse the score line
  vector<string> parts;
  rmap::split_whitespace(score_line, parts);
  if (parts.size() % rmap::alphabet_size != 0)
    throw RMAPException("bad format:\n" + score_line);
  
  // check the read width
  if (read_width == 0) read_width = parts.size()/rmap::alphabet_size;
  else if (parts.size()/rmap::alphabet_size < read_width)
    throw RMAPException("Incorrect read width");
  else parts.resize(read_width*rmap::alphabet_size);
  if (read_count == 0)
    BisulfiteFastReadWC::set_read_width(read_width);

  // convert to numerical values
  vector<vector<double> > error_probs(read_width, vector<double>(rmap::alphabet_size));
  for (size_t i = 0; i < read_width; ++i)
    for (size_t j = 0; j < rmap::alphabet_size; ++j)
      error_probs[i][j] = atof(parts[i*rmap::alphabet_size + j].c_str());
  
  // convert to probability
  size_t bad_count = 0;
  for (size_t i = 0; i < read_width; ++i) {
    for (size_t j = 0; j < rmap::alphabet_size; ++j)
      error_probs[i][j] = 
	quality_score_to_error_probability(score_format, error_probs[i][j]);
    bad_count += (*min_element(error_probs[i].begin(), error_probs[i].end()) > 0.995);
  }
  
  const bool good_read = (bad_count <= max_diffs);
  if (good_read) {
    fast_reads.push_back(BisulfiteFastReadWC(error_probs));
    string read;
    for (size_t i = 0; i < read_width; ++i)
      read += int2base(min_element(error_probs[i].begin(), 
				   error_probs[i].end()) - error_probs[i].begin());
    read_words.push_back(get_read_word_bs(read));
    read_index.push_back(read_count);
  }
  ++read_count;
}


void
load_reads_from_prb_file(const string &filename, const size_t max_diffs,
			 size_t &read_width, vector<BisulfiteFastReadWC> &fast_reads,
			 vector<size_t> &read_words, vector<size_t> &read_index) {

  FASTQScoreType score_format = FASTQ_Solexa;

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
      const string score_line(buffer);
      check_and_add(score_format, max_diffs, score_line, read_width, 
		    fast_reads, read_words, read_index, read_count);
    }
    in.peek();
  }
  if (fast_reads.empty())
    throw RMAPException("no high-quality reads in file:\"" + filename + "\"");
}
 
 
static void
check_and_add(const FASTQScoreType score_format, const size_t max_diffs,
	      const string &score_line, size_t &read_width, 
	      vector<BisulfiteFastReadQuality> &fast_reads, vector<size_t> &read_words, 
	      vector<size_t> &read_index, size_t &read_count) {
  
  // parse the score line
  vector<string> parts;
  rmap::split_whitespace(score_line, parts);
  if (parts.size() % rmap::alphabet_size != 0)
    throw RMAPException("bad format:\n" + score_line);
  
  // check the read width
  if (read_width == 0) read_width = parts.size()/rmap::alphabet_size;
  else if (parts.size()/rmap::alphabet_size < read_width)
    throw RMAPException("Incorrect read width");
  else parts.resize(read_width*rmap::alphabet_size);
  if (read_count == 0)
    BisulfiteFastReadQuality::set_read_width(read_width);
  
  // convert to numerical values
  vector<vector<double> > error_probs(read_width, vector<double>(rmap::alphabet_size));
  for (size_t i = 0; i < read_width; ++i)
    for (size_t j = 0; j < rmap::alphabet_size; ++j)
      error_probs[i][j] = atof(parts[i*rmap::alphabet_size + j].c_str());
  
  // convert to probability
  size_t bad_count = 0;
  for (size_t i = 0; i < read_width; ++i) {
    for (size_t j = 0; j < rmap::alphabet_size; ++j)
      error_probs[i][j] = 
	quality_score_to_error_probability(score_format, error_probs[i][j]);
    bad_count += (*min_element(error_probs[i].begin(), error_probs[i].end()) > 0.5);
  }
  
  const bool good_read = (bad_count <= max_diffs);
  if (good_read) {
    fast_reads.push_back(BisulfiteFastReadQuality(error_probs));

    string read;
    for (size_t i = 0; i < read_width; ++i)
      read += int2base(min_element(error_probs[i].begin(), 
				   error_probs[i].end()) - error_probs[i].begin());
    read_words.push_back(get_read_word_bs(read));
    read_index.push_back(read_count);
  }
  ++read_count;
}


void
load_reads_from_prb_file(const string &filename, const size_t max_diffs,
			 size_t &read_width, vector<BisulfiteFastReadQuality> &fast_reads,
			 vector<size_t> &read_words, vector<size_t> &read_index) {
  FASTQScoreType score_format = FASTQ_Solexa;
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
      const string score_line(buffer);
      check_and_add(score_format, max_diffs, score_line, read_width, 
		    fast_reads, read_words, read_index, read_count);
    }
    in.peek();
  }
  if (fast_reads.empty())
    throw RMAPException("no high-quality reads in file:\"" + filename + "\"");
}
