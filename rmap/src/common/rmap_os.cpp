/*
 *    Part of RMAP software
 *
 *    Copyright (C) 2008 Cold Spring Harbor Laboratory, 
 *                       University of Southern California and
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

#include "rmap_os.hpp"
#include "rmap_utils.hpp"

#include <fstream>
#include <iostream>
#include <algorithm>
#include <map>
#include <cstring>

using std::string;
using std::vector;
using std::ios_base;


string
strip_path(string full_path) {
  size_t start = full_path.find_last_of('/');
  if (start == string::npos)
    start = 0;
  else ++start;
  return full_path.substr(start);
}

string
strip_path_and_suffix(string full_path) {
  size_t start = full_path.find_last_of('/');
  if (start == string::npos)
    start = 0;
  else ++start;
  size_t end = full_path.find_last_of('.');
  if (end == string::npos)
    end = full_path.length();
  return full_path.substr(start, end - start);
}

bool
isdir(const char *filename) {
  struct stat buffer;
  stat(filename, &buffer);
  return S_ISDIR(buffer.st_mode);
}


////////////////////////////////////////////////////////////////////////
// Stuff dealing with FASTA format sequence files

bool
is_valid_filename(const string name, const string& filename_suffix) {
  const string suffix(name.substr(name.find_last_of(".") + 1));
  return (suffix == filename_suffix);
}



string 
path_join(const string& a, const string& b) {
  return a + "/" + b;
}



void 
read_dir(const string& dirname, string filename_suffix,
	 vector<string> &filenames) {
  DIR *dir;
  if (!(dir = opendir(dirname.c_str())))
    throw RMAPException("could not open directory: " + dirname);
  
  errno = 0;
  struct dirent *ent;
  while ((ent = readdir(dir))) {
    if (is_valid_filename(ent->d_name, filename_suffix))
      filenames.push_back(path_join(dirname, string(ent->d_name)));
    errno = 0;
  }
  // check for some errors
  if (errno)
    throw RMAPException("error reading directory: " + dirname);
  if (filenames.empty())
    throw RMAPException("no valid files found in: " + dirname);
}



bool
is_sequence_line(const char *buffer) {
  return isvalid(buffer[0]);
}


void
parse_score_line(const char *buffer, vector<double> &scr) {
  vector<string> parts;
  rmap::split_whitespace(buffer, parts);
  for (size_t i = 0; i < parts.size(); ++i) {
    scr.push_back(atof(parts[i].c_str()));
  }
}


void
read_fastq_file(const char *filename, vector<string> &names, 
		vector<string> &sequences, vector<vector<double> > &scores) {
  
  static const size_t input_buffer_size = 1000000;
  
  std::ifstream in(filename);
  if (!in) {
    throw RMAPException("cannot open input file " + string(filename));
  }

  // @HANNIBAL_1_FC304RBAAXX_RD1:1:1:601:1775
  // GTTTCTTAAGACCGCCCCTACGGTGCTGGCGCTCGGCCTAATCCCATATATGTCACTTNGTGGATCAAGCA
  // +HANNIBAL_1_FC304RBAAXX_RD1:1:1:601:1775
  // 13 30 28 3 35 14 9 0 2 12 2 18 28 10 -1 -2 -3 8 -2 3 -1 19 8 4 13 15 -2 -0 11 2 -2 4 -1 -4 -2 1 16 15 24 19 14 26 3 16 16 8 8 10 2 3 5 4 6 5 2 3 12 9 -5 0 7 2 -1 1 10 4 3 2 1 0 4
  
  string s, name;
  vector<double> scr;
  bool first_line = true;
  while (!in.eof()) {
    char buffer[input_buffer_size + 1];
    in.getline(buffer, input_buffer_size);
    if (in.gcount() == static_cast<int>(input_buffer_size))
      throw RMAPException("Line in " + name + "\nexceeds max length: " +
			 toa(input_buffer_size));
    // correct for dos carriage returns before newlines
    if (buffer[strlen(buffer) - 1] == '\r')
      buffer[strlen(buffer) - 1] = '\0';
    if (buffer[0] == '@') {
      if (first_line == false && s.length() > 0) {
	names.push_back(name);
	sequences.push_back(s);
	scores.push_back(scr);
      }
      else first_line = false;
      name = buffer;
      name = name.substr(name.find_first_not_of("@ "));
      s = "";
      scr.clear();
    }
    else if (is_sequence_line(buffer))
      s += buffer;
    else {
      vector<double> curr_scr;
      parse_score_line(buffer, curr_scr);
      scr.insert(scr.end(), curr_scr.begin(), curr_scr.end());
    }
  }
  if (!first_line && s.length() > 0) {
    names.push_back(name);
    sequences.push_back(s);
    scores.push_back(scr);
  }
}

struct NotNewline {
  bool operator()(char c) {return c != '\n' && c != '\r';}
};

void
read_fasta_file(const char *filename, vector<string> &names, 
		vector<string> &sequences) {
  
  std::ifstream in(filename, std::ios::binary);
  if (!in) {
    throw RMAPException("cannot open input file " + string(filename));
  }

  static const size_t input_buffer_size = 1000000;
  
  string s, name;
  bool first_line = true;
  while (!in.eof()) {
    char buffer[input_buffer_size + 1];
    in.getline(buffer, input_buffer_size);
    if (in.gcount() == static_cast<int>(input_buffer_size))
      throw RMAPException("Line in " + name + "\nexceeds max length: " +
			  toa(input_buffer_size));
    // correct for dos carriage returns before newlines
    if (buffer[strlen(buffer) - 1] == '\r')
      buffer[strlen(buffer) - 1] = '\0';
    if (buffer[0] == '>') {
      if (first_line == false && s.length() > 0) {
	names.push_back(name);
	sequences.push_back(s);
      }
      else first_line = false;
      name = buffer;
      name = name.substr(name.find_first_not_of("> "));
      s = "";
    }
    else s += buffer;
  }
  if (!first_line && s.length() > 0) {
    names.push_back(name);
    sequences.push_back(s);
  }

//   size_t begin_pos = in.tellg();
//   in.seekg(0, std::ios_base::end);
//   size_t end_pos = in.tellg();
//   in.seekg(0, std::ios_base::beg);
  
//   size_t filesize = end_pos - begin_pos;
//   char *buffer = new char[filesize + 1];
  
//   in.read(buffer, filesize);
//   in.close();
  
//   NotNewline not_newline;
  
//   // find name starts and ends
//   vector<std::pair<size_t, size_t> > name_starts_ends;
//   bool in_name = false;
//   for (size_t i = 0; i < filesize; ++i) {
//     if (buffer[i] == '>' && !in_name) {
//       name_starts_ends.push_back(std::make_pair(i, static_cast<size_t>(0)));
//       in_name = true;
//     }
//     else if (in_name && !not_newline(buffer[i])) {
//       name_starts_ends.back().second = i;
//       in_name = false;
//     }
//   }
//   assert(name_starts_ends.back().second != 0);
  
//   if (name_starts_ends.size() == 0)
//     throw RMAPException("no sequences found in file: " + string(filename));

//   // resize sequences and names
//   names.resize(name_starts_ends.size());
//   sequences.resize(name_starts_ends.size());
  
//   names.front() = string(buffer + name_starts_ends.front().first + 1, // (+1 for '>')
// 			 buffer + name_starts_ends.front().second);
  
//   for (size_t i = 1; i < name_starts_ends.size(); ++i) {
//     rmap::copy_if(buffer + name_starts_ends[i - 1].second,
// 		  buffer + name_starts_ends[i].first,
// 		  back_inserter(sequences[i - 1]), not_newline);
//     names[i] = string(buffer + name_starts_ends[i].first + 1, // (+1 for '>')
// 		      buffer + name_starts_ends[i].second);
//   }
//   rmap::copy_if(buffer + name_starts_ends.back().second, buffer + filesize,
// 		back_inserter(sequences.back()), not_newline);
//   delete[] buffer;

}



void
read_filename_file(const char *filename, vector<string> &filenames) {
  
  static const size_t input_buffer_size = 1000000;
  
  std::ifstream in(filename);
  if (!in)
    throw RMAPException("cannot open input file " + string(filename));
  while (!in.eof()) {
    char buffer[input_buffer_size + 1];
    in.getline(buffer, input_buffer_size);
    if (in.gcount() == static_cast<int>(input_buffer_size))
      throw RMAPException("Line in " + string(filename) +
			  "\nexceeds max length: " +
			  toa(input_buffer_size));
    filenames.push_back(buffer);
    in.peek();
  }
}



size_t 
get_filesize(string filename) {
  std::ifstream f(filename.c_str());
  if (!f.good()) {return 0;}
  size_t begin_pos = f.tellg();
  f.seekg(0, ios_base::end);
  size_t end_pos = f.tellg();
  f.close();
  return end_pos - begin_pos;
}



string
basename(string filename) {
  const string s(filename.substr(0, filename.find_last_of(".")));
  const size_t final_slash = s.find_last_of("/");
  if (final_slash != string::npos)
    return s.substr(final_slash + 1);
  else return s;
}

static void
extract_regions_chrom(const string &chrom_name, const string &filename,
		      const vector<SimpleGenomicRegion> &regions, 
		      vector<string> &sequences) {
  
  std::ifstream in(filename.c_str());
  for (vector<SimpleGenomicRegion>::const_iterator i(regions.begin());
       i != regions.end(); ++i) {
    
    const size_t orig_start_pos = i->get_start();
    const size_t orig_end_pos = i->get_end();
    const size_t orig_region_size = orig_end_pos - orig_start_pos;
    
    const size_t start_pos = orig_start_pos;
    const size_t region_size = orig_region_size;
    assert(start_pos >= 0);
    
    in.seekg(start_pos);
    char buffer[region_size + 1];
    buffer[region_size] = '\0';
    in.read(buffer, region_size);
    sequences.push_back(buffer);
    std::transform(sequences.back().begin(), sequences.back().end(), 
		   sequences.back().begin(), std::ptr_fun(&toupper));
  }
  in.close();
}

static void
extract_regions_chrom(const string &chrom_name, const string &filename,
		      const vector<GenomicRegion> &regions, 
		      vector<string> &sequences) {
  
  std::ifstream in(filename.c_str());
  for (vector<GenomicRegion>::const_iterator i(regions.begin());
       i != regions.end(); ++i) {
    
    const size_t orig_start_pos = i->get_start();
    const size_t orig_end_pos = i->get_end();
    const size_t orig_region_size = orig_end_pos - orig_start_pos;
    
    const size_t start_pos = orig_start_pos;
    const size_t region_size = orig_region_size;
    assert(start_pos >= 0);
    
    in.seekg(start_pos);
    char buffer[region_size + 1];
    buffer[region_size] = '\0';
    in.read(buffer, region_size);
    sequences.push_back(buffer);
    std::transform(sequences.back().begin(), sequences.back().end(), 
		   sequences.back().begin(), std::ptr_fun(&toupper));
    assert(i->get_width() == sequences.back().length());
  }
  in.close();
}



void 
read_dir(const string& dirname, vector<string> &filenames) {
  DIR *dir;
  if (!(dir = opendir(dirname.c_str())))
    throw "could not open directory: " + dirname;
  
  errno = 0;
  struct dirent *ent;
  while ((ent = readdir(dir))) {
    filenames.push_back(path_join(dirname, string(ent->d_name)));
    errno = 0;
  }
  // check for some errors
  if (errno)
    throw "error reading directory: " + dirname;
  if (filenames.empty())
    throw "no valid files found in: " + dirname;
}

void
extract_regions(const string &dirname, 
		const vector<SimpleGenomicRegion> &regions_in, 
		vector<string> &sequences) {
  
  assert(check_sorted(regions_in));

  vector<string> filenames;
  read_dir(dirname, filenames);

  vector<vector<SimpleGenomicRegion> > regions;
  separate_chromosomes(regions_in, regions);
  
  std::map<string, size_t> chrom_regions_map;
  for (size_t i = 0; i < filenames.size(); ++i)
    chrom_regions_map[strip_path(filenames[i])] = i;
  
  for (size_t i = 0; i < regions.size(); ++i) {
    // get the right file
    const string chrom(regions[i].front().get_chrom());
    std::map<string, size_t>::const_iterator f_idx(chrom_regions_map.find(chrom));
    if (f_idx == chrom_regions_map.end())
      throw RMAPException("chrom not found:\t" + chrom);
    extract_regions_chrom(chrom, filenames[f_idx->second], regions[i], sequences);
  }
}

void
extract_regions(const string &dirname, 
		const vector<GenomicRegion> &regions_in, 
		vector<string> &sequences) {

  assert(check_sorted(regions_in));
  
  vector<string> filenames;
  read_dir(dirname, filenames);
  
  vector<vector<GenomicRegion> > regions;
  separate_chromosomes(regions_in, regions);
  
  std::map<string, size_t> chrom_regions_map;
  for (size_t i = 0; i < filenames.size(); ++i)
    chrom_regions_map[strip_path(filenames[i])] = i;
  
  for (size_t i = 0; i < regions.size(); ++i) {
    // get the right file
    const string chrom(regions[i].front().get_chrom());
    std::map<string, size_t>::const_iterator f_idx(chrom_regions_map.find(chrom));
    if (f_idx == chrom_regions_map.end())
      throw RMAPException("chrom not found:\t" + chrom);
    extract_regions_chrom(chrom, filenames[f_idx->second], regions[i], sequences);
  }
}



void
read_prb_file(string filename, vector<vector<vector<double> > > &scores) {
  static const size_t input_buffer_size = 1000000;
  scores.clear();
  std::ifstream in(filename.c_str());
  if (!in)
    throw RMAPException("cannot open input file " + filename);
  string s;
  size_t line_number = 0;
  while (!in.eof()) {
    ++line_number;
    char buffer[input_buffer_size + 1];
    in.getline(buffer, input_buffer_size);
    if (in.gcount() == static_cast<int>(input_buffer_size))
      throw RMAPException("Line in " + filename + 
			    "\nexceeds max length: " +
			    toa(input_buffer_size));
    if (buffer[strlen(buffer) - 1] == '\r')
      buffer[strlen(buffer) - 1] = '\0';
    
    vector<string> parts;
    rmap::split_whitespace(buffer, parts);
    if (parts.size() % rmap::alphabet_size != 0)
      throw RMAPException("Incorrect number of values on line " + toa(line_number) +
			    " in file " + filename);
    scores.push_back(vector<vector<double> >());
    for (size_t i = 0; i < parts.size(); i += rmap::alphabet_size) {
      scores.back().push_back(vector<double>());
      for (size_t j = 0; j < rmap::alphabet_size; ++j)
	scores.back().back().push_back(atof(parts[i + j].c_str()));
    }
    in.peek();
  }
}
