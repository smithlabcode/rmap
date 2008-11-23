/*    rmap: a program for mapping Solexa reads
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

#include <string>
#include <vector>
#include <iostream>
#include <iterator>
#include <fstream>
#include <cmath>

#include "FastRead.hpp"
#include "rmap_os.hpp"
#include "SeedMaker.hpp"
#include "MapResult.hpp"
#include "OptionParser.hpp"

#include <tr1/unordered_map>

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::ostream;
using std::ofstream;
using std::ostream_iterator;
using std::max;
using std::min;
using std::numeric_limits;
using std::ifstream;
using std::pair;
using std::make_pair;

using std::tr1::unordered_multimap;

typedef unordered_multimap<size_t, size_t> SeedHash;

void
get_read_matches(const size_t the_seed, const vector<string> &reads,
		 SeedHash &seed_hash) {
  for (size_t i = 0; i < reads.size(); ++i) {
    const size_t truncate_to = min(reads[i].length(), SeedMaker::max_seed_part);
    const string s(reads[i].substr(0, truncate_to));
    if (SeedMaker::valid_seed(the_seed, s)) {
      const size_t the_key = (the_seed & SeedMaker::make_read_word(s));
      seed_hash.insert(SeedHash::value_type(the_key, i));
    }
  }
}


void
map_reads(const string &chrom, const size_t chrom_id,
	  const size_t profile, const size_t read_width, 
	  const size_t max_diffs, const vector<FastRead> &fast_reads,
	  const SeedHash &seed_hash, const bool strand, vector<MapResult> &best_maps) {
  
  MASK_t bad_bases = rmap_bits::all_ones;
  MASK_t read_word = rmap_bits::all_zeros;
  FastRead fast_read;
  
  const size_t key_diff = read_width - min(read_width, SeedMaker::max_seed_part);
  const size_t chrom_size = chrom.size();
  
  size_t chrom_offset = 0;
  for (chrom_offset = 0; chrom_offset < key_diff; ++chrom_offset)
    fast_read.shift(base2int(chrom[chrom_offset]));
  
  for (; chrom_offset < chrom_size; ++chrom_offset) {
    const size_t base = base2int(chrom[chrom_offset]);
    const size_t key_base = base2int(chrom[chrom_offset - key_diff]);
    
    fast_read.shift(base);
    SeedMaker::update_bad_bases(key_base, bad_bases);
    SeedMaker::update_read_word(key_base, read_word);
    
    if ((bad_bases & profile) == 0) {
      std::pair<SeedHash::const_iterator, SeedHash::const_iterator>
	bucket(seed_hash.equal_range((read_word & profile)));
      if (bucket.first != bucket.second) {
	const SeedHash::const_iterator  limit(bucket.second);
	for (SeedHash::const_iterator to_test(bucket.first); to_test != limit; ++to_test) {
	  const size_t score = fast_reads[to_test->second].score(fast_read);
	  if (score <= max_diffs) {
	    const vector<MapResult>::iterator current(best_maps.begin() + to_test->second);
	    if (score <= current->score)
	      current->set(score, chrom_id, chrom_offset - read_width + 1, strand);
	  }
	}
      }
    }
  }
}

static char
to_base_symbol(char c) {return (isvalid(c)) ? toupper(c) : 'N';}

void
clean_reads(const size_t max_diffs, const size_t read_width,
	    const vector<string> &input_read_names,
	    vector<string> &reads, vector<vector<string> > &read_names) {
  vector<pair<string, string> > sorter;
  for (size_t i = 0; i < reads.size(); ++i) {
    if (reads[i].length() != read_width) {
      if (reads[i].length() < read_width)
	throw RMAPException("Incorrect read width");
      else reads[i] = reads[i].substr(0, read_width);
    }
    transform(reads[i].begin(), reads[i].end(), reads[i].begin(),
	      std::ptr_fun(&to_base_symbol));
    if (count(reads[i].begin(), reads[i].end(), 'N') <=
	static_cast<int>(max_diffs))
      sorter.push_back(make_pair(reads[i], input_read_names[i]));
  }
  sort(sorter.begin(), sorter.end());
  reads.clear();
  for (size_t i = 0; i < sorter.size(); ++i) {
    if (i == 0 || sorter[i - 1].first != sorter[i].first) {
      reads.push_back(sorter[i].first);
      read_names.push_back(vector<string>(1, sorter[i].second));
    }
    else read_names.back().push_back(sorter[i].second);
  }
}



void
sites_to_regions(const vector<string> &chrom, const vector<size_t> &chrom_sizes, 
		 const vector<string> &reads, const vector<vector<string> > &read_names, 
		 const vector<MapResult> &best_maps, const size_t max_diffs, 
		 vector<GenomicRegion> &hits) {
  for (size_t i = 0; i < best_maps.size(); ++i)
    if (best_maps[i].unique && best_maps[i].score <= max_diffs) {
      const size_t start = best_maps[i].strand ? best_maps[i].site : 
	chrom_sizes[best_maps[i].chrom] - best_maps[i].site - reads[i].length();
      const size_t end = start + reads[i].length();
      const size_t chrom_id = best_maps[i].chrom;
      const size_t score = best_maps[i].score;
      const char strand = (best_maps[i].strand) ? '+' : '-';
      for (size_t j = 0; j < read_names[i].size(); ++j)
	hits.push_back(GenomicRegion(chrom[chrom_id], start, end, 
				     read_names[i][j], score, strand));
    }
}



void
write_non_uniques(string filename, const vector<pair<string, size_t> > &ambigs) {
  ofstream out(filename.c_str());
  for (size_t i = 0; i < ambigs.size(); ++i)
    out << ambigs[i].first << "\t" << ambigs[i].second << endl;
  out.close();
}



void
eliminate_ambigs(const size_t max_mismatches,
		 vector<MapResult> &best_maps, vector<vector<string> > &read_names,
		 vector<string> &reads, vector<pair<string, size_t> > &ambigs,
		 vector<FastRead> &fast_reads) {
  size_t j = 0;
  for (size_t i = 0; i < best_maps.size(); ++i) {
    if (!best_maps[i].unique && best_maps[i].score <= max_mismatches)
      for (size_t j = 0; j < read_names[i].size(); ++j)
	ambigs.push_back(make_pair(read_names[i][j], best_maps[i].score));
    else {
      best_maps[j] = best_maps[i];
      read_names[j].swap(read_names[i]);
      reads[j].swap(reads[i]);
      fast_reads[j] = fast_reads[i];
      ++j;
    }
  }
  best_maps.erase(best_maps.begin() + j, best_maps.end());
  read_names.erase(read_names.begin() + j, read_names.end());
  reads.erase(reads.begin() + j, reads.end());
  fast_reads.erase(fast_reads.begin() + j, fast_reads.end());
}

int 
main(int argc, const char **argv) {
  try {
    
    string chrom_file;
    string filenames_file;
    string outfile;
    string ambiguous_file;
    string fasta_suffix = "fa";
    
    size_t n_seeds = 0;
    size_t seed_weight = 0;
    size_t read_width = 0;
    size_t max_mismatches = 0;
    
    bool VERBOSE = false;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse("rmap", "The rmap mapping tool for Solexa reads");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false , outfile);
    opt_parse.add_opt("chrom", 'c', "FASTA file or dir containing chromosome(s)", 
		      true , chrom_file);
    opt_parse.add_opt("suffix", 's', "suffix of FASTA files (assumes -c indicates dir)", 
		      false , fasta_suffix);
    opt_parse.add_opt("filenames", 'F', "file listing names of chromosome files", 
		      false , filenames_file);
    opt_parse.add_opt("seeds", 'S', "number of seeds", true , n_seeds);
    opt_parse.add_opt("hit", 'h', "weight of hit", true , seed_weight);
    opt_parse.add_opt("width", 'w', "width of reads", true , read_width);
    opt_parse.add_opt("mismatch", 'm', "maximum allowed mismatches", 
		      true , max_mismatches);
    opt_parse.add_opt("ambiguous", 'a', "file to write names of ambiguously mapped reads", 
		      false , ambiguous_file);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.empty()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string reads_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/
    
    if (read_width == 0) {
      cerr << "ERROR: must specify read width" << endl;
      return EXIT_FAILURE;
    }
    
    // Get the reads
    vector<string> input_read_names, reads;
    read_fasta_file(reads_file.c_str(), input_read_names, reads);
    
    if (VERBOSE)
      cerr << "TOTAL READS:    " << reads.size() << endl;
    
    // prepare the reads
    vector<vector<string> > read_names;
    clean_reads(max_mismatches, read_width, input_read_names, reads, read_names);
    if (VERBOSE)
      cerr << "READS AFTER QC: " << reads.size() << endl;
    
    // convert the reads into word pairs
    FastRead::set_read_width(read_width);
    vector<FastRead> fast_reads;
    for (size_t i = 0; i < reads.size(); ++i)
      fast_reads.push_back(FastRead(reads[i]));

    vector<string> chrom_files;
    if (!filenames_file.empty())
      read_filename_file(filenames_file.c_str(), chrom_files);
    else if (isdir(chrom_file.c_str())) 
      read_dir(chrom_file, fasta_suffix, chrom_files);
    else chrom_files.push_back(chrom_file);
    
    if (VERBOSE) {
      cerr << endl << "chromosome files found (approx size):" << endl;
      for (vector<string>::const_iterator i = chrom_files.begin();
	   i != chrom_files.end(); ++i)
	cerr << *i << " (" << roundf(get_filesize(*i)/1e06) << "Mbp)" << endl;
    }
    
    vector<string> chrom_names;
    for (size_t i = 0; i < chrom_files.size(); ++i)
      chrom_names.push_back(basename(chrom_files[i]));
    
    if (VERBOSE)
      cerr << endl << "scanning chromosomes:" << endl;

    vector<size_t> the_seeds;
    SeedMaker::first_last_seeds(min(read_width, SeedMaker::max_seed_part),
				n_seeds, seed_weight, the_seeds);
    
    if (VERBOSE) {
      cerr << "SEED STRUCTURES" << endl;
      for (size_t i = 0; i < the_seeds.size(); ++i)
	cerr << bits2string_masked(rmap_bits::all_ones, the_seeds[i]) << endl;
      cerr << endl;
    }
    vector<size_t> chrom_sizes(chrom_files.size());
    vector<pair<string, size_t> > ambigs;
    vector<MapResult> best_maps(reads.size(), MapResult(max_mismatches + 1));
    
    for (size_t j = 0; j < the_seeds.size() && !reads.empty(); ++j) {
      for (size_t i = 0; i < chrom_files.size() && !reads.empty(); ++i) {
	if (VERBOSE)
	  cerr << "[SEED:" << j + 1 << "/" 
	       << the_seeds.size() << "] [SEEDING] ";
	
	SeedHash seed_hash;
	get_read_matches(the_seeds[j], reads, seed_hash);
	
	vector<string> chrom_names, chrom;
	if (VERBOSE)
	  cerr << "[LOADING CHROM] ";
	read_fasta_file(chrom_files[i].c_str(), chrom_names, chrom);
	if (j == 0)
	  chrom_sizes[i] = chrom.front().length();
	
	if (VERBOSE)
	  cerr << "[SCANNING=" << chrom_names.front() << "] ";
	map_reads(chrom.front(), i, the_seeds[j], read_width, max_mismatches,
		  fast_reads, seed_hash, true, best_maps);
	revcomp_inplace(chrom.front());
	map_reads(chrom.front(), i, the_seeds[j], read_width, max_mismatches,
		  fast_reads, seed_hash, false, best_maps);
	if (VERBOSE)
	  cerr << "[CLEANING=" << chrom_names.front() << "] ";
	eliminate_ambigs(0, best_maps, read_names, reads, ambigs, fast_reads);
	if (VERBOSE)
	  cerr << "[AMBIG=" << ambigs.size() << "] " << endl;
      }
    }
    if (VERBOSE)
      cerr << endl << "[eliminating ambiguous reads...";
    eliminate_ambigs(max_mismatches, best_maps, read_names, reads, ambigs, fast_reads);
    if (!ambiguous_file.empty())
      write_non_uniques(ambiguous_file, ambigs);
    if (VERBOSE)
      cerr << "done]" << endl;
    
    // transform best matches into BED format
    vector<GenomicRegion> hits;
    sites_to_regions(chrom_names, chrom_sizes, reads, read_names, 
		     best_maps, max_mismatches, hits);
    
    // Output the results
    ostream* out = (!outfile.empty()) ? 
      new ofstream(outfile.c_str()) : &cout;
    copy(hits.begin(), hits.end(), ostream_iterator<GenomicRegion>(*out, "\n"));
    if (out != &cout) delete out;

  }
  catch (const RMAPException &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
