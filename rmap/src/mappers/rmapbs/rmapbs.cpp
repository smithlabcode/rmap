/*    rmapbs: a program for mapping bisulfite treated Solexa reads
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

#include <string>
#include <vector>
#include <iostream>
#include <iterator>
#include <fstream>
#include <cmath>
#include <tr1/unordered_map>

#include "BisulfiteFastRead.hpp"
#include "BisulfiteFastReadWC.hpp"
#include "BisulfiteFastReadQuality.hpp"
#include "rmap_os.hpp"
#include "SeedMaker.hpp"
#include "MapResult.hpp"
#include "OptionParser.hpp"
#include "load_bisulfite_reads.hpp"

using std::tr1::unordered_multimap;

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::ostream_iterator;
using std::max;
using std::min;
using std::numeric_limits;
using std::pair;
using std::make_pair;

enum { FASTA_FILE, FASTQ_FILE, FASTA_AND_PRB };
enum { RUN_MODE_MISMATCH, RUN_MODE_WILDCARD, RUN_MODE_WEIGHT_MATRIX };

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////
////  THIS STUFF DEALS WITH THE HASHING OF KEYS
////
////
template <typename T> struct SeedHash {
  typedef
  unordered_map<size_t,
		pair<typename vector<T>::const_iterator,
		     typename vector<T>::const_iterator> > type;
};
typedef unordered_multimap<size_t, size_t> SeedHashSorter;


static void
get_read_matches(const size_t the_seed, const vector<size_t> &read_words,
		 SeedHashSorter &sh_sorter) {
  for (size_t i = 0; i < read_words.size(); ++i)
    sh_sorter.insert(SeedHashSorter::value_type(the_seed & read_words[i], i));
}


template <class T> void
sort_by_key(const SeedHashSorter &sh, vector<T> &in) {
  vector<T> tmp(in.size(), in.front());
  size_t j = 0;
  for (SeedHashSorter::const_iterator i(sh.begin()); i != sh.end(); ++i, ++j)
    tmp[j] = in[i->second];
  in.swap(tmp);
}


template <class T> void
sort_by_key(SeedHashSorter &seed_hash, vector<MultiMapResult> &best_maps,
	    vector<size_t> &reads, vector<vector<size_t> > &read_index, 
	    vector<T> &fast_reads) {
  sort_by_key(seed_hash, best_maps);
  sort_by_key(seed_hash, reads);
  sort_by_key(seed_hash, read_index);
  sort_by_key(seed_hash, fast_reads);
  size_t j = 0;
  for (SeedHashSorter::iterator i(seed_hash.begin()); i != seed_hash.end(); ++i, ++j)
    i->second = j;
}


template <class T> void
build_seed_hash(const SeedHashSorter &sh_sorter, const vector<T> &fast_reads,
		typename SeedHash<T>::type &seed_hash) {
  typename vector<T>::const_iterator frb(fast_reads.begin());
  size_t prev_key = 0, prev_idx = 0, curr_idx = 0;
  for (SeedHashSorter::const_iterator shs(sh_sorter.begin()); 
       shs != sh_sorter.end(); ++shs) {
    curr_idx = shs->second;
    if (shs->first != prev_key) {
      seed_hash[prev_key] = make_pair(frb + prev_idx, frb + curr_idx);
      prev_key = shs->first;
      prev_idx = curr_idx;
    }
  }
  seed_hash[prev_key] = make_pair(frb + prev_idx, fast_reads.end());
}


static void
load_seeds(const bool VERBOSE, 
	   const size_t read_width, const size_t n_seeds, 
	   const size_t seed_weight, vector<size_t> &the_seeds) {
  
  SeedMaker::first_last_seeds(min(read_width, SeedMaker::max_seed_part),
			      n_seeds, seed_weight, the_seeds);
  if (VERBOSE) {
    cerr << endl << "SEED STRUCTURES:" << endl;
    for (size_t i = 0; i < the_seeds.size(); ++i)
      cerr << bits2string_masked(rmap_bits::all_ones, the_seeds[i]) << endl;
    cerr << endl;
  }
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////
////  WHERE THE ACTUAL MAPPING HAPPENS
////
////
template <class T> void
map_reads(const string &chrom, const size_t chrom_id,
	  const size_t profile, const size_t read_width, 
	  const size_t max_diffs, const vector<T> &fast_reads, 
	  //const
	  typename SeedHash<T>::type &seed_hash, 
	  const bool strand, vector<MultiMapResult> &best_maps) {
  
  MASK_t bad_bases = rmap_bits::all_ones;
  MASK_t read_word = rmap_bits::all_zeros;
  T fast_read;
  
  const size_t key_diff = read_width - min(read_width, SeedMaker::max_seed_part);
  const size_t chrom_size = chrom.size();
  
  size_t chrom_offset = 0;
  for (chrom_offset = 0; chrom_offset < key_diff; ++chrom_offset)
    fast_read.shift(base2int(chrom[chrom_offset]));
  
  for (; chrom_offset < read_width - 1; ++chrom_offset) {
    const size_t key_base = base2int_bs(chrom[chrom_offset - key_diff]);
    fast_read.shift(base2int(chrom[chrom_offset]));
    SeedMaker::update_bad_bases(key_base, bad_bases);
    SeedMaker::update_read_word(key_base, read_word);
  }
  
  for (; chrom_offset < chrom_size; ++chrom_offset) {
    const size_t key_base = base2int_bs(chrom[chrom_offset - key_diff]);
    fast_read.shift(base2int(chrom[chrom_offset]));
    SeedMaker::update_bad_bases(key_base, bad_bases);
    SeedMaker::update_read_word(key_base, read_word);
    
    if ((bad_bases & profile) == 0) {
      typename SeedHash<T>::type::const_iterator bucket(seed_hash.find(read_word & profile));
      if (bucket != seed_hash.end()) {
	const typename vector<T>::const_iterator limit(bucket->second.second);
	for (typename vector<T>::const_iterator to_test(bucket->second.first); 
	     to_test != limit; ++to_test) {
	  const size_t score = to_test->score(fast_read);
	  if (score <= max_diffs) {
	    const vector<MultiMapResult>::iterator 
	      current(best_maps.begin() + (to_test - fast_reads.begin()));
	    if (score <= current->score)
	      current->add(score, chrom_id, 
			   chrom_offset - read_width + 1, strand);
	  }
	}
      }
    }
  }
}


static void
treat_cpgs(string &chrom) {
  const size_t lim = chrom.length() - 1;
  for (size_t i = 0; i < lim; ++i) {
    if (chrom[i] == 'T' and chrom[i + 1] == 'G')
      chrom[i] = 'C';
  }
}


template <class T> void
iterate_over_seeds(const bool VERBOSE, 
		   const vector<size_t> &the_seeds, const vector<string> &chrom_files,
		   vector<size_t> &ambigs, 
		   vector<string> &chrom_names, vector<size_t> &chrom_sizes,
		   vector<T> &fast_reads,  vector<size_t> &read_words, 
		   vector<vector<size_t> > &read_index,
		   vector<MultiMapResult> &best_maps,
		   const size_t max_mismatches, const size_t read_width) {
  
  if (VERBOSE)
    cerr << "[SCANNING CHROMOSOMES]" << endl;
  
  for (size_t j = 0; j < the_seeds.size() && !fast_reads.empty(); ++j) {
    if (VERBOSE)
      cerr << "[SEED:" << j + 1 << "/" << the_seeds.size() << "] "
	   << "[FORMATTING READS]" << endl;
    
    SeedHashSorter sh_sorter;
    get_read_matches(the_seeds[j], read_words, sh_sorter);
    sort_by_key(sh_sorter, best_maps, read_words, read_index, fast_reads);
    typename SeedHash<T>::type seed_hash;
    build_seed_hash(sh_sorter, fast_reads, seed_hash);
    sh_sorter.clear();
    
    size_t prev_chrom_count = 0;
    for (size_t i = 0; i < chrom_files.size() && !fast_reads.empty(); ++i) {
      
      vector<string> tmp_chrom_names, chroms;
      if (VERBOSE)
	cerr << "[SEED:" << j + 1 << "/" << the_seeds.size() << "] "
	     << "[LOADING CHROM] ";
      read_fasta_file(chrom_files[i].c_str(), tmp_chrom_names, chroms);
      
      if (VERBOSE)
	cerr << "[SCANNING=" << tmp_chrom_names.front() << "] ";
      
      const clock_t start(clock());
      for (size_t k = 0; k < chroms.size(); ++k) {
	if (j == 0) {
	  chrom_sizes.push_back(chroms[k].length());
	  chrom_names.push_back(tmp_chrom_names[k]);
	}
	transform(chroms[k].begin(), chroms[k].end(), chroms[k].begin(), 
		  std::ptr_fun(&toupper));
	treat_cpgs(chroms[k]);
	
	if (VERBOSE)
	  cerr << "[SCANNING POS=" << tmp_chrom_names[k] << "] ";
	map_reads(chroms[k], prev_chrom_count + k, the_seeds[j], read_width, 
		  max_mismatches, fast_reads, seed_hash, true, best_maps);
	chroms[k].clear();
      }
      
      chroms.clear();
      tmp_chrom_names.clear();
      read_fasta_file(chrom_files[i].c_str(), tmp_chrom_names, chroms);
      
      for (size_t k = 0; k < chroms.size(); ++k) {
	revcomp_inplace(chroms[k]);
	transform(chroms[k].begin(), chroms[k].end(), chroms[k].begin(), 
		  std::ptr_fun(&toupper));
	treat_cpgs(chroms[k]);
	if (VERBOSE)
	  cerr << "[SCANNING NEG=" << tmp_chrom_names[k] << "] ";
	map_reads(chroms[k], prev_chrom_count + k, the_seeds[j], read_width, 
		  max_mismatches, fast_reads, seed_hash, false, best_maps);
	chroms[k].clear();
      }
      const clock_t end(clock());
      if (VERBOSE)
	cerr << "[" << static_cast<float>(end - start)/CLOCKS_PER_SEC << " SEC]" << endl;
      prev_chrom_count += chroms.size();
    }
    if (j == 0) {
      if (VERBOSE)
	cerr << "[CLEANING] ";
      eliminate_ambigs(0, best_maps, read_index, read_words, ambigs, fast_reads);
      if (VERBOSE)
	cerr << "[AMBIG=" << ambigs.size() << "] " << endl;
    }
  }
  
  if (VERBOSE)
    cerr << "[FINAL CLEANING] ";
  eliminate_ambigs(max_mismatches, best_maps, read_index, read_words, ambigs, fast_reads);
  if (VERBOSE)
    cerr << "[AMBIG=" << ambigs.size() << "] " << endl;
  
  fast_reads.clear();
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////
////  CODE TO PREPARE AND WRITE THE OUTPUT
////
////
static void
sites_to_regions(const size_t read_len,
		 const vector<string> &chrom, const vector<size_t> &chrom_sizes, 
		 const vector<vector<size_t> > &read_index, 
		 const vector<string> &read_names, vector<MultiMapResult> &bests, 
		 vector<GenomicRegion> &hits) {
  for (size_t i = 0; i < bests.size(); ++i)
    if (!bests[i].empty()) {
      bests[i].sort();
      for (size_t j = 0; j < bests[i].mr.size(); ++j)
	if (j == 0 || bests[i].mr[j - 1] < bests[i].mr[j]) {
	  const size_t chrom_id = bests[i].mr[j].chrom;
	  const size_t start = bests[i].mr[j].strand ? 
	    bests[i].mr[j].site : 
	    chrom_sizes[chrom_id] - bests[i].mr[j].site - read_len;
	  const size_t end = start + read_len;
	  const size_t score = bests[i].score;
	  const char strand = ((bests[i].mr[j].strand) ? '+' : '-');
	  for (size_t k = 0; k < read_index[i].size(); ++k)
	    hits.push_back(GenomicRegion(chrom[chrom_id], start, end,
					 read_names[read_index[i][k]], 
					 score, strand));
	}
    }
}


static void
write_non_uniques(string filename, const vector<size_t> &ambigs,
		  const vector<string> &read_names) {
  std::ofstream out(filename.c_str());
  for (size_t i = 0; i < ambigs.size(); ++i)
    out << read_names[ambigs[i]] << endl;
  out.close();
}


template <class T> void
eliminate_ambigs(const size_t max_mismatches, vector<MultiMapResult> &best_maps, 
		 vector<vector<size_t> > &read_index, vector<size_t> &reads, 
		 vector<size_t> &ambigs, vector<T> &fast_reads) {
  size_t j = 0;
  for (size_t i = 0; i < best_maps.size(); ++i) {
    best_maps[i].collapse();
    if (best_maps[i].ambiguous() && best_maps[i].score <= max_mismatches)
      for (size_t k = 0; k < read_index[i].size(); ++k)
	ambigs.push_back(read_index[i][k]);
    else {
      best_maps[j] = best_maps[i];
      read_index[j].swap(read_index[i]);
      reads[j] = reads[i];
      fast_reads[j] = fast_reads[i];
      ++j;
    }
  }
  best_maps.erase(best_maps.begin() + j, best_maps.end());
  read_index.erase(read_index.begin() + j, read_index.end());
  reads.erase(reads.begin() + j, reads.end());
  fast_reads.erase(fast_reads.begin() + j, fast_reads.end());
}


static void
identify_chromosomes(const bool VERBOSE,
		     const string filenames_file,
		     const string fasta_suffix,
		     const string chrom_file, 
		     vector<string> &chrom_files) {
  if (VERBOSE)
    cerr << "[IDENTIFYING CHROMS] ";
  if (!filenames_file.empty())
    read_filename_file(filenames_file.c_str(), chrom_files);
  else if (isdir(chrom_file.c_str())) 
    read_dir(chrom_file, fasta_suffix, chrom_files);
  else chrom_files.push_back(chrom_file);
  if (VERBOSE) {
    cerr << "[DONE]" << endl 
	 << "chromosome files found (approx size):" << endl;
    for (vector<string>::const_iterator i = chrom_files.begin();
	 i != chrom_files.end(); ++i)
      cerr << *i << " (" << roundf(get_filesize(*i)/1e06) << "Mbp)" << endl;
    cerr << endl;
  }
}


static void
load_read_names(string filename, vector<string> &names) {
  static const size_t INPUT_BUFFER_SIZE = 10000;
  std::ifstream in(filename.c_str(), std::ios::binary);
  if (!in)
    throw RMAPException("cannot open input file " + string(filename));
  while (!in.eof()) {
    char buffer[INPUT_BUFFER_SIZE + 1];
    in.getline(buffer, INPUT_BUFFER_SIZE);
    if (in.gcount() == static_cast<int>(INPUT_BUFFER_SIZE))
      throw RMAPException("Line in " + filename + "\nexceeds max length: " +
                          toa(INPUT_BUFFER_SIZE));
    if (buffer[0] == '>' || buffer[0] == '@') {
      names.push_back(buffer + 1);
      const size_t name_end = names.back().find_first_of(" \t");
      if (name_end != string::npos)
        names.back().erase(names.back().begin() + name_end,
                           names.back().end());
    }
  }
  in.close();
}


static size_t
get_input_mode(const bool VERBOSE, 
	       const string &reads_file, const string &prb_file) {
  size_t INPUT_MODE = FASTA_FILE;
  if (is_fastq(reads_file)) INPUT_MODE = FASTQ_FILE;
  if (!prb_file.empty()) INPUT_MODE = FASTA_AND_PRB;
  if (VERBOSE)
    cerr << "INPUT MODE: "
	 << ((INPUT_MODE == FASTA_FILE) ? 
	     "FASTA" : ((INPUT_MODE == FASTQ_FILE) ? 
			"FASTQ" : "FASTA+PRB")) << endl;
  return INPUT_MODE;
}  


static size_t
get_run_mode(const bool VERBOSE, const size_t INPUT_MODE, 
	     const bool WILDCARD, const bool QUALITY) {
  size_t RUN_MODE = RUN_MODE_MISMATCH;
  if (WILDCARD and QUALITY)
    throw RMAPException("wildcard and quality matching: mutually exclusive");
  if (WILDCARD) {
    if (INPUT_MODE == FASTA_FILE)
      throw RMAPException("quality score information "
			  "required to use wildcards");
    RUN_MODE = RUN_MODE_WILDCARD;
  }
  else if (INPUT_MODE == FASTA_AND_PRB || 
	   (INPUT_MODE == FASTQ_FILE && QUALITY))
    RUN_MODE = RUN_MODE_WEIGHT_MATRIX;
  
  if (VERBOSE)
    cerr << "MATCH MODE: "
	 << ((RUN_MODE == RUN_MODE_MISMATCH) ? 
	     "MISMATCH" : ((RUN_MODE == RUN_MODE_WILDCARD) ? 
			   "WILDCARD" : "WEIGHT-MATRIX")) << endl;
  return RUN_MODE;
}  


static void
load_reads(const bool VERBOSE, const size_t INPUT_MODE, const size_t RUN_MODE,
	   const size_t max_mismatches,
	   const string &reads_file, const string &prb_file,
	   vector<BisulfiteFastRead> &fast_reads, vector<BisulfiteFastReadWC> &fast_reads_wc,
	   vector<BisulfiteFastReadQuality> &fast_reads_q,
	   vector<vector<size_t> > &read_index_in, vector<size_t> &read_words,
	   size_t &read_width) {

  vector<size_t> read_index;
  //////////////////////////////////////////////////////////////
  // LOAD THE READS (AS SEQUENCES OR PROBABILITIES) FROM DISK
  if (VERBOSE) cerr << "[LOADING READ SEQUENCES] ";
  vector<string> reads;
  if (INPUT_MODE == FASTQ_FILE) {
    if (RUN_MODE == RUN_MODE_WILDCARD)
      load_reads_from_fastq_file(reads_file, max_mismatches, read_width,
				 fast_reads_wc, read_words, read_index);
    else if (RUN_MODE == RUN_MODE_WEIGHT_MATRIX)
      load_reads_from_fastq_file(reads_file, max_mismatches, read_width,
				 fast_reads_q, read_words, read_index);
    else
      load_reads_from_fastq_file(reads_file, max_mismatches, read_width,
				 fast_reads, read_words, read_index);
  }
  else if (INPUT_MODE == FASTA_AND_PRB) {
    if (RUN_MODE == RUN_MODE_WILDCARD)
      load_reads_from_prb_file(prb_file, max_mismatches, read_width,
			       fast_reads_wc, read_words, read_index);
    else load_reads_from_prb_file(prb_file, max_mismatches, read_width,
				  fast_reads_q, read_words, read_index);
  }
  else load_reads_from_fasta_file(reads_file, max_mismatches, read_width,
				  fast_reads, read_words, read_index);
  
  for (size_t i = 0; i < read_index.size(); ++i)
    read_index_in.push_back(vector<size_t>(1, read_index[i]));
  if (VERBOSE)
    cerr << "[DONE]" << endl
	 << "TOTAL HQ READS: " << read_index.size() << endl;
}


int 
main(int argc, const char **argv) {
  try {
    
    string chrom_file;
    string filenames_file;
    string outfile;
    string prb_file;
    string ambiguous_file;
    string fasta_suffix = "fa";
    
    size_t n_seeds = 3;
    size_t seed_weight = 11;
    size_t read_width = 0;
    size_t max_mismatches = 10;
    size_t max_mappings = 1;
    
    bool VERBOSE = false;
    bool WILDCARD = false;
    bool QUALITY = false;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse("rmapbs", "The rmapbs mapping tool for Solexa reads"
			   " following bisulfite treatment",
			   "<fast[a/q]-reads-file>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false , outfile);
    opt_parse.add_opt("chrom", 'c', "FASTA file or dir containing chromosome(s)", 
		      false , chrom_file);
    opt_parse.add_opt("suffix", 's', "suffix of FASTA files "
		      "(assumes -c indicates dir)", false , fasta_suffix);
    opt_parse.add_opt("filenames", 'F', "file listing names of "
		      "chromosome files", false , filenames_file);
    opt_parse.add_opt("prb", 'p', "file with quality scores (prb format)", 
		      false, prb_file);
    opt_parse.add_opt("seeds", 'S', "number of seeds", false , n_seeds);
    opt_parse.add_opt("hit", 'h', "weight of hit", false , seed_weight);
    opt_parse.add_opt("width", 'w', "width of reads", false, read_width);
    opt_parse.add_opt("mismatch", 'm', "maximum allowed mismatches", 
		      false , max_mismatches);
    opt_parse.add_opt("ambiguous", 'a', "file to write names of ambiguously "
		      "mapped reads", false , ambiguous_file);
    opt_parse.add_opt("max-map", 'M', "maximum allowed mappings for a read", 
		      false, max_mappings);
    opt_parse.add_opt("wc", 'W', "wildcard matching based on quality scores", 
		      false, WILDCARD);
    opt_parse.add_opt("qual", 'Q', "use quality scores (input must be FASTQ)", 
		      false, QUALITY);
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
    if (chrom_file.empty() && filenames_file.empty()) {
      cerr << "must specify chroms file/dir or filenames file" << endl;
      return EXIT_FAILURE;
    }
    if (chrom_file.empty() && filenames_file.empty()) {
      cerr << "must specify chroms file/dir or filenames file" << endl;
      return EXIT_FAILURE;
    }
    const string reads_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/

    
    //////////////////////////////////////////////////////////////
    //  CHECK HOW QUALITY SCORES ARE USED
    //
    const size_t INPUT_MODE = get_input_mode(VERBOSE, reads_file, prb_file);
    const size_t RUN_MODE = get_run_mode(VERBOSE, INPUT_MODE, WILDCARD, QUALITY);
    
    //////////////////////////////////////////////////////////////
    //  DETERMINE WHICH CHROMOSOMES WILL USED IN MAPPING
    //
    vector<string> chrom_files;
    identify_chromosomes(VERBOSE, filenames_file, fasta_suffix, chrom_file, chrom_files);
    
    
    //////////////////////////////////////////////////////////////
    // OBTAIN THE READS
    // 
    vector<BisulfiteFastRead> fast_reads;
    vector<BisulfiteFastReadWC> fast_reads_wc;
    vector<BisulfiteFastReadQuality> fast_reads_q;
    vector<vector<size_t> > read_index;
    vector<size_t> read_words;
    load_reads(VERBOSE, INPUT_MODE, RUN_MODE, max_mismatches, 
	       reads_file, prb_file, fast_reads, fast_reads_wc,
	       fast_reads_q, read_index, read_words, 
	       read_width);
    
    double max_match_score = max_mismatches*BisulfiteFastReadQuality::get_scaler();
    
    //////////////////////////////////////////////////////////////
    // INITIALIZE THE SEED STRUCTURES
    //
    vector<size_t> the_seeds;
    load_seeds(VERBOSE, read_width, n_seeds, seed_weight, the_seeds);
    
    
    //////////////////////////////////////////////////////////////
    // INITIALIZE THE STRUCTURES THAT HOLD THE RESULTS
    //
    MultiMapResult::init(max_mappings);
    vector<MultiMapResult> 
      best_maps(read_words.size(), MultiMapResult((RUN_MODE == RUN_MODE_WEIGHT_MATRIX) ?
						  max_match_score : max_mismatches));
    
    
    //////////////////////////////////////////////////////////////
    // THIS IS WHERE THE ACTUAL MAPPING HAPPENS
    //
    vector<size_t> chrom_sizes;
    vector<string> chrom_names;
    vector<size_t> ambigs;
    if (RUN_MODE == RUN_MODE_MISMATCH)
      iterate_over_seeds(VERBOSE, the_seeds, chrom_files, ambigs, 
			 chrom_names, chrom_sizes,
			 fast_reads, // USE REGULAR FAST READS
			 read_words, read_index,
			 best_maps, max_mismatches, read_width);
    if (RUN_MODE == RUN_MODE_WILDCARD)
      iterate_over_seeds(VERBOSE, the_seeds, chrom_files, ambigs, 
			 chrom_names, chrom_sizes,
			 fast_reads_wc, // USE FAST READS FOR WILDCARD MATCHING
			 read_words, read_index,
			 best_maps, max_mismatches, read_width);
    if (RUN_MODE == RUN_MODE_WEIGHT_MATRIX)
      iterate_over_seeds(VERBOSE, the_seeds, chrom_files, ambigs, 
			 chrom_names, chrom_sizes,
			 fast_reads_q, // USE FAST READS FOR QUALITY MATCHING
			 read_words, read_index,
			 best_maps, max_match_score, read_width);
    
    
    //////////////////////////////////////////////////////////////
    // LOAD THE NAMES OF READS AGAIN (THEY WILL BE NEEDED)
    //
    vector<string> read_names;
    load_read_names(reads_file, read_names);
    
    
    //////////////////////////////////////////////////////////////
    // IF IDENTITIES OF AMBIGUOUS READS ARE DESIRED, WRITE THEM
    //
    if (!ambiguous_file.empty()) {
      if (VERBOSE)
	cerr << "[WRITING AMBIGS] ";
      write_non_uniques(ambiguous_file, ambigs, read_names);
      if (VERBOSE)
	cerr << "[DONE]" << endl
	     << "TOTAL AMBIGS: " << ambigs.size() << endl;
      ambigs.clear();
    }    
    
    //////////////////////////////////////////////////////////////
    // TRANSFORM THE RESULT STRUCTURES INTO BED FORMAT FOR OUTPUT
    //
    // First make sure the chrom names don't have spaces (cause
    // problems for later processing)
    for (size_t i = 0; i < chrom_names.size(); ++i) {
      const size_t chr_name_end = chrom_names[i].find_first_of(" \t");
      if (chr_name_end != string::npos)
	chrom_names[i].erase(chr_name_end);
    }
    vector<GenomicRegion> hits;
    sites_to_regions(read_width, chrom_names, chrom_sizes, read_index, read_names,
		     best_maps, hits);
    
    
    //////////////////////////////////////////////////////////////
    // CORRECT SCORES IF WEIGHT MATRIX MATCHING USED
    //
    if (RUN_MODE == RUN_MODE_WEIGHT_MATRIX)
      for (size_t i = 0; i < hits.size(); ++i)
	hits[i].set_score(BisulfiteFastReadQuality::value_to_quality(hits[i].get_score()));
    
    
    //////////////////////////////////////////////////////////////
    // OUTPUT THE RESULTS
    //
    std::ostream* out = (!outfile.empty()) ? 
      new std::ofstream(outfile.c_str()) : &cout;
    copy(hits.begin(), hits.end(), ostream_iterator<GenomicRegion>(*out, "\n"));
    if (out != &cout) delete out;
    if (VERBOSE)
      cerr << "[DONE] " << endl
	   << "TOTAL READS MAPPED: " << hits.size() << endl;
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
