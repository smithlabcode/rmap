/*    rmappe: a program for mapping paired-end Solexa reads
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

#include "FastRead.hpp"
#include "FastReadWC.hpp"
#include "FastReadQuality.hpp"
#include "rmap_os.hpp"
#include "SeedMaker.hpp"
#include "MapResultPE.hpp"
#include "OptionParser.hpp"
#include "load_paired_end_reads.hpp"

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
typedef vector<pair<size_t, unsigned int> > SeedHashSorter;
// typedef unordered_multimap<size_t, size_t> SeedHashSorter;


static void
get_read_matches(const size_t the_seed, const vector<size_t> &read_words,
		 SeedHashSorter &sh_sorter) {
  const size_t lim = read_words.size();
  sh_sorter.resize(read_words.size());
  for (size_t i = 0; i < lim; ++i)
    sh_sorter[i] = make_pair(the_seed & read_words[i], i);
  sort(sh_sorter.begin(), sh_sorter.end());
}


template <class T> void
sort_by_key(const SeedHashSorter &sh, vector<T> &in) {
  vector<T> tmp(in.size(), T());
  size_t j = 0;
  for (SeedHashSorter::const_iterator i(sh.begin()); i != sh.end(); ++i, ++j)
    tmp[j] = in[i->second];
  in.swap(tmp);
}


template <class T> void
sort_by_key(SeedHashSorter &sh_sorter, vector<MultiMapResultPE> &best_maps,
	    vector<size_t> &reads, vector<size_t> &read_index, 
	    vector<T> &fast_reads_left, vector<T> &fast_reads_right) {
  sort_by_key(sh_sorter, best_maps);
  sort_by_key(sh_sorter, reads);
  sort_by_key(sh_sorter, read_index);
  sort_by_key(sh_sorter, fast_reads_left);
  sort_by_key(sh_sorter, fast_reads_right);
  size_t j = 0;
  for (SeedHashSorter::iterator i(sh_sorter.begin()); i != sh_sorter.end(); ++i, ++j)
    i->second = j;
}

template <class T> void
sort_by_key(SeedHashSorter &sh_sorter, vector<MultiMapResult> &best_maps,
	    vector<size_t> &reads, vector<size_t> &read_index, 
	    vector<T> &fast_reads) {
  sort_by_key(sh_sorter, best_maps);
  sort_by_key(sh_sorter, reads);
  sort_by_key(sh_sorter, read_index);
  sort_by_key(sh_sorter, fast_reads);
  size_t j = 0;
  for (SeedHashSorter::iterator i(sh_sorter.begin()); i != sh_sorter.end(); ++i, ++j)
    i->second = j;
}


template <class T> void
build_seed_hash(const SeedHashSorter &sh_sorter, const vector<T> &fast_reads,
		typename SeedHash<T>::type &seed_hash) {
  typename vector<T>::const_iterator frb(fast_reads.begin());
  size_t prev_key = 0, prev_idx = 0, curr_idx = 0;
  seed_hash.clear();
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
load_seeds(const bool VERBOSE, const bool FASTER_MODE,
	   const size_t read_width, const size_t n_seeds, 
	   const size_t seed_weight, vector<size_t> &the_seeds) {
  if (FASTER_MODE) {
    // 0b0000111111001100001111110011000011111100110000111111001111111111
    the_seeds.push_back(1138354285449573375ul);
    // 0b1111110011000011111100110000111111001100001111110011000000111100
    the_seeds.push_back(18213668567193169980ul);
    // 0b0011111100110000111111001100001111110011000011111100110000001111
    the_seeds.push_back(4553417141798292495ul);
    // 0b1100001111110011000011111100110000111111001100001111110000110011
    the_seeds.push_back(14119646626644556851ul);
    // 0b0011000011111100110000111111001100001111110011000011111111110000
    the_seeds.push_back(3529911656661139440ul);
    // 0b1100110000111111001100001111110011000011111100110000111111001100
    the_seeds.push_back(14717535969447448524ul);
    // 0b1111001100001111110011000011111100110000111111001100001111000011
    the_seeds.push_back(17514442047644025795ul);
  }
  else SeedMaker::first_last_seeds(min(read_width, SeedMaker::max_seed_part),
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
map_reads(const string &chrom, const size_t chrom_id, const size_t profile, 
	  const size_t read_width, const size_t max_diffs, 
	  const size_t min_sep, size_t max_sep, 
	  const vector<T> &reads_left, const vector<T> &reads_right,
	  //const
	  typename SeedHash<T>::type &seed_hash, const bool strand, 
	  vector<MultiMapResultPE> &best_maps) {
  
  MASK_t bad_bases = rmap_bits::all_ones;
  MASK_t read_word = rmap_bits::all_zeros;
  T fast_read;
  
  const size_t chrom_size = chrom.size();
  vector<T> possible_lefts(max_sep + 1);
  
  size_t chrom_offset = 0;
  while (chrom_offset < min(chrom_size, read_width - 1)) {
    const size_t key_base = base2int(chrom[chrom_offset++]);
    fast_read.shift(key_base);
    SeedMaker::update_bad_bases(key_base, bad_bases);
    SeedMaker::update_read_word(key_base, read_word);
  }
  
  const string::const_iterator chrom_lim(chrom.end());
  for (string::const_iterator chrom_pos(chrom.begin() + chrom_offset); 
       chrom_pos < chrom_lim; ++chrom_pos) {
    const size_t key_base = *chrom_pos;
    fast_read.shift(key_base);
    SeedMaker::update_bad_bases(key_base, bad_bases);
    SeedMaker::update_read_word(key_base, read_word);
    if ((bad_bases & profile) == 0) {
      possible_lefts[distance(chrom.begin(), chrom_pos) % max_sep] = fast_read;
      typename SeedHash<T>::type::const_iterator bucket(seed_hash.find(read_word & profile));
      if (bucket != seed_hash.end()) {
	pair<typename vector<T>::const_iterator, 
	  typename vector<T>::const_iterator> tmp(bucket->second);
	const typename vector<T>::const_iterator limit(tmp.second);
	for (typename vector<T>::const_iterator to_test(tmp.first); 
	     to_test != limit; ++to_test) {
	  const size_t score = to_test->score(fast_read);
	  if (score <= max_diffs) {
	    chrom_offset = distance(chrom.begin(), chrom_pos);
	    const size_t idx = distance(reads_right.begin(), to_test);
	    const size_t lookback_limit = (chrom_offset > min_sep) ? 
	      chrom_offset - min_sep : chrom_offset;
	    for (size_t i = (chrom_offset > max_sep) ? 
		   chrom_offset - max_sep : 0; i <= lookback_limit; ++i) {
	      const size_t pair_score = score + 
		reads_left[idx].score(possible_lefts[i % max_sep]);
	      const vector<MultiMapResultPE>::iterator current(best_maps.begin() + idx);
	      if (pair_score <= current->score) {
		const size_t left_start = i - read_width + 1;
		const size_t right_start = distance(chrom.begin(), chrom_pos) - read_width + 1;
		current->add(pair_score, chrom_id, left_start, right_start, strand);
	      }
	    }
	  }
	}
      }
    }
  }
}


unsigned char
b2i(char c) {
  switch(c) {
  case 'A' : return 0;
  case 'C' : return 1;
  case 'G' : return 2;
  case 'T' : return 3;
  }
  return 4;
}


unsigned char
comp(char c) {
  switch(c) {
  case 0 : return 3;
  case 1 : return 2;
  case 2 : return 1;
  case 3 : return 0;
  }
  return 4;
}


template <class T> void
iterate_over_seeds(const bool VERBOSE,
		   const vector<size_t> &the_seeds, const vector<string> &chrom_files,
		   vector<size_t> &ambigs, 
		   vector<string> &chrom_names, vector<size_t> &chrom_sizes,
		   vector<T> &fast_reads_left, vector<T> &fast_reads_right,
		   vector<size_t> &read_words, 
		   vector<size_t> &read_index,
		   vector<MultiMapResultPE> &best_maps,
		   const size_t max_mismatches, const size_t read_width,
		   const size_t min_sep, const size_t max_sep) {

  
  if (VERBOSE)
    cerr << "[SCANNING CHROMOSOMES]" << endl;
  
  for (size_t j = 0; j < the_seeds.size() && !fast_reads_left.empty(); ++j) {
    if (VERBOSE)
      cerr << "[SEED:" << j + 1 << "/" << the_seeds.size() << "] "
	   << "[FORMATTING READS]" << endl;
    
    SeedHashSorter sh_sorter;
    get_read_matches(the_seeds[j], read_words, sh_sorter);
    sort_by_key(sh_sorter, best_maps, read_words, read_index, 
		fast_reads_left, fast_reads_right);
    typename SeedHash<T>::type seed_hash;
    build_seed_hash(sh_sorter, fast_reads_right, seed_hash);
    sh_sorter.clear();
    
    size_t prev_chrom_count = 0;
    for (size_t i = 0; i < chrom_files.size() && !fast_reads_left.empty(); ++i) {
      
      if (VERBOSE)
	cerr << "[SEED:" << j + 1 << "/" << the_seeds.size() << "] "
	     << "[LOADING CHROM] ";
      vector<string> tmp_chrom_names, chroms;
      read_fasta_file(chrom_files[i].c_str(), tmp_chrom_names, chroms);
      
      for (size_t k = 0; k < chroms.size(); ++k) {
	if (j == 0) {
	  chrom_sizes.push_back(chroms[k].length());
	  chrom_names.push_back(tmp_chrom_names[k]);
	}
	transform(chroms[k].begin(), chroms[k].end(), chroms[k].begin(), 
		  std::ptr_fun(&toupper));
	transform(chroms[k].begin(), chroms[k].end(), chroms[k].begin(), 
		  std::ptr_fun(&b2i));
	if (VERBOSE)
	  cerr << "[SCANNING=" << tmp_chrom_names[k] << "] ";
	const clock_t start(clock());
	map_reads(chroms[k], prev_chrom_count + k, 
		  the_seeds[j], read_width, max_mismatches,
		  min_sep, max_sep, fast_reads_left, fast_reads_right,
		  seed_hash, true, best_maps);
	transform(chroms[k].begin(), chroms[k].end(), chroms[k].begin(), 
		  std::ptr_fun(&comp));
	reverse(chroms[k].begin(), chroms[k].end());
	map_reads(chroms[k], prev_chrom_count + k, 
		  the_seeds[j], read_width, max_mismatches,
		  min_sep, max_sep, fast_reads_left, fast_reads_right,
		  seed_hash, false, best_maps);
	const clock_t end(clock());
	if (VERBOSE)
	  cerr << "[" << static_cast<float>(end - start)/CLOCKS_PER_SEC << " SEC]" << endl;
	string().swap(chroms[k]);
      }
      prev_chrom_count += chroms.size();
    }
    if (j == 0) {
      if (VERBOSE)
	cerr << "[CLEANING] ";
      eliminate_ambigs(0, best_maps, read_index, read_words, ambigs, 
		       fast_reads_left, fast_reads_right);
      if (VERBOSE)
	cerr << "[AMBIG=" << ambigs.size() << "] " << endl;
    }
  }

  if (VERBOSE)
    cerr << "[FINAL CLEANING] ";
  // Need to multiply max_mismatches by 2 for paired-end!!
  eliminate_ambigs(2*max_mismatches, best_maps, read_index, read_words, ambigs, 
		   fast_reads_left, fast_reads_right);
  if (VERBOSE)
    cerr << "[AMBIG=" << ambigs.size() << "] " << endl;
  
  fast_reads_left.clear();
  fast_reads_right.clear();
}


template <class T> void
map_reads_se(const string &chrom, const size_t chrom_id,
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
    const size_t key_base = base2int(chrom[chrom_offset - key_diff]);
    fast_read.shift(base2int(chrom[chrom_offset]));
    SeedMaker::update_bad_bases(key_base, bad_bases);
    SeedMaker::update_read_word(key_base, read_word);
  }
  
  for (; chrom_offset < chrom_size; ++chrom_offset) {
    const size_t key_base = base2int(chrom[chrom_offset - key_diff]);
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
	  if (score <= max_diffs)
	    (best_maps.begin() + distance(fast_reads.begin(), to_test))->
	      add(score, chrom_id, chrom_offset - read_width + 1, strand);
	}
      }
    }
  }
}


static void
get_best_map(MultiMapResult &a, MultiMapResult &b, 
	     const vector<size_t> &chrom_sizes, 
	     const size_t max_sep, const size_t min_sep,
	     const size_t read_len,
	     MultiMapResultPE &c) {
  if (!a.mr.empty() && !b.mr.empty()) {
    a.sort();
    b.sort();
    for (size_t i = 0; i < a.mr.size(); ++i) {
      if (i == 0 || a.mr[i - 1] < a.mr[i]) {
	const size_t chrom_a = a.mr[i].chrom;
	const size_t left_start = 
	  a.mr[i].strand ? a.mr[i].site : 
	  chrom_sizes[chrom_a] - a.mr[i].site - read_len;
	for (size_t j = 0; j < b.mr.size(); ++j)
	  if (j == 0 || b.mr[j - 1] < b.mr[j]) {
	    const size_t chrom_b = b.mr[j].chrom;
	    if (chrom_a == chrom_b && a.mr[i].strand != b.mr[j].strand) {
	      const size_t right_start = 
		b.mr[j].strand ? b.mr[j].site : 
		chrom_sizes[chrom_b] - b.mr[j].site - read_len;
	      if ((max(left_start, right_start) - 
		   min(left_start, right_start)) <= max_sep &&
		  (max(left_start, right_start) - 
		   min(left_start, right_start)) >= min_sep)
		c.add(a.score + b.score, chrom_a, left_start, 
		      right_start, a.mr[i].strand);
	    }
	  }
      }
    }
  }
}


template <class T> void
iterate_over_seeds_cheat(const bool VERBOSE, const vector<size_t> &the_seeds, 
			 const vector<string> &chrom_files, vector<size_t> &ambigs, 
			 vector<string> &chrom_names, vector<size_t> &chrom_sizes,
			 vector<T> &fast_reads, vector<T> &fast_reads_right,
			 vector<size_t> &read_words, vector<size_t> &read_words_r, 
			 vector<size_t> &read_index, vector<size_t> &read_index_r,
			 vector<MultiMapResultPE> &best_maps_in,
			 const size_t max_mismatches, const size_t read_width,
			 const size_t min_sep, const size_t max_sep) {
  
  const size_t n_reads = read_words.size();
  fast_reads.insert(fast_reads.end(), fast_reads_right.begin(), fast_reads_right.end());
  vector<T>().swap(fast_reads_right);
  read_words.insert(read_words.end(), read_words_r.begin(), read_words_r.end());
  vector<size_t>().swap(read_words_r);
  transform(read_index_r.begin(), read_index_r.end(), read_index_r.begin(),
	    std::bind2nd(std::plus<size_t>(), n_reads));
  read_index.insert(read_index.end(), read_index_r.begin(), read_index_r.end());
  vector<size_t>().swap(read_index_r);
  
  if (VERBOSE)
    cerr << "[SCANNING CHROMOSOMES]" << endl;
  
  vector<MultiMapResult> best_maps(2*n_reads, 2*max_mismatches);
  for (size_t j = 0; j < the_seeds.size() && !fast_reads.empty(); ++j) {
    if (VERBOSE)
      cerr << "[SEED:" << j + 1 << "/" << the_seeds.size() << "] "
	   << "[FORMATTING READS]" << endl;
    
    SeedHashSorter sh_sorter;
    get_read_matches(the_seeds[j], read_words, sh_sorter);
    sort_by_key(sh_sorter, best_maps, read_words, read_index, fast_reads);
    typename SeedHash<T>::type seed_hash;
    build_seed_hash(sh_sorter, fast_reads, seed_hash);
    SeedHashSorter().swap(sh_sorter);
    
    size_t prev_chrom_count = 0;
    for (size_t i = 0; i < chrom_files.size() && !fast_reads.empty(); ++i) {
      
      if (VERBOSE)
	cerr << "[SEED:" << j + 1 << "/" << the_seeds.size() << "] "
	     << "[LOADING CHROM] ";
      vector<string> tmp_chrom_names, chroms;
      read_fasta_file(chrom_files[i].c_str(), tmp_chrom_names, chroms);
      
      for (size_t k = 0; k < chroms.size(); ++k) {
	if (j == 0) {
	  chrom_sizes.push_back(chroms[k].length());
	  chrom_names.push_back(tmp_chrom_names[k]);
	}
	transform(chroms[k].begin(), chroms[k].end(), chroms[k].begin(), 
		  std::ptr_fun(&toupper));
	
	if (VERBOSE)
	  cerr << "[SCANNING=" << tmp_chrom_names[k] << "] ";
	
	const clock_t start(clock());
	map_reads_se(chroms[k], prev_chrom_count + k, 
		     the_seeds[j], read_width, max_mismatches,
		     fast_reads, seed_hash, true, best_maps);
	revcomp_inplace(chroms[k]);
	map_reads_se(chroms[k], prev_chrom_count + k, 
		     the_seeds[j], read_width, max_mismatches,
		     fast_reads, seed_hash, false, best_maps);
	const clock_t end(clock());
	if (VERBOSE)
	  cerr << "[" << static_cast<float>(end - start)/CLOCKS_PER_SEC << " SEC]" << endl;
	chroms[k].clear();
      }
      prev_chrom_count += chroms.size();
    }
    for (size_t i = 0; i < best_maps.size(); ++i)
      best_maps[i].collapse();
  }

  vector<T>().swap(fast_reads);
  vector<size_t>().swap(read_words);
  
  vector<pair<unsigned, unsigned> > sorter;
  for (size_t i = 0; i < read_index.size(); ++i)
    sorter.push_back(make_pair(read_index[i], i));
  sort(sorter.begin(), sorter.end());
  
  vector<MultiMapResult> best_maps2(2*n_reads, max_mismatches);
  for (size_t i = 0; i < sorter.size(); ++i) {
    read_index[i] = sorter[i].first;
    best_maps2[i] = best_maps[sorter[i].second];
  }  
  best_maps.swap(best_maps2);
  vector<MultiMapResult>().swap(best_maps2);
  for (size_t i = 0, j = n_reads; i < n_reads; ++i, ++j) {
    best_maps_in.push_back(MultiMapResultPE(2*max_mismatches + 1));
    get_best_map(best_maps[i], best_maps[j], chrom_sizes, 
		 max_sep, min_sep, read_width, best_maps_in[i]);
  }
  size_t j = 0;
  for (size_t i = 0; i < best_maps_in.size(); ++i) {
    best_maps_in[i].collapse();
    if (best_maps_in[i].ambiguous())
      ambigs.push_back(read_index[i]);
    else {
      best_maps_in[j] = best_maps_in[i];
      read_index[j] = read_index[i];
      ++j;
    }
  }
  
  best_maps_in.erase(best_maps_in.begin() + j, best_maps_in.end());
  read_index.erase(read_index.begin() + j, read_index.end());
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////
////  CODE TO PREPARE AND WRITE THE OUTPUT
////
////
static void
sites_to_regions(const bool VERBOSE, const size_t RUN_MODE, const string outfile, 
		 const size_t read_len, const size_t max_mappings,
		 const vector<string> &chrom, const vector<size_t> &chrom_sizes, 
 		 const vector<size_t> &read_index, vector<string> &read_names, 
		 vector<MultiMapResultPE> &bests) {
  
  static const string LEFT_TAG("_L");
  static const string RIGHT_TAG("_R");
  
  if (VERBOSE)
    cerr << "[WRITING OUTPUT]" << endl;
  std::ostream* out = (!outfile.empty()) ? 
    new std::ofstream(outfile.c_str()) : &cout;
  size_t n_reads_mapped = 0;
  
  for (size_t i = 0; i < bests.size(); ++i)
    if (!bests[i].mr.empty()) {
      bests[i].sort();
      if (bests[i].mr.size() <= max_mappings) {
	for (size_t j = 0; j < bests[i].mr.size(); ++j)
	  if (j == 0 || bests[i].mr[j - 1] < bests[i].mr[j]) {
	    const size_t chrom_id = bests[i].mr[j].chrom;
	    const size_t left_start = 
	      bests[i].mr[j].strand ? bests[i].mr[j].site : 
	      chrom_sizes[chrom_id] - bests[i].mr[j].site2 - read_len;
	    const size_t right_start = 
	      bests[i].mr[j].strand ? bests[i].mr[j].site2 :
	      chrom_sizes[chrom_id] - bests[i].mr[j].site - read_len;
	    const double score = (RUN_MODE == RUN_MODE_WEIGHT_MATRIX) ?
	      FastReadQuality::value_to_quality(bests[i].score) :
	      bests[i].score;
	    const char strand = ((bests[i].mr[j].strand) ? '+' : '-');
	    const vector<string>::iterator 
	      rn(read_names.begin() + read_index[i]);
	    *out << GenomicRegion(chrom[chrom_id], left_start, 
				  left_start + read_len, *rn + LEFT_TAG, 
				  score, strand) << '\n'
		 << GenomicRegion(chrom[chrom_id], right_start,
				  right_start + read_len, *rn + RIGHT_TAG, 
				  score, strand) << '\n';
	    string().swap(*rn);
	    n_reads_mapped++;
	  }
      }
    }
  if (out != &cout) delete out;
  if (VERBOSE)
    cerr << "[DONE] " << endl
	 << "TOTAL READS MAPPED: " << n_reads_mapped << endl;
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
eliminate_ambigs(const size_t max_mismatches, vector<MultiMapResultPE> &best_maps, 
		 vector<size_t> &read_index, vector<size_t> &reads, 
		 vector<size_t> &ambigs, vector<T> &fast_reads_left, 
		 vector<T> &fast_reads_right) {
  size_t j = 0;
  for (size_t i = 0; i < best_maps.size(); ++i) {
    best_maps[i].collapse();
    if (best_maps[i].ambiguous() && best_maps[i].score <= max_mismatches)
      ambigs.push_back(read_index[i]);
    else {
      best_maps[j] = best_maps[i];
      read_index[j] = read_index[i];
      reads[j] = reads[i];
      fast_reads_left[j] = fast_reads_left[i];
      fast_reads_right[j] = fast_reads_right[i];
      ++j;
    }
  }
  best_maps.erase(best_maps.begin() + j, best_maps.end());
  read_index.erase(read_index.begin() + j, read_index.end());
  reads.erase(reads.begin() + j, reads.end());
  fast_reads_left.erase(fast_reads_left.begin() + j, fast_reads_left.end());
  fast_reads_right.erase(fast_reads_right.begin() + j, fast_reads_right.end());
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
 	   vector<FastRead> &fast_reads_left, 
 	   vector<FastRead> &fast_reads_right, 
 	   vector<FastReadWC> &fast_reads_wc_left,
 	   vector<FastReadWC> &fast_reads_wc_right,
 	   vector<FastReadQuality> &fast_reads_q_left,
 	   vector<FastReadQuality> &fast_reads_q_right,
	   vector<size_t> &read_index, vector<size_t> &read_words,
	   size_t &read_width) {

  //////////////////////////////////////////////////////////////
  // LOAD THE READS (AS SEQUENCES OR PROBABILITIES) FROM DISK
  if (VERBOSE) cerr << "[LOADING READ SEQUENCES] ";
  vector<string> reads;
  if (INPUT_MODE == FASTQ_FILE) {
    if (RUN_MODE == RUN_MODE_WILDCARD)
      load_reads_from_fastq_file(reads_file, max_mismatches, read_width,
				 fast_reads_wc_left, fast_reads_wc_right, 
				 read_words, read_index);
    else if (RUN_MODE == RUN_MODE_WEIGHT_MATRIX)
      load_reads_from_fastq_file(reads_file, max_mismatches, read_width,
				 fast_reads_q_left, fast_reads_q_right, 
				 read_words, read_index);
    else
      load_reads_from_fastq_file(reads_file, max_mismatches, read_width,
				 fast_reads_left, fast_reads_right, 
				 read_words, read_index);
  }
  else if (INPUT_MODE == FASTA_AND_PRB) {
    if (RUN_MODE == RUN_MODE_WILDCARD)
      load_reads_from_prb_file(prb_file, max_mismatches, read_width,
			       fast_reads_wc_left, fast_reads_wc_right,
			       read_words, read_index);
    else load_reads_from_prb_file(prb_file, max_mismatches, read_width,
				  fast_reads_q_left, fast_reads_q_right, 
				  read_words, read_index);
  }
  else load_reads_from_fasta_file(reads_file, max_mismatches, read_width,
				  fast_reads_left, fast_reads_right, 
				  read_words, read_index);
  if (VERBOSE)
    cerr << "[DONE]" << endl
	 << "TOTAL HQ READS: " << read_index.size() << endl;
}


static void
load_reads_cheat(const bool VERBOSE, const size_t INPUT_MODE, const size_t RUN_MODE,
		 const size_t max_mismatches,
		 const string &reads_file, const string &prb_file,
		 vector<FastRead> &fast_reads_left, 
		 vector<FastRead> &fast_reads_right, 
		 vector<FastReadWC> &fast_reads_wc_left,
		 vector<FastReadWC> &fast_reads_wc_right,
		 vector<FastReadQuality> &fast_reads_q_left,
		 vector<FastReadQuality> &fast_reads_q_right,
		 vector<size_t> &read_index_l, 
		 vector<size_t> &read_index_r, 
		 vector<size_t> &read_words_l,
		 vector<size_t> &read_words_r,
		 size_t &read_width) {

  //////////////////////////////////////////////////////////////
  // LOAD THE READS (AS SEQUENCES OR PROBABILITIES) FROM DISK
  if (VERBOSE) cerr << "[LOADING READ SEQUENCES] ";
  vector<string> reads;
  if (INPUT_MODE == FASTQ_FILE) {
    if (RUN_MODE == RUN_MODE_WILDCARD)
      load_reads_from_fastq_file_cheat(reads_file, max_mismatches, read_width,
				       fast_reads_wc_left, fast_reads_wc_right, 
				       read_words_l, read_words_r, read_index_l, read_index_r);
    else if (RUN_MODE == RUN_MODE_WEIGHT_MATRIX)
      load_reads_from_fastq_file_cheat(reads_file, max_mismatches, read_width,
				       fast_reads_q_left, fast_reads_q_right, 
				       read_words_l, read_words_r, read_index_l, read_index_r);
    else
      load_reads_from_fastq_file_cheat(reads_file, max_mismatches, read_width,
				       fast_reads_left, fast_reads_right, 
				       read_words_l, read_words_r, read_index_l, read_index_r);
  }
  else if (INPUT_MODE == FASTA_AND_PRB) {
    if (RUN_MODE == RUN_MODE_WILDCARD)
      load_reads_from_prb_file_cheat(prb_file, max_mismatches, read_width,
				     fast_reads_wc_left, fast_reads_wc_right,
				     read_words_l, read_words_r, read_index_l, read_index_r);
    else load_reads_from_prb_file_cheat(prb_file, max_mismatches, read_width,
					fast_reads_q_left, fast_reads_q_right, 
					read_words_l, read_words_r, read_index_l, read_index_r);
  }
  else load_reads_from_fasta_file_cheat(reads_file, max_mismatches, read_width,
					fast_reads_left, fast_reads_right, 
					read_words_l, read_words_r, read_index_l, read_index_r);
  if (VERBOSE)
    cerr << "[DONE]" << endl
	 << "TOTAL HQ READS: " << read_index_l.size() << endl;
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
    
    size_t max_sep = 200;
    size_t min_sep = 0;
    
    bool VERBOSE = false;
    bool WILDCARD = false;
    bool QUALITY = false;
    bool ACCURATE = false;
    bool FASTER_MODE = false;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse("rmappe", "The rmappe mapping tool for "
			   "paired-end Solexa reads",
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
    opt_parse.add_opt("min-sep", '\0', "min separation between ends", 
		      false, min_sep);
    opt_parse.add_opt("max-sep", '\0', "max separation between ends", 
		      false, max_sep);
    opt_parse.add_opt("max-map", 'M', "maximum allowed mappings for a read", 
		      false, max_mappings);
    opt_parse.add_opt("accurate", 'A', "most accurate mode (slow)", 
		      false, ACCURATE);
    opt_parse.add_opt("wc", 'W', "wildcard matching based on quality scores", 
		      false, WILDCARD);
    opt_parse.add_opt("qual", 'Q', "use quality scores (input must be FASTQ)", 
		      false, QUALITY);
    opt_parse.add_opt("faster", 'f', "faster seeds (sensitive to 2 mismatches)", 
		      false, FASTER_MODE);
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
    vector<FastRead> fast_reads_left, fast_reads_right;
    vector<FastReadWC> fast_reads_wc_left, fast_reads_wc_right;
    vector<FastReadQuality> fast_reads_q_left, fast_reads_q_right;
    
    vector<size_t> read_index, read_words;
    vector<size_t> read_index_right;
    vector<size_t> read_words_right;

    
    if (!ACCURATE) {
      load_reads_cheat(VERBOSE, INPUT_MODE, RUN_MODE, max_mismatches, 
		       reads_file, prb_file, 
		       fast_reads_left, fast_reads_right, 
		       fast_reads_wc_left, fast_reads_wc_right,
		       fast_reads_q_left, fast_reads_q_right, 
		       read_index, read_index_right,
		       read_words, read_words_right, read_width);
    }
    else {
      load_reads(VERBOSE, INPUT_MODE, RUN_MODE, max_mismatches, 
		 reads_file, prb_file, 
		 fast_reads_left, fast_reads_right, 
		 fast_reads_wc_left, fast_reads_wc_right,
		 fast_reads_q_left, fast_reads_q_right, 
		 read_index, read_words, 
		 read_width);
    }
    
    
    double max_match_score = max_mismatches*FastReadQuality::get_scaler();
    
    //////////////////////////////////////////////////////////////
    // INITIALIZE THE SEED STRUCTURES
    //
    vector<size_t> the_seeds;
    load_seeds(VERBOSE, FASTER_MODE, read_width, 
	       n_seeds, seed_weight, the_seeds);
    
    
    //////////////////////////////////////////////////////////////
    // INITIALIZE THE STRUCTURES THAT HOLD THE RESULTS
    //
    if (!ACCURATE)
      MultiMapResult::init(max(2*max_mappings, 100ul));
    MultiMapResultPE::init(max_mappings);
    
    vector<MultiMapResultPE> best_maps;
    if (ACCURATE)
      best_maps.resize(read_words.size(), MultiMapResultPE((RUN_MODE == RUN_MODE_WEIGHT_MATRIX) ?
							   2*max_match_score : 2*max_mismatches));
    
    //////////////////////////////////////////////////////////////
    // THIS IS WHERE THE ACTUAL MAPPING HAPPENS
    //
    vector<size_t> chrom_sizes;
    vector<string> chrom_names;
    vector<size_t> ambigs;
    if (ACCURATE) {
      if (RUN_MODE == RUN_MODE_MISMATCH)
	iterate_over_seeds(VERBOSE, the_seeds, chrom_files, ambigs, 
			   chrom_names, chrom_sizes,
			   // USE REGULAR FAST READS
			   fast_reads_left, fast_reads_right, 
			   read_words, read_index,
			   best_maps, max_mismatches, read_width,
			   min_sep, max_sep);
      if (RUN_MODE == RUN_MODE_WILDCARD)
	iterate_over_seeds(VERBOSE, the_seeds, chrom_files, ambigs, 
			   chrom_names, chrom_sizes,
			   // USE FAST READS FOR WILDCARD MATCHING
			   fast_reads_wc_left, fast_reads_wc_right,
			   read_words, read_index,
			   best_maps, max_mismatches, read_width,
			   min_sep, max_sep);
      if (RUN_MODE == RUN_MODE_WEIGHT_MATRIX)
	iterate_over_seeds(VERBOSE, the_seeds, chrom_files, ambigs, 
			   chrom_names, chrom_sizes,
			   // USE FAST READS FOR QUALITY MATCHING
			   fast_reads_q_left, fast_reads_q_right,
			   read_words, read_index,
			   best_maps, max_match_score, read_width,
			   min_sep, max_sep);
    }
    else {
      if (RUN_MODE == RUN_MODE_MISMATCH)
	iterate_over_seeds_cheat(VERBOSE, the_seeds, chrom_files, ambigs, 
				 chrom_names, chrom_sizes,
				 // USE REGULAR FAST READS
				 fast_reads_left, fast_reads_right,
				 read_words, read_words_right, 
				 read_index, read_index_right, 
				 best_maps, max_mismatches, read_width,
				 min_sep, max_sep);
      if (RUN_MODE == RUN_MODE_WILDCARD)
	iterate_over_seeds_cheat(VERBOSE, the_seeds, chrom_files, ambigs, 
				 chrom_names, chrom_sizes,
				 // USE FAST READS FOR WILDCARD MATCHING
				 fast_reads_wc_left, fast_reads_wc_right, 
				 read_words, read_words_right,
				 read_index, read_index_right,
				 best_maps, max_mismatches, read_width,
				 min_sep, max_sep);
      if (RUN_MODE == RUN_MODE_WEIGHT_MATRIX)
	iterate_over_seeds_cheat(VERBOSE, the_seeds, chrom_files, ambigs, 
				 chrom_names, chrom_sizes,
				 // USE FAST READS FOR QUALITY MATCHING
				 fast_reads_q_left, fast_reads_q_right,
				 read_words, read_words_right,
				 read_index, read_index_right,
				 best_maps, max_match_score, read_width,
				 min_sep, max_sep);
    }
    
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
      vector<size_t>().swap(ambigs);
    }
    ambiguous_file.clear();
    
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
    
    sites_to_regions(VERBOSE, RUN_MODE, outfile,
		     read_width, max_mappings, chrom_names, 
		     chrom_sizes, read_index, read_names, best_maps);
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
