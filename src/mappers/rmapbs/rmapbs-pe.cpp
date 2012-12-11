/*    rmapbs: a program for mapping bisulfite treated Solexa reads
 *
 *    Copyright (C) 2012 University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Andrew D. Smith, Song Qiang
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

// #define NDEBUG

#include <string>
#include <vector>
#include <iostream>
#include <iterator>
#include <fstream>
#include <cmath>
#include <unistd.h>
#include <functional>
#include <tr1/unordered_map>

#include "FastRead.hpp"
#include "FastReadWC.hpp"
#include "FastReadQuality.hpp"
#include "smithlab_os.hpp"
#include "SeedMaker.hpp"
#include "MapResult.hpp"
#include "OptionParser.hpp"
#include "load_reads.hpp"
#include "clip_adaptor_from_reads.hpp"
#include "MappedRead.hpp"

using std::tr1::unordered_map;

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
using std::ptr_fun;

enum { FASTA_FILE, FASTQ_FILE, FASTA_AND_PRB };
enum { RUN_MODE_MISMATCH, RUN_MODE_WILDCARD, RUN_MODE_WEIGHT_MATRIX };

static void
bisulfite_treatment(bool AG_WC, size_t &x) 
{
  if (AG_WC)
    x = ((x & 0x5555555555555555ul) | 
	 (((x & 0x5555555555555555ul) << 1) & 
	  (x & 0xAAAAAAAAAAAAAAAAul)));
  else x |= ((x & 0x5555555555555555ul) << 1);
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////
////  THIS STUFF DEALS WITH THE HASHING OF KEYS
////
////
template <typename T> struct SeedHash 
{
  typedef
  unordered_map<size_t,
		pair<typename vector<T>::const_iterator,
		     typename vector<T>::const_iterator> > type;
};
typedef vector<pair<size_t, unsigned int> > SeedHashSorter;


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

inline size_t
bisulfite_treatment_tc(char c) {
  switch(c) {
  case 1 : return 3;
  default : return c;
  }
}

struct regular_score {
  regular_score() 
  {}
  template <class T, class U> size_t score_tc(T a, U b) const {
    return a->score_tc(b);
  }
  template <class T, class U> size_t score_ag(T a, U b) const {
    return a->score_ag(b);
  }
};

struct wildcard_score {
  wildcard_score() {}
  template <class T, class U> size_t score_tc(T a, U b) const {
    return a->score_tc_wild_n(b);
  }
  template <class T, class U> size_t score_ag(T a, U b) const {
    return a->score_ag_wild_n(b);
  }
};


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////
////  WHERE THE ACTUAL MAPPING HAPPENS
////
////
template <class T, class U> void
map_reads_tc(const U &specialized_score,
             const string &chrom, const size_t chrom_id,
             const size_t profile, const size_t read_width, 
             const size_t max_diffs, const vector<T> &fast_reads, 
             //const
             typename SeedHash<T>::type &seed_hash, 
             const bool strand, vector<MultiMapResult> &best_maps) 
{
  
  MASK_t bad_bases = rmap_bits::all_ones;
  MASK_t read_word = rmap_bits::all_zeros;
  T fast_read;
  
  const size_t key_diff = read_width - min(read_width, SeedMaker::max_seed_part);
  const size_t chrom_size = chrom.size();
  
  size_t chrom_offset = 0;
  while(chrom_offset < min(chrom_size, key_diff))
    fast_read.shift(chrom[chrom_offset++]);
  
  while (chrom_offset < min(chrom_size, read_width - 1)) 
    {
      const size_t key_base = chrom[chrom_offset - key_diff];
      fast_read.shift(chrom[chrom_offset++]);
      SeedMaker::update_bad_bases(key_base, bad_bases);
      SeedMaker::update_read_word(key_base, read_word);
    }
  
  string::const_iterator key_pos(chrom.begin() + (chrom_offset - key_diff));
  const string::const_iterator chrom_lim(chrom.end());
  for (string::const_iterator chrom_pos(chrom.begin() + chrom_offset); 
       chrom_pos < chrom_lim; ++chrom_pos, ++key_pos) 
    {
      const size_t key_base = bisulfite_treatment_tc(*key_pos);
      fast_read.shift(*chrom_pos);
      SeedMaker::update_bad_bases(key_base, bad_bases);
      SeedMaker::update_read_word(key_base, read_word);
      if ((bad_bases & profile) == 0) 
        {
	  typename SeedHash<T>::type::const_iterator 
	    bucket(seed_hash.find(read_word & profile));
	  if (bucket != seed_hash.end()) 
            {
	      pair<typename vector<T>::const_iterator, 
		typename vector<T>::const_iterator> tmp(bucket->second);
	      const typename vector<T>::const_iterator limit(tmp.second);
	      for (typename vector<T>::const_iterator to_test(tmp.first); 
		   to_test != limit; ++to_test) 
                {
		  // const size_t score = to_test->score_tc(fast_read);
		  const size_t score = specialized_score.score_tc(to_test, fast_read);
		  if (score <= max_diffs) 
                    {
		      const vector<MultiMapResult>::iterator 
			current(best_maps.begin() + (to_test - fast_reads.begin()));
		      if (score <= current->score)
			current->add(score, chrom_id, distance(chrom.begin(), chrom_pos)
				     - read_width + 1, strand);
                    }
                }
	    }
	}
    }
}


inline size_t
bisulfite_treatment_ag(char c) 
{
  switch(c) 
    {
    case 2 : return 0;
    default : return c;
    }
}


template <class T, class U> void
map_reads_ag(const U &specialized_score,
             const string &chrom, const size_t chrom_id,
             const size_t profile, const size_t read_width, 
             const size_t max_diffs, const vector<T> &fast_reads, 
             //const
             typename SeedHash<T>::type &seed_hash, 
             const bool strand, vector<MultiMapResult> &best_maps) 
{
  
  MASK_t bad_bases = rmap_bits::all_ones;
  MASK_t read_word = rmap_bits::all_zeros;
  T fast_read;
  
  const size_t key_diff = read_width - min(read_width, SeedMaker::max_seed_part);
  const size_t chrom_size = chrom.size();
  
  size_t chrom_offset = 0;
  while(chrom_offset < min(chrom_size, key_diff))
    fast_read.shift(chrom[chrom_offset++]);

  while (chrom_offset < min(chrom_size, read_width - 1)) 
    {
      const size_t key_base = chrom[chrom_offset - key_diff];
      fast_read.shift(chrom[chrom_offset++]);
      SeedMaker::update_bad_bases(key_base, bad_bases);
      SeedMaker::update_read_word(key_base, read_word);
    }
  
  string::const_iterator key_pos(chrom.begin() + (chrom_offset - key_diff));
  const string::const_iterator chrom_lim(chrom.end());
  for (string::const_iterator chrom_pos(chrom.begin() + chrom_offset); 
       chrom_pos < chrom_lim; ++chrom_pos, ++key_pos) 
    {
      const size_t key_base = bisulfite_treatment_ag(*key_pos);
      fast_read.shift(*chrom_pos);
      SeedMaker::update_bad_bases(key_base, bad_bases);
      SeedMaker::update_read_word(key_base, read_word);
      if ((bad_bases & profile) == 0) 
        {
	  typename SeedHash<T>::type::const_iterator 
	    bucket(seed_hash.find(read_word & profile));
	  if (bucket != seed_hash.end()) 
            {
	      pair<typename vector<T>::const_iterator, 
		typename vector<T>::const_iterator> tmp(bucket->second);
	      const typename vector<T>::const_iterator limit(tmp.second);
	      for (typename vector<T>::const_iterator to_test(tmp.first); 
		   to_test != limit; ++to_test) 
                {
		  // const size_t score = to_test->score_ag(fast_read);
		  const size_t score = specialized_score.score_ag(to_test, fast_read);
		  if (score <= max_diffs) 
                    {
		      const vector<MultiMapResult>::iterator 
			current(best_maps.begin() + distance(fast_reads.begin(), to_test));
		      if (score <= current->score)
			current->add(score, chrom_id, distance(chrom.begin(), chrom_pos)
				     - read_width + 1, strand);
                    }
                }
	    }
	}
    }
}


static void
treat_cpgs(const bool AG_WILDCARD, string &chrom) 
{
  const size_t lim = chrom.length() - (!AG_WILDCARD);
  if (AG_WILDCARD) 
    {
      for (size_t i = 1; i < lim; ++i)
	if (chrom[i] == 0 and chrom[i - 1] == 1) chrom[i] = 2;
    }
  else for (size_t i = 0; i < lim; ++i)
	 if (chrom[i] == 3 and chrom[i + 1] == 2) chrom[i] = 1;
}


static void
get_read_matches(const size_t the_seed, const vector<size_t> &read_words,
                 SeedHashSorter &sh_sorter) 
{
  const size_t lim = read_words.size();
  sh_sorter.resize(read_words.size());
  for (size_t i = 0; i < lim; ++i)
    sh_sorter[i] = make_pair(the_seed & read_words[i], i);
  sort(sh_sorter.begin(), sh_sorter.end());
}


template <class T> void
sort_by_key(const SeedHashSorter &sh, vector<T> &in) 
{
  vector<T> tmp(in);
  size_t j = 0;
  for (SeedHashSorter::const_iterator i(sh.begin()); i != sh.end(); ++i, ++j)
    in[j] = tmp[i->second];
}


template <class T> void
sort_by_key(SeedHashSorter &sh_sorter, vector<MultiMapResult> &best_maps,
            vector<size_t> &reads, vector<unsigned int> &read_index, 
            vector<T> &fast_reads) 
{
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
                typename SeedHash<T>::type &seed_hash) 
{
  seed_hash.clear();
  typename vector<T>::const_iterator frb(fast_reads.begin());
  size_t prev_key = 0, prev_idx = 0, curr_idx = 0;
  for (SeedHashSorter::const_iterator shs(sh_sorter.begin()); 
       shs != sh_sorter.end(); ++shs) 
    {
      curr_idx = shs->second;
      if (shs->first != prev_key) 
        {
	  seed_hash[prev_key] = make_pair(frb + prev_idx, frb + curr_idx);
	  prev_key = shs->first;
	  prev_idx = curr_idx;
        }
    }
  seed_hash[prev_key] = make_pair(frb + prev_idx, fast_reads.end());
}


template <class T> static void
resort_reads(const size_t the_seed,
             vector<T> &fast_reads,  vector<size_t> &read_words, 
             vector<unsigned int> &read_index,
             vector<MultiMapResult> &best_maps,
             typename SeedHash<T>::type &seed_hash) 
{
  seed_hash.clear();
  SeedHashSorter sh_sorter;
  get_read_matches(the_seed, read_words, sh_sorter);
  sort_by_key(sh_sorter, best_maps, read_words, read_index, fast_reads);
  build_seed_hash(sh_sorter, fast_reads, seed_hash);
}


unsigned char
b2i(char c) 
{
  switch(c) 
    {
    case 'A' : return 0;
    case 'C' : return 1;
    case 'G' : return 2;
    case 'T' : return 3;
    default : return 4;
    }
}


unsigned char
b2i_rc(char c) 
{
  switch(c) 
    {
    case 'A' : return 3;
    case 'C' : return 2;
    case 'G' : return 1;
    case 'T' : return 0;
    default : return 4;
    }
}


template <class T, class U> void
iterate_over_seeds(const bool VERBOSE, const bool AG_WILDCARD, 
                   const bool ALLOW_METH_BIAS,
                   const U &specialized_score,
                   const vector<size_t> &the_seeds, 
                   const vector<string> &chrom_files,
                   vector<pair<unsigned int, unsigned int> > &ambigs, 
                   vector<string> &chrom_names, 
                   vector<size_t> &chrom_sizes,
                   vector<T> &fast_reads,  vector<size_t> &read_words, 
                   vector<unsigned int> &read_index,
                   vector<MultiMapResult> &best_maps,
                   const size_t max_mismatches, const size_t read_width) 
{
  
  if (VERBOSE)
    cerr << "[SCANNING CHROMOSOMES]" << endl;
  
  for (size_t j = 0; j < the_seeds.size() && !fast_reads.empty(); ++j) 
    {
      if (VERBOSE)
	cerr << "[SEED:" << j + 1 << "/" << the_seeds.size() << "] "
	     << "[FORMATTING READS]" << endl;
    
      typename SeedHash<T>::type seed_hash;
      resort_reads(the_seeds[j], fast_reads, read_words, 
		   read_index, best_maps, seed_hash);
    
      size_t prev_chrom_count = 0;
      for (size_t i = 0; i < chrom_files.size() && !fast_reads.empty(); ++i) 
        {
      
	  vector<string> tmp_chrom_names, chroms;
	  if (VERBOSE)
	    cerr << "[SEED:" << j + 1 << "/" << the_seeds.size() << "] "
		 << "[LOADING CHROM] ";
	  read_fasta_file(chrom_files[i].c_str(), tmp_chrom_names, chroms);
      
	  if (VERBOSE)
	    cerr << "[SCANNING=" << tmp_chrom_names.front() << "] ";
      
	  const clock_t start(clock());
	  for (size_t k = 0; k < chroms.size(); ++k) 
            {
	      if (j == 0) 
                {
		  chrom_sizes.push_back(chroms[k].length());
		  chrom_names.push_back(tmp_chrom_names[k]);
                }
	      transform(chroms[k].begin(), chroms[k].end(), chroms[k].begin(), 
			ptr_fun(&toupper));
	      string tmp_chrom(chroms[k]);
	      transform(chroms[k].begin(), chroms[k].end(), tmp_chrom.begin(),
			ptr_fun(&b2i));
	      if (!ALLOW_METH_BIAS)
		treat_cpgs(AG_WILDCARD, tmp_chrom);
	      if (AG_WILDCARD)
		map_reads_ag(specialized_score, 
			     tmp_chrom, prev_chrom_count + k, the_seeds[j], read_width, 
			     max_mismatches, fast_reads, seed_hash, true, best_maps);
	      else map_reads_tc(specialized_score, 
				tmp_chrom, prev_chrom_count + k, the_seeds[j], read_width, 
				max_mismatches, fast_reads, seed_hash, true, best_maps);
	      string().swap(tmp_chrom);
	
	      std::reverse(chroms[k].begin(), chroms[k].end());
	      transform(chroms[k].begin(), chroms[k].end(), chroms[k].begin(),
			ptr_fun(&b2i_rc));
	      if (!ALLOW_METH_BIAS)
		treat_cpgs(AG_WILDCARD, chroms[k]);
	      if (AG_WILDCARD)
		map_reads_ag(specialized_score, 
			     chroms[k], prev_chrom_count + k, the_seeds[j], read_width, 
			     max_mismatches, fast_reads, seed_hash, false, best_maps);
	      else map_reads_tc(specialized_score, 
				chroms[k], prev_chrom_count + k, the_seeds[j], read_width,
				max_mismatches, fast_reads, seed_hash, false, best_maps);
	      string().swap(chroms[k]);
            }
	  const clock_t end(clock());
	  if (VERBOSE)
	    cerr << "[" << static_cast<float>(end - start)/
	      CLOCKS_PER_SEC << " SEC]" << endl;
	  if (j == 0 && (i + 1) < chrom_files.size()) 
            {
	      if (VERBOSE)
		cerr << "[CLEANING] ";
	      eliminate_ambigs(0, the_seeds[j], best_maps, read_index, 
			       read_words, ambigs, fast_reads, seed_hash);
	      if (VERBOSE)
		cerr << "[AMBIG=" << ambigs.size() << "] " << endl;
            }
	  prev_chrom_count += chroms.size();
        }
      if (j == 0) 
        {
	  if (VERBOSE)
	    cerr << "[CLEANING] ";
	  eliminate_ambigs(1, best_maps, read_index, read_words, ambigs, fast_reads);
	  if (VERBOSE)
	    cerr << "[AMBIG=" << ambigs.size() << "] " << endl;
        } else {
          for (size_t x = 0; x < best_maps.size(); ++x)
              best_maps[x].collapse();          
      }
    }
  if (VERBOSE)
    cerr << "[FINAL CLEANING] ";
  eliminate_ambigs(max_mismatches, best_maps, read_index, 
		   read_words, ambigs, fast_reads);
  if (VERBOSE)
    cerr << "[AMBIG=" << ambigs.size() << "] " << endl;
  vector<T>().swap(fast_reads);
  vector<size_t>().swap(read_words);
}

template <class T> void
eliminate_ambigs(const size_t max_mismatches, const size_t the_seed,
                 vector<MultiMapResult> &best_maps, 
                 vector<unsigned int> &read_index, vector<size_t> &read_words, 
                 vector<pair<unsigned int, unsigned int> > &ambigs, vector<T> &fast_reads,
                 typename SeedHash<T>::type &seed_hash) 
{
  size_t prev_idx = 0, j = 0;
  size_t prev_key = 0;
  seed_hash.clear();
  typename vector<T>::const_iterator frb(fast_reads.begin());
  for (size_t i = 0; i < best_maps.size(); ++i) 
    {
      best_maps[i].collapse();
      if (best_maps[i].ambiguous() && best_maps[i].score <= max_mismatches)
	ambigs.push_back(make_pair(read_index[i], best_maps[i].score));
      else 
        {
	  best_maps[j] = best_maps[i];
	  read_index[j] = read_index[i];
	  fast_reads[j] = fast_reads[i];
	  read_words[j] = read_words[i];
	  const size_t key = (read_words[j] & the_seed);
	  if (j > 0 && key != prev_key) 
            {
	      seed_hash[prev_key] = make_pair(frb + prev_idx, frb + j);
	      prev_idx = j;
            }
	  ++j;
	  prev_key = key;
        }
    }
  seed_hash[prev_key] = make_pair(frb + prev_idx, frb + j);
  best_maps.erase(best_maps.begin() + j, best_maps.end());
  // This below should work but doesn't... Is there a bug elsewhere?
  // vector<MultiMapResult>(best_maps).swap(best_maps);
  read_index.erase(read_index.begin() + j, read_index.end());
  vector<unsigned int>(read_index).swap(read_index);
  read_words.erase(read_words.begin() + j, read_words.end());
  vector<size_t>(read_words).swap(read_words);
  fast_reads.erase(fast_reads.begin() + j, fast_reads.end());
}


template <class T> void
eliminate_ambigs(const size_t max_mismatches, vector<MultiMapResult> &best_maps,
                 vector<unsigned int> &read_index, vector<size_t> &reads, 
                 vector<pair<unsigned int, unsigned int> > &ambigs, 
		 vector<T> &fast_reads) {
  size_t j = 0;
  for (size_t i = 0; i < best_maps.size(); ++i) {
    best_maps[i].collapse();
    if (best_maps[i].ambiguous() && best_maps[i].score <= max_mismatches)
      ambigs.push_back(make_pair(read_index[i], best_maps[i].score));
    else {
      best_maps[j] = best_maps[i];
      read_index[j] = read_index[i];
      reads[j] = reads[i];
      fast_reads[j] = fast_reads[i];
      ++j;
    }
  }
  best_maps.erase(best_maps.begin() + j, best_maps.end());
  vector<MultiMapResult>(best_maps).swap(best_maps);
  read_index.erase(read_index.begin() + j, read_index.end());
  vector<unsigned int>(read_index).swap(read_index);
  reads.erase(reads.begin() + j, reads.end());
  vector<size_t>(reads).swap(reads);
  fast_reads.erase(fast_reads.begin() + j, fast_reads.end());
  vector<T>(fast_reads).swap(fast_reads);
}


static void
identify_chromosomes(const bool VERBOSE,
                     const string filenames_file,
                     const string fasta_suffix,
                     const string chrom_file, 
                     vector<string> &chrom_files) 
{
  if (VERBOSE)
    cerr << "[IDENTIFYING CHROMS] ";
  if (!filenames_file.empty())
    read_filename_file(filenames_file.c_str(), chrom_files);
  else if (isdir(chrom_file.c_str())) 
    read_dir(chrom_file, fasta_suffix, chrom_files);
  else chrom_files.push_back(chrom_file);
  if (VERBOSE) 
    {
      cerr << "[DONE]" << endl 
	   << "chromosome files found (approx size):" << endl;
      for (vector<string>::const_iterator i = chrom_files.begin();
	   i != chrom_files.end(); ++i)
	cerr << *i << " (" << roundf(get_filesize(*i)/1e06) << "Mbp)" << endl;
      cerr << endl;
    }
}

static size_t
get_input_mode(const bool VERBOSE, const string &reads_file) 
{
  size_t INPUT_MODE = FASTA_FILE;
  if (is_fastq(reads_file)) INPUT_MODE = FASTQ_FILE;
  if (VERBOSE)
    cerr << "INPUT MODE: "
	 << ((INPUT_MODE == FASTA_FILE) ? 
	     "FASTA" : "FASTQ") << endl;
  return INPUT_MODE;
}


static size_t
get_run_mode(const bool VERBOSE, const size_t INPUT_MODE, 
             const bool WILDCARD, const bool QUALITY) 
{
  size_t RUN_MODE = RUN_MODE_MISMATCH;
  if (WILDCARD and QUALITY)
    throw SMITHLABException("wildcard and quality matching: mutually exclusive");
  if (WILDCARD) 
    {
      if (INPUT_MODE == FASTA_FILE)
	throw SMITHLABException("quality score information "
				"required to use wildcards");
      RUN_MODE = RUN_MODE_WILDCARD;
    }
  else if (INPUT_MODE == FASTQ_FILE && QUALITY)
    RUN_MODE = RUN_MODE_WEIGHT_MATRIX;
  if (VERBOSE)
    cerr << "MATCH MODE: "
	 << ((RUN_MODE == RUN_MODE_MISMATCH) ? 
	     "MISMATCH" : ((RUN_MODE == RUN_MODE_WILDCARD) ? 
			   "WILDCARD" : "WEIGHT-MATRIX")) << endl;
  return RUN_MODE;
}  


static void
load_reads(const bool VERBOSE, const size_t INPUT_MODE, 
           const size_t RUN_MODE, const bool AG_WILDCARD,
           const size_t max_mismatches, const string &adaptor,
           const string &reads_file,
	   const size_t read_start_index, const size_t n_reads_to_process,
	   vector<FastRead> &fast_reads, vector<FastReadWC> &fast_reads_wc,
           vector<FastReadQuality> &fast_reads_q,
           vector<unsigned int> &read_index, vector<size_t> &read_words,
           size_t &read_width) 
{

  //////////////////////////////////////////////////////////////
  // LOAD THE READS (AS SEQUENCES OR PROBABILITIES) FROM DISK
  if (VERBOSE) cerr << "[LOADING READ SEQUENCES] ";
  vector<string> reads;
  if (INPUT_MODE == FASTQ_FILE) 
    {
      if (RUN_MODE == RUN_MODE_WILDCARD)
	load_reads_from_fastq_file(reads_file, read_start_index, 
				   n_reads_to_process,
				   adaptor, max_mismatches, read_width,
				   fast_reads_wc, read_words, read_index);
      else if (RUN_MODE == RUN_MODE_WEIGHT_MATRIX)
	load_reads_from_fastq_file(reads_file, read_start_index, 
				   n_reads_to_process, adaptor, max_mismatches, read_width,
				   fast_reads_q, read_words, read_index);
      else
	load_reads_from_fastq_file(reads_file, read_start_index, 
				   n_reads_to_process, adaptor, max_mismatches, read_width,
				   fast_reads, read_words, read_index);
    }
  else load_reads_from_fasta_file(reads_file, read_start_index, 
				  n_reads_to_process, adaptor, max_mismatches, read_width,
				  fast_reads, read_words, read_index);
  for (size_t i = 0; i < fast_reads_wc.size(); ++i)
    fast_reads_wc[i].bisulfite_treatment(AG_WILDCARD);
  for (size_t i = 0; i < fast_reads_q.size(); ++i)
    fast_reads_q[i].bisulfite_treatment(AG_WILDCARD);
  for (size_t i = 0; i < read_words.size(); ++i)
    bisulfite_treatment(AG_WILDCARD, read_words[i]);
  if (VERBOSE)
    cerr << "[DONE]" << endl
	 << "TOTAL HQ READS: " << read_index.size() << endl
	 << "READ WIDTH: " << read_width << endl;
}


template<class T>
struct indexed_best_less 
{
  bool operator()(const pair<unsigned int, T> &a,
		  const pair<unsigned int, T> &b) const 
  {
    return a.first < b.first;
  }
};


// static void
// invert_bests_list(vector<unsigned int> &read_index, 
//                   vector<MultiMapResult> &bests) 
// {
//   vector<pair<unsigned int, MultiMapResult> > sorter;
//   for (size_t i = 0; i < bests.size(); ++i)
//     sorter.push_back(make_pair(read_index[i], bests[i]));
//   sort(sorter.begin(), sorter.end(), indexed_best_less<MultiMapResult>());
//   for (size_t i = 0; i < sorter.size(); ++i) 
//     {
//       read_index[i] = sorter[i].first;
//       bests[i] = sorter[i].second;
//     }
// }

static void
clean_and_invert_bests_list(vector<unsigned int> &read_index,
                            vector<MultiMapResult> &best_maps,
                            vector<vector<MultiMapResult>::iterator > &best_itrs)
{
  // clean those reads that have no good mapping location
  size_t j = 0;
  for (size_t i = 0; i < read_index.size(); ++i)
    if (best_maps[i].mr.size() > 0) {
      std::swap(read_index[i], read_index[j]);
      std::swap(best_maps[i], best_maps[j]);
      ++j;
    }
  read_index.erase(read_index.begin() + j, read_index.end());
  best_maps.erase(best_maps.begin() + j, best_maps.end());
  
  vector<pair<unsigned int, vector<MultiMapResult>::iterator > > sorter;
  for (size_t i = 0; i < best_maps.size(); ++i) 
    sorter.push_back(make_pair(read_index[i], best_maps.begin() + i));

  sort(sorter.begin(), sorter.end(),
       indexed_best_less<vector<MultiMapResult>::iterator>());

  best_itrs.resize(best_maps.size());
  for (size_t i = 0; i < sorter.size(); ++i) {
    read_index[i] = sorter[i].first;
    best_itrs[i] = sorter[i].second;
  }
}

static void
extract_adaptors(const string &adaptor, string & T_adaptor, string & A_adaptor)
{
  const vector<string> adaptors(smithlab::split(adaptor, ":"));
  if (adaptors.size() == 0)
    T_adaptor = A_adaptor = "";
  else if (adaptors.size() == 1)
    T_adaptor = A_adaptor = adaptors.front();
  else if (adaptors.size() == 2)
    {
      T_adaptor = adaptors[0];
      A_adaptor = adaptors[1];
    }
  else
    throw SMITHLABException("Invalid adaptors: " + adaptor + "\n"
			    + "Format: T_adaptor:A_adaptor");
}

static void
map_reads(const string &chrom_file, const string &filenames_file,
	  const string &fasta_suffix, const string &reads_file,
	  const string &adaptor_sequence, 
	  const size_t read_start_index, const size_t n_reads_to_process,
	  const size_t n_seeds, const size_t seed_weight, size_t read_width,
	  size_t max_mismatches, const double wildcard_cutoff,
	  const bool FASTER_MODE, const bool QUALITY, const bool ALLOW_METH_BIAS,
	  const bool WILDCARD, const bool WILD_N_MODE, const bool ORIGINAL_OUTPUT,
	  const bool VERBOSE,
	  const bool AG_WILDCARD, 
	  vector<unsigned int> &read_index,
	  vector<MultiMapResult> &best_maps)
{
  FastReadWC::set_cutoff(wildcard_cutoff);
  if (WILDCARD && wildcard_cutoff != numeric_limits<double>::max() &&
      (wildcard_cutoff > 1.0 || wildcard_cutoff < 0)) 
    throw SMITHLABException("wildcard cutoff must be in [0, 1]");
    
  //////////////////////////////////////////////////////////////
  //  CHECK HOW QUALITY SCORES ARE USED
  //
  const size_t INPUT_MODE = get_input_mode(VERBOSE, reads_file);
  if (INPUT_MODE == FASTA_FILE && !ORIGINAL_OUTPUT)
    throw SMITHLABException("when reads are given as FASTA, "
			    "original output format is required");
    
  const size_t RUN_MODE = get_run_mode(VERBOSE, INPUT_MODE, WILDCARD, QUALITY);
    
  vector<string> chrom_files;
  identify_chromosomes(VERBOSE, filenames_file, fasta_suffix,
		       chrom_file, chrom_files);

  // extract adatpers
  string T_adaptor, A_adaptor;
  extract_adaptors(adaptor_sequence, T_adaptor, A_adaptor);
  const string adaptor = AG_WILDCARD ? A_adaptor : T_adaptor;
    
  //////////////////////////////////////////////////////////////
  // OBTAIN THE READS
  // 
  vector<FastRead> fast_reads;
  vector<FastReadWC> fast_reads_wc;
  vector<FastReadQuality> fast_reads_q;
  vector<size_t> read_words;
  load_reads(VERBOSE, INPUT_MODE, RUN_MODE, AG_WILDCARD,
	     max_mismatches, adaptor, reads_file, 
	     read_start_index, n_reads_to_process,
	     fast_reads, fast_reads_wc, fast_reads_q, 
	     read_index, read_words, read_width);
  
  
  if (max_mismatches == numeric_limits<size_t>::max())
    max_mismatches = static_cast<size_t>(0.07*read_width);
      
  if (VERBOSE)
    cerr << "MAX MISMATCHES: " << max_mismatches << endl;
	
  double max_match_score = max_mismatches*FastReadQuality::get_scaler();

  /////////////////////////////////////////////////////////////
  // INITIALIZE THE SEED STRUCTURES
  //
  vector<size_t> the_seeds;
  load_seeds(VERBOSE, FASTER_MODE,
	     read_width, n_seeds, seed_weight, the_seeds);

  //////////////////////////////////////////////////////////////
  // THIS IS WHERE THE ACTUAL MAPPING HAPPENS
  //
  vector<size_t> chrom_sizes;
  vector<string> chrom_names;
  vector<pair<unsigned int, unsigned int> > ambigs;
    
  const wildcard_score specialized_score;

  best_maps = vector<MultiMapResult>(
				     read_words.size(), 
				     MultiMapResult((RUN_MODE == RUN_MODE_WEIGHT_MATRIX) ?
						    static_cast<size_t>(max_match_score) : 
						    max_mismatches));

  iterate_over_seeds(VERBOSE, AG_WILDCARD, ALLOW_METH_BIAS,
		     specialized_score,
		     the_seeds, chrom_files, ambigs, 
		     chrom_names, chrom_sizes,
		     fast_reads, // USE REGULAR FAST READS
		     read_words, read_index,
		     best_maps, max_mismatches, read_width);

  // First make sure the chrom names don't have spaces (cause
  // problems for later processing)
  for (size_t i = 0; i < chrom_names.size(); ++i) 
    {
      const size_t chr_name_end = chrom_names[i].find_first_of(" \t");
      if (chr_name_end != string::npos)
	chrom_names[i].erase(chr_name_end);
    }
  assert(chrom_names.size() == chrom_sizes.size());
  assert(MapResult::chrom_names.empty()
	 || std::equal(MapResult::chrom_names.begin(),
		       MapResult::chrom_names.end(), chrom_names.begin()));
  
  MapResult::chrom_names = chrom_names;
  MapResult::chrom_sizes = chrom_sizes;
}

static bool
fast_forward_fastq(std::istream &in, const size_t n_entries)
{
  const size_t INPUT_BUFFER_SIZE = 10000;
  char buffer[INPUT_BUFFER_SIZE + 1];

  // fast forward 
  for (size_t i = 0; i < n_entries && in.good(); ++i)
    for (size_t line_count = 0; line_count < 4 && in.good(); ++line_count)
      in.getline(buffer, INPUT_BUFFER_SIZE);
  return in;
}

static bool
get_entry_fastq(std::istream &in,
                const size_t curr_idx,  // the first unread entry index
                const size_t target_idx, // the entry index to be read
                string &name, string &seq, string &scr)
{
  static const size_t INPUT_BUFFER_SIZE = 10000;
  static char buffer[INPUT_BUFFER_SIZE + 1];

  // fast forward 
  for (size_t i = curr_idx; i < target_idx && in.good(); ++i)
    for (size_t line_count = 0; line_count < 4 && in.good(); ++line_count)
      in.getline(buffer, INPUT_BUFFER_SIZE);

  // read and parse the target read
  for (size_t line_count = 0; line_count < 4 && in.good(); ++line_count)
  {
    in.getline(buffer, INPUT_BUFFER_SIZE);
    if (in.gcount() == static_cast<int>(INPUT_BUFFER_SIZE))
      throw SMITHLABException("Line exceeds max length: " +
                              toa(INPUT_BUFFER_SIZE));
    
    if (line_count % 4 == 0) 
    {
      name = string(buffer + 1);
      const size_t name_end = name.find_first_of(" \t");
      if (name_end != string::npos)
        name.erase(name.begin() + name_end, name.end());
    }
    else if (line_count % 4 == 1) 
    {
      seq = string(buffer);
    }
    else if (line_count % 4 == 3) 
    {
      scr = string(buffer);
    }
  }
  
  return in;
}

inline static bool
same_read(const size_t suffix_len, 
	  const string &sa, const string &sb) {
  return std::equal(sa.begin(), sa.end() - suffix_len, sb.begin());
}

inline static bool
name_smaller(const size_t suffix_len, 
	     const string &sa, const string &sb) {
  return std::lexicographical_compare(sa.begin(), sa.end() - suffix_len, 
				      sb.begin(), sb.end() - suffix_len);
}

static MappedRead
MapResult_to_MappedRead(const MapResult &r, const string &name,
                        const string &seq, const string &scr,
                        const size_t mismatch)
{
  const size_t chrom_id = r.chrom;
  assert(chrom_id < MapResult::chrom_sizes.size());
  const size_t read_len = seq.length();
  size_t start = r.strand ? r.site : 
    MapResult::chrom_sizes[chrom_id] - r.site - read_len;
  size_t end = start + read_len;
  if (start >= end)
    throw SMITHLABException("invalid mapped reads\n" + 
			    r.tostring() + "\t" + seq + "\t" + scr);
  const char strand = ((r.strand) ? '+' : '-');
  GenomicRegion region(MapResult::chrom_names[chrom_id], start, end,
		       name, mismatch, strand);
  MappedRead mr;
  mr.r = region;
  mr.seq = seq;
  mr.scr = scr;
  return mr;
}

static bool
merge_mates(const size_t suffix_len,
	    const size_t MAX_SEGMENT_LENGTH,
	    const MappedRead &one, const MappedRead &two,
            MappedRead &merged) {

  if (!(one.r.same_chrom(two.r) && one.r.get_strand() == two.r.get_strand()))
    return false;
  
  merged = one;
  size_t start_one = 0, end_one = 0, start_two = 0, 
    end_two = 0, start_overlap = 0, end_overlap = 0;
  int len = 0;
  if (merged.r.pos_strand()) {
    start_overlap = std::max(one.r.get_start(), two.r.get_start());
    end_overlap = std::min(one.r.get_end(), two.r.get_end());
    start_one = one.r.get_start();
    end_one = std::min(start_overlap, one.r.get_end());
    start_two = std::max(end_overlap, two.r.get_start());
    end_two = two.r.get_end();
    len = end_two - start_one;
    merged.r.set_start(start_one);
    merged.r.set_end(end_two);
  }
  else { // if (merged.r.neg_strand())
    start_overlap = std::max(one.r.get_start(), two.r.get_start());
    end_overlap = std::min(one.r.get_end(), two.r.get_end());
    start_one = std::max(end_overlap, one.r.get_start());
    end_one = one.r.get_end();
    start_two = two.r.get_start();
    end_two = std::min(start_overlap, two.r.get_end());
    len = end_one - start_two;
    merged.r.set_start(start_two);
    merged.r.set_end(end_one);
  }
  
  assert(end_one >= start_one && end_two >= start_two);
  assert(start_overlap >= end_overlap || 
	 static_cast<size_t>(len) == end_one - start_one
	 + end_two - start_two + end_overlap - start_overlap);
  
  merged.r.set_score(one.r.get_score() + two.r.get_score());
  
  if (len > 0 && len <= static_cast<int>(MAX_SEGMENT_LENGTH)) {
    merged.seq = string(len, 'N');
    merged.scr = string(len, 'B');
    const string name(one.r.get_name());
    merged.r.set_name("FRAG:" + name.substr(0, name.size() - suffix_len));
    
    // lim_one ios the offset within the merged sequence where the
    // overlapping portion begins
    const size_t lim_one = end_one - start_one;
    copy(one.seq.begin(), one.seq.begin() + lim_one, merged.seq.begin());
    copy(one.scr.begin(), one.scr.begin() + lim_one, merged.scr.begin());
    const size_t lim_two = end_two - start_two;
    copy(two.seq.end() - lim_two, two.seq.end(), merged.seq.end() - lim_two);
    copy(two.scr.end() - lim_two, two.scr.end(), merged.scr.end() - lim_two);
    
    // deal with overlapping part
    if (start_overlap < end_overlap) {

      const int info_one = one.seq.length() - 
	count(one.seq.begin(), one.seq.end(), 'N') - 
	static_cast<int>(one.r.get_score());
      const int info_two = two.seq.length() - 
	count(two.seq.begin(), two.seq.end(), 'N') - 
	static_cast<int>(two.r.get_score());

      // use the mate with the most information ("info") to fill in
      // the overlapping portion
      if (info_one >= info_two) {
	const size_t source_start = merged.r.pos_strand() ?
	  (start_overlap - one.r.get_start()) : (one.r.get_end() - end_overlap);
	const size_t source_end = merged.r.pos_strand() ?
	  (end_overlap -  one.r.get_start()) : (one.r.get_end() - start_overlap);
	copy(one.seq.begin() + source_start, one.seq.begin() + source_end,
	     merged.seq.begin() + lim_one);
	copy(one.scr.begin() + source_start, one.scr.begin() + source_end,
	     merged.scr.begin() + lim_one);
      }
      else { // if (info_two > info_one)
	const size_t source_start = merged.r.pos_strand() ?
	  (start_overlap - two.r.get_start()) : (two.r.get_end() - end_overlap);
	const size_t source_end = merged.r.pos_strand() ?
	  (end_overlap -  two.r.get_start()) : (two.r.get_end() - start_overlap);
	copy(two.seq.begin() + source_start, two.seq.begin() + source_end,
	     merged.seq.begin() + lim_one);
	copy(two.scr.begin() + source_start, two.scr.begin() + source_end,
	     merged.scr.begin() + lim_one);
      }
    }
    return true;
  }
  return false;
}


class SwitchStrand : public std::unary_function<MapResult, void>
{
public:
  SwitchStrand(const size_t len) : read_len(len) {}
  void operator()(MapResult &r) const
  {
    r.strand = !r.strand;
    r.site = MapResult::chrom_sizes[r.chrom] - r.site - read_len;
  }
private:
  size_t read_len;
};


static void
purge_remaining(const bool A_RICH,
		const string &reads_file, 
		const vector<unsigned int> &read_index,
		vector<MultiMapResult> &best_maps,
        vector<vector<MultiMapResult>::iterator > &best_itrs,
		string &adaptor, size_t read_idx, size_t curr_idx,
        string &name, string &seq, string &scr, string &prev_name,
		std::ifstream &in, std::ofstream &out) {
  
  while (curr_idx < read_index.size()) {
    const vector<MultiMapResult>::const_iterator bm(best_itrs[curr_idx]);
    if (bm->mr.size() == 1) {
      const size_t curr_idx_val = read_index[curr_idx];
      
      if (read_idx != curr_idx_val) {
        if (get_entry_fastq(in, read_idx + 1, curr_idx_val, name, seq, scr)) {
          if (!adaptor.empty())
            clip_adaptor_from_read(adaptor, MIN_ADAPTOR_MATCH_SCORE, seq);
          read_idx = curr_idx_val;
          if (A_RICH) {
            revcomp_inplace(seq);
            std::reverse(scr.begin(), scr.end());
          }
        }
        else throw SMITHLABException("Error reading " + reads_file);
      }
      
      out << MapResult_to_MappedRead(bm->mr.front(), name, 
				     seq, scr, bm->score) << endl;
    }
    ++curr_idx;
  }
}

static void
clip_mates(const string &T_reads_file, vector<unsigned int> &t_read_index,
           vector<MultiMapResult> &t_best_maps,
           vector<vector<MultiMapResult>::iterator > &t_best_itrs,
           const string &A_reads_file, vector<unsigned int> &a_read_index,
           vector<MultiMapResult> &a_best_maps, 
           vector<vector<MultiMapResult>::iterator > &a_best_itrs,
           const size_t read_start_index,
           const size_t suffix_len, const size_t MAX_SEGMENT_LENGTH,
           const string &adaptor_sequence, const string &outfile, 
           const bool VERBOSE) {

// extract adatpers
  string T_adaptor, A_adaptor;
  extract_adaptors(adaptor_sequence, T_adaptor, A_adaptor);
  
  std::ifstream t_in(T_reads_file.c_str());
  if (!(t_in.good() && fast_forward_fastq(t_in, read_start_index)))
    throw SMITHLABException("cannot open input file " + T_reads_file);
  
  std::ifstream a_in(A_reads_file.c_str());
  if (!(a_in.good() && fast_forward_fastq(a_in, read_start_index)))
    throw SMITHLABException("cannot open input file " + A_reads_file);

  std::ofstream out(outfile.c_str());
  if (!out)
    throw SMITHLABException("cannot open output file " + outfile);
  
  size_t t_read_idx = 0, t_curr_idx = 0;
  string t_name, t_seq, t_scr, prev_t_name;
  if (get_entry_fastq(t_in, 0, 0, t_name, t_seq, t_scr)) {
    if (!T_adaptor.empty())
      clip_adaptor_from_read(T_adaptor, MIN_ADAPTOR_MATCH_SCORE, t_seq);
  }
  else throw SMITHLABException("Error reading " + T_reads_file);

  size_t a_read_idx = 0, a_curr_idx = 0;
  string a_name, a_seq, a_scr, prev_a_name;
  if (get_entry_fastq(a_in, 0, 0, a_name, a_seq, a_scr)) {
    if (!A_adaptor.empty())
      clip_adaptor_from_read(A_adaptor, MIN_ADAPTOR_MATCH_SCORE, a_seq);
    revcomp_inplace(a_seq);
    std::reverse(a_scr.begin(), a_scr.end());
  }
  else throw SMITHLABException("Error reading " + A_reads_file);
  
  /* Looping over both the T-rich and A-rich reads in parallel
   */
  while (t_curr_idx < t_read_index.size() && a_curr_idx < a_read_index.size()) {
    
    // Move to the appropriate T-rich read
    const size_t t_curr_idx_val = t_read_index[t_curr_idx];
    if (t_read_idx != t_curr_idx_val) {
      if (get_entry_fastq(t_in, t_read_idx + 1, t_curr_idx_val,
                          t_name, t_seq, t_scr)) {
        if (!T_adaptor.empty())
          clip_adaptor_from_read(T_adaptor, MIN_ADAPTOR_MATCH_SCORE, t_seq);
        t_read_idx = t_curr_idx_val;
      }
      else throw SMITHLABException("Error reading " + T_reads_file);
    }
    
    // Move to the appropriate A-rich read
    const size_t a_curr_idx_val = a_read_index[a_curr_idx];
    if (a_read_idx != a_curr_idx_val) {
      if (get_entry_fastq(a_in,  a_read_idx + 1, a_curr_idx_val,
                          a_name, a_seq, a_scr)) {
        if (!A_adaptor.empty())
          clip_adaptor_from_read(A_adaptor, MIN_ADAPTOR_MATCH_SCORE, a_seq);
        revcomp_inplace(a_seq);
        std::reverse(a_scr.begin(), a_scr.end());
        a_read_idx = a_curr_idx_val;
      }
      else throw SMITHLABException("Error reading " + T_reads_file);
    }
    
    const vector<MultiMapResult>::const_iterator tbm(t_best_itrs[t_curr_idx]);
    const vector<MultiMapResult>::iterator abm(a_best_itrs[a_curr_idx]);
    
    /*
      (1) THE A-RICH AND T-RICH READS CORRESPOND
     */
    if (t_curr_idx_val == a_curr_idx_val) {

      assert(same_read(suffix_len, t_name, a_name));
      std::for_each(abm->mr.begin(), abm->mr.end(), SwitchStrand(a_seq.size()));
      
      const size_t t_map_count = tbm->mr.size();
      const size_t a_map_count = abm->mr.size();

      /// TODO: in this loop we should only be CHECKING if the merging
      /// can happen, not actually doing the merging
      vector<MappedRead> valid_frags;
      for (size_t i = 0; i < t_map_count && valid_frags.size() <= 1; ++i) {
	const MappedRead 
	  one(MapResult_to_MappedRead(tbm->mr[i], t_name, t_seq, t_scr, tbm->score));
	
	for (size_t j = 0; j < a_map_count && valid_frags.size() <= 1; ++j) {
	  const MappedRead 
	    two(MapResult_to_MappedRead(abm->mr[j], a_name, a_seq, a_scr, abm->score));
	  
	  MappedRead merged;
	  if (merge_mates(suffix_len, MAX_SEGMENT_LENGTH, one, two, merged))
	    valid_frags.push_back(merged);
	}
      }
      
      if (valid_frags.size() == 1)
          out << valid_frags.front() << endl;
      else {
	if (t_map_count == 1)
	  out << MapResult_to_MappedRead(tbm->mr.front(), t_name, t_seq, t_scr, 
					 tbm->score) << endl;
	if (a_map_count == 1)
	  out << MapResult_to_MappedRead(abm->mr.front(), a_name, a_seq, a_scr, 
					 abm->score) << endl;
      }
      ++t_curr_idx;
      ++a_curr_idx;
    }
    
    /*
      (2) THE T-RICH READ COMES FIRST
     */
    else if (t_curr_idx_val < a_curr_idx_val) {
      if (tbm->mr.size() == 1)
	out << MapResult_to_MappedRead(tbm->mr.front(), t_name, t_seq, t_scr, 
				       tbm->score) << endl;
      ++t_curr_idx;
    }
    /*
      (3) THE A-RICH READ COMES FIRST
     */
    else { // if (a_curr_idx_val < t_curr_idx_val)
      if (abm->mr.size() == 1) {
	const SwitchStrand switcher(a_seq.size());
	std::for_each(abm->mr.begin(), abm->mr.end(), switcher);
	out << MapResult_to_MappedRead(abm->mr.front(), a_name, a_seq, a_scr, 
				       abm->score) << endl;
      }
      ++a_curr_idx;
    }
  }
  
  purge_remaining(false, T_reads_file, t_read_index, t_best_maps, t_best_itrs,
    	  T_adaptor, t_read_idx, t_curr_idx,
          t_name, t_seq, t_scr, prev_t_name, t_in, out);
  purge_remaining(true, A_reads_file, a_read_index, a_best_maps, a_best_itrs,
    	  A_adaptor, a_read_idx, a_curr_idx,
    	  a_name, a_seq, a_scr, prev_a_name, a_in, out);
}


int 
main(int argc, const char **argv) 
{
  try 
    {
      string chrom_file;
      string filenames_file;
      string outfile("/dev/stdout");
      string ambiguous_file;
      string fasta_suffix = "fa";
      string adaptor_sequence;
    
      size_t n_seeds = 3;
      size_t seed_weight = 11;
      size_t read_width = 0;
      size_t max_mismatches = numeric_limits<size_t>::max();
      size_t max_mappings = 500;
      double wildcard_cutoff = numeric_limits<double>::max();

      size_t MAX_SEGMENT_LENGTH = 500;
      size_t suffix_len = 1;
        
      bool VERBOSE = false;
      bool FASTER_MODE = true;
      bool QUALITY = false;
      bool ALLOW_METH_BIAS = false;
      bool WILDCARD = false;

      bool WILD_N_MODE = true;
      bool ORIGINAL_OUTPUT = false;

      size_t read_start_index = 0;
      size_t n_reads_to_process = std::numeric_limits<size_t>::max();
      
      /****************** COMMAND LINE OPTIONS ********************/
      OptionParser opt_parse(strip_path(argv[0]),
			     "The rmapbs mapping tool for Solexa reads"
			     " following bisulfite treatment",
			     "<T-rich-reads> <A-rich-reads>" );
      opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
			false , outfile);
      opt_parse.add_opt("chrom", 'c', "FASTA file or dir containing chromosome(s)", 
			false , chrom_file);
      opt_parse.add_opt("start", 'T', "index of first read to map", 
			false , read_start_index);
      opt_parse.add_opt("number", 'N', "number of reads to map", 
			false , n_reads_to_process);
      opt_parse.add_opt("suffix", 's', "suffix of FASTA files "
			"(assumes -c indicates dir)", false , fasta_suffix);
      opt_parse.add_opt("filenames", 'F', "file listing names of "
			"chromosome files", false , filenames_file);
      opt_parse.add_opt("seeds", 'S', "number of seeds", false , n_seeds);
      opt_parse.add_opt("hit", 'h', "weight of hit", false , seed_weight);
      opt_parse.add_opt("width", 'w', "width of reads", false, read_width);
      opt_parse.add_opt("mismatch", 'm', "maximum allowed mismatches", 
			false , max_mismatches);
      opt_parse.add_opt("max-map", 'M',
			"maximum allowed mappings for SE read (default "
			+ smithlab::toa(max_mappings) + ")", 
			false, max_mappings);
      opt_parse.add_opt("qual", 'Q', "use quality scores (input must be FASTQ)", 
			false, QUALITY);
      opt_parse.add_opt("bias", 'B', "allow CpG non-conversion to assist", 
			false, ALLOW_METH_BIAS);
      opt_parse.add_opt("faster-off", 'f', "turn off faster "
			"seeds (more sensitive)", false, FASTER_MODE);
      opt_parse.add_opt("clip", 'C', "clip the specified adaptor", 
			false, adaptor_sequence);
      opt_parse.add_opt("max-seg-len", 'L', "MAX DNA segment length", 
			false, MAX_SEGMENT_LENGTH);
      opt_parse.add_opt("suffix-len", '\0', "Suffix length of reads name", 
			false, suffix_len);
      opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
      opt_parse.add_opt("n-mismatch", 'y', "make N in read "
			"cause mismatch",  false, WILD_N_MODE);
      vector<string> leftover_args;
      opt_parse.parse(argc, argv, leftover_args);
      if (argc == 1 || opt_parse.help_requested()) 
        {
	  cerr << opt_parse.help_message() << endl;
	  return EXIT_SUCCESS;
        }
      if (opt_parse.about_requested()) 
        {
	  cerr << opt_parse.about_message() << endl;
	  return EXIT_SUCCESS;
        }
      if (opt_parse.option_missing()) 
        {
	  cerr << opt_parse.option_missing_message() << endl;
	  return EXIT_SUCCESS;
        }
      if (leftover_args.empty()) 
        {
	  cerr << opt_parse.help_message() << endl;
	  return EXIT_SUCCESS;
        }
      if (chrom_file.empty() && filenames_file.empty()) 
        {
	  cerr << "must specify chroms file/dir or filenames file" << endl;
	  return EXIT_FAILURE;
        }
      if (leftover_args.size() != 2)
        {
	  cerr << "must gives T-rich read file and A-rich read file" << endl;
	  cerr << opt_parse.help_message() << endl;
	  return EXIT_FAILURE;
        }
        
      const string T_reads_file = leftover_args.front();
      const string A_reads_file = leftover_args.back();
      
      // the index user gives is 1-based, rmap internal representation
      // is 0-based
      read_start_index = read_start_index == 0 ? 0 : read_start_index - 1;
      
      /****************** END COMMAND LINE OPTIONS *****************/
      //////////////////////////////////////////////////////////////


      if (VERBOSE)
        cerr << "[MAPPING T-RICH MATES] " << T_reads_file << endl;
      vector<unsigned int> t_read_index;
      MultiMapResult::init(max_mappings); //  max_mappings >= 500
      vector<MultiMapResult> t_best_maps;
      const static bool TC_WILDCARD = false;
      map_reads(chrom_file, filenames_file, fasta_suffix,
		T_reads_file, adaptor_sequence, 
		read_start_index, n_reads_to_process,
		n_seeds, seed_weight,
		read_width, max_mismatches, wildcard_cutoff,
		FASTER_MODE, QUALITY, ALLOW_METH_BIAS,
		WILDCARD, WILD_N_MODE, ORIGINAL_OUTPUT, VERBOSE,
		TC_WILDCARD, t_read_index, t_best_maps);

      vector<vector<MultiMapResult>::iterator > t_best_itrs;
      clean_and_invert_bests_list(t_read_index, t_best_maps, t_best_itrs);
      if (VERBOSE)
        cerr << "[T-RICH CANDIDATES: ] " << t_read_index.size() << endl;
      
      if (VERBOSE)
	cerr << endl << "[MAPPING A-RICH MATES] " << A_reads_file << endl;
      vector<unsigned int> a_read_index;
      MultiMapResult::init(max_mappings); //  max_mappings >= 500
      vector<MultiMapResult> a_best_maps;
      const static bool AG_WILDCARD = true;
      map_reads(chrom_file, filenames_file, fasta_suffix,
		A_reads_file, adaptor_sequence, 
		read_start_index, n_reads_to_process,
		n_seeds, seed_weight,
		read_width, max_mismatches, wildcard_cutoff,
		FASTER_MODE, QUALITY, ALLOW_METH_BIAS,
		WILDCARD, WILD_N_MODE, ORIGINAL_OUTPUT, VERBOSE,
		AG_WILDCARD, a_read_index, a_best_maps);

      vector<vector<MultiMapResult>::iterator > a_best_itrs;
      clean_and_invert_bests_list(a_read_index, a_best_maps, a_best_itrs);

      if (VERBOSE)
        cerr << "[A-RICH CANDIDATES: ] " << a_read_index.size() << endl;

      if (VERBOSE)
	cerr << "[JOINING MATES]" << endl;
      clip_mates(
        T_reads_file, t_read_index, t_best_maps, t_best_itrs,
        A_reads_file, a_read_index, a_best_maps, a_best_itrs, 
        read_start_index,        
        suffix_len, MAX_SEGMENT_LENGTH, adaptor_sequence, 
        outfile, VERBOSE);
      if (VERBOSE)
	cerr << "[MAPPING DONE]" << endl;
    }
  catch (const SMITHLABException &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
