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
get_read_words(const vector<string> &reads, vector<size_t> &read_words) {
  read_words.resize(reads.size());
  for (size_t i = 0; i < reads.size(); ++i) {
    const size_t truncate_to = min(reads[i].length(), SeedMaker::max_seed_part);
    // Need to replace the Ns because otherwise they will destroy the
    // conversion from DNA to integers. Could replace with random
    // bases, but everyone hates non-deterministic programs.
    string s;
    replace_copy(reads[i].end() - truncate_to, reads[i].end(),
		 std::back_inserter(s), 'N', 'A');
    read_words[i] = SeedMaker::make_read_word(s);
  }
}


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
sort_by_key(SeedHashSorter &sh_sorter, vector<MultiMapResultPE> &best_maps,
	    vector<size_t> &reads, vector<vector<size_t> > &read_index, 
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
  seed_hash[prev_key] = make_pair(frb + prev_idx, frb + curr_idx);
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
map_reads(const string &chrom, const size_t chrom_id, const size_t profile, 
	  const size_t read_width, const size_t max_diffs, 
	  const size_t min_sep, const size_t max_sep, 
	  const vector<T> &reads_left, const vector<T> &reads_right,
	  const typename SeedHash<T>::type &seed_hash, const bool strand, 
	  vector<MultiMapResultPE> &best_maps) {
  
  MASK_t bad_bases = rmap_bits::all_ones;
  MASK_t read_word = rmap_bits::all_zeros;
  T fast_read;
  
  const size_t chrom_size = chrom.size();
  
  vector<T> possible_lefts(max_sep);
  
  size_t chrom_offset = 0;
  for (; chrom_offset < read_width - 1; ++chrom_offset) {
    const size_t key_base = base2int(chrom[chrom_offset]);
    fast_read.shift(key_base);
    SeedMaker::update_bad_bases(key_base, bad_bases);
    SeedMaker::update_read_word(key_base, read_word);
  }
  
  for (; chrom_offset < chrom_size; ++chrom_offset) {
    const size_t key_base = base2int(chrom[chrom_offset]);
    fast_read.shift(key_base);
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
	    const size_t idx = (to_test - reads_right.begin());
	    const size_t lookback_limit = chrom_offset - min_sep;
	    for (size_t i = (chrom_offset > max_sep) ? 
		   chrom_offset - max_sep : 0; i < lookback_limit; ++i) {
	      const size_t pair_score = score + 
		reads_left[idx].score(possible_lefts[i % max_sep]);
	      const vector<MultiMapResultPE>::iterator current(best_maps.begin() + idx);
	      if (pair_score <= current->score) {
		const size_t left_start = i - read_width + 1;
		const size_t right_start = chrom_offset - read_width + 1;
		current->add(pair_score, chrom_id, left_start, right_start, strand);
	      }
	    }
	  }
	}
      }
      possible_lefts[chrom_offset % max_sep] = fast_read;
    }
  }
}


template <class T> void
iterate_over_seeds(const bool VERBOSE, 
		   const vector<size_t> &the_seeds, const vector<string> &chrom_files,
		   vector<size_t> &ambigs, 
		   vector<string> &chrom_names, vector<size_t> &chrom_sizes,
		   vector<T> &fast_reads_left, vector<T> &fast_reads_right,
		   vector<size_t> &read_words, 
		   vector<vector<size_t> > &read_index,
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
    
    for (size_t i = 0; i < chrom_files.size() && !fast_reads_left.empty(); ++i) {
      
      vector<string> tmp_chrom_names, chroms;
      if (VERBOSE)
	cerr << "[SEED:" << j + 1 << "/" << the_seeds.size() << "] "
	     << "[LOADING CHROM] ";
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

	map_reads(chroms[k], i, the_seeds[j], read_width, max_mismatches,
		  min_sep, max_sep, fast_reads_left, fast_reads_right,
		  seed_hash, true, best_maps);
	revcomp_inplace(chroms[k]);
	map_reads(chroms[k], i, the_seeds[j], read_width, max_mismatches,
		  min_sep, max_sep, fast_reads_left, fast_reads_right,
		  seed_hash, false, best_maps);
	const clock_t end(clock());
	
	if (VERBOSE)
	  cerr << "[" << static_cast<float>(end - start)/CLOCKS_PER_SEC << " SEC] ";
      }
      if (VERBOSE) cerr << endl;
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


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////
////  CODE FOR VALIDATING THE READS TO ENSURE QUALITY OF MAPPED READS
////
////
// NOTE: BAD POSITIONS ARE THOSE FOR WHICH THE HIGHEST QUALITY SCORE
// IS LESS THAN OR EQUAL TO -5 AND AT MOST 20 SUCH POSITIONS ARE
// ALLOWED IN A GOOD READ.
static bool
good_read(const vector<vector<double> > &read) {
  size_t bad_count = 0;
  for (size_t i = 0; i < read.size(); ++i)
    bad_count += (*std::max_element(read[i].begin(), read[i].end()) <= -5.0);
  return (bad_count <= 20);
}


static void
clean_reads(const size_t min_match_score, const size_t read_width,
	    const vector<size_t> &input_read_index,
	    vector<string> &reads_left, 
	    vector<string> &reads_right, 
	    vector<vector<vector<double> > > &scores_left,
	    vector<vector<vector<double> > > &scores_right,
	    vector<vector<size_t> > &read_index) {
  vector<vector<vector<double> > > good_scores_left, good_scores_right;
  vector<string> good_reads_left, good_reads_right;
  for (size_t i = 0; i < reads_left.size(); ++i) {
    if (good_read(scores_left[i]) && good_read(scores_right[i])) {
      if (reads_left[i].length() > read_width) {
	reads_left[i] = reads_left[i].substr(0, read_width);
	scores_left[i].erase(scores_left[i].begin() + read_width, scores_left[i].end());
	reads_right[i] = reads_right[i].substr(0, read_width);
	scores_right[i].erase(scores_right[i].begin() + read_width, scores_right[i].end());
      }
      read_index.push_back(vector<size_t>(1, input_read_index[i]));
      good_reads_left.push_back(reads_left[i]);
      good_scores_left.push_back(scores_left[i]);
      good_reads_right.push_back(reads_right[i]);
      good_scores_right.push_back(scores_right[i]);
    }
  }
  reads_left.swap(good_reads_left);
  scores_left.swap(good_scores_left);
  
  reads_right.swap(good_reads_right);
  scores_right.swap(good_scores_right);
}


static char
to_base_symbol(char c) {return (isvalid(c)) ? toupper(c) : 'N';}


static void
clean_reads(const size_t max_diffs, const size_t read_width,
	    const vector<size_t> &input_read_index,
	    vector<string> &reads_left, 
	    vector<string> &reads_right, 
	    vector<vector<size_t> > &read_index) {
  vector<pair<pair<string, string>, size_t> > sorter;
  for (size_t i = 0; i < reads_left.size(); ++i) {
    transform(reads_left[i].begin(), reads_left[i].end(), reads_left[i].begin(),
	      std::ptr_fun(&to_base_symbol));
    transform(reads_right[i].begin(), reads_right[i].end(), reads_right[i].begin(),
	      std::ptr_fun(&to_base_symbol));
    if (count(reads_left[i].begin(), reads_left[i].end(), 'N') <=
	static_cast<int>(max_diffs) &&
	count(reads_right[i].begin(), reads_right[i].end(), 'N') <=
	static_cast<int>(max_diffs))
      sorter.push_back(make_pair(make_pair(reads_left[i], reads_right[i]),
				 input_read_index[i]));
  }
  sort(sorter.begin(), sorter.end());
  reads_left.clear();
  reads_right.clear();
  for (size_t i = 0; i < sorter.size(); ++i) {
    if (i == 0 || sorter[i - 1].first != sorter[i].first) {
      reads_left.push_back(sorter[i].first.first);
      reads_right.push_back(sorter[i].first.second);
      read_index.push_back(vector<size_t>(1, sorter[i].second));
    }
    else read_index.back().push_back(sorter[i].second);
  }
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
		 const vector<string> &read_names, vector<MultiMapResultPE> &bests, 
		 vector<GenomicRegion> &hits) {
  
  static const string LEFT_TAG("_L");
  static const string RIGHT_TAG("_R");

  for (size_t i = 0; i < bests.size(); ++i)
    if (!bests[i].mr.empty()) {
      bests[i].sort();
      for (size_t j = 0; j < bests[i].mr.size(); ++j)
	if (j == 0 || bests[i].mr[j - 1] < bests[i].mr[j]) {
	  const size_t chrom_id = bests[i].mr[j].chrom;
	  const size_t left_start = 
	    bests[i].mr[j].strand ? bests[i].mr[j].site : 
	    chrom_sizes[chrom_id] - bests[i].mr[j].site2 - read_len;
	  const size_t right_start = 
	    bests[i].mr[j].strand ? bests[i].mr[j].site2 :
	    chrom_sizes[chrom_id] - bests[i].mr[j].site - read_len;
	  const size_t score = bests[i].score;
	  const char strand = ((bests[i].mr[j].strand) ? '+' : '-');
	  for (size_t k = 0; k < read_index[i].size(); ++k) {
	    const string read_name(read_names[read_index[i][k]]);
	    hits.push_back(GenomicRegion(chrom[chrom_id], left_start, 
					 left_start + read_len,
					 read_name + LEFT_TAG, 
					 score, strand));
	    hits.push_back(GenomicRegion(chrom[chrom_id], right_start,
					 right_start + read_len,
					 read_name + RIGHT_TAG, 
					 score, strand));
	  }
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
eliminate_ambigs(const size_t max_mismatches, vector<MultiMapResultPE> &best_maps, 
		 vector<vector<size_t> > &read_index, vector<size_t> &reads, 
		 vector<size_t> &ambigs, vector<T> &fast_reads_left, 
		 vector<T> &fast_reads_right) {
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
process_score_data(const size_t read_width, const size_t max_mismatches,
		   double &max_quality_score, double &max_match_score,
		   vector<vector<vector<double> > > &left,
		   vector<vector<vector<double> > > &right) {
  
  double min_quality_score = numeric_limits<double>::max();
  for (size_t i = 0; i < left.size(); ++i)
    for (size_t j = 0; j < left[i].size(); ++j) {
      max_quality_score =
	max(max_quality_score, 
	    max(*max_element(left[i][j].begin(), left[i][j].end()),
		*max_element(right[i][j].begin(), right[i][j].end())));
      min_quality_score =
	min(min_quality_score, 
	    min(*min_element(left[i][j].begin(), left[i][j].end()),
		*min_element(right[i][j].begin(), right[i][j].end())));
    }
  
  FastReadQuality::set_read_properties(read_width, 0,
				       max_quality_score - 
				       min_quality_score);
  FastReadWC::set_read_properties(read_width, 0, max_quality_score - 
 				  min_quality_score);
  max_match_score = max_mismatches*
    FastReadQuality::quality_to_value(max_quality_score - min_quality_score);
  
  for (size_t i = 0; i < left.size(); ++i) {
    for (size_t j = 0; j < left[i].size(); ++j) {
      const double max_score_left = 
	*max_element(left[i][j].begin(), left[i][j].end());
      const double max_score_right = 
	*max_element(right[i][j].begin(), right[i][j].end());
      for (size_t k = 0; k < left[i][j].size(); ++k) {
	left[i][j][k] = max_score_left - left[i][j][k];
	right[i][j][k] = max_score_right - right[i][j][k];
      }
    }
  }
}


static void
fastq_to_prb(const vector<string> &reads,
	     const vector<vector<double> > &fastq_scores, 
	     vector<vector<vector<double> > > &scores) {
  typedef vector<double> vv;
  scores.resize(reads.size());
  for (size_t i = 0; i < reads.size(); ++i) {
    for (size_t j = 0; j < reads[i].length(); ++j) {
      const double score = 
	solexa_to_error_probability(fastq_scores[i][j]);
      const double other_scores = 
	error_probability_to_solexa(1 - score);
      scores[i].push_back(vv(rmap::alphabet_size, other_scores));
      if (isvalid(reads[i][j]))
	scores[i].back()[base2int(reads[i][j])] = score;
    }
  }
}


// static void
// set_read_width(const vector<string> &reads, size_t &read_width) {
//   if (read_width == 0) {
//     read_width = reads.front().size();
//     for (size_t i = 1; i < reads.size(); ++i)
//       read_width = min(read_width, reads[i].length());
//   }
// }


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
load_read_names(const size_t INPUT_MODE, 
		string reads_file, vector<string> &read_names) {
  vector<string> reads;
  vector<vector<double> > scores;
  if (INPUT_MODE == FASTQ_FILE)
    read_fastq_file(reads_file.c_str(), read_names,
		    reads, scores);
  else read_fasta_file(reads_file.c_str(), read_names, reads);
  reads.clear();
  scores.clear();
  for (vector<string>::iterator i(read_names.begin()); 
       i != read_names.end(); ++i)
    i->erase(i->begin() + i->find_first_of(" \t"), i->end());
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
get_run_mode(const bool VERBOSE, const size_t INPUT_MODE, const bool WILDCARD) {
  size_t RUN_MODE = RUN_MODE_MISMATCH;
  if (WILDCARD) {
    if (INPUT_MODE == FASTA_FILE)
      throw RMAPException("quality score information "
			  "required to use wildcards");
    RUN_MODE = RUN_MODE_WILDCARD;
  }
  else if (INPUT_MODE == FASTA_AND_PRB || INPUT_MODE == FASTQ_FILE)
    RUN_MODE = RUN_MODE_WEIGHT_MATRIX;
  
  if (VERBOSE)
    cerr << "MATCH MODE: "
	 << ((RUN_MODE == RUN_MODE_MISMATCH) ? 
	     "MISMATCH" : ((RUN_MODE == RUN_MODE_WILDCARD) ? 
			   "WILDCARD" : "WEIGHT-MATRIX")) << endl;
  return RUN_MODE;
}  



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//////
////// FOR PAIRED-END READS
static void
split_scores(const size_t read_width,
	     vector<vector<vector<double> > > &scores, 
	     vector<vector<vector<double> > > &left,
	     vector<vector<vector<double> > > &right) {
  typedef vector<vector<double> > vv;
  for (size_t i = 0; i < scores.size(); ++i) {
    left.push_back(vv(scores[i].begin(), scores[i].begin() + read_width));
    right.push_back(vv(scores[i].end() - read_width, scores[i].end()));
    reverse(right.back().begin(), right.back().end());
    for (size_t j = 0; j < read_width; ++j)
      reverse(right.back()[j].begin(), right.back()[j].end());
    scores[i].clear();
  }
  scores.clear();
}


static void
split_reads(const size_t read_width, vector<string> &reads,
	    vector<string> &left, vector<string> &right) {
  for (size_t i = 0; i < reads.size(); ++i) {
    left.push_back(reads[i].substr(0, read_width));
    right.push_back(revcomp(reads[i].substr(reads[i].length() - read_width)));
  }
  reads.clear();
}


static size_t
determine_read_width(const vector<string> &reads, size_t read_width) {
  size_t both_ends = reads.front().size();
  for (size_t i = 1; i < reads.size(); ++i)
    if (both_ends != reads[i].length())
      throw RMAPException("paired-end reads must have uniform width");
  if (read_width == 0)
    read_width = std::floor(both_ends/2);
  else if (read_width > both_ends/2)
    read_width = both_ends - read_width;
  return read_width;
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
	   vector<vector<size_t> > &read_index, vector<size_t> &read_words,
	   size_t &read_width,
	   double &max_match_score) {

  //////////////////////////////////////////////////////////////
  // LOAD THE READS (AS SEQUENCES OR PROBABILITIES) FROM DISK
  if (VERBOSE) cerr << "[LOADING READ SEQUENCES] ";
  vector<string> reads, read_names;
  vector<vector<vector<double> > > scores;

  if (INPUT_MODE == FASTQ_FILE) {
    vector<vector<double> > fastq_scores;
    read_fastq_file(reads_file.c_str(), read_names, reads, fastq_scores);
    read_names.clear();
    fastq_to_prb(reads, fastq_scores, scores);
  }
  else if (INPUT_MODE == FASTA_AND_PRB) {
    read_fasta_file(reads_file.c_str(), read_names, reads);
    // set_read_width(reads, read_width);
    read_names.clear();
    read_prb_file(prb_file.c_str(), scores);
    if (scores.size() != reads.size())
      throw RMAPException("different number of reads in prb and reads files");
  }
  else {
    read_fasta_file(reads_file.c_str(), read_names, reads);
    if (reads.empty())
      throw RMAPException("empty reads file:\"" + reads_file + "\"");
    read_names.clear();
  }
  if (VERBOSE)
    cerr << "[DONE]" << endl
	 << "TOTAL READS: " << reads.size() << endl;


  //////////////////////////////////////////////////////////////
  // SET THE READ WIDTH (IMPORTANT PARAMETER!)
  read_width = determine_read_width(reads, read_width);
  // set_read_width(reads, read_width);
  if (VERBOSE) cerr << "READ WIDTH: " << read_width << endl;
  
  
  vector<string> reads_left, reads_right;
  split_reads(read_width, reads, reads_left, reads_right);
  reads.clear();
  vector<vector<vector<double> > > scores_left, scores_right;
  if (INPUT_MODE > 0) {
    split_scores(read_width, scores, scores_left, scores_right);
    scores.clear();
  }
  
  //////////////////////////////////////////////////////////////
  // BUILD THE INDEX FOR THE READ NAMES, AND CHECK THE READS FOR
  // QUALITY
  if (VERBOSE) cerr << "[CLEANING READS] ";
  vector<size_t> dummy_read_index(reads_left.size(), 1ul);
  partial_sum(dummy_read_index.begin(), dummy_read_index.end(),
	      dummy_read_index.begin());
  transform(dummy_read_index.begin(), dummy_read_index.end(),
	    dummy_read_index.begin(), std::bind2nd(std::minus<size_t>(), 1ul));
  if (INPUT_MODE > 0)
    clean_reads(max_mismatches, read_width, dummy_read_index, 
		reads_left, reads_right, scores_left, scores_right, read_index);
  else clean_reads(max_mismatches, read_width, dummy_read_index, 
		   reads_left, reads_right, read_index);
  dummy_read_index.clear();
  if (VERBOSE)
    cerr << "[DONE]" << endl
	 << "TOTAL READS (HQ): " << read_index.size() << endl;
  
  //////////////////////////////////////////////////////////////
  // CONVERT THE READS INTO FASTREADS
  if (VERBOSE) cerr << "[BUILDING FAST READS] ";
  double max_quality_score = 0;
  if (RUN_MODE == RUN_MODE_WILDCARD) {
    process_score_data(read_width, max_mismatches,
		       max_quality_score, max_match_score, 
		       scores_left, scores_right);
    for (size_t i = 0; i < reads_left.size(); ++i) {
      fast_reads_wc_left.push_back(FastReadWC(scores_left[i]));
      fast_reads_wc_right.push_back(FastReadWC(scores_right[i]));
      scores_left[i].clear();
      scores_right[i].clear();
    }
    scores_left.clear();
    scores_right.clear();
  }
  else if (RUN_MODE == RUN_MODE_WEIGHT_MATRIX) {
    process_score_data(read_width, max_mismatches,
		       max_quality_score, max_match_score, 
		       scores_left, scores_right);
    for (size_t i = 0; i < reads_left.size(); ++i) {
      fast_reads_q_left.push_back(FastReadQuality(scores_left[i]));
      fast_reads_q_right.push_back(FastReadQuality(scores_right[i]));
      scores_left[i].clear();
      scores_right[i].clear();
    }
    scores_left.clear();
    scores_right.clear();
  }
  else {
    FastRead::set_read_width(read_width);
    for (size_t i = 0; i < reads_left.size(); ++i) {
      fast_reads_left.push_back(FastRead(reads_left[i]));
      fast_reads_right.push_back(FastRead(reads_right[i]));
    }
  }
  if (VERBOSE) cerr << "[DONE]" << endl;
  
  //////////////////////////////////////////////////////////////
  // GET THE 2-BIT FORMAT OF THE FIRST PART OF EACH READ
  get_read_words(reads_right, read_words);
  reads_left.clear();
  reads_right.clear();
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
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse("rmappe", "The rmappe mapping tool for "
			   "paired-end Solexa reads",
			   "<fast[a/q]-reads-file>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false , outfile);
    opt_parse.add_opt("chrom", 'c', "FASTA file or dir containing chromosome(s)", 
		      true , chrom_file);
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
    opt_parse.add_opt("wc", 'W', "wildcard matching based on quality scores", 
		      false, WILDCARD);
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

    
    //////////////////////////////////////////////////////////////
    //  CHECK HOW QUALITY SCORES ARE USED
    //
    const size_t INPUT_MODE = get_input_mode(VERBOSE, reads_file, prb_file);
    const size_t RUN_MODE = get_run_mode(VERBOSE, INPUT_MODE, WILDCARD);
    
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
    vector<vector<size_t> > read_index;
    vector<size_t> read_words;
    double max_match_score = 0;
    load_reads(VERBOSE, INPUT_MODE, RUN_MODE, max_mismatches, 
	       reads_file, prb_file, 
	       fast_reads_left, fast_reads_right, 
	       fast_reads_wc_left, fast_reads_wc_right,
	       fast_reads_q_left, fast_reads_q_right, 
	       read_index, read_words, read_width, max_match_score);
    
    
    //////////////////////////////////////////////////////////////
    // INITIALIZE THE SEED STRUCTURES
    //
    vector<size_t> the_seeds;
    load_seeds(VERBOSE, read_width, n_seeds, seed_weight, the_seeds);
    
    
    //////////////////////////////////////////////////////////////
    // INITIALIZE THE STRUCTURES THAT HOLD THE RESULTS
    //
    MultiMapResultPE::init(max_mappings);
    vector<MultiMapResultPE> 
      best_maps(read_words.size(), MultiMapResultPE((RUN_MODE == RUN_MODE_WEIGHT_MATRIX) ?
						    2*max_match_score : 2*max_mismatches));
    

    //////////////////////////////////////////////////////////////
    // THIS IS WHERE THE ACTUAL MAPPING HAPPENS
    //
    vector<size_t> chrom_sizes;
    vector<string> chrom_names;
    vector<size_t> ambigs;
    if (RUN_MODE == RUN_MODE_MISMATCH)
      iterate_over_seeds(VERBOSE, the_seeds, chrom_files, ambigs, 
			 chrom_names, chrom_sizes,
			 fast_reads_left, fast_reads_right, // USE REGULAR FAST READS
			 read_words, read_index,
			 best_maps, max_mismatches, read_width,
			 min_sep, max_sep);
    if (RUN_MODE == RUN_MODE_WILDCARD)
      iterate_over_seeds(VERBOSE, the_seeds, chrom_files, ambigs, 
			 chrom_names, chrom_sizes,
			 fast_reads_wc_left, fast_reads_wc_right, // USE FAST READS FOR WILDCARD MATCHING
			 read_words, read_index,
			 best_maps, max_mismatches, read_width,
			 min_sep, max_sep);
    if (RUN_MODE == RUN_MODE_WEIGHT_MATRIX)
      iterate_over_seeds(VERBOSE, the_seeds, chrom_files, ambigs, 
			 chrom_names, chrom_sizes,
			 fast_reads_q_left, fast_reads_q_right, // USE FAST READS FOR QUALITY MATCHING
			 read_words, read_index,
			 best_maps, max_match_score, read_width,
			 min_sep, max_sep);
    
    
    //////////////////////////////////////////////////////////////
    // LOAD THE NAMES OF READS AGAIN (THEY WILL BE NEEDED)
    //
    vector<string> read_names;
    load_read_names(INPUT_MODE, reads_file, read_names);


    //////////////////////////////////////////////////////////////
    // IF IDENTITIES OF AMBIGUOUS READS ARE DESIRED, WRITE THEM
    //
    if (!ambiguous_file.empty())
      write_non_uniques(ambiguous_file, ambigs, read_names);

    
    //////////////////////////////////////////////////////////////
    // TRANSFORM THE RESULT STRUCTURES INTO BED FORMAT FOR OUTPUT
    //
    vector<GenomicRegion> hits;
    sites_to_regions(read_width, chrom_names, chrom_sizes, read_index, read_names,
		     best_maps, hits);


    //////////////////////////////////////////////////////////////
    // CORRECT SCORES IF WEIGHT MATRIX MATCHING USED
    //
    if (RUN_MODE == RUN_MODE_WEIGHT_MATRIX)
      for (size_t i = 0; i < hits.size(); ++i)
	hits[i].set_score(FastReadQuality::value_to_quality(hits[i].get_score()));
    
    
    //////////////////////////////////////////////////////////////
    // OUTPUT THE RESULTS
    //
    std::ostream* out = (!outfile.empty()) ? 
      new std::ofstream(outfile.c_str()) : &cout;
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
