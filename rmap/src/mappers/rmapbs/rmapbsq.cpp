/*    rmapbsq: a program for mapping bisulfite treated Solexa reads
 *    using quality scores
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

#include <string>
#include <vector>
#include <iostream>
#include <iterator>
#include <fstream>
#include <cmath>

#include "BisulfiteFastReadQuality.hpp"
#include <rmap_os.hpp>
#include <SeedMaker.hpp>
#include <MapResult.hpp>
#include "OptionParser.hpp"

#include <tr1/unordered_map>

using std::tr1::unordered_multimap;

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

typedef unordered_multimap<size_t, size_t> SeedHash;
void
get_read_matches(const size_t the_seed, const vector<string> &reads,
		 SeedHash &seed_hash) {
  for (size_t i = 0; i < reads.size(); ++i) {
    const size_t truncate_to = min(reads[i].length(), SeedMaker::max_seed_part);
    string s(reads[i].substr(0, truncate_to));
    if (SeedMaker::valid_seed(the_seed, s)) {
      std::replace(s.begin(), s.end(), 'C', 'T');
      const size_t the_key = (the_seed & SeedMaker::make_read_word(s));
      seed_hash.insert(SeedHash::value_type(the_key, i));
    }
  }
}

void
map_reads(const string &chrom, const size_t chrom_id,
	  const size_t profile, const size_t read_width, const size_t min_match_score, 
	  const vector<BisulfiteFastReadQuality> &fast_reads, 
	  const SeedHash &seed_hash, vector<MapResult> &best_maps, bool strand) {
  
  MASK_t bad_bases = rmap_bits::all_ones;
  MASK_t read_word = 0;
  BisulfiteFastReadQuality fast_read;
  
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
      std::pair<SeedHash::const_iterator, SeedHash::const_iterator>
	bucket(seed_hash.equal_range((read_word & profile)));
      if (bucket.first != bucket.second) {
	const SeedHash::const_iterator  limit(bucket.second);
	for (SeedHash::const_iterator to_test(bucket.first); to_test != limit; ++to_test) {
	  const size_t score = fast_reads[to_test->second].score(fast_read);
	  if (score <= min_match_score) {
	    const vector<MapResult>::iterator current(best_maps.begin() + to_test->second);
	    if (score <= current->score)
	      current->set(score, chrom_id, chrom_offset - read_width + 1, strand);
	  }
	}
      }
    }
  }
}


bool
good_read(const vector<vector<double> > &read) {
  size_t bad_count = 0;
  for (size_t i = 0; i < read.size(); ++i)
    bad_count += (*std::max_element(read[i].begin(), read[i].end()) <= 0);
  return (bad_count <= 2);
}

void
clean_reads(const size_t min_match_score, const size_t read_width,
	    const vector<string> &input_read_names,
	    vector<string> &reads, vector<vector<vector<double> > > &scores,
	    vector<vector<string> > &read_names) {
  vector<vector<vector<double> > > good_scores;
  vector<string> good_reads;
  for (size_t i = 0; i < reads.size(); ++i) {
    if (good_read(scores[i])) {
      if (reads[i].length() > read_width) {
	reads[i] = reads[i].substr(0, read_width);
	scores[i].erase(scores[i].begin() + read_width, scores[i].end());
      }
      read_names.push_back(vector<string>(1, input_read_names[i]));
      good_reads.push_back(reads[i]);
      good_scores.push_back(scores[i]);
    }
  }
  reads.swap(good_reads);
  scores.swap(good_scores);
}

void
sites_to_regions(const vector<string> &chrom, const vector<size_t> &chrom_sizes, 
		 const vector<string> &reads, const vector<vector<string> > &read_names, 
		 const vector<MapResult> &best_maps, const size_t min_match_score, 
		 const double max_quality_score,
		 vector<GenomicRegion> &hits) {
  for (size_t i = 0; i < best_maps.size(); ++i)
    if (best_maps[i].unique && best_maps[i].score <= min_match_score) {
      const size_t start = best_maps[i].strand ? best_maps[i].site : 
	chrom_sizes[best_maps[i].chrom] - best_maps[i].site - reads[i].length();
      const size_t end = start + reads[i].length();
      const size_t chrom_id = best_maps[i].chrom;
      const size_t score = best_maps[i].score;
      const char strand = (best_maps[i].strand) ? '+' : '-';
      for (size_t j = 0; j < read_names[i].size(); ++j)
	hits.push_back(GenomicRegion(chrom[chrom_id], start, end, 
				     read_names[i][j],
				     BisulfiteFastReadQuality::value_to_quality(score), 
				     strand));
    }
}

void
write_non_uniques(string filename, const vector<pair<string, double> > &ambigs) {
  ofstream out(filename.c_str());
  for (size_t i = 0; i < ambigs.size(); ++i)
    out << ambigs[i].first << "\t" << ambigs[i].second << endl;
  out.close();
}

void
eliminate_ambigs(const size_t min_match_score,
		 vector<MapResult> &best_maps, vector<vector<string> > &read_names,
		 vector<string> &reads, vector<pair<string, double> > &ambigs,
		 vector<BisulfiteFastReadQuality> &fast_reads) {
  size_t j = 0;
  for (size_t i = 0; i < best_maps.size(); ++i) {
    if (!best_maps[i].unique && best_maps[i].score <= min_match_score)
      for (size_t j = 0; j < read_names[i].size(); ++j)
	ambigs.push_back(make_pair(read_names[i][j], 
				   BisulfiteFastReadQuality::value_to_quality(best_maps[i].score)));
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

void
process_score_data(const vector<vector<vector<double> > > &scores, 
		   const size_t read_width, const size_t max_mismatches,
		   double &max_quality_score, double &max_match_score) {
  double min_quality_score = numeric_limits<double>::max();
  for (size_t i = 0; i < scores.size(); ++i)
    for (size_t j = 0; j < scores[i].size(); ++j) {
      max_quality_score = max(max_quality_score, 
			      *max_element(scores[i][j].begin(), 
					   scores[i][j].end()));
      min_quality_score = min(min_quality_score, 
			      *min_element(scores[i][j].begin(), 
					   scores[i][j].end()));
    }
  
  BisulfiteFastReadQuality::set_read_properties(read_width, 0, max_quality_score - min_quality_score);
  const double max_final_score = max_mismatches*(max_quality_score - min_quality_score);
  max_match_score = max_mismatches*
    BisulfiteFastReadQuality::quality_to_value(max_final_score/max_mismatches);
}

void
adjust_prb_matrices(vector<vector<vector<double> > > &scores) {
  for (size_t i = 0; i < scores.size(); ++i) {
    for (size_t j = 0; j < scores[i].size(); ++j) {
      const double max_score = *max_element(scores[i][j].begin(), scores[i][j].end());
      for (size_t k = 0; k < scores[i][j].size(); ++k)
	scores[i][j][k] = max_score - scores[i][j][k];
      scores[i][j][base2int('C')] = min(scores[i][j][base2int('C')],
					scores[i][j][base2int('T')]);
    }
  }
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
    
    size_t n_seeds = 0;
    size_t seed_weight = 0;
    size_t read_width = 0;
    size_t max_mismatches = 0;
    
    bool VERBOSE = false;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse("rmapbsq", "The rmapbsq mapping tool for bisulfite "
			   "treated Solexa reads making use of base-call "
			   "quality scores");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false , outfile);
    opt_parse.add_opt("chrom", 'c', "FASTA file or dir containing chromosome(s)",
		      true , chrom_file);
    opt_parse.add_opt("suffix", 's', "suffix of FASTA files (assumes -c "
		      "indicates dir)", false , fasta_suffix);
    opt_parse.add_opt("filenames", 'F', "file listing names of chromosome files",
		      false , filenames_file);
    opt_parse.add_opt("prb", 'p', "file with quality scores (prb format)", true, prb_file);
    opt_parse.add_opt("seeds", 'S', "number of seeds", true , n_seeds);
    opt_parse.add_opt("hit", 'h', "weight of hit", true , seed_weight);
    opt_parse.add_opt("width", 'w', "width of reads", true , read_width);
    opt_parse.add_opt("mismatch", 'm', "maximum allowed mismatches", true , max_mismatches);
    opt_parse.add_opt("ambiguous", 'a', "file to write names of ambiguously "
		      "mapped reads", false , ambiguous_file);
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

    if (read_width == 0)
      cerr << "ERROR: must specify read width" << endl;
    
    // Get the reads
    vector<string> input_read_names, reads;
    read_fasta_file(reads_file.c_str(), input_read_names, reads);

    vector<vector<vector<double> > > scores;
    read_prb_file(prb_file, scores);

    if (VERBOSE)
      cerr << "TOTAL READS:    " << reads.size() << endl;

    // prepare the reads
    vector<vector<string> > read_names;
    clean_reads(max_mismatches, read_width, input_read_names, reads, scores, read_names);
    if (VERBOSE)
      cerr << "READS AFTER QC: " << reads.size() << endl;
    
    double max_match_score = 0, max_quality_score = 0;
    process_score_data(scores, read_width, max_mismatches, 
		       max_quality_score, max_match_score);
    adjust_prb_matrices(scores);
    
    if (VERBOSE) {
      cerr << "MAX_MATCH_SCORE=" << max_match_score << endl;
      cerr << "MAX_QUALITY_SCORE=" << max_quality_score << endl;
    }
    
    if (scores.size() != reads.size()) {
      cerr << "different number of reads in prb and reads files" << endl;
      return EXIT_FAILURE;
    }

    // convert the reads into word pairs
    // BisulfiteFastReadQuality::set_read_width(read_width);
    vector<BisulfiteFastReadQuality> fast_reads;
    for (size_t i = 0; i < reads.size(); ++i) {
      fast_reads.push_back(BisulfiteFastReadQuality(scores[i]));
      scores[i].clear();
    }
    scores.clear();
    
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
    vector<pair<string, double> > ambigs;
    vector<MapResult> best_maps(reads.size(), MapResult(max_match_score + 1));
    
    for (size_t j = 0; j < the_seeds.size(); ++j) {
      
      for (size_t i = 0; i < chrom_files.size() && !reads.empty(); ++i) {
	if (VERBOSE)
	  cerr << "[SEED:" << j + 1 << "/" 
	       << the_seeds.size() << "] [SEEDING] ";
	
	SeedHash seed_hash;
	get_read_matches(the_seeds[j], reads, seed_hash);
	
	vector<string> chrom_names, chrom;
	if (VERBOSE) cerr << "[LOADING_CHROM] ";
	read_fasta_file(chrom_files[i].c_str(), chrom_names, chrom);
	if (j == 0)
	  chrom_sizes[i] = chrom.front().length();
	
	transform(chrom.front().begin(), chrom.front().end(), chrom.front().begin(), 
		  std::ptr_fun(&toupper));
	if (VERBOSE) cerr << "[SCANNING=" << chrom_names.front() << "] ";
	map_reads(chrom.front(), i, the_seeds[j], read_width, max_match_score,
		  fast_reads, seed_hash, best_maps, true);
	
	revcomp_inplace(chrom.front());
 	map_reads(chrom.front(), i, the_seeds[j], read_width, max_match_score, 
		  fast_reads, seed_hash, best_maps, false);
	
	if (VERBOSE) cerr << "[CLEANING=" << chrom_names.front() << "] ";
	eliminate_ambigs(0, best_maps, read_names, reads, ambigs, fast_reads);
	if (VERBOSE) cerr << "[AMBIG=" << ambigs.size() << "] " << endl;
      }
    }
    
    if (VERBOSE)
      cerr << endl << "[eliminating ambiguous reads...";

    eliminate_ambigs(max_match_score, best_maps, read_names, reads, ambigs, fast_reads);
    
    if (!ambiguous_file.empty())
      write_non_uniques(ambiguous_file, ambigs);
    if (VERBOSE)
      cerr << " done]" << endl;
    
    // Transform best matches into BED format
    vector<GenomicRegion> hits;
    sites_to_regions(chrom_names, chrom_sizes, reads, read_names, 
		     best_maps, max_match_score, max_quality_score, hits);
    
    // Output the results
    ostream* out = (!outfile.empty()) ? new ofstream(outfile.c_str()) : &cout;
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
