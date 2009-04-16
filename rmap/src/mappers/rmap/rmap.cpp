/*    rmap: a program for mapping Solexa reads
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

#include "FastRead.hpp"
#include "FastReadQuality.hpp"
#include "rmap_os.hpp"
#include "SeedMaker.hpp"
#include "MapResult.hpp"
#include "OptionParser.hpp"

#include <tr1/unordered_map>

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

typedef unordered_multimap<size_t, size_t> SeedHash;
static void
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


template <class T> void
map_reads(const string &chrom, const size_t chrom_id,
	  const size_t profile, const size_t read_width, 
	  const size_t max_diffs, const vector<T> &fast_reads,
	  const SeedHash &seed_hash, const bool strand, 
	  vector<MultiMapResult> &best_maps) {
  
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
      std::pair<SeedHash::const_iterator, SeedHash::const_iterator>
	bucket(seed_hash.equal_range((read_word & profile)));
      if (bucket.first != bucket.second) {
	const SeedHash::const_iterator  limit(bucket.second);
	for (SeedHash::const_iterator to_test(bucket.first); 
	     to_test != limit; ++to_test) {
	  const size_t score = fast_reads[to_test->second].score(fast_read);
	  if (score <= max_diffs) {
	    const vector<MultiMapResult>::iterator 
	      current(best_maps.begin() + to_test->second);
	    if (score <= current->score)
	      current->add(score, chrom_id, 
			   chrom_offset - read_width + 1, strand);
	  }
	}
      }
    }
  }
}


static bool
good_read(const vector<vector<double> > &read) {
  size_t bad_count = 0;
  for (size_t i = 0; i < read.size(); ++i)
    bad_count += (*std::max_element(read[i].begin(), read[i].end()) <= 0);
  return (bad_count <= 2);
}


static void
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


static void
sites_to_regions(const vector<string> &chrom, const vector<size_t> &chrom_sizes, 
		 const vector<string> &reads, const vector<vector<string> > &read_names, 
		 vector<MultiMapResult> &bests, const size_t max_diffs, 
		 vector<GenomicRegion> &hits) {
  for (size_t i = 0; i < bests.size(); ++i)
    if (!bests[i].mr.empty()) {
      sort(bests[i].mr.begin(), bests[i].mr.end());
      for (size_t j = 0; j < bests[i].mr.size(); ++j)
	if (j == 0 || bests[i].mr[j - 1] < bests[i].mr[j]) {
	  const size_t chrom_id = bests[i].mr[j].chrom;
	  const size_t start = bests[i].mr[j].strand ? 
	    bests[i].mr[j].site : 
	    chrom_sizes[chrom_id] - bests[i].mr[j].site - reads[i].length();
	  const size_t end = start + reads[i].length();
	  const size_t score = bests[i].mr[j].score;
	  const char strand = ((bests[i].mr[j].strand) ? '+' : '-');
	  for (size_t k = 0; k < read_names[i].size(); ++k)
	    hits.push_back(GenomicRegion(chrom[chrom_id], start, end,
					 read_names[i][k], score, strand));
	}
    }
}


static char
to_base_symbol(char c) {return (isvalid(c)) ? toupper(c) : 'N';}


static void
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


static void
sites_to_regions(const vector<string> &chrom, const vector<size_t> &chrom_sizes, 
		 const vector<string> &reads, const vector<vector<string> > &read_names, 
		 vector<MultiMapResult> &bests, const size_t min_match_score, 
		 const double max_quality_score,
		 vector<GenomicRegion> &hits) {
  for (size_t i = 0; i < bests.size(); ++i)
    if (!bests[i].mr.empty()) {
      sort(bests[i].mr.begin(), bests[i].mr.end());
      for (size_t j = 0; j < bests[i].mr.size(); ++j)
	if (j == 0 || bests[i].mr[j - 1] < bests[i].mr[j]) {
	  const size_t chrom_id = bests[i].mr[j].chrom;
	  const size_t start = bests[i].mr[j].strand ? 
	    bests[i].mr[j].site : 
	    chrom_sizes[chrom_id] - bests[i].mr[j].site - reads[i].length();
	  const size_t end = start + reads[i].length();
	  const size_t score = bests[i].mr[j].score;
	  const char strand = ((bests[i].mr[j].strand) ? '+' : '-');
	  for (size_t k = 0; k < read_names[i].size(); ++k)
	    hits.push_back(GenomicRegion(chrom[chrom_id], start, end,
					 read_names[i][k], 
					 FastReadQuality::value_to_quality(score), 
					 strand));
	}
    }
}


static void
write_non_uniques(string filename, const vector<pair<string, size_t> > &ambigs) {
  std::ofstream out(filename.c_str());
  for (size_t i = 0; i < ambigs.size(); ++i)
    out << ambigs[i].first << "\t" << ambigs[i].second << endl;
  out.close();
}


template <class T> void
eliminate_ambigs(const size_t max_mismatches,
		 vector<MultiMapResult> &best_maps, vector<vector<string> > &read_names,
		 vector<string> &reads, vector<pair<string, size_t> > &ambigs,
		 vector<T> &fast_reads) {

  size_t j = 0;
  for (size_t i = 0; i < best_maps.size(); ++i) {
    best_maps[i].collapse();
    if (best_maps[i].ambiguous() && best_maps[i].score <= max_mismatches)
      for (size_t k = 0; k < read_names[i].size(); ++k)
	ambigs.push_back(make_pair(read_names[i][k], best_maps[i].score));
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


static void
process_score_data(const size_t read_width, const size_t max_mismatches,
		   double &max_quality_score, double &max_match_score,
		   vector<vector<vector<double> > > &scores) {
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
  
  FastReadQuality::set_read_properties(read_width, 0, max_quality_score - 
				       min_quality_score);
  
  max_match_score = max_mismatches*
    FastReadQuality::quality_to_value(max_quality_score - min_quality_score);
  
  for (size_t i = 0; i < scores.size(); ++i) {
    for (size_t j = 0; j < scores[i].size(); ++j) {
      const double max_score = *max_element(scores[i][j].begin(), 
					    scores[i][j].end());
      for (size_t k = 0; k < scores[i][j].size(); ++k)
	scores[i][j][k] = max_score - scores[i][j][k];
    }
  }
}


static void
fastq_to_prb(const vector<string> &reads,
	     const vector<vector<double> > &fastq_scores, 
	     vector<vector<vector<double> > > &scores) {
  typedef vector<double> vv;
  for (size_t i = 0; i < reads.size(); ++i) {
    scores.push_back(vector<vv>(reads[i].length(),
				vv(rmap::alphabet_size, -40)));
    for (size_t j = 0; j < reads[i].size(); ++j)
      scores[i][j][base2int(reads[i][j])] = 
	10.0*log(1 + pow(10.0, (fastq_scores[i][j] - 64)/10.0))/log(10.0);
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
    size_t max_mappings = 1;
    
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
    opt_parse.add_opt("prb", 'p', "file with quality scores (prb format)", 
		      false, prb_file);
    opt_parse.add_opt("seeds", 'S', "number of seeds", true , n_seeds);
    opt_parse.add_opt("hit", 'h', "weight of hit", true , seed_weight);
    opt_parse.add_opt("width", 'w', "width of reads", false, read_width);
    opt_parse.add_opt("mismatch", 'm', "maximum allowed mismatches", 
		      true , max_mismatches);
    opt_parse.add_opt("ambiguous", 'a', "file to write names of ambiguously mapped reads", 
		      false , ambiguous_file);
    opt_parse.add_opt("max-map", 'M', "maximum allowed mappings for a read", 
		      false, max_mappings);
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

    const bool FASTQ_READS = (is_fastq(reads_file));
    const bool USING_QUALITY = (!prb_file.empty() || FASTQ_READS);
    
    if (VERBOSE && USING_QUALITY)
      cerr << "USING QUALITY SCORES" << endl;
    
    //// GET THE READS
    if (VERBOSE)
      cerr << "READING DATA" << endl;
    vector<string> input_read_names, reads;
    vector<vector<vector<double> > > scores;
    if (FASTQ_READS) {
      vector<vector<double> > fastq_scores;
      read_fastq_file(reads_file.c_str(), input_read_names, reads,
		      fastq_scores);
      fastq_to_prb(reads, fastq_scores, scores);
    }
    else {
      read_fasta_file(reads_file.c_str(), input_read_names, reads);
      if (USING_QUALITY) {
	if (VERBOSE)
	  cerr << "READING QUALITY SCORES" << endl;
	read_prb_file(prb_file.c_str(), scores);
	if (scores.size() != reads.size())
	  throw RMAPException("different number of reads in prb and reads files");
      }
    }
    
    if (VERBOSE)
      cerr << "TOTAL READS:    " << reads.size() << endl;
    
    if (read_width == 0) {
      read_width = reads.front().size();
      for (size_t i = 1; i < reads.size(); ++i)
	read_width = min(read_width, reads[i].length());
    }
    
    if (VERBOSE)
      cerr << "READ WIDTH:     " << read_width << endl;
    
    // prepare the reads
    if (VERBOSE)
      cerr << "CHECKING READ QUALITY" << endl;
    vector<vector<string> > read_names;
    if (USING_QUALITY)
      clean_reads(max_mismatches, read_width, input_read_names, 
		  reads, scores, read_names);
    else clean_reads(max_mismatches, read_width, input_read_names, 
		     reads, read_names);
    if (VERBOSE)
      cerr << "READS AFTER QC: " << reads.size() << endl;
    
    // convert the reads into FastReads
    vector<FastRead> fast_reads;
    vector<FastReadQuality> fast_reads_q;
    double max_match_score = 0, max_quality_score = 0;
    if (USING_QUALITY) {
      process_score_data(read_width, max_mismatches,
			 max_quality_score, max_match_score, scores);
      if (VERBOSE)
	cerr << "MAX_MATCH_SCORE=" << max_match_score << endl
	     << "MAX_QUALITY_SCORE=" << max_quality_score << endl;
      for (size_t i = 0; i < reads.size(); ++i) {
	fast_reads_q.push_back(FastReadQuality(scores[i]));
	scores[i].clear();
      }
      scores.clear();
    }
    else {
      FastRead::set_read_width(read_width);
      for (size_t i = 0; i < reads.size(); ++i)
	fast_reads.push_back(FastRead(reads[i]));
    }
    
    if (VERBOSE)
      cerr << "IDENTIFYING CHROMS" << endl;
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

    MultiMapResult::init(max_mappings);
    vector<MultiMapResult> best_maps(reads.size(), 
				     MultiMapResult((USING_QUALITY) ? 
						    max_match_score : 
						    max_mismatches));
    
    if (VERBOSE)
      cerr << endl << "scanning chromosomes:" << endl;
    
    for (size_t j = 0; j < the_seeds.size() && !reads.empty(); ++j) {
      for (size_t i = 0; i < chrom_files.size() && !reads.empty(); ++i) {
	if (VERBOSE)
	  cerr << "[SEED:" << j + 1 << "/" << the_seeds.size() << "] [SEEDING] ";
	
	SeedHash seed_hash;
	get_read_matches(the_seeds[j], reads, seed_hash);
	
	vector<string> chrom_names, chrom;
	if (VERBOSE) cerr << "[LOADING CHROM] ";
	read_fasta_file(chrom_files[i].c_str(), chrom_names, chrom);
	if (j == 0)
	  chrom_sizes[i] = chrom.front().length();
	
	transform(chrom.front().begin(), chrom.front().end(), chrom.front().begin(), 
		  std::ptr_fun(&toupper));
	
	if (VERBOSE)
	  cerr << "[SCANNING=" << chrom_names.front() << "] ";
	
	if (USING_QUALITY)
	  map_reads(chrom.front(), i, the_seeds[j], read_width, max_match_score,
		    fast_reads_q, seed_hash, true, best_maps);
	else map_reads(chrom.front(), i, the_seeds[j], read_width, max_mismatches,
		       fast_reads, seed_hash, true, best_maps);
	
	revcomp_inplace(chrom.front());
	
	if (USING_QUALITY)
	  map_reads(chrom.front(), i, the_seeds[j], read_width, max_match_score,
		    fast_reads_q, seed_hash, false, best_maps);
	else map_reads(chrom.front(), i, the_seeds[j], read_width, max_mismatches,
		       fast_reads, seed_hash, false, best_maps);
	
	if (VERBOSE)
	  cerr << "[CLEANING=" << chrom_names.front() << "] ";
	if (USING_QUALITY)
	  eliminate_ambigs(0, best_maps, read_names, reads, ambigs, fast_reads_q);
	else
	  eliminate_ambigs(0, best_maps, read_names, reads, ambigs, fast_reads);
	if (VERBOSE)
	  cerr << "[AMBIG=" << ambigs.size() << "] " << endl;
      }
    }

    if (VERBOSE)
      cerr << endl << "[eliminating ambiguous reads...";
    
    if (USING_QUALITY)
      eliminate_ambigs(max_match_score, best_maps, read_names, 
		       reads, ambigs, fast_reads_q);
    else eliminate_ambigs(max_mismatches, best_maps, read_names, 
			  reads, ambigs, fast_reads);
    
    if (!ambiguous_file.empty())
      write_non_uniques(ambiguous_file, ambigs);
    if (VERBOSE)
      cerr << "done]" << endl;
    
    // transform best matches into BED format
    vector<GenomicRegion> hits;

    if (USING_QUALITY)
      sites_to_regions(chrom_names, chrom_sizes, reads, read_names,
		       best_maps, max_match_score, max_quality_score, hits);
    else sites_to_regions(chrom_names, chrom_sizes, reads, read_names,
			  best_maps, max_mismatches, hits);
    
    // Output the results
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
