/*    bedoverlap: a program for finding intervals overlapping members
 *    of another set of intervals
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

#include <fstream>

#include "OptionParser.hpp"
#include "rmap_utils.hpp"
#include "GenomicRegion.hpp"

using std::ofstream;
using std::string;
using std::vector;
using std::ostream;
using std::cout;
using std::endl;
using std::cerr;
using std::numeric_limits;


void
filter_scores(const float lower_bound, const float upper_bound,
	      vector<GenomicRegion> &regions) {
  vector<GenomicRegion> new_regions;
  new_regions.reserve(regions.size()/2);
  for (size_t i = 0; i < regions.size(); ++i) {
    const double score(regions[i].get_score());
    if (score >= lower_bound && score <= upper_bound)
      new_regions.push_back(regions[i]);
  }
  regions.swap(new_regions);
}


void
sift_single_chrom(const vector<GenomicRegion> &other_regions,
		  const bool exclude, const bool full_containment,
		  const vector<GenomicRegion> &regions,
		  vector<GenomicRegion> &good_regions) {
  
  typedef vector<GenomicRegion>::const_iterator region_itr;
  
  if (exclude) {
    if (full_containment) {
      for (size_t i = 0; i < regions.size(); ++i) {
	region_itr closest(find_closest(other_regions, regions[i]));
	if (!closest->contains(regions[i]))
	  good_regions.push_back(regions[i]);
      }
    }
    else {
      for (size_t i = 0; i < regions.size(); ++i) {
	region_itr closest(find_closest(other_regions, regions[i]));
	if (!closest->overlaps(regions[i]))
	  good_regions.push_back(regions[i]);
      }
    }
  }
  else {
    if (full_containment) {
      for (size_t i = 0; i < regions.size(); ++i) {
	region_itr closest(find_closest(other_regions, regions[i]));
	if (closest->contains(regions[i]))
	  good_regions.push_back(regions[i]);
      }
    }
    else {
      for (size_t i = 0; i < regions.size(); ++i) {
	region_itr closest(find_closest(other_regions, regions[i]));
	if (closest->overlaps(regions[i]))
	  good_regions.push_back(regions[i]);
      }
    }
  }
}


void
sift(const vector<GenomicRegion> &other_regions,
     const bool exclude, const bool full_containment,
     vector<GenomicRegion> &regions) {
  
  vector<vector<GenomicRegion> > other_regions_by_chrom;
  separate_chromosomes(other_regions, other_regions_by_chrom);
  
  unordered_map<string, size_t> chrom_lookup;
  for (size_t i = 0; i < other_regions_by_chrom.size(); ++i) 
    chrom_lookup[other_regions_by_chrom[i].front().get_chrom()] = i;
  const vector<GenomicRegion> dummy;
  
  vector<vector<GenomicRegion> > regions_by_chrom;
  separate_chromosomes(regions, regions_by_chrom);
  regions.clear();
  
  vector<GenomicRegion> good_regions;
  for (size_t i = 0; i < regions_by_chrom.size(); ++i) {
    const unordered_map<string, size_t>::const_iterator j = 
      chrom_lookup.find(regions_by_chrom[i].front().get_chrom());
    if (j != chrom_lookup.end()) {
      sift_single_chrom(other_regions_by_chrom[j->second],
			exclude, full_containment, regions_by_chrom[i], 
			good_regions);
    }
    else if (exclude)
      good_regions.insert(good_regions.end(), regions_by_chrom[i].begin(),
			  regions_by_chrom[i].end());
  }
  regions.swap(good_regions);
}


int
main(int argc, const char **argv) {

  /* FILES */
  string target_regions_file;
  string outfile;
  
  bool exclude = false;
  bool full_containment = false;
  bool VERBOSE = false;

  float upper_bound = numeric_limits<float>::max();
  float lower_bound = -numeric_limits<float>::max();

  /****************** GET COMMAND LINE ARGUMENTS ***************************/
  OptionParser opt_parse("mapsifter", "a program for sifting through mapped "
			 "read locations (given in BED format)",
			 "<bed-format-file>");
  opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		    false , outfile);
  opt_parse.add_opt("verbose", 'v', "print more run info", 
		    false , VERBOSE);
  opt_parse.add_opt("exclude", 'e', "exclude contained", 
		    false , exclude);
  opt_parse.add_opt("upper", 'u', "upper bound on scores", 
		    false , upper_bound);
  opt_parse.add_opt("lower", 'l', "lower bound on scores", 
		    false , lower_bound);
  opt_parse.add_opt("target", 't', "target regions file", 
		    false , target_regions_file);
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
  if (leftover_args.size() != 1) {
    cerr << opt_parse.help_message() << endl;
    return EXIT_SUCCESS;
  }
  const string input_file_name = leftover_args.back();
  /**********************************************************************/
  
  try {
    
    if (VERBOSE && full_containment)
      cerr << "WARNING: full containment may not "
	   << "work if target regions overlap" << endl;
    
    vector<GenomicRegion> regions;
    ReadBEDFile(input_file_name, regions);
    
    if (!target_regions_file.empty()) {
      if (!check_sorted(regions)) {
	cerr << "ERROR: regions in \"" << input_file_name
	     << "\" not sorted" << endl;
	return EXIT_FAILURE;
      }
      vector<GenomicRegion> other_regions;
      ReadBEDFile(target_regions_file, other_regions);
      if (!check_sorted(other_regions)) {
	cerr << "ERROR: regions in \"" << target_regions_file
	     << "\" not sorted" << endl;
	return EXIT_FAILURE;
      }
      sift(other_regions, exclude, full_containment, regions);
    }
    
    if (lower_bound != -numeric_limits<float>::max() ||
	upper_bound != numeric_limits<float>::max())
      filter_scores(lower_bound, upper_bound, regions);
    
    ostream* out = (!outfile.empty()) ? 
      new ofstream(outfile.c_str()) : &cout;
    copy(regions.begin(), regions.end(), 
	 std::ostream_iterator<GenomicRegion>(*out, "\n"));
    if (out != &cout) delete out;
    
  }
  catch (RMAPException &e) {
    cerr << "ERROR:\t" << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
