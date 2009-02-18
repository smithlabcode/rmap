/*    bedoverlap: a program for finding intervals overlapping members
 *    of another set of intervals
 *
 *    Copyright (C) 2009 Cold Spring Harbor Laboratory, 
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

#include <fstream>

#include "OptionParser.hpp"
#include "rmap_utils.hpp"
#include "GenomicRegion.hpp"

using std::string;
using std::vector;
using std::ostream;
using std::cout;
using std::endl;
using std::cerr;

typedef GenomicRegion* GenomicRegionPointer;

struct region_pointer_less {
  bool operator()(const GenomicRegionPointer a, 
		  const GenomicRegionPointer b) const {
    return (*a) < (*b);
  }
};

void
sort_regions(vector<GenomicRegion> &regions) {
  vector<GenomicRegionPointer> sorter;
  for (vector<GenomicRegion>::iterator i = regions.begin(); 
       i != regions.end(); ++i) sorter.push_back(&(*i));
  sort(sorter.begin(), sorter.end(), region_pointer_less());
  
  vector<GenomicRegion> r;
  r.reserve(regions.size());
  for (vector<GenomicRegionPointer>::const_iterator i(sorter.begin());
       i != sorter.end(); ++i)
    r.push_back(*(*i));
  r.swap(regions);
}

int main(int argc, const char **argv) {

  try {
    /* FILES */
    string outfile;
    
    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse("bedoverlap", "a program for finding intervals "
			   "overlapping members of another set of intervals");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false , outfile);
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
    if (leftover_args.size() != 2) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string containing_regions_file = leftover_args.front();
    const string input_file_name = leftover_args.back();
    /**********************************************************************/

    vector<GenomicRegion> big_regions;
    ReadBEDFile(containing_regions_file, big_regions);
    
    if (!check_sorted(big_regions)) {
      cerr << "ERROR: regions in \"" << containing_regions_file
	   << "\" not sorted" << endl;
      return EXIT_FAILURE;
    }
    
    vector<GenomicRegion> regions;
    ReadBEDFile(input_file_name, regions);
    
    if (!check_sorted(regions)) {
      cerr << "ERROR: input regions not sorted" << endl;
      return EXIT_FAILURE;
    }
    
    vector<GenomicRegion> contained;
    genomic_region_intersection(big_regions, regions, contained);
    regions.clear();
    regions.swap(contained);
    
    ostream* out = (outfile.empty()) ? &cout : new std::ofstream(outfile.c_str());
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
