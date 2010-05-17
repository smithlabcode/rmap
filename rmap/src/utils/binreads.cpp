/*    binreads: 
 *
 *    Copyright (C) 2010 University of Southern California and
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


static void
BinReads(const vector<GenomicRegion> &reads,
	 const size_t region_start, const size_t region_end,
	 const size_t bin_size, std::ostream &out) {
  const string chrom_name(reads.front().get_chrom());
  size_t read_idx = 0;
  for (size_t i = region_start; i < region_end; i += bin_size) {
    size_t counts = 0;
    while (read_idx < reads.size() && reads[read_idx].get_start() < i + bin_size) {
      if (reads[read_idx].get_start() >= i)
	++counts;
      ++read_idx;
    }
    if (counts > 0)
      out << chrom_name << '\t' << i << '\t' << i + bin_size << '\t' << counts << '\n';
  }
}


int
main(int argc, const char **argv) {

  /* FILES */
  string outfile;
  string chroms_file_name;
  bool VERBOSE = false;
  size_t bin_size = 100;

  /****************** GET COMMAND LINE ARGUMENTS ***************************/
  OptionParser opt_parse("binreads", "",
			 "<bed-format-file>");
  opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		    false , outfile);
  opt_parse.add_opt("chrom", 'c', "chrom sizes file", 
		    true , chroms_file_name);
  opt_parse.add_opt("verbose", 'v', "print more run info", 
		    false , VERBOSE);
  opt_parse.add_opt("bin", 'b', "size of bins", false , bin_size);
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
    
    vector<GenomicRegion> regions;
    ReadBEDFile(input_file_name, regions);
    if (!check_sorted(regions)) {
      cerr << "ERROR: regions in \"" << input_file_name
	   << "\" not sorted" << endl;
      return EXIT_FAILURE;
    }
    
    vector<GenomicRegion> chroms;
    ReadBEDFile(chroms_file_name, chroms);
    unordered_map<string, size_t> chrom_sizes;
    for (size_t i = 0; i < chroms.size(); ++i)
      chrom_sizes[chroms[i].get_chrom()] = chroms[i].get_width();
    
    vector<vector<GenomicRegion> > separated_by_chrom;
    separate_chromosomes(regions, separated_by_chrom);

    ostream* out = (!outfile.empty()) ? 
      new ofstream(outfile.c_str()) : &cout;
    for (size_t i = 0; i < separated_by_chrom.size(); ++i) {
      const string chrom_name(separated_by_chrom[i].front().get_chrom());
      if (VERBOSE)
	cerr << "[BINNING=" << chrom_name << "]" << endl;
      const size_t chrom_size = chrom_sizes[separated_by_chrom[i].front().get_chrom()];
      BinReads(separated_by_chrom[i], 0, chrom_size, bin_size, *out);
    }
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
