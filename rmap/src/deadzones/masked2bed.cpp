/*    masked2bed: program for obtaining the locations of deadzones
 *    from a genome masked at dead bases
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

#include "rmap_utils.hpp"
#include "rmap_os.hpp"
#include "GenomicRegion.hpp"
#include "OptionParser.hpp"

using std::ofstream;
using std::string;
using std::vector;
using std::ostream;
using std::cout;
using std::endl;
using std::cerr;

int
main(int argc, const char **argv) {

  try {
    
    string outfile;
    bool VERBOSE = false;
    char mask_char = 'N';
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse("masked2bed", "program for obtaining the locations of "
			   "deadzones from a genome masked at dead bases",
			   "<fasta-chrom-file>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false, outfile);
    opt_parse.add_opt("char", 'c', "masking character", false, mask_char);
    opt_parse.add_opt("verbose", 'v', "print more run information", 
		      false, VERBOSE);
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
    const string seqfile = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/
    
    vector<string> sequences, names;
    read_fasta_file(seqfile.c_str(), names, sequences);
    
    ostream* out = (!outfile.empty()) ? new ofstream(outfile.c_str()) : &cout;
    for (size_t i = 0; i < sequences.size(); ++i) {
      bool in_ambig = false;
      size_t ambig_start = 0;
      
      string chrom;
      size_t start = 0, end = 0;
      parse_region_name(names[i], chrom, start, end);
      
      for (size_t j = 0; j < sequences[i].length(); ++j) {
	if (sequences[i][j] == mask_char) {
	  if (!in_ambig) {
	    in_ambig = true;
	    ambig_start = j;
	  }
	}
	else {
	  if (in_ambig) {
	    *out << SimpleGenomicRegion(chrom, start + 
					ambig_start, start + j) << endl;
	    in_ambig = false;
	  }
	}
      }
    }
    if (out != &cout) delete out;
  }
  catch (RMAPException &e) {
    cerr << "ERROR: " << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
