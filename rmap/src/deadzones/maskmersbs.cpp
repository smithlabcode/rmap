/*    maskmersbs: program for masking sequences at k-mers that appear in
 *    a given bisulfite treated file (special binary format)
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
#include "OptionParser.hpp"
#include "deadzone_utils.hpp"
#include "bisulfite_utils.hpp"

#include <tr1/unordered_set>

using std::tr1::unordered_set;

using std::ofstream;
using std::ostringstream;
using std::ostream_iterator;
using std::string;
using std::vector;
using std::ostream;
using std::cout;
using std::endl;
using std::cerr;

void
read_ambigs(string filename, unordered_set<long_index> &ambigs) {
  std::ifstream fin(filename.c_str(), std::ios::binary);
  if (!fin) {
    const string message = "ERROR:\n"
      "cannot open table file: " + filename;
    throw RMAPException(message);
  }
  
  size_t begin_pos = fin.tellg();
  fin.seekg(0, std::ios_base::end);
  size_t end_pos = fin.tellg();
  fin.seekg(0, std::ios_base::beg);
  const size_t filesize = end_pos - begin_pos;
  
  vector<long_index> tmp(filesize/sizeof(long_index));
  fin.read((char*)&tmp.front(), filesize);
  
  ambigs.insert(tmp.begin(), tmp.end());
  fin.close();
}

int
main(int argc, const char **argv) {

  try {

    // Parameter variables
    size_t kmer = 0;
    string seqfile;
    string outfile;
    char mask_char = 'N';
    bool VERBOSE = false;
    bool REVCOMP = false;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse("maskmersbs", "program for masking sequences "
			   " at k-mers that appear in a given bisulfite "
			   "treated file (special binary format) ",
			   "<kmers-binary-format>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false, outfile);
    opt_parse.add_opt("char", 'c', "masking character", false, mask_char);
    opt_parse.add_opt("kmer", 'k', "Width of k-mers", true, kmer);
    opt_parse.add_opt("seq", 's', "File of sequences to mask", true, seqfile);
    opt_parse.add_opt("rev", 'r', "mask reverse complement", 
		      false, REVCOMP);
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
    const vector<string> filenames(leftover_args);
    /****************** END COMMAND LINE OPTIONS *****************/
    
    if (VERBOSE)
      cerr << seqfile << endl;
    vector<string> names;
    vector<string> sequences;
    read_fasta_file(seqfile.c_str(), names, sequences);
    
    vector<string> to_mask(sequences);

    for (size_t f = 0; f < filenames.size(); ++f) {
      const string current_filename(filenames[f]);
      unordered_set<long_index> ambigs;
      read_ambigs(current_filename, ambigs);
      if (VERBOSE)
	cerr << "\r" << current_filename << "\t" << f + 1 << "/" 
	     << filenames.size() << endl;
      
      for (size_t i = 0; i < sequences.size(); ++i) {
	string treated(sequences[i]);
	if (REVCOMP) {
	  revcomp_inplace(treated);
	  bisulfite_treatment(treated);
	  const string::const_iterator lim(treated.end() - kmer + 1);
	  const string::const_iterator seq_start(treated.begin());
	  for (string::const_iterator j(treated.begin()); j < lim; ++j) {
	    bool good = true;
	    for (size_t k = 0; k < static_cast<size_t>(kmer) && good; ++k)
	      if (!isvalid(*(j + k))) {
		good = false;
		j += k;
	      }
	    if (good) {
	      const long_index key(j, j + kmer);
	      if (ambigs.find(key) != ambigs.end())
		to_mask[i][lim - j - 1] = mask_char;
	    }
	  }
	}
	else {
	  bisulfite_treatment(treated);
	  const string::const_iterator lim(treated.end() - kmer + 1);
	  const string::const_iterator seq_start(treated.begin());
	  for (string::const_iterator j(treated.begin()); j < lim; ++j) {
	    bool good = true;
	    for (size_t k = 0; k < static_cast<size_t>(kmer) && good; ++k)
	      if (!isvalid(*(j + k))) {
		good = false;
		j += k;
	      }
	    if (good) {
	      const long_index key(j, j + kmer);
	      if (ambigs.find(key) != ambigs.end())
		to_mask[i][j - seq_start] = mask_char;
	    }
	  }
	}
      }
    }

    ostream* out = (!outfile.empty()) ? new ofstream(outfile.c_str()) : &cout;
    for (size_t i = 0; i < to_mask.size(); ++i) {
      *out << ">" << names[i] << endl;
      for (size_t j = 0; j < to_mask[i].length(); j += 50)
	*out << to_mask[i].substr(j, 50) << endl;
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

    
