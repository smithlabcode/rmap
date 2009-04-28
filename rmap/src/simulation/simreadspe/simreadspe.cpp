/*    simreadspe: a program for simulating paired-end Solexa reads to
 *    test rmappe
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
#include "RNG.hpp"
#include "sim_utils.hpp"
#include "OptionParser.hpp"

#include <string>
#include <vector>
#include <iostream>
#include <fstream>

using std::ostream_iterator;
using std::numeric_limits;
using std::string;
using std::vector;
using std::endl;
using std::cerr;

void
simreadspe(const Runif &rng,
	   const size_t n_reads, const size_t read_width, 
	   const size_t separation, const size_t max_errors,
	   const string &name, const string &sequence,
	   vector<string> &read_names, vector<string> &reads) {
  
  const size_t lim = sequence.length() - separation - read_width + 1;
  
  for (size_t i = 0; i < n_reads; ++i) {
    bool found_one = false;
    while (!found_one) {          
      const size_t start = rng.runif(0ul, lim);      
      bool valid = true;
      for (size_t j = start; j < start + read_width && valid; ++j)
	valid = isvalid(sequence[j]);
      
      for (size_t j = start + separation; 
	   j < start + separation + read_width && valid; ++j)
	valid = isvalid(sequence[j]);
      
      if (valid) {
	string seq_left(sequence.substr(start, read_width));
	string seq_right(sequence.substr(start + separation, 
					 read_width));
	const bool rc = (rng.runif(0.0,1.0) > 0.5);
	if (rc) {
	  seq_left = revcomp(seq_left);	 
	  seq_right = revcomp(seq_right);
	  seq_left.swap(seq_right);
	}
	
	string error_log_left;
	add_sequencing_errors(rng, max_errors, seq_left, error_log_left);
	
	string error_log_right;
	add_sequencing_errors(rng, max_errors, seq_right, error_log_right);
	
	const string read_name(name + ":" + toa(start) + "-" + 
			       toa(start + separation + read_width) + "_" + toa(!rc) + "_" +
			       error_log_left + "_" + 
			       toa(count(error_log_left.begin(), 
					 error_log_left.end(), '1')) + "_" +
			       error_log_right + "_" + 
			       toa(count(error_log_right.begin(), 
					 error_log_right.end(), '1')));
	reads.push_back(seq_left + seq_right);
	read_names.push_back(read_name);
	found_one = true;
      }
    }
  }
}

int
main(int argc, const char **argv) {

  try {
    
    string seqfile;
    string outfile;
    size_t n_reads = 1000;
    size_t read_width = 25;
    size_t max_errors = 0;
    size_t random_number_seed = numeric_limits<size_t>::max();
    size_t separation = 200;

    bool VERBOSE = false;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse("simreadspe", "program for generating simulated "
			   "paired-end reads", "<fasta-chrom-files>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false , outfile);
    opt_parse.add_opt("reads", 'n', "number of reads to simulate", 
		      false, n_reads);
    opt_parse.add_opt("sep", 's', "separation between ends", 
		      false, separation);
    opt_parse.add_opt("width", 'w', "width of reads to simulate", 
		      false, read_width);
    opt_parse.add_opt("err", 'e', "maximum number of simulated sequencing errors", 
		      false, max_errors);
    opt_parse.add_opt("seed", 'S', "random number seed", false, random_number_seed);
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
    vector<string> filenames(leftover_args);
    /****************** END COMMAND LINE OPTIONS *****************/

    const Runif rng(random_number_seed);

    vector<string> reads, read_names;
    vector<size_t> filesizes;
    double total = 0;
    for (size_t i = 0; i < filenames.size(); ++i) {
      filesizes.push_back(get_filesize(filenames[i]));
      total += filesizes.back();
    }
    vector<size_t> samples;
    for (size_t i = 0; i < filesizes.size(); ++i)
      samples.push_back(n_reads*filesizes[i]/total);
    
    for (size_t i = 0; i < filenames.size(); ++i) {
      if (VERBOSE)
	cerr << filenames[i] << endl;
      
      vector<string> names, sequences;
      read_fasta_file(filenames[i].c_str(), names, sequences);
      for (size_t j = 0; j < names.size(); ++j) {
	const size_t offset = names[j].find(':');
	const string name(names[j].substr(0, offset));
	simreadspe(rng, samples[i], read_width, separation,
		   max_errors, name, sequences[j], read_names, reads);
      }
    }
    
    std::ofstream out(outfile.c_str());
    for (size_t i = 0; i < reads.size(); ++i)
      out << ">" << read_names[i] << endl << reads[i] << endl;
    out.close();
  }      
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  catch (RMAPException &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::exception &e) {
    cerr << "ERROR: " << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
