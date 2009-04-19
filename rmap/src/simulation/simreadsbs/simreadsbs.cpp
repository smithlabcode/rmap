/*    simreadsbs: a program for simulating Solexa reads to test rmapbs
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

#include "rmap_utils.hpp"
#include "rmap_os.hpp"
#include "sim_utils.hpp"
#include "bisulfite_utils.hpp"
#include "RNG.hpp"
#include "OptionParser.hpp"

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <set>

using std::ofstream;
using std::ostream_iterator;
using std::numeric_limits;
using std::string;
using std::vector;
using std::ostream;
using std::endl;
using std::cerr;
using std::ptr_fun;


void
simreadsbs(const Runif rng,
	   const size_t n_reads, const size_t read_width, 
	   const size_t max_errors, 
	   const double bs_rate, const double meth_rate,
	   const string &name, const string &sequence,
	   vector<string> &read_names, vector<string> &reads,
	   vector<string> &reads_bs) {

  const size_t lim = sequence.length() - read_width + 1;

  std::set<size_t> used;

  for (size_t i = 0; i < n_reads; ++i) {
    bool found_one = false;
    while (!found_one) {
      const size_t start = rng.runif(0ul, lim);
      
      bool valid = true;
      for (size_t j = start; j < start + read_width && valid; ++j)
	valid = isvalid(sequence[j]);
      
      if (valid) {
	string seq(sequence.substr(start, read_width));
	const bool rc = (rng.runif(0.0,1.0) > 0.5);
	if (rc) seq = revcomp(seq);
	
	transform(seq.begin(), seq.end(), seq.begin(), std::ptr_fun(&toupper));
	string error_log;
	add_sequencing_errors(rng, max_errors, seq, error_log);
	
	reads.push_back(seq);
	bisulfite_treatment(rng, seq, bs_rate, meth_rate);
	const string read_name(name + ":" + toa(start) + "-" + 
			       toa(start + read_width) + "_" + toa(!rc) + "_" +
			       error_log + "_" + toa(count(error_log.begin(), 
							   error_log.end(), '1')));
	reads_bs.push_back(seq);
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
    string unconverted_outfile;
    size_t n_reads = 1000;
    size_t read_width = 25;
    size_t max_errors = 0;
    double meth_rate = 0.0;
    double bs_rate = 1.0;
    size_t random_number_seed = -numeric_limits<size_t>::max();
    
    bool VERBOSE = false;
    
    /****************** COMMAND LINE OPTIONS ********************/
    static OptionParser opt_parse("simreadsbs",
				  "program for generating simulated reads",
				  "<fasta-chrom-files>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", false , outfile);
    opt_parse.add_opt("reads", 'n', "number of reads to simulate", false , n_reads);
    opt_parse.add_opt("width", 'w', "width of reads to simulate", false, read_width);
    opt_parse.add_opt("err", 'e', "maximum number of simulated sequencing errors", false, max_errors);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    opt_parse.add_opt("unconv", 'u', "output file for unconverted reads", false, unconverted_outfile);
    opt_parse.add_opt("meth", 'm', "rate of CpG methylation", false, meth_rate); 
    opt_parse.add_opt("bs", 'b', "rate of bisulfite conversion", false, bs_rate);
    opt_parse.add_opt("seed", 'S', "random number seed", false, random_number_seed);
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
    /**********************************************************************/

    const Runif rng(random_number_seed);

    vector<string> reads, reads_bs, read_names;

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
	simreadsbs(rng, samples[i], read_width, max_errors, bs_rate, meth_rate, 
		   name, sequences[j], read_names, reads, reads_bs);
      }
    }
    
    if (!unconverted_outfile.empty()) {
      ofstream out(unconverted_outfile.c_str());
      for (size_t i = 0; i < reads.size(); ++i)
	out << ">" << read_names[i] << endl << reads[i] << endl;
      out.close();
    }
    
    ofstream out(outfile.c_str());
    for (size_t i = 0; i < reads_bs.size(); ++i)
      out << ">" << read_names[i] << endl << reads_bs[i] << endl;
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
