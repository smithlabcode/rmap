/*    simreads: a program for simulating Solexa reads to test smithlab
 *
 *    Copyright (C) 2008 University of Southern California and
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

#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "RNG.hpp"
#include "sim_utils.hpp"
#include "bisulfite_utils.hpp"
#include "OptionParser.hpp"
#include "QualityScore.hpp"

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

using std::ofstream;
using std::ostream_iterator;
using std::numeric_limits;
using std::string;
using std::vector;
using std::ostream;
using std::endl;
using std::cerr;
using std::max;
using std::min;



static void
prb_to_fastq(const string &read,
	     const vector<vector<double> > &prb,
	     string &fastq) {
  for (size_t i = 0; i < read.length(); ++i) {
    const double score = *max_element(prb[i].begin(), prb[i].end());
    fastq += solexa_to_quality_character(score + 5);
  }
}


static void
write_read_fastq(ostream &out, 
		 const string &read,
		 const string &read_name,
		 const vector<vector<double> > &probs) {
  string fastq;
  prb_to_fastq(read, probs, fastq);
  out << '@' << read_name << '\n'
      << read << "\n+" << read_name << '\n'
      << fastq << '\n';
}
  

static void
write_read_fasta(ostream &out, const string &read, const string &read_name) {
  out << '>' << read_name << '\n' << read << '\n';
}

static void
write_read_prb(ofstream &out, 
	       const vector<vector<double> > &probs) {
  for (size_t j = 0; j < probs.size(); ++j) {
    if (j != 0) out << '\t';
    out << probs[j].front();
    for (size_t k = 1; k < probs[j].size(); ++k)
      out << '\t' << probs[j][k];
  }
  out << '\n';
}

static void
get_error_log(const string &seq, const string &called_seq, 
	      string &error_log) {
  error_log = string(seq.length(), '0');
  for (size_t i = 0; i < seq.length(); ++i)
    if (seq[i] != called_seq[i])
      error_log[i] = '1';
}



static void
simreads(const bool FASTQ_OUTPUT,
	 const string &outfile,
	 const string &prb_file,
	 const Runif &rng,
	 const size_t n_reads, const size_t read_width, const size_t max_errors, 
	 const string &name, const string &sequence,
	 vector<string> &read_names, vector<string> &reads,
	 vector<vector<vector<double> > > &probs) {

  ostream *out = (outfile.empty()) ? &std::cout : 
    new ofstream(outfile.c_str(), std::ios::app);
  ofstream *prb = (prb_file.empty()) ? 0 : 
    new ofstream(prb_file.c_str(), std::ios::app);
  
  const size_t lim = sequence.length() - read_width + 1;
  for (size_t i = 0; i < n_reads; ++i) {
    
    // Get a valid starting point
    size_t start = numeric_limits<size_t>::max();
    while (start == numeric_limits<size_t>::max()) {
      start = rng.runif(0ul, lim);
      for (size_t j = 0; j < read_width && 
	     start != numeric_limits<size_t>::max(); ++j)
	if (!isvalid(sequence[start + j]))
	  start = numeric_limits<size_t>::max();
    }
    
    // extract the sequence (and decide one revcomp)
    string seq(sequence.substr(start, read_width));
    const bool rc = (rng.runif(0.0,1.0) > 0.5);
    if (rc) seq = revcomp(seq);

    vector<vector<double> > matrix;
    sequence_to_consensus_matrix(seq, matrix);
    add_sequencing_errors(rng, max_errors, matrix);
    
    string called_seq;
    call_bases_solexa(matrix, called_seq);
    
    // Sample errors
    string error_log;
    get_error_log(seq, called_seq, error_log);
    
    for (size_t j = 0; j < matrix.size(); ++j)
      for (size_t k = 0; k < matrix[j].size(); ++k)
	matrix[j][k] = round(error_probability_to_solexa(1.0 - matrix[j][k]));
    
    // Make name
    size_t actual_mismatches = count(error_log.begin(), error_log.end(), '1');
    const string read_name(name + ":" + toa(start) + "-" + 
			   toa(start + read_width) + "_" + toa(!rc) + "_" +
			   error_log + "_" + toa(actual_mismatches));
    
    assert(matrix.size() == read_width);
    if (FASTQ_OUTPUT)
      write_read_fastq(*out, called_seq, read_name, matrix);
    else {
      if (prb != 0)
	write_read_prb(*prb, matrix);
      write_read_fasta(*out, called_seq, read_name);
    }
  }
  if (out != &std::cout) delete out;
  if (prb) delete prb;
}



int
main(int argc, const char **argv) {

  try {
    
    string prb_file;
    string seqfile;
    string outfile;
    size_t n_reads = 1000;
    size_t read_width = 25;
    size_t max_errors = 0;
    size_t random_number_seed = numeric_limits<size_t>::max();
    
    bool VERBOSE = false;
    bool FASTQ_OUTPUT = false;

    /****************** COMMAND LINE OPTIONS ********************/
    static OptionParser opt_parse("simreads",
				  "program for generating simulated reads"
				  " with simulated quality scores",
				  "<fasta-chrom-files>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)",
		      false , outfile);
    opt_parse.add_opt("reads", 'n', "number of reads to simulate", 
		      false, n_reads);
    opt_parse.add_opt("width", 'w', "width of reads to simulate", 
		      false, read_width);
    opt_parse.add_opt("err", 'e', "maximum number of simulated sequencing errors",
		      false, max_errors);
    opt_parse.add_opt("verbose", 'v', "print more run info", 
		      false, VERBOSE);
    opt_parse.add_opt("fastq", 'q', "write FASTQ format reads", 
		      false, FASTQ_OUTPUT);
    opt_parse.add_opt("prob", 'p', "prb output file", 
		      false, prb_file);
    opt_parse.add_opt("seed", 'S', "random number seed", 
		      false, random_number_seed);
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
    
    if (FASTQ_OUTPUT && !prb_file.empty())
      throw SMITHLABException("fastq output is incompatible "
			  "with specifying a prb file");

    const Runif rng(random_number_seed);
    
    vector<string> reads, read_names;
    vector<vector<vector<double> > > probs;
    
    vector<size_t> filesizes;
    double total = 0;
    for (size_t i = 0; i < filenames.size(); ++i) {
      filesizes.push_back(get_filesize(filenames[i]));
      total += filesizes.back();
    }

    if (VERBOSE)
      cerr << "[OBTAINING READ SAMPLE DISTRIBUTION]";
    vector<size_t> samples;
    for (size_t i = 0; i < filesizes.size(); ++i)
      samples.push_back(n_reads*filesizes[i]/total);
    if (VERBOSE) cerr << "[DONE]" << endl;
    
    if (!outfile.empty())
      ofstream out(outfile.c_str());
    
    if (!prb_file.empty())
      ofstream prb(prb_file.c_str());
    
    for (size_t i = 0; i < filenames.size(); ++i) {
      if (isdir(filenames[i].c_str()))
	throw SMITHLABException("\"" + filenames[i] + 
			    "\" not a FASTA format sequence file?");
      if (VERBOSE)
	cerr << "[LOADING=" << filenames[i] << "]";
      vector<string> names, sequences;
      read_fasta_file(filenames[i].c_str(), names, sequences);
      if (VERBOSE) cerr << "[DONE]" << endl;

      double sub_total = 0;
      for (size_t k = 0; k < sequences.size(); ++k)
	sub_total += sequences[i].length();
      
      vector<size_t> sub_samples;
      for (size_t k = 0; k < sequences.size(); ++k)
	sub_samples.push_back(samples[i]*sequences[i].length()/sub_total);
      
      for (size_t j = 0; j < names.size(); ++j) {
	if (VERBOSE)
	  cerr << "[PROCESSING CHROM=" << names[j] << "]" << endl;
	const size_t offset = names[j].find(':');
	const string name(names[j].substr(0, offset));
	
	transform(sequences[j].begin(), sequences[j].end(), 
		  sequences[j].begin(), std::ptr_fun(&toupper));
	
	simreads(FASTQ_OUTPUT, outfile, prb_file,
		 rng, sub_samples[i], read_width, max_errors, 
		 name, sequences[j], read_names, reads, probs);
      }
    }
  }      
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  catch (SMITHLABException &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::exception &e) {
    cerr << "ERROR: " << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
