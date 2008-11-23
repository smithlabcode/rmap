/*    simreadsbsq: a program for simulating Solexa reads to test rmapbsq
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
#include <cmath>

using std::ofstream;
using std::ostream_iterator;
using std::numeric_limits;
using std::string;
using std::vector;
using std::ostream;
using std::endl;
using std::cerr;
using std::ptr_fun;
using std::max;
using std::min;

bool
finite_matrix(const vector<vector<double> > &matrix) {
  for (size_t i = 0; i < matrix.size(); ++i) {
    for (size_t j = 0; j < matrix[i].size(); ++j)
      if (!finite(matrix[i][j])) return false;
  }
  return true;
}

void
get_max_and_min_quality_scores(const vector<vector<double> > &scores, 
			       double &min_quality_score, 
			       double &max_quality_score) {
  for (size_t i = 0; i < scores.size(); ++i) {
    max_quality_score = max(max_quality_score, 
			    *max_element(scores[i].begin(), 
					 scores[i].end()));
    min_quality_score = min(min_quality_score, 
			    *min_element(scores[i].begin(), 
					 scores[i].end()));
  }
}


void
seq_to_quality_scores(const string &seq, vector<vector<double> > &quality_scores) {
  quality_scores.resize(seq.length(), vector<double>(rmap::alphabet_size, 0));
  for (size_t i = 0; i < seq.length(); ++i)
    quality_scores[i][base2int(seq[i])] = 1;
}

double
score_against_sequence_bs(const vector<vector<double> > &quality_scores, const string &seq) {
  double score = 0;
  for (size_t i = 0; i < seq.length(); ++i) {
    if (toupper(seq[i]) == 'C')
      score += std::min(quality_scores[i][1], quality_scores[i][3]);
    else 
      score += quality_scores[i][base2int(seq[i])];
  }
  return score;
}

void
simreads_bs(const Runif &rng, 
	    const size_t n_reads, 
	    const size_t read_width, 
	    const size_t max_errors, 
	    const double meth_rate,
	    
	    const string &name, 
	    const string &sequence,
	    
	    vector<string> &read_names, 
	    vector<string> &reads,
	    vector<string> &reads_bs,
	    vector<vector<vector<double> > > &probs) {
  
  const size_t lim = sequence.length() - read_width + 1;

  std::set<size_t> used;

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
    transform(seq.begin(), seq.end(), seq.begin(), std::ptr_fun(&toupper));
    const bool rc = (rng.runif(0.0,1.0) > 0.5);
    if (rc) seq = revcomp(seq);
    
    // Do the bisulfite treatment
    string seq_bs(seq);
    bisulfite_treatment(rng, seq_bs, meth_rate);
    
    // Sample errors
    vector<vector<double> > errors;
    generate_sequencing_errors(rng, read_width, max_errors, errors);
    
    vector<vector<double> > quality_scores;
    seq_to_quality_scores(seq_bs, quality_scores);
    
    add_sequencing_errors(errors, quality_scores);
    
    prob_to_quality_scores_solexa(quality_scores, quality_scores);
    
    double max_quality_score = 0;
    double min_quality_score = numeric_limits<double>::max();
    get_max_and_min_quality_scores(quality_scores, min_quality_score, max_quality_score);

    assert(finite_matrix(quality_scores));

    adjust_seq_using_matrix(quality_scores, seq_bs);
    
    vector<vector<double> > score_matrix;
    complement_score_matrix(quality_scores, max_quality_score, score_matrix);
    
    // Make name
    const string read_name(name + ":" + toa(start) + "-" + 
			   toa(start + read_width) + "_" + toa(!rc) + "_" +
			   toa(score_against_sequence_bs(score_matrix, seq)/
			       (max_quality_score - min_quality_score)));
    
    // Push back what was sampled
    reads.push_back(seq);
    reads_bs.push_back(seq_bs);
    read_names.push_back(read_name);
    probs.push_back(quality_scores);
  }
}

int
main(int argc, const char **argv) {

  try {
    
    string prb_file;
    string seqfile;
    string outfile;
    string unconverted_outfile;
    size_t n_reads = 1000;
    size_t read_width = 25;
    size_t max_errors = 0;
    double meth_rate = 0;
    size_t random_number_seed = -numeric_limits<size_t>::max();

    bool VERBOSE = false;
    
    /****************** COMMAND LINE OPTIONS ********************/
    static OptionParser 
      opt_parse("simreadsbsq",
		"program for generating simulated "
		"bisulfite treated reads with "
		"simulated quality scores");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", false , outfile);
    opt_parse.add_opt("reads", 'n', "number of reads to simulate", false , n_reads);
    opt_parse.add_opt("width", 'w', "width of reads to simulate", false, read_width);
    opt_parse.add_opt("err", 'e', "maximum number of simulated sequencing errors", false, max_errors);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    opt_parse.add_opt("prob", 'p', "prb output file", true, prb_file);
    opt_parse.add_opt("unconv", 'u', "output file for unconverted reads", false, unconverted_outfile);
    opt_parse.add_opt("meth", 'm', "rate of CpG methylation", false, meth_rate);
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
    /****************** END COMMAND LINE OPTIONS *****************/

    const Runif rng(random_number_seed);

    vector<string> reads, reads_bs, read_names;
    vector<vector<vector<double> > > probs;

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
	simreads_bs(rng, samples[i], read_width, max_errors, meth_rate, 
		    name, sequences[j], read_names, reads, reads_bs, probs);
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

    out.open(prb_file.c_str());
    for (size_t i = 0; i < probs.size(); ++i) {
      for (size_t j = 0; j < probs[i].size(); ++j) {
	copy(probs[i][j].begin(), probs[i][j].end(),
	     ostream_iterator<double>(out, "\t"));
      }
      out << endl;
    }
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
