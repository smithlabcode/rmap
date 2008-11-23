/*
 *    Part of RMAP software
 *
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

#include "sim_utils.hpp"
#include "rmap_utils.hpp"

#include <numeric>
#include <algorithm>
#include <set>
#include <cmath>

using std::string;
using std::vector;
using std::accumulate;
using std::min;
using std::log;

using std::cerr;
using std::endl;

void
add_sequencing_errors(const Runif &rng, const double max_errors,
		      string &seq, string &error_log) {
  size_t no_err,temp;
  temp = max_errors+1;
  error_log = string(seq.length(), '0');
  std::set<size_t> errors;
  no_err = double(rng.runif(1ul, temp) + 0.7);
  for (size_t i = 0; i < no_err; ++i) {
    size_t error_pos = rng.runif(0ul, seq.length());
    while (errors.find(error_pos) != errors.end())
	error_pos = rng.runif(0ul, seq.length());
    errors.insert(error_pos);
    error_log[error_pos] = '1';

    size_t c = rng.runif(0ul, rmap::alphabet_size);
    while (c == base2int(seq[error_pos]))
      c = rng.runif(0ul, rmap::alphabet_size);
    
    seq[error_pos] = int2base(c);
  }

}


void
generate_sequencing_errors(const Runif &rng, 
			   const size_t read_width,const double total_error, 
			   vector<vector<double> > &errors) {
  
  errors.resize(read_width, vector<double>(rmap::alphabet_size, 0));
  double remaining_error = total_error;
  while (remaining_error > 0) {
    const size_t error_pos = rng.runif(0ul, read_width);
    size_t error_base = rng.runif(0ul, rmap::alphabet_size);
    
    const double error_amount = min(min(1.0 - errors[error_pos][error_base], remaining_error),
				    rng.runif(0.0, 1.0));
    
    remaining_error -= error_amount;
    errors[error_pos][error_base] += error_amount;
  }
}

void
add_sequencing_errors(const vector<vector<double> > &errors, vector<vector<double> > &prb) {
  for (size_t i = 0; i < prb.size(); ++i) {
    size_t base = max_element(prb[i].begin(), prb[i].end()) - prb[i].begin();
    const double sum = accumulate(errors[i].begin(), errors[i].end(), 0.0);
    prb[i][base] -= sum;
    prb[i][base] = std::max(0.0, prb[i][base]);
    transform(prb[i].begin(), prb[i].end(), 
	      errors[i].begin(), prb[i].begin(),
	      std::plus<double>());
  }
}

void
adjust_seq_using_matrix(const vector<vector<double> > &prb, string &seq) {
  assert(prb.size() == seq.length());
  for (size_t i = 0; i < prb.size(); ++i) {
    size_t base = max_element(prb[i].begin(), prb[i].end()) - prb[i].begin();
    seq[i] = int2base(base);
  }
}

void
prob_to_quality_scores_solexa(const vector<vector<double> > &prb, vector<vector<double> > &quality) {
  quality = prb;
  for (size_t i = 0; i < prb.size(); ++i) {
    std::transform(prb[i].begin(), prb[i].end(),
		   quality[i].begin(), std::bind2nd(std::plus<double>(), 1e-3));
//     cerr << quality[i][0] << "\t" << prb[i][0] << endl;
//     assert(quality[i][0] > 0);
    const double column_sum = accumulate(quality[i].begin(), quality[i].end(), 0.0);
//     cerr << column_sum << endl;
//     assert(quality[i][0] > 0);
//     cerr << quality[i][0] << endl;
    std::transform(quality[i].begin(), quality[i].end(),
		   quality[i].begin(), std::bind2nd(std::divides<double>(), column_sum));
//     cerr << quality[i][0] << endl;
//     assert(quality[i][0] > 0);
  }
  for (size_t i = 0; i < quality.size(); ++i)
    for (size_t j = 0; j < rmap::alphabet_size; ++j) {
      assert(quality[i][j] > 0);
      quality[i][j] = 10*(log(quality[i][j]) - log(1 - quality[i][j]))/log(10);
    }
}

void
add_sequencing_errors(const Runif &rng, const double max_errors,
		      string &seq, vector<vector<double> > &quality_scores) {
  
  // first make the pwm:
  quality_scores.resize(seq.length(), vector<double>(rmap::alphabet_size, 0.0));
  for (size_t i = 0; i < seq.length(); ++i)
    quality_scores[i][base2int(seq[i])] = 1.0;
  
  double total_error = max_errors;
  while (total_error > 0) {
    // sample an error position:
    const size_t error_pos = rng.runif(0ul, seq.length());
    // sample an error amount:
    double remaining_freq = min(quality_scores[error_pos][base2int(seq[error_pos])], total_error);
    const double error_amount = min(min(rng.runif(0.0, 1.0), max_errors), 
				    remaining_freq);
    size_t error_base = base2int(seq[error_pos]);
    while (error_base == base2int(seq[error_pos]))
      error_base = rng.runif(0ul, rmap::alphabet_size);
    
    quality_scores[error_pos][base2int(seq[error_pos])] -= error_amount;
    quality_scores[error_pos][error_base] += error_amount;
    
    total_error -= error_amount;
  }

  for (size_t i = 0; i < quality_scores.size(); ++i) {
    std::transform(quality_scores[i].begin(), quality_scores[i].end(),
		   quality_scores[i].begin(), std::bind2nd(std::plus<double>(), 1e-3));
    const double column_sum = accumulate(quality_scores[i].begin(),
					 quality_scores[i].end(), 0.0);
    std::transform(quality_scores[i].begin(), quality_scores[i].end(),
		   quality_scores[i].begin(), std::bind2nd(std::divides<double>(), column_sum));
  }
  
  for (size_t i = 0; i < quality_scores.size(); ++i)
    for (size_t j = 0; j < rmap::alphabet_size; ++j)
      quality_scores[i][j] = 10*(log(quality_scores[i][j]) - log(1 - quality_scores[i][j]))/log(10);
}

void
complement_score_matrix(const vector<vector<double> > &matrix,
			const double max_quality_score,
			vector<vector<double> > &scores) {
  scores = matrix;
  for (size_t i = 0; i < scores.size(); ++i) {
    //     for (size_t j = 0; j < scores[i].size(); ++j)
    //       assert(finite(scores[i][j]));
//    const size_t max_idx = max_element(scores[i].begin(), scores[i].end())
      //    - scores[i].begin();
    //const double max_score = scores[i][max_idx];
    for (size_t j = 0; j < scores[i].size(); ++j) {
      scores[i][j] = max_quality_score - scores[i][j];
    }
  }
}
