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

#ifndef SIM_UTILS_HPP
#define SIM_UTILS_HPP

#include <string>
#include <vector>
#include "RNG.hpp"

void
add_sequencing_errors(const Runif &rng, const double n_errors,
		      std::string &seq, std::string &error_log);

void
add_sequencing_errors(const Runif &rng, const double n_errors,
		      std::string &seq, std::vector<std::vector<double> > &quality_scores);

void
generate_sequencing_errors(const Runif &rng, 
			   const size_t read_width,const double total_error, 
			   std::vector<std::vector<double> > &errors);

void
add_sequencing_errors(const std::vector<std::vector<double> > &errors,
		      std::vector<std::vector<double> > &prb);


void
adjust_seq_using_matrix(const std::vector<std::vector<double> > &prb, std::string &seq);


void
prob_to_quality_scores_solexa(const std::vector<std::vector<double> > &prb, 
			      std::vector<std::vector<double> > &quality);


void
complement_score_matrix(const std::vector<std::vector<double> > &matrix,
			const double max_quality_score,
			std::vector<std::vector<double> > &scores);

#endif


