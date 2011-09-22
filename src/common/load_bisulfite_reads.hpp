/*
 *    Part of RMAP software
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

#ifndef LOAD_BISULFITE_READS_HPP
#define LOAD_BISULFITE_READS_HPP

#include "rmap_utils.hpp"
#include "BisulfiteFastRead.hpp"
#include "BisulfiteFastReadQuality.hpp"
#include "BisulfiteFastReadWC.hpp"

#include <vector>
#include <string>


void
load_reads_from_fasta_file(const std::string &filename, const size_t max_diffs,
			   size_t &read_width, std::vector<BisulfiteFastRead> &fast_reads,
			   std::vector<size_t> &read_words, std::vector<size_t> &read_index);

void
load_reads_from_fastq_file(const std::string &filename, const size_t max_diffs,
			   size_t &read_width, std::vector<BisulfiteFastRead> &fast_reads,
			   std::vector<size_t> &read_words, std::vector<size_t> &read_index);

void
load_reads_from_fastq_file(const std::string &filename, const size_t max_diffs,
			   size_t &read_width, std::vector<BisulfiteFastReadWC> &fast_reads,
			   std::vector<size_t> &read_words, std::vector<size_t> &read_index);

void
load_reads_from_fastq_file(const std::string &filename, const size_t max_diffs,
			   size_t &read_width, std::vector<BisulfiteFastReadQuality> &fast_reads,
			   std::vector<size_t> &read_words, std::vector<size_t> &read_index);

void
load_reads_from_prb_file(const std::string &filename, const size_t max_diffs,
			 size_t &read_width, std::vector<BisulfiteFastReadWC> &fast_reads,
			 std::vector<size_t> &read_words, std::vector<size_t> &read_index);

void
load_reads_from_prb_file(const std::string &filename, const size_t max_diffs,
			 size_t &read_width, std::vector<BisulfiteFastReadQuality> &fast_reads,
			 std::vector<size_t> &read_words, std::vector<size_t> &read_index);


#endif
