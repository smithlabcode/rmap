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

#include "BisulfiteFastReadWC.hpp"

#include <cmath>
#include <iomanip>

using std::string;
using std::vector;
using std::cerr;
using std::endl;

size_t BisulfiteFastReadWC::score_mask = 0;
size_t BisulfiteFastReadWC::segments = 0;
size_t BisulfiteFastReadWC::read_width = 0;
size_t BisulfiteFastReadWC::right_most_bit = 0;

////////////////////////////////////////////////////////////////////////
// WORD PAIR

size_t
BisulfiteFastReadWC::Words::quality_to_value(double quality) {
  return quality > 0.995;
}

BisulfiteFastReadWC::Words::Words(const vector<vector<double> > &s,
				  bool AG_WILDCARD) : 
  a_vec(0), c_vec(0), g_vec(0), t_vec(0) {
  const vector<vector<double> >::const_iterator limit = s.end();
  for (vector<vector<double> >::const_iterator i(s.begin()); i != limit; ++i) {
    a_vec = ((a_vec << 1) + (quality_to_value((*i)[0])));
    if (AG_WILDCARD) {
      c_vec = ((c_vec << 1) + (quality_to_value((*i)[1])));
      g_vec = ((g_vec << 1) + (std::min(quality_to_value((*i)[2]), quality_to_value((*i)[0]))));
    }
    else {
      c_vec = ((c_vec << 1) + (std::min(quality_to_value((*i)[1]), quality_to_value((*i)[3]))));
      g_vec = ((g_vec << 1) + (quality_to_value((*i)[2])));
    }
    t_vec = ((t_vec << 1) + (quality_to_value((*i)[3])));
  }
  if (s.size() < rmap_bits::word_size) {
    const size_t additional_shift = (rmap_bits::word_size - s.size());
    a_vec <<= additional_shift;
    c_vec <<= additional_shift;
    g_vec <<= additional_shift;
    t_vec <<= additional_shift;
  }
}

BisulfiteFastReadWC::Words::Words(string::const_iterator i, const string::const_iterator limit) :
  a_vec(0), c_vec(0), g_vec(0), t_vec(0) {
  const size_t length = limit - i;
  for (; i < limit; ++i) {
    a_vec = ((a_vec << 1) | contains_a(base2int(*i)));
    c_vec = ((c_vec << 1) | contains_c(base2int(*i)));
    g_vec = ((g_vec << 1) | contains_g(base2int(*i)));
    t_vec = ((t_vec << 1) | contains_t(base2int(*i)));
  }
  if (length < rmap_bits::word_size) {
    const size_t additional_shift = (rmap_bits::word_size - length);
    a_vec <<= additional_shift;
    c_vec <<= additional_shift;
    g_vec <<= additional_shift;
    t_vec <<= additional_shift;
  }
}


size_t
BisulfiteFastReadWC::Words::get_val(MASK_t mask, MASK_t base_vec, size_t pos) {
  // 00 -> A, 01 -> C, 10 -> G, 11 -> T
  const MASK_t selector = (BisulfiteFastReadWC::low_val_bits << (pos - 1));
  return (((mask & base_vec) & selector) >> (pos - 1));
}


string
BisulfiteFastReadWC::Words::bits2string(size_t mask, size_t bits) {
  string s;
  size_t selector = rmap_bits::high_bit;
  for (size_t i = 0; i < rmap_bits::word_size; ++i) {
    s += (selector & bits & mask) ? '1' : '0';
    selector >>= 1;
  }
  return s;
}


string
BisulfiteFastReadWC::Words::tostring_bits(size_t mask) const {
  return (bits2string(mask, a_vec) + "\n" +
	  bits2string(mask, c_vec) + "\n" +
	  bits2string(mask, g_vec) + "\n" +
	  bits2string(mask, t_vec) + "\n");

}


////////////////////////////////////////////////////////////////////////
// FAST READ

void
BisulfiteFastReadWC::set_read_width(const size_t rw) {
  read_width = rw;
  segments = static_cast<size_t>(std::ceil(rw/static_cast<double>(segment_size))) - 1;
  right_most_bit = (rmap_bits::word_size - (rw % rmap_bits::word_size));
  score_mask = (rmap_bits::all_ones << right_most_bit);
}


BisulfiteFastReadWC::BisulfiteFastReadWC(const vector<vector<double> > &s,
					 bool AG_WILDCARD) {
  words.resize(segments + 1);
  for (size_t i = 0; i < segments; ++i) {
    const vector<vector<double> > 
      this_seg(s.begin() + i*segment_size, s.begin() + (i + 1)*segment_size);
    words[i] = Words(this_seg, AG_WILDCARD);
  }
  const vector<vector<double> > this_seg(s.begin() + segments*segment_size, s.end());
  words.back() = Words(this_seg, AG_WILDCARD);
}


string
BisulfiteFastReadWC::tostring_bits() const {
  std::ostringstream ss;
  for (size_t i = 0; i < segments; ++i)
    ss << words[i].tostring_bits(rmap_bits::all_ones) << std::endl;
  ss << words.back().tostring_bits(score_mask);
  return ss.str();
}
