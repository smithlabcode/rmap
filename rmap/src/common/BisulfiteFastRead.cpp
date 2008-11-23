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

#include "BisulfiteFastRead.hpp"

#include <cmath>

using std::string;
using std::cerr;
using std::endl;

size_t BisulfiteFastRead::score_mask = 0;
size_t BisulfiteFastRead::segments = 0;
size_t BisulfiteFastRead::read_width = 0;
size_t BisulfiteFastRead::right_most_bit = 0;

////////////////////////////////////////////////////////////////////////
// WORD PAIR

BisulfiteFastRead::WordPair::WordPair(const string &s) : 
  a_vec(0), c_vec(0), g_vec(0), t_vec(0), bads(rmap_bits::all_ones) {
  const string::const_iterator limit = s.end();
  for (string::const_iterator i(s.begin()); i != limit; ++i) {
    a_vec = ((a_vec << 1) + (contains_a(*i)));
    c_vec = ((c_vec << 1) + (contains_c_bs(*i)));
    g_vec = ((g_vec << 1) + (contains_g(*i)));
    t_vec = ((t_vec << 1) + (contains_t(*i)));
    bads  = ((bads  << 1) + (contains_n(*i)));
  }
  if (s.length() < rmap_bits::word_size) {
    const size_t additional_shift = (rmap_bits::word_size - s.length());
    a_vec <<= additional_shift;
    c_vec <<= additional_shift;
    g_vec <<= additional_shift;
    t_vec <<= additional_shift;
    bads  <<= additional_shift;
    bads += ((1 << additional_shift) - 1);
  }
}


char
BisulfiteFastRead::WordPair::get_char(size_t mask, size_t pos) const {
  // 00 -> A, 01 -> C, 10 -> G, 11 -> T
  const MASK_t selector = (rmap_bits::low_bit << (pos - 1));
  
  if ((mask & bads) & selector)
    return 'N';
  
  if ((mask & a_vec) & selector) {
    return 'A';
  }
  if (((mask & c_vec) & selector) && !((mask & t_vec) & selector)) {
    return 'C';
  }
  if ((mask & g_vec) & selector) {
    return 'G';
  }
  if ((mask & t_vec) & selector) {
    return 'T';
  }
  return 'N';
}


string
BisulfiteFastRead::WordPair::bits2string(size_t mask, size_t bits) {
  string s;
  size_t selector = rmap_bits::high_bit;
  for (size_t i = 0; i < rmap_bits::word_size; ++i) {
    s += (selector & bits & mask) ? '1' : '0';
    selector >>= 1;
  }
  return s;
}


string
BisulfiteFastRead::WordPair::tostring_bits(size_t mask) const {
  const string s(bits2string(mask, a_vec) + "\n" +
		 bits2string(mask, c_vec) + "\n" +
		 bits2string(mask, g_vec) + "\n" +
		 bits2string(mask, t_vec) + "\n" +
		 bits2string(rmap_bits::all_ones, bads) + "\n");
  string seq;
  for (size_t i = rmap_bits::word_size; i > 0; --i)
    seq += get_char(mask, i);
  return s + seq;
}

string
BisulfiteFastRead::WordPair::tostring_bases(size_t mask) const {
  string seq;
  for (size_t i = rmap_bits::word_size; i > 0; --i)
    seq += get_char(mask, i);
  return seq;
}

void
BisulfiteFastRead::WordPair::shift(const BisulfiteFastRead::WordPair &other) {
  a_vec <<= 1;
  c_vec <<= 1;
  g_vec <<= 1;
  t_vec <<= 1;
  bads  <<= 1;

  a_vec += ((other.a_vec & rmap_bits::high_bit) != 0);
  c_vec += ((other.c_vec & rmap_bits::high_bit) != 0);
  g_vec += ((other.g_vec & rmap_bits::high_bit) != 0);
  t_vec += ((other.t_vec & rmap_bits::high_bit) != 0);
  bads  += ((other.bads  & rmap_bits::high_bit) != 0);
}

////////////////////////////////////////////////////////////////////////
// FAST READ

void
BisulfiteFastRead::set_read_width(const size_t m) {
  read_width = m;
  score_mask = (rmap_bits::low_bit << (m % rmap_bits::word_size)) - 1;
  right_most_bit = (rmap_bits::word_size - (m % rmap_bits::word_size));
  score_mask <<= right_most_bit;
  segments = size_t(std::ceil(m/static_cast<float>(rmap_bits::word_size)));
}

BisulfiteFastRead::BisulfiteFastRead(const std::string &s) {
  for (size_t i = 0; i < segments - 1; ++i) {
    const string this_seg(s.substr(i*rmap_bits::word_size, rmap_bits::word_size));
    wp.push_back(WordPair(this_seg));
  }
  wp.push_back(WordPair(s.substr((segments - 1)*rmap_bits::word_size)));
}

string
BisulfiteFastRead::tostring_bases() const {
  std::ostringstream ss;
  for (size_t i = 0; i < wp.size() - 1; ++i)
    ss << wp[i].tostring_bases(rmap_bits::all_ones);
  ss << wp.back().tostring_bases(score_mask);
  return ss.str();
}

string
BisulfiteFastRead::tostring_bits() const {
  std::ostringstream ss;
  for (size_t i = 0; i < wp.size() - 1; ++i)
    ss << wp[i].tostring_bits(rmap_bits::all_ones) << endl;
  ss << wp.back().tostring_bits(score_mask) << endl;
  return ss.str();
}
