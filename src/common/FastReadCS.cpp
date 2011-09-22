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

#include "FastReadCS.hpp"

#include <cmath>

using std::string;

////////////////////////////////////////////////////////////////////////
// WORD PAIR

using std::cerr;
using std::endl;

FastReadCS::WordPair::WordPair(const string &s) : 
  upper(0), lower(0), bads(static_cast<size_t>(-1)) {
  string::const_iterator i = s.begin();
  const string::const_iterator limit = s.end();
  while (i != limit) {
    const char c = base2int(*i) & static_cast<size_t>(3);
    upper = ((upper << 1) + get_upper(c));
    lower = ((lower << 1) + get_lower(c));
    bads  = ((bads << 1) + get_bads(*i));
    ++i;
  }
  if (s.length() < rmap_bits::word_size) {
    const size_t additional_shift = (rmap_bits::word_size - s.length());
    upper <<= additional_shift;
    lower <<= additional_shift;
    bads <<= additional_shift;
    bads += ((1ul << additional_shift) - 1);
  }
}


void
FastReadCS::WordPair::to_color_space(const WordPair &wp) {
  upper = (upper ^ ((upper >> 1) | (wp.upper << 63)));
  lower = (lower ^ ((lower >> 1) | (wp.lower << 63)));
}

void
FastReadCS::WordPair::to_color_space(const size_t prev) {
  upper = (upper ^ ((upper >> 1) | (get_upper(prev) << 63)));
  lower = (lower ^ ((lower >> 1) | (get_lower(prev) << 63)));
}


void
FastReadCS::to_color_space(const size_t prev) {
  for (size_t i = segments; i > 0; --i)
    wp[i].to_color_space(wp[i - 1]);
  wp[0].to_color_space(prev);
}


char
FastReadCS::WordPair::get_char(size_t mask, size_t pos) const {
  // 00 -> A, 01 -> C, 10 -> G, 11 -> T
  const size_t selector = (rmap_bits::low_bit << (pos - 1));
  if ((mask & bads) & selector) return 'N';
  const bool upper_bit = ((mask & upper) & selector);
  const bool lower_bit = ((mask & lower) & selector);
  if (upper_bit) return (lower_bit) ? 'T' : 'G';
  else return (lower_bit) ? 'C' : 'A';
}


string
FastReadCS::WordPair::bits2string(size_t mask, size_t bits) {
  string s;
  size_t selector = rmap_bits::high_bit;
  for (size_t i = 0; i < rmap_bits::word_size; ++i) {
    s += (selector & bits & mask) ? '1' : '0';
    selector >>= 1;
  }
  return s;
}


string
FastReadCS::WordPair::tostring_bits(size_t mask) const {
  const string s(bits2string(mask, upper) + "\n" +
		 bits2string(mask, lower) + "\n" + 
		 bits2string(rmap_bits::all_ones, bads) + "\n");
  string seq;
  for (size_t i = rmap_bits::word_size; i > 0; --i)
    seq += get_char(mask, i);
  return s + seq;
}

string
FastReadCS::WordPair::tostring_bases(size_t mask) const {
  string seq;
  for (size_t i = rmap_bits::word_size; i > 0; --i)
    seq += get_char(mask, i);
  return seq;
}

////////////////////////////////////////////////////////////////////////
// FAST READ

size_t FastReadCS::score_mask = 0;
size_t FastReadCS::segments = 0;
size_t FastReadCS::read_width = 0;
size_t FastReadCS::right_most_bit = 0;

void
FastReadCS::set_read_width(const size_t m) {
  read_width = m;
  score_mask = (rmap_bits::low_bit << (m % rmap_bits::word_size)) - 1;
  right_most_bit = (rmap_bits::word_size - (m % rmap_bits::word_size));
  score_mask <<= right_most_bit;
  segments = size_t(std::ceil(m/static_cast<float>(rmap_bits::word_size))) - 1;
}

FastReadCS::FastReadCS(const std::string &s_) {
  assert(s_.length() > 0);
  wp.resize(segments + 1);
  for (size_t i = 0; i < segments; ++i) {    
    const string this_seg(s_.substr(i*rmap_bits::word_size, rmap_bits::word_size));
    wp[i] = WordPair(this_seg);
  }
  wp[segments] = WordPair(s_.substr(segments*rmap_bits::word_size));
}

FastReadCS::FastReadCS(std::string::const_iterator a,
		   const std::string::const_iterator b) {
  wp.resize(segments + 1);
  for (size_t i = 0; i < segments; ++i) {
    const string this_seg(a + i*rmap_bits::word_size, 
			  a + (i + 1)*rmap_bits::word_size);
    wp[i] = WordPair(this_seg);
  }
  wp[segments] = WordPair(string(a + segments*rmap_bits::word_size, b));
}

string
FastReadCS::tostring_bases() const {
  std::ostringstream ss;
  for (size_t i = 0; i < segments; ++i)
    ss << wp[i].tostring_bases(rmap_bits::all_ones);
  ss << wp[segments].tostring_bases(score_mask);
  return ss.str();
}

string
FastReadCS::tostring_bits() const {
  std::ostringstream ss;
  for (size_t i = 0; i < segments; ++i)
    ss << wp[i].tostring_bits(rmap_bits::all_ones) << std::endl;
  ss << wp[segments].tostring_bits(score_mask) << std::endl;
  return ss.str();
}
