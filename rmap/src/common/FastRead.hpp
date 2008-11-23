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

#ifndef FAST_READ_HPP
#define FAST_READ_HPP

#include <string>
#include <vector>
#include <iostream>
#include <ostream>
#include "rmap_utils.hpp"

class FastRead {
public:
  FastRead(const std::string &s);
  FastRead(std::string::const_iterator a,
	   const std::string::const_iterator b);
  FastRead() {
    wp.resize(segments);
    wp.back() = WordPair((rmap_bits::all_ones << right_most_bit));
  }
  std::string tostring_bases() const;
  std::string tostring_bits() const;
  size_t score(const FastRead &other) const;
  void shift(const size_t i);
  static void set_read_width(const size_t m);
private:
  struct WordPair {
  public:
    WordPair() : upper(0), lower(0), bads(rmap_bits::all_ones) {}
    WordPair(size_t bads_mask) : upper(0), lower(0), bads(bads_mask) {}
    WordPair(const std::string &s);
    char get_char(size_t mask, size_t pos) const;
    void shift(const size_t i);
    void shift(const size_t i, const size_t shifter);
    void shift(const WordPair &other);
    size_t score(const WordPair &other, size_t mask) const;
    std::string tostring_bases(size_t mask) const;
    std::string tostring_bits(size_t mask) const;
  private:
    size_t upper;
    size_t lower;
    size_t bads;
    
    static size_t get_upper(const size_t i) {return i > 1;}
    static size_t get_lower(const size_t i) {return (i % 2);}
    static size_t get_bads(char c) {return (c == 4);}
    static std::string bits2string(size_t mask, size_t bits);
  };

  std::vector<WordPair> wp;
  static size_t score_mask;
  static size_t segments;
  static size_t read_width;
  static size_t right_most_bit;
};

inline 
std::ostream& 
operator<<(std::ostream& s, const FastRead& fr) {
  return s << fr.tostring_bits();
}

inline void
FastRead::WordPair::shift(const size_t i) {
  upper = ((upper << 1) + (i > 1));
  lower = ((lower << 1) + (i % 2));
  bads  = ((bads  << 1) + (i == 4));
}

inline void
FastRead::WordPair::shift(const size_t i, const size_t shifter) {
  upper = ((upper << 1) + (static_cast<size_t>((i > 1))  << shifter));
  lower = ((lower << 1) + (static_cast<size_t>((i % 2))  << shifter));
  bads  = ((bads << 1)  + (static_cast<size_t>((i == 4)) << shifter));
}

inline void
FastRead::shift(const size_t idx) {
  for (size_t i = 0; i < wp.size() - 1; ++i)
    wp[i].shift(wp[i + 1]);
  wp.back().shift(idx, right_most_bit);
}

inline size_t
FastRead::score(const FastRead &other) const {
  size_t ss = 0;
  for (size_t i = 0; i < wp.size() - 1; ++i)
    ss += wp[i].score(other.wp[i], rmap_bits::all_ones);
  return ss + wp.back().score(other.wp.back(), score_mask);
}

inline size_t
FastRead::WordPair::score(const FastRead::WordPair &other, const size_t score_mask) const {
  register size_t bits = ((other.upper ^ upper) | 
			  (other.lower ^ lower) | other.bads | bads) & score_mask;
  bits = ((bits & 0xAAAAAAAAAAAAAAAAul) >> 1)  + (bits & 0x5555555555555555ul);
  bits = ((bits & 0xCCCCCCCCCCCCCCCCul) >> 2)  + (bits & 0x3333333333333333ul);
  bits = ((bits & 0xF0F0F0F0F0F0F0F0ul) >> 4)  + (bits & 0x0F0F0F0F0F0F0F0Ful);
  bits = ((bits & 0xFF00FF00FF00FF00ul) >> 8)  + (bits & 0x00FF00FF00FF00FFul);
  bits = ((bits & 0xFFFF0000FFFF0000ul) >> 16) + (bits & 0x0000FFFF0000FFFFul);
  return ((bits & 0xFFFFFFFF00000000ul) >> 32) + (bits & 0x00000000FFFFFFFFul);
}

inline void
FastRead::WordPair::shift(const FastRead::WordPair &other) {
  upper = ((upper << 1) | static_cast<size_t>((other.upper & rmap_bits::high_bit) != 0));
  lower = ((lower << 1) | static_cast<size_t>((other.lower & rmap_bits::high_bit) != 0));
  bads  = ((bads << 1)  | static_cast<size_t>((other.bads  & rmap_bits::high_bit) != 0));
}

#endif
