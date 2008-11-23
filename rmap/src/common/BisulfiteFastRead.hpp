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

#ifndef BISULFITE_FAST_READ_HPP
#define BISULFITE_FAST_READ_HPP

#include <string>
#include <vector>
#include <iostream>
#include <ostream>
#include "rmap_utils.hpp"

class BisulfiteFastRead {
public:
  BisulfiteFastRead(const std::string &s);
  BisulfiteFastRead() {
    wp.resize(segments);
    wp.back() = WordPair((rmap_bits::all_ones << right_most_bit));
  }
  std::string tostring_bases() const;
  std::string tostring_bits() const;
  size_t score(const BisulfiteFastRead &other) const;
  void shift(const size_t i);
  static void set_read_width(const size_t m);
private:
  struct WordPair {
  public:
    WordPair() : a_vec(0), c_vec(0), 
		 g_vec(0), t_vec(0), bads(rmap_bits::all_ones) {}
    WordPair(size_t bads_mask) : a_vec(0), c_vec(0), 
				 g_vec(0), t_vec(0), bads(bads_mask) {}
    WordPair(const std::string &s);
    char get_char(size_t mask, size_t pos) const;
    void shift(const size_t i);
    void shift(const size_t i, const size_t shifter);
    void shift(const WordPair &other);
    size_t score(const WordPair &other, size_t mask) const;
    std::string tostring_bases(size_t mask) const;
    std::string tostring_bits(size_t mask) const;
  private:

    size_t a_vec;
    size_t c_vec;
    size_t g_vec;
    size_t t_vec;
    size_t bads;
    
    static size_t contains_a(char c) {return (c == 'a' || c == 'A');}
    static size_t contains_c(char c) {return (c == 'c' || c == 'C');}
    static size_t contains_c_bs(char c) {
      return contains_c(c) || contains_t(c);}
    static size_t contains_g(char c) {return (c == 'g' || c == 'G');}
    static size_t contains_g_bs(char c) {
      return contains_a(c) || contains_g(c);}
    static size_t contains_t(char c) {return (c == 't' || c == 'T');}
    static size_t contains_n(char c) {return (c == 'n' || c == 'N');}

    static size_t contains_a(size_t c) {return (c == 0);}
    static size_t contains_c(size_t c) {return (c == 1);}
    static size_t contains_c_bs(size_t c) {
      return contains_c(c) || contains_t(c);}
    static size_t contains_g(size_t c) {return (c == 2);}
    static size_t contains_g_bs(size_t c) {
      return contains_a(c) || contains_g(c);}
    static size_t contains_t(size_t c) {return (c == 3);}
    static size_t contains_n(size_t c) {return (c == 4);}
    
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
operator<<(std::ostream& s, const BisulfiteFastRead& fr) {
  return s << fr.tostring_bits();
}

inline void
BisulfiteFastRead::WordPair::shift(const size_t i) {
  a_vec = ((a_vec << 1) + (contains_a(i)));
  c_vec = ((c_vec << 1) + (contains_c(i)));
  g_vec = ((g_vec << 1) + (contains_g(i)));
  t_vec = ((t_vec << 1) + (contains_t(i)));
  bads =  ((bads  << 1) + (contains_n(i)));
}

inline void
BisulfiteFastRead::WordPair::shift(const size_t i, const size_t shifter) {
  a_vec = ((a_vec << 1) + (contains_a(i) << shifter));
  c_vec = ((c_vec << 1) + (contains_c(i) << shifter));
  g_vec = ((g_vec << 1) + (contains_g(i) << shifter));
  t_vec = ((t_vec << 1) + (contains_t(i) << shifter));
  bads =  ((bads  << 1) + (contains_n(i) << shifter));
}

inline void
BisulfiteFastRead::shift(const size_t idx) {
  for (size_t i = 0; i < wp.size() - 1; ++i)
    wp[i].shift(wp[i + 1]);
  wp.back().shift(idx, right_most_bit);
}

inline size_t
BisulfiteFastRead::score(const BisulfiteFastRead &other) const {
  size_t ss = 0;
  for (size_t i = 0; i < wp.size() - 1; ++i)
    ss += wp[i].score(other.wp[i], rmap_bits::all_ones);
  return ss + wp.back().score(other.wp.back(), score_mask);
}

inline size_t
BisulfiteFastRead::WordPair::score(const BisulfiteFastRead::WordPair &other, const size_t score_mask) const {
  register size_t bits = (~((other.a_vec & a_vec) | 
			    (other.c_vec & c_vec) | 
			    (other.g_vec & g_vec) | 
			    (other.t_vec & t_vec)) | other.bads | bads) & score_mask;
  bits = ((bits & 0xAAAAAAAAAAAAAAAAul) >> 1)  + (bits & 0x5555555555555555ul);
  bits = ((bits & 0xCCCCCCCCCCCCCCCCul) >> 2)  + (bits & 0x3333333333333333ul);
  bits = ((bits & 0xF0F0F0F0F0F0F0F0ul) >> 4)  + (bits & 0x0F0F0F0F0F0F0F0Ful);
  bits = ((bits & 0xFF00FF00FF00FF00ul) >> 8)  + (bits & 0x00FF00FF00FF00FFul);
  bits = ((bits & 0xFFFF0000FFFF0000ul) >> 16) + (bits & 0x0000FFFF0000FFFFul);
  return ((bits & 0xFFFFFFFF00000000ul) >> 32) + (bits & 0x00000000FFFFFFFFul);
}

#endif
