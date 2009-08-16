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

#ifndef BISULFITE_FAST_READ_QUALITY_HPP
#define BISULFITE_FAST_READ_QUALITY_HPP

#include <string>
#include <vector>
#include <iostream>
#include <ostream>
#include "rmap_utils.hpp"

class BisulfiteFastReadQuality {
public:
  BisulfiteFastReadQuality(const std::vector<std::vector<double> > &s,
			   bool AG_WILDCARD = false);
  BisulfiteFastReadQuality() {words.resize(segments);}
  std::string tostring_values() const;
  std::string tostring_bits() const;
  std::string tostring() const {return tostring_values() + "\n" + tostring_bits();}
  size_t score(const BisulfiteFastReadQuality &other) const;
  size_t score_ag(const BisulfiteFastReadQuality &other) const {return score(other);}
  void shift(const size_t i);
  static void set_read_width(const size_t rw);
  
  static double value_to_quality(size_t val);
  static size_t quality_to_value(double val);

  static double get_scaler() {return scaler;}
  
private:
  struct Words {
  public:
    Words() : a_vec(0), c_vec(0), g_vec(0), t_vec(0) {}
    Words(const std::vector<std::vector<double> > &s);
    Words(std::string::const_iterator i, const std::string::const_iterator limit);
    char get_char(size_t mask, size_t pos) const;
    void shift_last(const size_t i);
    void shift(const Words &other);

    static size_t get_val(MASK_t mask, MASK_t base_vec, size_t pos);

    size_t score(const Words &other, size_t mask) const;
    std::string tostring_values(size_t mask) const;
    std::string tostring_bits(size_t mask) const;
  private:
    
    size_t a_vec;
    size_t c_vec;
    size_t g_vec;
    size_t t_vec;
    
    static size_t contains_a(size_t c) {return (c == 0) ? low_val_bits : 0;}
    static size_t contains_c(size_t c) {return (c == 1) ? low_val_bits : 0;}
    static size_t contains_g(size_t c) {return (c == 2) ? low_val_bits : 0;}
    static size_t contains_t(size_t c) {return (c == 3) ? low_val_bits : 0;}

    static size_t quality_to_value(double quality);
    static double value_to_quality(size_t val);
    
    static std::string bits2string(size_t mask, size_t bits);
  };
  
  std::vector<Words> words;

  static size_t score_mask;
  static size_t segments;
  static size_t read_width;
  static size_t right_most_bit;
  
  static double scaler;
  static double min_quality_score;
  static double max_quality_score;
  
  static const size_t n_val_bits = 4;
  static const size_t segment_size = 16;
  static const size_t high_val_bits = 0xF000000000000000;
  static const size_t low_val_bits = 0x000000000000000F;
  static const size_t high_val_bits_to_low_val_bits_shift = 60;
  
};

inline 
std::ostream& 
operator<<(std::ostream& s, const BisulfiteFastReadQuality& fr) {
  return s << fr.tostring();
}

inline void
BisulfiteFastReadQuality::Words::shift(const BisulfiteFastReadQuality::Words &other) {
  a_vec <<= n_val_bits;
  c_vec <<= n_val_bits;
  g_vec <<= n_val_bits;
  t_vec <<= n_val_bits;
  a_vec |= ((other.a_vec & BisulfiteFastReadQuality::high_val_bits) >> high_val_bits_to_low_val_bits_shift);
  c_vec |= ((other.c_vec & BisulfiteFastReadQuality::high_val_bits) >> high_val_bits_to_low_val_bits_shift);
  g_vec |= ((other.g_vec & BisulfiteFastReadQuality::high_val_bits) >> high_val_bits_to_low_val_bits_shift);
  t_vec |= ((other.t_vec & BisulfiteFastReadQuality::high_val_bits) >> high_val_bits_to_low_val_bits_shift);
}

inline void
BisulfiteFastReadQuality::Words::shift_last(const size_t i) {
  a_vec = ((a_vec << n_val_bits) + (contains_a(i) << right_most_bit));
  c_vec = ((c_vec << n_val_bits) + (contains_c(i) << right_most_bit));
  g_vec = ((g_vec << n_val_bits) + (contains_g(i) << right_most_bit));
  t_vec = ((t_vec << n_val_bits) + (contains_t(i) << right_most_bit));
}

inline void
BisulfiteFastReadQuality::shift(const size_t idx) {
  std::vector<Words>::iterator i(words.begin());
  std::vector<Words>::const_iterator j(i + 1);
  const std::vector<Words>::const_iterator lim(words.end());
  for (; j < lim; ++i, ++j)
    i->shift(*j);
  i->shift_last(idx);
}


inline size_t
BisulfiteFastReadQuality::score(const BisulfiteFastReadQuality &other) const {
  size_t ss = 0;
  std::vector<Words>::const_iterator i(words.begin());
  std::vector<Words>::const_iterator j(other.words.begin());
  const std::vector<Words>::const_iterator lim(words.end() - 1);
  for (; i < lim; ++i, ++j)
    ss += i->score(*j, rmap_bits::all_ones);
  return ss + i->score(*j, score_mask);
}


inline size_t
BisulfiteFastReadQuality::Words::score(const BisulfiteFastReadQuality::Words &other, const size_t score_mask) const {
  register size_t bits = ((other.a_vec & a_vec) | 
			  (other.c_vec & c_vec) | 
			  (other.g_vec & g_vec) | 
			  (other.t_vec & t_vec)) & score_mask;
  //   bits = ((bits & 0xAAAAAAAAAAAAAAAA) >> 1)  + (bits & 0x5555555555555555);
  //   bits = ((bits & 0xCCCCCCCCCCCCCCCC) >> 2)  + (bits & 0x3333333333333333);
  bits = ((bits & 0xF0F0F0F0F0F0F0F0) >> 4)  + (bits & 0x0F0F0F0F0F0F0F0F);
  bits = ((bits & 0xFF00FF00FF00FF00) >> 8)  + (bits & 0x00FF00FF00FF00FF);
  bits = ((bits & 0xFFFF0000FFFF0000) >> 16) + (bits & 0x0000FFFF0000FFFF);
  return ((bits & 0xFFFFFFFF00000000) >> 32) + (bits & 0x00000000FFFFFFFF);
}

#endif
