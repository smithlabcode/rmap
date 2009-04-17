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

#include "rmap_utils.hpp"
#include "deadzone_utils.hpp"
#include <sstream>
#include <algorithm>
#include <vector>

using std::string;
using std::vector;
using std::min;

size_t long_index::index_size = 0;

static const size_t THIRTY_TWO = static_cast<size_t>(32);
static const size_t SIXTY_FOUR = static_cast<size_t>(64);

long_index::long_index(const string &s, bool revcomp) {
  const size_t len = s.length();
  if (revcomp) {
    first = mer2i_rc(s.end() - min(THIRTY_TWO, len), s.end());
    if (len > THIRTY_TWO)
      second = mer2i_rc(s.begin(), s.end() - min(THIRTY_TWO, len));
    else second = 0;
  }
  else {
    first = mer2i(s.begin(), s.begin() + min(THIRTY_TWO, len));
    if (len > THIRTY_TWO)
      second = mer2i(s.begin() + THIRTY_TWO, s.end());
    else second = 0;
  }
}

long_index::long_index(string::const_iterator a, 
		       string::const_iterator b, bool revcomp) {
  const size_t len = b - a;
  if (revcomp) {
    first = mer2i_rc(b - min(THIRTY_TWO, len), b);
    if (len > THIRTY_TWO)
      second = mer2i_rc(a, b - min(THIRTY_TWO, len));
    else second = 0;
  }
  else {
    first = mer2i(a, a + min(THIRTY_TWO, len));
    if (len > THIRTY_TWO)
      second = mer2i(a + THIRTY_TWO, b);
    else second = 0;
  }
}

static string 
bits2string(size_t positions, size_t bits) {
  static const size_t high_bit = (static_cast<size_t>(1) << (SIXTY_FOUR - 1));
  static const size_t word_size = SIXTY_FOUR;
  std::string s;
  size_t selector = high_bit;
  for (size_t i = 0; i < word_size; ++i) {
    s += (selector & bits) ? '1' : '0';
    selector >>= 1;
  }
  return s.substr(s.length() - positions);
}

string
long_index::tostring() const {
  std::ostringstream ss;
  ss << bits2string(THIRTY_TWO, first) << bits2string(THIRTY_TWO, second);
  return ss.str();
}

string
long_index::tostring_bases() const {
  std::ostringstream ss;
  ss << i2mer(min(THIRTY_TWO, index_size), first);
  if (index_size > THIRTY_TWO)
    ss << i2mer(index_size - THIRTY_TWO, second);
  return ss.str();
}

inline static int
base2int_bs_l(char b) {
  static const int b2i_bs_size = 20;
  static const int b2i_bs[] = {
  //A, b, C, d, e, f, G, h, i, j, k, l, m, n, o, p, q, r, s, T
    0,-1, 3,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 3
  };
  b = std::toupper(b);
  if (b - 'A' >= 0 && b - 'A' < b2i_bs_size)
    return b2i_bs[b - 'A'];
  else return -1;
}

inline static int
base2int_bs_rc_l(char b) {
  static const int b2i_bs_size = 20;
  static const int b2i_bs_rc[] = {
  //A, b, C, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, T
    3,-1, 2,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0
  };
  b = std::toupper(b);
  if (b - 'A' >= 0 && b - 'A' < b2i_bs_size)
    return b2i_bs_rc[b - 'A'];
  else return -1;
}


inline static size_t
mer2i_bs(const std::string::const_iterator a, std::string::const_iterator b) {
  size_t multiplier = 1, index = 0;
  do {
    --b;
    index += base2int_bs_l(*b)*multiplier;
    multiplier *= rmap::alphabet_size;
  } while (b > a);
  return index;
}

inline static size_t
mer2i_bs_rc(std::string::const_iterator a, const std::string::const_iterator b) {
  size_t multiplier = 1, index = 0;
  do {
    index += base2int_bs_rc_l(*a)*multiplier;
    multiplier *= rmap::alphabet_size;
  } while (++a < b);
  return index;
}

long_index
get_long_index_bs(string::const_iterator a, string::const_iterator b, 
		  bool revcomp) {
  const size_t len = b - a;
  size_t first = 0, second = 0;
  if (revcomp) {
    first = mer2i_bs_rc(b - min(THIRTY_TWO, len), b);
    if (len > THIRTY_TWO)
      second = mer2i_bs_rc(a, b - min(THIRTY_TWO, len));
  }
  else {
    first = mer2i_bs(a, a + min(THIRTY_TWO, len));
    if (len > THIRTY_TWO)
      second = mer2i_bs(a + THIRTY_TWO, b);
  }
  return long_index(first, second);
}



inline static size_t
mer2i_bs_meth(const std::string::const_iterator a, std::string::const_iterator b) {
  size_t multiplier = 1, index = 0;
  --b;
  index += base2int(*b)*multiplier;
  multiplier *= rmap::alphabet_size;
  do {
    --b;
    index += (*b == 'C' && *(b + 1) == 'G') ?
      base2int(*b)*multiplier : base2int_bs_l(*b)*multiplier;
    multiplier *= rmap::alphabet_size;
  } while (b > a);
  return index;
}

inline static size_t
mer2i_bs_meth_rc(std::string::const_iterator a, const std::string::const_iterator b) {
  size_t multiplier = 1, index = 0;
  const std::string::const_iterator c(b - 1);
  do {
    index += (*a == 'G' && *(a - 1) == 'C') ? base2int_rc(*a)*multiplier :
      base2int_bs_rc_l(*a)*multiplier;
    multiplier *= rmap::alphabet_size;
  } while (++a < c);
  index += base2int_rc(*a)*multiplier;
  return index;
}

long_index
get_long_index_bs_meth(string::const_iterator a, string::const_iterator b, 
		       bool revcomp) {
  const size_t len = b - a;
  size_t first = 0, second = 0;
  if (revcomp) {
    first = mer2i_bs_meth_rc(b - min(THIRTY_TWO, len), b);
    if (len > THIRTY_TWO)
      second = mer2i_bs_meth_rc(a, b - min(THIRTY_TWO, len));
  }
  else {
    first = mer2i_bs_meth(a, a + min(THIRTY_TWO, len));
    if (len > THIRTY_TWO)
      second = mer2i_bs_meth(a + THIRTY_TWO, b);
  }
  return long_index(first, second);
}
