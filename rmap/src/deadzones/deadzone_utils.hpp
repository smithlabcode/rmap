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

#ifndef DEAD_ZONE_UTILS
#define DEAD_ZONE_UTILS

#include <cstdlib>
#include <string>
#include <ostream>

#include <tr1/unordered_map>
#include "rmap_utils.hpp"

struct long_index {
  long_index() : first(0), second(0) {}
  long_index(const size_t a, const size_t b) : first(a), second(b) {}
  long_index(std::string::const_iterator a,
	     std::string::const_iterator b, bool revcomp = false);
  long_index(const std::string &s, bool revcomp = false);
  size_t first;
  size_t second;
  std::string tostring() const;
  std::string tostring_bases() const;
  static void set_index_size(size_t i) {index_size = i;}
  static size_t index_size;
  
  bool operator==(const long_index &rhs) const {
    return first == rhs.first && second == rhs.second;}
  bool operator!=(const long_index &rhs) const {
    return first != rhs.first || second != rhs.second;}
  bool operator<(const long_index &rhs) const {
    return first < rhs.first || (first == rhs.first && second < rhs.second);}
};

namespace std {
  namespace tr1 {
    template<> struct hash<long_index> {
      size_t operator()(long_index __x) const { 
	hash<size_t> h;
	return h(__x.first) ^ h(__x.second);
      }
    };
  };
};

inline std::ostream& 
operator<<(std::ostream& s, const long_index& li) {return s << li.tostring();}

long_index
get_long_index_bs(std::string::const_iterator a, std::string::const_iterator b, 
		  bool revcomp = false);

long_index
get_long_index_bs_meth(std::string::const_iterator a, 
		       std::string::const_iterator b, 
		       bool revcomp = false);

#endif
