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

#ifndef MAP_RESULT_HPP
#define MAP_RESULT_HPP

#include <numeric>

struct MapResult {
  MapResult(size_t scr,
	   size_t chr = std::numeric_limits<size_t>::max(),
	   size_t ste = std::numeric_limits<size_t>::max(),
	   size_t str = true, size_t unq = true) :
    chrom(chr), site(ste), score(scr), strand(str), unique(unq) {}
  void set(size_t scr, size_t chr, size_t ste, size_t str);
  unsigned chrom  : 15;
  unsigned site   : 32;
  unsigned score  : 15;
  unsigned strand : 1;
  unsigned unique : 1;
  bool operator<(const MapResult& rhs) const {
    return (chrom < rhs.chrom ||
            (chrom == rhs.chrom && site < rhs.site));
  }
  bool operator==(const MapResult& rhs) const {
    return chrom == rhs.chrom && site == rhs.site;
  }
};

inline void 
MapResult::set(size_t scr, size_t chr, size_t ste, size_t str) {
  unique = (scr < score || (site == ste && unique && chrom == chr));
  chrom = chr;
  site = ste;
  score = scr;
  strand = str;
}

struct MultiMapResult {
  void add(size_t scr, size_t chr, size_t ste, size_t str) {
    if (scr < score) {
      mr.clear();
      score = scr;
    }
    if (mr.size() <= max_count)
      mr.push_back(MapResult(scr, chr, ste, str));
  }
  std::vector<MapResult> mr;
  size_t score;
  static size_t max_count;
  bool is_unique() const {return mr.size() <= max_count;}
  void collapse() {
    sort(mr.begin(), mr.end());
    mr.erase(std::unique(mr.begin(), mr.end()), mr.end());
  }
};

size_t MultiMapResult::max_count;

#endif
