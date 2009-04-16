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

#ifndef MAP_RESULT_PE_HPP
#define MAP_RESULT_PE_HPP

#include <numeric>

struct MapResultPE {
  MapResultPE(size_t scr,
	      size_t chr = std::numeric_limits<size_t>::max(),
	      size_t ste = std::numeric_limits<size_t>::max(),
	      size_t st2 = std::numeric_limits<size_t>::max(),
	      size_t str = true, size_t unq = true) :
    chrom(chr), site(ste), site2(st2), score(scr), strand(str), unique(unq) {}
  void set(size_t scr, size_t chr, size_t ste, size_t st2, size_t str);
  unsigned chrom  : 15;
  unsigned site   : 32;
  unsigned site2  : 32;
  unsigned score  : 15;
  unsigned strand : 1;
  unsigned unique : 1;
  bool operator<(const MapResultPE& rhs) const {
    return (chrom < rhs.chrom ||
            (chrom == rhs.chrom && site < rhs.site) ||
            (chrom == rhs.chrom && site == rhs.site && site2 < rhs.site2));
  }
  bool operator==(const MapResultPE& rhs) const {
    return chrom == rhs.chrom && site == rhs.site && site2 == rhs.site2;
  }
};

inline void 
MapResultPE::set(size_t scr, size_t chr, size_t ste, size_t st2, size_t str) {
  unique = (scr < score || (site == ste && st2 == site2 && chrom == chr && unique));
  chrom = chr;
  site = ste;
  site2 = st2;
  score = scr;
  strand = str;
}

std::ostream&
operator<<(std::ostream &os, const MapResultPE &mrp) {
  return os << "CHROM=" << mrp.chrom << "\t" 
	    << "SITE =" << mrp.site << "\t"
	    << "SITE2=" << mrp.site2 << "\t"
	    << "SCORE=" << mrp.score << "\t"
	    << "STRAN=" << mrp.strand << "\t"
	    << "UNIQ =" << mrp.unique;
}

struct MultiMapResultPE {
  MultiMapResultPE(size_t scr) : score(scr) {}
  void add(size_t scr, size_t chr, size_t ste, size_t st2, size_t str) {
    if (scr < score) {
      mr.clear();
      score = scr;
    }
    if (mr.size() <= twice_max_count)
      mr.push_back(MapResultPE(scr, chr, ste, st2, str));
  }
  std::vector<MapResultPE> mr;
  size_t score;
  bool ambiguous() const {return mr.size() > max_count;}
  void collapse() {
    sort(mr.begin(), mr.end());
    mr.erase(std::unique(mr.begin(), mr.end()), mr.end());
  }
  static size_t max_count;
  static size_t twice_max_count;
  static void init(const size_t mc) {
    max_count = mc;
    twice_max_count = 2*mc;
  }
};

size_t MultiMapResultPE::max_count;
size_t MultiMapResultPE::twice_max_count;

#endif
