/*
 *    Part of SMITHLAB software
 *
 *    Copyright (C) 2010 University of Southern California and
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

#include "clip_adaptor_from_reads.hpp"

using std::string;
using std::min;

size_t
similarity(const string &s, const size_t pos, const string &adaptor) {
  const size_t lim = min(min(s.length() - pos, adaptor.length()), 
			 MIN_ADAPTOR_MATCH_SCORE + static_cast<size_t>(3));
  size_t count = 0;
  for (size_t i = 0; i < lim; ++i)
    count += (s[pos + i] == adaptor[i]);
  return count;
}

size_t 
clip_adaptor_from_read(const string &adaptor, 
		       const size_t min_match_score, string &s) {
  const size_t lim = s.length() - min_match_score + 1;
  for (size_t i = 0; i < lim; ++i) {
    const size_t score = similarity(s, i, adaptor);
    if (score >= min_match_score) {
      fill(s.begin() + i, s.end(), 'N');
      return s.length() - i;
    }
  }
  return 0;
}
