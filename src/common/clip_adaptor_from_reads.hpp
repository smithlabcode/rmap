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

#ifndef CLIP_ADAPTOR_FROM_READS_HPP
#define CLIP_ADAPTOR_FROM_READS_HPP

#include "smithlab_utils.hpp"

#include <string>
#include <limits>


const size_t MIN_ADAPTOR_MATCH_SCORE = 10;

size_t 
clip_adaptor_from_read(const std::string &adaptor, 
		       const size_t min_match_score, std::string &s);

#endif
