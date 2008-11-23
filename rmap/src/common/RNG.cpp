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

#include "RNG.hpp"
#include <algorithm>
#include<iostream>
using std::cerr;
using std::max;

size_t Runif :: instance_count = 0;

/* most significant w-r bits */
static const unsigned long UPPER_MASK = 0x80000000UL;   

/* least significant r bits */
static const unsigned long LOWER_MASK = 0x7fffffffUL;

static const size_t N = 624;   

static const size_t M = 397;

/*signed long int 
  static const LCG(unsigned long int x)
  {
  return( ((69069 * x) + 1) &0xffffffffUL);
  }*/

unsigned long int 
static const MAGIC(unsigned long int y) {
  return((((y)&0x1) ? 0x9908b0dfUL : 0));
}


void 
set_rng(vector<unsigned long>&randm,size_t &x ,unsigned long int seed) {
  size_t i;
  if(seed ==0)
    seed = 4357; 
  randm[0] = seed & 0xffffffffUL;
  for (i = 1; i < N; i++) {
    randm[i] =
      (1812433253UL * (randm[i-1] ^ (randm[i-1] >> 30)) + i);
    
    randm[i] &= 0xffffffffUL;
  }
  x = i;
}

unsigned long int
pre_compute(vector<unsigned long>&mt,size_t &x) {
  size_t kk;
  if (x >= N) {   /* generate N words at one time */
    for (kk = 0; kk < N - M; kk++) {
      unsigned long y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
      mt[kk] = mt[kk + M] ^ (y >> 1) ^ MAGIC(y);      
    }
    for (; kk<N - 1; kk++) {
      unsigned long y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
      mt[kk] = mt[kk + (M - N)] ^ (y >> 1) ^ MAGIC(y);
    }
    { 
      unsigned long y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
      mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ MAGIC(y);
    }
    x = 0;
  }
  /* Tempering */
  unsigned long int  k = mt[x];
  k ^= (k >> 11);
  k ^= (k << 7) & 0x9d2c5680UL;
  k ^= (k << 15) & 0xefc60000UL;
  k ^= (k >> 18);
  x++; 
  return k;
}

unsigned long int 
rng_integer(vector<unsigned long>&randm,size_t &x,unsigned int n) {
  unsigned long int offset = 0;
  unsigned long int range =  0xffffffffUL - offset;
  unsigned long int scale;
  unsigned long int k;

  if (n > range || n == 0) 
    {
      cerr<<"invalid n, either 0 or exceeds maximum value of generator";
            
    }

  scale = range / n;
  do
    {
      k = (pre_compute(randm,x) - offset) / scale;
    }
  while (k >= n);
   
  return k;
}

double
rng_double(vector<unsigned long>&randm,size_t &x) 
{
  return pre_compute(randm,x)/ 4294967296.0 ;
}


Runif::Runif(size_t seed) : randm(N,0) {
  srand((seed == std::numeric_limits<size_t>::max()) ? 
	((time(0) + getpid())*Runif::getCount()) : seed);
  instance_count ++ ;
  set_rng(randm,x,rand());
}

size_t Runif::getCount()
{
  return (instance_count);
}
Runif::~Runif() {
}

int
Runif::runif(int min_val, int max_val) const {
  return min_val + rng_integer(randm,x,max_val - min_val);
}

size_t
Runif::runif(size_t min_val, size_t max_val) const {
  return min_val  + rng_integer(randm,x,max_val - min_val);
}

double 
Runif::runif(double min_val, double max_val) const {
  return min_val + rng_double(randm,x)*(max_val - min_val);
}
