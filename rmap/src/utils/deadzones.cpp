/*    deadzones: A program for identifying genomic deadzones
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
#include "rmap_os.hpp"
#include "OptionParser.hpp"

#include <numeric>
#include <cmath>

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::cerr;

class IndexLess {
public:  
  IndexLess(const size_t k, const string &s) : 
    kmer(k), itr(s.begin()) {}
  bool operator()(size_t a, size_t b) const {
    const string::const_iterator lim(itr + a + kmer);
    string::const_iterator a_itr(itr + a), b_itr(itr + b);
    while (a_itr < lim && *(a_itr) == *(b_itr)) {
      ++a_itr; ++b_itr;
    }
    return (*a_itr < *b_itr && a_itr < lim);
  }
private:
  size_t kmer;
  const string::const_iterator itr;
};

template <class In> bool
lexico_equal(In first, In last, In first2) {
  while (first != last)
    if (*first++ != *first2++) return false;
  return true;
}

static void
sort_index(const bool VERBOSE, const size_t kmer, const string &prefix,
	   const string &seq, vector<size_t> &ambigs) {
  
  if (VERBOSE) cerr << "[BUILDING INDEX] ";
  vector<size_t> index;
  const string::const_iterator lim(seq.end() - kmer + 1);
  for (string::const_iterator j = seq.begin(); j != lim; ++j)
    if (lexico_equal(prefix.begin(), prefix.end(), j))
      index.push_back(j - seq.begin());
  
  if (!index.empty()) {
  
    if (VERBOSE) cerr << "[SORTING INDEX] ";
    IndexLess index_less(kmer, seq);
    sort(index.begin(), index.end(), index_less);
  
    if (VERBOSE) cerr << "[FINDING DEADS] ";
    const size_t len = seq.length();
    const string::const_iterator start(seq.begin());
    const string::const_iterator end(start + kmer);
    size_t prev = index.front();
    bool prev_inserted = false;
    for (size_t i = 1; i < index.size(); ++i) {
      const size_t curr = index[i];
      if (lexico_equal(start + prev, end + prev, start + curr) && 
	  prev + curr != len) {
	if (!prev_inserted) 
	  ambigs.push_back(prev);
	ambigs.push_back(curr);
	prev_inserted = true;
      
      }
      else prev_inserted = false;
      prev = curr;
    }
  }
  else if (VERBOSE) cerr << "[EMPTY INDEX] ";
}


static void
sort_index(const bool VERBOSE, const bool BISULFITE,
	   const size_t kmer, const size_t prefix_len,
	   const string &seq, vector<size_t> &ambigs) {

  static const float DENOM = CLOCKS_PER_SEC;

  const size_t n_prefix = 
    static_cast<size_t>(pow(rmap::alphabet_size, prefix_len));
  for (size_t i = 0; i < n_prefix; ++i) {
    const string prefix(i2mer(prefix_len, i));
    if (!BISULFITE || prefix.find('C') == string::npos) {
      const clock_t start(clock());
      if (VERBOSE) cerr << "[PREFIX=" << prefix << "] ";
      sort_index(VERBOSE, kmer, prefix, seq, ambigs);
      const clock_t end(clock());
      if (VERBOSE)
	cerr << "[" << (end - start)/DENOM << " SEC] [DONE]" << endl;
    }
  }
}


template <class T> void
merge_contiguous(vector<T> &dead) {
  sort(dead.begin(), dead.end());
  collapse(dead);
  size_t j = 0;
  for (size_t i = 1; i < dead.size(); ++i) {
    if (dead[j].get_end() == dead[i].get_start())
      dead[j].set_end(dead[i].get_end());
    else {
      ++j;
      dead[j].swap(dead[i]);
    }
  }
  ++j;
  dead.erase(dead.begin() + j, dead.end());
}


static void
get_dead(const size_t kmer, 
	 const vector<size_t> &seqoffsets, const vector<size_t> &seqlens, 
	 const vector<string> &chrom_names, vector<size_t> &ambigs, 
	 vector<SimpleGenomicRegion> &dead) {

  const size_t max_offset = seqoffsets.back();
  
  for (size_t i = 0; i < ambigs.size(); ++i) {
    if (ambigs[i] >= max_offset)
      ambigs[i] = 2*max_offset - ambigs[i] - kmer;
    assert(ambigs[i] < max_offset);
  }
  
  sort(ambigs.begin(), ambigs.end());
  ambigs.erase(std::unique(ambigs.begin(), ambigs.end()), ambigs.end());
  
  size_t j = 0;
  for (size_t i = 0; i < ambigs.size(); ++i) {
    while (ambigs[i] >= seqoffsets[j]) ++j;
    const size_t offset = ambigs[i] - (seqoffsets[j] - seqlens[j]);
    dead.push_back(SimpleGenomicRegion(chrom_names[j], offset, offset + 1));
  }
  merge_contiguous(dead);
}


static void
get_dead_bs(const size_t kmer, 
	    const vector<size_t> &seqoffsets, const vector<size_t> &seqlens, 
	    const vector<string> &chrom_names, vector<size_t> &ambigs, 
	    vector<GenomicRegion> &dead) {

  sort(ambigs.begin(), ambigs.end());
  
  const size_t max_offset = seqoffsets.back();
  
  // Do the positive strand bisulfite deadzones
  assert(!ambigs.empty());
  size_t i = 0, j = 0;
  for (; i < ambigs.size() && ambigs[i] < max_offset; ++i) {
    while (ambigs[i] >= seqoffsets[j]) ++j;
    const size_t offset = ambigs[i] - (seqoffsets[j] - seqlens[j]);
    dead.push_back(GenomicRegion(chrom_names[j], offset, offset + 1, "X", 0, '+'));
  }
  merge_contiguous(dead);
  
  // Move the negative strand deadzones into the first portion of the
  // vector and correct their indexes.
  for (j = 0; i < ambigs.size(); ++i)
    ambigs[j++] = 2*max_offset - ambigs[i] - kmer;
  ambigs.erase(ambigs.begin() + j, ambigs.end());
  reverse(ambigs.begin(), ambigs.end());
  
  // Do the negative strand bisulfite deadzones
  vector<GenomicRegion> dead_neg;
  for (i = 0, j = 0; i < ambigs.size(); ++i) {
    while (ambigs[i] >= seqoffsets[j])
      ++j;
    const size_t offset = ambigs[i] - (seqoffsets[j] - seqlens[j]);
    dead_neg.push_back(GenomicRegion(chrom_names[j], offset, offset + 1, "X", 0, '-'));
  }
  merge_contiguous(dead_neg);
  dead.insert(dead.end(), dead_neg.begin(), dead_neg.end());
  dead_neg.clear();
}


// This function appends the reverse complement in a space efficient way
static void
append_revcomp(string &long_seq) {
  const size_t seqlen = long_seq.length();
  long_seq.resize(2*seqlen);
  copy(long_seq.begin(), long_seq.begin() + seqlen, long_seq.begin() + seqlen);
  revcomp_inplace(long_seq.begin() + seqlen, long_seq.end());
}


int
main(int argc, const char **argv) {
  
  try {
    
    // Parameter variables
    size_t kmer = 0;
    size_t prefix_len = 0;
    string outfile;
  
    bool VERBOSE = false;
    bool BISULFITE = false;
  
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse("deadzones", "program for finding deadzones",
			   "<fasta-chrom-files>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      true, outfile);
    opt_parse.add_opt("kmer", 'k', "Width of k-mers", true, kmer);
    opt_parse.add_opt("prefix", 'p', "prefix length", true, prefix_len);
    opt_parse.add_opt("bisulfite", 'B', "get bisulfite deadzones", 
		      false, BISULFITE);
    opt_parse.add_opt("verbose", 'v', "print more run information", 
		      false, VERBOSE);
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.empty()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    vector<string> seqfiles = leftover_args;
    /****************** END COMMAND LINE OPTIONS *****************/
    
    string long_seq;
    vector<size_t> seqlens, seqoffsets;
    vector<string> chrom_names;
    
    if (VERBOSE)
      cerr << "[READING SEQUENCE FILES]" << endl;
    for (size_t i = 0; i < seqfiles.size(); ++i) {
      vector<string> names, sequences;
      read_fasta_file(seqfiles[i].c_str(), names, sequences);
      for (size_t j = 0; j < sequences.size(); ++j) {
	seqlens.push_back(sequences[j].length());
	long_seq += sequences[j];
	seqoffsets.push_back(long_seq.length());
	chrom_names.push_back(names[j]);
      }
      if (VERBOSE)
	cerr << seqfiles[i] << "\t(SEQS: " << names.size() << ")" << endl;
    }
    transform(long_seq.begin(), long_seq.end(), long_seq.begin(),
	      std::ptr_fun(&::toupper));
    
    if (VERBOSE)
      cerr << "[PREPARING CONCATENATED SEQUENCE]" << endl;
    append_revcomp(long_seq);
    
    if (BISULFITE)
      replace(long_seq.begin(), long_seq.end(), 'C', 'T');
    
    if (VERBOSE)
      cerr << "[IDENTIFYING AMBIGUOUS INDEXES]" << endl;
    vector<size_t> ambigs;
    sort_index(VERBOSE, BISULFITE, kmer, prefix_len, long_seq, ambigs);
    long_seq.clear();
    
    if (BISULFITE) {
      if (VERBOSE)
	cerr << "[PREPARING BS DEADZONES]" << endl;
      vector<GenomicRegion> dead;
      get_dead_bs(kmer, seqoffsets, seqlens, chrom_names, ambigs, dead);
      if (VERBOSE)
	cerr << "[WRITING OUTPUT]" << endl;
      std::ofstream out(outfile.c_str());
      copy(dead.begin(), dead.end(), 
	   std::ostream_iterator<GenomicRegion>(out, "\n"));
      out.close();
    }
    else {
      if (VERBOSE)
	cerr << "[PREPARING DEADZONES]" << endl;
      vector<SimpleGenomicRegion> dead;
      get_dead(kmer, seqoffsets, seqlens, chrom_names, ambigs, dead);
      if (VERBOSE)
	cerr << "[WRITING OUTPUT]" << endl;
      std::ofstream out(outfile.c_str());
      copy(dead.begin(), dead.end(), 
	   std::ostream_iterator<SimpleGenomicRegion>(out, "\n"));
      out.close();
    }
    
  }
  catch (RMAPException &e) {
    cerr << "ERROR: " << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
