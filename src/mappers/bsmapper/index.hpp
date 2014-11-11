/*
 * This file builds index (hash table) for the reference genome.
 * */

#ifndef INDEX_H_
#define INDEX_H_

#include "util.hpp"
#include "reference.hpp"

#include <utility>
#include <string>

using std::pair;
using std::string;
using std::make_pair;

/*
 * HashTable stores positions in the genome for each k-mer
 */
typedef struct {
  /* counter is an array which indicates the start position of k-mers in index array.
   * So the hash key is the k-mer (transfered to integer), and the hash value is the
   * start position for the partiular k-mer in the index array */
  uint32_t* counter;

  /* index array stores positions for all k-mer */
  uint32_t* index;

  /* counter_size is the size of counter array */
  uint32_t counter_size;

  /* index_size is the size of index array */
  uint32_t index_size;
} HashTable;

class BuildIndex {
 public:
  BuildIndex(const string& _index_file, const Genome* _genome,
             HashTable* _hash_table)
      : index_file(_index_file),
        genome(_genome),
        hash_table(_hash_table) {
    TIME_INFO(CountBucketSize(), "Count Bucket Size");
    TIME_INFO(HashToBucket(), "Hash to Bucket");
    WriteIndex();
  }
  ~BuildIndex() {
    free(hash_table->counter);
    free(hash_table->index);
  }
 private:
  /* count how many a particular k-mer in the genome */
  void CountBucketSize();

  /* put all the k-mer in to the corresponding index array */
  void HashToBucket();

  /* write index for mapping */
  void WriteIndex();

  string index_file;
  const Genome* genome;
  HashTable* hash_table;
};

class ReadIndex {
 public:
  ReadIndex(const string& _index_file, Genome* _genome, HashTable* _hash_table);
  ~ReadIndex() {
    free(genome->chrom_seqs);
    free(hash_table->counter);
    free(hash_table->index);
  }

 private:
  string index_file;
  Genome* genome;
  HashTable* hash_table;
};

#endif /* INDEX_H_ */
