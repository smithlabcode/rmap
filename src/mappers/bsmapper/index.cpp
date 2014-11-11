/* The detail description of each function please refer to head file */
#include "index.hpp"

void BuildIndex::CountBucketSize() {
  cerr << "[COUNT HASH BUCKET SIZE]" << endl;
  hash_table->counter_size = 1 << (2 * HASHLEN);
  MEMORY_ALLOCATE_CHECK(
      hash_table->counter = (uint32_t * ) malloc(
          (hash_table->counter_size + 1) * sizeof(uint32_t)));
  memset(hash_table->counter, 0,
         (hash_table->counter_size + 1) * sizeof(uint32_t));

  uint32_t size = 0, hashValue = 0;
  for (uint16_t i = 0; i < genome->num_of_chroms; ++i) {
    cerr << "[" << i + 1 << "/" << genome->num_of_chroms << "]";
    if (genome->chrom_sizes[i] + genome->chrom_start_pos[i] < HASHLEN)
      continue;
    size = genome->chrom_sizes[i] + genome->chrom_start_pos[i] - HASHLEN;
    for (uint32_t j = genome->chrom_start_pos[i]; j <= size; ++j) {
      hashValue = getHashValue(&(genome->chrom_seqs[j]));
      hash_table->counter[hashValue]++;
    }
  }
  cerr << endl;

  for (uint32_t i = 1; i <= hash_table->counter_size; ++i) {
    hash_table->counter[i] += hash_table->counter[i - 1];
  }
  hash_table->index_size = hash_table->counter[hash_table->counter_size];
  cerr << "[THE SIZE OF INDEX ARRAY IS " << hash_table->index_size << "]"
       << endl;

  for (uint32_t i = hash_table->counter_size - 1; i >= 1; --i) {
    hash_table->counter[i] = hash_table->counter[i - 1];
  }
  hash_table->counter[0] = 0;
}

void BuildIndex::HashToBucket() {
  cerr << "[HASH TO BUCKET]" << endl;
  cerr << "[THE MEMORY OF HASH TABLE INDEX ARRAY IS "
       << sizeof(uint32_t) * (hash_table->index_size / GB) << " GB]" << endl;

  MEMORY_ALLOCATE_CHECK(
      hash_table->index = (uint32_t * ) malloc(
          sizeof(uint32_t) * hash_table->index_size));
  uint32_t size = 0, hashValue = 0;
  for (uint16_t i = 0; i < genome->num_of_chroms; ++i) {
    cerr << "[" << i + 1 << "/" << genome->num_of_chroms << "]";
    if (genome->chrom_sizes[i] + genome->chrom_start_pos[i] < HASHLEN)
      continue;
    size = genome->chrom_sizes[i] + genome->chrom_start_pos[i] - HASHLEN;
    for (uint32_t j = genome->chrom_start_pos[i]; j <= size; ++j) {
      hashValue = getHashValue(&(genome->chrom_seqs[j]));
      hash_table->index[hash_table->counter[hashValue]++] = j;
    }
  }
  cerr << endl;

  for (uint32_t i = hash_table->counter_size - 1; i >= 1; --i) {
    hash_table->counter[i] = hash_table->counter[i - 1];
  }
  hash_table->counter[0] = 0;
}

void BuildIndex::WriteIndex() {
  FILE * fout = fopen(index_file.c_str(), "wb");
  cerr << "[WRITTING INDEX TO " << index_file << "]" << endl;
  uint32_t chrom_name_len;
  fwrite(&genome->num_of_chroms, sizeof(uint16_t), 1, fout);
  for (uint16_t i = 0; i < genome->num_of_chroms; ++i) {
    chrom_name_len = genome->chrom_names[i].size();
    if (chrom_name_len > 255)
      chrom_name_len = 255;
    fwrite(&chrom_name_len, sizeof(uint32_t), 1, fout);
    fwrite(&(genome->chrom_names[i][0]), sizeof(char), chrom_name_len, fout);
    fwrite(&(genome->chrom_sizes[i]), sizeof(uint32_t), 1, fout);
    fwrite(&(genome->chrom_start_pos[i]), sizeof(uint32_t), 1, fout);
  }
  fwrite(genome->chrom_seqs, sizeof(char), genome->all_chroms_len, fout);

  fwrite(&(hash_table->counter_size), sizeof(uint32_t), 1, fout);
  fwrite(&(hash_table->index_size), sizeof(uint32_t), 1, fout);
  fwrite(hash_table->counter, sizeof(uint32_t), hash_table->counter_size + 1,
         fout);
  fwrite(hash_table->index, sizeof(uint32_t), hash_table->index_size, fout);

  fclose(fout);
}

ReadIndex::ReadIndex(const string& _index_file, Genome* _genome,
                     HashTable* _hash_table)
    : index_file(_index_file),
      genome(_genome),
      hash_table(_hash_table) {

  cerr << "READING INDEX FROM " << index_file << endl;
  FILE * fin = fopen(index_file.c_str(), "rb");
  FILE_OPEN_CHECK(fin);
  FREAD_CHECK(fread(&genome->num_of_chroms, sizeof(uint16_t), 1, fin), 1);
  uint32_t tmp;
  char chrom_name[256];
  genome->all_chroms_len = 0;
  for (uint16_t i = 0; i < genome->num_of_chroms; ++i) {
    FREAD_CHECK(fread(&tmp, sizeof(uint32_t), 1, fin), 1);
    FREAD_CHECK(fread(chrom_name, sizeof(char), tmp, fin), tmp);
    chrom_name[tmp] = 0;
    genome->chrom_names.push_back(chrom_name);
    FREAD_CHECK(fread(&tmp, sizeof(uint32_t), 1, fin), 1);
    genome->chrom_sizes.push_back(tmp);
    FREAD_CHECK(fread(&tmp, sizeof(uint32_t), 1, fin), 1);
    genome->chrom_start_pos.push_back(tmp);
    genome->all_chroms_len += genome->chrom_sizes[i];
  }
  MEMORY_ALLOCATE_CHECK(
      genome->chrom_seqs = (char * ) malloc(
          sizeof(char) * (genome->all_chroms_len + 1)));
  FREAD_CHECK(
      fread(genome->chrom_seqs, sizeof(char), genome->all_chroms_len, fin),
      genome->all_chroms_len);
  genome->chrom_seqs[genome->all_chroms_len] = 0;

  FREAD_CHECK(fread(&(hash_table->counter_size), sizeof(uint32_t), 1, fin), 1);
  FREAD_CHECK(fread(&(hash_table->index_size), sizeof(uint32_t), 1, fin), 1);
  MEMORY_ALLOCATE_CHECK(
      hash_table->counter = (uint32_t * ) malloc(
          sizeof(uint32_t) * (hash_table->counter_size + 1)));
  FREAD_CHECK(
      fread(hash_table->counter, sizeof(uint32_t), hash_table->counter_size + 1,
            fin),
      hash_table->counter_size + 1);

  MEMORY_ALLOCATE_CHECK(
      hash_table->index = (uint32_t * ) malloc(
          sizeof(uint32_t) * hash_table->index_size));
  FREAD_CHECK(
      fread(hash_table->index, sizeof(uint32_t), hash_table->index_size, fin),
      hash_table->index_size);
  fclose(fin);
}
