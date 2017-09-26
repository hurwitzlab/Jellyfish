/*  This file is part of Jellyfish.

    Jellyfish is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Jellyfish is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Jellyfish.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <map>
#include <string>
#include <vector>

#include <jellyfish/err.hpp>
#include <jellyfish/thread_exec.hpp>
#include <jellyfish/file_header.hpp>
#include <jellyfish/stream_manager.hpp>
#include <jellyfish/whole_sequence_parser.hpp>
#include <jellyfish/mer_dna_bloom_counter.hpp>
#include <jellyfish/jellyfish.hpp>
#include "sequence_mers.hpp"

namespace err = jellyfish::err;

using jellyfish::mer_dna;
using jellyfish::mer_dna_bloom_counter;
typedef jellyfish::whole_sequence_parser<jellyfish::stream_manager<char**> > sequence_parser;

int findMode(std::vector<int> arr)
{
    std::map<int, int> modeMap;
    for (unsigned int i = 0; i < arr.size(); ++i) {
        ++modeMap[arr[i]];
    }

    auto x = std::max_element(modeMap.begin(), modeMap.end(),
        [](const std::pair<int, int>& a, const std::pair<int, int>& b) {
        return a.second < b.second; });

    return x->first;
}

template<typename PathIterator, typename Database>
void query_from_sequence(int min_mode, PathIterator file_begin, PathIterator file_end, const Database& db,
                         bool canonical) {
  jellyfish::stream_manager<PathIterator> streams(file_begin, file_end);
  sequence_parser                         parser(4, 100, 1, streams);
  sequence_mers                           mers(canonical);
  const sequence_mers                     mers_end(canonical);

  while(true) {
    sequence_parser::job j(parser);
    if(j.is_empty()) break;
    for(size_t i = 0; i < j->nb_filled; ++i) {
      mers = j->data[i].seq;
      std::vector<int> kmer_counts;
      if(mers != mers_end) {
        kmer_counts.push_back(db.check(*mers));
        ++mers;
      }
      for( ; mers != mers_end; ++mers) {
        kmer_counts.push_back(db.check(*mers));
        //std::cout << " " << db.check(*mers);
      }
      int mode = findMode(kmer_counts);
      if (mode >= min_mode) {
        std::cout << ">" << j->data[i].header << "\n";
        std::cout << "mode: " << mode << std::endl;
        for (auto i = kmer_counts.begin(); i != kmer_counts.end(); ++i) {
          std::cout << " " << *i;
        }
        std::cout << std::endl;
      }
    }
  }
}

int main(int argc, char *argv[])
{
  if(argc < 4)
    err::die(err::msg() << "Usage: " << argv[0] << " min_mode db.jf file.fa [...]");

  int min_mode = std::stoi(argv[1]);
  std::cout << min_mode << std::endl;
  std::ifstream in(argv[2], std::ios::in|std::ios::binary);
  jellyfish::file_header header(in);
  if(!in.good())
    err::die(err::msg() << "Failed to parse header of file '" << argv[2] << "'");
  mer_dna::k(header.key_len() / 2);
  if(header.format() == "bloomcounter") {
    jellyfish::hash_pair<mer_dna> fns(header.matrix(1), header.matrix(2));
    mer_dna_bloom_counter filter(header.size(), header.nb_hashes(), in, fns);
    if(!in.good())
      err::die("Bloom filter file is truncated");
    in.close();
    query_from_sequence(min_mode, argv + 3, argv + argc, filter, header.canonical());
  } else if(header.format() == binary_dumper::format) {
    jellyfish::mapped_file binary_map(argv[2]);
    binary_query bq(binary_map.base() + header.offset(), header.key_len(), header.counter_len(), header.matrix(),
                    header.size() - 1, binary_map.length() - header.offset());
    query_from_sequence(min_mode, argv + 3, argv + argc, bq, header.canonical());
  } else {
    err::die(err::msg() << "Unsupported format '" << header.format() << "'. Must be a bloom counter or binary list.");
  }

  return 0;
}
