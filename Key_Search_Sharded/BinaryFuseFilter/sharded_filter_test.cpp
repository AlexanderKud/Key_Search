#include "sharded_filter.hpp"
#include "binaryfusefilter.h"
#include "filter.hpp"
#include "page.hpp"
#include <iostream>
#include <cstdint>
#include <filesystem>
#include <span>
#include <stdexcept>
#include <vector>

using namespace std;

int main () {
/*
vector<uint64_t> numbers1 = { 3, 1, 2, 6, 4, 5};
vector<uint64_t> numbers2 = { 4, 5, 6};

binfuse::filter8 tiny_low(numbers1);
//binfuse::filter8 tiny_high(numbers2);
//tiny_low = numbers1);
*/

binfuse::filter8 tiny_low(std::vector<std::uint64_t>{
    // note the MSB is clear on all below
    0x0000000000000000,
    0x0000000000000001, // order is not important
    0x7000000000000002,
});

binfuse::filter8 tiny_high(std::vector<std::uint64_t>{
    // note the MSB is set on all below
    0x8000000000000000,
    0x9000000000000001, // order is not important
    0xA000000000000002,
});

binfuse::sharded_filter8_sink sink("sharded_filter8_numbers.bin", 1); // one bit sharding, ie 2 shards

sink.add_shard(tiny_low, 0);  // specify the prefix for each shard
sink.add_shard(tiny_high, 1); // order of adding is not important

cout << "Sink shards: " << sink.shards() << endl;

binfuse::sharded_filter8_source source("sharded_filter8_numbers.bin", 1);

    // verify all entries
cout << source.contains(0x0000000000000000) << endl;
cout << source.contains(0x6000000000000001) << endl;
cout << source.contains(0x7000000000000002) << endl;
cout << source.contains(0x8000000000000000) << endl;
cout << source.contains(0x9000000000000001) << endl;
cout << source.contains(0xB000000000000002) << endl;

std::filesystem::remove("sharded_filter8_numbers.bin");
/*
binfuse::sharded_filter8_source sharded_source;
cout << "Sharded source shards: " << sharded_source.shards() << endl;

binfuse::filter8 tiny_low(std::vector<std::uint64_t>{
    // note the MSB is clear on all below
    0x0000000000000000,
    0x0000000000000001, // order is not important
    0x0000000000000002,
});

binfuse::filter8 tiny_high(std::vector<std::uint64_t>{
    // note the MSB is set on all below
    0x8000000000000000,
    0x8000000000000001, // order is not important
    0x8000000000000002,
});

binfuse::sharded_filter8_sink sink("sharded_filter8_tiny.bin", 1); // one bit sharding, ie 2 shards

sink.add_shard(tiny_low, 0);  // specify the prefix for each shard
sink.add_shard(tiny_high, 1); // order of adding is not important

cout << "Sink shards: " << sink.shards() << endl;
*/
// now reopen the filter as a "source"
/*
binfuse::sharded_filter8_source source("sharded_filter8_tiny.bin", 1);

    // verify all entries
cout << source.contains(0x0000000000000003) << endl;
cout << source.contains(0x0000000000000001) << endl;
cout << source.contains(0x0000000000000002) << endl;
cout << source.contains(0x8000000000000000) << endl;
cout << source.contains(0x8000000000000001) << endl;
cout << source.contains(0x8000000000000002) << endl;
*/
//cout << "Source shards: " << source.shards() << endl;
// cleanup
//std::filesystem::remove("sharded_filter8_tiny.bin");

}
  