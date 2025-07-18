<pre>
Requires C++ Boost Library to be installed. <a href="https://www.boost.org">https://www.boost.org</a>
Requires C/C++ GMP Library to be installed. <a href="https://gmplib.org">https://gmplib.org</a>
Requires C/C++ OpenMP Library to be installed. <a href="https://www.openmp.org">https://www.openmp.org</a>
You can find the right package manager commands on the Internet for your Linux Distro.

prep_bloom.cpp
- batch addition
- batch inversion
- calculating just x coordinate for the batch - 1
- calculating x,y for the last of the batch entry (used as the next startPoint)
- bloom add only x coordinate uint64_t bits64[3] part

Concept code. Can be made faster with different variations.

key_search.cpp
- batch addition/subtraction
- batch inversion
- center of the partition as the starting point
- calculating just x coordinate for the batch - 1
- calculating x,y for the last of the batch entry (used as the next startPoint)
- bloom add only x coordinate uint64_t bits64[3] part

alexander@alexander-home:~/Documents/Key_Search$ ./prep_bloom
[20:54:44] Range Start: 60 bits
[20:54:44] Range End  : 61 bits
[20:54:44] Block Width: 2^30
[20:54:44] Search Pub : 02bfc4c349b81e87f2f5597252877bf074cfcec7372ff7f55264a916ae5f7b82f1
[20:54:44] Stride_sum written to file
[20:54:44] Creating bloomfilter image
[21:26:36] Writing image to bloom.bf
[21:26:55] Elapsed time: (0)hours (32)minutes (11)seconds

alexander@alexander-home:~/Documents/Key_Search$ ./key_search
[21:39:23] S_table generated
[21:39:23] Range Start: 60 bits
[21:39:23] Range End  : 61 bits
[21:39:23] Block Width: 2^30
[21:39:23] Search Pub : 02bfc4c349b81e87f2f5597252877bf074cfcec7372ff7f55264a916ae5f7b82f1
[21:39:23] Loading Bloomfilter image
[21:39:28] Key Search in progress...
[21:40:04] BloomFilter Hit (+) -> Success
[21:40:04] Private key: 1754855292826764224
[21:40:04] Elapsed time: (0)hours (0)minutes (36)seconds
</pre>
