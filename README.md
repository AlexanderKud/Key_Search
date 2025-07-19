<pre>
Requires C++ Boost Library to be installed. <a href="https://www.boost.org">https://www.boost.org</a>
Requires C/C++ GMP Library to be installed. <a href="https://gmplib.org">https://gmplib.org</a>
Requires C/C++ OpenMP Library to be installed. <a href="https://www.openmp.org">https://www.openmp.org</a>
You can find the right package manager commands on the Internet for your Linux Distro.

prep_bloom.cpp
- batch addition
- batch inversion
- calculating x coordinate for the batch - 1
- calculating x,y for the last of the batch entry (used as the next startPoint)
- bloom add x coordinate uint64_t bits64[3] part

key_search.cpp
- batch addition/subtraction
- batch inversion
- center of the partition and start/end of partition as the starting points
- calculating x coordinate for the batch - 1
- calculating x,y for the last of the batch entry (used as the next startPoint)
- bloom check x coordinate uint64_t bits64[3] part

Kali Linux XFCE
┌──(alexander㉿DESKTOP-ZMYIN9A)-[~/Documents/Key_Search]
└─$ ./prep_bloom
[23:33:10] Range Start: 63 bits
[23:33:10] Range End  : 64 bits
[23:33:10] Block Width: 2^32
[23:33:10] Search Pub : 02072a2a4dcdcd70f98b22900a08807258f2c4397044d493016bb86cff2fc35d57
[23:33:10] Stride_sum written to file
[23:33:10] Creating bloomfilter image
[00:20:35] Writing image to bloom.bf
[00:24:00] Elapsed time: (0)hours (50)minutes (51)seconds
                                                                                                            
┌──(alexander㉿DESKTOP-ZMYIN9A))-[~/Documents/Key_Search]
└─$ ./key_search
[00:30:25] S_table generated
[00:30:25] Range Start: 63 bits
[00:30:25] Range End  : 64 bits
[00:30:25] Block Width: 2^32
[00:30:25] Search Pub : 02072a2a4dcdcd70f98b22900a08807258f2c4397044d493016bb86cff2fc35d57
[00:30:25] Loading Bloomfilter image
[00:32:34] Key Search in progress...
[00:33:20] BloomFilter Hit (C-) -> False Positive
[00:34:07] BloomFilter Hit (S+) -> False Positive
[00:34:16] BloomFilter Hit (S-) -> Success
[00:34:16] Private key: 9743243711422511270
[00:34:16] Elapsed time: (0)hours (1)minutes (42)seconds
</pre>
