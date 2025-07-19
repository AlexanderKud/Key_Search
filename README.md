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

key_search.cpp
- batch addition/subtraction
- batch inversion
- center of the partition and start/end of partition as the starting points
- calculating just x coordinate for the batch - 1
- calculating x,y for the last of the batch entry (used as the next startPoint)
- bloom add only x coordinate uint64_t bits64[3] part

Kali Linux XFCE
┌──(alexander㉿DESKTOP-ZMYIN9A)-[~/Documents/Key_Search]
└─$ ./prep_bloom
[22:17:40] Range Start: 62 bits
[22:17:40] Range End  : 63 bits
[22:17:40] Block Width: 2^32
[22:17:40] Search Pub : 02fa8ea2106e6ffb35c3aca7a01d17dd6ef5f2ed82a9386419f59d32b4ea81cd8c
[22:17:40] Stride_sum written to file
[22:17:40] Creating bloomfilter image
[23:21:05] Writing image to bloom.bf
[23:24:57] Elapsed time: (1)hours (7)minutes (19)seconds
                                                                                                            
┌──(alexander㉿DESKTOP-ZMYIN9A))-[~/Documents/Key_Search]
└─$ ./key_search
[23:25:46] S_table generated
[23:25:46] Range Start: 62 bits
[23:25:46] Range End  : 63 bits
[23:25:46] Block Width: 2^32
[23:25:46] Search Pub : 02fa8ea2106e6ffb35c3aca7a01d17dd6ef5f2ed82a9386419f59d32b4ea81cd8c
[23:25:46] Loading Bloomfilter image
[23:27:55] Key Search in progress...
[23:28:38] BloomFilter Hit (S-) -> Success
[23:28:38] Private key: 7759623732764795293
[23:28:38] Elapsed time: (0)hours (0)minutes (43)seconds

</pre>
