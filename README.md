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
[21:38:20] Range Start: 61 bits
[21:38:20] Range End  : 62 bits
[21:38:20] Block Width: 2^31
[21:38:20] Search Pub : 02034eb1de50492bf0c26998f8dacfb4a8580601de3f2a9b4ced52bcc79a32ebfc
[21:38:20] Stride_sum written to file
[21:38:20] Creating bloomfilter image
[21:59:27] Writing image to bloom.bf
[22:00:45] Elapsed time: (0)hours (22)minutes (25)seconds
                                                                                                            
┌──(alexander㉿DESKTOP-ZMYIN9A))-[~/Documents/Key_Search]
└─$ ./key_search
[22:01:03] S_table generated
[22:01:03] Range Start: 61 bits
[22:01:03] Range End  : 62 bits
[22:01:03] Block Width: 2^31
[22:01:03] Search Pub : 02034eb1de50492bf0c26998f8dacfb4a8580601de3f2a9b4ced52bcc79a32ebfc
[22:01:03] Loading Bloomfilter image
[22:01:09] Key Search in progress...
[22:01:34] BloomFilter Hit (C+) -> False Positive
[22:01:39] BloomFilter Hit (S+) -> False Positive
[22:01:58] BloomFilter Hit (C-) -> Success
[22:01:58] Private key: 4415445441513132075
[22:01:58] Elapsed time: (0)hours (0)minutes (49)seconds

</pre>
