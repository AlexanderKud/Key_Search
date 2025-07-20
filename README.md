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
- center of the partition and start/end of partition as the starting points(cross meet)
- calculating x coordinate for the batch - 1
- calculating x,y for the last of the batch entry (used as the next startPoint)
- bloom check x coordinate uint64_t bits64[3] part

Kali Linux XFCE
┌──(alexander㉿DESKTOP-ZMYIN9A)-[~/Documents/Key_Search]
└─$ ./prep_bloom
[07:30:33] Range Start: 65 bits
[07:30:33] Range End  : 66 bits
[07:30:33] Block Width: 2^32
[07:30:33] Search Pub : 02259761f884f0c33054386939b2cece126ae1a4b1ab510212543a10f308945f06
[07:30:33] Stride_sum written to file
[07:30:33] Creating bloomfilter image
[08:17:13] Writing image to bloom.bf
[08:20:44] Elapsed time: (0)hours (50)minutes (11)seconds
                                                                                                            
┌──(alexander㉿DESKTOP-ZMYIN9A))-[~/Documents/Key_Search]
└─$ ./key_search
08:20:49] S_table generated
[08:20:50] Range Start: 65 bits
[08:20:50] Range End  : 66 bits
[08:20:50] Block Width: 2^32
[08:20:50] Search Pub : 02259761f884f0c33054386939b2cece126ae1a4b1ab510212543a10f308945f06
[08:20:50] Loading Bloomfilter image
[08:23:08] Key Search in progress...
[08:23:26] BloomFilter Hit (S+) -> False Positive
[08:24:00] BloomFilter Hit (C-) -> False Positive
[08:25:06] BloomFilter Hit (S+) -> False Positive
[08:26:11] BloomFilter Hit (C+) -> Success
[08:26:11] Private key: 38726004909084954807
[08:26:11] Elapsed time: (0)hours (3)minutes (3)seconds
</pre>
