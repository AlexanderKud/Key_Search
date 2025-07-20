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
[06:29:44] Range Start: 64 bits
[06:29:44] Range End  : 65 bits
[06:29:44] Block Width: 2^32
[06:29:44] Search Pub : 026d9eb6591ef95d0b26235b981f024d9019ab51b1fa10543910953cb6be201fac
[06:29:44] Stride_sum written to file
[06:29:44] Creating bloomfilter image
[07:15:48] Writing image to bloom.bf
[07:20:03] Elapsed time: (0)hours (50)minutes (18)seconds
                                                                                                            
┌──(alexander㉿DESKTOP-ZMYIN9A))-[~/Documents/Key_Search]
└─$ ./key_search
[07:20:39] S_table generated
[07:20:39] Range Start: 64 bits
[07:20:39] Range End  : 65 bits
[07:20:39] Block Width: 2^32
[07:20:39] Search Pub : 026d9eb6591ef95d0b26235b981f024d9019ab51b1fa10543910953cb6be201fac
[07:20:39] Loading Bloomfilter image
[07:22:49] Key Search in progress...
[07:24:10] BloomFilter Hit (C+) -> Success
[07:24:10] Private key: 33764838787567902188
[07:24:10] Elapsed time: (0)hours (1)minutes (20)seconds
</pre>
