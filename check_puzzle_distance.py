import multiprocessing as mp
import math

def get_range(number):
    start = end = 0
    for i in range(256):
        if S_table[i] > number:
            end = int(math.log2(S_table[i]))
            start = end - 1
            break
    print(f'Number ({number}) is within 2^{start}...2^{end} bits')

num_threads = mp.cpu_count() * 4

S_table = []
pk = 1
for k in range(256): 
    S_table.append(pk) 
    pk *= 2

batch_size = 1024
puzzle = 199976667976342049
start_bits = 57
end_bits   = 58  
range_start = S_table[start_bits]
range_end   = S_table[end_bits]
stride_bits = 28
stride      = S_table[stride_bits] # 2^stride_bits

block_width = range_start // num_threads
range_nums = []
range_nums.append(range_start)
walker = range_start
for i in range(num_threads):
    walker += block_width
    range_nums.append(walker)

print(f'Threads     : {num_threads}')
print(f'Range start : [{start_bits} bits] {range_start}')
print(f'Range end   : [{end_bits} bits] {range_end}')
print(f'Stride size : [{stride_bits} bits] {stride}')
print(f'Block width : {block_width}')
get_range(block_width)
start_num = 0
end_num   = 0
index_blk = 0
for idx, num in enumerate(range_nums):
    if puzzle > num:
        index_blk = idx
        start_num = range_nums[idx]
        end_num   = range_nums[idx + 1]
print()
print(f'Puzzle -> {puzzle}')
print(f'Found in block #{index_blk + 1} of {num_threads} [{start_num} - {end_num}]')
dist_s = puzzle - start_num
dist_e = end_num - puzzle
mid = start_num + (block_width // 2)
print(f'Block middle: {mid}')
if puzzle < mid:
   print(f'Puzzle is in the lower block half\n')
else:
   print(f'Puzzle is in the higher block half\n')
   
print(f'Distance from block start : {dist_s}')
get_range(dist_s)
print(f'Takes : {dist_s // stride} strides of size 2^{stride_bits} to solve')
print(f'Takes : {(dist_s // stride) // batch_size} batches of size {batch_size} to solve')
print()
print(f'Distance from block end   : {dist_e}')
get_range(dist_e)
print(f'Takes : {dist_e // stride} strides of size 2^{stride_bits} to solve')
print(f'Takes : {(dist_e // stride) // batch_size} batches of size {batch_size} to solve')
print()
dist_m = 0
if puzzle < mid:
    dist_m = mid - puzzle
    print(f'Distance from middle  [-] : {dist_m}')
    get_range(dist_m)
    print(f'Takes : {dist_m // stride} strides of size 2^{stride_bits} to solve')
    print(f'Takes : {(dist_m // stride) // batch_size} batches of size {batch_size} to solve')
else:
    dist_m = puzzle - mid
    print(f'Distance from middle  [+] : {dist_m}')
    get_range(dist_m)
    print(f'Takes : {dist_m // stride} strides of size 2^{stride_bits} to solve')
    print(f'Takes : {(dist_m // stride) // batch_size} batches of size {batch_size} to solve')
print()
min_val = min(dist_s, dist_e, dist_m)
if min_val == dist_s:
    print(f'Distance from block start is the shortest')
elif min_val == dist_e:
    print(f'Distance from block end is the shortest')
else:
    print(f'Distance from middle of the block is the shortest')

'''
[06:01:34] BloomFilter Hit (C+) -> Success
[06:01:34] Private key: 970436974005023690481
'''
