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

puzzle = 1103873984953507439627945351144005829577
start_bits = 129
end_bits   = 130  
range_start = S_table[start_bits]
range_end   = S_table[end_bits]

stride = range_start // num_threads
range_nums = []
range_nums.append(range_start)
walker = range_start
for i in range(num_threads):
    walker += stride
    range_nums.append(walker)

print(f'Threads     : {num_threads}')
print(f'Range start : {range_start}')
print(f'Range end   : {range_end}')
print(f'Block width : {stride}')
start_num = 0
end_num   = 0
for idx, num in enumerate(range_nums):
    if puzzle > num:
        start_num = range_nums[idx]
        end_num   = range_nums[idx + 1]
print()
print(f'Puzzle -> {puzzle}')
print(f'Found in block [{start_num} - {end_num}]')
print(f'Distance from start : {puzzle - start_num}')
print(f'Distance from end   : {end_num - puzzle}')
get_range(puzzle - start_num)
get_range(end_num - puzzle)

