#include "BinaryFuseFilter/sharded_filter.hpp"
#include "BinaryFuseFilter/binaryfusefilter.h"
#include "BinaryFuseFilter/filter.hpp"
#include "BinaryFuseFilter/page.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <thread>
#include <algorithm>
#include <span>
#include <atomic>


#include "secp256k1/SECP256k1.h"
#include "secp256k1/Int.h"
#include "secp256k1/IntGroup.h"
#include "util/util.h"

using namespace std;

const int POINTS_BATCH_SIZE = 1024; // Batch addition with batch inversion using IntGroup class

void unique_sort(std::vector<uint64_t>& vec) {
    std::sort(vec.begin(), vec.end());
    auto last = std::unique(vec.begin(), vec.end());
    vec.erase(last, vec.end());
}

int ready_flag = 1;
int vec_clear_flag = 1;
const int shard_bits = 4;
const int shards_count = 1 << shard_bits;

void process_shard(std::vector<uint64_t>& shard_nums, binfuse::sharded_filter16_sink &sink, int loop_id) {

    ::ready_flag = 0;

    binfuse::filter16 shard_fil(shard_nums);

    ::vec_clear_flag = 0;

    shard_nums.clear();

    ::vec_clear_flag = 1;
            
    print_time(); cout << "Shard " << (loop_id + 1) << '/' << shards_count << '\n';

    sink.add_shard(shard_fil, loop_id);

    ::ready_flag = 1;
}

auto main() -> int {
    
    auto chrono_start = std::chrono::high_resolution_clock::now();    // starting the timer
    Secp256K1* secp256k1 = new Secp256K1(); secp256k1->Init(); // initializing secp256k1 context

    std::filesystem::remove("sharded_filter.bin");

    uint64_t range_start, range_end, block_width; // block_width = number of elements in the bloomfilter and a stride size to walk the range
    string temp, search_pub;
    ifstream inFile("settings.txt");
    getline(inFile, temp); range_start = std::stoull(temp);
    getline(inFile, temp); range_end = std::stoull(temp);
    getline(inFile, temp); block_width = std::stoull(temp);
    getline(inFile, temp); search_pub = trim(temp);
    inFile.close();
    
    print_time(); cout << "Range Start: " << range_start << " bits" << endl;
    print_time(); cout << "Range End  : " << range_end << " bits" << endl;
    print_time(); cout << "Block Width: 2^" << block_width << endl;
    print_time(); cout << "Search Pub : " << search_pub << endl;

    ofstream outFile;
    outFile.open("stride_sum.txt", ios::out);
    outFile << "0" << '\n';
    outFile.close();

    print_time(); cout << "Stride_sum written to file" << endl;

    Point target_point = secp256k1->ParsePublicKeyHex(search_pub);
 
    uint64_t n_elements = pow(2, block_width);  // number of elements == 2^block_width
    uint64_t keysPerThread = n_elements / shards_count; // elements per thread
    Int add_key; add_key.SetInt64(keysPerThread);
    Point Add_Point = secp256k1->ScalarMultiplication(&add_key); // helper point to calculate the starting points
    
    Point addPoints[POINTS_BATCH_SIZE]; // array for the batch addition points(1G .. 1024G)
    Point batch_Add = secp256k1->DoublePoint(secp256k1->G); // 2G
    addPoints[0] = secp256k1->G; // 1G
    addPoints[1] = batch_Add;    // 2G
    for (int i = 2; i < POINTS_BATCH_SIZE; i++) // filling in the batch addition points array with points from(3G .. 1024G)
    {
        batch_Add = secp256k1->AddPoints(batch_Add, secp256k1->G);
        addPoints[i] = batch_Add;
    }
    
    int nbBatch = keysPerThread / POINTS_BATCH_SIZE; // number of batches for the single thread

    binfuse::sharded_filter16_sink sink("sharded_filter.bin", shard_bits);

    vector<uint64_t> shard_nums;
    
    auto binary_fuse_filter = [&]() {

        for (int loop_id = 0; loop_id < shards_count; loop_id++) {

            Point P = secp256k1->SubtractPoints(target_point, secp256k1->G);
            Point starting_points[shards_count];
            for (int i = 0; i < shards_count; i++) { // calculating the starting points 
                starting_points[i] = P;
                P = secp256k1->AddPoints(P, Add_Point);
            }

            vector<uint64_t> vec[shards_count];

            auto process_chunk = [&](Point start_point, int threadId) { // function for a thread
                
                IntGroup modGroup(POINTS_BATCH_SIZE); // group of deltaX (x1 - x2) set for batch inversion
                Int deltaX[POINTS_BATCH_SIZE]; // here we store (x1 - x2) batch that will be inverted for later multiplication
                modGroup.Set(deltaX); // assign array deltaX to modGroup for batch inversion (JLP one time Set)
                Int pointBatchX[POINTS_BATCH_SIZE]; // X coordinates of the batch
                Int pointBatchY[POINTS_BATCH_SIZE]; // Y coordinates of the batch
                Int deltaY, slope; // values to store the results of points addition formula
                          
                Point startPoint = start_point; // start point
                
                for (int s = 0; s < nbBatch; s++) {
                    
                    for (int i = 0; i < POINTS_BATCH_SIZE; i++) { // we compute (x1 - x2) for each entry of the entire batch
                        deltaX[i].ModSub(&startPoint.x, &addPoints[i].x); // insert each entry into the deltaX array
                    }
        
                    modGroup.ModInv();    // doing batch inversion
                    
                    int i;
                    for (i = 0; i < POINTS_BATCH_SIZE - 1; i++) { // follow points addition formula logic
                        
                        deltaY.ModSub(&startPoint.y, &addPoints[i].y);
                        slope.ModMulK1(&deltaY, &deltaX[i]); // deltaX already inverted for each entry of the batch

                        pointBatchX[i].ModSquareK1(&slope); // calculating just X coordinate
                        pointBatchX[i].ModSub(&pointBatchX[i], &startPoint.x);
                        pointBatchX[i].ModSub(&pointBatchX[i], &addPoints[i].x);

                    }
                    
                    deltaY.ModSub(&startPoint.y, &addPoints[i].y); // calculating X,Y coordinates for the last of the batch entry (used also as the next startPoint)
                    slope.ModMulK1(&deltaY, &deltaX[i]);

                    pointBatchX[i].ModSquareK1(&slope);
                    pointBatchX[i].ModSub(&pointBatchX[i], &startPoint.x);
                    pointBatchX[i].ModSub(&pointBatchX[i], &addPoints[i].x);
                        
                    pointBatchY[i].ModSub(&startPoint.x, &pointBatchX[i]);
                    pointBatchY[i].ModMulK1(&slope, &pointBatchY[i]);
                    pointBatchY[i].ModSub(&pointBatchY[i], &startPoint.y);
                    
                    for (int i = 0; i < POINTS_BATCH_SIZE; i++) { // inserting all batch points X coordinates into the bloomfilter
                        if (((pointBatchX[i].bits64[3] >> 60) & 0xF) == loop_id) vec[threadId].push_back(pointBatchX[i].bits64[3]);
                    }
                    
                    startPoint.x.Set(&pointBatchX[POINTS_BATCH_SIZE - 1]); // setting the new startPoint for the next batch iteration
                    startPoint.y.Set(&pointBatchY[POINTS_BATCH_SIZE - 1]);

                }
                
            };
            
            std::thread myThreads[shards_count]; // launching threads
            for (int i = 0; i < shards_count; i++) {
                myThreads[i] = std::thread(process_chunk, starting_points[i], i);
            }

            for (int i = 0; i < shards_count; i++) {
                myThreads[i].join(); // waiting for threads to finish
            }

            do {} while(!vec_clear_flag);

            for (int i = 0; i < shards_count; i++) {
                shard_nums.insert(shard_nums.end(), vec[i].begin(), vec[i].end());
            }

            do {} while(!ready_flag);
            
            std::thread shard_thread(process_shard, std::ref(shard_nums), std::ref(sink), loop_id);
            if((shards_count - loop_id) != 1) {
                shard_thread.detach();
            }
            else {
                shard_thread.join();
            }

        }       
    };

    std::thread binary_fuse_filter_thread(binary_fuse_filter);
    
    print_time(); cout << "Creating Binary Fuse Filter shards" << '\n';
    
    binary_fuse_filter_thread.join();
    
    print_elapsed_time(chrono_start);

}
