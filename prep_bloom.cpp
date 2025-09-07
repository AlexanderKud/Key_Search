#include <iostream>
#include <fstream>
#include <vector>
#include <thread>
#include <omp.h>

#include "secp256k1/SECP256k1.h"
#include "secp256k1/Int.h"
#include "secp256k1/IntGroup.h"
#include "bloom/filter.hpp"
#include "util/util.h"

using namespace std;
using filter = boost::bloom::filter<boost::uint64_t, 32>; // bloomfilter settings

const double error = 0.0000000001; // errror rate for bloomfilter
const int n_cores = 2; //number of processing cores
const int POINTS_BATCH_SIZE = 1024; // Batch addition with batch inversion using IntGroup class

auto main() -> int {
    
    auto chrono_start = std::chrono::high_resolution_clock::now();    // starting the timer
    Secp256K1* secp256k1 = new Secp256K1(); secp256k1->Init(); // initializing secp256k1 context
    
    omp_lock_t lock;
    omp_init_lock(&lock);
   
    std::remove("stride_sum.txt"); // remove previous settings and bloom files
    std::remove("bloom.bf");

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
    uint64_t keysPerThread = n_elements / n_cores; // elements per thread
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
    
    auto bloom_create = [&]() {
        
        string bloomfile = "bloom.bf";
        Point P = secp256k1->SubtractPoints(target_point, secp256k1->G);
        Point starting_points[n_cores];
        for (int i = 0; i < n_cores; i++) { // calculating the starting points 
            starting_points[i] = P;
            P = secp256k1->AddPoints(P, Add_Point);
        }

        filter bf(n_elements, error);
        
        auto process_chunk = [&](Point start_point) { // function for a thread
            
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
                
                omp_set_lock(&lock);
                for (int i = 0; i < POINTS_BATCH_SIZE; i++) { // inserting all batch points X coordinates into the bloomfilter
                    bf.insert(pointBatchX[i].bits64[3]);
                }
                omp_unset_lock(&lock);
                
                startPoint.x.Set(&pointBatchX[POINTS_BATCH_SIZE - 1]); // setting the new startPoint for the next batch iteration
                startPoint.y.Set(&pointBatchY[POINTS_BATCH_SIZE - 1]);

            }
            
        };
        
        std::thread myThreads[n_cores]; // launching threads
        for (int i = 0; i < n_cores; i++) {
            myThreads[i] = std::thread(process_chunk, starting_points[i]);
        }

        for (int i = 0; i < n_cores; i++) {
            myThreads[i].join(); // waiting for threads to finish
        }
        
        omp_destroy_lock(&lock);
        
        print_time(); cout << "Writing image to bloom.bf" << '\n';
        std::ofstream out(bloomfile, std::ios::binary);
        std::size_t c1= bf.capacity();
        out.write((const char*) &c1, sizeof(c1)); // save capacity (bits)
        boost::span<const unsigned char> s1 = bf.array();
        out.write((const char*) s1.data(), s1.size()); // save array
        out.close();
    };


    std::thread bloom_thread(bloom_create);
    
    print_time(); cout << "Creating bloomfilter image" << '\n';
    
    bloom_thread.join();
    
    print_elapsed_time(chrono_start);

}
