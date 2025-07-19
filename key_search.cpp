#include <iostream>
#include <fstream>
#include <vector>
#include <thread>

#include "secp256k1/SECP256k1.h"
#include "secp256k1/Int.h"
#include "secp256k1/IntGroup.h"
#include "bloom/filter.hpp"
#include "util/util.h"

using namespace std;
using filter = boost::bloom::filter<boost::uint64_t, 32>;

const int cpuCores = std::thread::hardware_concurrency() * 4; // number of processing cores
const int POINTS_BATCH_SIZE = 1024; // Batch addition / Batch inversion (IntGroup.h class)

auto main() -> int {

    Secp256K1* secp256k1 = new Secp256K1(); secp256k1->Init(); // initialize secp256k1 context
    
    Int pk; pk.SetInt32(1); // generating power of two values (2^0..2^256) table
    uint64_t mult = 2;
    vector<Int> S_table;
    for (int i = 0; i < 256; i++)
    {
        S_table.push_back(pk);
        pk.Mult(mult);
    }
    print_time(); cout << "S_table generated" << endl;

    uint64_t range_start, range_end, block_width; // block_width = number of elements in the bloomfilter and a stride size to walk the range
    string temp, search_pub;
    ifstream inFile("settings.txt"); // load setiings from file
    getline(inFile, temp); range_start = std::stoull(temp);
    getline(inFile, temp); range_end = std::stoull(temp);
    getline(inFile, temp); block_width = std::stoull(temp);
    getline(inFile, temp); search_pub = trim(temp);
    inFile.close();
    
    print_time(); cout << "Range Start: " << range_start << " bits" << endl;
    print_time(); cout << "Range End  : " << range_end << " bits" << endl;
    print_time(); cout << "Block Width: 2^" << block_width << endl;
    print_time(); cout << "Search Pub : " << search_pub << endl;

    uint64_t stride_bits = pow(2, block_width);
    
    string bloomfile = "bloom.bf";
    filter bf;
    
    print_time(); cout << "Loading Bloomfilter image" << endl;
    
    std::ifstream in1(bloomfile, std::ios::binary);
    std::size_t c1;
    in1.read((char*) &c1, sizeof(c1));
    bf.reset(c1); // restore capacity
    boost::span<unsigned char> s1 = bf.array();
    in1.read((char*) s1.data(), s1.size()); // load array
    in1.close();
    
    auto pow10_nums = break_down_to_pow10(uint64_t(pow(2, block_width))); // decomposing the 2^block_width to the power of ten values
    vector<Point> pow10_points;                                           // to get the offset from the target point
    Int pow_key;
    for (auto& n : pow10_nums) { // calculating points corresponding to the decomposition components
        pow_key.SetInt64(n);     
        pow10_points.push_back(secp256k1->ScalarMultiplication(&pow_key));
    }
    
    auto chrono_start = std::chrono::high_resolution_clock::now();
    
    auto key_search = [&]() {

        int save_counter = 0; 
        string temp;
        Point stride_point, calc_point;
        Int stride_sum, stride;
        ifstream inFile("stride_sum.txt");
        getline(inFile, temp);
        stride_sum.SetBase10(trim(temp).data());
        inFile.close();

        stride.SetInt64(stride_bits);
        stride_point = secp256k1->ScalarMultiplication(&stride);
        
        Int two, int_Cores, range_Start, range_End, partition_Size, center_Num;
        two.SetInt32(2);
        int_Cores.SetInt32(cpuCores);
        range_Start.Set(&S_table[range_start]);
        range_End.Set(&S_table[range_end]);
        partition_Size.floor_Div(&S_table[range_start], &int_Cores);
        center_Num.floor_Div(&partition_Size, &two);

        vector<Int> range_Nums;
        for (int i = 0; i < cpuCores; i++) {
            range_Nums.push_back(range_Start);
            range_Start.Add(&partition_Size);
        }
        range_Nums.push_back(range_End);

        Int center_Int;
        vector<Int> center_Nums;
        for (int i = 0; i < cpuCores; i++) {
            center_Int.Add(&range_Nums[i], &center_Num);
            center_Nums.push_back(center_Int);
        }
        
        vector<Point> center_Points;
        for (int i = 0; i < cpuCores; i++) {
            center_Points.push_back(secp256k1->ScalarMultiplication(&center_Nums[i]));
        }
        
        vector<Point> sideway_Points;
        for (int i = 0; i < cpuCores + 1; i++) {
            sideway_Points.push_back(secp256k1->ScalarMultiplication(&range_Nums[i]));
        }

        vector<Point> pos_Points;
        vector<Point> neg_Points;
        vector<Point> pos_PointsSw;
        vector<Point> neg_PointsSw;
        if (stride_sum.IsZero()) {
            for (int i = 0; i < cpuCores; i++) {
                pos_Points.push_back(center_Points[i]);
                neg_Points.push_back(center_Points[i]);
                pos_PointsSw.push_back(sideway_Points[i]);
                neg_PointsSw.push_back(sideway_Points[i + 1]);
            }
        }
        else {
            Point offset_PosPoint = secp256k1->ScalarMultiplication(&stride_sum);
            Point offset_NegPoint = offset_PosPoint; offset_NegPoint.y.ModNeg();
            for (int i = 0; i < cpuCores; i++) {
                pos_Points.push_back(secp256k1->AddPoints(center_Points[i], offset_PosPoint));
                neg_Points.push_back(secp256k1->AddPoints(center_Points[i], offset_NegPoint));
                pos_PointsSw.push_back(secp256k1->AddPoints(sideway_Points[i], offset_PosPoint));
                neg_PointsSw.push_back(secp256k1->AddPoints(sideway_Points[i + 1], offset_NegPoint));
            }
        }

        Point addPointsPos[POINTS_BATCH_SIZE]; // array for batch addition points positive
        Point addPointsNeg[POINTS_BATCH_SIZE]; // array for batch addition points negative      
        Point batch_Add = secp256k1->DoublePoint(stride_point);
        addPointsPos[0] = stride_point;
        addPointsNeg[0] = stride_point; addPointsNeg[0].y.ModNeg();
        addPointsPos[1] = batch_Add;
        addPointsNeg[1] = batch_Add; addPointsNeg[1].y.ModNeg();
        for (int i = 2; i < POINTS_BATCH_SIZE; i++)
        {
            batch_Add = secp256k1->AddPoints(batch_Add, stride_point);
            addPointsPos[i] = batch_Add;
            addPointsNeg[i] = batch_Add;
            addPointsNeg[i].y.ModNeg();
        }
        // center_key_search
        auto center_key_search = [&](Point posStartP, Point negStartP, int threadIdx, Int center_Num, Int stride_Sum) {

            Int stride_sum; stride_sum.Set(&stride_Sum);
            Int center_num; center_num.Set(&center_Num);
            Int Int_steps, Int_temp, privkey;
            int index, count;
            uint64_t steps;
            vector<uint64_t> privkey_num;
            
            IntGroup modGroupPos(POINTS_BATCH_SIZE); // group of deltaX (x1 - x2) set for batch inversion
            IntGroup modGroupNeg(POINTS_BATCH_SIZE); // group of deltaX (x1 - x2) set for batch inversion
            Int deltaXPos[POINTS_BATCH_SIZE]; // here we store (x1 - x2) batch that will be inverted for later multiplication
            Int deltaXNeg[POINTS_BATCH_SIZE]; // here we store (x1 - x2) batch that will be inverted for later multiplication
            modGroupPos.Set(deltaXPos); // assign array deltaX to modGroup for batch inversion (JLP way set it once)
            modGroupNeg.Set(deltaXNeg); // assign array deltaX to modGroup for batch inversion (JLP way set it once)
            Int pointBatchXPos[POINTS_BATCH_SIZE]; // X coordinates of the batch
            Int pointBatchYPos[POINTS_BATCH_SIZE]; // Y coordinates of the batch
            Int pointBatchXNeg[POINTS_BATCH_SIZE]; // X coordinates of the batch
            Int pointBatchYNeg[POINTS_BATCH_SIZE]; // Y coordinates of the batch
            Int deltaYPos, deltaYNeg; // values to store the results of points addition formula
            Int slopePos[POINTS_BATCH_SIZE];
            Int slopeNeg[POINTS_BATCH_SIZE];
            
            Point startPointPos = posStartP; // start point positive
            Point startPointNeg = negStartP; // start point negative
            Point BloomP;
        
            Int batch_stride, batch_index;
            batch_stride.Mult(&stride, uint64_t(POINTS_BATCH_SIZE));

            while (true) {

                for (int i = 0; i < POINTS_BATCH_SIZE; i++) { // we compute (x1 - x2) for each entry of the entire batch
                    deltaXPos[i].ModSub(&startPointPos.x, &addPointsPos[i].x);
                    deltaXNeg[i].ModSub(&startPointNeg.x, &addPointsNeg[i].x);
                }
    
                modGroupPos.ModInv();
                modGroupNeg.ModInv();
                
                int i;
                for (i = 0; i < POINTS_BATCH_SIZE - 1; i++) { // follow points addition formula logic
                    
                    deltaYPos.ModSub(&startPointPos.y, &addPointsPos[i].y);
                    slopePos[i].ModMulK1(&deltaYPos, &deltaXPos[i]); // deltaX already inverted for each entry of the batch
                    
                    deltaYNeg.ModSub(&startPointNeg.y, &addPointsNeg[i].y);
                    slopeNeg[i].ModMulK1(&deltaYNeg, &deltaXNeg[i]); // deltaX already inverted for each entry of the batch

                    pointBatchXPos[i].ModSquareK1(&slopePos[i]); // computing just x coordinate for the (batch_size - 1)
                    pointBatchXPos[i].ModSub(&pointBatchXPos[i], &startPointPos.x);
                    pointBatchXPos[i].ModSub(&pointBatchXPos[i], &addPointsPos[i].x);
                    
                    pointBatchXNeg[i].ModSquareK1(&slopeNeg[i]); // computing just x coordinate for the (batch_size - 1)
                    pointBatchXNeg[i].ModSub(&pointBatchXNeg[i], &startPointNeg.x);
                    pointBatchXNeg[i].ModSub(&pointBatchXNeg[i], &addPointsNeg[i].x);
                    
                }
                
                deltaYPos.ModSub(&startPointPos.y, &addPointsPos[i].y);
                slopePos[i].ModMulK1(&deltaYPos, &deltaXPos[i]);
                
                deltaYNeg.ModSub(&startPointNeg.y, &addPointsNeg[i].y);
                slopeNeg[i].ModMulK1(&deltaYNeg, &deltaXNeg[i]);

                pointBatchXPos[i].ModSquareK1(&slopePos[i]);
                pointBatchXPos[i].ModSub(&pointBatchXPos[i], &startPointPos.x);
                pointBatchXPos[i].ModSub(&pointBatchXPos[i], &addPointsPos[i].x);
                
                pointBatchXNeg[i].ModSquareK1(&slopeNeg[i]);
                pointBatchXNeg[i].ModSub(&pointBatchXNeg[i], &startPointNeg.x);
                pointBatchXNeg[i].ModSub(&pointBatchXNeg[i], &addPointsNeg[i].x);
                    
                pointBatchYPos[i].ModSub(&startPointPos.x, &pointBatchXPos[i]);
                pointBatchYPos[i].ModMulK1(&slopePos[i], &pointBatchYPos[i]);
                pointBatchYPos[i].ModSub(&pointBatchYPos[i], &startPointPos.y);
                
                pointBatchYNeg[i].ModSub(&startPointNeg.x, &pointBatchXNeg[i]);
                pointBatchYNeg[i].ModMulK1(&slopeNeg[i], &pointBatchYNeg[i]);
                pointBatchYNeg[i].ModSub(&pointBatchYNeg[i], &startPointNeg.y);

                for (int i = 0; i < POINTS_BATCH_SIZE; i++) {
                    // check positive
                    if (bf.may_contain(pointBatchXPos[i].bits64[3])) {
                        
                        print_time(); cout << "BloomFilter Hit (C+) ";
                        
                        BloomP.x.Set(&pointBatchXPos[i]);
                        BloomP.y.ModSub(&startPointPos.x, &pointBatchXPos[i]);
                        BloomP.y.ModMulK1(&slopePos[i], &BloomP.y);
                        BloomP.y.ModSub(&BloomP.y, &startPointPos.y);
                        
                        privkey_num.clear();
                        index = 0;
                        for (auto& p : pow10_points) {
                            count = 0;
                            while (bf.may_contain(BloomP.x.bits64[3])) {
                                BloomP = secp256k1->SubtractPoints(BloomP, p);
                                count += 1;
                            }
                            privkey_num.push_back(pow10_nums[index] * (count - 1));
                            BloomP = secp256k1->AddPoints(BloomP, p);
                            index += 1;
                        }
                        
                        steps = 0;
                        for (auto& n : privkey_num) { steps += n; }
                        Int_steps.SetInt64(steps);
                        batch_index.Mult(&stride, uint64_t(i + 1));
                        Int_temp.Add(&stride_sum, &batch_index);
                        center_num.Add(&Int_temp); 
                        privkey.Sub(&center_num, &Int_steps);
                        calc_point = secp256k1->ScalarMultiplication(&privkey);
                        
                        if (secp256k1->GetPublicKeyHex(calc_point) == search_pub) {
                            cout << "-> Success" << endl;
                            print_time(); cout << "Private key: " << privkey.GetBase10() << endl;
                            ofstream outFile;
                            outFile.open("found.txt", ios::app);
                            outFile << privkey.GetBase10() << '\n';
                            outFile.close();
                            print_elapsed_time(chrono_start);
                            exit(0);
                        }
                        
                        cout << "-> False Positive" << endl;
                    }
                    // check negative
                    if (bf.may_contain(pointBatchXNeg[i].bits64[3])) {
                        
                        print_time(); cout << "BloomFilter Hit (C-) ";
                        
                        BloomP.x.Set(&pointBatchXNeg[i]);
                        BloomP.y.ModSub(&startPointNeg.x, &pointBatchXNeg[i]);
                        BloomP.y.ModMulK1(&slopeNeg[i], &BloomP.y);
                        BloomP.y.ModSub(&BloomP.y, &startPointNeg.y);
                        
                        privkey_num.clear();
                        index = 0;
                        for (auto& p : pow10_points) {
                            count = 0;
                            while (bf.may_contain(BloomP.x.bits64[3])) {
                                BloomP = secp256k1->SubtractPoints(BloomP, p);
                                count += 1;
                            }
                            privkey_num.push_back(pow10_nums[index] * (count - 1));
                            BloomP = secp256k1->AddPoints(BloomP, p);
                            index += 1;
                        }
                        
                        steps = 0;
                        for (auto& n : privkey_num) { steps += n; }
                        Int_steps.SetInt64(steps);
                        batch_index.Mult(&stride, uint64_t(i + 1));
                        Int_temp.Add(&stride_sum, &batch_index);
                        center_num.Sub(&Int_temp); 
                        privkey.Sub(&center_num, &Int_steps);
                        calc_point = secp256k1->ScalarMultiplication(&privkey);
                        
                        if (secp256k1->GetPublicKeyHex(calc_point) == search_pub) {
                            cout << "-> Success" << endl;
                            print_time(); cout << "Private key: " << privkey.GetBase10() << endl;
                            ofstream outFile;
                            outFile.open("found.txt", ios::app);
                            outFile << privkey.GetBase10() << '\n';
                            outFile.close();
                            print_elapsed_time(chrono_start);
                            exit(0);
                        }
                        
                        cout << "-> False Positive" << endl;
                    }
                }
                
                startPointPos.x.Set(&pointBatchXPos[POINTS_BATCH_SIZE - 1]); // last batch entry as the new startPoint for the next batch iteration
                startPointPos.y.Set(&pointBatchYPos[POINTS_BATCH_SIZE - 1]);
                
                startPointNeg.x.Set(&pointBatchXNeg[POINTS_BATCH_SIZE - 1]); // last batch entry as the new startPoint for the next batch iteration
                startPointNeg.y.Set(&pointBatchYNeg[POINTS_BATCH_SIZE - 1]);
                
                stride_sum.Add(&batch_stride);
                    
                if (threadIdx == 0) {
                    save_counter += 1;
                    if (save_counter % 60000 == 0) {
                        ofstream outFile;
                        outFile.open("stride_sum.txt");
                        outFile << stride_sum.GetBase10() << '\n';
                        outFile.close();
                        save_counter = 0;
                        print_time(); cout << "Progress written to stride_sum.txt" << endl;
                    }
                }
            } // while (true) loop end curly brace
        }; // center_key_search
        // sideway_key_search
        auto sideway_key_search = [&](Point posStartP, Point negStartP, int threadIdx, Int range_NumPos, Int range_NumNeg, Int stride_Sum) {

            Int stride_sum; stride_sum.Set(&stride_Sum);
            Int range_num_pos; range_num_pos.Set(&range_NumPos);
            Int range_num_neg; range_num_neg.Set(&range_NumNeg);
            Int Int_steps, Int_temp, privkey;
            int index, count;
            uint64_t steps;
            vector<uint64_t> privkey_num;
            
            IntGroup modGroupPos(POINTS_BATCH_SIZE); // group of deltaX (x1 - x2) set for batch inversion
            IntGroup modGroupNeg(POINTS_BATCH_SIZE); // group of deltaX (x1 - x2) set for batch inversion
            Int deltaXPos[POINTS_BATCH_SIZE]; // here we store (x1 - x2) batch that will be inverted for later multiplication
            Int deltaXNeg[POINTS_BATCH_SIZE]; // here we store (x1 - x2) batch that will be inverted for later multiplication
            modGroupPos.Set(deltaXPos); // assign array deltaX to modGroup for batch inversion (JLP way set it once)
            modGroupNeg.Set(deltaXNeg); // assign array deltaX to modGroup for batch inversion (JLP way set it once)
            Int pointBatchXPos[POINTS_BATCH_SIZE]; // X coordinates of the batch
            Int pointBatchYPos[POINTS_BATCH_SIZE]; // Y coordinates of the batch
            Int pointBatchXNeg[POINTS_BATCH_SIZE]; // X coordinates of the batch
            Int pointBatchYNeg[POINTS_BATCH_SIZE]; // Y coordinates of the batch
            Int deltaYPos, deltaYNeg; // values to store the results of points addition formula
            Int slopePos[POINTS_BATCH_SIZE];
            Int slopeNeg[POINTS_BATCH_SIZE];
            
            Point startPointPos = posStartP; // start point positive
            Point startPointNeg = negStartP; // start point negative
            Point BloomP;
        
            Int batch_stride, batch_index;
            batch_stride.Mult(&stride, uint64_t(POINTS_BATCH_SIZE));

            while (true) {

                for (int i = 0; i < POINTS_BATCH_SIZE; i++) { // we compute (x1 - x2) for each entry of the entire batch
                    deltaXPos[i].ModSub(&startPointPos.x, &addPointsPos[i].x);
                    deltaXNeg[i].ModSub(&startPointNeg.x, &addPointsNeg[i].x);
                }
    
                modGroupPos.ModInv();
                modGroupNeg.ModInv();
                
                int i;
                for (i = 0; i < POINTS_BATCH_SIZE - 1; i++) { // follow points addition formula logic
                    
                    deltaYPos.ModSub(&startPointPos.y, &addPointsPos[i].y);
                    slopePos[i].ModMulK1(&deltaYPos, &deltaXPos[i]); // deltaX already inverted for each entry of the batch
                    
                    deltaYNeg.ModSub(&startPointNeg.y, &addPointsNeg[i].y);
                    slopeNeg[i].ModMulK1(&deltaYNeg, &deltaXNeg[i]); // deltaX already inverted for each entry of the batch

                    pointBatchXPos[i].ModSquareK1(&slopePos[i]); // computing just x coordinate for the (batch_size - 1)
                    pointBatchXPos[i].ModSub(&pointBatchXPos[i], &startPointPos.x);
                    pointBatchXPos[i].ModSub(&pointBatchXPos[i], &addPointsPos[i].x);
                    
                    pointBatchXNeg[i].ModSquareK1(&slopeNeg[i]); // computing just x coordinate for the (batch_size - 1)
                    pointBatchXNeg[i].ModSub(&pointBatchXNeg[i], &startPointNeg.x);
                    pointBatchXNeg[i].ModSub(&pointBatchXNeg[i], &addPointsNeg[i].x);
                    
                }
                
                deltaYPos.ModSub(&startPointPos.y, &addPointsPos[i].y);
                slopePos[i].ModMulK1(&deltaYPos, &deltaXPos[i]);
                
                deltaYNeg.ModSub(&startPointNeg.y, &addPointsNeg[i].y);
                slopeNeg[i].ModMulK1(&deltaYNeg, &deltaXNeg[i]);

                pointBatchXPos[i].ModSquareK1(&slopePos[i]);
                pointBatchXPos[i].ModSub(&pointBatchXPos[i], &startPointPos.x);
                pointBatchXPos[i].ModSub(&pointBatchXPos[i], &addPointsPos[i].x);
                
                pointBatchXNeg[i].ModSquareK1(&slopeNeg[i]);
                pointBatchXNeg[i].ModSub(&pointBatchXNeg[i], &startPointNeg.x);
                pointBatchXNeg[i].ModSub(&pointBatchXNeg[i], &addPointsNeg[i].x);
                    
                pointBatchYPos[i].ModSub(&startPointPos.x, &pointBatchXPos[i]);
                pointBatchYPos[i].ModMulK1(&slopePos[i], &pointBatchYPos[i]);
                pointBatchYPos[i].ModSub(&pointBatchYPos[i], &startPointPos.y);
                
                pointBatchYNeg[i].ModSub(&startPointNeg.x, &pointBatchXNeg[i]);
                pointBatchYNeg[i].ModMulK1(&slopeNeg[i], &pointBatchYNeg[i]);
                pointBatchYNeg[i].ModSub(&pointBatchYNeg[i], &startPointNeg.y);

                for (int i = 0; i < POINTS_BATCH_SIZE; i++) {
                    // check positive
                    if (bf.may_contain(pointBatchXPos[i].bits64[3])) {
                        
                        print_time(); cout << "BloomFilter Hit (S+) ";
                        
                        BloomP.x.Set(&pointBatchXPos[i]);
                        BloomP.y.ModSub(&startPointPos.x, &pointBatchXPos[i]);
                        BloomP.y.ModMulK1(&slopePos[i], &BloomP.y);
                        BloomP.y.ModSub(&BloomP.y, &startPointPos.y);
                        
                        privkey_num.clear();
                        index = 0;
                        for (auto& p : pow10_points) {
                            count = 0;
                            while (bf.may_contain(BloomP.x.bits64[3])) {
                                BloomP = secp256k1->SubtractPoints(BloomP, p);
                                count += 1;
                            }
                            privkey_num.push_back(pow10_nums[index] * (count - 1));
                            BloomP = secp256k1->AddPoints(BloomP, p);
                            index += 1;
                        }
                        
                        steps = 0;
                        for (auto& n : privkey_num) { steps += n; }
                        Int_steps.SetInt64(steps);
                        batch_index.Mult(&stride, uint64_t(i + 1));
                        Int_temp.Add(&stride_sum, &batch_index);
                        range_num_pos.Add(&Int_temp); 
                        privkey.Sub(&range_num_pos, &Int_steps);
                        calc_point = secp256k1->ScalarMultiplication(&privkey);
                        
                        if (secp256k1->GetPublicKeyHex(calc_point) == search_pub) {
                            cout << "-> Success" << endl;
                            print_time(); cout << "Private key: " << privkey.GetBase10() << endl;
                            ofstream outFile;
                            outFile.open("found.txt", ios::app);
                            outFile << privkey.GetBase10() << '\n';
                            outFile.close();
                            print_elapsed_time(chrono_start);
                            exit(0);
                        }
                        
                        cout << "-> False Positive" << endl;
                    }
                    // check negative
                    if (bf.may_contain(pointBatchXNeg[i].bits64[3])) {
                        
                        print_time(); cout << "BloomFilter Hit (S-) ";
                        
                        BloomP.x.Set(&pointBatchXNeg[i]);
                        BloomP.y.ModSub(&startPointNeg.x, &pointBatchXNeg[i]);
                        BloomP.y.ModMulK1(&slopeNeg[i], &BloomP.y);
                        BloomP.y.ModSub(&BloomP.y, &startPointNeg.y);
                        
                        privkey_num.clear();
                        index = 0;
                        for (auto& p : pow10_points) {
                            count = 0;
                            while (bf.may_contain(BloomP.x.bits64[3])) {
                                BloomP = secp256k1->SubtractPoints(BloomP, p);
                                count += 1;
                            }
                            privkey_num.push_back(pow10_nums[index] * (count - 1));
                            BloomP = secp256k1->AddPoints(BloomP, p);
                            index += 1;
                        }
                        
                        steps = 0;
                        for (auto& n : privkey_num) { steps += n; }
                        Int_steps.SetInt64(steps);
                        batch_index.Mult(&stride, uint64_t(i + 1));
                        Int_temp.Add(&stride_sum, &batch_index);
                        range_num_neg.Sub(&Int_temp); 
                        privkey.Sub(&range_num_neg, &Int_steps);
                        calc_point = secp256k1->ScalarMultiplication(&privkey);
                        
                        if (secp256k1->GetPublicKeyHex(calc_point) == search_pub) {
                            cout << "-> Success" << endl;
                            print_time(); cout << "Private key: " << privkey.GetBase10() << endl;
                            ofstream outFile;
                            outFile.open("found.txt", ios::app);
                            outFile << privkey.GetBase10() << '\n';
                            outFile.close();
                            print_elapsed_time(chrono_start);
                            exit(0);
                        }
                        
                        cout << "-> False Positive" << endl;
                    }
                }
                
                startPointPos.x.Set(&pointBatchXPos[POINTS_BATCH_SIZE - 1]); // last batch entry as the new startPoint for the next batch iteration
                startPointPos.y.Set(&pointBatchYPos[POINTS_BATCH_SIZE - 1]);
                
                startPointNeg.x.Set(&pointBatchXNeg[POINTS_BATCH_SIZE - 1]); // last batch entry as the new startPoint for the next batch iteration
                startPointNeg.y.Set(&pointBatchYNeg[POINTS_BATCH_SIZE - 1]);
                
                stride_sum.Add(&batch_stride);

            } // while (true) loop end curly brace
        }; // sideway_key_search
        
        std::thread keySearch_CThreads[cpuCores];
        std::thread keySearch_SThreads[cpuCores];
        for (int i = 0; i < cpuCores; i++) {
            keySearch_CThreads[i] = std::thread(center_key_search, pos_Points[i], neg_Points[i], i, center_Nums[i], stride_sum);
            keySearch_SThreads[i] = std::thread(sideway_key_search, pos_PointsSw[i], neg_PointsSw[i], i, range_Nums[i], range_Nums[i + 1], stride_sum);
        }

        for (int i = 0; i < cpuCores; i++) {
            keySearch_CThreads[i].join();
            keySearch_SThreads[i].join();
        }
    };
    
    print_time(); cout << "Key Search in progress..." << endl;
    
    std::thread thread(key_search);
    
    thread.join();
}
