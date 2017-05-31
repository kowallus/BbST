#ifndef BBSTCON_H
#define BBSTCON_H

#include <vector>
#include "common.h"

using namespace std;

enum sortingAlg_enum
{
    csort = 'q',
    stdsort = 's',
    kxradixsort = 'r',
    ompparallelsort = 'p',
    pssparallelsort = 'i'
};

class BbSTcon {
public:
    BbSTcon(std::vector<t_value , std::allocator<t_value>> &valuesArray,
              std::vector<t_array_size, std::allocator<t_array_size>> &queries,
              t_array_size* resultLoc,
              sortingAlg_enum sortingAlg, int noOfThreads);

    void solve();

    void verify();

    size_t memUsageInBytes();

private:
    int D = 32;
    //const int K = 512;
    sortingAlg_enum sortingAlg;

    vector<t_value> valuesArray;
    vector<t_array_size> queries;
    t_array_size* resultLoc;

    vector<t_value> verifyVal;
    vector<t_array_size> verifyLoc;

    t_array_size_2x* bounds = 0;
    t_array_size* queries2ContractedIdx = 0;
    t_value* contractedVal = 0;
    t_array_size* contractedLoc = 0;
    t_value* blocksVal2D = 0;
    t_array_size*  blocksLoc2D = 0;

    void getUniqueBoundsSorted();
    void getContractedMins();
    void getBlocksMins();
    t_array_size getRangeMinLoc(const t_array_size &begContIdx, const t_array_size &endContIdx);
    t_array_size scanContractedMinIdx(const t_array_size &begContIdx, const t_array_size &endContIdx);

    void cleanup();
};


#endif //BBSTCON_H
