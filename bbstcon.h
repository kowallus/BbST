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
    BbSTcon(sortingAlg_enum sortingAlg, int kExp);

    void rmqBatch(const t_value* valuesArray, const t_array_size n, const vector<t_array_size> &queries, t_array_size *resultLoc);

    size_t memUsageInBytes();

private:
    t_array_size blocksCount;
    int k, kExp, D;
    sortingAlg_enum sortingAlg;

    const t_value *valuesArray;
    t_array_size n;
    size_t q;

    vector<t_value> verifyVal;
    vector<t_array_size> verifyLoc;

    t_array_size_2x* bounds = 0;
    t_array_size* queries2ContractedIdx = 0;
    t_value* contractedVal = 0;
    t_array_size* contractedLoc = 0;
    t_value* blocksVal2D = 0;
    t_array_size*  blocksLoc2D = 0;

    void getUniqueBoundsSorted(const vector<t_array_size> &queries);
    void getContractedMins();
    void getBlocksMins();
    t_array_size getContractedRMQ(const t_array_size &begContIdx, const t_array_size &endContIdx);
    t_array_size scanContractedMinIdx(const t_array_size &begContIdx, const t_array_size &endContIdx);

    void cleanup();
};


#endif //BBSTCON_H
