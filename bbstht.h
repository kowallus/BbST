#ifndef BBST_BBSTHT_H
#define BBST_BBSTHT_H

#include <vector>
#include "common.h"
#include "hybtempl.h"

using namespace std;

class BbSTht {
public:

#ifdef MINI_BLOCKS
    BbSTht(const vector<t_value> &valuesArray, int kExp, int miniKExp, RMQAPI* secondaryRMQ);
#else
    BbSTht(const vector<t_value> &valuesArray, int kExp, RMQAPI* secondaryRMQ);
#endif
    void rmqBatch(const vector<t_array_size> &queries, t_array_size *resultLoc);

    t_array_size rmq(const t_array_size &begIdx, const t_array_size &endIdx);

    virtual ~BbSTht();

    size_t memUsageInBytes();

private:
    RMQAPI* secondaryRMQ;

    t_array_size blocksCount;
    int k, kExp, D, miniK, miniKExp;

    vector<t_value> verifyVal;
    vector<t_array_size> verifyLoc;

    t_value* blocksVal2D = 0;
    t_array_size*  blocksLoc2D = 0;

    t_array_size miniBlocksCount;
    int miniBlocksInBlock;
    uint8_t* miniBlocksLoc = 0;
    t_value* miniBlocksVal = 0;

    void getBlocksMinsBase(const vector<t_value> &valuesArray);
    void getBlocksSparseTable();

    inline t_array_size miniScanMinIdx(const t_array_size &begIdx, const t_array_size &endIdx);

    void cleanup();

};

#endif //BBST_BBSTHT_H
