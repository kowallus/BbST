#ifndef BBST_BBSTQHT_H
#define BBST_BBSTQHT_H

#include <vector>
#include "common.h"
#include "hybtempl.h"

using namespace std;

template<typename t_qvalue, int max_qvalue>
class BbSTqht {
public:

#ifdef MINI_BLOCKS
    BbSTqht(const vector<t_value> &valuesArray, int kExp, int miniKExp, RMQAPI* secondaryRMQ);
#else
    BbSTqht(const vector<t_value> &valuesArray, int kExp, RMQAPI* secondaryRMQ);
#endif

    void rmqBatch(const vector<t_array_size> &queries, t_array_size *resultLoc);

    t_array_size rmq(const t_array_size &begIdx, const t_array_size &endIdx);

    virtual ~BbSTqht();

    size_t memUsageInBytes();

private:
    RMQAPI* secondaryRMQ;

    t_array_size blocksCount;
    int k, kExp, D, miniK, miniKExp;

    t_qvalue* blocksQVal2D = 0;
    t_array_size*  blocksLoc2D = 0;

    t_array_size miniBlocksCount;
    int miniBlocksInBlock;
    uint8_t* miniBlocksLoc = 0;
    t_qvalue* miniBlocksQVal = 0;

    void prepareMinTables(const vector<t_value> &valuesArray);
    void prepareBlocksSparseTable(vector<t_value> &blocksVal2D, const t_value &minMinVal, const t_value &maxMinVal);

    inline t_array_size miniScanMinIdx(const t_array_size &begIdx, const t_array_size &endIdx);

    void cleanup();

};

#include "bbstqht.hpp"

#endif //BBST_BBSTQHT_H
