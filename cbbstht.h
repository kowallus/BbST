#ifndef BBST_BBSTQHT_H
#define BBST_BBSTQHT_H

#include <vector>
#include "common.h"
#include "hybtempl.h"

using namespace std;

template<typename t_qvalue, int max_qvalue>
class CBbSTht {
public:

#ifdef MINI_BLOCKS
    CBbSTht(const vector<t_value> &valuesArray, int kExp, int miniKExp, RMQAPI* secondaryRMQ);
#else
    CBbSTht(const vector<t_value> &valuesArray, int kExp, RMQAPI* secondaryRMQ);
#endif

    void rmqBatch(const vector<t_array_size> &queries, t_array_size *resultLoc);

    t_array_size rmq(const t_array_size &begIdx, const t_array_size &endIdx);

    virtual ~CBbSTht();

    size_t memUsageInBytes();

private:
    RMQAPI* secondaryRMQ;

    t_array_size blocksCount;
    int k, kExp, BD, D, miniK, miniKExp;

    uint8_t* baseBlocksValLoc2D = 0; // full ST info (location + value) for base layers (0, 9, 18, etc.)
    uint8_t* blocksRelativeLoc2D = 0; // location of minima in a closest lower base layer; for layers (1--8, 10--17, etc.)

    t_value minMinVal, maxMinVal;
    t_array_size miniBlocksCount;
    int miniBlocksInBlock;
    uint8_t* miniBlocksLoc = 0;
    t_qvalue* miniBlocksQVal = 0;

    void prepareMinTables(const vector<t_value> &valuesArray);
    void prepareBlocksSparseTable(vector<t_value> &tempBlocksVal, vector<t_array_size> &tempBlocksLoc);
    inline t_qvalue quantizeValue(const t_value value);

    inline t_array_size miniScanMinIdx(const t_array_size &begIdx, const t_array_size &endIdx);

    void cleanup();

};

#include "cbbstht.hpp"

#endif //BBST_BBSTQHT_H
