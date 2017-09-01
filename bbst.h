#ifndef SBRMA2_NOC_H
#define SBRMA2_NOC_H

#include <vector>
#include "common.h"

using namespace std;

class BbST {
public:
    void solve();

#ifdef MINI_BLOCKS
    BbST(vector<t_value> valuesArray, t_array_size *resultLoc, int kExp, int miniKExp);
    BbST(vector<t_value> valuesArray, vector<t_array_size> queries, t_array_size *resultLoc, int kExp, int miniKExp);
#else 
    BbST(vector<t_value> valuesArray, t_array_size *resultLoc, int kExp);
    BbST(vector<t_value> valuesArray, vector<t_array_size> queries, t_array_size *resultLoc, int kExp);
#endif
    void prepare();
    void solve(vector<t_array_size> queries);

    void verify();

    virtual ~BbST();

    size_t memUsageInBytes();

private:
    t_array_size blocksCount;
    int k, kExp, D, miniK, miniKExp;

    vector<t_value> valuesArray;
    vector<t_array_size> queries;
    t_array_size* resultLoc;

    vector<t_value> verifyVal;
    vector<t_array_size> verifyLoc;

    t_value* blocksVal2D = 0;
    t_array_size*  blocksLoc2D = 0;

    t_array_size miniBlocksCount;
    int miniBlocksInBlock;
    uint8_t* miniBlocksLoc = 0;

    void getBlocksMinsBase();
    void getBlocksSparseTable();
    t_array_size getRangeMinLoc(const t_array_size &begIdx, const t_array_size &endIdx);
    t_array_size scanMinIdx(const t_array_size &begIdx, const t_array_size &endIdx);
    inline t_array_size rawScanMinIdx(const t_array_size &begIdx, const t_array_size &endIdx);
    inline t_array_size miniScanMinIdx(const t_array_size &begIdx, const t_array_size &endIdx);

    bool batchMode = false;
    void cleanup();

};


#endif //SBRMA2_NOC_H