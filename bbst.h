#ifndef SBRMA2_NOC_H
#define SBRMA2_NOC_H

#include <vector>
#include "common.h"

using namespace std;

class BbST {
public:
    BbST(vector<t_value> valuesArray, vector<t_array_size> queries, t_array_size *resultLoc, int kExp);
    void solve();

    BbST(vector<t_value> valuesArray, t_array_size *resultLoc, int kExp);
    void solve(vector<t_array_size> queries);

    void verify();

    virtual ~BbST();

    size_t memUsageInBytes();

private:
    t_array_size blocksCount;
    int k, kExp, D;
    
    vector<t_value> valuesArray;
    vector<t_array_size> queries;
    t_array_size* resultLoc;

    vector<t_value> verifyVal;
    vector<t_array_size> verifyLoc;

    t_value* blocksVal2D = 0;
    t_array_size*  blocksLoc2D = 0;

    void getBlocksMinsBase();
    void getBlocksSparseTable();
    t_array_size getRangeMinLoc(const t_array_size &begIdx, const t_array_size &endIdx);
    t_array_size scanMinIdx(const t_array_size &begIdx, const t_array_size &endIdx);

    bool batchMode = false;
    void cleanup();
};


#endif //SBRMA2_NOC_H