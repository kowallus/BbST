#ifndef SBRMA2_NOC_H
#define SBRMA2_NOC_H

#include <vector>
#include "common.h"

using namespace std;

class BbST {
public:
    BbST(vector<t_value> valuesArray, vector<t_array_size> queries, t_array_size *resultLoc, int k);

    void solve();

    void verify();

    size_t memUsageInBytes();

private:
    t_array_size blocksCount;
    int k, D;
    
    vector<t_value> valuesArray;
    vector<t_array_size> queries;
    t_array_size* resultLoc;

    vector<t_value> verifyVal;
    vector<t_array_size> verifyLoc;

    t_value* blocksVal2D;
    t_array_size*  blocksLoc2D;

    void getBlocksMinsBase();
    void getBlocksSparseTable();
    t_array_size getRangeMinLoc(const t_array_size &begIdx, const t_array_size &endIdx);
    t_array_size scanMinIdx(const t_array_size &begIdx, const t_array_size &endIdx);

    void cleanup();
};


#endif //SBRMA2_NOC_H