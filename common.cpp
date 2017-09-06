#include "common.h"

void verify(const vector<t_value> &valuesArray, const vector<t_array_size> &queries, t_array_size *resultLoc) {
    cout << "Solution verification..." << std::endl;
    unsigned int i;
    const int q = queries.size() / 2;
    vector<t_value> verifyVal(q);
    vector<t_array_size> verifyLoc(q);
    const t_value *const vaPtr = &valuesArray[0];
    const t_array_size *const qPtr = &queries[0];
    for(i = 0; i < q; i++) {
        const t_value *const begPtr = vaPtr + qPtr[2*i];
        const t_value *const endPtr = vaPtr + qPtr[2*i + 1] + 1;
        const t_value *const minValPtr = std::min_element(begPtr, endPtr);
        verifyVal[i] = *minValPtr;
        verifyLoc[i] = minValPtr - vaPtr;
        if (resultLoc[i] != MAX_T_ARRAYSIZE && verifyLoc[i] != resultLoc[i]) {
            cout << "Error: " << i << " query (" << queries[i * 2] << ", " << queries[i * 2 + 1] << ") - expected "
                 << verifyLoc[i] << " is " << resultLoc[i] << std::endl;
        }
    }
}