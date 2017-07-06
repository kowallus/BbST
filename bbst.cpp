#include <algorithm>
#include <iostream>
#include <numeric>
#include "bbst.h"
#include <omp.h>

BbST::BbST(vector<t_value> valuesArray, vector<t_array_size> queries, t_array_size *resultLoc, int kExp) {
    this->valuesArray = valuesArray;
    this->queries = queries;
    this->resultLoc = resultLoc;
    this->kExp = kExp;
    this->k = 1 << kExp;
}

void BbST::solve() {
    getBlocksMinsBase();
    getBlocksSparseTable();
    #pragma omp parallel for
    for (int i = 0; i < queries.size(); i = i + 2) {
        if (queries[i] == queries[i + 1]) {
            this->resultLoc[i / 2] = queries[i];
            continue;
        }
        this->resultLoc[i / 2] = getRangeMinLoc(queries[i], queries[i + 1]);
    }
    cleanup();
/**/
}

void BbST::verify() {
    cout << "Solution verification..." << std::endl;
    unsigned int i;
    const int q = queries.size() / 2;
    this->verifyLoc.resize(q);
    this->verifyVal.resize(q);
    t_value *const vaPtr = &this->valuesArray[0];
    t_array_size *const qPtr = &this->queries[0];
    for(i = 0; i < q; i++) {
        t_value *const begPtr = vaPtr + qPtr[2*i];
        t_value *const endPtr = vaPtr + qPtr[2*i + 1] + 1;
        t_value *const minValPtr = std::min_element(begPtr, endPtr);
        verifyVal[i] = *minValPtr;
        verifyLoc[i] = minValPtr - vaPtr;
        if (verifyLoc[i] != resultLoc[i]) {
            cout << "Error: " << i << " query (" << queries[i * 2] << ", " << queries[i * 2 + 1] << ") - expected "
                 << verifyLoc[i] << " is " << resultLoc[i] << std::endl;
        }
    }
}

void BbST::getBlocksMinsBase() {
    this->blocksCount = (valuesArray.size() + k - 1) >> kExp;
    this->D = 32 - __builtin_clz(blocksCount);
    const t_array_size blocksSize = blocksCount * D;
    blocksVal2D = new t_value[blocksSize];
    blocksLoc2D = new t_array_size[blocksSize];
    #pragma omp parallel for
    for (t_array_size i = 0; i < blocksCount - 1; i++) {
        auto minPtr = std::min_element(&valuesArray[i << kExp], &valuesArray[(i + 1) << kExp]);
        blocksVal2D[i] = *minPtr;
        blocksLoc2D[i] = minPtr - &valuesArray[0];
    }
    auto minPtr = std::min_element(&valuesArray[(blocksCount - 1) << kExp], &(*valuesArray.end()));
    blocksVal2D[blocksCount - 1] = *minPtr;
    blocksLoc2D[blocksCount - 1] = minPtr - &valuesArray[0];
}

void BbST::getBlocksSparseTable() {
    for(t_array_size e = 1, step = 1; e < D; ++e, step <<= 1) {
        for (t_array_size i = 0; i < blocksCount; i++) {
            t_array_size minIdx = i;
            const t_array_size e0offset = (e - 1) * blocksCount;
            if (i + step < blocksCount && blocksVal2D[(i + step) + e0offset] < blocksVal2D[i + e0offset]) {
                minIdx = i + step;
            }
            const t_array_size e1offset = e * blocksCount;
            blocksVal2D[i + e1offset] = blocksVal2D[minIdx + e0offset];
            blocksLoc2D[i + e1offset] = blocksLoc2D[minIdx + e0offset];
        }
    }
}

t_array_size BbST::getRangeMinLoc(const t_array_size &begIdx, const t_array_size &endIdx) {
    t_array_size result = -1;
    const t_array_size begCompIdx = begIdx >> kExp;
    const t_array_size endCompIdx = endIdx >> kExp;
    if (endCompIdx - begCompIdx <= 1) {
        return scanMinIdx(begIdx, endIdx);
    }

    t_array_size kBlockCount = endCompIdx - begCompIdx - 1;
    t_array_size e = 31 - __builtin_clz(kBlockCount);
    t_array_size step = 1 << e;
    t_value minVal = blocksVal2D[(begCompIdx + 1) + e * blocksCount];
    result = blocksLoc2D[(begCompIdx + 1) + e * blocksCount];
    t_array_size endShiftCompIdx = endCompIdx - step;
    if (endShiftCompIdx != begCompIdx + 1) {
        t_value temp = blocksVal2D[(endShiftCompIdx) + e * blocksCount];
        if (temp < minVal) {
            minVal = temp;
            result = blocksLoc2D[(endShiftCompIdx) + e * blocksCount];
        }
    }

    if (blocksVal2D[begCompIdx] <= minVal && begIdx != (begCompIdx + 1) << kExp) {
        t_array_size minIdx = scanMinIdx(begIdx, ((begCompIdx + 1) << kExp) - 1);
        if (valuesArray[minIdx] <= minVal) {
            minVal = valuesArray[minIdx];
            result = minIdx;
        }
    }
    if (blocksVal2D[endCompIdx] < minVal) {
       t_array_size minIdx = scanMinIdx(endCompIdx << kExp, endIdx);
        if (valuesArray[minIdx] < minVal) {
            result = minIdx;
        }
    }
    return result;
}

t_array_size BbST::scanMinIdx(const t_array_size &begIdx, const t_array_size &endIdx) {
    t_array_size minValIdx = begIdx;
    for(t_array_size i = begIdx + 1; i <= endIdx; i++) {
        if (valuesArray[i] < valuesArray[minValIdx]) {
            minValIdx = i;
        }
    }
    return minValIdx;
}

void BbST::cleanup() {
    delete[] this->blocksLoc2D;
    delete[] this->blocksVal2D;
}

size_t BbST::memUsageInBytes() {
    const t_array_size blocksCount = ((valuesArray.size() - 1 + k - 1) >> kExp);
    D = 32 - __builtin_clz(blocksCount);
    const t_array_size blocksSize = blocksCount * D;
    const size_t bytes = blocksSize * (sizeof(t_value) + sizeof(t_array_size));
    return bytes;
}
