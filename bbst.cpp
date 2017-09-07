#include <algorithm>
#include <iostream>
#include <numeric>
#include "bbst.h"
#include <omp.h>

#ifdef MINI_BLOCKS
BbST::BbST(int kExp, int miniKExp) {
    this->miniKExp = miniKExp;
    this->miniK = 1 << miniKExp;
#else
BbST::BbST(int kExp) {
#endif
    this->kExp = kExp;
    this->k = 1 << kExp;
}

BbST::~BbST() {
    if (batchMode) cleanup();
}

void BbST::rmqBatch(const t_value* valuesArray, const t_array_size n, const vector<t_array_size> &queries, t_array_size *resultLoc) {
    this->valuesArray = valuesArray;
    this->n = n;
    getBlocksMinsBase();
    getBlocksSparseTable();
    this->rmqBatch(queries, resultLoc);
    cleanup();
/**/
}

#ifdef MINI_BLOCKS
BbST::BbST(const t_value* valuesArray, const t_array_size n, int kExp, int miniKExp) {
    this->miniKExp = miniKExp;
    this->miniK = 1 << miniKExp;
#else
BbST::BbST(const t_value* valuesArray, const t_array_size n, int kExp) {
#endif
    this->valuesArray = valuesArray;
    this->n = n;
    this->kExp = kExp;
    this->k = 1 << kExp;
    this->batchMode = true;
    getBlocksMinsBase();
    getBlocksSparseTable();
}

void BbST::rmqBatch(const vector<t_array_size> &queries, t_array_size *resultLoc) {
    #pragma omp parallel for
    for (int i = 0; i < queries.size(); i = i + 2) {
        resultLoc[i / 2] = rmq(queries[i], queries[i + 1]);
    }
}

void BbST::getBlocksMinsBase() {
#ifdef MINI_BLOCKS
    this->miniBlocksCount = (n + miniK - 1) >> miniKExp;
    this->miniBlocksLoc = new uint8_t[miniBlocksCount];
    this->miniBlocksInBlock = k / miniK;
#endif
    this->blocksCount = (n + k - 1) >> kExp;
    this->D = 32 - __builtin_clz(blocksCount);
    const t_array_size blocksSize = blocksCount * D;
    blocksVal2D = new t_value[blocksSize];
    blocksLoc2D = new t_array_size[blocksSize];
    #pragma omp parallel for
    for (t_array_size i = 0; i < blocksCount - 1; i++) {
#ifdef MINI_BLOCKS
        t_array_size miniI = i << (kExp - miniKExp);
        for (t_array_size j = 0; j < miniBlocksInBlock; j++, miniI++) {
            auto miniMinPtr = std::min_element(&valuesArray[miniI << miniKExp], &valuesArray[(miniI + 1) << miniKExp]);
            miniBlocksLoc[miniI] = miniMinPtr - &valuesArray[miniI << miniKExp];
        }
#endif
        auto minPtr = std::min_element(&valuesArray[i << kExp], &valuesArray[(i + 1) << kExp]);
        blocksVal2D[i] = *minPtr;
        blocksLoc2D[i] = minPtr - &valuesArray[0];
    }
#ifdef MINI_BLOCKS
    t_array_size miniI = (blocksCount - 1) << (kExp - miniKExp);
    for (; ((miniI + 1) << miniKExp) < n; miniI++) {
        auto miniMinPtr = std::min_element(&valuesArray[miniI << miniKExp], &valuesArray[(miniI + 1) << miniKExp]);
        miniBlocksLoc[miniI] = miniMinPtr - &valuesArray[miniI << miniKExp];
    }
    auto miniMinPtr = std::min_element(&valuesArray[miniI << miniKExp], &valuesArray[n]);
    miniBlocksLoc[miniI] = miniMinPtr - &valuesArray[miniI << miniKExp];
#endif
    auto minPtr = std::min_element(&valuesArray[(blocksCount - 1) << kExp], &valuesArray[n]);
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

t_array_size BbST::rmq(const t_array_size &begIdx, const t_array_size &endIdx) {
    if (begIdx == endIdx)
        return begIdx;
    t_array_size result = -1;
    const t_array_size begCompIdx = begIdx >> kExp;
    const t_array_size endCompIdx = endIdx >> kExp;
    const t_array_size firstBlockMinLoc = blocksLoc2D[begCompIdx];
    if (endCompIdx == begCompIdx) {
        if (begIdx <= firstBlockMinLoc && firstBlockMinLoc <= endIdx)
            return firstBlockMinLoc;
        return scanMinIdx(begIdx, endIdx);
    }
    t_value minVal = MAX_T_VALUE;
    if (endCompIdx - begCompIdx > 1) {
        t_array_size kBlockCount = endCompIdx - begCompIdx - 1;
        t_array_size e = 31 - __builtin_clz(kBlockCount);
        t_array_size step = 1 << e;
        minVal = blocksVal2D[(begCompIdx + 1) + e * blocksCount];
        result = blocksLoc2D[(begCompIdx + 1) + e * blocksCount];
        t_array_size endShiftCompIdx = endCompIdx - step;
        if (endShiftCompIdx != begCompIdx + 1) {
            t_value temp = blocksVal2D[(endShiftCompIdx) + e * blocksCount];
            if (temp < minVal) {
                minVal = temp;
                result = blocksLoc2D[(endShiftCompIdx) + e * blocksCount];
            }
        }
    }
#ifndef WORST_CASE
    if (blocksVal2D[begCompIdx] <= minVal && begIdx != (begCompIdx + 1) << kExp) {
        if (firstBlockMinLoc >= begIdx) {
            minVal = blocksVal2D[begCompIdx];
            result = firstBlockMinLoc;
        } else {
            t_array_size minIdx = scanMinIdx(begIdx, ((begCompIdx + 1) << kExp) - 1);
            if (valuesArray[minIdx] <= minVal) {
                minVal = valuesArray[minIdx];
                result = minIdx;
            }
        }
    }
    if (blocksVal2D[endCompIdx] < minVal) {
        t_array_size lastBlockMinLoc = blocksLoc2D[endCompIdx];
        if (lastBlockMinLoc <= endIdx) {
            result = lastBlockMinLoc;
        } else {
            t_array_size minIdx = scanMinIdx(endCompIdx << kExp, endIdx);
            if (valuesArray[minIdx] < minVal) {
                result = minIdx;
            }
        }
    }
#else
    t_array_size minIdx = scanMinIdx(begIdx, ((begCompIdx + 1) << kExp) - 1);
    if (valuesArray[minIdx] <= minVal) {
        minVal = valuesArray[minIdx];
        result = minIdx;
    }
    minIdx = scanMinIdx(endCompIdx << kExp, endIdx);
    if (valuesArray[minIdx] < minVal) {
        result = minIdx;
    }
#endif
    return result;
}

inline t_array_size BbST::rawScanMinIdx(const t_array_size &begIdx, const t_array_size &endIdx) {
    t_array_size minValIdx = begIdx;
    for(t_array_size i = begIdx + 1; i <= endIdx; i++) {
        if (valuesArray[i] < valuesArray[minValIdx]) {
            minValIdx = i;
        }
    }
    return minValIdx;
}

inline t_array_size BbST::miniScanMinIdx(const t_array_size &begIdx, const t_array_size &endIdx) {
    t_array_size result = -1;
    const t_array_size begMiniIdx = begIdx >> miniKExp;
    const t_array_size endMiniIdx = endIdx >> miniKExp;
    const t_array_size firstMiniBlockMinLoc = (begMiniIdx << miniKExp) + miniBlocksLoc[begMiniIdx];
    if (endMiniIdx == begMiniIdx) {
        if (begIdx <= firstMiniBlockMinLoc && firstMiniBlockMinLoc <= endIdx)
            return firstMiniBlockMinLoc;
        return rawScanMinIdx(begIdx, endIdx);
    }
    t_value minVal = MAX_T_VALUE;
    if (endMiniIdx - begMiniIdx > 1) {
        for(t_array_size i = begMiniIdx + 1; i <= endMiniIdx - 1; i++) {
            t_array_size tempLoc = (i << miniKExp) + miniBlocksLoc[i];
            t_value tempVal = valuesArray[tempLoc];
            if (tempVal < minVal) {
                minVal = tempVal;
                result = tempLoc;
            }
        }
    }
#ifndef WORST_CASE
    t_value tempVal = valuesArray[firstMiniBlockMinLoc];
    if (tempVal <= minVal && begIdx != (begMiniIdx + 1) << miniKExp) {
        if (firstMiniBlockMinLoc >= begIdx) {
            minVal = tempVal;
            result = firstMiniBlockMinLoc;
        } else {
            t_array_size minIdx = rawScanMinIdx(begIdx, ((begMiniIdx + 1) << miniKExp) - 1);
            if (valuesArray[minIdx] <= minVal) {
                minVal = valuesArray[minIdx];
                result = minIdx;
            }
        }
    }
    t_array_size lastBlockMinLoc = (endMiniIdx << miniKExp) + miniBlocksLoc[endMiniIdx];
    tempVal = valuesArray[lastBlockMinLoc];
    if (tempVal < minVal) {
        if (lastBlockMinLoc <= endIdx) {
            result = lastBlockMinLoc;
        } else {
            t_array_size minIdx = rawScanMinIdx(endMiniIdx << miniKExp, endIdx);
            if (valuesArray[minIdx] < minVal) {
                result = minIdx;
            }
        }
    }
#else
    t_array_size minIdx = rawScanMinIdx(begIdx, ((begMiniIdx + 1) << miniKExp) - 1);
    if (valuesArray[minIdx] <= minVal) {
        minVal = valuesArray[minIdx];
        result = minIdx;
    }
    minIdx = rawScanMinIdx(endMiniIdx << miniKExp, endIdx);
    if (valuesArray[minIdx] < minVal) {
        result = minIdx;
    }
#endif
    return result;
}

t_array_size BbST::scanMinIdx(const t_array_size &begIdx, const t_array_size &endIdx) {
#ifdef MINI_BLOCKS
    return miniScanMinIdx(begIdx, endIdx);
#else
    return rawScanMinIdx(begIdx, endIdx);
#endif
}

void BbST::cleanup() {
    delete[] this->blocksLoc2D;
    delete[] this->blocksVal2D;
#ifdef MINI_BLOCKS
    delete[] this->miniBlocksLoc;
#endif
}

size_t BbST::memUsageInBytes() {
    D = 32 - __builtin_clz(blocksCount);
    const t_array_size blocksSize = blocksCount * D;
    size_t bytes = blocksSize * (sizeof(t_value) + sizeof(t_array_size));
#ifdef MINI_BLOCKS
    bytes += miniBlocksCount;
#endif
    return bytes;
}

