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
    t_value* miniBlocksVal= new t_value[miniBlocksCount];
    this->miniBlocksInBlock = k / miniK;
#endif
    this->blocksCount = (n + k - 1) >> kExp;
    this->D = 32 - __builtin_clz(blocksCount);
    const t_array_size blocksSize = blocksCount * D;
    blocksVal2D = new t_value[blocksSize];
    blocksLoc2D = new t_array_size[blocksSize];

#ifdef MINI_BLOCKS
    #pragma omp parallel for
    for (t_array_size miniI = 0; miniI < this->miniBlocksCount - 1; miniI++) {
        auto miniMinPtr = std::min_element(&valuesArray[miniI << miniKExp], &valuesArray[(miniI + 1) << miniKExp]);
        miniBlocksVal[miniI] = *miniMinPtr;
        miniBlocksLoc[miniI] = miniMinPtr - &valuesArray[miniI << miniKExp];
    }
    int delKExp = kExp - miniKExp;
    #pragma omp parallel for
    for (t_array_size i = 0; i < blocksCount - 1; i++) {
        auto minPtr = std::min_element(&miniBlocksVal[i << delKExp], &miniBlocksVal[(i + 1) << delKExp]);
        blocksVal2D[i] = *minPtr;
        t_array_size miniI = minPtr - &miniBlocksVal[0];
        blocksLoc2D[i] = (miniI << miniKExp) + miniBlocksLoc[miniI];
    }
#else
    #pragma omp parallel for
    for (t_array_size i = 0; i < blocksCount - 1; i++) {
        auto minPtr = std::min_element(&valuesArray[i << kExp], &valuesArray[(i + 1) << kExp]);
        blocksVal2D[i] = *minPtr;
        blocksLoc2D[i] = minPtr - &valuesArray[0];
    }
#endif
#ifdef MINI_BLOCKS
    t_array_size miniI = (blocksCount - 1) << (kExp - miniKExp);
    for (; miniI < this->miniBlocksCount - 1; miniI++) {
        auto miniMinPtr = std::min_element(&valuesArray[miniI << miniKExp], &valuesArray[(miniI + 1) << miniKExp]);
        miniBlocksVal[miniI] = *miniMinPtr;
        miniBlocksLoc[miniI] = miniMinPtr - &valuesArray[miniI << miniKExp];
    }
    auto miniMinPtr = std::min_element(&valuesArray[miniI << miniKExp], &valuesArray[n]);
    miniBlocksLoc[miniI] = miniMinPtr - &valuesArray[miniI << miniKExp];
    miniBlocksVal[miniI] = *miniMinPtr;
    auto minPtr = std::min_element(&miniBlocksVal[blocksCount - 1 << delKExp], &miniBlocksVal[miniBlocksCount]);
    blocksVal2D[blocksCount - 1] = *minPtr;
    miniI = minPtr - &miniBlocksVal[0];
    blocksLoc2D[blocksCount - 1] = (miniI << miniKExp) + miniBlocksLoc[miniI];
    delete miniBlocksVal;
#else
    auto minPtr = std::min_element(&valuesArray[(blocksCount - 1) << kExp], &valuesArray[n]);
    blocksVal2D[blocksCount - 1] = *minPtr;
    blocksLoc2D[blocksCount - 1] = minPtr - &valuesArray[0];
#endif
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
    if (begIdx == endIdx) {
        return begIdx;
    }
    t_array_size result = MAX_T_ARRAYSIZE;
    const t_array_size begCompIdx = begIdx >> kExp;
    const t_array_size endCompIdx = endIdx >> kExp;
#ifdef START_FROM_NARROW_RANGES
    if (endCompIdx == begCompIdx) {
        const t_array_size result = blocksLoc2D[begCompIdx];
        if (begIdx <= result && result <= endIdx)
            return result;
        t_value minVal = MAX_T_VALUE;
        return scanMinIdx(begIdx, endIdx, miniSmallerThen, true);
    }
#endif
    const t_array_size kBlockCount = endCompIdx - begCompIdx; // actual kBlock count is +1
    const t_array_size e = kBlockCount?(31 - __builtin_clz(kBlockCount)):0;
    const t_array_size step = 1 << e;
    const t_array_size endShiftCompIdx = endCompIdx - step + 1;
    t_value leftMin = blocksVal2D[begCompIdx + e * blocksCount];
    t_value rightMin = blocksVal2D[endShiftCompIdx + e * blocksCount];
    bool minOnTheLeft = leftMin <= rightMin;
    result = blocksLoc2D[(minOnTheLeft?begCompIdx:endShiftCompIdx) + e * blocksCount];
#ifndef WORST_CASE
    if (begIdx <= result && result <= endIdx)
        return result;
#endif
    t_value minVal = MAX_T_VALUE;
    if (kBlockCount <= 1) {
        return scanMinIdx(begIdx, endIdx, minVal, true);
    }
#ifndef WORST_CASE
    result = blocksLoc2D[begCompIdx + e * blocksCount];
    if (result >= begIdx) {
        minVal = leftMin;
    } else
#endif
    {
        const t_array_size inner2DbegIdx = begCompIdx + 1 + (e - (step == kBlockCount)) * blocksCount;
        minVal = blocksVal2D[inner2DbegIdx];
        result = blocksLoc2D[inner2DbegIdx];
        const t_array_size minIdx = scanMinIdx(begIdx, ((begCompIdx + 1) << kExp) - 1, minVal, true);
        if (minIdx != MAX_T_ARRAYSIZE)
            result = minIdx;
    }

#ifndef WORST_CASE
    t_array_size tempLoc;
    if (rightMin < minVal)
        if ((tempLoc = blocksLoc2D[endShiftCompIdx + e * blocksCount]) <= endIdx) {
            return tempLoc;
        } else
#endif
        {
            const t_array_size inner2DEndShiftIdx = (endCompIdx - (step >> (step == kBlockCount)) + ((e - (step == kBlockCount)) * blocksCount));
            t_value tempVal = blocksVal2D[inner2DEndShiftIdx];
            if (tempVal < minVal) {
                minVal = tempVal;
                result = blocksLoc2D[inner2DEndShiftIdx];
            }
            const t_array_size minIdx = scanMinIdx(endCompIdx << kExp, endIdx, minVal, false);
            if (minIdx != MAX_T_ARRAYSIZE)
                return minIdx;
        }

    return result;
}

inline t_array_size BbST::rawScanMinIdx(const t_array_size &begIdx, const t_array_size &endIdx, t_value& minVal, bool smallerOrEqual) {
    t_array_size minValIdx = begIdx;
    for(t_array_size i = begIdx + 1; i <= endIdx; i++) {
        if (valuesArray[i] < valuesArray[minValIdx]) {
            minValIdx = i;
        }
    }
    if (smallerOrEqual?valuesArray[minValIdx]<=minVal:valuesArray[minValIdx]<minVal) {
        minVal = valuesArray[minValIdx];
        return minValIdx;
    } else
        return MAX_T_ARRAYSIZE;
}

inline t_array_size BbST::miniScanMinIdx(const t_array_size &begIdx, const t_array_size &endIdx, t_value& minVal, bool smallerOrEqual) {
    const t_array_size begMiniIdx = begIdx >> miniKExp;
    const t_array_size endMiniIdx = endIdx >> miniKExp;
    const t_array_size firstMiniBlockMinLoc = (begMiniIdx << miniKExp) + miniBlocksLoc[begMiniIdx];
    if (endMiniIdx == begMiniIdx) {
#ifndef WORST_CASE
        if (begIdx <= firstMiniBlockMinLoc && firstMiniBlockMinLoc <= endIdx) {
            if (smallerOrEqual?valuesArray[firstMiniBlockMinLoc]<=minVal:valuesArray[firstMiniBlockMinLoc]<minVal) {
                minVal = valuesArray[firstMiniBlockMinLoc];
                return firstMiniBlockMinLoc;
            } else
                return MAX_T_ARRAYSIZE;
        }
#endif
        return rawScanMinIdx(begIdx, endIdx, minVal, smallerOrEqual);
    }
    t_array_size result = MAX_T_ARRAYSIZE;
    if (endMiniIdx - begMiniIdx > 1) {
        t_value innerMinVal = minVal;
        for(t_array_size i = endMiniIdx - 1; i > begMiniIdx; i--) {
            t_array_size tempLoc = (i << miniKExp) + miniBlocksLoc[i];
            t_value tempVal = valuesArray[tempLoc];
            if (tempVal <= innerMinVal) {
                innerMinVal = tempVal;
                result = tempLoc;
            }
        }
        if (innerMinVal < minVal) {
            smallerOrEqual = true;
            minVal = innerMinVal;
        } else if (result != MAX_T_ARRAYSIZE && !smallerOrEqual)
            result = MAX_T_ARRAYSIZE;
    }
#ifndef WORST_CASE
    t_value tempVal = valuesArray[firstMiniBlockMinLoc];
    if ((smallerOrEqual?tempVal <= minVal:tempVal < minVal) && begIdx != (begMiniIdx + 1) << miniKExp)
        if (firstMiniBlockMinLoc >= begIdx) {
            smallerOrEqual = false;
            minVal = tempVal;
            result = firstMiniBlockMinLoc;
        } else
#endif
        {
            t_array_size minIdx = rawScanMinIdx(begIdx, ((begMiniIdx + 1) << miniKExp) - 1, minVal, smallerOrEqual);
            if (minIdx != MAX_T_ARRAYSIZE) {
                smallerOrEqual = false;
                result = minIdx;
            }
        }

#ifndef WORST_CASE
    t_array_size lastBlockMinLoc = (endMiniIdx << miniKExp) + miniBlocksLoc[endMiniIdx];
    tempVal = valuesArray[lastBlockMinLoc];
    if ((smallerOrEqual && result > lastBlockMinLoc)?tempVal <= minVal:tempVal < minVal)
        if (lastBlockMinLoc <= endIdx) {
            minVal = tempVal;
            return lastBlockMinLoc;
        } else
#endif
        {
            t_array_size minIdx = rawScanMinIdx(endMiniIdx << miniKExp, endIdx, minVal, result == MAX_T_ARRAYSIZE && smallerOrEqual);
            if (minIdx != MAX_T_ARRAYSIZE) {
                return minIdx;
            }
        }

    return result;
}

t_array_size BbST::scanMinIdx(const t_array_size &begIdx, const t_array_size &endIdx, t_value& minVal, bool smallerOrEqual) {
#ifdef MINI_BLOCKS
    return miniScanMinIdx(begIdx, endIdx, minVal, smallerOrEqual);
#else
    return rawScanMinIdx(begIdx, endIdx, minVal, smallerOrEqual);
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

