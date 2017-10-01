#include <algorithm>
#include <iostream>
#include <numeric>

#include "bbstx.h"

#include <omp.h>

BbSTx::~BbSTx() {
    cleanup();
}

#ifdef MINI_BLOCKS
BbSTx::BbSTx(const vector<t_value> &valuesArray, int kExp, int miniKExp, RMQAPI* secondaryRMQ) {
    this->miniKExp = miniKExp;
    this->miniK = 1 << miniKExp;
#else
BbSTx::BbSTx(const vector<t_value> &valuesArray, int kExp, RMQAPI* secondaryRMQ) {
#endif
    this->kExp = kExp;
    this->k = 1 << kExp;
    this->secondaryRMQ = secondaryRMQ;
    getBlocksMinsBase(valuesArray);
    getBlocksSparseTable();
}

void BbSTx::rmqBatch(const vector<t_array_size> &queries, t_array_size *resultLoc) {
    for (int i = 0; i < queries.size(); i = i + 2) {
        resultLoc[i / 2] = rmq(queries[i], queries[i + 1]);
    }
}

void BbSTx::getBlocksMinsBase(const vector<t_value> &valuesArray) {
#ifdef MINI_BLOCKS
    this->miniBlocksCount = (valuesArray.size() + miniK - 1) >> miniKExp;
    this->miniBlocksLoc = new uint8_t[miniBlocksCount];
    this->miniBlocksVal = new t_value[miniBlocksCount];
    this->miniBlocksInBlock = k / miniK;
#endif
    this->blocksCount = (valuesArray.size() + k - 1) >> kExp;
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
            miniBlocksVal[miniI] = *miniMinPtr;
        }
#endif
        auto minPtr = std::min_element(&valuesArray[i << kExp], &valuesArray[(i + 1) << kExp]);
        blocksVal2D[i] = *minPtr;
        blocksLoc2D[i] = minPtr - &valuesArray[0];
    }
#ifdef MINI_BLOCKS
    t_array_size miniI = (blocksCount - 1) << (kExp - miniKExp);
    for (; ((miniI + 1) << miniKExp) < valuesArray.size(); miniI++) {
        auto miniMinPtr = std::min_element(&valuesArray[miniI << miniKExp], &valuesArray[(miniI + 1) << miniKExp]);
        miniBlocksLoc[miniI] = miniMinPtr - &valuesArray[miniI << miniKExp];
        miniBlocksVal[miniI] = *miniMinPtr;
    }
    auto miniMinPtr = std::min_element(&valuesArray[miniI << miniKExp], &(*valuesArray.end()));
    miniBlocksLoc[miniI] = miniMinPtr - &valuesArray[miniI << miniKExp];
    miniBlocksVal[miniI] = *miniMinPtr;
#endif
    auto minPtr = std::min_element(&valuesArray[(blocksCount - 1) << kExp], &(*valuesArray.end()));
    blocksVal2D[blocksCount - 1] = *minPtr;
    blocksLoc2D[blocksCount - 1] = minPtr - &valuesArray[0];
}

void BbSTx::getBlocksSparseTable() {
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

t_array_size BbSTx::rmq(const t_array_size &begIdx, const t_array_size &endIdx) {
    if (begIdx == endIdx) {
        return begIdx;
    }
    t_array_size result = MAX_T_ARRAYSIZE;
    const t_array_size begCompIdx = begIdx >> kExp;
    const t_array_size endCompIdx = endIdx >> kExp;
/* OPTIMIZATION FOR NARROW RANGES
    if (endCompIdx == begCompIdx) {
        const t_array_size result = blocksLoc2D[begCompIdx];
        if (begIdx <= result && result <= endIdx)
            return result;
#ifndef MINI_BLOCKS
        return secondaryRMQ->rmq(begIdx, endIdx);
#else
        t_value miniNotSmallerThen = MAX_T_VALUE;
        t_array_size minIdx = miniScanMinIdx(begIdx, endIdx, miniNotSmallerThen);
        if (minIdx == MAX_T_ARRAYSIZE) {
            return secondaryRMQ->rmq(begIdx, endIdx);
        }
        return minIdx;
#endif
    }
*/
    const t_array_size kBlockCount = endCompIdx - begCompIdx; // actual kBlock count is +1
    const t_array_size e = kBlockCount?(31 - __builtin_clz(kBlockCount)):0;
    const t_array_size step = 1 << e;
    const t_array_size endShiftCompIdx = endCompIdx - step + 1;
    t_value leftMin = blocksVal2D[begCompIdx + e * blocksCount];
    t_value rightMin = blocksVal2D[endShiftCompIdx + e * blocksCount];
    bool minOnTheLeft = leftMin <= rightMin;
    result = blocksLoc2D[(minOnTheLeft?begCompIdx:endShiftCompIdx) + e * blocksCount];
    if (begIdx <= result && result <= endIdx)
        return result;
#ifndef MINI_BLOCKS
    return secondaryRMQ->rmq(begIdx, endIdx);
#else
    t_value miniNotSmallerThen = MAX_T_VALUE;
    if (kBlockCount <= 1) {
        const t_array_size minIdx = miniScanMinIdx(begIdx, endIdx, miniNotSmallerThen);
//        cout << "(" << begIdx << ", " << endIdx << ") " << minIdx << " " << miniNotSmallerThen << endl;

        if (minIdx == MAX_T_ARRAYSIZE)
            return secondaryRMQ->rmq(begIdx, endIdx);
        return minIdx;
    }
    bool uncertainMini = false;
    result = blocksLoc2D[begCompIdx + e * blocksCount];
    t_value minVal;
    if (result >= begIdx) {
        minVal = leftMin;
    } else {
        const t_array_size inner2DbegIdx = begCompIdx + 1 + (e - (step == kBlockCount)) * blocksCount;
        minVal = blocksVal2D[inner2DbegIdx];
        result = blocksLoc2D[inner2DbegIdx];
        const t_array_size minIdx = miniScanMinIdx(begIdx, ((begCompIdx + 1) << kExp) - 1, miniNotSmallerThen);
        if (minIdx == MAX_T_ARRAYSIZE) {
            if (miniNotSmallerThen <= minVal) {
                uncertainMini = true;
                minVal = miniNotSmallerThen;
            }
        } else {
            const t_value tempVal = miniBlocksVal[minIdx >> miniKExp];
            if (tempVal <= minVal) {
                minVal = tempVal;
                result = minIdx;
            }
        }
    }

    if (rightMin < minVal) {
        t_array_size tempLoc = blocksLoc2D[endShiftCompIdx + e * blocksCount];
        if (tempLoc <= endIdx) {
            return tempLoc;
        } else {
            const t_array_size inner2DEndShiftIdx = (endCompIdx - (step >> (step == kBlockCount)) + ((e - (step == kBlockCount)) * blocksCount));
            t_value tempVal = blocksVal2D[inner2DEndShiftIdx];
            if (tempVal < minVal) {
                uncertainMini = false;
                minVal = tempVal;
                result = blocksLoc2D[inner2DEndShiftIdx];
            }
            const t_array_size minIdx = miniScanMinIdx(endCompIdx << kExp, endIdx, miniNotSmallerThen);
            if (minIdx == MAX_T_ARRAYSIZE) {
                if (miniNotSmallerThen < minVal)
                    return secondaryRMQ->rmq(begIdx, endIdx);
            } else if (miniBlocksVal[minIdx >> miniKExp] < minVal) {
                return minIdx;
            }
        }
    }

    if (uncertainMini)
        return secondaryRMQ->rmq(begIdx, endIdx);

    return result;

#endif
}

inline t_array_size BbSTx::miniScanMinIdx(const t_array_size &begIdx, const t_array_size &endIdx, t_value &notSmallerThan) {
    t_array_size result = MAX_T_ARRAYSIZE;
    const t_array_size begMiniIdx = begIdx >> miniKExp;
    const t_array_size endMiniIdx = endIdx >> miniKExp;
    if (endMiniIdx == begMiniIdx) {
        const t_array_size firstMiniBlockMinLoc = (begMiniIdx << miniKExp) + miniBlocksLoc[begMiniIdx];
        if (begIdx <= firstMiniBlockMinLoc && firstMiniBlockMinLoc <= endIdx)
            return firstMiniBlockMinLoc;
        notSmallerThan = miniBlocksVal[begMiniIdx];
        return MAX_T_ARRAYSIZE;
    }
    t_value minVal = MAX_T_VALUE;
    if (endMiniIdx - begMiniIdx > 1) {
        for(t_array_size i = begMiniIdx + 1; i <= endMiniIdx - 1; i++) {
            t_value tempVal = miniBlocksVal[i];
            if (tempVal < minVal) {
                minVal = tempVal;
                result = (i << miniKExp) + miniBlocksLoc[i];
            }
        }
    }
    bool uncertainMini = false;
    t_value tempVal = miniBlocksVal[begMiniIdx];
    if (tempVal <= minVal && begIdx != (begMiniIdx + 1) << miniKExp) {
        const t_array_size firstMiniBlockMinLoc = (begMiniIdx << miniKExp) + miniBlocksLoc[begMiniIdx];
        minVal = tempVal;
        if (firstMiniBlockMinLoc >= begIdx) {
            result = firstMiniBlockMinLoc;
        } else {
            uncertainMini = true;
        }
    }
    tempVal = miniBlocksVal[endMiniIdx];
    if (tempVal < minVal) {
        t_array_size lastBlockMinLoc = (endMiniIdx << miniKExp) + miniBlocksLoc[endMiniIdx];
        if (lastBlockMinLoc <= endIdx) {
            return lastBlockMinLoc;
        } else {
            notSmallerThan = tempVal;
            return MAX_T_ARRAYSIZE;
        }
    }

    if (uncertainMini) {
        notSmallerThan = minVal;
        return MAX_T_ARRAYSIZE;
    }
    return result;
}

void BbSTx::cleanup() {
    delete[] this->blocksLoc2D;
    delete[] this->blocksVal2D;
#ifdef MINI_BLOCKS
    delete[] this->miniBlocksLoc;
    delete[] this->miniBlocksVal;
#endif
}

size_t BbSTx::memUsageInBytes() {
    const t_array_size blocksSize = blocksCount * D;
    size_t bytes = blocksSize * (sizeof(t_value) + sizeof(t_array_size));
#ifdef MINI_BLOCKS
    bytes += miniBlocksCount * (sizeof(uint8_t) + sizeof(t_value));
#endif
    return bytes + secondaryRMQ->memUsageInBytes();
}

