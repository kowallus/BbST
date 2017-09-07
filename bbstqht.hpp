#include <algorithm>
#include <iostream>
#include <numeric>
#include "math.h"

#include <omp.h>

template<typename t_qvalue, int max_qvalue> BbSTqht<t_qvalue, max_qvalue>::~BbSTqht() {
    cleanup();
}

#ifdef MINI_BLOCKS
template<typename t_qvalue, int max_qvalue> BbSTqht<t_qvalue, max_qvalue>::BbSTqht(const vector<t_value> &valuesArray, int kExp, int miniKExp, RMQAPI* secondaryRMQ) {
    this->miniKExp = miniKExp;
    this->miniK = 1 << miniKExp;
#else
template<typename t_qvalue, int max_qvalue> BbSTqht<t_qvalue, max_qvalue>::BbSTqht(const vector<t_value> &valuesArray, int kExp, RMQAPI* secondaryRMQ) {
#endif
    this->kExp = kExp;
    this->k = 1 << kExp;
    this->secondaryRMQ = secondaryRMQ;
    prepareMinTables(valuesArray);
}

template<typename t_qvalue, int max_qvalue> void BbSTqht<t_qvalue, max_qvalue>::rmqBatch(const vector<t_array_size> &queries, t_array_size *resultLoc) {
    for (int i = 0; i < queries.size(); i = i + 2) {
        resultLoc[i / 2] = rmq(queries[i], queries[i + 1]);
    }
}

template<typename t_qvalue> inline t_qvalue quantizeValue(const t_value value, const t_qvalue max_qvalue, const t_value &minMinVal, const t_value &maxMinVal) {
//    double valMinMaxRatio = 1 - ((((double) maxMinVal - value)) / (((double) maxMinVal - minMinVal)));
//    double valMinMaxRatio = 1 - ((((double) maxMinVal - value)*(maxMinVal - value)) / (((double) maxMinVal - minMinVal)*(maxMinVal - minMinVal)));
//    double valMinMaxRatio = 1 - ((((double) maxMinVal - value)*(maxMinVal - value)*(maxMinVal - value)) / (((double) maxMinVal - minMinVal)*(maxMinVal - minMinVal)*(maxMinVal - minMinVal)));
//    double ratioExp = 5;
//    double valMinMaxRatio = 1 - (pow((double) maxMinVal - value, ratioExp) / (pow((double) maxMinVal - minMinVal, ratioExp)));
    double a = maxMinVal - value;
    a *= a;
    a *= a * a;
    double b = maxMinVal - minMinVal;
    b *= b;
    b *= b * b;
    double valMinMaxRatio = 1 - a / b;
    return (max_qvalue - 1) * valMinMaxRatio;
}

template<typename t_qvalue, int max_qvalue> void BbSTqht<t_qvalue, max_qvalue>::prepareMinTables(const vector<t_value> &valuesArray) {
#ifdef MINI_BLOCKS
    this->miniBlocksCount = (valuesArray.size() + miniK - 1) >> miniKExp;
    this->miniBlocksLoc = new uint8_t[miniBlocksCount];
    vector<t_value> miniBlocksVal(miniBlocksCount);
    this->miniBlocksQVal = new t_qvalue[miniBlocksCount];
    this->miniBlocksInBlock = k / miniK;
#endif
    this->blocksCount = (valuesArray.size() + k - 1) >> kExp;
    this->D = 32 - __builtin_clz(blocksCount);
    const t_array_size blocksSize = blocksCount * D;
    vector<t_value> blocksVal2D(blocksSize);
    blocksQVal2D = new t_qvalue[blocksSize];
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
    t_value minMinVal = blocksVal2D[0];
    t_value maxMinVal = blocksVal2D[0];
    for (t_array_size i = 1; i < blocksCount; i++)
        if (minMinVal > blocksVal2D[i])
            minMinVal = blocksVal2D[i];
        else if (maxMinVal < blocksVal2D[i])
            maxMinVal = blocksVal2D[i];
#ifdef MINI_BLOCKS
    for (t_array_size miniI = 0; miniI < miniBlocksCount; miniI++) {
        if (miniBlocksVal[miniI] >= maxMinVal)
            miniBlocksQVal[miniI] = (max_qvalue - 1);
        else
            miniBlocksQVal[miniI] = quantizeValue(miniBlocksVal[miniI], max_qvalue, minMinVal, maxMinVal);
    }
#endif
    for (t_array_size i = 0; i < blocksCount; i++) {
        blocksQVal2D[i] = quantizeValue(blocksVal2D[i], max_qvalue, minMinVal, maxMinVal);
    }
    prepareBlocksSparseTable(blocksVal2D, minMinVal, maxMinVal);
}

template<typename t_qvalue, int max_qvalue> void BbSTqht<t_qvalue, max_qvalue>::prepareBlocksSparseTable(vector<t_value> &blocksVal2D, const t_value &minMinVal, const t_value &maxMinVal) {
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
            blocksQVal2D[i + e1offset] = quantizeValue(blocksVal2D[i + e1offset], max_qvalue, minMinVal, maxMinVal);
        }
    }
}

template<typename t_qvalue, int max_qvalue> t_array_size BbSTqht<t_qvalue, max_qvalue>::rmq(const t_array_size &begIdx, const t_array_size &endIdx) {
    if (begIdx == endIdx) {
        return begIdx;
    }
    t_array_size result = -1;
    const t_array_size begCompIdx = begIdx >> kExp;
    const t_array_size endCompIdx = endIdx >> kExp;
    if (endCompIdx == begCompIdx) {
        const t_array_size firstBlockMinLoc = blocksLoc2D[begCompIdx];
        if (begIdx <= firstBlockMinLoc && firstBlockMinLoc <= endIdx)
            return firstBlockMinLoc;
#ifdef MINI_BLOCKS
        t_array_size minIdx = miniScanMinIdx(begIdx, endIdx);
        if (minIdx == MAX_T_ARRAYSIZE) {
            return secondaryRMQ->rmq(begIdx, endIdx);
        }
        return minIdx;
#else
        return secondaryRMQ->rmq(begIdx, endIdx);
#endif
    }
    t_qvalue minQVal = max_qvalue;
    bool uncertainResult = false;
    if (endCompIdx - begCompIdx > 1) {
        t_array_size kBlockCount = endCompIdx - begCompIdx - 1;
        t_array_size e = 31 - __builtin_clz(kBlockCount);
        t_array_size step = 1 << e;
        minQVal = blocksQVal2D[(begCompIdx + 1) + e * blocksCount];
        result = blocksLoc2D[(begCompIdx + 1) + e * blocksCount];
        t_array_size endShiftCompIdx = endCompIdx - step;
        if (endShiftCompIdx != begCompIdx + 1) {
            t_value temp = blocksQVal2D[(endShiftCompIdx) + e * blocksCount];
            if (temp < minQVal) {
                minQVal = temp;
                result = blocksLoc2D[(endShiftCompIdx) + e * blocksCount];
            } else {
                if (temp == minQVal && result != blocksLoc2D[(endShiftCompIdx) + e * blocksCount])
                    uncertainResult = true;
            }
        }
    }

    if (begIdx != (begCompIdx + 1) << kExp) {
        if (blocksQVal2D[begCompIdx] == minQVal) {
#ifdef MINI_BLOCKS
            const t_array_size firstBlockMinLoc = blocksLoc2D[begCompIdx];
            if (firstBlockMinLoc >= begIdx) {
                uncertainResult = true;
            } else {
                t_array_size minIdx = miniScanMinIdx(begIdx, ((begCompIdx + 1) << kExp) - 1);
                if (minIdx == MAX_T_ARRAYSIZE || miniBlocksQVal[minIdx >> miniKExp] == minQVal) {
                    uncertainResult = true;
                }
            }
#else
                uncertainResult = true;
#endif
        } else if (blocksQVal2D[begCompIdx] < minQVal) {
            const t_array_size firstBlockMinLoc = blocksLoc2D[begCompIdx];
            if (firstBlockMinLoc >= begIdx) {
                uncertainResult = false;
                minQVal = blocksQVal2D[begCompIdx];
                result = firstBlockMinLoc;
            } else {
#ifdef MINI_BLOCKS
                t_array_size minIdx = miniScanMinIdx(begIdx, ((begCompIdx + 1) << kExp) - 1);
                if (minIdx == MAX_T_ARRAYSIZE) {
                    return secondaryRMQ->rmq(begIdx, endIdx);
                }
                t_value tempQVal = miniBlocksQVal[minIdx >> miniKExp];
                if (tempQVal < minQVal) {
                    uncertainResult = false;
                    minQVal = tempQVal;
                    result = minIdx;
                } else if (tempQVal == minQVal) {
                        uncertainResult = true;
                }
#else
                return secondaryRMQ->rmq(begIdx, endIdx);
#endif
            }
        }
    }
    if (blocksQVal2D[endCompIdx] == minQVal) {
#ifdef MINI_BLOCKS
        if (blocksLoc2D[endCompIdx] <= endIdx) {
            return secondaryRMQ->rmq(begIdx, endIdx);
        } else {
            t_array_size minIdx = miniScanMinIdx(endCompIdx << kExp, endIdx);
            if (minIdx == MAX_T_ARRAYSIZE || miniBlocksQVal[minIdx >> miniKExp] == minQVal) {
                return secondaryRMQ->rmq(begIdx, endIdx);
            }
        }
#else
        return secondaryRMQ->rmq(begIdx, endIdx);
#endif
    } else if (blocksQVal2D[endCompIdx] < minQVal) {
        t_array_size lastBlockMinLoc = blocksLoc2D[endCompIdx];
        if (lastBlockMinLoc <= endIdx) {
            return lastBlockMinLoc;
        } else {
#ifdef MINI_BLOCKS
            t_array_size minIdx = miniScanMinIdx(endCompIdx << kExp, endIdx);
            if (minIdx == MAX_T_ARRAYSIZE) {
                return secondaryRMQ->rmq(begIdx, endIdx);
            }
            if (miniBlocksQVal[minIdx >> miniKExp] < minQVal) {
                return minIdx;
            } else if (miniBlocksQVal[minIdx >> miniKExp] == minQVal) {
                return secondaryRMQ->rmq(begIdx, endIdx);
            }
#else
            return secondaryRMQ->rmq(begIdx, endIdx);
#endif
        }
    }
    if (uncertainResult) {
        return secondaryRMQ->rmq(begIdx, endIdx);
    }

    return result;
}

template<typename t_qvalue, int max_qvalue> inline t_array_size BbSTqht<t_qvalue, max_qvalue>::miniScanMinIdx(const t_array_size &begIdx, const t_array_size &endIdx) {
    t_array_size result = -1;
    const t_array_size begMiniIdx = begIdx >> miniKExp;
    const t_array_size endMiniIdx = endIdx >> miniKExp;
    if (endMiniIdx == begMiniIdx) {
        const t_array_size firstMiniBlockMinLoc = (begMiniIdx << miniKExp) + miniBlocksLoc[begMiniIdx];
        if (begIdx <= firstMiniBlockMinLoc && firstMiniBlockMinLoc <= endIdx)
            return firstMiniBlockMinLoc;
        return MAX_T_ARRAYSIZE;
    }

    bool uncertainResult = false;
    t_qvalue minQVal = max_qvalue;
    if (endMiniIdx - begMiniIdx > 1) {
        for(t_array_size i = begMiniIdx + 1; i <= endMiniIdx - 1; i++) {
            t_qvalue tempQVal = miniBlocksQVal[i];
            if (tempQVal < minQVal) {
                uncertainResult = false;
                minQVal = tempQVal;
                result = (i << miniKExp) + miniBlocksLoc[i];
            } else if (tempQVal == minQVal) {
                uncertainResult = true;
            }
        }
    }

    t_value tempQVal = miniBlocksQVal[begMiniIdx];
    if (tempQVal <= minQVal && begIdx != (begMiniIdx + 1) << miniKExp) {
        const t_array_size firstMiniBlockMinLoc = (begMiniIdx << miniKExp) + miniBlocksLoc[begMiniIdx];
        minQVal = tempQVal;
        result = firstMiniBlockMinLoc;
        uncertainResult = (tempQVal == minQVal) || (firstMiniBlockMinLoc < begIdx);
    }

    tempQVal =  miniBlocksQVal[endMiniIdx];
    if (tempQVal == minQVal)
        return MAX_T_ARRAYSIZE;
    if (tempQVal < minQVal) {
        t_array_size lastBlockMinLoc = (endMiniIdx << miniKExp) + miniBlocksLoc[endMiniIdx];
        if (lastBlockMinLoc <= endIdx) {
            return lastBlockMinLoc;
        } else {
            return MAX_T_ARRAYSIZE;
        }
    }
    if (uncertainResult)
        return MAX_T_ARRAYSIZE;

    return result;
}

template<typename t_qvalue, int max_qvalue> void BbSTqht<t_qvalue, max_qvalue>::cleanup() {
    delete[] this->blocksLoc2D;
    delete[] this->blocksQVal2D;
#ifdef MINI_BLOCKS
    delete[] this->miniBlocksLoc;
    delete[] this->miniBlocksQVal;
#endif
}

template<typename t_qvalue, int max_qvalue> size_t BbSTqht<t_qvalue, max_qvalue>::memUsageInBytes() {
    D = 32 - __builtin_clz(blocksCount);
    const t_array_size blocksSize = blocksCount * D;
    size_t bytes = blocksSize * (sizeof(t_qvalue) + sizeof(t_array_size));
#ifdef MINI_BLOCKS
    bytes += miniBlocksCount * (sizeof(uint8_t) + sizeof(t_qvalue));
#endif
    return bytes;
}