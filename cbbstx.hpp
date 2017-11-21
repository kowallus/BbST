#include <algorithm>
#include <iostream>
#include <numeric>
#include "math.h"
#include "cbbstx.h"

#include <omp.h>

template<typename t_qvalue, int max_qvalue> CBbSTx<t_qvalue, max_qvalue>::~CBbSTx() {
    cleanup();
}

#ifdef MINI_BLOCKS
template<typename t_qvalue, int max_qvalue> CBbSTx<t_qvalue, max_qvalue>::CBbSTx(const vector<t_value> &valuesArray, int kExp, int miniKExp, RMQAPI* secondaryRMQ) {
    this->miniKExp = miniKExp;
    this->miniK = 1 << miniKExp;
#else
template<typename t_qvalue, int max_qvalue> CBbSTx<t_qvalue, max_qvalue>::CBbSTx(const vector<t_value> &valuesArray, int kExp, RMQAPI* secondaryRMQ) {
#endif
    this->kExp = kExp;
    this->k = 1 << kExp;
    this->secondaryRMQ = secondaryRMQ;
    prepareMinTables(valuesArray);
}

template<typename t_qvalue, int max_qvalue> void CBbSTx<t_qvalue, max_qvalue>::rmqBatch(const vector<t_array_size> &queries, t_array_size *resultLoc) {
    for (int i = 0; i < queries.size(); i = i + 2) {
        resultLoc[i / 2] = rmq(queries[i], queries[i + 1]);
    }
}

template<typename t_qvalue, int max_qvalue>  inline t_qvalue CBbSTx<t_qvalue, max_qvalue>::quantizeValue(const t_value value) {
//    double valMinMaxRatio = 1 - ((((double) maxMinVal - value)) / (((double) maxMinVal - minMinVal)));
//    double valMinMaxRatio = 1 - ((((double) maxMinVal - value)*(maxMinVal - value)) / (((double) maxMinVal - minMinVal)*(maxMinVal - minMinVal)));
//    double valMinMaxRatio = 1 - ((((double) maxMinVal - value)*(maxMinVal - value)*(maxMinVal - value)) / (((double) maxMinVal - minMinVal)*(maxMinVal - minMinVal)*(maxMinVal - minMinVal)));
//    double ratioExp = 5;
//    double valMinMaxRatio = 1 - (pow((double) maxMinVal - value, ratioExp) / (pow((double) maxMinVal - minMinVal, ratioExp)));
    double a = maxMinVal - value;
    a *= a * a;
    a *= a * a;
    a *= a * a;
    double b = maxMinVal - minMinVal;
    b *= b * b;
    b *= b * b;
    b *= b * b;
    double valMinMaxRatio = 1 - a / b;
    return (max_qvalue - 1) * valMinMaxRatio;
}

template<typename t_qvalue, int max_qvalue> void CBbSTx<t_qvalue, max_qvalue>::prepareMinTables(const vector<t_value> &valuesArray) {
#ifdef MINI_BLOCKS
    this->miniBlocksCount = (valuesArray.size() + miniK - 1) >> miniKExp;
    this->miniBlocksLoc = new uint8_t[miniBlocksCount];
    vector<t_value> miniBlocksVal(miniBlocksCount);
    this->miniBlocksQVal = new t_qvalue[miniBlocksCount];
    this->miniBlocksInBlock = k / miniK;
#endif
    this->blocksCount = (valuesArray.size() + k - 1) >> kExp;
    this->D = 32 - __builtin_clz(blocksCount);
    this->BD = 1 + ((D - 1) / 9);
    const t_array_size baseBlocksSize = blocksCount * BD;
    const t_array_size relativeBlocksSize = blocksCount * (D - BD);
    vector<t_value> tempBlocksVal(blocksCount);
    vector<t_array_size> tempBlocksLoc(blocksCount);
    this->baseBlocksValLoc2D = new uint8_t[baseBlocksSize * VALUE_AND_LOCATION_BYTES];
    this->blocksRelativeLoc2D = new uint8_t[relativeBlocksSize];
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
        t_value* valLocPtr = (t_value*) (baseBlocksValLoc2D + VALUE_AND_LOCATION_BYTES * i);
        *valLocPtr++ = tempBlocksVal[i] = *minPtr;
        *(t_array_size*) valLocPtr = tempBlocksLoc[i] = minPtr - &valuesArray[0];
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
    t_value* valLocPtr = (t_value*) (baseBlocksValLoc2D + VALUE_AND_LOCATION_BYTES * (blocksCount - 1));
    *valLocPtr++ = tempBlocksVal[blocksCount - 1] = *minPtr;
    *(t_array_size*) valLocPtr = tempBlocksLoc[blocksCount - 1] = minPtr - &valuesArray[0];
#ifdef MINI_BLOCKS
    minMinVal = miniBlocksVal[0];
    maxMinVal = miniBlocksVal[0];
    for (t_array_size i = 1; i < miniBlocksCount; i++)
        if (minMinVal > miniBlocksVal[i])
            minMinVal = miniBlocksVal[i];
        else if (maxMinVal < miniBlocksVal[i])
            maxMinVal = miniBlocksVal[i];
    for (t_array_size miniI = 0; miniI < miniBlocksCount; miniI++) {
        if (miniBlocksVal[miniI] >= maxMinVal)
            miniBlocksQVal[miniI] = (max_qvalue - 1);
        else
            miniBlocksQVal[miniI] = quantizeValue(miniBlocksVal[miniI]);
    }
#endif

    prepareBlocksSparseTable(tempBlocksVal, tempBlocksLoc);
}

template<typename t_qvalue, int max_qvalue> void CBbSTx<t_qvalue, max_qvalue>::prepareBlocksSparseTable(vector<t_value> &tempBlocksVal, vector<t_array_size> &tempBlocksLoc) {
    t_array_size eBDoffset = 0;
    t_array_size eRDoffset = 0;
    for(t_array_size e = 1, step = 1; e < D; ++e, step <<= 1) {
        uint8_t relFlag = e % 9;
        if (relFlag == 1)
            eBDoffset += blocksCount * VALUE_AND_LOCATION_BYTES;
        else
            eRDoffset += blocksCount;
        if (relFlag) {
            for (t_array_size i = 0; i < blocksCount; i++) {
                uint8_t tempLoc = 0;
                t_array_size minIdx = i;
                if (i + step < blocksCount && tempBlocksVal[i + step] < tempBlocksVal[i]) {
                    tempLoc = 1 << (relFlag - 1);
                    minIdx = i + step;
                    tempBlocksVal[i] = tempBlocksVal[i + step];
                    tempBlocksLoc[i] = tempBlocksLoc[i + step];
                }

                if (relFlag - 1)
                    tempLoc += blocksRelativeLoc2D[eRDoffset - blocksCount + minIdx];
                blocksRelativeLoc2D[eRDoffset + i] = tempLoc;
            }
        } else {
            for (t_array_size i = 0; i < blocksCount; i++) {
                if (i + step < blocksCount && tempBlocksVal[i + step] < tempBlocksVal[i]) {
                    tempBlocksVal[i] = tempBlocksVal[i + step];
                    tempBlocksLoc[i] = tempBlocksLoc[i + step];
                }
                t_value* valLocPtr = (t_value*) (baseBlocksValLoc2D + eBDoffset + VALUE_AND_LOCATION_BYTES * i);
                *valLocPtr++ = tempBlocksVal[i];
                *(t_array_size*) valLocPtr = tempBlocksLoc[i];
            }
        }
    }
}

template<typename t_qvalue, int max_qvalue> t_array_size CBbSTx<t_qvalue, max_qvalue>::rmq(const t_array_size &begIdx, const t_array_size &endIdx) {
    if (begIdx == endIdx) {
        return begIdx;
    }
    t_array_size result = MAX_T_ARRAYSIZE;
    const t_array_size begCompIdx = begIdx >> kExp;
    const t_array_size endCompIdx = endIdx >> kExp;
#ifdef START_FROM_NARROW_RANGES
    if (endCompIdx == begCompIdx) {
        const t_array_size firstBlockMinLoc = *((t_array_size *) (baseBlocksValLoc2D +
                                                                  begCompIdx * VALUE_AND_LOCATION_BYTES +
                                                                  LOCATION_OFFSET_IN_BYTES));
        if (begIdx <= firstBlockMinLoc && firstBlockMinLoc <= endIdx)
            return firstBlockMinLoc;
#ifndef MINI_BLOCKS
        return secondaryRMQ->rmq(begIdx, endIdx);
#else
        t_qvalue qMiniNotSmallerThen = max_qvalue;
        t_array_size minIdx = miniScanMinIdx(begIdx, endIdx, qMiniNotSmallerThen);
        if (minIdx == MAX_T_ARRAYSIZE) {
            return secondaryRMQ->rmq(begIdx, endIdx);
        }
        return minIdx;
#endif
    }
#endif
    const t_array_size kBlockCount = endCompIdx - begCompIdx; // actual kBlock count is +1
    const t_array_size e = kBlockCount?(31 - __builtin_clz(kBlockCount)):0;
    const t_array_size step = 1 << e;
    const uint8_t eBD = e / 9;
    const uint8_t relFlag = e % 9;
    t_array_size baseLocBegIdx = begCompIdx;
    t_array_size baseLocEndShiftIdx = endCompIdx - step + 1;
    if (relFlag) {
        const t_array_size eRDoffset = (e - eBD - 1) * blocksCount;
        baseLocBegIdx += ((t_array_size) blocksRelativeLoc2D[eRDoffset + baseLocBegIdx]) << (eBD * 9);
        baseLocEndShiftIdx += ((t_array_size) blocksRelativeLoc2D[eRDoffset + baseLocEndShiftIdx]) << (eBD * 9);
    }
    const t_array_size eBDoffset = eBD * blocksCount * VALUE_AND_LOCATION_BYTES;
    t_value *leftMinValLocPtr = (t_value *) (baseBlocksValLoc2D + eBDoffset + VALUE_AND_LOCATION_BYTES * baseLocBegIdx);
    t_value *rightMinValLocPtr = (t_value *) (baseBlocksValLoc2D + eBDoffset + VALUE_AND_LOCATION_BYTES * baseLocEndShiftIdx);
    bool minOnTheLeft = (*leftMinValLocPtr) <= (*rightMinValLocPtr);
    result = *((t_array_size*) ((minOnTheLeft?leftMinValLocPtr:rightMinValLocPtr) + 1));
    if (begIdx <= result && result <= endIdx)
        return result;
#ifndef MINI_BLOCKS
    return secondaryRMQ->rmq(begIdx, endIdx);
#else
    t_qvalue qMiniNotSmallerThen = max_qvalue;
    if (kBlockCount == 1) {
        if (*(t_array_size*) (leftMinValLocPtr + 1) >= begIdx) {
            t_qvalue minQVal = quantizeValue(*leftMinValLocPtr);
            t_array_size minIdx = miniScanMinIdx(endCompIdx << kExp, endIdx, qMiniNotSmallerThen);
            if (minIdx == MAX_T_ARRAYSIZE) {
                if (qMiniNotSmallerThen > minQVal)
                    return *(t_array_size *) (leftMinValLocPtr + 1);
                return secondaryRMQ->rmq(begIdx, endIdx);
            } else {
                if (miniBlocksQVal[minIdx >> miniKExp] < minQVal)
                    return minIdx;
                if (miniBlocksQVal[minIdx >> miniKExp] > minQVal)
                    return *(t_array_size *) (leftMinValLocPtr + 1);
                return secondaryRMQ->rmq(begIdx, endIdx);
            }
        } else if (*(t_array_size*) (rightMinValLocPtr + 1) <= endIdx) {
            t_qvalue minQVal = quantizeValue(*rightMinValLocPtr);
            t_array_size minIdx = miniScanMinIdx(begIdx, ((begCompIdx + 1) << kExp) - 1, qMiniNotSmallerThen);
            if (minIdx == MAX_T_ARRAYSIZE) {
                if (qMiniNotSmallerThen > minQVal)
                    return *(t_array_size *) (rightMinValLocPtr + 1);
                return secondaryRMQ->rmq(begIdx, endIdx);
            } else {
                if (miniBlocksQVal[minIdx >> miniKExp] < minQVal)
                    return minIdx;
                if (miniBlocksQVal[minIdx >> miniKExp] > minQVal)
                    return *(t_array_size *) (rightMinValLocPtr + 1);
                return secondaryRMQ->rmq(begIdx, endIdx);
            }
        }
    }
    if (kBlockCount <= 1) {
        t_array_size minIdx = miniScanMinIdx(begIdx, endIdx, qMiniNotSmallerThen);
        if (minIdx == MAX_T_ARRAYSIZE)
            return secondaryRMQ->rmq(begIdx, endIdx);
        return minIdx;
    }

    bool quantizedMode = false;
    bool uncertainResult = false;
    t_qvalue minQVal = max_qvalue;

    result = *((t_array_size*) (leftMinValLocPtr + 1));
    t_value minVal;
    if (result >= begIdx) {
        minVal = *leftMinValLocPtr;
    } else {
        t_array_size innerBaseLocBegIdx = begCompIdx + 1;
        t_array_size eBDoffset2;
        if (step == kBlockCount) {
            if (relFlag - 1) {
                const t_array_size eRDoffset = (e - (eBD - relFlag?0:1) - 2) * blocksCount;
                innerBaseLocBegIdx += ((t_array_size) blocksRelativeLoc2D[eRDoffset + innerBaseLocBegIdx]) << ((eBD - relFlag?0:1) * 9);
            }
            eBDoffset2 = (eBD - relFlag?0:1) * blocksCount * VALUE_AND_LOCATION_BYTES;
        } else {
            if (relFlag) {
                const t_array_size eRDoffset = (e - eBD - 1) * blocksCount;
                innerBaseLocBegIdx += ((t_array_size) blocksRelativeLoc2D[eRDoffset + innerBaseLocBegIdx]) << (eBD * 9);
            }
            eBDoffset2 = eBD * blocksCount * VALUE_AND_LOCATION_BYTES;
        }
        t_value *innerLeftMinValLocPtr = (t_value *) (baseBlocksValLoc2D + eBDoffset2 + VALUE_AND_LOCATION_BYTES * innerBaseLocBegIdx);

        minVal = *(innerLeftMinValLocPtr++);
        result = *((t_array_size*) innerLeftMinValLocPtr);
        const t_array_size minIdx = miniScanMinIdx(begIdx, ((begCompIdx + 1) << kExp) - 1, qMiniNotSmallerThen);
        minQVal = quantizeValue(minVal);
        if (minIdx == MAX_T_ARRAYSIZE) {
            if (qMiniNotSmallerThen <= minQVal) {
                quantizedMode = true;
                uncertainResult = true;
                minQVal = qMiniNotSmallerThen;
            }
        } else {
            t_value tempQVal = miniBlocksQVal[minIdx >> miniKExp];
            if (tempQVal < minQVal) {
                quantizedMode = true;
                minQVal = tempQVal;
                result = minIdx;
            } else if (tempQVal == minQVal) {
                quantizedMode = true;
                uncertainResult = true;
            }
        }
    }

    t_qvalue rightQVal;
    if (quantizedMode?(rightQVal = quantizeValue(*rightMinValLocPtr)) <= minQVal:(*rightMinValLocPtr < minVal)) {
        t_array_size tempLoc = *((t_array_size*) (rightMinValLocPtr + 1));
        if (tempLoc <= endIdx) {
            if (quantizedMode && (rightQVal == minQVal)) {
                uncertainResult = true;
            } else {
                return tempLoc;
            }
        }
        t_array_size innerBaseLocEndShiftIdx = endCompIdx - (step >> (step == kBlockCount));
        t_array_size eBDoffset2;
        if (step == kBlockCount) {
            if (relFlag - 1) {
                const t_array_size eRDoffset = (e - (eBD - relFlag?0:1) - 2) * blocksCount;
                innerBaseLocEndShiftIdx += ((t_array_size) blocksRelativeLoc2D[eRDoffset + innerBaseLocEndShiftIdx]) << ((eBD - relFlag?0:1) * 9);
            }
            eBDoffset2 = (eBD - relFlag?0:1) * blocksCount * VALUE_AND_LOCATION_BYTES;
        } else {
            if (relFlag) {
                const t_array_size eRDoffset = (e - eBD - 1) * blocksCount;
                innerBaseLocEndShiftIdx += ((t_array_size) blocksRelativeLoc2D[eRDoffset + innerBaseLocEndShiftIdx]) << (eBD * 9);
            }
            eBDoffset2 = eBD * blocksCount * VALUE_AND_LOCATION_BYTES;
        }
        t_value *innerRightMinValLocPtr = (t_value *) (baseBlocksValLoc2D + eBDoffset2 + VALUE_AND_LOCATION_BYTES * innerBaseLocEndShiftIdx);
        if (!quantizedMode) {
            if (*innerRightMinValLocPtr < minVal) {
                minVal = *innerRightMinValLocPtr;
                result = *((t_array_size *) (innerRightMinValLocPtr + 1));
            }
            rightQVal = quantizeValue(*rightMinValLocPtr);
        } else {
            t_qvalue tempQVal = quantizeValue(*innerRightMinValLocPtr);
            if (tempQVal < minQVal) {
                uncertainResult = false;
                minQVal = tempQVal;
                result = *((t_array_size *) (innerRightMinValLocPtr + 1));
            } else if (tempQVal == minQVal) {
                uncertainResult = true;
            }
        }
        minQVal = quantizedMode?minQVal:quantizeValue(minVal);
        quantizedMode = true;
    }
    if (quantizedMode && *((t_array_size*) (rightMinValLocPtr + 1)) > endIdx && rightQVal <= minQVal) {
        t_array_size minIdx = miniScanMinIdx(endCompIdx << kExp, endIdx, qMiniNotSmallerThen);
        if (minIdx == MAX_T_ARRAYSIZE) {
            if (qMiniNotSmallerThen <= minQVal)
                return secondaryRMQ->rmq(begIdx, endIdx);
        } else {
            if (miniBlocksQVal[minIdx >> miniKExp] < minQVal) {
                return minIdx;
            } else if (miniBlocksQVal[minIdx >> miniKExp] == minQVal) {
                return secondaryRMQ->rmq(begIdx, endIdx);
            }
        }
    }
    if (uncertainResult) {
        return secondaryRMQ->rmq(begIdx, endIdx);
    }
    return result;
#endif
}

template<typename t_qvalue, int max_qvalue> inline t_array_size CBbSTx<t_qvalue, max_qvalue>::miniScanMinIdx(const t_array_size &begIdx, const t_array_size &endIdx, t_qvalue &qNotSmallerThan) {
    t_array_size result = -1;
    const t_array_size begMiniIdx = begIdx >> miniKExp;
    const t_array_size endMiniIdx = endIdx >> miniKExp;
    if (endMiniIdx == begMiniIdx) {
        const t_array_size firstMiniBlockMinLoc = (begMiniIdx << miniKExp) + miniBlocksLoc[begMiniIdx];
        if (begIdx <= firstMiniBlockMinLoc && firstMiniBlockMinLoc <= endIdx)
            return firstMiniBlockMinLoc;
        qNotSmallerThan = miniBlocksQVal[begMiniIdx];
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
        uncertainResult = (tempQVal == minQVal) || (firstMiniBlockMinLoc < begIdx);
        minQVal = tempQVal;
        result = firstMiniBlockMinLoc;
    }

    tempQVal =  miniBlocksQVal[endMiniIdx];
    if (tempQVal == minQVal) {
        qNotSmallerThan = miniBlocksQVal[begMiniIdx]<tempQVal?:tempQVal;
        return MAX_T_ARRAYSIZE;
    }
    if (tempQVal < minQVal) {
        t_array_size lastBlockMinLoc = (endMiniIdx << miniKExp) + miniBlocksLoc[endMiniIdx];
        if (lastBlockMinLoc <= endIdx) {
            return lastBlockMinLoc;
        } else {
            qNotSmallerThan = tempQVal;
            return MAX_T_ARRAYSIZE;
        }
    }
    if (uncertainResult) {
        qNotSmallerThan = minQVal;
        return MAX_T_ARRAYSIZE;
    }

    return result;
}

template<typename t_qvalue, int max_qvalue> void CBbSTx<t_qvalue, max_qvalue>::cleanup() {
    delete[] this->blocksRelativeLoc2D;
    delete[] this->baseBlocksValLoc2D;
#ifdef MINI_BLOCKS
    delete[] this->miniBlocksLoc;
    delete[] this->miniBlocksQVal;
#endif
}

template<typename t_qvalue, int max_qvalue> size_t CBbSTx<t_qvalue, max_qvalue>::memUsageInBytes() {
    size_t bytes = blocksCount * (D - BD) * (sizeof(uint8_t));
    bytes += blocksCount * BD *(sizeof(t_value) + sizeof(t_array_size));
#ifdef MINI_BLOCKS
    bytes += miniBlocksCount * (sizeof(uint8_t) + sizeof(t_qvalue));
#endif
    return bytes + secondaryRMQ->memUsageInBytes();
}