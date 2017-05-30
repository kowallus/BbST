#include <algorithm>
#include <iostream>
#include <numeric>
#include "bbstcon.h"

#include <omp.h>
#include <parallel/algorithm>
#include "utils/kxsort.h"
#include "utils/parallel_stable_sort.h"

BbSTcon::BbSTcon(std::vector<t_value, std::allocator<t_value>> &valuesArray,
                     std::vector<t_array_size, std::allocator<t_array_size>> &queries,
                     t_array_size* resultLoc,
                     sortingAlg_enum sortingAlg, int noOfThreads) {
    this->valuesArray = valuesArray;
    this->queries = queries;
    this->resultLoc = resultLoc;
    this->sortingAlg = sortingAlg;
    omp_set_num_threads(noOfThreads);
}

void BbSTcon::solve() {
    getUniqueBoundsSorted();
    getContractedMins();
    getBlocksMins();
    for(int i = 0; i < queries.size(); i = i + 2) {
        if (queries[i] == queries[i+1]) {
            this->resultLoc[i/2] = queries[i];
            continue;
        }
        this->resultLoc[i/2] = getRangeMinLoc(queries2ContractedIdx[i], queries2ContractedIdx[i+1]);
    }
    cleanup();
/**/
}

void BbSTcon::verify() {
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

struct RadixTraitsBounds {
    static const int nBytes = sizeof(t_array_size);
    int kth_byte(const t_array_size_2x &x, int k) {
        return *((t_array_size*) &x) >> (k * 8) & 0xFF;
    }
    bool compare(const t_array_size_2x &a, const t_array_size_2x &b) {
        return *((t_array_size*) &a) < *((t_array_size*) &b);
    }
};

void BbSTcon::getUniqueBoundsSorted() {
    bounds = new t_array_size_2x[queries.size()];
    for (t_array_size i = 0; i < queries.size(); i++) {
        *((t_array_size *) (&bounds[i])) = queries[i];
        *(((t_array_size *) (&bounds[i])) + 1) = i;
    }
    switch (sortingAlg) {
        case stdsort : std::sort(bounds, bounds + queries.size(), [](const t_array_size_2x& a, const t_array_size_2x& b) -> bool { return *((t_array_size*) &a) < *((t_array_size*) &b); });
            break;
        case csort : qsort((void*) &bounds[0], queries.size(), sizeof(unsigned long long int), [](const void* a, const void* b) -> int { return *((int32_t*) a) - *((t_array_size*) b); });
            break;
        case kxradixsort : kx::radix_sort(bounds, bounds + queries.size(), RadixTraitsBounds());
            break;
        case pssparallelsort : pss::parallel_stable_sort(bounds, bounds + queries.size(), [](const t_array_size_2x& a, const t_array_size_2x& b) -> bool { return *((t_array_size*) &a) < *((t_array_size*) &b); });
            break;
        case ompparallelsort  : __gnu_parallel::sort(bounds, bounds + queries.size(), [](const t_array_size_2x& a, const t_array_size_2x& b) -> bool { return *((t_array_size*) &a) < *((t_array_size*) &b); });
//        case ompparallelsort  : __gnu_parallel::sort(bounds, bounds + queries.size(), [](const t_array_size_2x& a, const t_array_size_2x& b) -> bool { return *((t_array_size*) &a) < *((t_array_size*) &b); }, __gnu_parallel::multiway_mergesort_tag());
            break;
    }
    queries2ContractedIdx = new t_array_size[queries.size()];
    for(t_array_size i = 0; i < queries.size(); i++) {
        queries2ContractedIdx[*(((t_array_size*) &bounds[i]) + 1)] = i;
    }
}

void BbSTcon::getContractedMins() {
    contractedVal = new t_value[queries.size() - 1];
    contractedLoc = new t_array_size[queries.size() - 1];

    #pragma omp parallel for
    for(int i = 1;  i < queries.size(); i++) {
        t_array_size minIdx = *((t_array_size*) (bounds + i - 1));
        const t_array_size endIdx = *((t_array_size*) (bounds + i));

        const int *const minPtr = std::min_element(&valuesArray[minIdx], &valuesArray[endIdx + 1]);
        contractedVal[i-1] = *minPtr;
        contractedLoc[i-1] = minPtr - &valuesArray[0];

//        t_value minVal = valuesArray[minIdx];
//        for (t_array_size idx = minIdx + 1; idx <= endIdx; ++idx)
//            if (valuesArray[idx] < minVal) {
//                minVal = valuesArray[idx];
//                minIdx = idx;
//            }
//        contractedVal[i-1] = minVal;
//        contractedLoc[i-1] = minIdx;

    }
}

void BbSTcon::getBlocksMins() {
    const t_array_size blocksCount = ((queries.size() - 1 + K - 1)/ K);
    D = 32 - __builtin_clz(blocksCount);
    const t_array_size blocksSize = blocksCount * D;
    blocksVal2D = new t_value[blocksSize];
    blocksLoc2D = new t_array_size[blocksSize];
    t_array_size i = 0;
    t_value* ptr;
    t_value* endPtr = contractedVal + (queries.size() - 1);
    for(ptr = contractedVal; ptr < endPtr - K; ptr = ptr + K) {
        auto minPtr = std::min_element(ptr, ptr + K);
        blocksVal2D[i * D] = *minPtr;
        blocksLoc2D[i++ * D] = contractedLoc[minPtr - contractedVal];
    }
    auto minPtr = std::min_element(ptr, endPtr);
    blocksVal2D[i * D] = *minPtr;
    blocksLoc2D[i * D] = contractedLoc[minPtr - contractedVal];
    for(t_array_size e = 1, step = 1; e < D; ++e, step <<= 1) {
        for (i = 0; i < blocksCount; i++) {
            t_array_size minIdx = i;
            if (i + step < blocksCount && blocksVal2D[(i + step) * D + e - 1] < blocksVal2D[i * D + e - 1]) {
                minIdx = i + step;
            }
            blocksVal2D[i * D + e] = blocksVal2D[minIdx * D + e - 1];
            blocksLoc2D[i * D + e] = blocksLoc2D[minIdx * D + e - 1];
        }
    }
}

t_array_size BbSTcon::getRangeMinLoc(const t_array_size &begContIdx, const t_array_size &endContIdx) {
    t_array_size result = -1;
    const t_array_size begCompIdx = begContIdx / K;
    const t_array_size endCompIdx = endContIdx / K;
    if (endCompIdx - begCompIdx <= 1) {
        return contractedLoc[scanContractedMinIdx(begContIdx, endContIdx)];
    }

    t_array_size kBlockCount = endCompIdx - begCompIdx - 1;
    t_array_size e = 31 - __builtin_clz(kBlockCount);
    t_array_size step = 1 << e;
    t_value minVal = blocksVal2D[(begCompIdx + 1) * D + e];
    result = blocksLoc2D[(begCompIdx + 1) * D + e];
    t_array_size endShiftCompIdx = endCompIdx - step;
    if (endShiftCompIdx != begCompIdx + 1) {
        t_value temp = blocksVal2D[(endShiftCompIdx) * D + e];
        if (temp < minVal) {
            minVal = temp;
            result = blocksLoc2D[(endShiftCompIdx) * D + e];
        }
    }

    if (blocksVal2D[begCompIdx * D + 0] <= minVal && begContIdx != (begCompIdx + 1) * K) {
        t_array_size contMinIdx = scanContractedMinIdx(begContIdx, (begCompIdx + 1) * K);
        if (contractedVal[contMinIdx] <= minVal) {
            minVal = contractedVal[contMinIdx];
            result = contractedLoc[contMinIdx];
        }
    }
    if (blocksVal2D[endCompIdx * D + 0] < minVal && endCompIdx * K != endContIdx) {
       t_array_size contMinIdx = scanContractedMinIdx(endCompIdx * K, endContIdx);
        if (contractedVal[contMinIdx] < minVal) {
//            minVal = contractedVal[contMinIdx];
            result = contractedLoc[contMinIdx];
        }
    }
    return result;
}

t_array_size BbSTcon::scanContractedMinIdx(const t_array_size &begContIdx, const t_array_size &endContIdx) {
    t_array_size minValIdx = begContIdx;
    for(t_array_size i = begContIdx + 1; i < endContIdx; i++) {
        if (contractedVal[i] < contractedVal[minValIdx]) {
            minValIdx = i;
        }
    }
    return minValIdx;
}

void BbSTcon::cleanup() {
    delete[] this->bounds;
    delete[] this->blocksLoc2D;
    delete[] this->blocksVal2D;
    delete[] this->contractedLoc;
    delete[] this->contractedVal;
    delete[] this->queries2ContractedIdx;
}
