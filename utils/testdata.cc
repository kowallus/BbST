#include "testdata.h"

std::mt19937 randgenerator;

void getRandomValues(vector<t_value> &data, const t_value modulo) {
     randgenerator.seed(randgenerator.default_seed);
     for(t_array_size i = 0; i < data.size(); i++) {
         data[i] = randgenerator();
         if (modulo > 0)
             data[i] %= modulo;
     }
}

void getPermutationOfRange(vector<t_value> &data) {
    t_value n = data.size();
    randgenerator.seed(randgenerator.default_seed);

    for (t_value i= 0; i < n; i++ )
        data[i] = i;

    for ( t_value i = 0; i < (n/2); i++ )
    {
        t_value a = randgenerator() % n;
        t_value b = randgenerator() % n;
        t_value temp = data[a];
        data[a] = data[b];
        data[b] = temp;
    }
}

void getPseudoMonotonicValues(vector<t_value> &data, t_value delta, bool decreasing) {
    t_value n = data.size();
    randgenerator.seed(randgenerator.default_seed);
    for(t_array_size i = 0; i < data.size(); i++) {
        data[i] = i + (randgenerator() % (2 * delta + 1)) - delta;
    }
    if (decreasing)
        for(t_array_size i = 0; i < data.size(); i++) {
            data[i] = n - data[i];
        }
}

void getRandomRangeQueries(vector<pair<t_array_size, t_array_size>> &queries, const t_array_size array_size, const t_array_size max_range_size) {
    randgenerator.seed(randgenerator.default_seed);
    for(long long int i = 0; i < queries.size(); i++) {
        const t_array_size randA = randgenerator() % (array_size);
        t_array_size maxB = randA + max_range_size;
        if (maxB >= array_size) {
            maxB = array_size - 1;
        }
        t_array_size minB = 0;
        if (max_range_size < randA) {
            minB = randA - max_range_size;
        }
        const t_array_size randB = minB + randgenerator() % (maxB - minB + 1);
        if (randA < randB) {
            queries[i].first =  randA;
            queries[i].second = randB;
        } else {
            queries[i].first =  randB;
            queries[i].second = randA;
        }
    }
}

vector<t_array_size> flattenQueries(const vector<pair<t_array_size, t_array_size>> &queriesPairs, const t_array_size queries_count) {
    vector<t_array_size> result;
    result.reserve(queries_count * 2);
    for(long long int i = 0; i < queries_count; i++) {
        result.push_back(queriesPairs[i % queriesPairs.size()].first);
        result.push_back(queriesPairs[i % queriesPairs.size()].second);
    }
    return result;
}

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
                 << verifyLoc[i] << " is " << resultLoc[i]
                 << " (vals:" << valuesArray[verifyLoc[i]] << " and " << valuesArray[resultLoc[i]] << ")" << std::endl;
        }
    }
}

void cleanCache() {
    const int ccSize = 10*1024*1024; // Allocate 10M. Set much larger then L2
    const int ccRepeats = 0xff;
    char *cc = (char *)malloc(ccSize);
    for (int cci = 0; cci < ccRepeats; cci++)
        for (int ccj = 0; ccj < ccSize; ccj++)
            cc[ccj] = cci+ccj;
    free(cc);
}
        