#include "testdata.h"

std::mt19937 randgenerator;

void getRandomValues(vector<t_value> &data, t_value modulo) {
     randgenerator.seed(randgenerator.default_seed);
     for(long long int i = 0; i < data.size(); i++) {
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

void getRandomRangeQueries(vector<pair<t_array_size, t_array_size>> &queries, t_value array_size) {
    randgenerator.seed(randgenerator.default_seed);
    for(long long int i = 0; i < queries.size(); i++) {
        const unsigned int randA = randgenerator() % (array_size);
        const unsigned int randB = randgenerator() % (array_size);
        if (randA < randB) {
            queries[i].first =  randA;
            queries[i].second = randB;
        } else {
            queries[i].first =  randB;
            queries[i].second = randA;
        }
    }
}

vector<t_array_size> flattenQueries(vector<pair<t_array_size, t_array_size>> queriesPairs, t_array_size queries_count) {
    vector<t_array_size> result;
    result.reserve(queries_count * 2);
    for(long long int i = 0; i < queries_count; i++) {
        result.push_back(queriesPairs[i % queriesPairs.size()].first);
        result.push_back(queriesPairs[i % queriesPairs.size()].second);
    }
    return result;
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
        