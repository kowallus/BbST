#ifndef TESTDATA_H
#define TESTDATA_H

#include <vector>
#include "../common.h"

using namespace std;

void getRandomValues(vector<t_value> &data, const t_value modulo = 0);
void getPermutationOfRange(vector<t_value> &data);
void getPseudoMonotonicValues(vector<t_value> &data, t_value delta, bool decreasing);

void getRandomRangeQueries(vector<pair<t_array_size, t_array_size>> &queries, const t_array_size array_size, const t_array_size max_range_size);

vector<t_array_size> flattenQueries(const vector<pair<t_array_size, t_array_size>> &queriesPairs, const t_array_size queries_count);

void verify(const vector<t_value> &valuesArray, const vector<t_array_size> &queries, t_array_size *resultLoc);

void cleanCache();


#endif /* TESTDATA_H */

