#ifndef TESTDATA_H
#define TESTDATA_H

#include <vector>
#include "../common.h"

using namespace std;

void getRandomValues(vector<t_value> &data, t_value modulo = 0);
void getPermutationOfRange(vector<t_value> &data);

void getRandomRangeQueries(vector<pair<t_array_size, t_array_size>> &queries, t_array_size array_size, t_array_size max_range_size);

vector<t_array_size> flattenQueries(vector<pair<t_array_size, t_array_size>> vector, t_array_size queries_count);

void cleanCache();


#endif /* TESTDATA_H */

