#ifndef COMMON_H
#define COMMON_H

#include <limits.h>
#include <vector>
#include "bits/stdc++.h"

#define RANDOM_DATA

#define MAX_T_VALUE INT32_MAX
#define MAX_T_ARRAYSIZE UINT32_MAX
typedef int32_t t_value;

typedef unsigned int t_array_size;
typedef long long int t_array_size_2x;

using namespace std;

void verify(const vector<t_value> &valuesArray, const vector<t_array_size> &queries, t_array_size *resultLoc);

#endif /* COMMON_H */

