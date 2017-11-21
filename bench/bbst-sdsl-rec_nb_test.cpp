#include <iostream>
#include <algorithm>

#include "../utils/testdata.h"
#include "../utils/timer.h"
#include "../includes/sdsl/rmq_succinct_rec.hpp"
#ifdef QUANTIZED
#include "../cbbstx.h"
#else
#include "../bbstx.h"
#endif
#include <unistd.h>
#include <omp.h>

class CompetitorRMQ: public RMQAPI {
private:
    typedef sdsl::rmq_succinct_rec<> rmqStruct;
    rmqStruct *rmqImpl;
public:

    CompetitorRMQ(const t_value* valuesArray, const t_array_size n) {
        typedef sdsl::int_vector<32> int_vector;
        int_vector intVector(n);
        for (t_array_size i = 0; i < n; i++) {
            intVector[i] = (uint32_t)((int64_t) valuesArray[i]) - INT32_MIN;
        }
        rmqImpl = new rmqStruct(&intVector);
    }

    ~CompetitorRMQ() {
        delete rmqImpl;
    }

    t_array_size rmq(const t_array_size &begIdx, const t_array_size &endIdx) {
        return rmqImpl->operator()(begIdx, endIdx);
    }

    size_t memUsageInBytes() {
        return sdsl::size_in_bytes(*rmqImpl);
    }
};

string rmqName = "BbST-sdsl-REC";

int main(int argc, char**argv) {

#ifdef QUANTIZED
    rmqName = string("c") + rmqName;
    fstream fout(rmqName + "_nb_res.txt", ios::out | ios::binary | ios::app);
#else
    fstream fout(rmqName + "_nb_res.txt", ios::out | ios::binary | ios::app);
#endif
    ChronoStopWatch timer;
    bool verbose = true;
    bool verification = false;
    int kExp = 14;
    int noOfThreads = 1;
    int opt; // current option
    int repeats = 1;
    t_array_size max_range = 0;

    while ((opt = getopt(argc, argv, "k:t:r:m:vq?")) != -1) {
        switch (opt) {
            case 'q':
                verbose = false;
                break;
            case 'v':
                verification = true;
                break;
            case 'k':
                kExp = atoi(optarg);
                if (kExp < 0 || kExp > 24) {
                    fprintf(stderr, "%s: Expected 24>=k>=0\n", argv[0]);
                    fprintf(stderr, "try '%s -?' for more information\n", argv[0]);
                    exit(EXIT_FAILURE);
                }
                break;
            case 't':
                noOfThreads = atoi(optarg);
                if (noOfThreads <= 0) {
                    fprintf(stderr, "%s: Expected noOfThreads >=1\n", argv[0]);
                    fprintf(stderr, "try '%s -?' for more information\n", argv[0]);
                    exit(EXIT_FAILURE);
                }
                break;
            case 'r':
                repeats = atoi(optarg);
                if (repeats <= 0) {
                    fprintf(stderr, "%s: Expected number of repeats >=1\n", argv[0]);
                    fprintf(stderr, "try '%s -?' for more information\n", argv[0]);
                    exit(EXIT_FAILURE);
                }
                break;
            case 'm':
                max_range = atoi(optarg);
                if (max_range <= 0) {
                    fprintf(stderr, "%s: Expected maximum size of a range>=1\n", argv[0]);
                    fprintf(stderr, "try '%s -?' for more information\n", argv[0]);
                    exit(EXIT_FAILURE);
                }
                break;
            case '?':
            default: /* '?' */
                fprintf(stderr, "Usage: %s [-k block size power of 2 exponent] [-t noOfThreads] [-v] [-q] n q\n\n",
                        argv[0]);
                fprintf(stderr, "-k [24>=k>=0] \n-t [noOfThreads>=1] \n-v verify results (extremely slow)\n-q quiet output (only parameters)\n\n");
                exit(EXIT_FAILURE);
        }
    }

    if (optind > (argc - 2)) {
        fprintf(stderr, "%s: Expected 2 arguments after options (found %d)\n", argv[0], argc-optind);
        fprintf(stderr, "try '%s -?' for more information\n", argv[0]);

        exit(EXIT_FAILURE);
    }

    t_array_size n = atoi(argv[optind++]);
    t_array_size q = atoi(argv[optind]);
    if (max_range == 0) {
        max_range = n;
    }

    if (verbose) cout << "Generation of values..." << std::endl;
    vector<t_value> valuesArray(n);
#ifdef RANDOM_DATA
    getRandomValues(valuesArray, MAX_T_VALUE / 4);
#else
    getPermutationOfRange(valuesArray);
#endif

    if (verbose) cout << "Generation of queries..." << std::endl;
    vector<pair<t_array_size, t_array_size>> queriesPairs(q);

    getRandomRangeQueries(queriesPairs, n, max_range);

    vector<t_array_size> queries = flattenQueries(queriesPairs, q);
    t_array_size* resultLoc = new t_array_size[queries.size() / 2];

    if (verbose) cout << "Building "<< rmqName << "... " << std::endl;
    timer.startTimer();
    CompetitorRMQ rmqIdx(&valuesArray[0], valuesArray.size());
#ifdef QUANTIZED
    CBbSTx<uint8_t, 255> solver(valuesArray, kExp, &rmqIdx);
#else
    BbSTx solver(valuesArray, kExp, &rmqIdx);
#endif

    timer.stopTimer();
    double buildTime = timer.getElapsedTime();
    if (verbose) cout << "Solving... " << std::endl;

    omp_set_num_threads(noOfThreads);
    vector<double> times;
    for(int i = 0; i < repeats; i++) {
        cleanCache();
        timer.startTimer();
        solver.rmqBatch(queries, resultLoc);
        timer.stopTimer();
        times.push_back(timer.getElapsedTime());
    }
    std::sort(times.begin(), times.end());
    double nanoqcoef = 1000000000.0 / q;
    double maxQueryTime = times[repeats - 1] * nanoqcoef ;
    double medianQueryTime = times[times.size()/2] * nanoqcoef;
    double minQueryTime = times[0] * nanoqcoef;
    if (verbose) cout << "query time [ns]; n; q; m; size [KB]; k; noOfThreads; BbST build time [s]; max/min time [ns]" << std::endl;
    cout << medianQueryTime << "\t" << valuesArray.size() << "\t" << (queries.size() / 2) << "\t" << max_range
         << "\t" << (solver.memUsageInBytes() / 1000) << "\t" << (1 << kExp) << "\t" << noOfThreads
         << "\t" << buildTime << "\t" << maxQueryTime << "\t" << minQueryTime << "\t" << std::endl;
    fout << medianQueryTime << "\t" << valuesArray.size() << "\t" << (queries.size() / 2) << "\t" << max_range <<
         "\t" << (solver.memUsageInBytes() / 1000) << "\t" << (1 << kExp) << "\t" << noOfThreads <<
         "\t" << buildTime << "\t" << maxQueryTime << "\t" << minQueryTime << "\t" << std::endl;
    if (verification) verify(valuesArray, queries, resultLoc);

    if (verbose) cout << "The end..." << std::endl;
    return 0;
}