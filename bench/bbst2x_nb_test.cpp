#include <iostream>
#include <algorithm>

#include "../utils/testdata.h"
#include "../utils/timer.h"
#include "../bbstx.h"

#include <unistd.h>
#include <omp.h>

int main(int argc, char**argv) {

    fstream fout("BbST2x_nb_res.txt", ios::out | ios::binary | ios::app);

    ChronoStopWatch timer;
    bool verbose = true;
    bool verification = false;
    int kExp = 14;
    int miniKExp = 7;
    int noOfThreads = 1;
    int opt; // current option
    int repeats = 1;
    t_array_size max_range = 0;

    while ((opt = getopt(argc, argv, "k:l:t:r:m:vq?")) != -1) {
        switch (opt) {
            case 'q':
                verbose = false;
                break;
            case 'v':
                verification = true;
                break;
            case 'k':
                kExp = atoi(optarg);
                if (kExp < 1 || kExp > 24) {
                    fprintf(stderr, "%s: Expected 24>=k>=1\n", argv[0]);
                    fprintf(stderr, "try '%s -?' for more information\n", argv[0]);
                    exit(EXIT_FAILURE);
                }
                break;
            case 'l':
                miniKExp = atoi(optarg);
                if (miniKExp < 0 || miniKExp > 8) {
                    fprintf(stderr, "%s: Expected 8>=l>=0\n", argv[0]);
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
                fprintf(stderr, "Usage: %s [-k block size power of 2 exponent] [-l miniblock size power of 2 exponent] [-t noOfThreads] [-v] [-q] n q\n\n",
                        argv[0]);
                fprintf(stderr, "-k [24>=k>=1] \n-l [8>=l>=0] \n-t [noOfThreads>=1] \n-v verify results (extremely slow)\n-q quiet output (only parameters)\n\n");
                exit(EXIT_FAILURE);
        }
    }

    if (optind > (argc - 2)) {
        fprintf(stderr, "%s: Expected 2 arguments after options (found %d)\n", argv[0], argc-optind);
        fprintf(stderr, "try '%s -?' for more information\n", argv[0]);

        exit(EXIT_FAILURE);
    }

    if (kExp <= miniKExp) {
        fprintf(stderr, "%s: k block size must be greater then miniblock size (k=%d, l=%d) \n", argv[0], kExp, miniKExp);

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

    if (verbose) cout << "Building sBbST2... " << std::endl;
    timer.startTimer();
    RMQCounter rmqCounter;
    BbSTx solver(valuesArray, kExp, miniKExp, &rmqCounter);
    timer.stopTimer();
    double buildTime = timer.getElapsedTime();
    if (verbose) cout << "Solving... " << std::endl;

    omp_set_num_threads(noOfThreads);
    vector<double> times;
    for(int i = 0; i < repeats; i++) {
        cleanCache();
        rmqCounter.resetCounter();
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
    double successRate = 100 - (100.0 * ((double) rmqCounter.getRMQCount()) / q);
    if (verbose) cout << "query time [ns]; successRate [%]; n; q; m; size [KB]; k; miniK; noOfThreads; BbST build time [s]; max/min time [ns]" << std::endl;
    cout << medianQueryTime << "\t" << successRate << "\t" << valuesArray.size() << "\t" << (queries.size() / 2) << "\t" << max_range
         << "\t" << (solver.memUsageInBytes() / 1000) << "\t" << (1 << kExp) << "\t" << (1 << miniKExp) << "\t" << noOfThreads
         << "\t" << buildTime << "\t" << maxQueryTime << "\t" << minQueryTime << "\t" << std::endl;
    fout << medianQueryTime << "\t" << successRate << "\t" << valuesArray.size() << "\t" << (queries.size() / 2) << "\t" << max_range <<
         "\t" << (solver.memUsageInBytes() / 1000) << "\t" << (1 << kExp) << "\t" << (1 << miniKExp) << "\t" << noOfThreads <<
         "\t" << buildTime << "\t" << maxQueryTime << "\t" << minQueryTime << "\t" << std::endl;
    if (verification) verify(valuesArray, queries, resultLoc);

    if (verbose) cout << "The end..." << std::endl;
    return 0;
}