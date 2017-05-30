#include <iostream>
#include <algorithm>

#include "../utils/testdata.h"
#include "../utils/timer.h"
#include "../bbstcon.h"

#include <unistd.h>

int main(int argc, char**argv) {

    fstream fout("bbstcon_res.txt", ios::out | ios::binary | ios::app);

    ChronoStopWatch timer;
    bool verbose = true;
    bool verification = false;
    int noOfThreads = 1;
    sortingAlg_enum sortingAlg = kxradixsort;
    int opt; // current option
    int repeats = 1;

    while ((opt = getopt(argc, argv, "t:s:r:vq?")) != -1) {
        switch (opt) {
            case 'q':
                verbose = false;
                break;
            case 'v':
                verification = true;
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
            case 's':
                sortingAlg = (sortingAlg_enum) optarg[0];
                if (sortingAlg != csort && sortingAlg != ompparallelsort && sortingAlg != pssparallelsort &&
                        sortingAlg != stdsort && sortingAlg != kxradixsort) {
                    fprintf(stderr, "%s: Unknown sorting algorithm option.\n", argv[0]);
                    fprintf(stderr, "try '%s -?' for more information\n", argv[0]);
                    exit(EXIT_FAILURE);
                }
                break;
            case '?':
            default: /* '?' */
                fprintf(stderr, "Usage: %s [-t noOfThreads] [-s sortingAlgorithm] [-v] [-q] n q\n\n",
                        argv[0]);
                fprintf(stderr, "-t [noOfThreads>=1] \n-s [q-quicksort;s-stdsort;r-kxradixsort;i-psspparallelsort;p-ompparallelsort;\n-v verify results (extremely slow)\n-q quiet output (only parameters)\n\n");
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

    if (verbose) cout << "Generation of values..." << std::endl;
    vector<t_value> valuesArray(n);
#ifdef RANDOM_DATA
    getRandomValues(valuesArray, MAX_T_VALUE / 4);
#else
    getPermutationOfRange(valuesArray);
#endif

    if (verbose) cout << "Generation of queries..." << std::endl;
    vector<pair<t_array_size, t_array_size>> queriesPairs(q);
    getRandomRangeQueries(queriesPairs, n);
    vector<t_array_size> queries = flattenQueries(queriesPairs, q);
    t_array_size* resultLoc = new t_array_size[queries.size() / 2];

    BbSTcon solver(valuesArray, queries, resultLoc, sortingAlg, noOfThreads);
    if (verbose) cout << "Solving... ";

    vector<double> times;
    for(int i = 0; i < repeats; i++) {
        if (i > 0) {
            cleanCache();
        }
        timer.startTimer();
        solver.solve();
        timer.stopTimer();
        times.push_back(timer.getElapsedTime());
    }
    std::sort(times.begin(), times.end());
    double medianTime = times[times.size()/2];
    if (verbose) cout << "elapsed time [s] n q sorting noOfThreads K maxTime minTime:\t";
    cout << medianTime << "\t" << valuesArray.size() << "\t" << (queries.size() / 2) <<
        "\t" << (char) sortingAlg << "\t" << noOfThreads << "\t" << K <<
        "\t" << times[repeats - 1] << "\t" << times[0] << "\t" << std::endl;
    fout << medianTime << "\t" << valuesArray.size() << "\t" << (queries.size() / 2) <<
        "\t" << (char) sortingAlg << "\t" << noOfThreads << "\t" << K <<
        "\t" << times[repeats - 1] << "\t" << times[0] << "\t" << std::endl;
    if (verification) solver.verify();

    if (verbose) cout << "The end..." << std::endl;
    return 0;
}

