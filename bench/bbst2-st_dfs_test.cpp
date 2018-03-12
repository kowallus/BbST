#include <iostream>
#include <algorithm>
#include "sdsl/suffix_arrays.hpp"
#include "sdsl/construct_lcp.hpp"

#include "../utils/testdata.h"
#include "../utils/timer.h"
#ifdef SDSL_REC
#include "../includes/sdsl/rmq_succinct_rec.hpp"
    #ifdef QUANTIZED
#include "../cbbstx.h"
    #else
#include "../bbstx.h"
    #endif
#else
#include "../bbst.h"
#endif
#include <unistd.h>
#include <omp.h>

#ifdef SDSL_REC
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
    #ifdef MINI_BLOCKS
string rmqName = "BbST2-sdsl-REC";
    #else
string rmqName = "BbST-sdsl-REC";
    #endif
#else
    #ifdef MINI_BLOCKS
string rmqName = "BbST2";
    #else
string rmqName = "BbST";
    #endif
#endif

using namespace sdsl;

typedef map<string, void (*)(cache_config&)> tMSFP;

using ll = long long;
using ival = pair<ll,ll>;

struct state {
    ll i,j,l;
    state(ll i, ll j, ll l) : i(i), j(j), l(l) { }
};

void construct_lcp(cache_config& test_config, string& test_file) {
    tMSFP lcp_function;
    lcp_function["bwt_based"] = &construct_lcp_bwt_based;
    lcp_function["bwt_based2"] = &construct_lcp_bwt_based2;
    lcp_function["PHI"] = &construct_lcp_PHI<8>;
    lcp_function["semi_extern_PHI"] = &construct_lcp_semi_extern_PHI;
    lcp_function["go"] = &construct_lcp_go;
    lcp_function["goPHI"] = &construct_lcp_goPHI;

    {
        cout << "Load text..." << endl;
        int_vector<8> text;
        load_vector_from_file(text,test_file,1);
        append_zero_symbol(text);
        store_to_cache(text, conf::KEY_TEXT, test_config);

        cout << "Construct Suffix Array..." << endl;
        int_vector<> sa(text.size(), 0,bits::hi(text.size())+1);
        algorithm::calculate_sa((const unsigned char*)text.data(), text.size(), sa);
        store_to_cache(sa,conf::KEY_SA,test_config);
    }

    {
        cout << "Construct LCP Array..." << endl;
        construct_lcp_PHI<8>(test_config);
    }

}

#ifdef SDSL_REC
    #ifdef QUANTIZED
size_t traverseSuffixTree(CBbSTx<uint8_t, 255>& rmq, int_vector<> lcp) {
    #else
size_t traverseSuffixTree(BbSTx& rmq, int_vector<> lcp) {
    #endif
#else
size_t traverseSuffixTree(BbST& rmq, int_vector<> lcp) {
#endif
    size_t num_queries = 0;
    size_t N = lcp.size();
    stack<state> s; s.emplace(0,N-1,0);

    while(!s.empty()) {
        state cur = s.top(); s.pop();
//                    if(cur.i != cur.j) cout << "Internal Node: (" << cur.i << "," << cur.j << ") - " << cur.l << endl;
        if(cur.i == cur.j) {
//                           cout << "Leaf: " << cur.i << endl;
            continue;
        }

        ll cur_i = cur.i;
        while(cur_i < cur.j) {
            num_queries++;
            size_t min_i = rmq.rmq(cur_i+1,cur.j);
//                           cout << cur_i << " " << cur.j << " " << min_i << endl;
            ll ii = cur_i; ll jj = cur.j-1;
            if(lcp[min_i] == cur.l && min_i < cur.j) jj = min_i-1;
            else if(lcp[min_i] != cur.l) jj = cur.j;
            if(ii+1 <= jj) {
                num_queries++;
                size_t l_idx = rmq.rmq(ii+1,jj);
                s.emplace(ii,jj,lcp[l_idx]);
            }
            else if(ii == jj) s.emplace(ii,ii,lcp[ii]);
            cur_i = jj+1;
        }

    }

    return num_queries;
}

int main(int argc, char**argv) {

#ifdef QUANTIZED
    rmqName = string("c") + rmqName;
#endif
    ChronoStopWatch timer;
    bool verbose = true;
    int kExp = 14;
#ifdef MINI_BLOCKS
    int miniKExp = 7;
#endif
    int noOfThreads = 1;
    int opt; // current option
#ifdef MINI_BLOCKS
    while ((opt = getopt(argc, argv, "k:l:t:q?")) != -1) {
#else
    while ((opt = getopt(argc, argv, "k:t:q?")) != -1) {
#endif
        switch (opt) {
            case 'q':
                verbose = false;
                break;
            case 'k':
                kExp = atoi(optarg);
                if (kExp < 1 || kExp > 24) {
                    fprintf(stderr, "%s: Expected 24>=k>=1\n", argv[0]);
                    fprintf(stderr, "try '%s -?' for more information\n", argv[0]);
                    exit(EXIT_FAILURE);
                }
                break;
#ifdef MINI_BLOCKS
            case 'l':
                miniKExp = atoi(optarg);
                if (miniKExp < 0 || miniKExp > 8) {
                    fprintf(stderr, "%s: Expected 8>=l>=0\n", argv[0]);
                    fprintf(stderr, "try '%s -?' for more information\n", argv[0]);
                    exit(EXIT_FAILURE);
                }
                break;
#endif
            case 't':
                noOfThreads = atoi(optarg);
                if (noOfThreads <= 0) {
                    fprintf(stderr, "%s: Expected noOfThreads >=1\n", argv[0]);
                    fprintf(stderr, "try '%s -?' for more information\n", argv[0]);
                    exit(EXIT_FAILURE);
                }
                break;
            case '?':
            default: /* '?' */
#ifdef MINI_BLOCKS
                fprintf(stderr, "Usage: %s [-k block size power of 2 exponent] [-l miniblock size power of 2 exponent] [-t noOfThreads] [-q]\n\n",
                        argv[0]);
                fprintf(stderr, "-k [24>=k>=1] \n-l [8>=l>=0] \n-t [noOfThreads>=1] \n-q quiet output (only parameters)\n\n");
#else
                fprintf(stderr, "Usage: %s [-k block size power of 2 exponent] [-t noOfThreads] [-v] [-q] n q\n\n",
                        argv[0]);
                fprintf(stderr, "-k [24>=k>=0] \n-t [noOfThreads>=1] \n-v verify results (extremely slow)\n-q quiet output (only parameters)\n\n");
#endif
                exit(EXIT_FAILURE);
        }
    }

    if (optind > (argc - 2)) {
        fprintf(stderr, "%s: Expected 2 argument after options (found %d)\n", argv[0], argc-optind);
        fprintf(stderr, "try '%s -?' for more information\n", argv[0]);

        exit(EXIT_FAILURE);
    }
#ifdef MINI_BLOCKS
    if (kExp <= miniKExp) {
        fprintf(stderr, "%s: k block size must be greater then miniblock size (k=%d, l=%d) \n", argv[0], kExp, miniKExp);

        exit(EXIT_FAILURE);
    }
#endif

    string test_file   = argv[optind];
    string temp_dir    = argv[optind+1];
    string test_id     = test_file.substr(test_file.find_last_of("/\\") + 1);

    cache_config test_config = cache_config(false, temp_dir, test_id);

    string lcp_file = cache_file_name(conf::KEY_LCP, test_config);
    int_vector<> lcp;
    if(!load_from_file(lcp,lcp_file)) {
        construct_lcp(test_config,test_file);
        load_from_file(lcp, lcp_file);
    }

    t_array_size N = lcp.size();
    vector<t_value> B(N);
    for(size_t i = 0; i < N; ++i) {
        B[i] = lcp[i];
    }
    if (verbose) cout << "Building "<< rmqName << "... " << std::endl;
    timer.startTimer();
#ifdef SDSL_REC
    #ifndef COUNT_2ND_RMQ
    CompetitorRMQ rmqIdx(&B[0], N);
    #else
    CompetitorRMQ sndRmqIdx(&B[0], N);
    RMQCounterDecorator rmqIdx(&sndRmqIdx);
    #endif
    #ifdef QUANTIZED
        #ifdef MINI_BLOCKS
    CBbSTx<uint8_t, 255> solver(B, kExp, miniKExp, &rmqIdx);
        #else
    CBbSTx<uint8_t, 255> solver(B, kExp, &rmqIdx);
        #endif
    #else
        #ifdef MINI_BLOCKS
    BbSTx solver(B, kExp, miniKExp, &rmqIdx);
        #else
    BbSTx solver(B, kExp, &rmqIdx);
        #endif
    #endif
#else
    #ifdef MINI_BLOCKS
    BbST solver(&B[0], N, kExp, miniKExp);
    #else
    BbST solver(&B[0], N, kExp);
    #endif
#endif
    timer.stopTimer();
    double buildTime = timer.getElapsedTime();
    if (verbose) cout << "Start Suffix-Tree Traversion ... " << std::endl;
    omp_set_num_threads(noOfThreads);

    timer.startTimer();
    size_t num_queries = traverseSuffixTree(solver,lcp);
    timer.stopTimer();

    double dfsTime = timer.getElapsedTime();
    double nanoqcoef = 1000000000.0 / num_queries;
    double queryTime = dfsTime * nanoqcoef;
    fstream fout(rmqName + "_st_dfs_res.txt", ios::out | ios::binary | ios::app);
#ifdef MINI_BLOCKS
    if (verbose) cout << "query time [ns]; N; q; dataset; size [KB]; k; miniK; noOfThreads; BbST build time [s]; dfs time [s]";
#else
    if (verbose) cout << "query time [ns]; N; q; dataset; size [KB]; k; noOfThreads; BbST build time [s]; dfs time [s]";
#endif
#ifdef COUNT_2ND_RMQ
    if (verbose) cout << "; 2nd-ary rmq query counter; log avaraged range";
#endif
    if (verbose) cout << std::endl;
#ifdef MINI_BLOCKS
    cout << queryTime << "\t" << N << "\t" << num_queries << "\t" << test_file
         << "\t" << (solver.memUsageInBytes() / 1000) << "\t" << (1 << kExp) << "\t" << (1 << miniKExp) << "\t" << noOfThreads
         << "\t" << buildTime << "\t" << dfsTime;
    fout << queryTime << "\t" << N << "\t" << num_queries << "\t" << test_file
         << "\t" << (solver.memUsageInBytes() / 1000) << "\t" << (1 << kExp) << "\t" << (1 << miniKExp) << "\t" << noOfThreads
         << "\t" << buildTime << "\t" << dfsTime;
#else
    cout << queryTime << "\t" << N << "\t" << num_queries << "\t" << test_file
         << "\t" << (solver.memUsageInBytes() / 1000) << "\t" << (1 << kExp) << "\t" << noOfThreads
         << "\t" << buildTime << "\t" << dfsTime;
    fout << queryTime << "\t" << N << "\t" << num_queries << "\t" << test_file
         << "\t" << (solver.memUsageInBytes() / 1000) << "\t" << (1 << kExp) << "\t" << noOfThreads
         << "\t" << buildTime << "\t" << dfsTime;
#endif
#ifdef COUNT_2ND_RMQ
    cout << "\t" << rmqIdx.getRMQCount() << "\t" << rmqIdx.getLogAvaregedM();
    fout << "\t" << rmqIdx.getRMQCount() << "\t" << rmqIdx.getLogAvaregedM();
#endif
    cout << std::endl;
    fout << std::endl;

    if (verbose) cout << "The end..." << std::endl;

    return 0;
}