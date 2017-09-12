cmake_minimum_required(VERSION 3.7)
project(BbST)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2 -fopenmp")

set(COMMON_SOURCE_FILES
        common.h
        utils/testdata.cc
        utils/testdata.h
        utils/timer.cpp
        utils/timer.h)

set(BBSTCON_SOURCE_FILES
        ${COMMON_SOURCE_FILES}
        bench/bbstcon_test.cpp
        bbstcon.cpp
        bbstcon.h
        utils/kxsort.h
        utils/parallel_stable_sort.h
        utils/pss_common.h)

set(BBST_SOURCE_FILES
        ${COMMON_SOURCE_FILES}
        bbst.cpp
        bbst.h)

set(BBSTHT_SOURCE_FILES
        ${COMMON_SOURCE_FILES}
        bbstht.cpp
        bbstht.h)

set(CBBSTHT_SOURCE_FILES
        ${COMMON_SOURCE_FILES}
        cbbstht.h
        cbbstht.hpp)

set(HFERRADA_RMQ_SOURCE_FILES
        RMQRMM64.cpp
        includes/RMQRMM64.h
        includes/Basic_rmq.cpp
        includes/Basic_rmq.h)

add_executable(bbstcon ${BBSTCON_SOURCE_FILES})
add_executable(bbst bench/bbst_test.cpp ${BBST_SOURCE_FILES})
add_executable(bbst_nb bench/bbst_nb_test.cpp ${BBST_SOURCE_FILES})
add_executable(bbst2_nb bench/bbst2_nb_test.cpp ${BBST_SOURCE_FILES})
target_compile_definitions(bbst2_nb PUBLIC "-DMINI_BLOCKS")
add_executable(bbst2_nb_wc bench/bbst2_nb_test.cpp ${BBST_SOURCE_FILES})
target_compile_definitions(bbst2_nb_wc PUBLIC "-DMINI_BLOCKS -DWORST_CASE")
add_executable(bbstht_nb bench/bbstht_nb_test.cpp ${BBSTHT_SOURCE_FILES})
add_executable(bbst2ht_nb bench/bbst2ht_nb_test.cpp ${BBSTHT_SOURCE_FILES})
target_compile_definitions(bbst2ht_nb PUBLIC "-DMINI_BLOCKS")
add_executable(cbbstht_nb bench/cbbstht_nb_test.cpp ${CBBSTHT_SOURCE_FILES})
add_executable(cbbst2ht_nb bench/cbbst2ht_nb_test.cpp ${CBBSTHT_SOURCE_FILES})
target_compile_definitions(cbbst2ht_nb PUBLIC "-DMINI_BLOCKS")
add_executable(bbst-fn_nb bench/bbst-fn_nb_test.cpp ${BBSTHT_SOURCE_FILES} ${HFERRADA_RMQ_SOURCE_FILES})
add_executable(cbbst-fn_nb bench/bbst-fn_nb_test.cpp ${CBBSTHT_SOURCE_FILES} ${HFERRADA_RMQ_SOURCE_FILES})
target_compile_definitions(cbbst-fn_nb PUBLIC "-DQUANTIZED")
add_executable(bbst2-fn_nb bench/bbst2-fn_nb_test.cpp ${BBSTHT_SOURCE_FILES} ${HFERRADA_RMQ_SOURCE_FILES})
target_compile_definitions(bbst2-fn_nb PUBLIC "-DMINI_BLOCKS")
add_executable(cbbst2-fn_nb bench/bbst2-fn_nb_test.cpp ${CBBSTHT_SOURCE_FILES} ${HFERRADA_RMQ_SOURCE_FILES})
target_compile_definitions(cbbst2-fn_nb PUBLIC "-DMINI_BLOCKS -DQUANTIZED")
