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

add_executable(bbstcon ${BBSTCON_SOURCE_FILES})
add_executable(bbst bench/bbst_test.cpp ${BBST_SOURCE_FILES})
add_executable(bbst_nb bench/bbst_nb_test.cpp ${BBST_SOURCE_FILES})
add_executable(bbst2_nb bench/bbst2_nb_test.cpp ${BBST_SOURCE_FILES})
target_compile_definitions(bbst2_nb PUBLIC "-DMINI_BLOCKS")
add_executable(bbst_nb_wc bench/bbst_nb_test.cpp ${BBST_SOURCE_FILES})
target_compile_definitions(bbst_nb_wc PUBLIC "-DWORST_CASE")
add_executable(bbst2_nb_wc bench/bbst2_nb_test.cpp ${BBST_SOURCE_FILES})
target_compile_definitions(bbst2_nb_wc PUBLIC "-DMINI_BLOCKS -DWORST_CASE")