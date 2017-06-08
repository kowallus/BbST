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
        bench/bbst_test.cpp
        bbst.cpp
        bbst.h)

add_executable(bbstcon ${BBSTCON_SOURCE_FILES})
add_executable(bbst ${BBST_SOURCE_FILES})