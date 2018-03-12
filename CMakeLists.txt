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
        bbstx.cpp
        bbstx.h)

set(CBBSTX_SOURCE_FILES
        ${COMMON_SOURCE_FILES}
        cbbstx.h
        cbbstx.hpp)

set(HFERRADA_RMQ_SOURCE_FILES
        RMQRMM64.cpp
        includes/RMQRMM64.h
        includes/Basic_rmq.cpp
        includes/Basic_rmq.h)

set(SDSL-LITE_FILES
        sdsl/memory_management.cpp
        sdsl/ram_fs.cpp
        sdsl/bits.cpp
        sdsl/util.cpp
        sdsl/bp_support_algorithm.cpp
        includes/sdsl/rmq_support_sparse_table.hpp
        includes/sdsl/rmq_succinct_sada.hpp
        includes/sdsl/rmq_succinct_sct.hpp
        includes/sdsl/rmq_support.hpp)

add_executable(bbstcon ${BBSTCON_SOURCE_FILES})
add_executable(bbst bench/bbst_test.cpp ${BBST_SOURCE_FILES})
add_executable(bbst2 bench/bbst2_test.cpp ${BBST_SOURCE_FILES})
target_compile_definitions(bbst2 PUBLIC "-DMINI_BLOCKS")
add_executable(bbst_nb bench/bbst_nb_test.cpp ${BBST_SOURCE_FILES})
add_executable(bbst_pm_nb bench/bbst_nb_test.cpp ${BBST_SOURCE_FILES})
target_compile_definitions(bbst_pm_nb PUBLIC "-DPSEUDO_MONO")
add_executable(bbst2_nb bench/bbst2_nb_test.cpp ${BBST_SOURCE_FILES})
target_compile_definitions(bbst2_nb PUBLIC "-DMINI_BLOCKS")
add_executable(bbst2_pm_nb bench/bbst2_nb_test.cpp ${BBST_SOURCE_FILES})
target_compile_definitions(bbst2_pm_nb PUBLIC "-DMINI_BLOCKS -DPSEUDO_MONO")
add_executable(bbst_nb_wc bench/bbst_nb_test.cpp ${BBST_SOURCE_FILES})
target_compile_definitions(bbst_nb_wc PUBLIC "-DWORST_CASE")
add_executable(bbst2_nb_wc bench/bbst2_nb_test.cpp ${BBST_SOURCE_FILES})
target_compile_definitions(bbst2_nb_wc PUBLIC "-DMINI_BLOCKS -DWORST_CASE")
add_executable(bbstx_nb bench/bbstx_nb_test.cpp ${BBSTHT_SOURCE_FILES})
add_executable(bbst2x_nb bench/bbst2x_nb_test.cpp ${BBSTHT_SOURCE_FILES})
target_compile_definitions(bbst2x_nb PUBLIC "-DMINI_BLOCKS")
add_executable(cbbstx_nb bench/cbbstx_nb_test.cpp ${CBBSTX_SOURCE_FILES})
add_executable(cbbst2x_nb bench/cbbst2x_nb_test.cpp ${CBBSTX_SOURCE_FILES})
target_compile_definitions(cbbst2x_nb PUBLIC "-DMINI_BLOCKS")
add_executable(bbst-bp_nb bench/bbst-bp_nb_test.cpp ${BBSTHT_SOURCE_FILES} ${HFERRADA_RMQ_SOURCE_FILES})
add_executable(cbbst-bp_nb bench/bbst-bp_nb_test.cpp ${CBBSTX_SOURCE_FILES} ${HFERRADA_RMQ_SOURCE_FILES})
target_compile_definitions(cbbst-bp_nb PUBLIC "-DQUANTIZED")
add_executable(bbst2-bp_nb bench/bbst2-bp_nb_test.cpp ${BBSTHT_SOURCE_FILES} ${HFERRADA_RMQ_SOURCE_FILES})
target_compile_definitions(bbst2-bp_nb PUBLIC "-DMINI_BLOCKS")
add_executable(cbbst2-bp_nb bench/bbst2-bp_nb_test.cpp ${CBBSTX_SOURCE_FILES} ${HFERRADA_RMQ_SOURCE_FILES})
target_compile_definitions(cbbst2-bp_nb PUBLIC "-DMINI_BLOCKS -DQUANTIZED")

add_subdirectory(includes/sdsl/mman EXCLUDE_FROM_ALL)
include_directories(AFTER includes)

add_executable(bbst-sdsl-bp_nb bench/bbst-sdsl-bp_nb_test.cpp ${BBSTHT_SOURCE_FILES} ${SDSL-LITE_FILES})
add_executable(cbbst-sdsl-bp_nb bench/bbst-sdsl-bp_nb_test.cpp ${CBBSTX_SOURCE_FILES} ${SDSL-LITE_FILES})
target_link_libraries(bbst-sdsl-bp_nb PUBLIC mman)
target_link_libraries(cbbst-sdsl-bp_nb PUBLIC mman)
target_compile_definitions(cbbst-sdsl-bp_nb PUBLIC "-DQUANTIZED")
add_executable(bbst2-sdsl-bp_nb bench/bbst2-sdsl-bp_nb_test.cpp ${BBSTHT_SOURCE_FILES} ${SDSL-LITE_FILES})
add_executable(cbbst2-sdsl-bp_nb bench/bbst2-sdsl-bp_nb_test.cpp ${CBBSTX_SOURCE_FILES} ${SDSL-LITE_FILES})
target_link_libraries(bbst2-sdsl-bp_nb PUBLIC mman)
target_compile_definitions(bbst2-sdsl-bp_nb PUBLIC "-DMINI_BLOCKS")
target_link_libraries(cbbst2-sdsl-bp_nb PUBLIC mman)
target_compile_definitions(cbbst2-sdsl-bp_nb PUBLIC "-DMINI_BLOCKS -DQUANTIZED")

add_executable(bbst-sdsl-rec_nb bench/bbst-sdsl-rec_nb_test.cpp ${BBSTHT_SOURCE_FILES} ${SDSL-LITE_FILES})
add_executable(cbbst-sdsl-rec_nb bench/bbst-sdsl-rec_nb_test.cpp ${CBBSTX_SOURCE_FILES} ${SDSL-LITE_FILES})
target_link_libraries(bbst-sdsl-rec_nb PUBLIC mman)
target_link_libraries(cbbst-sdsl-rec_nb PUBLIC mman)
target_compile_definitions(cbbst-sdsl-rec_nb PUBLIC "-DQUANTIZED")
add_executable(bbst2-sdsl-rec_nb bench/bbst2-sdsl-rec_nb_test.cpp ${BBSTHT_SOURCE_FILES} ${SDSL-LITE_FILES})
add_executable(cbbst2-sdsl-rec_nb bench/bbst2-sdsl-rec_nb_test.cpp ${CBBSTX_SOURCE_FILES} ${SDSL-LITE_FILES})
target_link_libraries(bbst2-sdsl-rec_nb PUBLIC mman)
target_compile_definitions(bbst2-sdsl-rec_nb PUBLIC "-DMINI_BLOCKS")
target_link_libraries(cbbst2-sdsl-rec_nb PUBLIC mman)
target_compile_definitions(cbbst2-sdsl-rec_nb PUBLIC "-DMINI_BLOCKS -DQUANTIZED")