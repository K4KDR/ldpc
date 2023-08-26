cmake_minimum_required(VERSION 3.0)

add_executable(test_encoder test_encoder.cpp)
target_link_libraries(test_encoder libldpc)

add_executable(test_decoder test_decoder.cpp)
target_link_libraries(test_decoder libldpc)

add_executable(test_chain test_chain.cpp)
target_link_libraries(test_chain libldpc)

add_executable(test_benchmark test_benchmark.cpp)
target_link_libraries(test_benchmark libldpc)

add_executable(ber_simulation ber_simulation.cpp)
target_link_libraries(ber_simulation libldpc)


add_test(TestEncoder test_encoder)
add_test(TestDecoder test_decoder)
add_test(TestChain test_chain)
add_test(TestBenchmark test_benchmark)
