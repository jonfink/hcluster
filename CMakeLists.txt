cmake_minimum_required (VERSION 2.6)
project (hcluster)

add_library(hcluster SHARED src/hcluster.cc)

add_executable(test_hcluster src/test_hcluster.cc)
target_link_libraries(test_hcluster hcluster)
