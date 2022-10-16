cmake_minimum_required(VERSION 3.22)
project(sea_current LANGUAGES CXX)

include(GNUInstallDirs)
set(CMAKE_CXX_STANDARD 20)

add_library(sea_current INTERFACE)

# install(
#   TARGETS sea_current
#   # EXPORT install_targets
# )


# examples
add_executable(test examples/test.cpp)
target_link_libraries(test PRIVATE sea_current)
set_target_properties(test PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

install(FILES
  "${PROJECT_SOURCE_DIR}/sea_current.hpp"
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
