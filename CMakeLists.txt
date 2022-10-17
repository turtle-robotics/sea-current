cmake_minimum_required(VERSION 3.22)
project(sea_current LANGUAGES CXX)

include(GNUInstallDirs)
set(CMAKE_CXX_STANDARD 20)

add_library(sea_current INTERFACE)
target_include_directories(sea_current
  INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/examples>
    $<INSTALL_INTERFACE:include>
)

find_package(Eigen3 3.4)
target_link_libraries(sea_current INTERFACE
    Eigen3::Eigen
)

# tests
add_executable(test examples/test.cpp)
target_link_libraries(test PRIVATE sea_current)
set_target_properties(test PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

install(FILES
  "${PROJECT_SOURCE_DIR}/sea_current.hpp"
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
