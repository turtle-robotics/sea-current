cmake_minimum_required(VERSION 3.22)
project(sea_current LANGUAGES CXX)

include(GNUInstallDirs)
set(CMAKE_CXX_STANDARD 20)

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif ()
string(TOUPPER "${CMAKE_BUILD_TYPE}" BUILD_TYPE_UPPER)

if ("${BUILD_TYPE_UPPER}" STREQUAL "RELEASE")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build Type" FORCE)

  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
elseif ("${BUILD_TYPE_UPPER}" STREQUAL "DEBUG")
  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Build Type" FORCE)

  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")

  add_definitions(-DDEBUG)
endif ()

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

# examples
if (BUILD_EXAMPLES)
  add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/third-party/matplotlib-cpp")

  add_executable(test examples/test.cpp)
  target_link_libraries(test PRIVATE sea_current matplotlib_cpp)
  set_target_properties(test PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
endif ()

install(FILES
  "${PROJECT_SOURCE_DIR}/sea_current.hpp"
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
