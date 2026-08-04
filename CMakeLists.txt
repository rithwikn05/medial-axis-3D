cmake_minimum_required(VERSION 3.16)
project(medial_axis_3d LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(MEDIAL_AXIS_ENABLE_CUDA "Enable CUDA surface-query acceleration" OFF)

include(FetchContent)

find_package(polyscope CONFIG QUIET)
if(NOT TARGET polyscope)
  set(POLYSCOPE_BUILD_EXAMPLES OFF CACHE BOOL "Build Polyscope examples" FORCE)
  set(POLYSCOPE_BUILD_TESTS OFF CACHE BOOL "Build Polyscope tests" FORCE)
  set(POLYSCOPE_BUILD_SHARED OFF CACHE BOOL "Build Polyscope as shared" FORCE)
  FetchContent_Declare(
    polyscope
    GIT_REPOSITORY https://github.com/nmwsharp/polyscope.git
    GIT_TAG master
  )
  FetchContent_MakeAvailable(polyscope)
endif()

add_library(medial_axis_3d_core STATIC
  src/core/predicates.c
  src/core/surface_query_cpu.cpp
)
target_include_directories(medial_axis_3d_core PUBLIC src)

if(MEDIAL_AXIS_ENABLE_CUDA)
  include(CheckLanguage)
  check_language(CUDA)
  if(NOT CMAKE_CUDA_COMPILER)
    message(FATAL_ERROR
      "MEDIAL_AXIS_ENABLE_CUDA is ON, but a CUDA compiler was not found. "
      "Install a compatible CUDA Toolkit or add nvcc to PATH."
    )
  endif()

  enable_language(CUDA)
  add_library(medial_axis_3d_cuda STATIC
    src/cuda/surface_query_cuda.cu
  )
  target_include_directories(medial_axis_3d_cuda PUBLIC src)
  target_link_libraries(medial_axis_3d_cuda PUBLIC medial_axis_3d_core)
  target_compile_definitions(medial_axis_3d_cuda
    PUBLIC MEDIAL_AXIS_HAS_CUDA=1
  )
  set_target_properties(medial_axis_3d_cuda PROPERTIES
    CUDA_STANDARD 17
    CUDA_STANDARD_REQUIRED ON
  )
endif()

add_executable(test_runner tests/test_runner.cpp)
target_include_directories(test_runner PRIVATE src)
target_link_libraries(test_runner PRIVATE medial_axis_3d_core)

add_executable(polyscope_viewer src/viewer/polyscope_viewer.cpp)
target_include_directories(polyscope_viewer PRIVATE src)
target_link_libraries(polyscope_viewer PRIVATE medial_axis_3d_core polyscope)
if(TARGET medial_axis_3d_cuda)
  target_link_libraries(polyscope_viewer PRIVATE medial_axis_3d_cuda)
endif()

enable_testing()
add_test(NAME medial_axis_3d_tests COMMAND test_runner)
