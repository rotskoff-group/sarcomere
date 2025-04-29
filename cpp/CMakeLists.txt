cmake_minimum_required(VERSION 3.0.0)

# Force CMake to use the compiler from the Conda environment
if(DEFINED ENV{CONDA_PREFIX})
    set(CMAKE_PREFIX_PATH $ENV{CONDA_PREFIX} CACHE PATH "Path to Conda environment")
    set(CMAKE_C_COMPILER $ENV{CONDA_PREFIX}/bin/gcc CACHE FILEPATH "C compiler")
    set(CMAKE_CXX_COMPILER $ENV{CONDA_PREFIX}/bin/g++ CACHE FILEPATH "C++ compiler")
    set(CMAKE_LIBRARY_PATH $ENV{CONDA_PREFIX}/lib CACHE PATH "Library path")
    set(CMAKE_INCLUDE_PATH $ENV{CONDA_PREFIX}/include CACHE PATH "Include path")
endif()

project(sarcomere)

# Set source and include directories
set(SRC_DIR "${PROJECT_SOURCE_DIR}/src")
set(INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")

# Include the header files
include_directories(${INCLUDE_DIR})

# Find required libraries
find_package(GSL REQUIRED)
#find_package(GTest REQUIRED)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

find_package(HDF5 REQUIRED COMPONENTS C CXX)
find_package(OpenMP REQUIRED)
include_directories("${PROJECT_SOURCE_DIR}/autodiff")

# Set compiler flags
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS_DEBUG "-g")  # Add debug symbols

# Source files
set(SOURCES
    ${SRC_DIR}/sarcomere.cpp
    ${SRC_DIR}/geometry.cpp
    ${SRC_DIR}/interaction.cpp
    ${SRC_DIR}/utils.cpp
    ${SRC_DIR}/neighborlist.cpp
    ${SRC_DIR}/h5_utils.cpp
    ${SRC_DIR}/components.cpp
    ${SRC_DIR}/langevin.cpp


)



# FetchContent is built into CMake from version 3.11 onward.
include(FetchContent)

# Optionally disable Google Benchmark's internal tests.
set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "Disable benchmark tests" FORCE)

# Declare and download the Google Benchmark repository.
FetchContent_Declare(
  benchmark
  GIT_REPOSITORY https://github.com/google/benchmark.git
  GIT_TAG v1.7.0  # Specify the desired tag/commit.
)

# Make the benchmark content available (it will be built as part of your project).
FetchContent_MakeAvailable(benchmark)

#check if the benchmark.cpp file exists
if(EXISTS ${SRC_DIR}/benchmark.cpp)
    # Now create your target and link it against Google Benchmark.
    add_executable(my_benchmark ${SRC_DIR}/benchmark.cpp)
    target_link_libraries(my_benchmark PRIVATE my_library benchmark::benchmark)
endif()


# Create the library (STATIC or SHARED)
add_library(my_library STATIC ${SOURCES})  # Use SHARED for a dynamic library
target_include_directories(my_library PUBLIC ${INCLUDE_DIR})

# Link required libraries
target_link_libraries(my_library PUBLIC ${GSL_LIBRARIES} ${HDF5_LIBRARIES} Eigen3::Eigen OpenMP::OpenMP_CXX)

# Executables
add_executable(${PROJECT_NAME} ${SRC_DIR}/main.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE my_library)

# Include directories for required libraries
target_include_directories(my_library PUBLIC ${GSL_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS})
target_include_directories(${PROJECT_NAME} PRIVATE ${GSL_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS})

# Link GSL, HDF5, and OpenMP
target_link_libraries(${PROJECT_NAME} PRIVATE ${GSL_LIBRARIES} ${HDF5_LIBRARIES} OpenMP::OpenMP_CXX)

# Google Test (Uncomment if needed)
# include_directories(${GTEST_INCLUDE_DIRS})
# add_executable(tests ${SRC_DIR}/tests.cpp)
# target_link_libraries(tests PRIVATE my_library gtest_main)
# include(GoogleTest)
# gtest_discover_tests(tests)

# Define custom target to run Valgrind
add_custom_target(valgrind
    COMMAND valgrind --leak-check=full ./${PROJECT_NAME}
    DEPENDS ${PROJECT_NAME}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Running Valgrind on ${PROJECT_NAME}..."
)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
