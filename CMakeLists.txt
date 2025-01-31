# CMake project configuration for RaMA.
cmake_minimum_required (VERSION 3.1)

project(RaMA VERSION 1.2.0)
# Set C++17 as the project standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF) # Prefer -std=c++xx over -std=gnu++xx

# Ensure GCC version is at least 9 for full C++17 support
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9)
  message(FATAL_ERROR "GCC version must be at least 9 for full C++17 support.")
endif()

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

# Automatically detect and enable OpenMP support.
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  message(STATUS "OpenMP found")
  set(OPENMP ON CACHE BOOL "Enable OpenMP" FORCE)
else()
  message(STATUS "OpenMP not found")
endif()

# Option to build with M64 flag
option(USE_M64 "Build with M64 flag" OFF)
if(USE_M64)
  add_compile_definitions(M64)
endif()

# -DEXTRA_FLAGS="-mavx2"
# Add optimization flags if provided
set(EXTRA_FLAGS "" CACHE STRING "Additional optimization flags for C/C++ compiler")
if(EXTRA_FLAGS)
  message(STATUS "Adding optimization flags: ${EXTRA_FLAGS} for C/C++ compiler")
  add_compile_options(${EXTRA_FLAGS})
endif()

# Find and link Threads library.
find_package(Threads REQUIRED)
# Add subdirectory containing the WFA2-lib project.
add_subdirectory(Alignment/WFA2-lib)

# Include directories
include_directories(
  ${PROJECT_SOURCE_DIR}
  Anchor
  Alignment
  ArgParser
  ThreadPool
  Logging
  Utils
  Alignment/WFA2-lib/
)

set(SOURCE_FILES
      Anchor/anchor.h Anchor/gsacak.h Utils/kseq.h Logging/logging.h 
  Alignment/pairwise_alignment.h Anchor/rare_match.h ThreadPool/threadpool.h 
  Utils/utils.h Anchor/anchor.cpp Anchor/gsacak.c Logging/logging.cpp 
  Alignment/pairwise_alignment.cpp Anchor/rare_match.cpp Utils/utils.cpp 
  Anchor/RMQ.h Anchor/RMQ.cpp ArgParser/argparser.h
)

# Specify the target executable and its sources
add_executable(RaMA 
  RaMA.cpp 
  ${SOURCE_FILES}
)
target_compile_options(RaMA PRIVATE -O3)
target_link_libraries(RaMA PRIVATE Threads::Threads wfa2cpp wfa2_static)


# Create a new executable target for testing RMQ with test_RMQ.cpp
# add_executable(test_RMQ 
#  test_RMQ.cpp 
#  ${SOURCE_FILES}
#)

# target_compile_options(test_RMQ PRIVATE -O3)
# Add OpenMP flags if found for test_RMQ
#if(OpenMP_CXX_FOUND)
#  target_link_libraries(test_RMQ PRIVATE OpenMP::OpenMP_CXX)
#endif()

# Link necessary libraries for test_RMQ
#target_link_libraries(test_RMQ PRIVATE Threads::Threads wfa2cpp wfa2_static)

