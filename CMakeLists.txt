# CMake project configuration for RaMA.
cmake_minimum_required (VERSION 3.1)

project ("RaMA")
set(CMAKE_CXX_STANDARD 17)

# Check for GCC version
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9)
    message(FATAL_ERROR "GCC version must be at least 9 for full C++17 support.")
  endif()
endif()

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2")

# Automatically detect and enable OpenMP support.
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  message(STATUS "OpenMP found")
  set(OPENMP ON CACHE BOOL "Enable OpenMP" FORCE)
else()
  message(STATUS "OpenMP not found")
endif()

# -DUSE_M64=ON
option(USE_M64 "Build with M64 flag" OFF)

if(USE_M64)
  add_definitions(-DM64)
endif()

# -DEXTRA_FLAGS="-mavx2"
if(EXTRA_FLAGS)
  message(STATUS "Add optimization flags ${EXTRA_FLAGS} for C/C++ compiler")
  set(EXTRA_FLAGS ${EXTRA_FLAGS} CACHE STRING "Add optimization flags for C/C++ compiler" FORCE)
endif(EXTRA_FLAGS)

# Find and link Threads library.
find_package(Threads REQUIRED)
# Add subdirectory containing the WFA2-lib project.
add_subdirectory(Alignment/WFA2-lib)

# Include directories for the project.
include_directories(${PROJECT_SOURCE_DIR} Anchor Alignment ArgParser ThreadPool Logging Utils Alignment/WFA2-lib/)


# Define the executable for the project.
add_executable(RaMA "RaMA.cpp" "Anchor/anchor.h" "Anchor/gsacak.h" "Utils/kseq.h" "Logging/logging.h" 
"Alignment/pairwise_alignment.h" "Anchor/rare_match.h" "ThreadPool/threadpool.h" "Utils/utils.h" 
"Anchor/anchor.cpp" "Anchor/gsacak.c" "Logging/logging.cpp" "Alignment/pairwise_alignment.cpp" 
"Anchor/rare_match.cpp" "Utils/utils.cpp" "Anchor/RMQ.h" "Anchor/RMQ.cpp" "ArgParser/argparser.h")

# Set C++17 as the standard.
set_property(TARGET RaMA PROPERTY CXX_STANDARD 17)

target_link_libraries(RaMA PRIVATE Threads::Threads wfa2cpp wfa2_static)

