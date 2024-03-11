# CMakeList.txt: RaMA 的 CMake 项目，在此处包括源代码并定义
# 项目特定的逻辑。
#
cmake_minimum_required (VERSION 3.8)

# 如果支持，请为 MSVC 编译器启用热重载。
if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

project ("RaMA")

# list (APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/Alignment/seqan3/build_system")

# find_package (seqan3 3.0 REQUIRED)

# Search for zlib as a dependency for SeqAn.
# find_package (ZLIB)
# find_package(OpenMP REQUIRED)
# if(OpenMP_CXX_FOUND)
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
# endif()

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2")

# 启用OpenMP
set(OPENMP ON CACHE BOOL "Enable OpenMP" FORCE)
# 添加编译器优化标志
set(EXTRA_FLAGS "-mavx2" CACHE STRING "Add optimization flags for C/C++ compiler" FORCE)

# find_package(Threads REQUIRED)
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

# list(APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/Alignment/seqan/util/cmake")
# Load the SeqAn module and fail if not found.
# find_package (SeqAn REQUIRED)

add_subdirectory(Alignment/WFA2-lib)
# list(APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/Alignment/WFA2-lib")
# find_package(WFA2-lib)

include_directories(${PROJECT_SOURCE_DIR} Anchor Alignment Logging Utils Alignment/WFA2-lib/)

# add_definitions (${SEQAN_DEFINITIONS})

# 将源代码添加到此项目的可执行文件。
add_executable (RaMA "RaMA.cpp"  "Anchor/anchor.h" "config.h" "Anchor/gsacak.h" "Utils/kseq.h" "Logging/logging.h" "Alignment/pairwise_alignment.h" "Anchor/rare_match.h" "Utils/threadpool.h" "Utils/utils.h" "Anchor/anchor.cpp" "Anchor/gsacak.c" "Logging/logging.cpp" "Alignment/pairwise_alignment.cpp" "Anchor/rare_match.cpp" "Utils/utils.cpp")

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET RaMA PROPERTY CXX_STANDARD 17)
endif()

# 添加这两行来找到并链接pthread库
find_package(Threads REQUIRED)
target_link_libraries(RaMA PRIVATE Threads::Threads ${SEQAN_LIBRARIES} wfa2cpp wfa2_static)


# TODO: 如有需要，请添加测试并安装目标。
