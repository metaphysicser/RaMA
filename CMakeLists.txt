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

include_directories(include)


# 将源代码添加到此项目的可执行文件。
add_executable (RaMA "RaMA.cpp"  "include/anchor.h" "include/common.h" "include/gsacak.h" "include/kseq.h" "include/logging.h" "include/pairwise_alignment.h" "include/rare_match.h" "include/threadpool.h" "include/utils.h" "src/anchor.cpp" "src/gsacak.c" "src/logging.cpp" "src/pairwise_alignment.cpp" "src/rare_match.cpp" "src/utils.cpp")

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET RaMA PROPERTY CXX_STANDARD 17)
endif()

# 添加这两行来找到并链接pthread库
find_package(Threads REQUIRED)
target_link_libraries(RaMA PRIVATE Threads::Threads)


# TODO: 如有需要，请添加测试并安装目标。
