cmake_minimum_required(VERSION 3.10)

# 项目名称和版本
project(CFD VERSION 1.0)

# 设置 C 标准
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED True)

# 强制使用 GCC 编译器
set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)

# 指定源代码目录
set(SOURCE_DIR src)

# 添加可执行文件，包含 src 目录下的所有源文件
file(GLOB SOURCES "${SOURCE_DIR}/*.c")

add_executable(FVM2D ${SOURCES})

# 编译选项
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(FVM2D PRIVATE -O0 -g)
else()
    target_compile_options(FVM2D PRIVATE -O3)
endif()

# 链接数学库和 OpenMP
target_link_libraries(FVM2D m)

# 设置输出目录为当前目录下的 bin 文件夹
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

find_package(OpenMP)
if(OpenMP_C_FOUND)
    target_link_libraries(FVM2D OpenMP::OpenMP_C)
endif()

# 确保可执行文件输出目录设置
set_target_properties(FVM2D PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}
)
