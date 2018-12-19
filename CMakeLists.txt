cmake_minimum_required(VERSION 2.8.8)
set(PROJECT_NAME_STR sinoduplex)
project(${PROJECT_NAME_STR} C CXX)
# external dependencies to be discovered:
set(SAMTOOLS_DIR ${PROJECT_SOURCE_DIR}/thirdparty/samtools/samtools-1.3.1)
set(HTSLIB_DIR ${SAMTOOLS_DIR}/htslib-1.3.1)

if (APPLE)
    # For Macs, override the paths here
else ()
    # to compile in aliyun VM using g++-6.3.0
    include_directories(/mnt/pipeline-programs/gcc/gcc_6.3.0/include)
    set(CMAKE_CXX_COMPILER "/mnt/pipeline-programs/gcc/gcc_6.3.0/bin/g++")
    set(CMAKE_C_COMPILER "/mnt/pipeline-programs/gcc/gcc_6.3.0/bin/gcc")

#    include_directories(/usr/include/c++/4.8.5)
#
#    set(CMAKE_CXX_COMPILER "/usr/bin/g++")
#    set(CMAKE_C_COMPILER "/usr/bin/gcc")

    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath -Wl,/usr/local/lib64")

endif ()

# no version requirement at the moment
find_package(Threads REQUIRED)
find_package(ZLIB REQUIRED)

# Discover samtools and htslib
find_library(SAMTOOLS_BAM_LIB NAMES bam HINTS ${SAMTOOLS_DIR} STATIC IMPORTED)
find_library(HTS_LIB NAMES hts HINTS ${HTSLIB_DIR} STATIC IMPORTED)

# Compiler flags
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    # using Clang
    add_definitions(-Wall -g -O2 -std=c++11)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # using GCC
    add_definitions(-Wall -g -O2  -g -std=c++11)
    #add_definitions(-Wall -g -std=c++11 -fprofile-arcs -ftest-coverage)
endif ()

# Create the installdir.h include file
set(REFERENCE_DATA_MAKER "./make_ref_data.sh")
add_subdirectory(${PROJECT_SOURCE_DIR}/src)




