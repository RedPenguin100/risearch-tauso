cmake_minimum_required(VERSION 3.20)
project(risearch1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

option(RISEARCH_BUILD_TESTS "Build the gtest suite" ON)

# Path RIsearch searches for on-disk energy tables. RIsearch1 has its matrices
# compiled in and never reads this, but the source requires the macro to exist.
set(RISEARCH_MATPATH "${CMAKE_CURRENT_SOURCE_DIR}/../tables" CACHE STRING
    "Directory holding energy matrix files")

set(RISEARCH1_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src)
# Every source under src/, at any depth, so a new subdirectory needs no CMake change.
file(GLOB_RECURSE RISEARCH1_SOURCES CONFIGURE_DEPENDS ${RISEARCH1_SRC}/*.cpp)

# Same for the include path: every directory that holds a header goes on it, so
# "Matrix.h" and "math/Matrix.h" both resolve.
file(GLOB_RECURSE RISEARCH1_HEADERS CONFIGURE_DEPENDS
     ${RISEARCH1_SRC}/*.h ${RISEARCH1_SRC}/*.hpp)
set(RISEARCH1_INCLUDES ${RISEARCH1_SRC})
foreach(header ${RISEARCH1_HEADERS})
    get_filename_component(dir ${header} DIRECTORY)
    list(APPEND RISEARCH1_INCLUDES ${dir})
endforeach()
list(REMOVE_DUPLICATES RISEARCH1_INCLUDES)

# The usage banner reports this fork's version, not upstream's. Read it from
# pyproject.toml so there is one place to bump; fall back when RIsearch1/ is
# configured outside the repository.
set(RISEARCH_TAUSO_VERSION "dev")
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../pyproject.toml")
    file(READ "${CMAKE_CURRENT_SOURCE_DIR}/../pyproject.toml" _pyproject)
    # The trailing group takes a PEP 440 pre-release suffix -- 1.6.0b1, 1.6.0rc2 --
    # so a beta reports the version it was published under rather than falling
    # through to "dev".
    string(REGEX MATCH "\nversion = \"([0-9]+\.[0-9]+\.[0-9]+[0-9a-z]*)\"" _unused "${_pyproject}")
    if(CMAKE_MATCH_1)
        set(RISEARCH_TAUSO_VERSION "${CMAKE_MATCH_1}")
    endif()
endif()
message(STATUS "RIsearch1: risearch-tauso version ${RISEARCH_TAUSO_VERSION}")

set(RISEARCH1_DEFINES MATPATH="${RISEARCH_MATPATH}" RISVERSION=1
    RISEARCH_TAUSO_VERSION="${RISEARCH_TAUSO_VERSION}")

set(RISEARCH1_WARNINGS -Wall -Wextra -pedantic
    -Wdouble-promotion -Wcast-qual -Wextra-semi -Wcast-align -Wzero-as-null-pointer-constant
    -Wpointer-arith -Wnon-virtual-dtor -Wduplicated-branches -Wduplicated-cond
    -Woverloaded-virtual -Wlogical-op -Wnull-dereference -Wformat=2
    -fvisibility-inlines-hidden -fvisibility=hidden)

# Split by where the flag belongs: -ffunction-sections/-fdata-sections put every
# function and object in its own section at compile time, and --gc-sections drops
# the unreferenced ones at link time. One without the other does nothing.
set(RISEARCH1_COMPILE_OPTS -O3 -ffunction-sections -fdata-sections)

# Unroll where gcc judges it worthwhile. The kernels are written as intrinsics, so
# nothing here relies on the vectoriser, but the sweep and the traceback both run
# short bodies whose loop overhead is a real share of their time: unrolling is
# worth about 5% of a run. Choosing the factor by hand is worse than leaving it
# to gcc -- forcing two on the sweep's inner loop alone measures slower.
list(APPEND RISEARCH1_COMPILE_OPTS -funroll-loops)

# Never split a 256-bit unaligned access into two 128-bit ones. GCC 10 and older
# do that by default: it suited Sandy Bridge, whose load path was 128 bits wide
# internally, and AMD Bulldozer, which cracked every 256-bit op in two. Haswell
# widened the path and the split became two loads and a vinserti128 where one
# vmovdqu would do.
#
# The vector kernels here need AVX2, which no chip that wanted the split has --
# AVX2 arrived with Haswell -- so the split is a loss for this program on every
# machine that runs those kernels at all. It costs about 45% of a sweep. GCC 11
# changed the default; the flag is what carries that back to older toolchains,
# which is what the manylinux wheels are built with.
include(CheckCXXCompilerFlag)
foreach(flag -mno-avx256-split-unaligned-load -mno-avx256-split-unaligned-store)
    string(MAKE_C_IDENTIFIER "HAVE${flag}" flag_var)
    check_cxx_compiler_flag(${flag} ${flag_var})
    if(${flag_var})
        list(APPEND RISEARCH1_COMPILE_OPTS ${flag})
    endif()
endforeach()
set(RISEARCH1_LINK_OPTS -Wl,--gc-sections)

# LTO needs the same flag on both the compile and the link line. CMake's
# INTERPROCEDURAL_OPTIMIZATION sets both, and check_ipo_supported keeps the build
# working on toolchains that lack it.
include(CheckIPOSupported)
check_ipo_supported(RESULT RISEARCH1_IPO_SUPPORTED OUTPUT RISEARCH1_IPO_MESSAGE)
if(NOT RISEARCH1_IPO_SUPPORTED)
    message(STATUS "RIsearch1: building without LTO -- ${RISEARCH1_IPO_MESSAGE}")
endif()

# Everything the optimised targets need. Kept in one place so the binary the
# wheel ships and the binary the benchmarks time are built the same way.
function(risearch1_optimise target)
    target_compile_options(${target} PRIVATE ${RISEARCH1_COMPILE_OPTS})
    target_link_options(${target} PRIVATE ${RISEARCH1_LINK_OPTS})
    if(RISEARCH1_IPO_SUPPORTED)
        set_property(TARGET ${target} PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    endif()
endfunction()

# --- fmt, which formats the hit lines ---------------------------------------
# Fetched at a pinned tag rather than taken from wherever the machine happens to
# have one. A build environment on the PATH ahead of the system -- conda is the
# common one -- can offer an fmt built against a different runtime, and mixing it
# with this build's LTO fails at link time. Packagers who want to link the fmt
# they ship can turn the option on.
option(RISEARCH_USE_SYSTEM_FMT "Link the installed fmt instead of fetching one" OFF)
if(RISEARCH_USE_SYSTEM_FMT)
    find_package(fmt 9 REQUIRED)
else()
    include(FetchContent)
    FetchContent_Declare(fmt
        GIT_REPOSITORY https://github.com/fmtlib/fmt.git
        GIT_TAG 11.0.2
        GIT_SHALLOW TRUE)
    FetchContent_MakeAvailable(fmt)
endif()

# setup.py copies bin/RIsearch into the wheel, so both binaries land where the
# Makefile used to put them rather than in the build tree.
set(RISEARCH1_BIN ${CMAKE_CURRENT_SOURCE_DIR}/bin)

# --- the CLI, optimised: what the wheel ships -------------------------------
add_executable(RIsearch ${RISEARCH1_SOURCES})
target_include_directories(RIsearch PRIVATE ${RISEARCH1_INCLUDES})
target_compile_definitions(RIsearch PRIVATE ${RISEARCH1_DEFINES})
target_compile_options(RIsearch PRIVATE ${RISEARCH1_WARNINGS})
target_link_libraries(RIsearch PRIVATE fmt::fmt-header-only)
risearch1_optimise(RIsearch)
set_target_properties(RIsearch PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${RISEARCH1_BIN})

# --- the same program with the debug printing compiled in -------------------
add_executable(RIsearch_dbg ${RISEARCH1_SOURCES})
target_include_directories(RIsearch_dbg PRIVATE ${RISEARCH1_INCLUDES})
target_compile_definitions(RIsearch_dbg PRIVATE ${RISEARCH1_DEFINES} DEBUG VERBOSE=2)
target_compile_options(RIsearch_dbg PRIVATE ${RISEARCH1_WARNINGS} -g -O0)
target_link_libraries(RIsearch_dbg PRIVATE fmt::fmt-header-only)
set_target_properties(RIsearch_dbg PROPERTIES
    OUTPUT_NAME RIsearch.dbg
    RUNTIME_OUTPUT_DIRECTORY ${RISEARCH1_BIN})

# --- the same translation units as a library, for the tests -----------------
# main() is renamed so the test binary can link these under gtest's own main().
# The rename is PRIVATE so it applies only to these sources.
add_library(risearch1_core STATIC ${RISEARCH1_SOURCES})
target_include_directories(risearch1_core PUBLIC ${RISEARCH1_INCLUDES})
target_compile_definitions(risearch1_core
    PUBLIC ${RISEARCH1_DEFINES}
    PRIVATE main=risearch1_main)
target_compile_options(risearch1_core PRIVATE ${RISEARCH1_WARNINGS})
# The benchmarks time this library, so it is built like the shipped binary.
risearch1_optimise(risearch1_core)
target_link_libraries(risearch1_core PUBLIC m fmt::fmt-header-only)

if(RISEARCH_BUILD_TESTS)
    add_subdirectory(tests)
endif()
