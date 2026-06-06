cmake_minimum_required(VERSION 3.10)
project(face-off VERSION 1.0.0 LANGUAGES CXX)

OPTION(Debug "Set to ON to compile with debug symbols" OFF)
option(G "Set to ON to compile with optimisations and debug symbols" OFF)
option(NATIVE_COMPILATION "Set to OFF to compile without the -march=native flag. This may be required when compiling binaries to be used elsewhere" ON)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# The default behaviour is to compile with optimisations
set(CMAKE_BUILD_TYPE Release)

IF(Debug)
        SET(CMAKE_BUILD_TYPE "Debug")
        MESSAGE(STATUS "Compiling in Debug mode")
ENDIF(Debug)

if(G)
    set(CMAKE_BUILD_TYPE RelWithDebInfo)
    message(STATUS "Compiling in Release+Debug symbols mode")
endif(G)

if(NOT Debug)
    if(INTEL)
        ADD_DEFINITIONS(-fast)
    else()
        ADD_DEFINITIONS(-ffast-math)
        if(NATIVE_COMPILATION)
            include(CheckCXXCompilerFlag)
            unset(COMPILER_SUPPORTS_MARCH_NATIVE CACHE)
            CHECK_CXX_COMPILER_FLAG(-march=native COMPILER_SUPPORTS_MARCH_NATIVE)
            if(COMPILER_SUPPORTS_MARCH_NATIVE)
                ADD_DEFINITIONS(-march=native)
            else()
                message(STATUS "The '-march=native' flag is not supported, compiling without it")
            endif()
        else()
            message(STATUS "Compiling without the '-march=native' flag")
        endif()
    endif(INTEL)
endif()

add_definitions(-Wall -Wshadow -Wextra -Wno-unused-parameter)

include_directories(extern/cxxopts-3.3.1)

# Find fmt library
find_package(fmt REQUIRED)

# Add the executable
add_executable(face-off src/main.cpp)

# Link fmt library
target_link_libraries(face-off fmt::fmt)
