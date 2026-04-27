cmake_minimum_required(VERSION 3.16)

project(ldsc LANGUAGES CXX)

find_package(ZLIB REQUIRED)
find_package(BZip2 REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(ldsc
    src/ldsc.cpp
    src/common.cpp
    src/regression.cpp
    src/sumstats_commands.cpp
    src/munge.cpp
    src/annotation.cpp
    src/ld_score.cpp
    src/continuous.cpp
)

target_include_directories(ldsc PRIVATE include)
target_link_libraries(ldsc PRIVATE ZLIB::ZLIB BZip2::BZip2)

if(MSVC)
    target_compile_options(ldsc PRIVATE /W4)
else()
    target_compile_options(ldsc PRIVATE -Wall -Wextra -pedantic)
endif()
