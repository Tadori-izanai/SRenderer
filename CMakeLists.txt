cmake_minimum_required(VERSION 3.10)
project(TinyRenderer)

set(CMAKE_CXX_STANDARD 17)

add_executable(TinyRenderer main.cpp geometry.cpp geometry.h gl.cpp gl.h model.cpp model.h tgaimage.h tgaimage.cpp)

target_compile_options(TinyRenderer PRIVATE -O3)

# OpenMP
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message("OPENMP FOUND")
    target_link_libraries(TinyRenderer PUBLIC OpenMP::OpenMP_CXX)
endif()
