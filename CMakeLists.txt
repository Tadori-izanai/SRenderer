cmake_minimum_required(VERSION 3.10)
project(TinyRenderer)

set(CMAKE_CXX_STANDARD 17)

add_executable(TinyRenderer main.cpp geometry.cpp geometry.h gl.cpp gl.h model.cpp model.h tgaimage.h tgaimage.cpp)
