cmake_minimum_required(VERSION 2.8.9)

project(test)
set(SOURCES test.cpp polIMParton.cpp)

add_library(polIMParton SHARED polIMParton.cpp)
add_executable(polPDFtest ${SOURCES})
