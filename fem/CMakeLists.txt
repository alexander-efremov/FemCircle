cmake_minimum_required(VERSION 2.8)

project(fem)

add_definitions(-std=c++0x)

set(HEADER_FILES
    tecplot.h
    utils.h
    timer.h
    )

set(SOURCE_FILES
    tecplot.cpp
        consts.h solver.h solver.cpp)

add_library(fem STATIC ${SOURCE_FILES} ${HEADER_FILES})
