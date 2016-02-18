cmake_minimum_required(VERSION 2.8)

project(calendars)

add_definitions(-std=c++0x)

set(HEADER_FILES
    tecplot.h
    utils.h
    timer.h
    )

set(SOURCE_FILES
    tecplot.cpp
    )

add_library(calendars STATIC ${SOURCE_FILES} ${HEADER_FILES})
