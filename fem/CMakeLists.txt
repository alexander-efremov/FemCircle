cmake_minimum_required(VERSION 2.8)

project(fem)

add_definitions(-std=c++0x)

set(HEADER_FILES
    utils.h
    timer.h
    )

set(SOURCE_FILES
       consts.h solver1.h solver1.cpp solver2.h solver2.cpp consts.cpp common.h)

add_library(fem STATIC ${SOURCE_FILES} ${HEADER_FILES})
