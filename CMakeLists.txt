cmake_minimum_required(VERSION 3.26)
project(MES)

set(CMAKE_CXX_STANDARD 17)

add_executable(MES main.cpp
        wczytywanie.cpp
        ElementUniversal.cpp
        ElementUniversal.h
        Gauss.cpp
        Gauss.h
        struktury.cpp
        struktury.h
)
