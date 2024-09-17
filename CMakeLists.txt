cmake_minimum_required(VERSION 3.22)
project(capClone)

set(CMAKE_CXX_STANDARD 14)

add_executable(capClone main.cpp Instance.cpp Instance.h Utils.cpp Utils.h Antibody.cpp Antibody.h Search.cpp Search.h Config.h Stats.cpp Stats.h)
