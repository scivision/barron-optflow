cmake_minimum_required (VERSION 2.8.12)
project(optflow C)

add_compile_options(-march=native -Wall -Wpedantic -Werror=array-bounds)

include_directories(SOFTWARE)
link_libraries(m)

add_executable(flow2ps SOFTWARE/psflow.c)

add_executable(horn HORN/horn.c)

add_executable(lucas LUCAS/lucas.c)

add_executable(anandan ANANDAN/anandan.c)
