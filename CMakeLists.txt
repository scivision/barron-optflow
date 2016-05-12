cmake_minimum_required (VERSION 2.8.12)
project(optflow C Fortran)

add_compile_options(-march=native -Wall -Wpedantic -Werror=array-bounds)

include_directories(SOFTWARE)
link_libraries(m)

add_library(linpack SOFTWARE/svddouble.f)
target_compile_options(linpack PRIVATE -std=legacy)

add_executable(flow2ps SOFTWARE/psflow.c)

add_executable(horn HORN/horn.c)

add_executable(lucas LUCAS/lucas.c)

add_executable(anandan ANANDAN/anandan.c)

add_executable(fleet FLEET/fleet.c)
target_link_libraries(fleet linpack)

add_executable(singh SINGH/singh.c)
target_link_libraries(singh linpack)
