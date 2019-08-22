cmake_minimum_required (VERSION 3.1)
project(optflow C)
enable_testing()

include(cmake/compilers.cmake)

if(MSVC OR (WIN32 AND (CMAKE_C_COMPILER_ID STREQUAL Intel OR CMAKE_C_COMPILER_ID STREQUAL PGI)))
  set(math_lib)
else()
  set(math_lib m)
endif()

# if working Fortran compiler present, enable Fortran before searching for libraries
include(CheckLanguage)
check_language(Fortran)
if(CMAKE_Fortran_COMPILER)
  enable_language(Fortran)

  find_package(BLAS)

  if(BLAS_FOUND)
    add_library(linpack SOFTWARE/dsvdc.f)
    target_link_libraries(linpack PUBLIC ${BLAS_LIBRARIES})

    add_executable(fleet FLEET/fleet.c)
    target_link_libraries(fleet linpack ${math_lib})
    target_compile_options(fleet PRIVATE ${CFLAGS})
    add_test(NAME Fleet
      COMMAND fleet new2rastertreet. 2. 9 0.05
      ${CMAKE_SOURCE_DIR}/TESTDATA/TREE_DATA/TRANS
      ${CMAKE_SOURCE_DIR}/out ${CMAKE_SOURCE_DIR}/out)

    add_executable(singh SINGH/singh.c)
    target_link_libraries(singh linpack)
    target_compile_options(singh PRIVATE ${CFLAGS})
    target_include_directories(singh PRIVATE SOFTWARE)
    add_test(NAME Singh
      COMMAND singh new2binarytreet. 20
      ${CMAKE_SOURCE_DIR}/TESTDATA/TREE_DATA/TRANS
      ${CMAKE_SOURCE_DIR}/out -B 150 150
      -i 25 -n 2 -w 2 -N 4)
  endif()
endif()

add_executable(flow2ps SOFTWARE/psflow.c)
target_link_libraries(flow2ps ${math_lib})
target_compile_options(flow2ps PRIVATE ${CFLAGS})

add_executable(anandan ANANDAN/anandan.c)
target_link_libraries(anandan ${math_lib})
target_compile_options(anandan PRIVATE ${CFLAGS})
add_test(NAME Anandan
 COMMAND anandan
 ${CMAKE_SOURCE_DIR}/TESTDATA/TREE_DATA/TRANS
 ${CMAKE_SOURCE_DIR}/out/
 new2binarytreet. -N 9 10 -W 5 -L 3 -I 15)


add_executable(horn HORN/horn.c)
target_link_libraries(horn ${math_lib})
target_include_directories(horn PRIVATE SOFTWARE)
target_compile_options(horn PRIVATE ${CFLAGS})
add_test(NAME HornSchunck
  COMMAND horn new2binarytreet. 0.5 1.5 20 100
  ${CMAKE_SOURCE_DIR}/TESTDATA/TREE_DATA/TRANS
  ${CMAKE_SOURCE_DIR}/out/
  -B 150 150 -MH -T 5.0)


add_executable(lucas LUCAS/lucas.c)
target_link_libraries(lucas ${math_lib})
target_include_directories(lucas PRIVATE SOFTWARE)
add_test(NAME LucasKanade
  COMMAND lucas new2binarytreet. 1.5 20 1.0
  ${CMAKE_SOURCE_DIR}/TESTDATA/TREE_DATA/TRANS
  ${CMAKE_SOURCE_DIR}/out/
  -B 150 150)

add_executable(nagel NAGEL/nagel.c)
target_link_libraries(nagel ${math_lib})
target_compile_options(nagel PRIVATE ${CFLAGS})
#add_test(NAME Nagel
#  COMMAND nagel
#  ${CMAKE_SOURCE_DIR}/TESTDATA/TREE_DATA/TRANS
#  ${CMAKE_SOURCE_DIR}/out/
#  newbinarytreet. -SG 1.5 -TG 1.5 -M 20 -A 0.5 -I 100 -B 150 150 -F 5.0)

add_executable(uras URAS/uras.c)
target_link_libraries(uras ${math_lib})
target_compile_options(uras PRIVATE ${CFLAGS})
#add_test(NAME Uras
#  COMMAND uras
#  ${CMAKE_SOURCE_DIR}/TESTDATA/TREE_DATA/TRANS
#  ${CMAKE_SOURCE_DIR}/out/
#  new2binarytreet. -SG 3.0 -TG 1.5 -M 9 -R 1 -F 0.5)
