cmake_minimum_required (VERSION 3.1)
project(optflow C Fortran)
enable_testing()

include(cmake/compilers.cmake)

find_package(BLAS)

if(BLAS_FOUND)
  add_library(linpack SOFTWARE/dsvdc.f)
  target_link_libraries(linpack PUBLIC ${BLAS_LIBRARIES})
  target_compile_options(linpack PRIVATE -Wall)
  
  add_executable(fleet FLEET/fleet.c)
  target_link_libraries(fleet linpack m)
  target_compile_options(fleet PRIVATE ${CFLAGS})
  add_test(NAME Fleet
           COMMAND fleet new2rastertreet. 2. 9 0.05 
           ${CMAKE_SOURCE_DIR}/TESTDATA/TREE_DATA/TRANS 
           ${CMAKE_SOURCE_DIR}/out ${CMAKE_SOURCE_DIR}/out)
  
  add_executable(singh SINGH/singh.c)
  target_link_libraries(singh linpack)
  target_compile_options(singh PRIVATE ${CFLAGS})
  target_include_directories(singh PRIVATE SOFTWARE)
endif()

add_executable(flow2ps SOFTWARE/psflow.c)
target_link_libraries(flow2ps m)

add_executable(anandan ANANDAN/anandan.c)
target_link_libraries(anandan m)
target_compile_options(anandan PRIVATE ${CFLAGS})


add_executable(horn HORN/horn.c)
target_link_libraries(horn m)
target_include_directories(horn PRIVATE SOFTWARE)
add_test(NAME HornSchunck 
         COMMAND horn new2binarytreet. 0.5 1.5 20 100 
                 ${CMAKE_SOURCE_DIR}/TESTDATA/TREE_DATA/TRANS 
                 ${CMAKE_SOURCE_DIR}/out/ 
                 -B 150 150 -MH -T 5.0)

add_executable(lucas LUCAS/lucas.c)
target_link_libraries(lucas m)
target_include_directories(lucas PRIVATE SOFTWARE)
add_test(NAME LucasKanade
         COMMAND lucas new2binarytreet. 1.5 20 1.0 
                 ${CMAKE_SOURCE_DIR}/TESTDATA/TREE_DATA/TRANS 
                 ${CMAKE_SOURCE_DIR}/out/ 
                 -B 150 150)

add_executable(nagel NAGEL/nagel.c)
target_link_libraries(nagel m)
target_compile_options(nagel PRIVATE ${CFLAGS})

add_executable(uras URAS/uras.c)
target_link_libraries(uras m)
target_compile_options(uras PRIVATE ${CFLAGS})
