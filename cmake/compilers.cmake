set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

if(CMAKE_BUILD_TYPE STREQUAL Debug)
  add_compile_options(-g -O0)
else()
  add_compile_options(-O3)
endif()

if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  add_compile_options(-march=native -Werror=array-bounds -Wno-unused-variable)
  list(APPEND CFLAGS -Wno-unused-result -Wno-implicit-int)
endif()
