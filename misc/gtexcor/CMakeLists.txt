cmake_minimum_required(VERSION 3.10)
project(gtexcor)

set (CMAKE_STATIC_LINKER_FLAGS "-Wl,--as-needed.-lcurl")

include(ExternalProject)

# thanks for the macro: David Grayson
macro(use_cxx11)
    if (CMAKE_VERSION VERSION_LESS "3.1")
        if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
            set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
        endif ()
    else ()
        set (CMAKE_CXX_STANDARD 11)
    endif ()
endmacro(use_cxx11)

use_cxx11()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")

include("cmake/htslib.cmake")

add_executable(gtexcor
               src/gtexcor.cpp
               src/arg_parse.h
               src/arg_parse.cpp
               src/Multimap.hh
               src/gff.h
               src/gff.cpp
               src/GBase.h
               src/GBase.cpp
               src/gdna.h
               src/gdna.cpp
               src/codons.h
               src/codons.cpp
               src/GFaSeqGet.h
               src/GFaSeqGet.cpp
               src/GList.hh
               src/GVec.hh
               src/GHash.hh)
add_dependencies(gtexcor htslib)
target_link_libraries(gtexcor ${htslib_LIB} -lbz2 -llzma -pthread -lcrypto -lcurl -lz)
