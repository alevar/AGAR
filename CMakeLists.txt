cmake_minimum_required(VERSION 3.10)
project(agar)

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

#use_cxx11()

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")

include("cmake/htslib.cmake")
include("cmake/hisat2.cmake")

include_directories(src)

add_executable(agar_indexer
        src/arg_parse.cpp
        src/arg_parse.h
        src/codons.cpp
        src/codons.h
        src/FastaTools.cpp
        src/FastaTools.h
        src/GBase.cpp
        src/GBase.h
        src/gdna.cpp
        src/gdna.h
        src/GFaSeqGet.cpp
        src/GFaSeqGet.h
        src/gff.cpp
        src/gff.h
        src/GHash.hh
        src/GList.hh
        src/Indexer.cpp
        src/Indexer.h
        src/GVec.hh
        src/tokenize.cpp
        src/Multimap.hh)

add_dependencies(agar_indexer htslib)
target_link_libraries(agar_indexer ${htslib_LIB} -lbz2 -llzma -pthread -lcrypto -lcurl -lz)

add_executable(agar
        src/arg_parse.cpp
        src/arg_parse.h
        src/codons.cpp
        src/codons.h
        src/GBase.cpp
        src/GBase.h
        src/GHash.hh
        src/GList.hh
        src/GVec.hh
        src/Converter.cpp
        src/Converter.h
        src/tokenize.cpp
        src/tokenize.h
        src/agar.cpp
        src/Multimap.hh)

add_dependencies(agar htslib)
target_link_libraries(agar ${htslib_LIB} -lbz2 -llzma -pthread -lcrypto -lcurl -lz)

add_executable(agar_diff
        src/arg_parse.cpp
        src/arg_parse.h
        src/codons.cpp
        src/codons.h
        src/GBase.cpp
        src/GBase.h
        src/GHash.hh
        src/GList.hh
        src/GVec.hh
        src/Comparator.cpp
        src/Comparator.h
        src/tokenize.cpp
        src/tokenize.h
        src/agar_diff.cpp
        src/Multimap.hh)

add_dependencies(agar_diff htslib)
target_link_libraries(agar_diff ${htslib_LIB} -lbz2 -llzma -pthread -lcrypto -lcurl -lz)