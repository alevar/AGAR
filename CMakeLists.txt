cmake_minimum_required(VERSION 3.10)
project(trans2genome)

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

add_executable(trans2genome_indexer
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

add_dependencies(trans2genome_indexer htslib)
target_link_libraries(trans2genome_indexer ${htslib_LIB} -lbz2 -llzma -pthread -lcrypto -lcurl -lz)

add_executable(trans2genome
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
        src/trans2genome.cpp
        src/Multimap.hh)

add_dependencies(trans2genome htslib)
target_link_libraries(trans2genome ${htslib_LIB} -lbz2 -llzma -pthread -lcrypto -lcurl -lz)
