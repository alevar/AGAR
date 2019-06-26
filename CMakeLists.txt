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

#set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer -g")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=undefined -fno-omit-frame-pointer -g")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address -fno-omit-frame-pointer -g")

#include("cmake/zlib.cmake")
include("cmake/htslib.cmake")

include_directories(src)

add_executable(gtf_to_fasta
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
        src/GTFToFasta.cpp
        src/GTFToFasta.h
        src/GVec.hh
        src/tokenize.cpp)

add_executable(gtf_to_fasta_salmon
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
        src/GTFToFasta_salmon.cpp
        src/GTFToFasta_salmon.h
        src/GVec.hh
        src/tokenize.cpp
        src/Multimap.hh)

add_dependencies(gtf_to_fasta_salmon htslib)
target_link_libraries(gtf_to_fasta_salmon ${htslib_LIB} -lbz2 -llzma -pthread -lcrypto -lcurl -lz)

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
        src/map2gff.cpp
        src/map2gff.h
        src/tokenize.cpp
        src/tokenize.h
        src/trans2genome.cpp
        src/Multimap.hh)

add_dependencies(trans2genome htslib)
target_link_libraries(trans2genome ${htslib_LIB} -lbz2 -llzma -pthread -lcrypto -lcurl -lz)

add_executable(salmon2genome
        src/arg_parse.cpp
        src/arg_parse.h
        src/codons.cpp
        src/codons.h
        src/GBase.cpp
        src/GBase.h
        src/GHash.hh
        src/GList.hh
        src/GVec.hh
        src/map2gff_salmon.cpp
        src/map2gff_salmon.h
        src/tokenize.cpp
        src/tokenize.h
        src/salmon2genome.cpp
        src/Multimap.hh)

add_dependencies(salmon2genome htslib)
target_link_libraries(salmon2genome ${htslib_LIB} -lbz2 -llzma -pthread -lcrypto -lcurl -lz)
