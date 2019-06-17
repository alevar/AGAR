//
//  gtfToFasta.h
//  TopHat
//
//  Created by Harold Pimentel on 10/26/11.
//

#ifndef GTFToFasta_H
#define GTFToFasta_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <sstream>
#include <utility>
#include <vector>
#include <unordered_map>
#include <set>
#include <algorithm>

#include "gff.h"
#include "GFaSeqGet.h"
#include "FastaTools.h"
#include "arg_parse.h"
#include "Multimap.hh"

class GTFToFasta {
public:
    GTFToFasta(std::string gtf_fname, std::string genome_fname,const std::string& out_fname, int kmer_length,bool multi);
    ~GTFToFasta();
    void make_transcriptome();

    // for debugging
    void print_mapping();
private:
    GffReader gtfReader_;

    std::string gtf_fname_;
    std::string genome_fname_;
    FILE* gtf_fhandle_;

    bool multi = false; // multimapper resolution
    int kmerlen; // kmer length to index
    std::string out_fname; // base name for all files
    int topTransID = 0; // the highest transcript ID assigned for any transcript in the current transcriptome. This information is written to the info file

    std::ofstream *infofp;
    std::ofstream *tlst;
    std::ofstream *multimap;
    std::ofstream *uniquefp;
    std::ofstream *genefp;
    std::ofstream *out_file;

    void transcript_map();
    std::string get_exonic_sequence(GffObj& p_trans, FastaRecord& rec, std::string& coords);

    GTFToFasta() = default; // Don't want anyone calling the constructor w/o options

    // MULTIMAPPERS

    //convert set to unordered set
    std::map<std::vector<std::pair<int,int> >, std::pair<std::string,int> > kmer_coords; //genomic positions encountered and counts of how many transcripts share a given position on the genome if the position is not a multimapper. The second feature in the value pointed to by the key shows the first transcript for which a position is defined. If only a single transcript contains the position, than it is its transcript id written
    std::unordered_map<std::string,std::vector<std::vector<std::pair<int,int>>>> kmers; // kmers and their positions for efficient lookup

    std::unordered_map<int,const char* const> int_to_chr; // conversion back from intigers to chromosomes
    std::unordered_map<const char*,int> chr_to_int; // conversion table from chromosome id to int

    typedef std::map<std::string, std::vector< int >* > ContigTransMap;
    ContigTransMap contigTransMap_;
};

#endif
