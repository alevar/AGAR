//
//  gtfToFasta.h
//  TopHat
//
//  Created by Ales Varabyou on 07/26/19 based on gtf_to_fasta from Tophat written by Harold Pimentel
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
    void print_mmap();

    void save(bool multi, bool uniq);
private:
    GffReader gtfReader_;

    std::string gtf_fname_;
    std::string genome_fname_;
    FILE* gtf_fhandle_;

    bool multi = false; // multimapper resolution
    int kmerlen; // kmer length to index
    std::string out_fname; // base name for all files
    int topTransID = 0; // the highest transcript ID assigned for any transcript in the current transcriptome. This information is written to the info file

    std::string trans_fastaname;
    std::ofstream *infofp;
    std::ofstream *tlst;
    std::ofstream *multimap;
    std::ofstream *uniquefp;
    std::ofstream *genefp;

    void transcript_map();
    std::string get_exonic_sequence(GffObj& p_trans, FastaRecord& rec, std::string& coords);
    void add_to_geneMap(GffObj &p_trans);

    GTFToFasta() = default; // Don't want anyone calling the constructor w/o options

    // MULTIMAPPERS
    Multimap mmap;

    std::unordered_map<int,const char* const> int_to_chr; // conversion back from intigers to chromosomes
    std::unordered_map<const char*,int> chr_to_int; // conversion table from chromosome id to int

    typedef std::map<std::string, std::vector< int >* > ContigTransMap;
    ContigTransMap contigTransMap_;

    int curGeneID = 0;
    std::map<std::string,std::tuple<int,int,int,int>> geneMap; // stores minimum and maximum gene coordinates of transcripts in a given gene as well as strand and unique ID
    std::pair< std::map<
            std::string,
            std::tuple<int,int,int,int>
    >::iterator,bool> exists_cur_gene;
    std::map<std::string,std::tuple<int,int,int,int>>::iterator found_gene;
};

#endif
