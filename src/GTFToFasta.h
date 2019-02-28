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

// #include "common.h"
#include "gff.h"
#include "GFaSeqGet.h"
#include "FastaTools.h"
#include "arg_parse.h"

// make -C /home/sparrow/genomicTools/tophat/ && ~/genomicTools/tophat/src/gtf_to_fasta 76 ~/genomicTools/tophat/ann.gff ~/genomicData/hg38/hg38_p8.fa ~/genomicTools/tophat/ann.fa
class GTFToFasta {
public:
    GTFToFasta(std::string gtf_fname, std::string genome_fname);
    ~GTFToFasta();
    void make_transcriptome(std::string out_fname, int kmer_length);

    // for debugging
    void print_mapping();
private:
    GffReader gtfReader_;
    // The genome_fhandle_ isn't used anywhere after being
    // initialized, and double opening the fasta file means that pipes
    // cannot work.
    // GFaSeqGet genome_fhandle_;

    std::string gtf_fname_;
    std::string genome_fname_;

    FILE* gtf_fhandle_;

    //convert set to unordered set
    std::set<std::vector<std::pair<int,int> > > kmer_coords; // genomic positions encountered
    std::unordered_map<std::string,std::vector<std::vector<std::pair<int,int>>>> kmers; // kmers and their positions for efficient lookup

    std::unordered_map<int,const char* const> int_to_chr; // conversion back from intigers to chromosomes
    std::unordered_map<const char*,int> chr_to_int; // conversion table from chromosome id to int

    // "contig" => vector(index_of_gff_obj)
//    typedef std::map<const char* , std::vector< size_t >* > ContigTransMap;
    typedef std::map<std::string, std::vector< int >* > ContigTransMap;

    ContigTransMap contigTransMap_;

    // std::transform(seq.cbegin(), seq.cend(), seq.begin(), lambda);

    void transcript_map();
    std::string get_exonic_sequence(GffObj& p_trans, FastaRecord& rec, std::string& coords, int kmer_length, std::ofstream& multimap);

    GTFToFasta(); // Don't want anyone calling the constructor w/o options
};

#endif
