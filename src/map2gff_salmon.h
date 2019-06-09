//
// Created by sparrow on 6/6/19.
//

#ifndef TRANS2GENOME_MAP2GFF_SALMON_H
#define TRANS2GENOME_MAP2GFF_SALMON_H

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <cassert>
#include <cstdlib>
#include <vector>
#include <cstring>
#include <unordered_map>

#include <htslib/sam.h>

#include "GVec.hh"

#define BAM_CMATCH  0
#define BAM_CINS    1
#define BAM_CDEL    2
#define BAM_CREF_SKIP   3
#define BAM_CSOFT_CLIP  4
#define BAM_CHARD_CLIP  5
#define BAM_CPAD    6
#define BAM_CEQUAL  7
#define BAM_CDIFF   8
#define MAX_CIGARS  1024

struct GffTranscript: public GSeg {
    GVec<GSeg> exons;
    int numID; //numeric ID in tlst
    std::string gffID;
    std::string refID;
    char strand;
    GffTranscript():exons(1), numID(-1), gffID(),
                    refID(), strand(0) { }

    std::string& getRefName() {
        return refID;
    }
    explicit GffTranscript(const std::string& tline);
};

class Map2GFF_SALMON{
public:
    Map2GFF_SALMON(const std::string& tlstFP, const std::string& alFP, const int& threads);
    ~Map2GFF_SALMON() = default;

    void convert_coords(const std::string& outFP, const std::string& genome_header);

private:
    bool get_read_start(GVec<GSeg>& exon_list,size_t gff_start, size_t& genome_start,int& exon_idx);
    int convert_cigar(int i,int cur_intron_len,int miss_length,GSeg *next_exon,int match_length,GVec<GSeg>& exon_list,
                      int &num_cigars,int read_start,bam1_t* curAl,std::string &cigar_str,int cigars[MAX_CIGARS],std::vector<std::pair<int,int>> &coords);
    int merge_cigar(const std::vector<std::pair<int,int>> *cor,bam1_t *al, uint32_t *cur_cigar_full, int n_cigar);

    int numThreads=1;

    GPVec<GffTranscript> transcripts;
    std::unordered_map<std::string,GffTranscript*> tidx_to_t;

    std::ifstream tlststream;

    samFile *al;
    bam_hdr_t *al_hdr;
    bam1_t *curAl;
    samFile *genome_al;
    bam_hdr_t *genome_al_hdr;

    std::unordered_map<std::string,int> ref_to_id; // built from the genome header file

};


#endif //TRANS2GENOME_MAP2GFF_SALMON_H
