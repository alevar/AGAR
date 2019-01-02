#ifndef _MAP2GFF_H_
#define _MAP2GFF_H_

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

#include "include/htslib/sam.h"

#include "GVec.hh"

// #define BAM_CIGAR_SHIFT     4
// #define BAM_CIGAR_MASK      ((1<<BAM_CIGAR_SHIFT)-1)

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

struct MatePair{
    bam1_t* firstMate;
    bam1_t* secondMate;
    MatePair(bam1_t* al1,bam1_t* al2){
        this->firstMate=al1;
        this->secondMate=al2;
    }
    ~MatePair(){
        // std::cout<<"DELETING 1"<<std::endl;
        bam_destroy1(firstMate);
        bam_destroy1(secondMate);
    }
};

struct MateRead{
    bam1_t *al;
    int pos;
    std::string cigar;
    MateRead(bam1_t* curAl,int pos, std::string cigar){
        this->al=curAl;
        this->pos=pos;
        this->cigar=cigar;
    }
    ~MateRead(){
        // std::cout<<"DELETING 2"<<std::endl;
        bam_destroy1(al);
    }
};

struct Map_Al{
    bam1_t* read;
    std::string sub_key_own;  // p_trans->refID+"$"+std::to_string(read_start)
    std::string sub_key_mate; // mate_p_trans->refID+"$"+std::to_string(mate_read_start)
    int origRefID;
    Map_Al(bam1_t* curAl,std::string sub_key_own,std::string sub_key_mate,int origRefID){
        this->read=curAl;
        this->sub_key_own=sub_key_own;
        this->sub_key_mate=sub_key_mate;
        this->origRefID=origRefID;
    }
    ~Map_Al(){
        // std::cout<<"DELETING 3"<<std::endl;
        bam_destroy1(read);
    }
};

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
	GffTranscript(const std::string& tline);
};

class Map2GFF{
    public:
        Map2GFF(const std::string& gffFP, const std::string& alFP);
        ~Map2GFF();

        void convert_coords(const std::string& outFP, const std::string& genome_header);
    
    private:
        bool get_read_start(GVec<GSeg>& exon_list,size_t gff_start, size_t& genome_start,int& exon_idx);
        void add_to_group();
        int convert_cigar(int i,int cur_intron_len,int miss_length,GSeg *next_exon,int match_length,GVec<GSeg>& exon_list,
                          int &num_cigars,int read_start,bam1_t* curAl,std::string &cigar_str,int cigars[MAX_CIGARS]);

        GPVec<GffTranscript> transcripts;
        std::unordered_map<std::string,GffTranscript*> tidx_to_t;

        std::ifstream tlststream;

        samFile *al; 
        bam_hdr_t *al_hdr;
        bam1_t *curAl;
        samFile *genome_al;
        bam_hdr_t *genome_al_hdr;

        std::unordered_map<std::string,int> ref_to_id;
};

#endif
