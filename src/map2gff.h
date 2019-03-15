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
#include <map>
#include <set>
// #include <boost/config.hpp>
// #include <boost/container_hash/hash.hpp>

#include <htslib/sam.h>

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
    std::vector<std::pair<int,int>> coords;
    int pos;
    uint32_t cigar_full[MAX_CIGARS];
    int n_cigar;
    std::string cigar;
    MateRead(bam1_t* curAl,int pos, std::string cigar,std::vector<std::pair<int,int>> al_coords,uint32_t *cf,int nc){
        this->al=curAl;
        this->pos=pos;
        this->cigar=cigar;
        this->coords=al_coords;
        memcpy(this->cigar_full,cf,MAX_CIGARS);
        this->n_cigar=nc;
        // this->cigar_full=cf;
    }
    ~MateRead(){
        // std::cout<<"DELETING 2"<<std::endl;
        bam_destroy1(al);
        // delete[] cigar_full;
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

// typedef std::vector<std::pair<int,int> > Coords;

struct coord_hash {
    uint64_t operator()(const std::vector<std::pair<int,int> > &coords ) const
    {
        uint64_t resHash=1;
        for (auto const& coord: coords) {
            resHash=resHash^std::hash<uint64_t>()(coord.first)^std::hash<uint64_t>()(coord.second);
        }
        return resHash;
    }
};

typedef std::tuple<int,int,int> coord_range;
struct coord_cmp {
    bool operator()(coord_range const & a, coord_range const & b) {
        if(std::get<2>(a) < std::get<1>(b) || std::get<0>(a)!=std::get<0>(b)){
            return true;
        }
        return false;
    }
};

class Map2GFF{
    public:
        Map2GFF(const std::string& tlstFP, const std::string& alFP, const std::string& multiFP, const std::string& glstFP, const bool& multi_flag);
        ~Map2GFF();

        void convert_coords(const std::string& outFP, const std::string& genome_header);
    
    private:
        bool get_read_start(GVec<GSeg>& exon_list,size_t gff_start, size_t& genome_start,int& exon_idx);
        void add_to_group();
        int convert_cigar(int i,int cur_intron_len,int miss_length,GSeg *next_exon,int match_length,GVec<GSeg>& exon_list,
                          int &num_cigars,int read_start,bam1_t* curAl,std::string &cigar_str,int cigars[MAX_CIGARS],std::vector<std::pair<int,int>> &coords);
        void add_multimapper_pair(const std::vector<std::pair<int,int>> *cor1, const std::vector<std::pair<int,int>> *cor2, bam1_t *al1, bam1_t *al2,
                          uint32_t *cur_cigar_full1,uint32_t *cur_cigar_full2, int n_cigar1, int n_cigar2);
        int merge_cigar(const std::vector<std::pair<int,int>> *cor,bam1_t *al, uint32_t *cur_cigar_full, int n_cigar);
        
        // coord_range make_coord_range(int strand, int lower, int upper);

        bool multi_flag=false;

        GPVec<GffTranscript> transcripts;
        std::unordered_map<std::string,GffTranscript*> tidx_to_t;

        std::ifstream tlststream;
        std::ifstream multistream;
        std::ifstream glststream;

        samFile *al; 
        bam_hdr_t *al_hdr;
        bam1_t *curAl;
        samFile *genome_al;
        bam_hdr_t *genome_al_hdr;

        std::unordered_map<std::string,MateRead*> curReadGroup_tmp;
        std::unordered_map<std::string,MateRead*> curReadGroup_tmp_multi;
        // std::unordered_map<std::string,MatePair*> curReadGroup_paired;

        std::map<
            std::pair<
                std::vector<
                    std::pair<int,int>
                >, std::vector<
                    std::pair<int,int>
                >
            >,MatePair*> curReadGroup_paired; // key is two vectors of coordinate pairs - one vector for each mate
        std::pair<
            std::map<
                std::pair<
                    std::vector<
                        std::pair<int,int>
                    >, std::vector<
                        std::pair<int,int>
                    >
                >,MatePair*>::iterator,
            bool> exists_curReadGroup_paired; // key is two vectors of coordinate pairs - one vector for each mate

        std::map<std::vector<std::pair<int,int>>,bam1_t*> curReadGroup;
        std::pair<std::map<std::vector<std::pair<int,int>>,bam1_t*>::iterator,bool> exists_curReadGroup;
        std::string curReadName="";

        std::unordered_map<std::string,int> ref_to_id; // built from the genome header file

        std::map<std::vector<std::pair<int,int> >,std::vector<const std::vector<std::pair<int,int> >*> > multimappers;
        std::unordered_map< std::string,int > ref_to_id_mult; // contains a map from chromosome name to the unique ID
        std::unordered_map< int,std::string > id_to_ref_mult; // contains a map from unique chromosome ID to the chromosome name

        std::map<
                std::vector< std::pair< int,int> >,
                std::vector<const std::vector<std::pair<int,int> >*>
            >::iterator exists_al_coord; // iterator to find the key in the multimappers map based on the coordinates of the current alignment
        std::map<
                std::vector< std::pair< int,int> >,
                std::vector<const std::vector<std::pair<int,int> >*>
            >::iterator exists_al_coord_mate1; // iterator to find the key in the multimappers map based on the coordinates of the current alignment
        std::map<
                std::vector< std::pair< int,int> >,
                std::vector<const std::vector<std::pair<int,int> >*>
            >::iterator exists_al_coord_mate2; // iterator to find the key in the multimappers map based on the coordinates of the current alignment

        std::map<coord_range,int,coord_cmp> gene_coords;
        std::map<coord_range,int>::iterator exists_geneID_mate1;
        std::map<coord_range,int>::iterator exists_geneID_mate2;
};

#endif
