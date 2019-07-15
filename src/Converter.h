//
// Created by sparrow on 6/6/19.
//

#ifndef CONVERTER_H
#define CONVERTER_H

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <cstdlib>
#include <vector>
#include <cstring>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <math.h>

//#include <htslib/khash.h> // TODO: consider moving the implementation to a better hashmap (such as khash or something else)
#include <htslib/sam.h>

#include "GVec.hh"
#include "Multimap.hh"

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

#define MAX_EDITS 100

#include <cassert>

// TODO: for the final version it needs to scan through the first n reads to check for constraints
//    such as whether multimappers are reported, and whatever else
//    or the number of bases per read

// TODO: needs to be parallelized

// class for the unmapped reads
// evaluates reads and outputs them accordingly
// TODO: consider making this into a probabilistic lookup in order to reduce memory
//     alternatively, assuming a sorted output, can now simply reset when new readname appears
class UMAP{
public:
    UMAP() = default; // do not want anyone to call from outside
    UMAP(const std::string& outFP){ // the output file allows saving the output into a specific file as opposed to passing it to stdout
        this->outFP = outFP;
        std::string out_fname_r1(this->outFP);
        std::string out_fname_r2(this->outFP);
        std::string out_fname_s(this->outFP);
        out_fname_r1.append(".r1.fa");
        out_fname_r2.append(".r2.fa");
        out_fname_s.append(".s.fa");
        this->out_stream_r1 = new std::ofstream(out_fname_r1.c_str());
        this->out_stream_r2 = new std::ofstream(out_fname_r2.c_str());
        this->out_stream_s = new std::ofstream(out_fname_s.c_str());
    }
    ~UMAP(){
        if(!this->outFP.empty()){
            this->out_stream_r1->close();
            this->out_stream_r2->close();
            this->out_stream_s->close();
        }
        delete this->out_stream_r1;
        delete this->out_stream_r2;
        delete this->out_stream_s;
    }

    void set_outFP(const std::string& outFP){
        this->outFP = outFP;
        std::string out_fname_r1(this->outFP);
        std::string out_fname_r2(this->outFP);
        std::string out_fname_s(this->outFP);
        out_fname_r1.append(".r1.fa");
        out_fname_r2.append(".r2.fa");
        out_fname_s.append(".s.fa");
        this->out_stream_r1 = new std::ofstream(out_fname_r1.c_str());
        this->out_stream_r2 = new std::ofstream(out_fname_r2.c_str());
        this->out_stream_s = new std::ofstream(out_fname_s.c_str());
    }

    // this function verifies whether a read needs to be written to the output
    void insert(bam1_t* al){
//        std::cout<<this->unpaired1.size()<<"\t"<<this->unpaired2.size()<<"\t"<<this->ps1.size()<<"\t"<<this->ps2.size()<<std::endl;
        this->clean(al);
        // also need to check if the second pair does not
        if(this->pair_unmapped(al)){
            if(this->is1(al)){
                this->pse = this->ps2.find(bam_get_qname(al));
                if(this->pse != this->ps2.end()){ // already exists and is waiting to be returned
                    // in this case we can simply write out both records
                    this->outputFasta(al,*this->out_stream_r1);
                    this->outputFasta(this->pse->second,*this->out_stream_r2);
                    // remove from the stack
                    this->ps2.erase(this->pse);
                    // add the entry to unpaired
                    this->unpaired1.insert(bam_get_qname(al));
                }
                else{
                    if(this->unpaired1.insert(bam_get_qname(al)).second){ // successfully inserted a new key
                        // need to add to the paired_stack to be outputted when the second mate is found
                        this->ps1.insert(std::make_pair(bam_get_qname(al),al));
                    }
                }
            }
            else{
                this->pse = this->ps1.find(bam_get_qname(al));
                if(this->pse != this->ps1.end()){ // already exists and is waiting to be returned
                    // in this case we can simply write out both records
                    this->outputFasta(al, *this->out_stream_r2);
                    this->outputFasta(this->pse->second,*this->out_stream_r1);
                    // remove from the stack
                    this->ps1.erase(this->pse);
                    this->unpaired2.insert(bam_get_qname(al));
                }
                else{
                    if(this->unpaired2.insert(bam_get_qname(al)).second){ // successfully inserted a new key
                        // need to add to the paired_stack to be outputted when the second mate is found
                        this->ps2.insert(std::make_pair(bam_get_qname(al),al));
                    }
                }
            }
        }
        else if(this->is1(al)){
            if(this->unpaired1.insert(bam_get_qname(al)).second){ // did not previously exist
                this->outputFasta(al,*this->out_stream_s);
            }
        }
        else{
            if(this->unpaired2.insert(bam_get_qname(al)).second){ // did not previously exist
                this->outputFasta(al,*this->out_stream_s);
            }
        }
    }
private:
    std::string last_readname = ""; // used to erase entries for sorted-by-name input to optimize lookup and memory
    std::string outFP;
    std::ofstream *out_stream_r1,*out_stream_r2,*out_stream_s;
    // keep a map of the unmapped elements to directly output to hisat2
    std::unordered_set<std::string> unpaired1; // for first mates
    std::unordered_set<std::string> unpaired2; // for second mates
    std::pair<std::unordered_set<std::string>::iterator,bool> ue1,ue2; // an entry exists;

    std::unordered_map<std::string,bam1_t*> ps1,ps2; // stack for holding the paired end information
    std::unordered_map<std::string,bam1_t*>::iterator pse;

    // this method evaluates the state of the maps based on readnames and cleans it
    void clean(bam1_t *al){
        if(bam_get_qname(al) != this-> last_readname){ // if new read name
            // clean everything up
            this->unpaired1.clear();
            this->unpaired2.clear();
            this->ps1.clear();
            this->ps2.clear();
            this->last_readname = bam_get_qname(al);
        }
    }

    void outputFasta(bam1_t *al,std::ofstream& out_stream){
        // reports a fasta record from a given alignment record
        std::string read_name= bam_get_qname(al);
        out_stream<<">"<<read_name<<std::endl;
        // now deal with the sequence
        uint8_t * seq = bam_get_seq(al);
        for(int i=0;i<al->core.l_qseq;i++){
            out_stream<<seq_nt16_str[bam_seqi(seq,i)];
        }
        out_stream<<std::endl;
    }
    bool is1(bam1_t *al){return al->core.flag & 64;} // test if the current read is first in the pair
    bool is2(bam1_t *al){return al->core.flag & 128;} // test if the current read is second in the pair
    bool pair_unmapped(bam1_t *al){return (al->core.flag & 4) && (al->core.flag & 8);}; // test if bth reads in the pair are unmapped
};

// this structure defines all the things related to a genomic position of a transcriptomic alignment
struct ReadData{
    int read_start;
    int cigar[MAX_CIGARS];
    int num_cigars;
    int tid;
};

// this class defines the unique identification of a read mapping
class MapID{
public:
    explicit MapID(bam1_t *al){
        // populate object with all the info
        this->name = bam_get_qname(al);
        this->refid = al->core.tid; // assumes that tid == mtid
        if(al->core.flag & 0x40){
            this->pos1 = al->core.pos;
            this->pos2 = al->core.mpos;
        }
        else{
            this->pos1 = al->core.mpos;
            this->pos2 = al->core.pos;
        }
    }
    ~MapID()=default;

    bool operator==(const MapID& m) const{
        return this->name==m.get_name() &&
                this->refid==m.get_refid() &&
                this->pos1==m.get_pos1() &&
                this->pos2==m.get_pos2();
    }

    bool operator<(const MapID& m) const{
        return this->name<m.get_name() &&
               this->refid<m.get_refid() &&
               this->pos1<m.get_pos1() &&
               this->pos2<m.get_pos2();

    }

    std::string get_name() const{return this->name;}
    int get_refid() const{return this->refid;}
    int get_pos1() const{return this->pos1;}
    int get_pos2() const{return this->pos2;}
private:
    MapID()=default; // do not want to be called
    std::string name;
    int refid;
    int pos1,pos2; // positions of two mates
};

// hash function to be used for
namespace std {
    template<>
    struct hash<MapID> {
        size_t operator()(const MapID &k) const {
            return ((hash<string>()(k.get_name()) ^
                     (
                             (hash<int>()(k.get_refid()) ^ hash<int>()(k.get_pos1()) ^
                              hash<int>()(k.get_pos2())
                             ) << 1)) >> 1);
        }
    };
}

// this class facilitates efficient identification of mates from the same paired-end read
// when a read is added to the class, the pair is looked up and if found is reported back
// pairs are identified by the readname and respective positions (reference id and read start)
// once identified, the pair is immediately removed from the set
/*
* SAMPLE WORKFLOW
* 1. have a class to hold and release reads
*      - this class will do the mate resolution
*      - the class is needed since we can not guarantee without sorting by readname that reads will appear in the correct order
*          however, with high likelihood they will, and that can be exploited for efficiency
*
* 2. for a paired read we simply wait until the mate is found and then proceed to process them together
*
*/
class Pairs{
public:
    Pairs() = default;
    ~Pairs(){
        for(auto &v : this->mates){
            bam_destroy1(v.second);
        }
    };

    int add(bam1_t *al,bam1_t *mate){ // add read to the stack
        MapID m(al);
        bam1_t *al_dup = bam_init1();
        bam_copy1(al_dup,al);
        me = mates.insert(std::make_pair(m,al_dup));
        if(!me.second){ // entry previously existed - can report a pair and remove from the stack
            bam_copy1(mate,me.first->second);
            bam_destroy1(me.first->second);
            mates.erase(me.first);
            return 1;
        }
        return 0;
    }
private:
    std::unordered_map<MapID,bam1_t*> mates;
    std::pair<std::unordered_map<MapID,bam1_t*>::iterator,bool> me;
};

// this class describes the unique identifier of a genomic maping of a read
class ReadGenKey{
public:
    ReadGenKey() = default;
    ~ReadGenKey() = default;
    std::string name;
    int tid;
    int start;
    int mate_start;
    int cigar; // TODO: if this class is used, need to consider how to pass the cigar information
    int mate_cigar; // TODO: same here as above

    bool operator==(const ReadGenKey& m) const{
        return this->name==m.name &&
               this->tid==m.tid &&
               this->start==m.start &&
               this->cigar==m.cigar &&
               this->mate_start==m.mate_start &&
               this->mate_cigar==m.mate_cigar;
    }
};

// hash function to be used for
namespace std {
    template<>
    struct hash<ReadGenKey> {
        size_t operator()(const ReadGenKey &k) const {
            return ((hash<string>()(k.name) ^
                     (
                             (hash<int>()(k.tid) ^
                              hash<int>()(k.start) ^
                              hash<int>()(k.mate_start) ^
                              hash<int>()(k.cigar) ^
                              hash<int>()(k.mate_cigar)
                             ) << 1)) >> 1);
        }
    };
}

// the following class describes resolution of transcriptomic multimappers
// it notifies wheter a given mapping needs to be returned or not
class Collapser{
public:
    Collapser() = default;
    ~Collapser() = default;
    int add(bam1_t *curAl,size_t cigar_hash){ // add single read to the stack
        this->clean(curAl);
        ReadGenKey rgk;
        rgk.name = bam_get_qname(curAl);
        rgk.tid = curAl->core.tid;
        rgk.start = curAl->core.pos;
        rgk.cigar = cigar_hash;
        rgk.mate_start = curAl->core.mpos;
        rgk.mate_cigar = cigar_hash;
        gpe = this->genomic_positions.insert(rgk);
        if(!gpe.second){ // entry previously existed - can report a pair and remove from the stack
            return 0;
        }
        return 1;
    }
    int add(bam1_t *curAl,bam1_t *mateAl,size_t cigar_hash,size_t mate_cigar_hash){ // add pair to the stack
        this->clean(curAl);
        ReadGenKey rgk;
        rgk.name = bam_get_qname(curAl);
        rgk.tid = curAl->core.tid;
        if(curAl->core.flag & 0x40){ // first in a pair
            rgk.start = curAl->core.pos;
            rgk.cigar = cigar_hash;
            rgk.mate_start = mateAl->core.pos;
            rgk.mate_cigar = mate_cigar_hash;
        }
        else{
            rgk.start = mateAl->core.pos;
            rgk.cigar = mate_cigar_hash;
            rgk.mate_start = curAl->core.pos;
            rgk.mate_cigar = cigar_hash;
        }
        gpe = this->genomic_positions.insert(rgk);
        if(!gpe.second){ // entry previously existed - can report a pair and remove from the stack
            return 0;
        }
        return 1;
    }
private:
    std::string last_readname = "";
    std::unordered_set<ReadGenKey> genomic_positions;
    std::pair<std::unordered_set<ReadGenKey>::iterator,bool> gpe;

    // when a new read is detected then this method removes any old entries for memory and lookup efficiency
    void clean(bam1_t *al){
//        std::cout<<this->genomic_positions.size()<<std::endl;
        if(bam_get_qname(al) != this->last_readname){ // if new read is detected
//            std::cout<<"2"<<std::endl;
            // clear contents of the containers
            this->genomic_positions.clear();
            this->last_readname = bam_get_qname(al);
        }
    }
};

// this clss stores and evaluates information regading edit distances
class ErrorCheck{
public:
    ErrorCheck() = default;
    ~ErrorCheck() = default;

    void set_stdv(int stdv){
        this->stdv = stdv;
    }

    std::pair<double,double> get_thresh(){
        return std::make_pair(this->lower_bound,lower_bound_pair);
    }

    // TODO: values need to be recomputed not every turn
    bool add_read(bam1_t *al){ // returns true if read passes error-check; otherwise returns false; also appends reads to the distribution
        if(num_reads%500000 == 0){
            std::cout<<"1\t"<<lower_bound<<std::endl;
        }
        cur_nm = get_nm(al);
        num_reads++;
        sum_nm=sum_nm+cur_nm;
        mean_nm = sum_nm/num_reads;
        if(!ztest){ // compute based on standard deviation
            std2 = std::sqrt(mean_nm);
            lower_bound = this->mean_nm+(std2*this->stdv); // two standard deviations
            this->observed[cur_nm]++;
        }
        else{ // perform z test
            // TODO: implement z test here
        }

        return cur_nm<=lower_bound;
    }
    bool add_pair(bam1_t *al,bam1_t *mate){
        if(num_reads_pair%500000 == 0){
            std::cout<<"2\t"<<lower_bound_pair<<std::endl;
        }
        cur_nm_pair = get_nm(al)+get_nm(mate);
        num_reads_pair++;
        sum_nm_pair=sum_nm_pair+cur_nm_pair;
        mean_nm_pair = sum_nm_pair/num_reads_pair;
        if(!ztest){ // compute based on standard deviation
            std2_pair = std::sqrt(mean_nm_pair);
            lower_bound_pair = this->mean_nm_pair+(std2_pair*this->stdv_pair); // two standard deviations
            this->observed_pair[cur_nm_pair]++;
        }
        else{ // perform z test
            // TODO: implement z test here
        }

        return cur_nm_pair<=lower_bound_pair;
    }

private:
    int stdv = 3; // for standard deviation based outlier detection this parameter sets the number of standard deviations above the mean that the
    bool ztest = false;
    std::vector<uint32_t> observed = std::vector<uint32_t>(MAX_EDITS);
    int num_reads=0,sum_nm=0,cur_nm=0;
    float mean_nm,std2,lower_bound;

    int stdv_pair = 3; // for standard deviation based outlier detection this parameter sets the number of standard deviations above the mean that the
    std::vector<uint32_t> observed_pair = std::vector<uint32_t>(MAX_EDITS);
    int num_reads_pair=0,sum_nm_pair=0,cur_nm_pair=0;
    float mean_nm_pair,std2_pair,lower_bound_pair;

    uint32_t get_nm(bam1_t *al){
        uint8_t* ptr_nm_1=bam_aux_get(al,"NM");
        if(ptr_nm_1){
            return bam_aux2i(ptr_nm_1);
        }
        else{
            return md2nm();
        }
    }

    uint32_t md2nm() {
        std::cerr<<"Computing from md tag"<<std::endl;
        return 0; // TODO: add parsing of MD tag if NM is not unavailable
    }
};

class Converter{
public:
    Converter(const std::string& alFP,const std::string& outFP,const std::string& index_base,const int& threads,bool multi);
    ~Converter();

    void load_abundances(const std::string& abundFP); // add abundances to the transcripts

    void convert_coords();
    void precompute(int perc);
    void print_multimappers();
    void set_unaligned();
    void set_k1(){this->k1_mode = true;};
    void set_fraglen(int fraglen);
    void set_num_multi(int num_multi){this->mmap.set_num_multi(num_multi);};
    void set_all_multi(){this->mmap.set_all_multi();};
    void set_misalign();
    void set_stdv(int stdv){errorCheck.set_stdv(stdv);}

    void print_stats(){
        std::cerr<<"reads discarded as misalignments: "<<(this->num_err_discarded_pair*2)+this->num_err_discarded<<" at "<<errorCheck.get_thresh().first<<" singleton and "<<errorCheck.get_thresh().second<<" paired thresholds"<<std::endl;
        std::cerr<<"\tof which "<<this->num_err_discarded<<" were singles"<<std::endl;
        std::cerr<<"\tand "<<this->num_err_discarded_pair*2<<" were paired"<<std::endl;

        std::cerr<<"multimapped reads: "<<(this->num_multi_pair*2)+this->num_multi<<std::endl;
        std::cerr<<"\tof which "<<this->num_multi<<" were singles at multimapping rate of "<<double(num_multi_hits)/double(num_multi)<<std::endl;
        std::cerr<<"\tand "<<this->num_multi_pair*2<<" were paired at multimapping rate of "<<double(num_multi_hits_pair)/double(num_multi_pair)<<std::endl;
    }

private:
    // STATS
    int num_err_discarded = 0,num_err_discarded_pair = 0;

    // MISALIGNMENT METHODS AND DECLARATIONS
    ErrorCheck errorCheck;
    bool detect_misalign = false;
    bool evaluate_errors(bam1_t *curAl); // returns true if the read passes the error check
    bool evaluate_errors_pair(bam1_t *curAl,bam1_t *mate);
    void write_unaligned(bam1_t *curAl,std::ofstream &out_ss);
    void write_unaligned_pair(bam1_t *curAl,bam1_t *mate);
    std::string unal_r1_fname,unal_r2_fname,unal_s_fname;
    std::ofstream unal_r1,unal_r2,unal_s;
    int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 };
    int8_t *buf = NULL;
    size_t max_buf = 0;
    int perc_precomp = 10;
    int total_num_pair_al = 0,total_num_al = 0;
    int total_multi_detected;

    // INDEX METHODS
    void load_index(const std::string& index_base,bool multi);
    void load_info(const std::string& info_fname);
    void load_transcriptome(const std::string& tlst_fname); // add transcriptome positional information into the index
    void _load_transcriptome(std::ifstream& tlstfp,char* buffer);
    void load_genome_header(const std::string& genome_header_fname);
    void load_multi(const std::string& multiFP);

    // INDEX-SPECIFIC DECLARATIONS
    bool multi = false; // set to true if the multimapper index is loaded
    bool abund = false;
    bool k1_mode = false;
    bool unaligned_mode = false;
    int numThreads=1;
    int fraglen = false;
    int numTranscripts=0; // gtf_to_fasta returns a .info file with this information
    int maxLocID = 0; // what is the highest geneID assigned for the GFF by the index-builder
    int kmerlen; // kmer length used for index construction

    // INDEX DATA
    Multimap mmap;
    Loci loci;
    std::vector<GffTranscript> transcriptome;

    // ALIGNMENT SPECIFIC DECLARATIONS
    std::string alFP,outFP;

    samFile *al;
    bam_hdr_t *al_hdr;
    samFile *genome_al;
    bam_hdr_t *genome_al_hdr;
    samFile *outSAM;
    bam_hdr_t *outSAM_header;

    std::unordered_map<std::string,int> ref_to_id; // built from the genome header file

    UMAP umap;
    Pairs pairs;
    Collapser collapser;

    // ALIGNMENT METHODS
    int convert_cigar(int i,GSeg *next_exon,GVec<GSeg>& exon_list,int &num_cigars,int read_start,bam1_t* curAl,int cigars[MAX_CIGARS],Position& pos_obj);
    bool has_valid_mate(bam1_t *curAl);
    bool has_mate(bam1_t *curAl);
    bool get_read_start(GVec<GSeg>& exon_list,int32_t gff_start,int32_t& genome_start, int& exon_idx);
    void add_cigar(bam1_t *curAl,int num_cigars,int* cigars);
    void add_aux(bam1_t *curAl,char xs);
    void change_nh_flag(bam1_t *curAl,int nh);
    void fix_flag(bam1_t *curAl);
    int collapse_genomic(bam1_t *curAl,size_t cigar_hash);
    int collapse_genomic(bam1_t *curAl,bam1_t *mateAl,size_t cigar_hash,size_t mate_cigar_hash);
    void process_pair(bam1_t* curAl);
    void process_single(bam1_t* curAl);
    size_t process_read(bam1_t* curAl,Position& cur_pos,int cigars[MAX_CIGARS],int &num_cigars);
    void finish_read(bam1_t *curAl);
    void add_multi_tag(bam1_t* curAl);

    // Multimapper-related methods
    int evaluate_multimappers(bam1_t* curAl,Position& cur_pos,int cigars[MAX_CIGARS],int &num_cigars);
    int evaluate_multimappers_pair(bam1_t *curAl,bam1_t* curAl_mate,Position &cur_pos,Position &cur_pos_mate,
                                    int *cigars,int *cigars_mate,int &num_cigars,int &num_cigars_mate);
    int num_multi = 0,num_multi_pair = 0,num_multi_hits=0,num_multi_hits_pair=0;

    // various printers relevant only for debug
    void print_transcriptome();

};

#endif //CONVERTER_H
