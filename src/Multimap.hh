//
// Created by Ales Varabyou on 6/16/19.
//

#ifndef TRANS2GENOME_MULTIMAP_H
#define TRANS2GENOME_MULTIMAP_H

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <cstdlib>
#include <vector>
#include <cstring>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <random>
#include <ctime>
#include <algorithm>

#include "gff.h"
#include <htslib/sam.h>

class PID{
public:
    PID() = default;
    PID(double kp, double ki){
        // set tunings here
        if (kp<0 || ki<0){
            std::cerr<<"@ERROR::wrong PID tunings specified"<<std::endl;
            exit(1);
        };

        this->kp = kp;
        this->ki = ki;
        prev_input = 0;

    }
    ~PID() = default;

    double compute(double &setpoint,double &pv){
        /*Compute all the working error variables*/
        double error = setpoint - pv;
        this->I_term += error;

        /*Compute PID Output*/
        double res = (this->kp * error) + (this->I_term * this->ki);

        /*Remember some variables for next time*/
        this->prev_input = pv;
        return res;
    }

private:
    double kp, ki;                  // Tuning Parameters
    double prev_input = 0;

    double I_term = 0;

};

// this class holds a description of a single Locus
// the description includes:
// 1. effective length of the locus (excludes the intronic sequences)
// 2. number of reads
// Allows the following methods:
// 1. get abundance from reads and effective length
// 2. count reads and etc
class Locus{
public:
    Locus()=default;
    Locus(uint32_t elen){
        this->elen = elen;
    }
    ~Locus()=default;

    void set_elen(uint32_t elen){this->elen = elen;this->empty = false;}
    void set_start(uint32_t start){this->start=start;this->empty=false;}
    void set_end(uint32_t end){this->end=end;this->empty=false;}
    void set_strand(uint32_t strand){this->strand=strand,this->empty=false;}
    void set_id(uint32_t locid){this->locid=locid;this->empty=false;}
    void set_abund(float abund){this->abundance = abund;}

    uint32_t get_start(){return this->start;}
    float get_abund(){return this->abundance;}
    uint32_t get_end(){return this->end;}
    uint32_t get_locid(){return this->locid;}
    uint8_t get_strand(){return this->strand;}
    uint32_t get_elen(){return this->elen;}
    uint32_t get_nreads(){return this->n_reads;}
    bool is_empty(){return this->empty;}

    // return rpk difference which is used to increment the normalizing constant
    double inc(){ // increments the read count and recomputes the TPM
        this->n_reads++;
        double old_rpk = this->rpk;
        this->rpk = (double)this->n_reads/(double)this->elen;
        this->rpk_total = (double)(this->n_reads_multi+this->n_reads)/(double)this->elen;
        return this->rpk - old_rpk;
    }

    double inc_multi(){ // this function is used to keep track of how many multimapping reads have been allocated to this locus
        this->n_reads_multi++;
        double old_rpk_total = this->rpk_total;
        this->rpk_total = (double)(this->n_reads_multi+this->n_reads)/(double)this->elen;
        return this->rpk_total - old_rpk_total;
    }

    double get_rpk(){ // return the RPK
        return this->rpk;
    }

    double get_rpk_total(){
        return this->rpk_total;
    }

    void clear(){
        this->empty=true;
        this->elen=1;
        this->rpk=0.0;
        this->n_reads=0;
    }

    void print(){
        std::cout<<this->locid<<":"<<this->elen<<this->strand<<this->start<<" "<<this->end<<std::endl;
    }

    double compute_pid(double &setpoint,double &pv){
        return this->controller.compute(setpoint,pv);
    }
private:
    // PID
    PID controller = PID(3.0,1.0);

    // members for abundance estimation with unique reads
    double rpk = 0.0;
    uint32_t n_reads = 0;

    // members for abundance estimation with both unique and multimapping reads
    // these members are used to compute the error from the current state (unique+multi) to the target (unique only)
    double rpk_total = 0.0;
    uint32_t n_reads_multi = 0;

    // general purpose members
    uint32_t elen = 1;
    uint32_t start=MAX_INT,end = 0;
    uint8_t strand;
    uint32_t locid=0;
    float abundance = 0.0;
    bool empty = true;
};

class Gene: public Locus{
public:
    Gene() = default;
    Gene(uint32_t geneID,GffObj &p_trans){
        this->set_id(geneID);
        this->set_strand(p_trans.strand);
        this->add_transcript(p_trans);
    }

    // the function below creates an exon map for the effective length calculation
    // as well as updates the start and end of the current gene based on the exon coordinates
    void add_transcript(GffObj &p_trans){
        GList<GffExon>& exon_list = p_trans.exons;
        for (int i = 0; i < exon_list.Count(); ++i) {
            // add exon to the current exon-intron coordinates
            GffExon &cur_exon = *(exon_list.Get(i));
            this->exons.insert(std::make_pair(cur_exon.start,cur_exon.end));
            // update gene start and end coordinates
            if(cur_exon.start<this->get_start()){
                this->set_start(cur_exon.start);
            }
            if(cur_exon.end>this->get_end()){
                this->set_end(cur_exon.end);
            }
        }

    }
    void add_exon(uint32_t start,uint32_t end){
        this->exons.insert(std::make_pair(start,end));
    }
    uint32_t compute_elen(){
        uint32_t cur_elen = 0;
        uint32_t eo = 0;
        std::vector<std::pair<uint32_t,uint32_t> > stack;
        if(this->exons.size()==1){
            return (this->exons.begin()->second+1)-this->exons.begin()->first;
        }
        else if(this->exons.size()>1){
            stack.push_back(*(this->exons.begin()));
            for(auto it=this->exons.begin();it!=this->exons.end();it++){
                eo = exons_overlap(it,stack.back());
                if(eo){ // if the exons overlap, then top needs to be updated
                    stack.back().second = eo;
                }
                else{ // no overlap
                    // update the current effective length based on the previous entry in the stack
                    cur_elen += ((stack.back().second+1)-stack.back().first);
                    stack.push_back(*it);
                }
            }
            // update the current effective length based on the previous entry in the stack
            cur_elen += ((stack.back().second+1)-stack.back().first);
        }
        else{ // no exons in the data??? this is weird
            return cur_elen;
        }

        return cur_elen;
    }

private:
    uint32_t exons_overlap(std::set<std::pair<uint32_t,uint32_t> >::iterator e1,std::pair<uint32_t,uint32_t> e2){
        if(e2.second >= e1->first){ // overlap
            if(e2.second < e1->second){ // need to update the end
                return e1->second;
            }
            else{
                return e2.second;
            }
        }
        else{ // do not overlap
            return 0;
        }
    }
    std::set<std::pair<uint32_t,uint32_t> > exons; // guarantees correct lexicographical sorting
};

// this object holds the map of loci and different means of accessing the contents and modifying it
// for instance for the purpose of Abundance inference
class Loci{
public:
    explicit Loci(uint32_t nloc){
        this->nloc = nloc;
        this->loci = std::vector<Locus>(nloc); // Initialize the index to a given length
    };
    Loci() = default;
    ~Loci()=default;

    Locus& operator[] (const int index){
        return this->loci[index];
    }

    void set_precomputed_abundances(){this->precomputed_abundances = true;}

    void add_locus(uint32_t lid,uint32_t elen){
        this->loci[lid].set_elen(elen);
    }

    void add_read(const int index){
        double rpk_diff = this->loci[index].inc(); // update locus abundance
        this->scaling_factor += (rpk_diff/1000000.0); // update the scaling factor
        this->scaling_factor_total += (rpk_diff/1000000.0); // update the scaling factor
    }

    void add_read_multi(const int index){
        double rpk_diff_total = this->loci[index].inc_multi();
        this->scaling_factor_total += (rpk_diff_total/1000000.0);
    }

    // the following method is used to normalize the locus counts
    double get_abund(uint32_t locid) {
        return this->loci[locid].get_rpk() / this->scaling_factor;
    }

    double get_abund_total(uint32_t locid){
        return this->loci[locid].get_rpk_total()/this->scaling_factor_total;
    }

    double get_corrected(uint32_t locid,double setpoint, double pv){
        // perform PID computation for a given locus
        return this->loci[locid].compute_pid(setpoint,pv);
    }

    // load from a .glst file
    void load(std::string& locus_file){
        struct stat buffer{};
        if(stat (locus_file.c_str(), &buffer) != 0){ // if file does not exists
            std::cerr<<"@ERROR::Locus file not found: "<<locus_file<<std::endl;
            exit(1);
        }
        std::cerr<<"@LOG::loading locus data from: "<<locus_file<<std::endl;
        std::ifstream locfp(locus_file,std::ifstream::binary);
        if(locfp){
            locfp.seekg(0,locfp.end);
            int size = locfp.tellg();
            locfp.seekg(0,locfp.beg);
            char *buffer = new char[size];
            if(locfp.read(buffer,size)){
                this->_load(locfp,buffer);
            }
            delete[] buffer;
            std::cerr<<"@LOG::finished loading the locus data"<<std::endl;
        }
        else{
            std::cerr<<"@ERROR::failed to open the locus file"<<std::endl;
        }
        locfp.close();
    }

    void print(){
        for(auto &v : this->loci){
            v.print();
            std::cout<<this->scaling_factor<<"\t"<<v.get_nreads()<<"\t"<<v.get_rpk()<<"\t"<<this->get_abund(v.get_locid())<<std::endl;
        }
    }

private:
    std::vector<Locus> loci;
    uint32_t nloc;
    double scaling_factor = 1.0,scaling_factor_total = 1.0;
    bool precomputed_abundances = false;

    void _load(std::ifstream& locfp,char* buffer){
        int k = locfp.gcount();
        uint32_t locus=0,start=0,end=0,elen=0;
        enum Opt {LOCUS   = 0,
                  ELEN    = 1,
                  START   = 2,
                  END     = 3};
        uint32_t elem = Opt::LOCUS;
        Locus loc;
        for(int i=0;i<k;i++){
            switch(buffer[i]){
                case '\n':
                    loc.set_end(end);
                    this->loci[loc.get_locid()]=loc;
                    loc.clear();
                    elem = Opt::LOCUS;
                    end = 0;
                    break;
                case ':':
                    loc.set_id(locus);
                    elem = Opt::ELEN;
                    locus = 0;
                    break;
                case ' ':
                    loc.set_start(start);
                    elem = Opt::END;
                    start = 0;
                    break;
                case '-': case '+':
                    loc.set_elen(elen);
                    loc.set_strand((uint8_t)buffer[i]);
                    elem = Opt::START;
                    elen = 0;
                    break;
                case '0': case '1': case '2': case '3':
                case '4': case '5': case '6': case '7':
                case '8': case '9':
                    //add to the current integer
                    switch(elem){
                        case 0: //chromosome
                            locus = 10*locus + buffer[i] - '0';
                            break;
                        case 1:
                            elen = 10*elen + buffer[i] - '0';
                            break;
                        case 2:
                            start = 10*start + buffer[i] - '0';
                            break;
                        case 3:
                            end = 10*end + buffer[i] - '0';
                            break;
                        default:
                            std::cerr<<"@ERROR::should never happen"<<std::endl;
                            exit(1);
                    }
                    break;
                default:
                    std::cerr<<"@ERROR::unrecognized character"<<buffer[i]<<std::endl;
                    exit(1);
            }
        }
        if(!loc.is_empty()){ // no end of line character, so need to write the last entry
            loc.set_end(end);
            loc.set_elen((loc.get_end()+1)-loc.get_start());
            this->loci[loc.get_locid()]=loc;
        }
    }
};

class Position{
public:
    Position() = default;
    Position(uint32_t chr,uint32_t strand,uint32_t start,uint32_t locus,uint32_t transID){
        this->chr = chr;
        this->strand = strand;
        this->start = start;
        this->locus = locus;
        this->transID = transID;
    }
    Position(bam1_t* curAl){ // generate position from an alignment
        this->chr = (uint32_t)curAl->core.tid;
        this->start = (uint32_t)curAl->core.pos;
    }
    ~Position() = default;

    void add_move(uint32_t move){this->moves.push_back(move);this->num_elems++;}
    void set_chr(uint32_t chr){this->chr=chr;this->num_elems++;}
    void set_strand(uint32_t strand){this->strand=strand;this->num_elems++;}
    void set_start(uint32_t start){this->start=start;this->num_elems++;}
    void set_locus(uint32_t locus){this->locus=locus;this->num_elems++;}
    void set_trans(uint32_t trans){this->transID=trans;this->num_elems++;}

    std::string get_strg() const {
        std::string res;
        res.append(std::to_string(this->transID));
        res += ">";
        res.append(std::to_string(this->locus));
        res += "@";
        res.append(std::to_string(this->chr));
        res += this->strand;
        res.append(std::to_string(this->start));
        res += ':';
        for(auto &mit : this->moves){
            res.append(std::to_string(mit));
            res += ' ';
        }
        res.pop_back(); // removes the last whitespace
        return res;
    }

    bool operator==(const Position& m) const{
        return this->chr==m.chr &&
               this->strand==m.strand &&
               this->start==m.start &&
               this->locus==m.locus &&
               this->moves==m.moves;
    }

    bool operator>(const Position& m) const{
        return this->locus > m.locus;
    }

    bool operator<(const Position& m) const{
        return this->chr<m.chr ||
               this->strand<m.strand ||
               this->start<m.start ||
               this->locus<m.locus ||
               this->moves<m.moves;
    }

    void clear(){
        this->moves.clear();
        this->num_elems = 0;
    }

    uint8_t size(){return this->num_elems;}

    uint8_t num_elems = 0;
    uint32_t chr,strand,start,locus;
    uint32_t transID; // transID can be any transcript which can describe the position to which it belongs; it does not participate in the equality computation
    std::vector<uint32_t> moves; // simplified CIGAR describing the intron-exon coverage of the given kmer
};

inline void hash_combine(std::size_t& seed) { }

template <typename T, typename... Rest>
inline void hash_combine(std::size_t& seed, const T& v, Rest... rest) {
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
    hash_combine(seed, rest...);
}


// hash function to be used for
namespace std {
    template<>
    struct hash<Position> {
        size_t operator()(const Position &p) const {
            size_t ret = 0;
            hash_combine(ret,p.chr,p.strand,p.start,p.locus);
            return ret;
        }
    };
}

struct GffTranscript: public GSeg {
    GVec<GSeg> exons;
    uint32_t gffID;
    uint32_t refID;
    uint32_t geneID;
    uint32_t numID;
    uint8_t strand;
    uint32_t abundance; // only used if the salmon or kallisto quantifications are provided
    bool empty = true;
    GffTranscript():exons(1), numID(-1), gffID(),
                    refID(), strand(0), geneID(0),abundance(0) { }

    void set_gffID(uint32_t gffID){this->gffID=gffID;this->empty=false;}
    void set_refID(uint32_t refID){this->refID=refID;this->empty=false;}
    void set_geneID(uint32_t geneID){this->geneID=geneID;this->empty=false;}
    void set_numID(uint32_t numID){this->numID=numID;this->empty=false;}
    void set_strand(uint8_t strand){this->strand=strand;this->empty=false;}
    void add_exon(GSeg exon){ this->exons.Add(exon);this->empty=false;}

    uint32_t get_refID() {return refID;}
    uint32_t get_geneID(){return geneID;}

    void clear(){
        this->exons.Clear();
        this->empty = true;
    }
    void print(){
        std::cout<<this->gffID<<"\t"<<this->geneID<<"@"<<this->refID<<this->strand;
        for(int i=0;i<this->exons.Count();i++){
            std::cout<<this->exons[i].start<<"_"<<this->exons[i].end<<",";
        }
        std::cout<<std::endl;
    }
};

// this object holds information regarding the fragment length distribution
// this is then used to tell if a multimapping pair is valid
class Fragments{
public:
    Fragments() = default;
    ~Fragments() = default;

    void add_pair(bam1_t *al,bam1_t *mate){
        // calculate the distance from the first start - the mate start and add to the distribution
        cur_dist = get_dist(al,mate);
        min = std::min(this->min,cur_dist);
        max = std::max(this->max,cur_dist);
    }

    int get_min(){return this->min;}
    int get_max(){return this->max;}

    bool is_valid(bam1_t *al, bam1_t *mate){ // checks if the current pair passes the fragment length test // TODO: make better - currently will test if within bounds that have been observed - perhaps compute the possible ends of the distribution and evaluate the likelihood?
        cur_dist = get_dist(al,mate);
        return (cur_dist>=min && cur_dist<=max);
    }
private:
    int min = 0, max = MAX_INT;
    int cur_dist;
    uint32_t get_dist(bam1_t *al,bam1_t *mate){
        return std::abs(al->core.mpos - al->core.pos);
    }
};

// this class describes the multimappers in the index transcriptome
// and facilitates efficient storage and lookup
class Multimap{
public:
    Multimap() = default;
    ~Multimap() = default;

    void set_fraglen(int fraglen){
        this->fraglen = fraglen;
    }

    void set_num_multi(int num_multi){
        this->num_multi = num_multi;
    }

    void set_all_multi(){
        this->all_multi = true;
    }

    void set_loci(Loci* loci){
        this->loci = loci;
    }

    void set_transcriptome(std::vector<GffTranscript>* transcriptome){ // used to find loaded abundancs when available
        this->transcriptome = transcriptome;
    }

    void set_precomputed_abundances(){
        this->precomputed_abundances = true;
    }

    void add_frag(bam1_t* curAl,bam1_t* mate){
        this->frags.add_pair(curAl,mate);
    }

    int get_min_frag(){return this->frags.get_min();}
    int get_max_frag(){return this->frags.get_max();}

    Loci* loci;
    std::vector<GffTranscript>* transcriptome;

    int get_block_size(){
        int count = 0;
        while(this->ii->second){
            count++;
            this->ii++;
        }
        this->ii-=count; // reset multimappers to the beginning
        return count;
    }

    // copies the multimappers of the block pointed to by the current iterator
    void copy_current(std::vector<Position>& pos_res){
        while(this->ii->second){ // iterate until the end of the block
            pos_res.push_back(this->ii->first);
            this->ii++; // step to the next position
        }
        pos_res.push_back(this->ii->first);
    }

    // TODO: need to refactor the process_pos functions
    //     mainly would prefer having a single function that handles all cases and is no longer redundant
    int process_pos_precomp(Position& pos,Loci& loci,std::vector<Position>& pos_res){ // using precomputed abundance values only
        this->ltf = this->lookup_table.find(pos);
        if(this->ltf == this->lookup_table.end()){ // position is not a multimapper
            return 0;
        }
        else{ // multimappers exist - need to evaluate
            // if all mode - need to keep a bool and then just iterate through the entire block and output all values into a vector of positions which can be accessed by the converter and written to the alignment
            std::vector<double>abunds,res; // holds pid-corrected abundances
            double total = 0;

            this->ii = this->index.begin()+this->ltf->second; // get iterator to the start of a multimapping block in the array
            int cur_num_multi = get_block_size();
            if(this->all_multi){ // just output the entire block
                copy_current(pos_res);
                return cur_num_multi; // done
            }
            while(this->ii->second){ // iterate until the end of the block
                // first need to compute expected values
                // for this we need the uniq abundances which determine base likelihood
                // as well as current assignments of multimappers
                abunds.push_back(loci[this->ii->first.locus].get_abund());
                total += abunds.back();

                this->ii++; // step to the next position
            }
            abunds.push_back(loci[this->ii->first.locus].get_abund());
            total += abunds.back();

            for(auto &v : abunds){
                res.push_back(v/total);
            }

            // now need to make the decision which is best
            int pos_idx = this->get_likely(res,0.0,1.0);
            // now follow the iterator to get the actual position object which corresponds to the selected item
            pos_res.push_back((this->index.begin()+this->ltf->second+pos_idx)->first);
            return cur_num_multi;
        }
    }

    // given a position generated from an alignment this function searches for multimapper
    // and returns true if the position is not multimapping or false if it is multimapping
    // for multimappers, it also replaces the data in the position with a position selected by likelihood
    int process_pos(Position& pos,Loci& loci,std::vector<Position>& pos_res){
        this->ltf = this->lookup_table.find(pos);
        if(this->ltf == this->lookup_table.end()){ // position is not a multimapper
            return 0;
        }
        else{ // multimappers exist - need to evaluate
            std::vector<double> uniq,multi,tmp_res,res; // holds pid-corrected abundances
            double utotal=0,mtotal=0,rtotal=0;

            this->ii = this->index.begin()+this->ltf->second; // get iterator to the start of a multimapping block in the array
            int offset = this->ltf->second;
            int cur_num_multi = get_block_size();
            if(this->all_multi){ // just output the entire block
                copy_current(pos_res);
                return cur_num_multi; // done
            }
            while(this->ii->second){ // iterate until the end of the block
                // first need to compute expected values
                // for this we need the uniq abundances which determine base likelihood
                // as well as current assignments of multimappers
                uniq.push_back(loci.get_abund(this->ii->first.locus));
                utotal+=uniq.back();
                multi.push_back(loci.get_abund_total(this->ii->first.locus));
                mtotal+=multi.back();

                this->ii++; // step to the next position
            }
            uniq.push_back(loci.get_abund(this->ii->first.locus));
            utotal+=uniq.back();
            multi.push_back(loci.get_abund_total(this->ii->first.locus));
            mtotal+=multi.back();

            // get expected and compute PID
            utotal = (utotal==0)?1:utotal;
            mtotal = (mtotal==0)?1:mtotal;
            double min = (double)MAX_INT;
            for(int i=0; i<uniq.size();i++){
                double expected = (uniq[i]/utotal)*mtotal;
                double tmp = uniq[i]+loci.get_corrected(this->ii->first.locus,expected+uniq[i],multi[i]+uniq[i]);
                tmp = (tmp>0.0)?tmp:0.000001;
                tmp_res.push_back(tmp);
                min = (tmp<min)?tmp:min;
                rtotal+=tmp_res.back();
            }

            for(auto &v : tmp_res){
                res.push_back(v/rtotal);
            }

            // now need to make the decision which is best
            int pos_idx = this->get_likely(res,0.0,1.0);
            // now follow the iterator to get the actual position object which corresponds to the selected item
            pos_res.push_back((this->index.begin()+this->ltf->second+pos_idx)->first);
            return cur_num_multi;
        }
    }

    int process_pos_pair_precomp(Position& pos,Position& pos_mate,Loci& loci,std::vector<Position>& pos_res,std::vector<Position>& pos_res_mate){ // using precomputed abundance values only
        this->ltf = this->lookup_table.find(pos);
        if(this->ltf == this->lookup_table.end()){ // position is not a multimapper
            return 0;
        }
        this->ltf_mate = this->lookup_table.find(pos_mate);
        if(this->ltf_mate == this->lookup_table.end()){ // position is not a multimapper
            return 0;
        }

        else{ // multimappers exist - need to evaluate
            std::vector<double> abunds,res; // holds pid-corrected abundances
            double total=0;

            this->ii = this->index.begin()+this->ltf->second; // get iterator to the start of a multimapping block in the array
            this->ii_mate = this->index.begin()+this->ltf_mate->second; // get iterator to the start of a multimapping block in the array
            // here we rely on the sorted order of the multimappers in the index, which guarantees that multimappers on the same locus will appear close by
            int idx=0,idx_mate=0; // indices help us keep track of the positions within multimapper blocks
            int locus_index = 0; // indicates the position within block of the current locus
            std::vector<std::pair<int,int>> multi_pairs; // holds offsets of the correct paired positions
            while(true){
                this->ii_mate = this->index.begin()+this->ltf_mate->second+locus_index; // reset iterator to the beginning
                idx_mate = locus_index;
                bool locus_found = false; // did the inner loop find the locus of the outer loop
                while(true){
                    if(this->ii->first.locus == this->ii_mate->first.locus &&
                       std::abs((int)this->ii->first.start - (int)this->ii_mate->first.start) < this->fraglen){ // valid pair
                        multi_pairs.push_back(std::make_pair(idx,idx_mate));
                        locus_found = true;
                    }
                    if(this->ii->first.locus != this->ii_mate->first.locus){
                        if(locus_found){ // exceeded current locus - can now terminate
                            break;
                        }
                        else { // can reset locus index to skip iterations in future
                            locus_index++;
                        }
                    }
                    if(!this->ii_mate->second){
                        break;
                    }
                    idx_mate++;
                    this->ii_mate++;
                }
                if(!this->ii->second){
                    break;
                }
                idx++;
                this->ii++;
            }
            if(multi_pairs.empty()){
                return 0; // means no valid multimapping pairs were detected
            }
//            std::cerr<<"multi exists pair"<<std::endl;
            int cur_num_multi = multi_pairs.size();
            // now compute abundances for all these blocks
            this->ii = this->index.begin()+this->ltf->second; // return to the start of the multimapping block
            this->ii_mate = this->index.begin()+this->ltf_mate->second; // same for the mate
            if(this->all_multi){ // just output the entire block
                for(auto & mp: multi_pairs) { // iterate over all the detected pairs
                    pos_res.push_back((this->ii + mp.first)->first);
                    pos_res_mate.push_back((this->ii + mp.second)->first);
                }
                return cur_num_multi; // done
            }
            for(auto & mp: multi_pairs){
                // first need to compute expected values
                // for this we need the uniq abundances which determine base likelihood
                // as well as current assignments of multimappers
                abunds.push_back(loci[(this->ii+mp.first)->first.locus].get_abund());
                total += abunds.back();
            }

            for(auto &v : abunds){
                res.push_back(v/total);
            }

            // now need to make the decision which is best
            int pos_idx = this->get_likely(res,0.0,1.0);
            // now follow the iterator to get the actual position object which corresponds to the selected item

            int offset = multi_pairs[pos_idx].first; // get the actual offset
            int offset_mate = multi_pairs[pos_idx].second;
            pos_res.push_back((this->index.begin()+this->ltf->second+offset)->first);
            pos_res_mate.push_back((this->index.begin()+this->ltf_mate->second+offset_mate)->first);
            return cur_num_multi;
        }
    }

    int process_pos_pair(Position& pos,Position& pos_mate,Loci& loci,std::vector<Position>& pos_res,std::vector<Position>& pos_res_mate){
        this->ltf = this->lookup_table.find(pos);
        if(this->ltf == this->lookup_table.end()){ // position is not a multimapper
            return 0;
        }
        this->ltf_mate = this->lookup_table.find(pos_mate);
        if(this->ltf_mate == this->lookup_table.end()){ // position is not a multimapper
            return 0;
        }

        else{ // multimappers exist - need to evaluate
            std::vector<double> uniq,multi,tmp_res,res; // holds pid-corrected abundances
            double utotal=0,mtotal=0,rtotal=0;

            this->ii = this->index.begin()+this->ltf->second; // get iterator to the start of a multimapping block in the array
            this->ii_mate = this->index.begin()+this->ltf_mate->second; // get iterator to the start of a multimapping block in the array
            // here we rely on the sorted order of the multimappers in the index, which guarantees that multimappers on the same locus will appear close by
            int idx=0,idx_mate=0; // indices help us keep track of the positions within multimapper blocks
            int locus_index = 0; // indicates the position within block of the current locus
            std::vector<std::pair<int,int>> multi_pairs; // holds offsets of the correct paired positions
            while(true){
                this->ii_mate = this->index.begin()+this->ltf_mate->second+locus_index; // reset iterator to the beginning
                idx_mate = locus_index;
                bool locus_found = false; // did the inner loop find the locus of the outer loop
                while(true){
                    if(this->ii->first.locus == this->ii_mate->first.locus &&
                      std::abs((int)this->ii->first.start - (int)this->ii_mate->first.start) < this->fraglen){ // valid pair // TODO: set fragment length based on the precomputed distribution
                        multi_pairs.push_back(std::make_pair(idx,idx_mate));
                        locus_found = true;
                    }
                    if(this->ii->first.locus != this->ii_mate->first.locus){
                        if(locus_found){ // exceeded current locus - can now terminate
                            break;
                        }
                        else { // can reset locus index to skip iterations in future
                            locus_index++;
                        }
                    }
                    if(!this->ii_mate->second){
                        break;
                    }
                    idx_mate++;
                    this->ii_mate++;
                }
                if(!this->ii->second){
                    break;
                }
                idx++;
                this->ii++;
            }
            if(multi_pairs.empty()){
                return 0; // means no valid multimapping pairs were detected
            }
//            std::cerr<<"multi exists pair"<<std::endl;
            int cur_num_multi = multi_pairs.size();
            // now compute abundances for all these blocks
            this->ii = this->index.begin()+this->ltf->second; // return to the start of the multimapping block
            this->ii_mate = this->index.begin()+this->ltf_mate->second; // same for the mate
            if(this->all_multi){ // just output the entire block
                for(auto & mp: multi_pairs) { // iterate over all the detected pairs
                    pos_res.push_back((this->ii + mp.first)->first);
                    pos_res_mate.push_back((this->ii_mate + mp.second)->first);
                }
                return cur_num_multi; // done
            }
            for(auto & mp: multi_pairs){
                // first need to compute expected values
                // for this we need the uniq abundances which determine base likelihood
                // as well as current assignments of multimappers
                uniq.push_back(loci.get_abund((this->ii+mp.first)->first.locus));
                utotal+=uniq.back();
                multi.push_back(loci.get_abund_total((this->ii+mp.first)->first.locus));
                mtotal+=multi.back();
            }

            // get expected and compute PID
            utotal = (utotal==0)?1:utotal;
            mtotal = (mtotal==0)?1:mtotal;
            double min = (double)MAX_INT;
            for(int i=0; i<uniq.size();i++){
                double expected = (uniq[i]/utotal)*mtotal;
                double tmp = uniq[i]+loci.get_corrected(this->ii->first.locus,expected+uniq[i],multi[i]+uniq[i]);
                tmp = (tmp>0.0)?tmp:0.000001;
                tmp_res.push_back(tmp);
                min = (tmp<min)?tmp:min;
                rtotal+=tmp_res.back();
            }

            for(auto &v : tmp_res){
                res.push_back(v/rtotal);
            }

            // now need to make the decision which is best
            int pos_idx = this->get_likely(res,0.0,1.0);
            // now follow the iterator to get the actual position object which corresponds to the selected item

            int offset = multi_pairs[pos_idx].first; // get the actual offset
            int offset_mate = multi_pairs[pos_idx].second;
            pos_res.push_back((this->index.begin()+this->ltf->second+offset)->first);
            pos_res_mate.push_back((this->index.begin()+this->ltf_mate->second+offset_mate)->first);
            return cur_num_multi;
        }
    }

    // thanks to https://medium.com/@dimonasdf/hi-your-code-is-correct-but-in-pseudocode-it-should-be-fc1875cf9de3 for the idea
    int get_likely(std::vector<double>& abunds,double min, double max){
        double rv = (max - min) * ( (double)rand() / (double)RAND_MAX ) + min;

        for(int i=0;i<abunds.size();i++){
            rv = rv - abunds[i];
            if(rv<=0){
                return i;
            }
        }
        return abunds.size()-1;
    }

    void load(const std::string& input_file){
        srand(time(0));
        std::cerr<<"@LOG::loading multimapper data from: "<<input_file<<std::endl;
        std::ifstream infp(input_file,std::ifstream::binary);
        if(infp){
            infp.seekg(0,infp.end);
            int size = infp.tellg();
            infp.seekg(0,infp.beg);
            char *buffer = new char[size];
            if(infp.read(buffer,size)){
                this->_load(infp,buffer);
            }
            delete[] buffer;
            std::cerr<<"@LOG::finished loading the multimapper data"<<std::endl;
        }
        else{
            std::cerr<<"@LOG::failed to open the multimapping file"<<std::endl;
        }
        infp.close();
    }

    // print using Positions from the kmer_loadcoords
    void print(){
        for(auto &kv : this->kmer_coords){
            std::cout<<kv.first<<"\t";
            for(auto &cv : kv.second){
                std::cout<<cv.get_strg()<<"\t";
            }
            std::cout<<std::endl;
        }
    }

    // print from the loaded index
    void print_multimapers(){
        for(auto &mm : this->index){
            if(mm.second){
                std::cout<<mm.first.get_strg()<<"\t";
            }
            else{
                std::cout<<mm.first.get_strg()<<std::endl;
            }
        }
    }

    void set_kmerlen(int kmerlen){this->kmerlen = kmerlen;}

    void save_multimappers(std::string& outFP){
        std::cerr<<"@LOG::writing multimappers"<<std::endl;
        std::ofstream multi_ss(outFP.c_str());
        std::string res = "";
        for(auto &kv : this->kmer_coords){
            if(kv.second.size()>1) {
                // pre-sort coordinates here
                std::vector<Position> tmp;
                tmp.insert(tmp.end(), kv.second.begin(), kv.second.end());
                std::sort(tmp.begin(), tmp.end(),[](Position const& lhs, Position const& rhs) { return lhs > rhs; });

                for(auto &cv : tmp){
                    res.append(cv.get_strg());
                    res+='\t';
                }

                res.pop_back();
                res+='\n';
                multi_ss << res;
                multi_ss.flush();
                res = "";
            }
        }

        multi_ss.close();
        std::cerr<<"@LOG::done writing multimappers"<<std::endl;

    }

    void save_unique(std::string& outFP){
        std::ofstream uniq_ss(outFP.c_str());
        for(auto &kv : this->kmer_coords){
            if(kv.second.size()==1) {
                for (auto &cv : kv.second) {
                    uniq_ss << cv.get_strg() << "\t";
                }
                uniq_ss << std::endl;
            }
        }
        uniq_ss.close();
    }

    // given a full transcript sequence and the pointer to the transcript description (constituent exons)
    // get all kmers and coordinates and load them into the multimapper index
    void add_sequence(std::string& seq,GffObj &p_trans,uint32_t geneID,int transID){
        GList<GffExon>& exon_list = p_trans.exons;

        int el_pos = 0; // current position within the exon list
        int exon_pos = 0;
        std::string kmer = ""; // currently evaluated kmer
        GffExon *cur_exon = exon_list[el_pos];
        for(int i=0;i<seq.size()-this->kmerlen+1;i++){ // iterate over all kmers in the sequence
            kmer=seq.substr(i,this->kmerlen);

            Position p(p_trans.gseq_id,(uint32_t)p_trans.strand,cur_exon->start+exon_pos,geneID,transID); // initialize position and add information to it accordingly

            this->kce = this->kmer_coords.insert(std::make_pair(kmer,pcord{}));

            process_remaining(el_pos,exon_pos,p,exon_list,cur_exon); // adds moves to the position object

            this->pce = this->kce.first->second.insert(p); // insert new completed position now

            // reset_parameters
            exon_pos++; // increment and evaluate
            if(exon_pos+cur_exon->start > cur_exon->end){
                exon_pos = 0;
                el_pos++;
                cur_exon = exon_list[el_pos];
            }
        }
    }

private:
    // FRAGMENT LENGTH STUFF
    Fragments frags;

    // we know that the kmer spans an intron
    // this function appends the intron/exon information to moves of the position object
    void process_remaining(int& el_pos,int& exon_pos,Position& p,GList<GffExon>& exon_list,GffExon* cur_exon){
        int ee = cur_exon->end;
        int es = cur_exon->start;
        int el = (ee+1) - es; // length of the current exon
        if(el-exon_pos >= this->kmerlen){ // kmer fits well
            p.add_move(this->kmerlen);
            return;
        }
        else{
            int left = this->kmerlen - (el-exon_pos);
            p.add_move(el-exon_pos);

            GffExon *next_exon,*prev_exon;
            int nee,nes,nel,cee,ces,cel; // holders for current and next exons
            for(int i=el_pos+1;i<exon_list.Count();i++){ // iterate over the exons
                prev_exon = exon_list[i-1];
                ee = prev_exon->end,es = prev_exon->start,el = (ee+1)-es;
                next_exon = exon_list[i];
                nee = next_exon->end,nes = next_exon->start,nel = (nee+1)-nes;

                p.add_move(nes - ee); // adding the intron information to the moves
                if(nel >= left){ // remainder fits well
                    p.add_move(left);
                    return;
                }
                else{
                    left = left - nel;
                    p.add_move(nel);
                }
            }
            std::cerr<<"@ERROR::exon boundaries exceeded"<<std::endl;
            exit(1);
        }
    }

    struct coord_hash { // in case it is implemented as an unordered_map
        uint64_t operator()(const std::vector<uint32_t> &coords ) const{
            uint64_t resHash=1;
            for (auto const& coord: coords) {
                resHash=resHash^std::hash<uint64_t>()(coord);
            }
            return resHash;
        }
    };

    void _load(std::ifstream &infp,char *buffer){
        int k = infp.gcount();
        uint32_t cur=0,chr=0,strand=0,start=0,move=0,locus=0,trans=0;
        enum Opt {CHR   = 0,
                START = 1,
                MOVE  = 2,
                LOCUS = 3,
                TRANS = 4};
        uint32_t elem = Opt::TRANS;
        Position pos;
        // Multimappers are all stored on a single line, so instead of creating pointers for each, we can simply store two small integers, which describe which positions in the index characterize a multimapper
        uint32_t multi_block_start = 0;
        for(int i=0;i<k;i++){
            switch(buffer[i]){
                case '\n':
                    //end of line
                    pos.add_move(move);
                    this->index.push_back(std::make_pair(pos,false));
                    this->lte = this->lookup_table.insert(std::make_pair(pos,multi_block_start));
                    multi_block_start = this->index.size(); // end of line implies end of the multimapping block. Now set the start of the next block
                    pos.clear();
                    elem = Opt::TRANS;
                    move = 0;
                    break;
                case '\t':
                    // end of a full coordinate - write last integer
                    pos.add_move(move);
                    this->index.push_back(std::make_pair(pos,true));
                    this->lte = this->lookup_table.insert(std::make_pair(pos,multi_block_start));
                    pos.clear();
                    elem = Opt::TRANS;
                    move = 0;
                    break;
                case ' ':
                    //end of current int - write move
                    pos.add_move(move);
                    elem = Opt::MOVE;
                    move = 0;
                    break;
                case ':':
                    //end of current int - write start
                    pos.set_start(start);
                    elem = Opt::MOVE;
                    start = 0;
                    break;
                case '@':
                    // got the geneID
                    pos.set_locus(locus);
                    elem = Opt::CHR;
                    locus = 0;
                    break;
                case '>':
                    // got the sample transcript ID
                    pos.set_trans(trans);
                    elem = Opt::LOCUS;
                    trans = 0;
                    break;
                case '-': case '+':
                    // negative strand - write chromosome
                    pos.set_chr(chr);
                    pos.set_strand((uint32_t)buffer[i]);
                    elem = Opt::START;
                    chr = 0;
                    break;
                case '0': case '1': case '2': case '3':
                case '4': case '5': case '6': case '7':
                case '8': case '9':
                    //add to the current integer
                    switch(elem){
                        case 0: //chromosome
                            chr = 10*chr + buffer[i] - '0';
                            break;
                        case 1:
                            start = 10*start + buffer[i] - '0';
                            break;
                        case 2:
                            move = 10*move + buffer[i] - '0';
                            break;
                        case 3:
                            locus = 10*locus + buffer[i] - '0';
                            break;
                        case 4:
                            trans = 10*trans + buffer[i] - '0';
                            break;
                        default:
                            std::cerr<<"@ERROR::should never happen _load from Multimap with value: "<<elem<<std::endl;
                            exit(1);
                    }
                    break;
                default:
                    std::cerr<<"@ERROR:unrecognized character: "<<buffer[i]<<std::endl;
                    exit(1);
            }
        }
        if(pos.size() >= 4){ // no end of line character, so need to write the last entry
            this->index.push_back(std::make_pair(pos,false));
            this->lte = this->lookup_table.insert(std::make_pair(pos,multi_block_start));
        }
    }

    int kmerlen;
    int fraglen = 200000; // TODO: URGENT!
    int num_multi = 1;
    int all_multi = false;
    bool precomputed_abundances = false;

    typedef std::unordered_set<Position> pcord;
    std::pair<pcord::iterator,bool> pce;
    std::unordered_map<std::string,pcord> kmer_coords;
    std::pair<std::unordered_map<std::string,pcord>::iterator,bool> kce;

    // now time to write the unique kmers for transcripts
    std::unordered_map<std::string,int> uniq_cnt; // counts of unique kmers per transcript
    std::pair<std::unordered_map<std::string,int>::iterator,bool> ex_ucnt; // exists or not

    std::unordered_map<Position,uint32_t> lookup_table;
    std::pair<std::unordered_map<Position,uint32_t>::iterator,bool> lte; // entry exists in the lookup table
    std::unordered_map<Position,uint32_t>::iterator ltf,ltf_mate; // iterator to found entry or the end
    std::vector<std::pair<Position,bool>> index; // the actual position as well as whether the next position is related to the current
    std::vector<std::pair<Position,bool>>::iterator ii,ii_mate; // iterator within the index
    // every entry in the lookup table links to the first entry in the block of related multimappers

};

#endif //TRANS2GENOME_MULTIMAP_H
