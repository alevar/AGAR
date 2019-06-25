//
// Created by sparrow on 6/16/19.
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

#include "gff.h"
#include <htslib/sam.h>

class PID{
public:
    PID() = default;
    PID(double kp, double ki){
        // set tunings here
        if (kp<0 || ki<0){
            std::cerr<<"wrong PID tunings specified"<<std::endl;
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
        this->I_term += (this->ki * error);

        /*Compute PID Output*/
        double res = (this->kp * error) + this->I_term;

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

    uint32_t get_start(){return this->start;}
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
    PID controller = PID(2.0,1.0);

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
    double get_abund(uint32_t locid){
        return this->loci[locid].get_rpk()/this->scaling_factor;
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
            std::cerr<<"Locus file: "<<locus_file<<" is not found. Check that correct index is provided."<<std::endl;
            exit(1);
        }
        std::cerr<<"loading locus data from: "<<locus_file<<std::endl;
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
            std::cerr<<"finished loading the locus data"<<std::endl;
        }
        else{
            std::cerr<<"failed to open the locus file"<<std::endl;
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
                            std::cerr<<"should never happen"<<std::endl;
                            exit(1);
                    }
                    break;
                default:
                    std::cerr<<"unrecognized character"<<std::endl;
                    std::cerr<<buffer[i]<<std::endl;
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

// hash function to be used for
namespace std {
    template<>
    struct hash<Position> {
        size_t operator()(const Position &p) const {
            int hash = 1;
            hash ^= p.chr + 0x9e3779b9 + (hash<<6) + (hash>>2);
            hash ^= p.strand + 0x9e3779b9 + (hash<<6) + (hash>>2);
            hash ^= p.start + 0x9e3779b9 + (hash<<6) + (hash>>2);
            hash ^= p.locus + 0x9e3779b9 + (hash<<6) + (hash>>2);
            for(auto &v : p.moves){
                hash ^= v + 0x9e3779b9 + (hash<<6) + (hash>>2);
            }
            return hash;
        }
    };
}

// this class describes the multimappers in the index transcriptome
// and facilitates efficient storage and lookup
class Multimap{
public:
    Multimap() = default;
    ~Multimap() = default;

    // given a position generated from an alignment this function searches for multimapper
    // and returns true if the position is not multimapping or false if it is multimapping
    // for multimappers, it also replaces the data in the position with a position selected by likelihood
    bool process_pos(Position& pos,Loci& loci){
        this->ltf = this->lookup_table.find(pos);
        if(this->ltf == this->lookup_table.end()){ // position is not a multimapper
            return true;
        }
        else{ // multimappers exist - need to evaluate
            std::vector<double> uniq,multi,tmp_res,res; // holds pid-corrected abundances
            double utotal=0,mtotal=0,rtotal=0;

            this->ii = this->index.begin()+this->ltf->second; // get iterator to the start of a multimapping block in the array
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

            // now need to make the decision which is best
            int pos_idx = this->get_likely(tmp_res,min,rtotal);
            // now follow the iterator to get the actual position object which corresponds to the selected item
            pos = (this->index.begin()+this->ltf->second+pos_idx)->first;

            return false;
        }
    }

    // thanks to https://medium.com/@dimonasdf/hi-your-code-is-correct-but-in-pseudocode-it-should-be-fc1875cf9de3 for the idea
    int get_likely(std::vector<double>& abunds,double& min, double& max){
        srand(time(NULL));
        double rv = (max - min) * ( (double)rand() / (double)RAND_MAX ) + min;

        for(int i=0;i<abunds.size();i++){
            rv = rv - abunds[i];
            if(rv<=0){
                return i;
            }
        }
    }

    void load(const std::string& input_file){
        std::cerr<<"loading multimapper data from: "<<input_file<<std::endl;
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
            std::cerr<<"finished loading the multimapper data"<<std::endl;
        }
        else{
            std::cerr<<"failed to open the multimapping file"<<std::endl;
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
        std::cerr<<"NUMBER OF MULTIMAPPING POSITIONS: "<<this->index.size()<<std::endl;
        for(auto &mm : this->index){
            if(mm.second){
                std::cerr<<mm.first.get_strg()<<"\t";
            }
            else{
                std::cerr<<mm.first.get_strg()<<std::endl;
            }
        }
    }

    void set_kmerlen(int kmerlen){this->kmerlen = kmerlen;}

    void save_multimappers(std::string& outFP){
        std::string res = "";
        for(auto &kv : this->kmer_coords){
            if(kv.second.size()>1) {
                for (auto &cv : kv.second) {
                    res.append(cv.get_strg());
                    res+='\t';
                }
                res.pop_back();
                res+='\n';
            }
        }
        std::ofstream multi_ss(outFP.c_str());
        multi_ss << res;
        multi_ss.close();
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
            std::cerr<<"exon boundaries exceeded"<<std::endl;
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
//                    std::cerr<<"start: "<<start<<std::endl;
                    pos.set_start(start);
                    elem = Opt::MOVE;
                    start = 0;
                    break;
                case '@':
                    // got the geneID
//                    std::cerr<<"locus: "<<locus<<std::endl;
                    pos.set_locus(locus);
                    elem = Opt::CHR;
                    locus = 0;
                    break;
                case '>':
                    // got the sample transcript ID
//                    std::cerr<<"trans: "<<trans<<std::endl;
                    pos.set_trans(trans);
                    elem = Opt::LOCUS;
                    trans = 0;
                    break;
                case '-': case '+':
                    // negative strand - write chromosome
//                    std::cerr<<"chr: "<<chr<<std::endl;
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
                            std::cerr<<"should never happen _load from Multimap with value: "<<elem<<std::endl;
                            exit(1);
                    }
                    break;
                default:
                    std::cerr<<"unrecognized character"<<std::endl;
                    std::cerr<<buffer[i]<<std::endl;
                    exit(1);
            }
        }
        if(pos.size() >= 4){ // no end of line character, so need to write the last entry
            this->index.push_back(std::make_pair(pos,false));
            this->lte = this->lookup_table.insert(std::make_pair(pos,multi_block_start));
        }
    }

    int kmerlen;

    typedef std::unordered_set<Position> pcord;
    std::pair<pcord::iterator,bool> pce;
    std::unordered_map<std::string,pcord> kmer_coords;
    std::pair<std::unordered_map<std::string,pcord>::iterator,bool> kce;

    // now time to write the unique kmers for transcripts
    std::unordered_map<std::string,int> uniq_cnt; // counts of unique kmers per transcript
    std::pair<std::unordered_map<std::string,int>::iterator,bool> ex_ucnt; // exists or not

    std::unordered_map<Position,uint32_t> lookup_table;
    std::pair<std::unordered_map<Position,uint32_t>::iterator,bool> lte; // entry exists in the lookup table
    std::unordered_map<Position,uint32_t>::iterator ltf; // iterator to found entry or the end
    std::vector<std::pair<Position,bool>> index; // the actual position as well as whether the next position is related to the current
    std::vector<std::pair<Position,bool>>::iterator ii; // iterator within the index
    // every entry in the lookup table links to the first entry in the block of related multimappers

};

#endif //TRANS2GENOME_MULTIMAP_H
