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

#include "gff.h"

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

    void set_elen(uint32_t elen){
        this->elen = elen;
    }

    // return rpk difference which is used to increment the normalizing constant
    double inc(){ // increments the read count and recomputes the TPM
        this->n_reads++;
        double old_rpk = this->rpk;
        this->rpk = this->n_reads/this->elen;
        return this->rpk - old_rpk;
    }
    double get_rpk(){ // return the RPK
        return this->rpk;
    }
private:
    uint32_t elen;
    double rpk;
    uint32_t n_reads;
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

    void add_locus(uint32_t lid,uint32_t elen){
        this->loci[lid].set_elen(elen);
    }

    void add_read(uint32_t lid){
        double rpk_diff = this->loci[lid].inc(); // update locus abundance
        this->scaling_factor += rpk_diff/1000000; // update the scaling factor
    }

    // the following method is used to normalize the locus counts
    double get_abund(uint32_t locid){
        return this->loci[locid].get_rpk()/this->scaling_factor;
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

private:
    std::vector<Locus> loci;
    uint32_t nloc;
    double scaling_factor = 0;

    void _load(std::ifstream& locfp,char* buffer){
        int k = locfp.gcount();
        uint32_t cur=0,chr=0,strand=0,start=0,move=0,locus=0;
        enum Opt {CHR   = 0,
                START = 1,
                MOVE  = 2,
                LOCUS = 3};
//        uint32_t elem = Opt::CHR;
//        Position pos;
//        for(int i=0;i<k;i++){
//            switch(buffer[i]){
//                case '\n':
//                    //end of line
//                    pos.add_move(move);
//                    this->pce = this->index.insert(pos); // write the last entry
//                    pos.clear();
//                    elem = Opt::LOCUS;
//                    chr = 0;
//                    break;
//                case '\t':
//                    // end of a full coordinate - write last integer
//                    pos.add_move(move);
//                    this->pce = this->index.insert(pos); // write the last entry
//                    pos.clear();
//                    elem = Opt::LOCUS;
//                    chr = 0;
//                    break;
//                case ' ':
//                    //end of current int - write move
//                    pos.add_move(move);
//                    elem = Opt::MOVE;
//                    move = 0;
//                    break;
//                case ':':
//                    //end of current int - write start
//                    pos.set_start(start);
//                    elem = Opt::MOVE;
//                    move = 0;
//                    break;
//                case '@':
//                    // got the geneID
//                    pos.set_locus(locus);
//                    elem = Opt::CHR;
//                    locus = 0;
//                case '-': case '+':
//                    // negative strand - write chromosome
//                    pos.set_chr(chr);
//                    pos.set_strand((uint32_t)buffer[i]);
//                    elem = Opt::START;
//                    start = 0;
//                    break;
//                case '0': case '1': case '2': case '3':
//                case '4': case '5': case '6': case '7':
//                case '8': case '9':
//                    //add to the current integer
//                    switch(elem){
//                        case 0: //chromosome
//                            chr = 10*chr + buffer[i] - '0';
//                            break;
//                        case 1:
//                            start = 10*start + buffer[i] - '0';
//                            break;
//                        case 2:
//                            move = 10*move + buffer[i] - '0';
//                            break;
//                        case 3:
//                            locus = 10*locus + buffer[i] - '0';
//                            break;
//                        default:
//                            std::cerr<<"should never happen"<<std::endl;
//                            exit(1);
//                    }
//                    break;
//                default:
//                    std::cerr<<"unrecognized character"<<std::endl;
//                    std::cerr<<buffer[i]<<std::endl;
//                    exit(1);
//            }
//        }
//        if(pos.size() >= 4){ // no end of line character, so need to write the last entry
//            this->pce = this->index.insert(pos);
//        }
    }
};

class Position{
public:
    Position() = default;
    Position(uint32_t chr,uint32_t strand,uint32_t start,uint32_t locus){
        this->chr = chr;
        this->strand = strand;
        this->start = start;
        this->locus = locus;
    }
    ~Position() = default;

    void add_move(uint32_t move){this->moves.push_back(move);this->num_elems++;}
    void set_chr(uint32_t chr){this->chr=chr;this->num_elems++;}
    void set_strand(uint32_t strand){this->strand=strand;this->num_elems++;}
    void set_start(uint32_t start){this->start=start;this->num_elems++;}
    void set_locus(uint32_t locus){this->locus=locus;this->num_elems++;}

    std::string get_strg() const {
        std::string res;
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

    // print using Positions from the kmer_coords
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
        std::cerr<<this->index.size()<<std::endl;
        for(auto &mm : this->index){
            std::cerr<<mm.get_strg()<<std::endl;
        }
    }

    void set_kmerlen(int kmerlen){this->kmerlen = kmerlen;}

    void save_multimappers(std::string& outFP){
        std::ofstream multi_ss(outFP.c_str());
        for(auto &kv : this->kmer_coords){
            if(kv.second.size()>1) {
                for (auto &cv : kv.second) {
                    multi_ss << cv.get_strg() << "\t";
                }
                multi_ss << std::endl;
            }
        }
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
    void add_sequence(std::string& seq,GffObj &p_trans,uint32_t geneID){
        GList<GffExon>& exon_list = p_trans.exons;

        int el_pos = 0; // current position within the exon list
        int exon_pos = 0;
        std::string kmer = ""; // currently evaluated kmer
        GffExon *cur_exon = exon_list[el_pos];
        for(int i=0;i<seq.size()-this->kmerlen+1;i++){ // iterate over all kmers in the sequence
            kmer=seq.substr(i,this->kmerlen);

            Position p(p_trans.gseq_id,(uint32_t)p_trans.strand,cur_exon->start+exon_pos,geneID); // initialize position and add information to it accordingly

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
    // this function appends the intron/exon information to moves of of the position object
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

    // TODO: needs to handle reverse complement of a kmer
    // TODO: for efficient lookup pcord needs to contain pointers to related multimappers

    void _load(std::ifstream &infp,char *buffer){
        int k = infp.gcount();
        uint32_t cur=0,chr=0,strand=0,start=0,move=0,locus=0;
        enum Opt {CHR   = 0,
            START = 1,
            MOVE  = 2,
            LOCUS = 3};
        uint32_t elem = Opt::CHR;
        Position pos;
        for(int i=0;i<k;i++){
            switch(buffer[i]){
                case '\n':
                    //end of line
                    pos.add_move(move);
                    this->pce = this->index.insert(pos); // write the last entry
                    pos.clear();
                    elem = Opt::LOCUS;
                    chr = 0;
                    break;
                case '\t':
                    // end of a full coordinate - write last integer
                    pos.add_move(move);
                    this->pce = this->index.insert(pos); // write the last entry
                    pos.clear();
                    elem = Opt::LOCUS;
                    chr = 0;
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
                    move = 0;
                    break;
                case '@':
                    // got the geneID
                    pos.set_locus(locus);
                    elem = Opt::CHR;
                    locus = 0;
                case '-': case '+':
                    // negative strand - write chromosome
                    pos.set_chr(chr);
                    pos.set_strand((uint32_t)buffer[i]);
                    elem = Opt::START;
                    start = 0;
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
        if(pos.size() >= 4){ // no end of line character, so need to write the last entry
            this->pce = this->index.insert(pos);
        }
    }

    int kmerlen;

    typedef std::unordered_set<Position> pcord;
    pcord index;
    std::pair<pcord::iterator,bool> pce;
    std::unordered_map<std::string,pcord> kmer_coords;
    std::pair<std::unordered_map<std::string,pcord>::iterator,bool> kce;

    // now time to write the unique kmers for transcripts
    std::unordered_map<std::string,int> uniq_cnt; // counts of unique kmers per transcript
    std::pair<std::unordered_map<std::string,int>::iterator,bool> ex_ucnt; // exists or not

};


#endif //TRANS2GENOME_MULTIMAP_H
