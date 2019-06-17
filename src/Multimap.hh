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

class Position{
public:
    Position() = default;
    Position(uint32_t chr,uint32_t strand,uint32_t start){
        this->chr = chr;
        this->strand = strand;
        this->start = start;
    }
    ~Position() = default;

    void add_move(uint32_t move){
        this->moves.emplace_back(move);
    }

    std::string get_strg() const {
        std::string res;
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
               this->moves==m.moves;
    }

    uint32_t chr,strand,start;
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
                int k = infp.gcount();
                uint32_t cur=0,chr=0,strand=0,start=0,move=0;
                enum Opt {CHR   = 0,
                    START = 1,
                    MOVE  = 2};
                uint32_t elem = Opt::CHR;
                std::vector<uint32_t> coords;
                for(int i=0;i<k;i++){
                    switch(buffer[i]){
                        case '\n':
                            //end of line
                            coords.emplace_back(move);
                            this->ie = this->index.insert(coords); // write the last entry
                            coords.clear();
                            elem = Opt::CHR;
                            chr = 0;
                            break;
                        case '\t':
                            // end of a full coordinate - write lsat integer
                            coords.emplace_back(move);
                            this->ie = this->index.insert(coords); // write the last entry
                            coords.clear();
                            elem = Opt::CHR;
                            chr = 0;
                            break;
                        case ' ':
                            //end of current int - write move
                            coords.emplace_back(move);
                            elem = Opt::MOVE;
                            move = 0;
                            break;
                        case ':':
                            //end of current int - write start
                            coords.emplace_back(start);
                            elem = Opt::MOVE;
                            move = 0;
                            break;
                        case '-': case '+':
                            // negative strand - write chromosome
                            coords.emplace_back(chr);
                            coords.emplace_back((uint32_t)buffer[i]);
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
                if(coords.size() >= 4){ // no end of line character, so need to write the last entry
                    this->ie = this->index.insert(coords);
                }
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
            for(auto &c : mm){
                std::cerr<<c<<" ";
            }
            std::cerr<<std::endl;
        }
    }

    void set_kmerlen(int kmerlen){this->kmerlen = kmerlen;}

    void save_multimappers(std::string& outFP){
        for(auto &kc : this->kmer_coords){
            std::cout<<kc.first<<std::endl;
        }
    }

    // given a full transcript sequence and the pointer to the transcript description (constituent exons)
    // get all kmers and coordinates and load them into the multimapper index
    void add_sequence(std::string& seq,GffObj &p_trans){
        GList<GffExon>& exon_list = p_trans.exons;

        int el_pos = 0; // current position within the exon list
        int exon_pos = 0;
        std::string kmer = ""; // currently evaluated kmer
        GffExon *cur_exon = exon_list[el_pos];
        for(int i=0;i<seq.size()-this->kmerlen+1;i++){ // iterate over all kmers in the sequence
            kmer=seq.substr(i,this->kmerlen);

            Position p(p_trans.gseq_id,(uint32_t)p_trans.strand,cur_exon->start+exon_pos); // initialize position and add information to it accordingly

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

    typedef std::unordered_set<Position> pcord;
    std::pair<pcord::iterator,bool> pce;
    std::unordered_map<std::string,pcord> kmer_coords; // TODO: here the vector of positions needs to be reimplemented as a set to remove redundancy introduced by shared exons with identical genomic coordinates
    std::pair<std::unordered_map<std::string,pcord>::iterator,bool> kce;

    int kmerlen;

    struct coord_hash { // in case it is implemented as an unordered_map
        uint64_t operator()(const std::vector<uint32_t> &coords ) const{
            uint64_t resHash=1;
            for (auto const& coord: coords) {
                resHash=resHash^std::hash<uint64_t>()(coord);
            }
            return resHash;
        }
    };

    std::unordered_set<std::vector<uint32_t>,coord_hash> index;
    std::pair<std::unordered_set<std::vector<uint32_t>,coord_hash>::iterator,bool> ie; // iterator

    // TODO: deal with uniques here as well if needed
    // now time to write the unique kmers for transcripts
    std::unordered_map<std::string,int> uniq_cnt; // counts of unique kmers per transcript
    std::pair<std::unordered_map<std::string,int>::iterator,bool> ex_ucnt; // exists or not

};


#endif //TRANS2GENOME_MULTIMAP_H
