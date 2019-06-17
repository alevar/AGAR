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

// this class describes the multimappers in the index transcriptome
// and facilitates efficient storage and lookup
class Multimap{
public:
    Multimap() = default;
    ~Multimap() = default;

    void add_sequence(){

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

    void print_multimapers(){
        std::cerr<<this->index.size()<<std::endl;
        for(auto &mm : this->index){
            for(auto &c : mm){
                std::cerr<<c<<" ";
            }
            std::cerr<<std::endl;
        }
    }
private:

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
};


#endif //TRANS2GENOME_MULTIMAP_H
