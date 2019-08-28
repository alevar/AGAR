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
#include "htslib/sam.h"

#include "Multimap.hh"
#include "arg_parse.h"

#define CMATCH  0
#define CINS    1
#define CDEL    2
#define CREF_SKIP   3
#define CSOFT_CLIP  4
#define CHARD_CLIP  5
#define CPAD    6
#define CEQUAL  7
#define CDIFF   8

enum Opt {INPUT     = 'i',
          OUTPUT    = 'o',
          MULTIFP   = 'm',
          TGMAPFP   = 't',
          MULTILOC  = 'l',
          GLST      = 'g',
          TLST      = 's',
          INFO      = 'f'};

int get_tid(bam1_t *al){
    uint8_t* ptr_op_1=bam_aux_get(al,"OP");
    if(ptr_op_1){
        return bam_aux2i(ptr_op_1);
    }
    return -1;
}

char get_xs(bam1_t *al){
    uint8_t* ptr_xs_1=bam_aux_get(al,"XS");
    if(ptr_xs_1){
        return bam_aux2A(ptr_xs_1);
    }
    return -1;
}

int get_nh(bam1_t *al){
    uint8_t* ptr_nh_1=bam_aux_get(al,"NH");
    if(ptr_nh_1){
        return bam_aux2i(ptr_nh_1);
    }
    return -1;
}

bool is_multi(bam1_t *al){
    uint8_t* ptr_zz_1=bam_aux_get(al,"ZZ");
    if(ptr_zz_1){
        if(get_nh(al)>1){
            return true;
        }
        return false;
    }
    return false;
}

void set_xs(bam1_t *curAl,uint8_t xs){
    uint8_t* ptr=bam_aux_get(curAl,"XS");
    if(ptr){
        bam_aux_del(curAl,ptr);
    }
    if(xs == '+'){
        bam_aux_append(curAl,"XS",'A',1,(const unsigned char*)"+");
    }
    else if(xs == '-'){
        bam_aux_append(curAl,"XS",'A',1,(const unsigned char*)"-");
    }
    else{
        return;
    }
}

// loads data from the tgmap index into a vector, where index is the transcriptID and the value is the locusID
void load_tgmap(std::string tgmapfp,std::vector<int32_t>& tgmap){
    srand(time(0));
    std::cerr<<"@LOG::loading tgmap data from: "<<tgmapfp<<std::endl;
    std::ifstream infp(tgmapfp,std::ifstream::binary);
    if(infp){
        infp.seekg(0,infp.end);
        int size = infp.tellg();
        infp.seekg(0,infp.beg);
        char *buffer = new char[size];
        if(infp.read(buffer,size)){
            int k = infp.gcount();
            uint32_t transID=0,locID=0;
            enum Opt {TRANS  = 0,
                      LOC    = 1};
            uint32_t elem = Opt::TRANS;
            for(int i=0;i<k;i++){
                switch(buffer[i]){
                    case '\n':
                        //end of line
                        tgmap[transID] = locID;
                        transID = 0; locID = 0;
                        elem = Opt::TRANS;
                        break;
                    case '\t':
                        // end of a transcript
                        if(tgmap.size()<transID){ // resize
                            tgmap.resize(transID+100000);
                        }
                        elem = Opt::LOC;
                        break;
                    case '0': case '1': case '2': case '3':
                    case '4': case '5': case '6': case '7':
                    case '8': case '9':
                        //add to the current integer
                        switch(elem){
                            case 0: //transcript
                                transID = 10*transID + buffer[i] - '0';
                                break;
                            case 1: // locus
                                locID = 10*locID + buffer[i] - '0';
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
            if(elem == 1){ // no end of line character, so need to write the last entry
                tgmap[transID] = locID;
            }
        }
        delete[] buffer;
        std::cerr<<"@LOG::finished loading the tgmap data"<<std::endl;
    }
    else{
        std::cerr<<"@LOG::failed to open the tgmap file"<<std::endl;
        exit(-1);
    }
    infp.close();
}

struct INFO{
    int numTranscripts = 0;
    int maxLocID = 0;
    int kmerlen = 0;
};

void load_info(std::string infofp,struct INFO& info){
    // read file to get important stats
    struct stat buffer{};
    if(stat (infofp.c_str(), &buffer) != 0){ // if file does not exists
        std::cerr<<"@ERROR::Info file is not found: "<<infofp<<std::endl;
        exit(1);
    }
    // read file and save important info
    std::cerr<<"@LOG::Reading the info file: "<<infofp<<std::endl;
    std::string iline;
    std::ifstream infostream(infofp);

    std::getline(infostream,iline);
    info.numTranscripts = std::stoi(iline);

    std::getline(infostream,iline);
    info.maxLocID = std::stoi(iline);

    std::getline(infostream,iline);
    info.kmerlen = std::stoi(iline);

    infostream.close();
    std::cerr<<"@LOG::Loaded info data"<<std::endl;
}

// this function parses the .tlst file and inserts entries into the index
void _load_transcriptome(std::ifstream& tlstfp,char* buffer,std::vector<GffTranscript>& transcriptome,Loci& loci){
    int k = tlstfp.gcount();
    uint32_t tid=0,start=0,end=0,chr=0,locus=0,linenum=0;
    enum Opt {TID = 0,
        LOCUS = 1,
        CHR   = 2,
        START = 3,
        END   = 4};
    uint32_t elem = Opt::TID;
    GSeg exon;
    GffTranscript transcript;
    for(int i=0;i<k;i++){
        switch(buffer[i]){
            case '\n':
                // write last exon and copy the transcript over to the index
                exon.end = end;
                transcript.add_exon(exon);
                transcriptome[linenum] = transcript;
                loci[transcript.get_geneID()].set_chr(transcript.refID);
                linenum++;
                transcript.clear();
                elem = Opt::TID;
                end = 0;
                break;
            case '\t':
                // write tid
                transcript.set_gffID(tid);
                elem = Opt::LOCUS;
                tid = 0;
                break;
            case '_':
                exon.start = start;
                elem = Opt::END;
                start = 0;
                break;
            case ',':
                exon.end = end;
                transcript.add_exon(exon);
                elem = Opt::START;
                end = 0;
                break;
            case '@':
                transcript.set_geneID(locus);
                elem = Opt::CHR;
                locus = 0;
                break;
            case '-': case '+':
                transcript.set_refID(chr);
                transcript.set_strand((uint32_t)buffer[i]);
                elem = Opt::START;
                chr = 0;
                break;
            case '0': case '1': case '2': case '3':
            case '4': case '5': case '6': case '7':
            case '8': case '9':
                //add to the current integer
                switch(elem){
                    case 0: //chromosome
                        tid = 10*tid + buffer[i] - '0';
                        break;
                    case 1:
                        locus = 10*locus + buffer[i] - '0';
                        break;
                    case 2:
                        chr = 10*chr + buffer[i] - '0';
                        break;
                    case 3:
                        start = 10*start + buffer[i] - '0';
                        break;
                    case 4:
                        end = 10*end + buffer[i] - '0';
                        break;
                    default:
                        std::cerr<<"@ERROR::should never happen _load_transcriptome"<<std::endl;
                        exit(1);
                }
                break;
            default:
                std::cerr<<"@ERROR::unrecognized character: "<<buffer[i]<<std::endl;
                exit(1);
        }
    }
    if(!transcript.empty){
        transcript.add_exon(exon);
        transcriptome[linenum] = transcript;
    }
}

void load_tlst(std::string transfp,std::vector<GffTranscript>& transcriptome,Loci& loci){
    struct stat buffer{};
    if(stat (transfp.c_str(), &buffer) != 0){ // if file does not exists
        std::cerr<<"@ERROR::Locus file: "<<transfp<<" is not found. Check that correct index is provided."<<std::endl;
        exit(1);
    }
    std::cerr<<"@LOG::loading locus data from: "<<transfp<<std::endl;
    std::ifstream tlstfp(transfp,std::ifstream::binary);
    if(tlstfp){
        tlstfp.seekg(0,tlstfp.end);
        int size = tlstfp.tellg();
        tlstfp.seekg(0,tlstfp.beg);
        char *buffer = new char[size];
        if(tlstfp.read(buffer,size)){
            _load_transcriptome(tlstfp,buffer,transcriptome,loci);
        }
        delete[] buffer;
        std::cerr<<"@LOG::finished loading the locus data"<<std::endl;
    }
    else{
        std::cerr<<"@ERROR::failed to open the locus file"<<std::endl;
    }
    tlstfp.close();
}

// parse cigar and convert it to moves for the Position
int add_moves(bam1_t *al,Position& pos_obj){
    int soft_start = 0;
    uint16_t move = 0;
    for (uint8_t c=0;c<al->core.n_cigar;++c){
        uint32_t *cigar_full=bam_get_cigar(al);
        int opcode=bam_cigar_op(cigar_full[c]);
        int length=bam_cigar_oplen(cigar_full[c]);
        if (opcode==BAM_CINS){}
        else if (opcode==BAM_CDEL){
            move+=length;
        }
        else if (opcode==BAM_CSOFT_CLIP){
            if(c==0) {
                soft_start+=length;
            }
        }
        else if (opcode==BAM_CMATCH){
            move+=length;
        }
        else if (opcode==BAM_CREF_SKIP) {
            pos_obj.add_move(move);
            pos_obj.add_move(length+1);
            move = 0;
        }
        else{
            std::cerr<<"unknown opcode: "<<opcode<<" in "<<bam_get_qname(al)<<std::endl;
            exit(-1);
        }
    }
    if(move>0){
        pos_obj.add_move(move);
    }
    return soft_start;
}

// this function loads the multiloc file
void load_multiloc(std::string multilocfp,std::vector<std::vector<int32_t>>& multilocs){
    srand(time(0));
    std::cerr<<"@LOG::loading multiloc data from: "<<multilocfp<<std::endl;
    std::ifstream infp(multilocfp,std::ifstream::binary);
    bool ce = false;
    if(infp){
        infp.seekg(0,infp.end);
        int size = infp.tellg();
        infp.seekg(0,infp.beg);
        char *buffer = new char[size];
        if(infp.read(buffer,size)){
            int k = infp.gcount();
            uint32_t idx=0;
            int32_t multiloc=0;
            enum Opt {IDX  = 0, // index in the vector
                      LOC  = 1};
            uint32_t elem = Opt::IDX;
            for(int i=0;i<k;i++){
                switch(buffer[i]){
                    case '\n':
                        //end of line
                        multilocs[idx].push_back(multiloc);
                        idx = 0;
                        multiloc=0;
                        elem = Opt::IDX;
                        break;
                    case '\t':
                        // end of a transcript
                        if(multilocs.size()<idx){ // resize
                            multilocs.resize(idx+100000);
                        }
                        elem = Opt::LOC;
                        break;
                    case ',':
                        // new multi loc
                        multilocs[idx].push_back(multiloc);
                        multiloc = 0;
                        break;
                    case '0': case '1': case '2': case '3':
                    case '4': case '5': case '6': case '7':
                    case '8': case '9':
                        //add to the current integer
                        switch(elem){
                            case 0: // index
                                idx = 10*idx + buffer[i] - '0';
                                break;
                            case 1: // locus
                                multiloc = 10*multiloc + buffer[i] - '0';
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
            if(elem == 1){ // no end of line character, so need to write the last entry
                multilocs[idx].push_back(multiloc);
            }
        }
        delete[] buffer;
        std::cerr<<"@LOG::finished loading the multiloc data"<<std::endl;
    }
    else{
        std::cerr<<"@LOG::failed to open the multiloc file"<<std::endl;
        exit(-1);
    }
    infp.close();
}

void change_nh_flag(bam1_t *curAl,int nh){
    uint8_t* ptr_nh_1=bam_aux_get(curAl,"NH");
    if(ptr_nh_1){
        bam_aux_del(curAl,ptr_nh_1);
    }
    bam_aux_append(curAl,"NH",'i',4,(uint8_t*)&nh);
}

int main(int argc, char** argv) {

    ArgParse args("spoon help");
    args.add_string(Opt::INPUT,"input","","input alignment in BAM/CRAM format",true);
    args.add_string(Opt::OUTPUT,"output","","basename of the output files",true);
    args.add_string(Opt::MULTIFP,"multi","","path to the t2g/agar multimappers index",true);
    args.add_string(Opt::TGMAPFP,"tgmap","","path to the t2g/agar transcript to locus mapping",true);
    args.add_string(Opt::MULTILOC,"multiloc","","path to the file containing multimapping loci,",true);
    args.add_string(Opt::GLST,"glst","","path to the t2g/agar glst file,",true);
    args.add_string(Opt::TLST,"tlst","","path to the t2g/agar tlst file,",true);
    args.add_string(Opt::INFO,"info","","path to the t2g/agar info file,",true);

    if(argc == 1){
        std::cerr<<args.get_help()<<std::endl;
        exit(-1);
    }

    if(strcmp(argv[1],"--help")==0){
        std::cerr<<args.get_help()<<std::endl;
        exit(-1);
    }

    args.parse_args(argc,argv);

    std::cout<<"working with: "<<args.get_string(Opt::INPUT)<<std::endl;
    std::cout<<"writing output to: "<<args.get_string(Opt::OUTPUT)<<std::endl;

    std::vector<int32_t> tgmap;
    tgmap.reserve(100000);
    load_tgmap(args.get_string(Opt::TGMAPFP),tgmap);

    std::vector<std::vector<int32_t>> multilocs;
    multilocs.reserve(100000);
    multilocs.resize(100000);
    load_multiloc(args.get_string(Opt::MULTILOC),multilocs);

    Loci loci;
    std::string glst_fname = args.get_string(Opt::GLST);
    loci.load(glst_fname);

    struct INFO info;
    load_info(args.get_string(Opt::INFO),info);

    std::vector<GffTranscript> transcriptome;
    transcriptome.reserve(info.numTranscripts+1);
    load_tlst(args.get_string(Opt::TLST),transcriptome,loci);

    Multimap mmap;
    mmap.load(args.get_string(Opt::MULTIFP));

    std::string out_base(args.get_string(Opt::OUTPUT));

    samFile *al = hts_open(args.get_string(Opt::INPUT).c_str(),"r");
    bam_hdr_t *al_hdr = sam_hdr_read(al);

    samFile *outSAM=sam_open(args.get_string(Opt::OUTPUT).c_str(),"wb");
    bam_hdr_t *outSAM_header=bam_hdr_dup(al_hdr);
    int ret_hdr = sam_hdr_write(outSAM,outSAM_header);

    bam1_t *curAl = bam_init1(); // initialize the alignment record

    GffTranscript p_trans;

    bool skip_next = false;
    bam1_t * prevAl = bam_init1();
    bool prev_init = false;

    Position pos;

    while(sam_read1(al,al_hdr,curAl)>=0) {
        // first do correction of the wrong XS tag - this somehow requires us to get access to the transcript
        // we can do this potentially by looking up multimappers and once again finding the multimapper and seing if it has the correct xs tag and if not switching it up - only needed for the multimappers
        // this can be traced back through the op tag which specifies the original transcript to which the alignment was mapped
        if(is_multi(curAl)){
            int tid_orig = get_tid(curAl);
            if(tid_orig == -1){
                std::cerr<<"no OP aux tag found for the multimapping record (ZZ set)"<<std::endl;
                exit(-1);
            }

            p_trans = transcriptome[tid_orig];

            // get locus
            int locid = tgmap[p_trans.gffID];
            // iterate over multimappers:
            bool pos_strand=false,neg_strand=false;
            std::vector<uint32_t> mls;
            for(auto& ml : multilocs[locid]){
                // check if locus covers the current alignment
                if(curAl->core.tid == loci[ml].get_chr() &&
                   curAl->core.pos+1 >= loci[ml].get_start() && curAl->core.pos+1 <= loci[ml].get_end()){ // found correct locus
                    mls.push_back(ml);
                }
            }
            if(curAl->core.tid == loci[locid].get_chr() &&
               curAl->core.pos+1 >= loci[locid].get_start() && curAl->core.pos+1 <= loci[locid].get_end()){ // found correct locus
                mls.push_back(locid);
            }

            bool found_xs = false;
            // now having all valid multimappers we can proceed
            if(mls.size()==1){ // just set the new XS tag
                set_xs(curAl,loci[mls[0]].get_strand());
            }
            else if(mls.size()>1){ // need to consider the actual multimapping positions and figure out what to do with the xs...
                bool is_duplicate = false;
                std::set<int> multi_xs;
                std::set<int> multi_nh;
                for(auto& ml : mls){
                    pos.clear();
                    pos.set_chr(loci[ml].get_chr());
                    pos.set_strand(loci[ml].get_strand());
                    pos.set_locus(ml);
                    int soft_start = add_moves(curAl,pos);
                    pos.set_start(curAl->core.pos+1);
                    int new_nh = mmap.get_nh(pos);
                    if(new_nh == -1){
                        continue;
                    }
                    change_nh_flag(curAl,new_nh);
                    int multi_xs_res = mmap.get_multi_xs(pos);
                    if(multi_xs_res < 0){ // not valid
                        std::cerr<<"ERROR:: strand invalid: "<<bam_get_qname(curAl)<<std::endl;
                        exit(-1);
                    }
                    // what we can do here is the following:
                    // TODO: we can incorporate information about the mates start position
                    // TODO: and move skip_next block from the start of the loop to here
                    else{ // TODO: needs to be moved outside the for-loop
                        found_xs = true;
                        if(multi_xs_res == 0){
                            is_duplicate = true;
                            multi_nh.insert(new_nh);
                            break; // found a duplicate - can safely exit and decide what to do with the read
                        }
                        else{
                            multi_xs.insert(multi_xs_res);
                            multi_nh.insert(new_nh);
                        }
                    }
                }
                if(is_duplicate){
                    if(prev_init && // first alignment seed
                       std::strcmp(bam_get_qname(curAl),bam_get_qname(prevAl))==0 && // names match
                       curAl->core.pos == prevAl->core.pos && // postions also match
                       (curAl->core.flag & 0x40) == (prevAl->core.flag & 0x40)){ // and they belong to the same mate
                        if(!(curAl->core.flag & 0x100)){
                            std::cerr<<"ERROR::attempting to remove primary alignment: "<<bam_get_qname(curAl)<<"\t"<<curAl->core.pos<<std::endl;
                            exit(-1);
                        }
                        prevAl = bam_copy1(prevAl,curAl);
                        continue; // skip this read and continue to the next iteration of the while loop, thus skipping writing the current read to disk
                    }
                    else{ // check multi_xs and multi_nh and adjust accordingly and write to disk
                        prev_init = true;
                        prevAl = bam_copy1(prevAl,curAl);
                        set_xs(curAl,'0');
                        change_nh_flag(curAl,*multi_nh.begin());
                    }
                }
                else{
                    if(multi_xs.size()==1){
                        set_xs(curAl,*multi_xs.begin());
                        if(multi_nh.size()==1){
                            change_nh_flag(curAl,*multi_nh.begin());
                        }
                        else{
                            std::cerr<<"different NH tags for the read: "<<bam_get_qname(curAl)<<std::endl;
                            exit(-1);
                        }
                    }
                    else if(multi_xs.size()==0) { // did not find in the multimapper database - hence need to make corrections based on whatever else method

                    }
                    else{
                        std::cerr<<"different XS tags for the read: "<<bam_get_qname(curAl)<<std::endl;
                    }
                }
                //
                if(!found_xs){
                    std::cerr<<"ERROR::did not find XS in read: "<<bam_get_qname(curAl)<<" at position: "<<curAl->core.pos<<std::endl;
                    continue;
                }
            }
            else{
//                std::cerr<<"ERROR::NO MUlTIMAPPER FOUND: "<<bam_get_qname(curAl)<<std::endl;
            }
            int ret_val = sam_write1(outSAM, outSAM_header, curAl);
//            std::cout<<"XS after correction: "<<get_xs(curAl)<<std::endl;
//            std::cout<<"NH after correction: "<<get_nh(curAl)<<std::endl;
        }
        else{
            // not a multimapper so can be written as is
            int ret_val = sam_write1(outSAM, outSAM_header, curAl);
        }
    }
    bam_destroy1(curAl);
    bam_destroy1(prevAl);

    bam_hdr_destroy(al_hdr);
    sam_close(al);
    sam_close(outSAM);
    bam_hdr_destroy(outSAM_header);

    return 0;
}
