#include <iostream>
#include <math.h>
#include <sys/time.h>
#include <unistd.h>
#include <htslib/sam.h>
#include <string.h>
#include <sstream>

#include "arg_parse.h"

#include "GVec.hh"
#include "tokenize.h"
#include "map2gff.h"

// trans2genome -g ~/JHU/transcriptome/hisatTrans2Genome/errorData/me_mi/data/ann.gff -i ~/JHU/transcriptome/hisatTrans2Genome/errorData/me_mi/al_1/sample.gffread.bam -s ~/JHU/transcriptome/hisatTrans2Genome/errorData/me_mi/data/genomic_header.sam -o ./test.bam

enum Opt {TLST_FP   = 't',
          IN_AL     = 'i',
          OUT_AL    = 'o',
          GEN_HDR   = 's',
          MULTI     = 'm',
          GLST_FP   = 'g',
          MULTI_FLAG= 'f',
          THREADS   = 't'};

int main(int argc, char** argv) {

    ArgParse args("trans2genome");
    args.add_string(Opt::TLST_FP,"tlst","","",true);
    args.add_string(Opt::IN_AL,"input","","",true);
    args.add_string(Opt::OUT_AL,"output","","",true);
    args.add_string(Opt::GEN_HDR,"header","","",true);
    args.add_string(Opt::MULTI,"multi","","",false);
    args.add_string(Opt::GLST_FP,"glst","","",true);
    args.add_flag(Opt::MULTI_FLAG,"mf","",false);
    args.add_int(Opt::THREADS,"threads",1,"",false);
    
    args.parse_args(argc,argv);
    Map2GFF gffMapper(args.get_string(Opt::TLST_FP),args.get_string(Opt::IN_AL),args.get_string(Opt::MULTI),args.get_string(Opt::GLST_FP),args.get_flag(Opt::MULTI_FLAG),args.get_int(Opt::THREADS));
    gffMapper.convert_coords(args.get_string(Opt::OUT_AL),args.get_string(Opt::GEN_HDR));
    return 0;
}
