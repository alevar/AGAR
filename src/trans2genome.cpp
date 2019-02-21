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

enum Opt {GFF_FP    = 'g',
          IN_AL     = 'i',
          OUT_AL    = 'o',
          GEN_HDR   = 's',
          MULTI     = 'm'};

int main(int argc, char** argv) {

    ArgParse args("Map2GFF");
    args.add_string(Opt::GFF_FP,"gff","","");
    args.add_string(Opt::IN_AL,"input","","");
    args.add_string(Opt::OUT_AL,"output","","");
    args.add_string(Opt::GEN_HDR,"header","","");
    args.add_string(Opt::MULTI,"multi","","");
    
    args.parse_args(argc,argv);
    Map2GFF gffMapper(args.get_string(Opt::GFF_FP),args.get_string(Opt::IN_AL),args.get_string(Opt::MULTI));
    gffMapper.convert_coords(args.get_string(Opt::OUT_AL),args.get_string(Opt::GEN_HDR));
    return 0;
}
