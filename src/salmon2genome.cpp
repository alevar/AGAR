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
#include "map2gff_salmon.h"

// salmon2genome -g ~/JHU/transcriptome/hisatTrans2Genome/errorData/me_mi/data/ann.gff -i ~/JHU/transcriptome/hisatTrans2Genome/errorData/me_mi/al_1/sample.gffread.bam -s ~/JHU/transcriptome/hisatTrans2Genome/errorData/me_mi/data/genomic_header.sam -o ./test.bam

enum Opt {IN_AL     = 'i',
        OUT_AL    = 'o',
        THREADS   = 'p',
        INDEX     = 'x',
        MULTI     = 'm',
        ABUNDANCE = 'a',
        UNALIGNED = 'u'};

int main(int argc, char** argv) {

    ArgParse args("salmon2genome");
    args.add_string(Opt::IN_AL,"input","","input alignment SAM/BAM",true);
    args.add_string(Opt::OUT_AL,"output","","output file path (BAM)",true);
    args.add_int(Opt::THREADS,"threads",1,"number of threads (default = 1)",false);
    args.add_string(Opt::INDEX,"index","","path and basename of the index built by gtf_to_fasta",true);
    args.add_flag(Opt::MULTI,"multi","whether to search and evaluate multimappers",false);
    args.add_string(Opt::ABUNDANCE,"abund","","use abundances precomputed",false);
    args.add_flag(Opt::UNALIGNED,"unal","search for unaligned reads, extract from alignment into separate files",false);

    if(strcmp(argv[1],"--help")==0){
        std::cerr<<args.get_help()<<std::endl;
        exit(1);
    }

    args.parse_args(argc,argv);

    Map2GFF_SALMON gffMapper(args.get_string(Opt::IN_AL),args.get_string(Opt::OUT_AL),args.get_string(Opt::INDEX),args.get_int(Opt::THREADS),args.get_flag(Opt::MULTI));
    if(args.is_set(Opt::ABUNDANCE)){
        gffMapper.load_abundances(args.get_string(Opt::ABUNDANCE));
    }
    std::cerr<<"Begin Translating Coordinates"<<std::endl;
    gffMapper.convert_coords(args.get_flag(Opt::UNALIGNED)); // TODO: need to make sure that multimappers are evaluated if the index is provided
    return 0;
}
