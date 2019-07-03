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
#include "Converter.h"

// trans2genome -g ~/JHU/transcriptome/hisatTrans2Genome/errorData/me_mi/data/ann.gff -i ~/JHU/transcriptome/hisatTrans2Genome/errorData/me_mi/al_1/sample.gffread.bam -s ~/JHU/transcriptome/hisatTrans2Genome/errorData/me_mi/data/genomic_header.sam -o ./test.bam

enum Opt {IN_AL   = 'i',
        OUT_AL    = 'o',
        THREADS   = 'p',
        INDEX     = 'x',
        MULTI     = 'm',
        ABUNDANCE = 'a',
        UNALIGNED = 'u',
        UNIQ      = 'q',
        FRAGLEN   = 'f',
        ALL_MULTI = 'l',
        NUM_MULTI = 'k'};

int main(int argc, char** argv) {

    ArgParse args("salmon2genome");
    args.add_string(Opt::IN_AL,"input","","input alignment SAM/BAM",true);
    args.add_string(Opt::OUT_AL,"output","","output file path (BAM)",true);
    args.add_int(Opt::THREADS,"threads",1,"number of threads (default = 1)",false);
    args.add_string(Opt::INDEX,"index","","path and basename of the index built by gtf_to_fasta",true);
    args.add_flag(Opt::MULTI,"multi","whether to search and evaluate multimappers",false);
    args.add_string(Opt::ABUNDANCE,"abund","","use abundances precomputed. Only gene-level abundance estimates computed by salmon are supported at the moment",false); // TODO: now need to implement proper loading of external abundance estimates
    args.add_flag(Opt::UNALIGNED,"unal","search for unaligned reads, extract from alignment into separate files",false);
    args.add_flag(Opt::UNIQ,"uniq","input alignment contains only 1 mapping per read (no secondary alignments present such as in bowtie k1 mode)",false);
    args.add_int(Opt::FRAGLEN,"fraglen",200000,"fragment length of the paired reads",false);
    args.add_flag(Opt::ALL_MULTI,"all","whether to output all multimappers and not assign them based on likelihood. This flag negates -k",false); // TODO: needs to be implemented
    args.add_int(Opt::NUM_MULTI,"nmult",1,"The number of most likely multimappers to report",false); // TODO: needs to be implemented

    // TODO: compile the new version of hisat2 which does not miss the frequent multimappers and test wether it performs faster/better than the bowtie2 mode - does anything need ot be changed?

    // TODO: current branch work
    //       1. first implement discordnat pair handling -  namely, currently an alignment can be reported as discordant for two mates that map to different transcripts on the same locus, in which case, we should re-unite them and manage multimappers as a pair
    //       2. secondly implement the -a mode to output all multimappers as before and compare the results

    // TODO: implement the -k mode in which only a certain number of most frequent multimappers is reported (2,3,etc depending on the value set)

    // TODO: do any flags need to be fixed along with the auxilary tags?

    if(strcmp(argv[1],"--help")==0){
        std::cerr<<args.get_help()<<std::endl;
        exit(1);
    }

    args.parse_args(argc,argv);

    Converter converter(args.get_string(Opt::IN_AL),args.get_string(Opt::OUT_AL),args.get_string(Opt::INDEX),args.get_int(Opt::THREADS),args.get_flag(Opt::MULTI));
    if(args.is_set(Opt::ABUNDANCE)){
        converter.load_abundances(args.get_string(Opt::ABUNDANCE));
    }
    if(args.get_flag(Opt::UNALIGNED)){
        converter.set_unaligned();
    }
    if(args.get_flag(Opt::UNIQ)){
        converter.set_k1();
    }
    if(args.is_set(Opt::FRAGLEN)){
        converter.set_fraglen(Opt::FRAGLEN);
    }
    if(args.is_set(Opt::NUM_MULTI)){
        converter.set_num_multi(args.get_int(Opt::NUM_MULTI));
    }
    if(args.is_set(Opt::ALL_MULTI)){
        converter.set_all_multi();
    }
    std::cerr<<"Begin Translating Coordinates"<<std::endl;
    converter.convert_coords(); // TODO: need to make sure that multimappers are evaluated if the index is provided
    return 0;
}
