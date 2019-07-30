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

// TODO: log of events to see if anything goes wrong during gtex realignment

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
        NUM_MULTI = 'k',
        MISALIGN  = 's',
        PERCENT   = 'e',
        OUTLIER_STDV   = 't',
        OUTLIER_RAW   = 'r',
        NUM_READS_PRECOMP = 'n'};

int main(int argc, char** argv) {

    ArgParse args("salmon2genome");
    args.add_string(Opt::IN_AL,"input","","input alignment SAM/BAM",false);
    args.add_string(Opt::OUT_AL,"output","","output file path (BAM)",true);
    args.add_int(Opt::THREADS,"threads",1,"number of threads (default = 1)",false);
    args.add_string(Opt::INDEX,"index","","path and basename of the index built by gtf_to_fasta",true);
    args.add_flag(Opt::MULTI,"multi","whether to search and evaluate multimappers",false);
    args.add_string(Opt::ABUNDANCE,"abund","","use abundances precomputed. Only gene-level abundance estimates computed by salmon are supported at the moment",false);
    args.add_flag(Opt::UNALIGNED,"unal","search for unaligned reads, extract from alignment into separate files",false);
    args.add_flag(Opt::UNIQ,"uniq","input alignment contains only 1 mapping per read (no secondary alignments present such as in bowtie k1 mode)",false);
    args.add_int(Opt::FRAGLEN,"fraglen",200000,"fragment length of the paired reads",false);
    args.add_flag(Opt::ALL_MULTI,"all","whether to output all multimappers and not assign them based on likelihood. This flag negates -k",false);
    args.add_int(Opt::NUM_MULTI,"nmult",1,"The number of most likely multimappers to report",false);
    args.add_flag(Opt::MISALIGN,"mis","try to eliminate misaligned reads based on the error distribution",false);
    args.add_int(Opt::PERCENT,"perc",10,"percent of the reads to be evaluated when pre-loading data before parsing",false);
    args.add_int(Opt::NUM_READS_PRECOMP,"nrp",500000,"number of reads to precompute based on the streaming",false);
    args.add_int(Opt::OUTLIER_STDV,"outlier_stdv",5,"Poisson threshold as the number of standard deviations. Everythin above the threshold will be discarded as a misalignment",false);
    args.add_int(Opt::OUTLIER_RAW,"outlier_raw",10,"Edit distance threshold for misalignments. For paired alignments the threshold is applied to each read separately",false);

    // TODO: implement the -k mode in which only a certain number of most frequent multimappers is reported (2,3,etc depending on the value set)

    if(strcmp(argv[1],"--help")==0){
        std::cerr<<args.get_help()<<std::endl;
        exit(1);
    }

    args.parse_args(argc,argv);

    std::string inputAlFP = args.is_set(Opt::IN_AL)?args.get_string(Opt::IN_AL):"-";
    std::cout<<"input alignment: "<<inputAlFP<<std::endl;
    Converter converter(inputAlFP,args.get_string(Opt::OUT_AL),args.get_string(Opt::INDEX),args.get_int(Opt::THREADS),args.get_flag(Opt::MULTI));
    if(args.is_set(Opt::ABUNDANCE)){
        std::cerr<<"will be using abundance from "<<args.get_string(Opt::ABUNDANCE)<<std::endl;
        converter.load_abundances(args.get_string(Opt::ABUNDANCE));
    }
    if(args.get_flag(Opt::UNALIGNED)){
        std::cerr<<"will report unalgned reads"<<std::endl;
        converter.set_unaligned();
    }
    if(args.get_flag(Opt::UNIQ)){
        std::cerr<<"trusting the user that input is uniq"<<std::endl;
        converter.set_k1();
    }
    if(args.is_set(Opt::FRAGLEN)){
        std::cerr<<"fragment length is set to: "<<args.get_int(Opt::FRAGLEN)<<std::endl;
        converter.set_fraglen(args.get_int(Opt::FRAGLEN));
    }
    if(args.is_set(Opt::NUM_MULTI)){
        std::cerr<<"num multi enabled and set to "<<args.get_int(Opt::NUM_MULTI)<<std::endl;
        converter.set_num_multi(args.get_int(Opt::NUM_MULTI));
    }
    if(args.is_set(Opt::ALL_MULTI)){ // TODO: does this need to be converted into get_flag instead of is_set?
        std::cerr<<"all multi reporting enabled"<<std::endl;
        converter.set_all_multi();
    }
    if(args.get_flag(Opt::MISALIGN)){
        std::cerr<<"misalignment detection enabled"<<std::endl;
        converter.set_misalign();

        if(args.is_set(Opt::OUTLIER_STDV)){
            std::cerr<<"misalignment is set to be detected at "<<args.get_int(Opt::OUTLIER_STDV)<<std::endl;
            converter.set_stdv(args.get_int(Opt::OUTLIER_STDV));
        }
        if(args.is_set(Opt::OUTLIER_RAW)){
            std::cerr<<"misalignment is set to be detected at "<<args.get_int(Opt::OUTLIER_RAW)<<std::endl;
            converter.set_raw(args.get_int(Opt::OUTLIER_RAW));
        }
    }

    if(inputAlFP=="-"){
        std::cerr<<"Begin pre-loading data from stream"<<std::endl;
        converter.precompute_save(args.get_int(Opt::NUM_READS_PRECOMP));
        std::cerr<<"Done pre-loading data from stream"<<std::endl;
    }
    else{
        std::cerr<<"Begin pre-loading data"<<std::endl;
        converter.precompute(args.get_int(Opt::PERCENT));
        std::cerr<<"Done pre-loading data"<<std::endl;
    }

    std::cerr<<"Begin Translating Coordinates"<<std::endl;
    converter.convert_coords();
    std::cerr<<"Done Translating Coordinates"<<std::endl;
    if(inputAlFP=="-"){
        std::cerr<<"Begin Translating Coordinates of precomputed reads"<<std::endl;
        converter.convert_coords_precomp();
        std::cerr<<"Done Translating Coordinates of precomputed reads"<<std::endl;
    }

    converter.print_stats();

    return 0;
}
