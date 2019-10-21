#include <iostream>
#include <math.h>
#include <string.h>

#include "arg_parse.h"
#include "tokenize.h"
#include "Comparator.h"

enum Opt {TRANS_AL   = 't',
          GENOM_AL   = 'g',
          OUTPUT_AL  = 'o',
          THREADS   = 'p',
          INDEX     = 'x',
          MULTI     = 'm',
          ALL_MULTI = 'l',
          NUM_MULTI = 'k',
          ABUNDANCE = 'a',
          FRAGLEN   = 'f',
          PERCENT   = 'e',
          NUM_READS_PRECOMP = 'n'};

int main(int argc, char** argv) {

    ArgParse args("agar_diff");
    args.add_string(Opt::TRANS_AL,"trans","","Transcriptome alignment in SAM/BAM",true);
    args.add_string(Opt::GENOM_AL,"trans","","Genome alignment in SAM/BAM",true);
    args.add_string(Opt::OUTPUT_AL,"output","","output file path (BAM)",true);
    args.add_int(Opt::THREADS,"threads",1,"number of threads (default = 1)",false);
    args.add_string(Opt::INDEX,"index","","path and basename of the index built by gtf_to_fasta",true);
    args.add_flag(Opt::MULTI,"multi","whether to search and evaluate multimappers",false);
    args.add_flag(Opt::ALL_MULTI,"all","whether to output all multimappers and not assign them based on likelihood. This flag negates -k",false);
    args.add_string(Opt::ABUNDANCE,"abund","","use abundances precomputed. Only gene-level abundance estimates computed by salmon are supported at the moment",false);
    args.add_int(Opt::NUM_MULTI,"nmult",1,"The number of most likely multimappers to report",false);
    args.add_int(Opt::PERCENT,"perc",10,"percent of the reads to be evaluated when pre-loading data before parsing",false);
    args.add_int(Opt::NUM_READS_PRECOMP,"nrp",500000,"number of reads to precompute based on the streaming",false);
    args.add_int(Opt::FRAGLEN,"fraglen",200000,"fragment length of the paired reads",false);

    if(argc <= 1 || strcmp(argv[1],"--help")==0){
        std::cerr<<args.get_help()<<std::endl;
        exit(1);
    }

    args.parse_args(argc,argv);

    // first create the execution string
    std::string cl="agar_diff ";
    for (int i=0;i<argc;i++){
        if(i==0){
            cl+=argv[i];
        }
        else{
            cl+=" ";
            cl+=argv[i];
        }
    }

    std::string trans_al_fname = args.get_string(Opt::TRANS_AL);
    std::cerr<<"@LOG::transcriptome alignment: "<<trans_al_fname<<std::endl;
    std::string genome_al_fname = args.get_string(Opt::GENOM_AL);
    std::cerr<<"@LOG::genome alignment: "<<trans_al_fname<<std::endl;
    Comparator comparator(trans_al_fname,genome_al_fname,args.get_string(Opt::OUTPUT_AL),args.get_string(Opt::INDEX),args.get_int(Opt::THREADS),args.get_flag(Opt::MULTI),cl);
    if(args.is_set(Opt::ABUNDANCE)){
        std::cerr<<"@LOG::will be using abundance from "<<args.get_string(Opt::ABUNDANCE)<<std::endl;
        comparator.load_abundances(args.get_string(Opt::ABUNDANCE));
    }
    if(args.is_set(Opt::FRAGLEN)){
        std::cerr<<"@LOG::fragment length is set to: "<<args.get_int(Opt::FRAGLEN)<<std::endl;
        comparator.set_fraglen(args.get_int(Opt::FRAGLEN));
    }
    if(args.is_set(Opt::NUM_MULTI)){
        std::cerr<<"@LOG::num multi enabled and set to "<<args.get_int(Opt::NUM_MULTI)<<std::endl;
        comparator.set_num_multi(args.get_int(Opt::NUM_MULTI));
    }
    if(args.is_set(Opt::ALL_MULTI)){ // TODO: does this need to be converted into get_flag instead of is_set?
        std::cerr<<"@LOG::all multi reporting enabled"<<std::endl;
        comparator.set_all_multi();
    }
// TODO:    if(args.is_set(Opt::MULTI) && !args.is_set(Opt::ALL_MULTI) && !args.is_set(Opt::ABUNDANCE)){ // preload only when estimating likelihood of multimappers
//        std::cerr<<"@LOG::Begin pre-loading data"<<std::endl;
//        comparator.precompute(args.get_int(Opt::PERCENT));
//        std::cerr<<"@LOG::Done pre-loading data"<<std::endl;
//    }
    std::cerr<<"@LOG::Begin Translating Coordinates"<<std::endl;
    comparator.merge();
    std::cerr<<"@LOG::Done Translating Coordinates"<<std::endl;

    std::string out_abund_fname = args.get_string(Opt::OUTPUT_AL);
    out_abund_fname.append(".loc.abunds");
    comparator.print_abundances(out_abund_fname);

    comparator.print_stats();

    return 0;
}
