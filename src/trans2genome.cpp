#include <iostream>
#include <math.h>
#include <string.h>

#include "arg_parse.h"
#include "tokenize.h"
#include "Converter.h"

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
        NUM_READS_PRECOMP = 'n',
        NODISCORD = 'j',
        NOSINGLE = 'g'};

int main(int argc, char** argv) {

    ArgParse args("trans2genome");
    args.add_string(Opt::IN_AL,"input","","input alignment SAM/BAM",false);
    args.add_string(Opt::OUT_AL,"output","","output file path (BAM)",true);
    args.add_int(Opt::THREADS,"threads",1,"number of threads (default = 1)",false);
    args.add_string(Opt::INDEX,"index","","path and basename of the index built by gtf_to_fasta",true);
    args.add_flag(Opt::MULTI,"multi","whether to search and evaluate multimappers",false);
    args.add_string(Opt::ABUNDANCE,"abund","","use abundances precomputed. Only gene-level abundance estimates computed by salmon are supported at the moment",false);
    args.add_flag(Opt::UNIQ,"uniq","input alignment contains only 1 mapping per read (no secondary alignments present such as in bowtie k1 mode)",false);
    args.add_int(Opt::FRAGLEN,"fraglen",200000,"fragment length of the paired reads",false);
    args.add_flag(Opt::ALL_MULTI,"all","whether to output all multimappers and not assign them based on likelihood. This flag negates -k",false);
    args.add_int(Opt::NUM_MULTI,"nmult",1,"The number of most likely multimappers to report",false);
    args.add_flag(Opt::MISALIGN,"mis","try to eliminate misaligned reads based on the error distribution",false);
    args.add_int(Opt::PERCENT,"perc",10,"percent of the reads to be evaluated when pre-loading data before parsing",false);
    args.add_int(Opt::NUM_READS_PRECOMP,"nrp",500000,"number of reads to precompute based on the streaming",false);
    args.add_int(Opt::OUTLIER_STDV,"outlier_stdv",2,"Poisson threshold as the number of standard deviations. Everything above the threshold will be discarded as a misalignment",false);
    args.add_int(Opt::OUTLIER_RAW,"outlier_raw",3,"Edit distance threshold for misalignments. For paired alignments the threshold is applied to each read separately",false);
    args.add_string(Opt::UNALIGNED,"un","","output unaligned reads (singletons and pairs) to separate fastq.",false);
    args.add_flag(Opt::NODISCORD,"nodiscord","report all dicscordant pairs as unaligned. If --un is specified, the reads will be deposited into respective fasta/fastq files. If not, SAM information will be set to indicate unaligned",false);
    args.add_flag(Opt::NOSINGLE,"nosingle","report all mates for which the second mate is unaligned as unaligned. If --un is specified, the reads will be deposited into respective fasta/fastq files. If not, SAM information will be set to indicate unaligned",false);

    // TODO: implement the -k mode in which only a certain number of most frequent multimappers is reported (2,3,etc depending on the value set)

    if(argc <= 1 || strcmp(argv[1],"--help")==0){
        std::cerr<<args.get_help()<<std::endl;
        exit(1);
    }

    args.parse_args(argc,argv);

    // first create the execution string
    std::string cl="trans2genome ";
    for (int i=0;i<argc;i++){
        if(i==0){
            cl+=argv[i];
        }
        else{
            cl+=" ";
            cl+=argv[i];
        }
    }

    std::string inputAlFP = args.is_set(Opt::IN_AL)?args.get_string(Opt::IN_AL):"-";
    std::cerr<<"@LOG::input alignment: "<<inputAlFP<<std::endl;
    Converter converter(inputAlFP,args.get_string(Opt::OUT_AL),args.get_string(Opt::INDEX),args.get_int(Opt::THREADS),args.get_flag(Opt::MULTI),cl);
    if(args.is_set(Opt::ABUNDANCE)){
        std::cerr<<"@LOG::will be using abundance from "<<args.get_string(Opt::ABUNDANCE)<<std::endl;
        converter.load_abundances(args.get_string(Opt::ABUNDANCE));
    }
    if(args.is_set(Opt::UNALIGNED)){
        std::cerr<<"@LOG::will report unaligned reads"<<std::endl;
        converter.set_unaligned(args.get_string(Opt::UNALIGNED));
    }
    if(args.is_set(Opt::NODISCORD)){
        std::cerr<<"@LOG::will report discordant alignments as unaligned unless able to repair"<<std::endl;
        converter.set_discord_unaligned();
    }
    if(args.is_set(Opt::NOSINGLE)){
        std::cerr<<"@LOG::will report mates for which the mate is unaligned as unaligned"<<std::endl;
        converter.set_single_unaligned();
    }
    if(args.get_flag(Opt::UNIQ)){
        std::cerr<<"@LOG::trusting the user that input is uniq"<<std::endl;
        converter.set_k1();
    }
    if(args.is_set(Opt::FRAGLEN)){
        std::cerr<<"@LOG::fragment length is set to: "<<args.get_int(Opt::FRAGLEN)<<std::endl;
        converter.set_fraglen(args.get_int(Opt::FRAGLEN));
    }
    if(args.is_set(Opt::NUM_MULTI)){
        std::cerr<<"@LOG::num multi enabled and set to "<<args.get_int(Opt::NUM_MULTI)<<std::endl;
        converter.set_num_multi(args.get_int(Opt::NUM_MULTI));
    }
    if(args.is_set(Opt::ALL_MULTI)){ // TODO: does this need to be converted into get_flag instead of is_set?
        std::cerr<<"@LOG::all multi reporting enabled"<<std::endl;
        converter.set_all_multi();
    }
    if(args.get_flag(Opt::MISALIGN)){
        std::cerr<<"@LOG::misalignment detection enabled"<<std::endl;
        converter.set_misalign();

        if(args.is_set(Opt::OUTLIER_STDV)){
            std::cerr<<"@LOG::misalignment is set to be detected with standard deviation at "<<args.get_int(Opt::OUTLIER_STDV)<<std::endl;
            converter.set_stdv(args.get_int(Opt::OUTLIER_STDV));
        }
        else{
            std::cerr<<"@LOG::misalignment is set to be detected with raw threshold at "<<args.get_int(Opt::OUTLIER_RAW)<<std::endl;
            converter.set_raw(args.get_int(Opt::OUTLIER_RAW));
        }
    }

    if(inputAlFP=="-"){
        std::cerr<<"@LOG::Begin pre-loading data from stream"<<std::endl;
        converter.precompute_save(args.get_int(Opt::NUM_READS_PRECOMP));
        std::cerr<<"@LOG::Done pre-loading data from stream"<<std::endl;
    }
    else{
        std::cerr<<"@LOG::Begin pre-loading data"<<std::endl;
        converter.precompute(args.get_int(Opt::PERCENT));
        std::cerr<<"@LOG::Done pre-loading data"<<std::endl;
    }

    std::cerr<<"@LOG::Begin Translating Coordinates"<<std::endl;
    converter.convert_coords();
    std::cerr<<"@LOG::Done Translating Coordinates"<<std::endl;
    if(inputAlFP=="-"){
        std::cerr<<"@LOG::Begin Translating Coordinates of precomputed reads"<<std::endl;
        converter.convert_coords_precomp();
        std::cerr<<"@LOG::Done Translating Coordinates of precomputed reads"<<std::endl;
    }

    std::string out_abund_fname = args.get_string(Opt::OUT_AL);
    out_abund_fname.append(".loc.abunds");
    converter.print_abundances(out_abund_fname);

    converter.print_stats();

    return 0;
}
