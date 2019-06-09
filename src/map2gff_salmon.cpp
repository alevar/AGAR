//
// Created by varabyou on 6/6/19.
//

#include "map2gff_salmon.h"
#include "tokenize.h"

void tline_parserr(const std::string& tline, std::string add="") {
    std::cerr << "Error at parsing .tlst line " << add << ":"
              << std::endl << '\t' << tline << std::endl;
    exit(1);
}

GffTranscript::GffTranscript(const std::string& tline){

    std::istringstream f(tline);
    std::string token;
    std::vector<std::string> tokens;
    while (std::getline(f, token, ' ')) {
        tokens.push_back(token);
    }

    if (tokens.size()!=4) {
        tline_parserr(tline);
    }
    numID=atoi(tokens[0].c_str());
    gffID=tokens[1];
    refID=tokens[2];
    if (refID.length()<1) {
        tline_parserr(tline, "(refID empty)");
    }
    strand=refID[refID.length()-1];
    if (strand!='-' && strand!='+') {
        tline_parserr(tline, "(invalid strand)");
    }
    refID.erase(refID.length()-1);

    f.clear(); //to reset the std::getline() iterator
    f.str(tokens[3]);
    while (std::getline(f, token, ',')) {
        size_t sp_pos=token.find('-');
        if (sp_pos == std::string::npos) {
            std::string s("(invalid exon str: ");
            s+=token;s+=")";
            tline_parserr(tline, s);
        }
        std::string s_start=token.substr(0,sp_pos);
        std::string s_end=token.substr(sp_pos+1);
        GSeg exon(atoi(s_start.c_str()), atoi(s_end.c_str()));
        if (exon.start==0 || exon.end==0 || exon.end<exon.start) {
            std::string s("(invalid exon: ");
            s+=token;s+=")";
            tline_parserr(tline, s);
        }
        if (start==0 || start>exon.start){
            start=exon.start;
        }
        if (end==0 || end<exon.end){
            end=exon.end;
        }
        exons.Add(exon);
    } //while exons
}

void print_cigar(bam1_t *al){
    for (uint8_t c=0;c<al->core.n_cigar;++c){
        uint32_t *cigar_full=bam_get_cigar(al);
        int opcode=bam_cigar_op(cigar_full[c]);
        int length=bam_cigar_oplen(cigar_full[c]);
        std::cout<<length<<bam_cigar_opchr(opcode);
    }
    std::cout<<std::endl;
}

void print_aux(bam1_t *al) {
    uint8_t *s= bam_get_aux(al);
    uint8_t *sStop = s+bam_get_l_aux(al);
    std::cout<<*sStop<<std::endl;
}

Map2GFF_SALMON::Map2GFF_SALMON(const std::string& tlstFP, const std::string& alFP, const int& threads){
    this->numThreads=threads;

    tlststream.open(tlstFP.c_str(),std::ios::in);
    if (!tlststream.good()){
        std::cerr<<"FATAL: Couldn't open trascript data: "<<tlstFP<<std::endl;
        exit(1);
    }
    std::ios::sync_with_stdio(false);

    std::cout<<"Reading the transcript data: "<<tlstFP<<std::endl;
    std::string tline;
    while (std::getline(tlststream,tline)){
        if (tline.length()>4){
            auto* t=new GffTranscript(tline);
            transcripts.Add(t);
            tidx_to_t[std::to_string(t->numID)]=t;
        }
    }
    tlststream.close();
    std::cout<<"Loaded Transcript Data"<<std::endl;

    // now initialize the sam file
    al=hts_open(alFP.c_str(),"r");
    al_hdr=sam_hdr_read(al); // read the alignment header
    curAl=bam_init1(); // initialize the alignment
}

void Map2GFF_SALMON::convert_coords(const std::string& outFP, const std::string& genome_header){
    genome_al=hts_open(genome_header.c_str(),"r");
    genome_al_hdr=sam_hdr_read(genome_al);

    samFile *outSAM=sam_open(outFP.c_str(),"wb");
    bam_hdr_t *outSAM_header=bam_hdr_init();

    outSAM_header=bam_hdr_dup(genome_al_hdr);
    sam_hdr_write(outSAM,outSAM_header);
    for (int i=0;i<genome_al_hdr->n_targets;++i){
        ref_to_id[genome_al_hdr->target_name[i]]=i;
    }
    bam_hdr_destroy(outSAM_header);

    GffTranscript *p_trans=NULL;
    GffTranscript *mate_p_trans=NULL;

    bool start=true;

    while(sam_read1(al,al_hdr,curAl)>0) {
        std::string newReadName = bam_get_qname(curAl);
        bool unmapped = curAl->core.flag & 4;
        if (unmapped) {
            // this is where we can simply output the reads to the stdandard out to be accepted by hisat2 for realignment onto the genome
            continue;
        }
//        std::cout << newReadName << std::endl;
//        print_cigar(curAl);
//        std::cout << "========" << std::endl;
    }
}