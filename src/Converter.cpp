//
// Created by varabyou on 6/6/19.
//

#include <sys/stat.h>

#include "Converter.h"
#include "tokenize.h"

// TODO: report how many multimappers were detected in the run

// TODO: can be made better hash
size_t cigar2hash(const int cigars[MAX_CIGARS],int n_cigar){
    size_t hash = 0;
    for (uint8_t c=0;c<n_cigar;++c){
        hash = hash * 31 + std::hash<int>()(bam_cigar_op(cigars[c]));
        hash = hash * 31 + std::hash<int>()(bam_cigar_oplen(cigars[c]));
    }
    return hash;
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

Converter::Converter(const std::string& alFP,const std::string& outFP,const std::string& index_base,const int& threads,bool multi){
    this->numThreads=threads;
    this->outFP = outFP;
    this->alFP = alFP;

    // now load the index
    this->load_index(index_base,multi);

    // now initialize the sam file
    struct stat buffer{};
    if(stat (alFP.c_str(), &buffer) != 0){ // if file does not exists
        std::cerr<<"Alignment file: "<<alFP<<" is not found."<<std::endl;
        exit(1);
    }
    this->al=hts_open(alFP.c_str(),"r");
    this->al_hdr=sam_hdr_read(this->al); // read the current alignment header

    genome_al_hdr=sam_hdr_read(genome_al);

    this->outSAM=sam_open(outFP.c_str(),"wb");
    this->outSAM_header = bam_hdr_init();

    // TODO: needs to handle both BAM and SAM input
    this->outSAM_header=bam_hdr_dup(genome_al_hdr);
    int ret = sam_hdr_write(outSAM,this->outSAM_header);
    for (int i=0;i<genome_al_hdr->n_targets;++i){
        ref_to_id[genome_al_hdr->target_name[i]]=i;
    }
}

Converter::~Converter() {
//    bam_destroy1(this->curAl);
    bam_hdr_destroy(this->al_hdr);
    bam_hdr_destroy(this->genome_al_hdr);
    sam_close(this->al);
    sam_close(this->genome_al);
    sam_close(outSAM);
    bam_hdr_destroy(this->outSAM_header);

    this->unal_r1.close();
    this->unal_r2.close();
    this->unal_s.close();
}

void Converter::set_unaligned(){
    this->umap.set_outFP(this->outFP);
    this->unaligned_mode = true;
}

void Converter::set_k1(){
    this->k1_mode = true;
}

void Converter::set_fraglen(int fraglen) {
    this->fraglen = fraglen;
    if(!this->multi){
        return; // no need to set the argument if multimappers are not being evaluated
    }
    this->mmap.set_fraglen(fraglen);
}

void Converter::set_num_multi(int num_multi) {
    this->mmap.set_num_multi(num_multi);
}

void Converter::set_all_multi() {
    this->mmap.set_all_multi();
}

void Converter::set_misalign() {
    this->detect_misalign = true;

    // init output files
    this->unal_r1_fname = outFP;
    this->unal_r1_fname.append(".unal_r1.fastq");
    this->unal_r1.open(this->unal_r1_fname, std::ofstream::out | std::ofstream::trunc);

    this->unal_r2_fname = outFP;
    this->unal_r2_fname.append(".unal_r2.fastq");
    this->unal_r2.open(this->unal_r2_fname, std::ofstream::out | std::ofstream::trunc);

    this->unal_s_fname = outFP;
    this->unal_s_fname.append(".unal_s.fastq");
    this->unal_s.open(this->unal_s_fname, std::ofstream::out | std::ofstream::trunc);
}

void Converter::load_info(const std::string& info_fname){
    // read file to get important stats
    struct stat buffer{};
    if(stat (info_fname.c_str(), &buffer) != 0){ // if file does not exists
        std::cerr<<"Info file: "<<info_fname<<" is not found. Check that correct index is provided."<<std::endl;
        exit(1);
    }
    // read file and save important info
    std::cerr<<"Reading the info file: "<<info_fname<<std::endl;
    std::string iline;
    std::ifstream infostream(info_fname);

    std::getline(infostream,iline);
    this->numTranscripts = std::stoi(iline);

    std::getline(infostream,iline);
    this->maxLocID = std::stoi(iline);
    this->loci = Loci(maxLocID+1);

    std::getline(infostream,iline);
    this->kmerlen = std::stoi(iline);

    infostream.close();
    std::cerr<<"Loaded info data"<<std::endl;
}

// this function parses the .tlst file and inserts entries into the index
void Converter::_load_transcriptome(std::ifstream& tlstfp,char* buffer){
    int k = tlstfp.gcount();
    uint32_t tid=0,start=0,end=0,chr=0,locus=0;
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
                this->transcriptome[transcript.gffID] = transcript;
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
                        std::cerr<<"should never happen _load_transcriptome"<<std::endl;
                        exit(1);
                }
                break;
            default:
                std::cerr<<"unrecognized character"<<std::endl;
                std::cerr<<buffer[i]<<std::endl;
                exit(1);
        }
    }
    if(!transcript.empty){
        transcript.add_exon(exon);
        this->transcriptome[transcript.gffID] = transcript;
    }
}

void Converter::print_transcriptome(){
    std::cerr<<"printing transcriptome"<<std::endl;
    for(auto &t : this->transcriptome){
        t.print();
    }
    std::cerr<<"done printing transcriptome"<<std::endl;
}

void Converter::load_transcriptome(const std::string &tlst_fname) {
    // first initiate the transcriptome array to hold the transcripts
    this->transcriptome = std::vector<GffTranscript>(this->numTranscripts+1);

    struct stat buffer{};
    if(stat (tlst_fname.c_str(), &buffer) != 0){ // if file does not exists
        std::cerr<<"Locus file: "<<tlst_fname<<" is not found. Check that correct index is provided."<<std::endl;
        exit(1);
    }
    std::cerr<<"loading locus data from: "<<tlst_fname<<std::endl;
    std::ifstream tlstfp(tlst_fname,std::ifstream::binary);
    if(tlstfp){
        tlstfp.seekg(0,tlstfp.end);
        int size = tlstfp.tellg();
        tlstfp.seekg(0,tlstfp.beg);
        char *buffer = new char[size];
        if(tlstfp.read(buffer,size)){
            this->_load_transcriptome(tlstfp,buffer);
        }
        delete[] buffer;
        std::cerr<<"finished loading the locus data"<<std::endl;
    }
    else{
        std::cerr<<"failed to open the locus file"<<std::endl;
    }
    tlstfp.close();
}

void Converter::load_genome_header(const std::string& genome_header_fname){
    struct stat buffer{};
    if(stat (genome_header_fname.c_str(), &buffer) != 0){ // if file does not exists
        std::cerr<<"Genome Header file: "<<genome_header_fname<<" is not found. Check that correct index is provided."<<std::endl;
        exit(1);
    }
    this->genome_al=hts_open(genome_header_fname.c_str(),"r");
}

void Converter::load_index(const std::string& index_base,bool multi){
    // first verify that all required files are present
    this->multi = multi;

    // load information
    std::string info_fname(index_base);
    info_fname.append(".info");
    this->load_info(info_fname);

    // load transcriptome
    std::string tlst_fname(index_base);
    tlst_fname.append(".tlst");
    this->load_transcriptome(tlst_fname);

    // load the locus data
    std::string glst_fname(index_base);
    glst_fname.append(".glst");
    this->loci.load(glst_fname);

    // load genome header
    std::string genome_header_fname(index_base);
    genome_header_fname.append(".genome_header");
    this->load_genome_header(genome_header_fname);

    if(this->multi){
        std::string multi_fname(index_base);
        multi_fname.append(".multi");
        this->load_multi(multi_fname);
        // now we can also set pointers to the related objects
        this->mmap.set_loci(&(this->loci));
        this->mmap.set_transcriptome(&(this->transcriptome));
    }
}

// this function takes in the abundance estimation from salmon and augments the transcriptome index with the data
// TODO: this function needs to be written differently to handle not just the salmon input
void Converter::load_abundances(const std::string& abundFP){
    this->abund = true;
    this->loci.set_precomputed_abundances();
    std::ifstream abundstream;
    std::stringstream *linestream;
    abundstream.open(abundFP.c_str(),std::ios::in);
    if (!abundstream.good()){
        std::cerr<<"FATAL: Couldn't open transcript abundance file: "<<abundFP<<std::endl;
        exit(1);
    }
    std::ios::sync_with_stdio(false);

    std::cerr<<"Reading the transcript abundance file: "<<abundFP<<std::endl;
    std::string aline,col;
    std::getline(abundstream,aline); // skip the header from salmon
    int locid;
    float abundance;
    while (std::getline(abundstream,aline)) {
        // given a line we need to extract the name and the TPM
        linestream = new std::stringstream(aline);
        std::getline(*linestream,col,'\t');
        locid = std::atoi(col.c_str());
        // now need to get the abundance
        std::getline(*linestream,col,'\t'); // skip second column
        std::getline(*linestream,col,'\t'); // skip third column
        std::getline(*linestream,col,'\t'); // abundance here
        abundance = std::atof(col.c_str());
        this->loci[locid].set_abund(abundance);
        delete linestream;
    }
    abundstream.close();
    std::cerr<<"Loaded transcript abundance data"<<std::endl;
}

// this metod performs evaluation of errors detected by the aligner
// and decides whether the read needs to be processed or not
bool Converter::evaluate_errors(bam1_t *curAl){ // return true if the read passes the error check
    return this->errorCheck.add_read(curAl);
}

// this function cycles through a section of the bam file loads some preliminary data
// necessary for thecoorect parsing of the errors, multimappers, etc
void Converter::precompute(int perc){
    bam1_t *curAl = bam_init1(); // initialize the alignment record

    int counter=0;
    while(sam_read1(al,al_hdr,curAl)>0) { // only perfom if unaligned flag is set to true
        if(counter%perc==0){
            bool ret = Converter::evaluate_errors(curAl); // add to the error checks
        }
        counter++;
    }
    bam_destroy1(curAl);

    bam_hdr_destroy(this->al_hdr);
    sam_close(this->al);
    this->al=hts_open(this->alFP.c_str(),"r");
    this->al_hdr=sam_hdr_read(this->al); // read the current alignment header
}

// this function writes reads as unaligned
// implementation is mostly from the bam2fq
void Converter::write_unaligned(bam1_t *curAl){
    this->unal_s << "@" << bam_get_qname(curAl) << "\n";
    int8_t *buf;
    buf = NULL;
    size_t max_buf;
    max_buf = 0;
    int i;
    int32_t qlen = curAl->core.l_qseq;
    assert(qlen >= 0);
    uint8_t* seq;
    uint8_t* qual = bam_get_qual(curAl);
    const uint8_t *oq = bam_aux_get(curAl, "OQ");
    bool has_qual = (qual[0] != 0xff || (oq)); // test if there is quality

    if (max_buf < qlen + 1) {
        max_buf = qlen + 1;
        kroundup32(max_buf);
        buf = static_cast<int8_t *>(realloc(buf, max_buf));
        if (buf == NULL) {
            fprintf(stderr, "Out of memory");
            exit(1);
        }
    }

    buf[qlen] = '\0';
    seq = bam_get_seq(curAl);
    for (i = 0; i < qlen; ++i)
        buf[i] = bam_seqi(seq, i);
    if (curAl->core.flag & BAM_FREVERSE) { // reverse complement
        for (i = 0; i < qlen>>1; ++i) {
            int8_t t = seq_comp_table[buf[qlen - 1 - i]];
            buf[qlen - 1 - i] = seq_comp_table[buf[i]];
            buf[i] = t;
        }
        if (qlen&1) buf[i] = seq_comp_table[buf[i]];
    }
    for (i = 0; i < qlen; ++i) {
        buf[i] = seq_nt16_str[buf[i]];
    }
    this->unal_s << (char*)buf<<std::endl;

    if (has_qual) {
        // Write quality
        this->unal_s << "+\n";
        if (oq) memcpy(buf, oq + 1, qlen);
        else {
            for (i = 0; i < qlen; ++i)
                buf[i] = 33 + qual[i];
        }
        if (curAl->core.flag & BAM_FREVERSE) { // reverse
            for (i = 0; i < qlen>>1; ++i) {
                int8_t t = buf[qlen - 1 - i];
                buf[qlen - 1 - i] = buf[i];
                buf[i] = t;
            }
        }
        this->unal_s << (char*)buf<<std::endl;
    }
}

void Converter::convert_coords(){
    bam1_t *curAl = bam_init1(); // initialize the alignment record

    while(sam_read1(al,al_hdr,curAl)>0) { // only perfom if unaligned flag is set to true
        if (this->unaligned_mode && curAl->core.flag & 4) { // if read is unmapped
            // output unaligned reads for hisat2 to realign
            this->umap.insert(curAl);
            continue;
        }

        // TODO: need to output proper log files so that we can detect when something goes wrong when realigning GTEx

        // first evaluate the error-rate of the read
        bool ret = Converter::evaluate_errors(curAl);
        // if the read is aligned as part of a pair and the current mate does not pass the error check
        // it then needs to be processed as a singleton (if the other mate passes the error check)
        if(!ret){
            // need to write to the output fasta file the detected poorly aligned reads for realignment
            Converter::write_unaligned(curAl); // TODO: needs to properly handle pairs
            continue; // didn't pass the error check - continue to the next read
            // TODO: needs to remember the info for mated reads (or to be evaluated as part of the pair)
        }

        // otherwise we proceed to evaluate the reads accordingly
        // first check if belongs to a valid pair
        if(this->has_valid_mate(curAl)){ // belongs to a valid pair
            this->process_pair(curAl);
        }
//        else if(this->has_mate(curAl)){ // TODO: this is the case when mates are aligned to different transcripts on the same locus, or can be re-united
//            this-
//        }
        else{ // does not belong to a valid pair
            this->process_single(curAl);
        }
    }
    bam_destroy1(curAl);
}

bool Converter::has_valid_mate(bam1_t *curAl){
    return (curAl->core.flag & 0x1) && // belongs to a pair
           (curAl->core.flag & 0x2) && // mapped as a pair
          !(curAl->core.flag & 0x4) && // is mapped
          !(curAl->core.flag & 0x8); // mate is mapped
}

bool Converter::has_mate(bam1_t *curAl){ // TODO: needs to be used to properly handle pairs that are mapped discordantly
    return (curAl->core.flag & 0x1) && // belongs to a pair
           !(curAl->core.flag & 0x4) && // is mapped
           !(curAl->core.flag & 0x8); // mate is mapped
}

bool Converter::get_read_start(GVec<GSeg>& exon_list,int32_t gff_start,int32_t& genome_start, int& exon_idx){
    if(gff_start==-1){ // conforming to SAM specification where 0 pos PNEXT means it is not set (for bowtie means reads are not paired)
        genome_start = -1;
        return true;
    }
    const GSeg* cur_exon;
    size_t cur_intron_dist=0;
    size_t trans_start=exon_list[0].start;
    int trans_offset=0;
    for(int i=0;i<exon_list.Count();++i){
        cur_exon=&(exon_list[i]);
        trans_offset=trans_start+cur_intron_dist;

        if (gff_start>=cur_exon->start-trans_offset && gff_start <=cur_exon->end-trans_offset){
            genome_start =gff_start+trans_start+cur_intron_dist;
            exon_idx=i;
            return true;
        }
        else{
            if (i+1<exon_list.Count()){
                cur_intron_dist+=exon_list[i+1].start-cur_exon->end-1;
            }
            else{
                return false;
            }
        }
    }
    return false;
}

int Converter::convert_cigar(int i,GSeg *next_exon,GVec<GSeg>& exon_list,int &num_cigars,int read_start,bam1_t* curAl,int cigars[MAX_CIGARS],Position& pos_obj){
    int miss_length = 0,match_length = 0;
    int cur_total_pos=read_start; // same as cur_pos but includes the soft clipping bases
    int cur_pos=read_start;
    int move = 0;
    for (uint8_t c=0;c<curAl->core.n_cigar;++c){
        uint32_t *cigar_full=bam_get_cigar(curAl);
        int opcode=bam_cigar_op(cigar_full[c]);
        int length=bam_cigar_oplen(cigar_full[c]);

        if (opcode==BAM_CINS){
            cigars[num_cigars]=opcode|(length<<BAM_CIGAR_SHIFT);
            ++num_cigars;
        }
        if (opcode==BAM_CDEL){
            move+=length;
        }
        if (opcode==BAM_CSOFT_CLIP){
            move+=length;
            cigars[num_cigars]=opcode|(length<<BAM_CIGAR_SHIFT);
            ++num_cigars;
            cur_total_pos+=bam_cigar_oplen(cigars[num_cigars]);
        }
        if (opcode != BAM_CMATCH && opcode != BAM_CDEL){
            continue;
        }
        int remaining_length=length;
        for (;i<exon_list.Count();++i){
            GSeg& cur_exon=exon_list[i];
            if (cur_pos>=(int)cur_exon.start && cur_pos+remaining_length-1<=(int)cur_exon.end){
                move+=remaining_length;
                cigars[num_cigars]=opcode | (remaining_length <<BAM_CIGAR_SHIFT);
                ++num_cigars;
                cur_pos+=remaining_length;
                cur_total_pos+=remaining_length;
                break;
            }
            else if (cur_pos >= (int)cur_exon.start && cur_pos+remaining_length-1>(int)cur_exon.end){
                match_length=(int)cur_exon.end-cur_pos+1;
                if (match_length>0){
                    move+=match_length;
                    cigars[num_cigars]=opcode | (match_length << BAM_CIGAR_SHIFT);
                    ++num_cigars;
                }
                if (i+1>=exon_list.Count()){
                    return 0;
                }
                else{
                    next_exon=&(exon_list[i+1]);
                }
                miss_length=next_exon->start-cur_exon.end-1;

                // when an intron is encountered need to push the move, set the new move to zero and push an intron to the moves
                pos_obj.add_move(move);
                pos_obj.add_move(miss_length+1);
                move = 0;
                cigars[num_cigars]=BAM_CREF_SKIP | (miss_length <<BAM_CIGAR_SHIFT);
                ++num_cigars;

                cur_pos+=match_length+miss_length;
                cur_total_pos+=match_length+miss_length;

                remaining_length-=match_length;
                assert(cur_pos == (int)next_exon->start);
            }
        }
    }
    pos_obj.add_move(move);
    return 1;
}

void Converter::add_cigar(bam1_t *curAl,int num_cigars,int* cigars){
    int old_num_cigars = curAl->core.n_cigar;
    int data_len=curAl->l_data+4*(num_cigars-old_num_cigars);
    int m_data=std::max(data_len,(int)curAl->m_data);
    kroundup32(m_data);

    auto* data = (uint8_t*)calloc(m_data,1);

    int copy1_len = (uint8_t*)bam_get_cigar(curAl) - curAl->data;
    memcpy(data, curAl->data, copy1_len);

    int copy2_len = num_cigars * 4;
    memcpy(data + copy1_len, cigars, copy2_len);

    int copy3_len = curAl->l_data - copy1_len - (old_num_cigars * 4);
    memcpy(data + copy1_len + copy2_len, bam_get_seq(curAl), copy3_len);

    curAl->core.n_cigar = num_cigars;

    free(curAl->data);
    curAl->data = data;
    curAl->l_data = data_len;
    curAl->m_data = m_data;
}

void Converter::change_nh_flag(bam1_t *curAl,int nh){
    uint8_t* ptr_nh_1=bam_aux_get(curAl,"NH");
    if(ptr_nh_1){
        bam_aux_del(curAl,ptr_nh_1);
    }
    bam_aux_append(curAl,"NH",'i',4,(uint8_t*)&nh);
}

void Converter::add_aux(bam1_t *curAl,char xs){
    uint8_t* ptr=bam_aux_get(curAl,"XS");
    if(ptr){
        bam_aux_del(curAl,ptr);
    }
    if (xs=='-'){
        bam_aux_append(curAl,"XS",'A',1,(const unsigned char*)"-");
    }
    if (xs=='+'){
        bam_aux_append(curAl,"XS",'A',1,(const unsigned char*)"+");
    }

    // add optional tag with transcript ID
    uint8_t* ptr_op=bam_aux_get(curAl,"OP");
    if(ptr_op){
        bam_aux_del(curAl,ptr_op);
    }
    bam_aux_append(curAl,"OP",'i',4,(uint8_t*)&curAl->core.tid);
}

// this function modifies the flags for the alignment
void Converter::fix_flag(bam1_t *curAl){
    // fix secondary to primary since we only output one mapping for each read
    curAl->core.flag &= ~BAM_FSECONDARY;
}

// this function performs analysis of the genomic coordinates of a read
// and notifies whether the read needs to be outputted or not (1 or 0)
int Converter::collapse_genomic(bam1_t *curAl,size_t cigar_hash){
    // this function needs several things.
    // 1. read start and mate start postitions
    // 2. strand
    // 3. cigar string (can we pass it efficiently and keep it in the map efficiently???
    // that's why we did the conversion to coordinates before, so that we can use the coordinate string as a unique identifier

    // TODO: here we assume that the input file is sorted by name
    //     this is important for salmon mode in order to collapse reads
    //     however in bowtie mode, this will be unimportant

    return collapser.add(curAl,cigar_hash);
}

// this function performs analysis of the genomic coordinates of a read
// and notifies whether the read needs to be outputted or not (1 or 0)
int Converter::collapse_genomic(bam1_t *curAl, bam1_t *mateAl,size_t cigar_hash,size_t mate_cigar_hash){

    // this function needs several things.
    // 1. read start and mate start postitions
    // 2. strand
    // 3. cigar string (can we pass it efficiently and keep it in the map efficiently???
    // that's why we did the conversion to coordinates before, so that we can use the coordinate string as a unique identifier

    // TODO: here we assume that the input file is sorted by name
    //     this is important for salmon mode in order to collapse reads
    //     however in bowtie mode, this will be unimportant

    return collapser.add(curAl,mateAl,cigar_hash,mate_cigar_hash);
}

void Converter::process_pair(bam1_t *curAl) {
    // need a queue to hold and release pairs
    bam1_t* mate = bam_init1();
    int ret = this->pairs.add(curAl,mate);
    size_t cigar_hash,mate_cigar_hash;
    if(!ret){ // mate not found
        bam_destroy1(mate);
        return;
    }
    Position cur_pos,cur_pos_mate;
    int cigars[MAX_CIGARS];
    int num_cigars=0;
    int cigars_mate[MAX_CIGARS];
    int num_cigars_mate=0;
    cigar_hash = process_read(curAl,cur_pos,cigars,num_cigars);
    mate_cigar_hash = process_read(mate,cur_pos_mate,cigars_mate,num_cigars_mate);

    if(!this->k1_mode){
        int ret_val = collapse_genomic(curAl,mate,cigar_hash,mate_cigar_hash);
        if(!ret_val) {
            bam_destroy1(mate);
            return;
        }
    }

    this->evaluate_multimappers_pair(curAl,mate,cur_pos,cur_pos_mate,cigars,cigars_mate,num_cigars,num_cigars_mate); // also finishes the read

    bam_destroy1(mate);
}

void Converter::process_single(bam1_t *curAl){
    Position cur_pos;
    int cigars[MAX_CIGARS];
    int num_cigars=0;
    size_t cigar_hash = this->process_read(curAl,cur_pos,cigars,num_cigars);

    if(!this->k1_mode){
        int ret_val = collapse_genomic(curAl,cigar_hash);
        if(!ret_val) {
            return;
        }
    }

    this->evaluate_multimappers(curAl,cur_pos,cigars,num_cigars); // also finishes the read
}

void Converter::add_multi_tag(bam1_t* curAl){
    uint8_t* ptr_op=bam_aux_get(curAl,"ZZ");
    if(ptr_op){
        bam_aux_del(curAl,ptr_op);
    }
    bam_aux_append(curAl,"ZZ",'A',1,(const unsigned char*)"+");
}

void Converter::evaluate_multimappers_pair(bam1_t *curAl,bam1_t* curAl_mate,Position &cur_pos,Position &cur_pos_mate,
                                                int *cigars,int *cigars_mate,int &num_cigars,int &num_cigars_mate) {
    bool unique;
    std::vector<Position> res_pos,res_pos_mate; // holds the results of the multimapper evaluation
    if(!this->abund){ // compute abundance dynamically
        unique = this->mmap.process_pos_pair(cur_pos,cur_pos_mate,this->loci,res_pos,res_pos_mate);
    }
    else{ // compute abundance dynamically
        unique = this->mmap.process_pos_pair_precomp(cur_pos,cur_pos_mate,this->loci,res_pos,res_pos_mate);
    }
    if(unique){ // increment abundance
        if(!this->abund) {
            this->loci.add_read(cur_pos.locus);
            this->loci.add_read(cur_pos_mate.locus);
        }
        curAl->core.pos = cur_pos.start-1;
        curAl->core.mpos = cur_pos_mate.start-1;
        curAl_mate->core.pos = cur_pos_mate.start-1;
        curAl_mate->core.mpos = cur_pos.start-1;
        // add already computed cigar to the read
        add_cigar(curAl, num_cigars, cigars); // will be performed afterwards
        add_cigar(curAl_mate, num_cigars_mate, cigars_mate); // will be performed afterwards
        this->finish_read(curAl);
        this->finish_read(curAl_mate);
    }
    else{
        add_multi_tag(curAl);
        add_multi_tag(curAl_mate);
        change_nh_flag(curAl,res_pos.size());
        change_nh_flag(curAl_mate,res_pos.size());
        for(int pos_idx=0;pos_idx<res_pos.size();pos_idx++){
            // increment total abundances of the locus to which the new cur_pos belongs
            if(!this->abund) {
                this->loci.add_read_multi(res_pos[pos_idx].locus);
                this->loci.add_read_multi(res_pos_mate[pos_idx].locus);
            }

            // first get the transcript
            GffTranscript& new_trans = transcriptome[res_pos[pos_idx].transID];
            GVec<GSeg>& exon_list=new_trans.exons;
            GSeg *next_exon=nullptr;
            int32_t read_start=res_pos[pos_idx].start;
            int i=0;
            for(;i<exon_list.Count();i++){
                if(read_start>=exon_list[i].start && read_start<=exon_list[i].end){
                    break;
                }
            }

            GffTranscript& new_trans_mate = transcriptome[res_pos_mate[pos_idx].transID];
            GVec<GSeg>& exon_list_mate=new_trans_mate.exons;
            GSeg *next_exon_mate=nullptr;
            int32_t read_start_mate=res_pos_mate[pos_idx].start;
            int i_mate=0;
            for(;i_mate<exon_list_mate.Count();i_mate++){
                if(read_start_mate>=exon_list_mate[i_mate].start && read_start_mate<=exon_list_mate[i_mate].end){
                    break;
                }
            }

            curAl->core.pos = res_pos[pos_idx].start-1;
            curAl->core.mpos = res_pos_mate[pos_idx].start-1;
            curAl_mate->core.pos = res_pos[pos_idx].start-1;
            curAl_mate->core.mpos = res_pos[pos_idx].start-1;

            // reconvert the cur_pos into a read and output
            num_cigars = 0,num_cigars_mate = 0;
            memset(cigars, 0, MAX_CIGARS);
            memset(cigars_mate, 0, MAX_CIGARS);

            int ret_val = Converter::convert_cigar(i,next_exon,exon_list,num_cigars,read_start,curAl,cigars,res_pos[pos_idx]);
            if (!ret_val) {
                std::cerr << "Can not create a new cigar string for the single read evaluate_multimapper_pair1" << std::endl;
                exit(1);
            }
            add_cigar(curAl, num_cigars, cigars); // will be performed afterwards
            this->finish_read(curAl);

            int ret_val_mate = Converter::convert_cigar(i_mate,next_exon_mate,exon_list_mate,num_cigars_mate,read_start_mate,curAl_mate,cigars_mate,res_pos_mate[pos_idx]);
            if (!ret_val_mate) {
                std::cerr << "Can not create a new cigar string for the single read evaluate_multimapper_pair2" << std::endl;
                exit(1);
            }
            add_cigar(curAl_mate, num_cigars_mate, cigars_mate); // will be performed afterwards
            this->finish_read(curAl_mate);
        }
    }
}

void Converter::evaluate_multimappers(bam1_t* curAl,Position& cur_pos,int cigars[MAX_CIGARS],int &num_cigars){ // TODO: make sure only mapped reads are passed through this function
    bool unique;
    std::vector<Position> res_pos; // holds the results of the multimapper evaluation
    if(!this->abund){ // compute abundance dynamically
        unique = this->mmap.process_pos(cur_pos,this->loci,res_pos);
    }
    else{ // compute abundance dynamically
        unique = this->mmap.process_pos_precomp(cur_pos,this->loci,res_pos);
    }
    if(unique){ // increment abundance
        if(!this->abund){
            this->loci.add_read(cur_pos.locus);
        }
        curAl->core.pos = cur_pos.start-1;
        curAl->core.tid = cur_pos.chr;
        // add already computed cigar to the read
        add_cigar(curAl, num_cigars, cigars); // will be performed afterwards
        this->finish_read(curAl);
    }
    else{
        add_multi_tag(curAl);
        change_nh_flag(curAl,res_pos.size());
        for(auto &v : res_pos){
            // modify the NH tag to include however many multimappers were added

            // increment total abundances of the locus to which the new cur_pos belongs - performed for each locus if multiple positions are reported (all or -k)
            if(!this->abund) {
                this->loci.add_read_multi(v.locus);
            }

            // first get the transcript
            GffTranscript& new_trans = transcriptome[v.transID];
            GVec<GSeg>& exon_list=new_trans.exons;
            GSeg *next_exon=nullptr;
            int32_t read_start=v.start;
            int i=0;
            for(;i<exon_list.Count();i++){
                if(read_start>=exon_list[i].start && read_start<=exon_list[i].end){
                    break;
                }
            }

            curAl->core.pos = v.start-1; // assign new position
            curAl->core.tid = v.chr;

            // reconvert the cur_pos into a read and output
            num_cigars = 0;
            memset(cigars, 0, MAX_CIGARS);

            int ret_val = Converter::convert_cigar(i,next_exon,exon_list,num_cigars,read_start,curAl,cigars,v);
            if (!ret_val) {
                std::cerr << "Can not create a new cigar string for the single read evaluate_multimappers" << std::endl;
                exit(1);
            }

            add_cigar(curAl, num_cigars, cigars); // will be performed afterwards
            this->finish_read(curAl);
        }
    }
}

size_t Converter::process_read(bam1_t *curAl,Position& cur_pos,int cigars[MAX_CIGARS],int &num_cigars) {
    // let's deal with this case first since there is less stuff
    int target_name = atoi(al_hdr->target_name[curAl->core.tid]); // name of the transcript from the input alignment
    GffTranscript& p_trans = transcriptome[target_name]; // get the transcript
    GVec<GSeg>& exon_list=p_trans.exons; // get exons

    GSeg *next_exon=nullptr;
    int i=0;
    int32_t read_start=0;

    // first find the genomic read start
    bool ret_val = Converter::get_read_start(exon_list,curAl->core.pos,read_start,i);
    if(!ret_val){
        std::cerr<<"Can not get the genomic read start"<<std::endl;
        exit(1);
    }

    // now get mate read start
    int target_name_mate = atoi(al_hdr->target_name[curAl->core.mtid]); // name of the transcript from the input alignment
    GffTranscript& p_trans_mate = transcriptome[target_name_mate]; // get the transcript
    GVec<GSeg>& exon_list_mate=p_trans_mate.exons; // get exons
    int i_mate=0;
    int32_t read_start_mate=0;

    ret_val = Converter::get_read_start(exon_list_mate,curAl->core.mpos,read_start_mate,i_mate);
    if(!ret_val){
        std::cerr<<"Can not get the genomic read start of the mate"<<std::endl;
        exit(1);
    }

    cur_pos.set_chr(p_trans.get_refID());
    cur_pos.set_start(read_start);
    cur_pos.set_locus(p_trans.get_geneID());
    cur_pos.set_strand(p_trans.strand);

    // secondly build a new cigar string
    int cur_cigar_full[MAX_CIGARS];
    memcpy(cur_cigar_full, bam_get_cigar(curAl), curAl->core.n_cigar);

    ret_val = Converter::convert_cigar(i,next_exon,exon_list,num_cigars,read_start,curAl,cigars,cur_pos);
    if (!ret_val) {
        std::cerr << "Can not create a new cigar string for the single read from process_read" << std::endl;
        exit(1);
    }

    // before we change anything we need to write the aux data which depends on the current information in the record
    // write data to the optional tags
    add_aux(curAl, p_trans.strand);

    // fix flags
    fix_flag(curAl);

    // now need to deal with the rest
    curAl->core.tid = p_trans.refID; // assign new reference
    curAl->core.pos = read_start - 1; // assign new position

    // convert the mate information for single reads without a concordantly mapped mate
    // unless the read is not paired, in which case the paired information should stay the same
    if(curAl->core.mpos != -1 || curAl->core.mtid != -1){
        curAl->core.mtid = p_trans_mate.refID;
        curAl->core.mpos = read_start_mate - 1;
    }

    // assign new cigar string to the record

    return cigar2hash(cigars,num_cigars);
}

// write a read without a mate
void Converter::finish_read(bam1_t *curAl){
    int ret_val = sam_write1(this->outSAM, genome_al_hdr, curAl);
}

void Converter::load_multi(const std::string& multiFP){
    this->multi = true;
    this->mmap.load(multiFP);
}

void Converter::print_multimappers() {
    this->mmap.print_multimapers();
}