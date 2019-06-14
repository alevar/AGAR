//
// Created by varabyou on 6/6/19.
//

#include "map2gff_salmon.h"
#include "tokenize.h"

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

Map2GFF_SALMON::Map2GFF_SALMON(const std::string& tlstFP, const std::string& alFP,const std::string& abundFP,const std::string& genome_headerFP,const std::string& outFP,const int& threads, const int& num_trans){
    this->tlstFP = tlstFP;
    this->alFP   = alFP;
    this->abundFP= abundFP;
    this->infoFP = infoFP;
    this->genome_headerFP = genome_headerFP;
    this->outFP = outFP;

    this->umap.set_outFP(outFP);

    this->numThreads=threads;
    this->numTranscripts = num_trans;

    // now load the appropriate data
    this->load_transcriptome();
    this->load_abundances();

    // now initialize the sam file
    this->al=hts_open(this->alFP.c_str(),"r");
    this->al_hdr=sam_hdr_read(this->al); // read the current alignment header
//    this->curAl=bam_init1(); // initialize the alignment record

    genome_al=hts_open(this->genome_headerFP.c_str(),"r");
    genome_al_hdr=sam_hdr_read(genome_al);

    this->outSAM=sam_open(outFP.c_str(),"wb");
    this->outSAM_header = bam_hdr_init();

    // TODO: needs to handle both BAM and SAM input
    outSAM_header=bam_hdr_dup(genome_al_hdr);
    int ret = sam_hdr_write(outSAM,outSAM_header);
    for (int i=0;i<genome_al_hdr->n_targets;++i){
        ref_to_id[genome_al_hdr->target_name[i]]=i;
    }
}

Map2GFF_SALMON::~Map2GFF_SALMON() {
//    bam_destroy1(this->curAl);
    bam_hdr_destroy(this->al_hdr);
    bam_hdr_destroy(this->genome_al_hdr);
    sam_close(this->al);
    sam_close(this->genome_al);
    bam_hdr_destroy(outSAM_header);
    sam_close(outSAM);
}

// this function loads the transcriptome from a tlst file, keeping the name of the transcript, position as well as the ID of the transcript
void Map2GFF_SALMON::load_transcriptome(){
    std::ifstream tstream;
    tstream.open(this->tlstFP.c_str(),std::ios::in);
    if(!tstream.good()){
        std::cerr<<"FATAL: Couldn't open transcript data: "<<tlstFP<<std::endl;
        exit(1);
    }
    std::ios::sync_with_stdio(false);
    std::cerr<<"Reading the transcript data: "<<tlstFP<<std::endl;
    this->transcriptome = std::vector<GffTranscript>(this->numTranscripts);

    std::string tline;
    while(std::getline(tstream,tline)){
        GffTranscript t(tline);
        transcriptome[t.numID] = t;
    }
    std::cerr<<"Loaded the transcript data"<<std::endl;
}

// this function takes in the abundance estimation from salmon and augments the transcriptome index with the data
void Map2GFF_SALMON::load_abundances(){
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
    std::getline(abundstream,aline);
    int count=0;
    int tid;
    float abundance;
    while (std::getline(abundstream,aline)) {
        // given a line we need to extract the name and the TPM
        linestream = new std::stringstream(aline);
        std::getline(*linestream,col,'\t');
        tid = std::atoi(col.c_str());
        // now need to get the abundance
        std::getline(*linestream,col,'\t'); // skip second column
        std::getline(*linestream,col,'\t'); // skip third column
        std::getline(*linestream,col,'\t'); // abundance here
        abundance = std::atof(col.c_str());
        this->transcriptome[tid].abundance = abundance;
        delete linestream;
    }
    abundstream.close();
    std::cerr<<"Loaded transcript abundance data"<<std::endl;
}

void Map2GFF_SALMON::convert_coords(){
    GffTranscript *p_trans=NULL;
    GffTranscript *mate_p_trans=NULL;

    bool start=true;

    bam1_t *curAl = bam_init1(); // initialize the alignment record

    while(sam_read1(al,al_hdr,curAl)>0) {
        if (curAl->core.flag & 4) { // if read is unmapped
            // output unaligned reads for hisat2 to realign
            this->umap.insert(curAl);
            continue;
        }

        // otherwise we proceed to evaluate the reads accordingly
        std::string newReadName = bam_get_qname(curAl);

        // TODO: need to output proper log files so that we can detect when something goes wrong when realigning GTEx

        // first check if belongs to a valid pair
        if(this->has_valid_mate(curAl)){ // belongs to a valid pair
            this->process_pair(curAl);
        }
        else{ // does not belong to a valid pair
            this->process_single(curAl);
        }
    }
    bam_destroy1(curAl);
}

bool Map2GFF_SALMON::has_valid_mate(bam1_t *curAl){
    return (curAl->core.flag & 0x1) && // belongs to a pair
           (curAl->core.flag & 0x2) && // mapped as a pair
          !(curAl->core.flag & 0x4) && // is mapped
          !(curAl->core.flag & 0x8);   // mate is mapped
}

bool Map2GFF_SALMON::get_read_start(GVec<GSeg>& exon_list,size_t gff_start,size_t& genome_start, int& exon_idx){
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
int Map2GFF_SALMON::convert_cigar(int i,int cur_intron_len,int miss_length,GSeg *next_exon,int match_length,
        GVec<GSeg>& exon_list,int &num_cigars,int read_start,bam1_t* curAl,int cigars[MAX_CIGARS]){

    int cur_total_pos=read_start; // same as cur_pos but includes the soft clipping bases
    int cur_pos=read_start;
    for (uint8_t c=0;c<curAl->core.n_cigar;++c){
        uint32_t *cigar_full=bam_get_cigar(curAl);
        int opcode=bam_cigar_op(cigar_full[c]);
        int length=bam_cigar_oplen(cigar_full[c]);

        if (opcode==BAM_CINS){
            cigars[num_cigars]=opcode|(length<<BAM_CIGAR_SHIFT);
            ++num_cigars;
        }
        if (opcode==BAM_CSOFT_CLIP){
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
                cigars[num_cigars]=opcode | (remaining_length <<BAM_CIGAR_SHIFT);
                ++num_cigars;
                cur_pos+=remaining_length;
                cur_total_pos+=remaining_length;
                break;
            }
            else if (cur_pos >= (int)cur_exon.start && cur_pos+remaining_length-1>(int)cur_exon.end){
                match_length=(int)cur_exon.end-cur_pos+1;
                if (match_length>0){
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
                cur_intron_len+=miss_length;

                cigars[num_cigars]=BAM_CREF_SKIP | (miss_length <<BAM_CIGAR_SHIFT);
                ++num_cigars;

                cur_pos+=match_length+miss_length;
                cur_total_pos+=match_length+miss_length;

                remaining_length-=match_length;
                assert(cur_pos == (int)next_exon->start);
            }
        }
    }
    return 1;
}

// TODO: store all of the computed information about genomic positions within a structure which can be easily created and passed as a single argument to any function
void Map2GFF_SALMON::add_cigar(bam1_t *curAl,int num_cigars,int* cigars){
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

void Map2GFF_SALMON::add_aux(bam1_t *curAl,char xs){
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

    // NH tag should always stay at 1 since we only output one mappng even for multimappers
    // due to the likelihood evaluation
    int nh=1;
    uint8_t* ptr_nh_1=bam_aux_get(curAl,"NH");
    if(ptr_nh_1){
        bam_aux_del(curAl,ptr_nh_1);
    }
    bam_aux_append(curAl,"NH",'i',4,(uint8_t*)&nh);
}

// this function modifies the flags for the alignment
void Map2GFF_SALMON::fix_flag(bam1_t *curAl){
    // fix secondary to primary since we only output one mapping for each read
    curAl->core.flag &= ~BAM_FSECONDARY;
}

// this function performs analysis of the genomic coordinates of a read
// and notifies whether the read needs to be outputted or not (1 or 0)
int Map2GFF_SALMON::collapse_genomic(bam1_t *curAl,size_t cigar_hash){
    // TODO: need to re-implement outputting of unique genomic coordinates per read
    //      consider whether cigar string information is needed

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
int Map2GFF_SALMON::collapse_genomic(bam1_t *curAl, bam1_t *mateAl,size_t cigar_hash,size_t mate_cigar_hash){
    // TODO: need to re-implement outputting of unique genomic coordinates per read
    //      consider whether cigar string information is needed

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

void Map2GFF_SALMON::process_pair(bam1_t *curAl) {
    // need a queue to hold and release pairs
    bam1_t* mate = bam_init1();
    int ret = this->pairs.add(curAl,mate);
    size_t cigar_hash,mate_cigar_hash;
    if(!ret){ // mate returned need to write information out
        return;
    }

    cigar_hash = process_read(curAl);
    mate_cigar_hash = process_read(mate);

    int ret_val = collapse_genomic(curAl,mate,cigar_hash,mate_cigar_hash);

    if(ret_val) {
        finish_read(curAl);
        finish_read(mate);
    }
}

void Map2GFF_SALMON::process_single(bam1_t *curAl){
    size_t cigar_hash = this->process_read(curAl);

    int ret_val = collapse_genomic(curAl,cigar_hash);
    if(ret_val) {
        this->finish_read(curAl);
    }
}

size_t Map2GFF_SALMON::process_read(bam1_t *curAl) {
    // let's deal with this case first since there is less stuff
    int target_name = atoi(al_hdr->target_name[curAl->core.tid]); // name of the transcript from the input alignment
    GffTranscript& p_trans = transcriptome[target_name]; // get the transcript
    GVec<GSeg>& exon_list=p_trans.exons; // get exons

    GSeg *next_exon=nullptr;
    int cigars[MAX_CIGARS];
    int match_length,miss_length,cur_intron_len=0,i=0,num_cigars=0;
    size_t read_start=0;

    // first find the genomic read start
    bool ret_val = Map2GFF_SALMON::get_read_start(exon_list,curAl->core.pos,read_start,i);
    if(!ret_val){
        std::cerr<<"Can not get the genomic read start"<<std::endl;
        exit(1);
    }

    // now get mate read start
    int match_length_mate,miss_length_mate,cur_intron_len_mate=0,i_mate=0,num_cigars_mate=0;
    size_t read_start_mate=0;
    ret_val = Map2GFF_SALMON::get_read_start(exon_list,curAl->core.mpos,read_start_mate,i_mate);
    if(!ret_val){
        std::cerr<<"Can not get the genomic read start of the mate"<<std::endl;
        exit(1);
    }

    // secondly build a new cigar string
    int cur_cigar_full[MAX_CIGARS];
    memcpy(cur_cigar_full, bam_get_cigar(curAl), curAl->core.n_cigar);

    ret_val = Map2GFF_SALMON::convert_cigar(i, cur_intron_len, miss_length, next_exon, match_length,
                                            exon_list, num_cigars, read_start, curAl, cigars);
    if (!ret_val) {
        std::cerr << "Can not create a new cigar string for the single read" << std::endl;
        exit(1);
    }

    // before we change anything we need to write the aux data which depends on the current information in the record
    // write data to the optional tags
    add_aux(curAl, p_trans.strand);

    // fix flags
    fix_flag(curAl);

    // now need to deal with the rest
    curAl->core.tid = ref_to_id[p_trans.refID]; // assign new reference
    curAl->core.pos = read_start - 1; // assign new position

    // convert the mate information for single reads without a concordantly mapped mate
    curAl->core.mtid = ref_to_id[p_trans.refID];
    curAl->core.mpos = read_start_mate - 1;

    // assign new cigar string to the record
    add_cigar(curAl, num_cigars, cigars);

    return cigar2hash(cigars,num_cigars);
}

// write a read without a mate
void Map2GFF_SALMON::finish_read(bam1_t *curAl){
    int ret_val = sam_write1(this->outSAM, genome_al_hdr, curAl);
}
