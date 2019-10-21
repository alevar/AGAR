//
// Created by varabyou on 6/6/19.
//

#include <sys/stat.h>

#include "Comparator.h"
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

std::string extract_pg(bam_hdr_t *al_hdr){
    std::string res_pg = "";

    std::string tmpstr(al_hdr->text,al_hdr->l_text); // length optional, but needed if there may be zero's in your data
    std::istringstream is(tmpstr);

    std::string line;
    while (getline(is,line)) {
        if(line.substr(0,3)=="@PG"){
            res_pg += line+"\n";
        }
    }
    return res_pg;
}

void Comparator::set_xs(bam1_t *curAl,uint8_t xs){
    uint8_t* ptr=bam_aux_get(curAl,"XS");
    if(ptr){
        bam_aux_del(curAl,ptr);
    }
    if(xs == '+'){
        bam_aux_append(curAl,"XS",'A',1,(const unsigned char*)"+");
    }
    else if(xs == '-'){
        bam_aux_append(curAl,"XS",'A',1,(const unsigned char*)"-");
    }
    else{
        return; // do not write an XS flag
    }
}

void Comparator::change_nh_flag(bam1_t *curAl,int nh){
    uint8_t* ptr_nh_1=bam_aux_get(curAl,"NH");
    if(ptr_nh_1){
        bam_aux_del(curAl,ptr_nh_1);
    }
    bam_aux_append(curAl,"NH",'i',4,(uint8_t*)&nh);
}

void Comparator::add_aux(bam1_t *curAl,char xs){
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
    int temp = atoi(al_t_hdr->target_name[curAl->core.tid]);
    bam_aux_append(curAl, "OP", 'i', 4, (uint8_t*)&temp);
}

// this function modifies the flags for the alignment
void Comparator::fix_flag(bam1_t *curAl){
    // fix secondary to primary since we only output one mapping for each read
    curAl->core.flag &= ~BAM_FSECONDARY;
}

// appends the transcript id to which the multimapping was done
void Comparator::add_multi_tag(bam1_t* curAl,int new_tid){
    uint8_t* ptr_zz=bam_aux_get(curAl,"ZZ");
    if(ptr_zz){
        bam_aux_del(curAl,ptr_zz);
    }
    bam_aux_append(curAl,"ZZ",'i',4,(uint8_t*)&new_tid);
}

Comparator::Comparator(const std::string& trans_al_fname,const std::string& genome_al_fname,const std::string& out_al_fname,const std::string& index_base,const int& threads,bool multi,std::string cl){
    this->num_threads=threads;
    this->out_al_fname = out_al_fname;
    this->trans_al_fname = trans_al_fname;
    this->genome_al_fname = genome_al_fname;

    // now load the index
    this->load_index(index_base,multi);

    // now initialize the sam file
    struct stat buffer{};
    if(stat (trans_al_fname.c_str(), &buffer) != 0 ){ // if file does not exists
        std::cerr<<"@ERROR::Transcriptome alignment file is not found: "<<trans_al_fname<<std::endl;
        exit(1);
    }
    if(stat (genome_al_fname.c_str(), &buffer) != 0 ){ // if file does not exists
        std::cerr<<"@ERROR::Genome alignment file is not found: "<<genome_al_fname<<std::endl;
        exit(1);
    }

    this->trans_al=hts_open(trans_al_fname.c_str(),"r");
    this->al_t_hdr=sam_hdr_read(this->trans_al); // read the current alignment header
    this->genome_al=hts_open(genome_al_fname.c_str(),"r");
    this->al_g_hdr=sam_hdr_read(this->genome_al); // read the current alignment header

    std::string prev_pgs = extract_pg(al_t_hdr);

    idx_al_hdr=sam_hdr_read(idx_al);

    this->merged_al=sam_open(out_al_fname.c_str(),"wb");
//    merge_headers(merged_al_hdr,al_g_hdr,idx_al_hdr);
    this->merged_al_hdr = bam_hdr_dup(idx_al_hdr);

    int ret = sam_hdr_write(merged_al,this->merged_al_hdr);
    for (int i=0;i<merged_al_hdr->n_targets;++i){
        ref_to_id[merged_al_hdr->target_name[i]]=i;
    }
}

Comparator::~Comparator() {
    bam_hdr_destroy(this->al_t_hdr);
    bam_hdr_destroy(this->al_g_hdr);
    bam_hdr_destroy(this->idx_al_hdr);
    bam_hdr_destroy(this->merged_al_hdr);
    sam_close(this->trans_al);
    sam_close(this->genome_al);
    sam_close(this->idx_al);
    sam_close(this->merged_al);

    free(buf);
}

void Comparator::merge_headers(bam_hdr_t* hdr_out,bam_hdr_t* hdr1,bam_hdr_t* hdr2){
    // set hdr1 as base keeping ids the same for the respective alignments (genomic in case of AGAR)
    hdr_out = bam_hdr_dup(hdr1);
    // append things that are in hdr2 but not in hdr1 to hdr_out
}

uint8_t Comparator::get_pair_status(bam1_t *al1, bam1_t *al2){
    if(al1->core.flag & 0x2 && al2->core.flag & 0x2){
        return STATUS::CONC;
    }
    else if(!(al1->core.flag & 0x4) && !(al2->core.flag & 0x4)){
        return STATUS::DISC;
    }
    else if(al1->core.flag & 0x4 && al2->core.flag & 0x4){
        return STATUS::UNALIGNED_PAIR;
    }
    else if(al1->core.flag & 0x4 && !(al2->core.flag & 0x4)){
        return STATUS::ALIGNED_MIXED_1; // mate was first
    }
    else if(!(al1->core.flag & 0x4) && al2->core.flag & 0x4){
        return STATUS::ALIGNED_MIXED_2;
    }
    else{
        std::cout<<"ERROR STATUS PAIR"<<std::endl;
        return STATUS::ERROR;
    }
}

uint8_t Comparator::get_read_status(bam1_t *al){
    if(al->core.flag & 0x2 && al->core.flag & 0x1){
        return STATUS_READ::CONC_R;
    }
    else if(!(al->core.flag & 0x4) && !(al->core.flag & 0x8)){
        return STATUS_READ::DISC_R;
    }
    else if(al->core.flag & 0x4 && al->core.flag & 0x8){
        return STATUS_READ::UNALIGNED_R;
    }
    else if((al->core.flag & 0x4 && !(al->core.flag & 0x8)) ||
            (al->core.flag & 0x8 && !(al->core.flag & 0x4))){
        return STATUS_READ::MIXED_R; // mate was first
    }
    else{
        std::cout<<"ERROR STATUS PAIR"<<std::endl;
        return STATUS::ERROR;
    }
}

uint32_t Comparator::md2nm(bam1_t* al) {
    return 0; // TODO: add parsing of MD tag if NM is not unavailable
}

uint32_t Comparator::get_nm(bam1_t* al){
    uint8_t* ptr_nm_1=bam_aux_get(al,"NM");
    if(ptr_nm_1){
        return bam_aux2i(ptr_nm_1);
    }
    else{
        return md2nm(al);
    }
}

// one two or none are aligned
int Comparator::get_num_aligned(bam1_t* curAl){
    int num_al = 0;
    if(!(curAl->core.flag & 0x4)){num_al++;}
    if(!(curAl->core.flag & 0x8)){num_al++;}
    return num_al;
}

bool Comparator::can_be_repaired(bam1_t* al1,bam1_t* al2){
    return al1->core.tid == al2->core.tid && std::abs(al1->core.pos - al2->core.pos) <= this->fraglen;
}

uint32_t Comparator::cmp_gt(bam1_t* al1_1,bam1_t* al1_2,bam1_t* al2_1,bam1_t* al2_2){
    if(std::strcmp(bam_get_qname(al1_1),bam_get_qname(al2_1)) == 0){ // same read
        int nm1 = get_nm(al1_1)+get_nm(al1_2);
        int nm2 = get_nm(al2_1)+get_nm(al2_2);
        int na1 = get_num_aligned(al1_1);
        int na2 = get_num_aligned(al2_1);
        // first look at unaligned vs aligned
        // TODO: perform an additional check based on whether alignment is concordant or not
        //   - this requires transcriptomic alignment to be checked if can be repaired if discordant
        if(na1==0 && na2==0){
            return STATUS_CMP_READ::UNALIGNED;
        }
        else if(!can_be_repaired(al2_1,al2_2)){
            return STATUS_CMP_READ::FIRST;
        }
        else if(na1>na2){
            return STATUS_CMP_READ::FIRST;
        }
        else if(na2>na1){
            return STATUS_CMP_READ::SECOND;
        }
        else if(nm1<nm2){
            return STATUS_CMP_READ::FIRST;
        }
        else if(nm2<nm1){
            return STATUS_CMP_READ::SECOND;
        }
        else{
            return STATUS_CMP_READ::BOTH;
        }
    }
    else{
        return STATUS_CMP_READ::ERR_CMP;
    }
}

uint32_t Comparator::cmp_gt(bam1_t* al1,bam1_t* al2){
    if(std::strcmp(bam_get_qname(al1),bam_get_qname(al2)) == 0){ // same read
        int nm1 = get_nm(al1);
        int nm2 = get_nm(al2);
        int na1 = get_num_aligned(al1);
        int na2 = get_num_aligned(al2);
        // first look at unaligned vs aligned
        if(na1==0 && na2==0){
            return STATUS_CMP_READ::UNALIGNED;
        }
        else if(na1>na2){
            return STATUS_CMP_READ::FIRST;
        }
        else if(na2>na1){
            return STATUS_CMP_READ::SECOND;
        }
        else if(nm1<nm2){
            return STATUS_CMP_READ::FIRST;
        }
        else if(nm2<nm1){
            return STATUS_CMP_READ::SECOND;
        }
        else{
            int read_status1 = get_read_status(al1);
            int read_status2 = get_read_status(al2);
            if(read_status1<read_status2){
                return STATUS_CMP_READ::FIRST;
            }
            else if(read_status2<read_status1){
                return STATUS_CMP_READ::SECOND;
            }
            else{
                return STATUS_CMP_READ::BOTH;
            }
        }
    }
    else{
        return STATUS_CMP_READ::ERR_CMP;
    }
}

int Comparator::skip_read(samFile* al,bam_hdr_t* al_hdr,bam1_t* curAl){
    std::string readname = bam_get_qname(curAl);
    int ret;
    while(true){
        ret = sam_read1(al,al_hdr,curAl);
        if(ret<0){return ret;}
        if(std::strcmp(bam_get_qname(curAl),readname.c_str())!=0){
            return ret;
        }
    }
}

void Comparator::change_tids(bam_hdr_t* al_hdr,bam1_t* curAl){
    if(!(curAl->core.flag & 0x4)){ // read is mapped
        this->r2i_it = this->ref_to_id.find(al_hdr->target_name[curAl->core.tid]);
        if(this->r2i_it == this->ref_to_id.end()){
            std::cerr<<"@ERROR:::reference not found: "<<al_hdr->target_name[curAl->core.tid]<<std::endl;
            exit(-1);
        }
        curAl->core.tid = this->r2i_it->second;
    }
    if(!(curAl->core.flag & 0x8)){ // mate is mapped
        this->r2i_it = this->ref_to_id.find(al_hdr->target_name[curAl->core.mtid]);
        if(this->r2i_it == this->ref_to_id.end()){
            std::cerr<<"@ERROR:::mate reference not found: "<<al_hdr->target_name[curAl->core.mtid]<<std::endl;
            exit(-1);
        }
        curAl->core.mtid = this->r2i_it->second;
    }
}

void Comparator::add_genome_tag(bam1_t* curAl){
    uint8_t* ptr_op=bam_aux_get(curAl,"ZG");
    if(ptr_op){
        bam_aux_del(curAl,ptr_op);
    }
    bam_aux_append(curAl,"ZG",'A',1,(const unsigned char*)"+");
}

void Comparator::add_both_tag(bam1_t* curAl){
    uint8_t* ptr_op=bam_aux_get(curAl,"ZB");
    if(ptr_op){
        bam_aux_del(curAl,ptr_op);
    }
    bam_aux_append(curAl,"ZB",'A',1,(const unsigned char*)"+");
}

int Comparator::load_remaining_genome_pair(samFile* genome_al,bam_hdr_t* al_g_hdr,bam1_t* al_g1,bam1_t* al_g2,std::vector<bam1_t*>& al_g_conc){
    change_tids(al_g_hdr,al_g1);
    change_tids(al_g_hdr,al_g2);
    add_genome_tag(al_g1);
    add_genome_tag(al_g2);
    al_g_conc.push_back(bam_dup1(al_g1));
    al_g_conc.push_back(bam_dup1(al_g2));
    std::string readname = bam_get_qname(al_g1);
    int ret;
    while(true){ // write genome alignments for the current read
        ret = sam_read1(genome_al,al_g_hdr,al_g1);
        if(ret<0){return ret;}
        if(std::strcmp(bam_get_qname(al_g1),readname.c_str()) == 0){
            change_tids(al_g_hdr,al_g1);
            add_genome_tag(al_g1);
            al_g_conc.push_back(bam_dup1(al_g1));
        }
        else{ // found next read
            return ret;
        }
    }
}

int Comparator::write_remaining_genome_pair(samFile* al,bam_hdr_t* al_hdr,bam1_t* curAl,bam1_t* mate){
    change_tids(al_hdr,curAl);
    change_tids(al_hdr,mate);
    add_genome_tag(curAl);
    add_genome_tag(mate);
    sam_write1(this->merged_al, merged_al_hdr, curAl);
    sam_write1(this->merged_al, merged_al_hdr, mate);
    std::string readname = bam_get_qname(curAl);
    int ret;
    while(true){ // write genome alignments for the current read
        ret = sam_read1(genome_al,al_g_hdr,curAl);
        if(ret<0){return ret;}
        if(std::strcmp(bam_get_qname(curAl),readname.c_str()) == 0){
            change_tids(al_hdr,curAl);
            add_genome_tag(curAl);
            sam_write1(this->merged_al, this->merged_al_hdr, curAl);
        }
        else{ // found next read
            return ret;
        }
    }
}

int Comparator::write_remaining_genome_read(samFile* al,bam_hdr_t* al_hdr,bam1_t* curAl){
    change_tids(al_hdr,curAl);
    add_genome_tag(curAl);
    sam_write1(this->merged_al, merged_al_hdr, curAl);
    std::string readname = bam_get_qname(curAl);
    int ret;
    while(true){ // write genome alignments for the current read
        ret = sam_read1(genome_al,al_g_hdr,curAl);
        if(ret<0){return ret;}
        if(std::strcmp(bam_get_qname(curAl),readname.c_str()) == 0){
            change_tids(al_hdr,curAl);
            add_genome_tag(curAl);
            sam_write1(this->merged_al, this->merged_al_hdr, curAl);
        }
        else{ // found next read
            return ret;
        }
    }
}

void Comparator::add_T2G_aux(bam1_t *curAl){
    // add optional tag with transcript ID
    uint8_t* ptr_op=bam_aux_get(curAl,"ZT");
    if(ptr_op){
        bam_aux_del(curAl,ptr_op);
    }
    bam_aux_append(curAl,"ZT",'A',1,(const unsigned char*)"+");
}

int Comparator::write_remaining_both_pair(samFile* al,bam_hdr_t* al_hdr,bam1_t* curAl,bam1_t* mate,std::vector<bam1_t*> al_g_conc){
    if(al_g_conc.size()%2!=0){
        std::cerr<<"@ERROR:::Number of pairs is not even: "<<bam_get_qname(curAl)<<std::endl;
        exit(-1);
    }
    int num_genome = al_g_conc.size()/2;

    add_both_tag(curAl);
    add_both_tag(mate);
    add_T2G_aux(curAl); // add T2G tag
    add_T2G_aux(mate); // add T2G tag

    bool check_mate = false;

    int cur_num_multi = process_status(mate,curAl,check_mate,al_hdr,num_genome);
    if(cur_num_multi==0){ // no multimappers
        cur_num_multi=1;
    }

    // now output the genome alignments correcting them with the respective NH
    for(auto & i : al_g_conc){
        change_nh_flag(i,num_genome+cur_num_multi);
        sam_write1(this->merged_al, this->merged_al_hdr, i);
        bam_destroy1(i);
    }

    std::string readname = bam_get_qname(curAl);
    int ret;
    while(true) { // only perfom if unaligned flag is set to true
        ret = sam_read1(al,al_hdr,curAl);
        if(ret<0){return ret;}
        if(std::strcmp(bam_get_qname(curAl),readname.c_str()) == 0){
            std::cerr<<"@ERROR:::this should never happen in write_remaining_both_pair"<<std::endl;
            exit(-1);
            add_both_tag(curAl);
            add_T2G_aux(curAl); // add T2G tag
            if(check_mate){
                int cur_num_multi = process_status(mate,curAl,check_mate,al_hdr,num_genome);
            }
            else{
                bam_copy1(mate,curAl);
                check_mate = true;
            }
        }
        else{ // found next read
            return ret;
        }
    }
}

int Comparator::write_remaining_trans_pair(samFile* al,bam_hdr_t* al_hdr,bam1_t* curAl,bam1_t* mate,bool both){
    if(both){
        add_both_tag(curAl);
        add_both_tag(mate);
    }
    add_T2G_aux(curAl); // add T2G tag
    add_T2G_aux(mate); // add T2G tag

    bool check_mate = false;

    process_status(mate,curAl,check_mate,al_hdr,0);

    std::string readname = bam_get_qname(curAl);
    int ret;
    while(true) { // only perfom if unaligned flag is set to true
        ret = sam_read1(al,al_hdr,curAl);
        if(ret<0){return ret;}
        if(std::strcmp(bam_get_qname(curAl),readname.c_str()) == 0){
            add_both_tag(curAl);
            add_T2G_aux(curAl); // add T2G tag
            if(check_mate){
                process_status(mate,curAl,check_mate,al_hdr,0);
            }
            else{
                bam_copy1(mate,curAl);
                check_mate = true;
            }
        }
        else{ // found next read
            return ret;
        }
    }
}

int Comparator::write_remaining_trans_read(samFile* al,bam_hdr_t* al_hdr,bam1_t* curAl,bam1_t* mate,bool both){
    if(both){
        add_both_tag(curAl);
    }
    bam_copy1(mate,curAl);
    bool check_mate = true;

    std::string readname = bam_get_qname(curAl);
    int ret;
    while(true) { // only perfom if unaligned flag is set to true
        ret = sam_read1(al,al_hdr,curAl);
        if(ret<0){return ret;}
        if(std::strcmp(bam_get_qname(curAl),readname.c_str()) == 0){
            if(both){
                add_both_tag(curAl);
            }
            add_T2G_aux(curAl); // add T2G tag
            if(check_mate){
                process_status(mate,curAl,check_mate,al_hdr,0);
            }
            else{
                bam_copy1(mate,curAl);
                check_mate = true;
            }
        }
        else{ // found next read
            return ret;
        }
    }
}

bool Comparator::repair(bam1_t *curAl,bam1_t *mate){
    // first check if two form a valid pair in genomic space
    if(curAl->core.tid ==  mate->core.tid &&
       std::abs(curAl->core.pos - mate->core.pos) <= this->fraglen){ // valid pair
//        std::cout<<this->fraglen<<"\t"<<std::abs(curAl->core.pos - mate->core.pos)<<std::endl;

        // fix secondary to primary since we only output one mapping for each read
        curAl->core.flag &= ~BAM_FPROPER_PAIR; // always false
        mate->core.flag &= ~BAM_FPROPER_PAIR; // always false
        curAl->core.mpos = mate->core.pos;
        mate->core.mpos = curAl->core.pos;

        return true;
    }
    return false;
}

int Comparator::_process_disc(bam1_t *curAl,bam1_t* mate,bam_hdr_t* al_hdr,int plus_nm){
    size_t cigar_hash,mate_cigar_hash;
    // Now that both mates of an alignment are found - first evaluate the error-rate of the read

    Position cur_pos,cur_pos_mate;
    int cigars[MAX_CIGARS];
    int num_cigars=0;
    int cigars_mate[MAX_CIGARS];
    int num_cigars_mate=0;
    cigar_hash = process_read(curAl,cur_pos,cigars,num_cigars,al_hdr);
    mate_cigar_hash = process_read(mate,cur_pos_mate,cigars_mate,num_cigars_mate,al_hdr);

    bool ret = repair(curAl,mate);
    if(ret){ // correct pair detected
        this->repaired++;

        int cur_num_multi = this->evaluate_multimappers_pair(curAl,mate,cur_pos,cur_pos_mate,cigars,cigars_mate,num_cigars,num_cigars_mate,plus_nm); // also finishes the read
        if(cur_num_multi>0){
            this->num_multi_pair++;
            this->num_multi_hits_pair=this->num_multi_hits_pair+cur_num_multi;
        }

        total_num_pair_al++;
        return cur_num_multi;
    }
    else{
        int cur_num_multi = this->evaluate_multimappers(curAl,cur_pos,cigars,num_cigars,plus_nm); // also finishes the read
        if(cur_num_multi>0){
            this->num_multi++;
            this->num_multi_hits=this->num_multi_hits+cur_num_multi;
        }

        cur_num_multi = this->evaluate_multimappers(mate,cur_pos_mate,cigars_mate,num_cigars_mate,plus_nm); // also finishes the read
        if(cur_num_multi>0){
            this->num_multi++;
            this->num_multi_hits=this->num_multi_hits+cur_num_multi;
        }

        total_num_disc_al++;
        return cur_num_multi;
    }
}

int Comparator::process_disc(bam1_t *curAl,bam1_t *mate,bam_hdr_t* al_hdr,int plus_nm) {
    return _process_disc(curAl,mate,al_hdr,plus_nm);
}

int Comparator::_process_pair(bam1_t *curAl,bam1_t* mate,bam_hdr_t* al_hdr,int plus_nm){
    size_t cigar_hash,mate_cigar_hash;

    Position cur_pos,cur_pos_mate;
    int cigars[MAX_CIGARS];
    int num_cigars=0;
    int cigars_mate[MAX_CIGARS];
    int num_cigars_mate=0;
    cigar_hash = process_read(curAl,cur_pos,cigars,num_cigars,al_hdr);
    mate_cigar_hash = process_read(mate,cur_pos_mate,cigars_mate,num_cigars_mate,al_hdr);

    if(std::strcmp(bam_get_qname(curAl),"read4234734/rna-XM_006723946.2;mate1:1923-2023;mate2:2322-2422")==0){
        std::cout<<"found"<<std::endl;
    }

    int cur_num_multi = this->evaluate_multimappers_pair(curAl,mate,cur_pos,cur_pos_mate,cigars,cigars_mate,num_cigars,num_cigars_mate,plus_nm); // also finishes the read
    if(cur_num_multi>0){
        this->num_multi_pair++;
        this->num_multi_hits_pair=this->num_multi_hits_pair+cur_num_multi;
    }

    total_num_pair_al++;

    return cur_num_multi;
}

int Comparator::process_pair(bam1_t *curAl,bam1_t *mate,bam_hdr_t* al_hdr,int plus_nm) {
    return _process_pair(curAl,mate,al_hdr,plus_nm);
}

size_t Comparator::process_read(bam1_t *curAl,Position& cur_pos,int cigars[MAX_CIGARS],int &num_cigars,bam_hdr_t* al_hdr) {
    // let's deal with this case first since there is less stuff
    int target_name = atoi(al_hdr->target_name[curAl->core.tid]); // name of the transcript from the input alignment
    GffTranscript& p_trans = transcriptome[target_name]; // get the transcript
    GVec<GSeg>& exon_list=p_trans.exons; // get exons

    GSeg *next_exon=nullptr;
    int i=0;
    int32_t read_start=0;

    // first find the genomic read start
    bool ret_val = this->get_read_start(exon_list,curAl->core.pos,read_start,i);
    if(!ret_val){
        std::cerr<<"@ERROR::Can not get the genomic read start"<<std::endl;
        exit(1);
    }

    // convert the mate information for single reads without a concordantly mapped mate
    // unless the read is not paired, in which case the paired information should stay the same
    if(curAl->core.mpos != -1 || curAl->core.mtid != -1){
        int target_name_mate = atoi(al_hdr->target_name[curAl->core.mtid]); // name of the transcript from the input alignment
        GffTranscript& p_trans_mate = transcriptome[target_name_mate]; // get the transcript
        GVec<GSeg>& exon_list_mate=p_trans_mate.exons; // get exons
        int i_mate=0;
        int32_t read_start_mate=0;

        ret_val = this->get_read_start(exon_list_mate,curAl->core.mpos,read_start_mate,i_mate);
        if(!ret_val){
            std::cerr<<"@ERROR::Can not get the genomic read start of the mate"<<std::endl;
            exit(1);
        }

        curAl->core.mtid = p_trans_mate.refID;
        curAl->core.mpos = read_start_mate - 1;
    }

    cur_pos.set_chr(p_trans.get_refID());
    cur_pos.set_start(read_start);
    cur_pos.set_locus(p_trans.get_geneID());
    cur_pos.set_strand(p_trans.strand);
    cur_pos.set_trans(target_name);

    // secondly build a new cigar string
    int cur_cigar_full[MAX_CIGARS];
    memcpy(cur_cigar_full, bam_get_cigar(curAl), curAl->core.n_cigar);

    ret_val = this->convert_cigar(i,next_exon,exon_list,num_cigars,read_start,curAl,cigars,cur_pos);
    if (!ret_val) {
        std::cerr << "@ERROR::Can not create a new cigar string for the single read from process_read" << std::endl;
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

    // assign new cigar string to the record

    return cigar2hash(cigars,num_cigars);
}

bool Comparator::get_read_start(GVec<GSeg>& exon_list,int32_t gff_start,int32_t& genome_start, int& exon_idx){
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

int Comparator::process_single(bam1_t *curAl,bam_hdr_t* al_hdr,int plus_nm){
    Position cur_pos;
    int cigars[MAX_CIGARS];
    int num_cigars=0;
    size_t cigar_hash = this->process_read(curAl,cur_pos,cigars,num_cigars,al_hdr);

    int cur_num_multi = this->evaluate_multimappers(curAl,cur_pos,cigars,num_cigars,plus_nm); // also finishes the read
    if(cur_num_multi>0){
        this->num_multi++;
        this->num_multi_hits=this->num_multi_hits+cur_num_multi;
    }
    total_num_al++;

    return cur_num_multi;
}

int  Comparator::process_single_mate(bam1_t *aligned, bam1_t *unaligned,bam_hdr_t* al_hdr,int plus_nm){
    int cur_num_multi = process_single(aligned,al_hdr);
    aligned->core.mpos = 0;
    aligned->core.mtid = -1;
    finish_read(aligned);
    unaligned->core.pos = 0;
    unaligned->core.tid = -1;
    unaligned->core.mpos = aligned->core.pos;
    unaligned->core.mtid = aligned->core.tid;
    finish_read(unaligned);
    return cur_num_multi;
}

void Comparator::add_cigar(bam1_t *curAl,int num_cigars,int* cigars){
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

int Comparator::evaluate_multimappers_pair(bam1_t *curAl,bam1_t* curAl_mate,Position &cur_pos,Position &cur_pos_mate,
                                          int *cigars,int *cigars_mate,int &num_cigars,int &num_cigars_mate,int plus_nm) {
    int unique;
    std::vector<Position> res_pos,res_pos_mate; // holds the results of the multimapper evaluation
    if(!this->abund){ // compute abundance dynamically
//        std::cout<<"\n================\n"<<std::endl;
//        std::cout<<bam_get_qname(curAl)<<std::endl;
        unique = this->mmap.process_pos_pair(cur_pos,cur_pos_mate,this->loci,res_pos,res_pos_mate);
    }
    else{ // rely on the salmon abundance
        unique = this->mmap.process_pos_pair_precomp(cur_pos,cur_pos_mate,this->loci,res_pos,res_pos_mate);
    }
//    std::cerr<<"eval multi_pair: "<<unique<<std::endl;
    if(res_pos.empty()){ // increment abundance
        if(!this->abund) { // TODO: only needed when the not reporting all multimappers but rather performing likelihood assignment
            this->loci.add_read(cur_pos.locus);
            this->loci.add_read(cur_pos_mate.locus);
        }
        curAl->core.pos = cur_pos.start-1;
        curAl->core.mpos = cur_pos_mate.start-1;
        curAl_mate->core.pos = cur_pos_mate.start-1;
        curAl_mate->core.mpos = cur_pos.start-1;
        set_xs(curAl,cur_pos.strand);
        set_xs(curAl_mate,cur_pos_mate.strand);
        // add already computed cigar to the read
        add_cigar(curAl, num_cigars, cigars); // will be performed afterwards
        add_cigar(curAl_mate, num_cigars_mate, cigars_mate); // will be performed afterwards
        change_nh_flag(curAl,plus_nm+1);
        change_nh_flag(curAl_mate,plus_nm+1);
        this->finish_read(curAl);
        this->finish_read(curAl_mate);
        return unique;

        // looking specifically at the cases where agar didn't get it right but hisat2 did and the wrong alignment was preferred
        // What else can we do then?
        // we can look at the reads that come from the known transcriptome which agar get's wrong (of which there should be virtually none)
        // to do so, we may want to look at the multi_dist test - this should hopefully give us information regarding the

        // what if we use genome fasta generated by salmontools to build hisat index???
    }
    else{
        change_nh_flag(curAl,res_pos.size()+plus_nm);
        change_nh_flag(curAl_mate,res_pos.size()+plus_nm);
        bool prim = true;
        for(int pos_idx=0;pos_idx<res_pos.size();pos_idx++){
            // increment total abundances of the locus to which the new cur_pos belongs
            if(!this->abund) {
                this->loci.add_read_multi(res_pos[pos_idx].locus);
                this->loci.add_read_multi(res_pos_mate[pos_idx].locus);
            }

            if(std::strcmp(bam_get_qname(curAl),"read4234734/rna-XM_006723946.2;mate1:1923-2023;mate2:2322-2422")==0){
                print_cigar(curAl);
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

            add_multi_tag(curAl,res_pos[pos_idx].transID);
            add_multi_tag(curAl_mate,res_pos_mate[pos_idx].transID);

            curAl->core.pos = res_pos[pos_idx].start-1;
            curAl->core.mpos = res_pos_mate[pos_idx].start-1;
            curAl->core.tid = res_pos[pos_idx].chr;
            curAl->core.mtid = res_pos_mate[pos_idx].chr;
            curAl_mate->core.pos = res_pos_mate[pos_idx].start-1;
            curAl_mate->core.mpos = res_pos[pos_idx].start-1;
            curAl_mate->core.tid = res_pos_mate[pos_idx].chr;
            curAl_mate->core.mtid = res_pos[pos_idx].chr;
            set_xs(curAl,res_pos[pos_idx].strand);
            set_xs(curAl_mate,res_pos_mate[pos_idx].strand);

            // reconvert the cur_pos into a read and output
            num_cigars = 0,num_cigars_mate = 0;
            memset(cigars, 0, MAX_CIGARS);
            memset(cigars_mate, 0, MAX_CIGARS);

            int ret_val = this->convert_cigar(i,next_exon,exon_list,num_cigars,read_start,curAl,cigars,res_pos[pos_idx]);
            if (!ret_val) {
                std::cerr << "@ERROR::Can not create a new cigar string for the single read evaluate_multimapper_pair1" << std::endl;
                exit(1);
            }
            add_cigar(curAl, num_cigars, cigars); // will be performed afterwards
            if(prim){ // set as primary alignment
                curAl->core.flag &= ~BAM_FSECONDARY;
            }
            else{ // set as secondary alignment
                curAl->core.flag |= BAM_FSECONDARY;
            }
            this->finish_read(curAl);

            if(std::strcmp(bam_get_qname(curAl),"read4234734/rna-XM_006723946.2;mate1:1923-2023;mate2:2322-2422")==0){
                print_cigar(curAl);
            }

            int ret_val_mate = this->convert_cigar(i_mate,next_exon_mate,exon_list_mate,num_cigars_mate,read_start_mate,curAl_mate,cigars_mate,res_pos_mate[pos_idx]);
            if (!ret_val_mate) {
                std::cerr << "@ERROR::Can not create a new cigar string for the single read evaluate_multimapper_pair2" << std::endl;
                exit(1);
            }
            add_cigar(curAl_mate, num_cigars_mate, cigars_mate); // will be performed afterwards
            if(prim){ // set as primary alignment
                curAl_mate->core.flag &= ~BAM_FSECONDARY;
            }
            else{ // set as secondary alignment
                curAl_mate->core.flag |= BAM_FSECONDARY;
            }
            this->finish_read(curAl_mate);

            if(std::strcmp(bam_get_qname(curAl_mate),"read4234734/rna-XM_006723946.2;mate1:1923-2023;mate2:2322-2422")==0){
                print_cigar(curAl_mate);
            }
            prim = false;
        }
        return res_pos.size();
    }
}

int Comparator::evaluate_multimappers(bam1_t* curAl,Position& cur_pos,int cigars[MAX_CIGARS],int &num_cigars,int plus_nm=0){
    int unique;
    std::vector<Position> res_pos; // holds the results of the multimapper evaluation
    if(!this->abund){ // compute abundance dynamically
        unique = this->mmap.process_pos(cur_pos,this->loci,res_pos);
    }
    else{ // compute abundance dynamically
        unique = this->mmap.process_pos_precomp(cur_pos,this->loci,res_pos);
    }
    if(res_pos.empty()){ // increment abundance
        if(!this->abund){
            this->loci.add_read(cur_pos.locus);
        }
        curAl->core.pos = cur_pos.start-1;
        curAl->core.tid = cur_pos.chr;
        set_xs(curAl,cur_pos.strand);
        // add already computed cigar to the read
        change_nh_flag(curAl,plus_nm+1);
        add_cigar(curAl, num_cigars, cigars); // will be performed afterwards
        this->finish_read(curAl);
        return unique;
    }
    else{
        change_nh_flag(curAl,res_pos.size()+plus_nm);
        bool prim = true;
        for(auto &v : res_pos){
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
            set_xs(curAl,v.strand);

            add_multi_tag(curAl,v.transID);

            // reconvert the cur_pos into a read and output
            num_cigars = 0;
            memset(cigars, 0, MAX_CIGARS);

            int ret_val = this->convert_cigar(i,next_exon,exon_list,num_cigars,read_start,curAl,cigars,v);
            if (!ret_val) {
                std::cerr << "@ERROR::Can not create a new cigar string for the single read evaluate_multimappers" << std::endl;
                exit(1);
            }

            add_cigar(curAl, num_cigars, cigars); // will be performed afterwards
            if(prim){ // set as primary alignment
                curAl->core.flag &= ~BAM_FSECONDARY;
            }
            else{ // set as secondary alignment
                curAl->core.flag |= BAM_FSECONDARY;
            }
            this->finish_read(curAl);
            prim = false;
        }
        return res_pos.size();
    }
}

int Comparator::convert_cigar(int i,GSeg *next_exon,GVec<GSeg>& exon_list,int &num_cigars,int read_start,bam1_t* curAl,int cigars[MAX_CIGARS],Position& pos_obj){
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

int Comparator::process_status(bam1_t *mate,bam1_t *curAl,bool &check_mate,bam_hdr_t* al_hdr,int plus_nm){
    int cur_num_multi;
    if(std::strcmp(bam_get_qname(curAl),bam_get_qname(mate)) == 0){ // is a pair
        uint8_t status = get_pair_status(curAl,mate);
        // do stuff with the pair based on the status
        switch(status){
            case STATUS::CONC:
                cur_num_multi = this->process_pair(curAl,mate,al_hdr,plus_nm);
                break;
            case STATUS::DISC:
                cur_num_multi = this->process_disc(curAl,mate,al_hdr,plus_nm);
                break;
            case STATUS::ALIGNED_MIXED_1:
                cur_num_multi = this->process_single_mate(mate,curAl,al_hdr,plus_nm);
                break;
            case STATUS::ALIGNED_MIXED_2:
                cur_num_multi = this->process_single_mate(curAl,mate,al_hdr,plus_nm);
                break;
            case STATUS::UNALIGNED_PAIR:
                std::cerr<<"@ERROR:::unaligned transcriptomic not supported now: "<<bam_get_qname(curAl)<<std::endl;
                exit(-1);
//                if (!this->unaligned_mode){finish_read(curAl);finish_read(mate);break;}
//                write_unaligned_pair(curAl,mate);
//                this->num_unal_paired++;
//                break;
            default:
                std::cout<<"@ERROR::SWITCH PAIR"<<std::endl;
                break;
        }
        check_mate = false;
    }
    else{ // is not a pair - process the first as a singleton and move the second to the mate
        uint8_t status = get_read_status(mate);
        // do stuff with the single read based on the status
        switch(status){
            case STATUS::ALIGNED_SINGLE:
                cur_num_multi = this->process_single(mate,al_hdr,plus_nm);
                break;
            case STATUS::UNALIGNED_SINGLE:
                std::cerr<<"@ERROR:::unaligned transcriptomic not supported now for single either"<<std::endl;
                exit(-1);
//                if (!this->unaligned_mode){finish_read(mate);break;}
//                write_unaligned(mate,this->unal_s);
//                this->num_unal_single++;
//                break;
            default:
                std::cout<<"@ERROR::SWITCH SINGLE"<<std::endl;
                break;
        }
        bam_copy1(mate,curAl);
        check_mate = true;
    }
    return cur_num_multi;
}

void Comparator::write_remaining_gt_read(bam1_t* al_g1,bam1_t* al_t1,bam1_t* al_g2,bam1_t* al_t2,int& ret_g,int& ret_t){
    // create set of coordinates
    // add genomic
    // add transcriptomic and multimappers
    // write all
}

void Comparator::merge(){
    std::cout<<"merging alignments"<<std::endl;
    // read genome and transcriptome files simultaneously - processing reads one group at a time (assumes the same ordering)
    bam1_t* al_t1 = bam_init1();
    bam1_t* al_t2 = bam_init1();
    bam1_t* al_g1 = bam_init1();
    bam1_t* al_g2 = bam_init1();

    std::vector<std::pair<bam1_t*,bam1_t*>> al_g_conc,al_g_disc;
    std::vector<bam1_t*> al_g_single = {};

    std::string readname_genome = "",readname_trans = "";

    int ret_t,ret_g;
    int gt_cmp_status;

    ret_g = sam_read1(genome_al,al_g_hdr,al_g1);
    ret_t = sam_read1(trans_al,al_t_hdr,al_t1);
    if(ret_g<0 || ret_t<0){
        bam_destroy1(al_t1);
        bam_destroy1(al_t2);
        bam_destroy1(al_g1);
        bam_destroy1(al_g2);
        return;
    }
    readname_genome = bam_get_qname(al_g1);
    readname_trans = bam_get_qname(al_t1);

    // TODO: process in pairs to compute more accurate nm estimate and perform better cmp_gt

    // genome on the top level, and transcriptome on the bottom level - get read from transcriptome - ASSUMES reads are in the same order
    while(true) { // only perfom if unaligned flag is set to true
        al_g_single.clear();
        // see if the mate needs to be loaded
//        if(al_g1->core.flag &0x2){ // read is paired and aligned as a pair - load mate
            ret_g = sam_read1(genome_al,al_g_hdr,al_g2);
            // make sure both mates are currently present
//            if(al_g1->core.pos != al_g2->core.mpos || al_g1->core.tid != al_g2->core.tid){
//                std::cerr<<"@ERROR:::detected concordant does not have valid matches"<<std::endl;
//                exit(-1);
//            }
            // load next transcriptomic read
            ret_t = sam_read1(trans_al,al_t_hdr,al_t2);
            if(std::strcmp(bam_get_qname(al_t1),bam_get_qname(al_t2))!=0){
                std::cerr<<"@ERROR:::transcriptomic read does not have expected mate"<<std::endl;
                exit(-1);
            }

            gt_cmp_status = cmp_gt(al_g1,al_g2,al_t1,al_t2);

            // PROCESS separately: concordant and discordant/single - discordant and single should be compared on a per-mate bases while concordant should be compared based on both reads
            if(gt_cmp_status == STATUS_CMP_READ::UNALIGNED){
                // skip to the next read in transcriptome
                ret_t = skip_read(trans_al,al_t_hdr,al_t1);

                // write remaining genomic reads from the same read group
                ret_g = write_remaining_genome_pair(genome_al,al_g_hdr,al_g1,al_g2);
            }
            else if(gt_cmp_status == STATUS_CMP_READ::FIRST){ // just output hisat2 alignments
                // skip to the next read in transcriptome
                ret_t = skip_read(trans_al,al_t_hdr,al_t1);

                // write remaining genomic reads from the same read group
                ret_g = write_remaining_genome_pair(genome_al,al_g_hdr,al_g1,al_g2);
            }
            else if(gt_cmp_status == STATUS_CMP_READ::SECOND){ // process transcriptomic and output TODO: transcriptomic can be discordant while in fact being concordant - potential solution is to simply ignore and process multimappers for hisat alignments
                // skip to the next read in the genome
                ret_g = skip_read(genome_al,al_g_hdr,al_g1);

                ret_t = write_remaining_trans_pair(trans_al,al_t_hdr,al_t1,al_t2,false);
            }
            else if(gt_cmp_status == STATUS_CMP_READ::BOTH){ // process both genomic and transcriptomic
                ret_g = skip_read(genome_al,al_g_hdr,al_g1);
//                ret_t = skip_read(trans_al,al_t_hdr,al_t1);

//                ret_g = write_remaining_genome_pair(genome_al,al_g_hdr,al_g1,al_g2);
                ret_t = write_remaining_trans_pair(trans_al,al_t_hdr,al_t1,al_t2,true);
//                ret_g = load_remaining_genome_pair(genome_al,al_g_hdr,al_g1,al_g2,al_g_single);
//                ret_t = write_remaining_both_pair(trans_al,al_t_hdr,al_t1,al_t2,al_g_single);

//            write_remaining_gt_read(al_g1,al_t1,al_g2,al_t2,ret_g,ret_t);
            }
            else if(gt_cmp_status == STATUS_CMP_READ::ERR_CMP){
                std::cerr<<"@ERROR:::should not be error in read comparison"<<std::endl;
            }
            else{
                std::cerr<<"@ERROR:::unrecognized error in read comparison: "<<bam_get_qname(al_g1)<<"\t"<<bam_get_qname(al_t1)<<std::endl;
            }

            // by now we should have advanced to the next read
            if(ret_g<0 && ret_t<0){
                break;
            }
            else if(ret_g<0 || ret_t<0){
                std::cerr<<"@ERROR:::Only one file is over"<<std::endl;
            }
            else{ // continue
                readname_genome = bam_get_qname(al_g1);
                readname_trans = bam_get_qname(al_t1);
            }
//        }
//        else{
//            gt_cmp_status = cmp_gt(al_g1,al_t1);
//
//            // PROCESS separately: concordant and discordant/single - discordant and single should be compared on a per-mate bases while concordant should be compared based on both reads
//            if(gt_cmp_status == STATUS_CMP_READ::UNALIGNED){
//                // skip to the next read in transcriptome
//                ret_t = skip_read(trans_al,al_t_hdr,al_t1);
//
//                // write remaining genomic reads from the same read group
//                ret_g = write_remaining_genome_read(genome_al,al_g_hdr,al_g1);
//            }
//            else if(gt_cmp_status == STATUS_CMP_READ::FIRST){ // just output hisat2 alignments
//                // skip to the next read in transcriptome
//                ret_t = skip_read(trans_al,al_t_hdr,al_t1);
//
//                // write remaining genomic reads from the same read group
//                ret_g = write_remaining_genome_read(genome_al,al_g_hdr,al_g1);
//            }
//            else if(gt_cmp_status == STATUS_CMP_READ::SECOND){ // process transcriptomic and output TODO: transcriptomic can be discordant while in fact being concordant - potential solution is to simply ignore and process multimappers for hisat alignments
//                // skip to the next read in the genome
//                ret_g = skip_read(genome_al,al_g_hdr,al_g1);
//
//                ret_t = write_remaining_trans_read(trans_al,al_t_hdr,al_t1,al_t2,false);
//            }
//            else if(gt_cmp_status == STATUS_CMP_READ::BOTH){ // process both genomic and transcriptomic
//            ret_g = skip_read(genome_al,al_g_hdr,al_g1);
//
////                ret_g = write_remaining_genome_read(genome_al,al_g_hdr,al_g1);
//                ret_t = write_remaining_trans_read(trans_al,al_t_hdr,al_t1,al_t2,true);
//
////            write_remaining_gt_read(al_g1,al_t1,al_g2,al_t2,ret_g,ret_t);
//            }
//            else if(gt_cmp_status == STATUS_CMP_READ::ERR_CMP){
//                std::cerr<<"@ERROR:::should not be error in read comparison"<<std::endl;
//            }
//            else{
//                std::cerr<<"@ERROR:::unrecognized error in read comparison: "<<bam_get_qname(al_g1)<<"\t"<<bam_get_qname(al_t1)<<std::endl;
//            }
//
//            // by now we should have advanced to the next read
//            if(ret_g<0 && ret_t<0){
//                break;
//            }
//            else if(ret_g<0 || ret_t<0){
//                std::cerr<<"@ERROR:::Only one file is over"<<std::endl;
//            }
//            else{ // continue
//                readname_genome = bam_get_qname(al_g1);
//                readname_trans = bam_get_qname(al_t1);
//            }
//        }
    }
    // TODO: consider prefering non-spliced transcriptomic to spliced genomic alignments potentially allowing a few additional mismatches

    bam_destroy1(al_t1);
    bam_destroy1(al_t2);
    bam_destroy1(al_g1);
    bam_destroy1(al_g2);
    std::cout<<"done merging alignments"<<std::endl;
}

void Comparator::set_fraglen(int fraglen) {
    this->fraglen = fraglen;
    if(!this->multi){
        return; // no need to set the argument if multimappers are not being evaluated
    }
    this->mmap.set_fraglen(fraglen);
}

void Comparator::load_info(const std::string& info_fname){
    // read file to get important stats
    struct stat buffer{};
    if(stat (info_fname.c_str(), &buffer) != 0){ // if file does not exists
        std::cerr<<"@ERROR::Info file is not found: "<<info_fname<<std::endl;
        exit(1);
    }
    // read file and save important info
    std::cerr<<"@LOG::Reading the info file: "<<info_fname<<std::endl;
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
    std::cerr<<"@LOG::Loaded info data"<<std::endl;
}

// this function parses the .tlst file and inserts entries into the index
void Comparator::_load_transcriptome(std::ifstream& tlstfp,char* buffer){
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
                        std::cerr<<"@ERROR::should never happen _load_transcriptome"<<std::endl;
                        exit(1);
                }
                break;
            default:
                std::cerr<<"@ERROR::unrecognized character: "<<buffer[i]<<std::endl;
                exit(1);
        }
    }
    if(!transcript.empty){
        transcript.add_exon(exon);
        this->transcriptome[transcript.gffID] = transcript;
    }
}

void Comparator::print_transcriptome(){
    std::cerr<<"@LOG::printing transcriptome"<<std::endl;
    for(auto &t : this->transcriptome){
        t.print();
    }
    std::cerr<<"@LOG::done printing transcriptome"<<std::endl;
}

void Comparator::load_transcriptome(const std::string &tlst_fname) {
    // first initiate the transcriptome array to hold the transcripts
    this->transcriptome = std::vector<GffTranscript>(this->numTranscripts+1);

    struct stat buffer{};
    if(stat (tlst_fname.c_str(), &buffer) != 0){ // if file does not exists
        std::cerr<<"@ERROR::Locus file: "<<tlst_fname<<" is not found. Check that correct index is provided."<<std::endl;
        exit(1);
    }
    std::cerr<<"@LOG::loading locus data from: "<<tlst_fname<<std::endl;
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
        std::cerr<<"@LOG::finished loading the locus data"<<std::endl;
    }
    else{
        std::cerr<<"@ERROR::failed to open the locus file"<<std::endl;
    }
    tlstfp.close();
}

void Comparator::load_genome_header(const std::string& genome_header_fname){
    struct stat buffer{};
    if(stat (genome_header_fname.c_str(), &buffer) != 0){ // if file does not exists
        std::cerr<<"@ERROR::Genome Header file: "<<genome_header_fname<<" is not found. Check that correct index is provided."<<std::endl;
        exit(1);
    }
    this->idx_al=hts_open(genome_header_fname.c_str(),"r");
}

void Comparator::load_index(const std::string& index_base,bool multi){
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
void Comparator::load_abundances(const std::string& abundFP){
    std::cout<<"using salmon abundances"<<std::endl;
    this->abund = true;
    this->loci.set_precomputed_abundances();
    std::ifstream abundstream;
    std::stringstream *linestream;
    abundstream.open(abundFP.c_str(),std::ios::in);
    if (!abundstream.good()){
        std::cerr<<"@ERROR::Couldn't open transcript abundance file: "<<abundFP<<std::endl;
        exit(1);
    }
    std::ios::sync_with_stdio(false);

    std::cerr<<"@LOG::Reading the transcript abundance file: "<<abundFP<<std::endl;
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
    std::cerr<<"@LOG::Loaded transcript abundance data"<<std::endl;
}

void Comparator::print_abundances(const std::string& out_abund_fname){
    std::ofstream out_abund_fp;
    out_abund_fp.open(out_abund_fname, std::ofstream::out | std::ofstream::trunc);
    if(this->loci.is_precomp_abund()){
        for(int i=0;i<this->loci.get_size();i++){
            out_abund_fp<<i<<"\t"<<this->loci[i].get_abund()<<std::endl;
        }
    }
    else{
        for(int i=0;i<this->loci.get_size();i++){
            out_abund_fp<<i<<"\t"<<this->loci.get_abund(i)<<std::endl;
        }
    }

    out_abund_fp.close();
}

void Comparator::load_multi(const std::string& multiFP){
    this->multi = true;
    this->mmap.load(multiFP);
}

void Comparator::print_multimappers() {
    this->mmap.print_multimapers();
}

void Comparator::print_stats(){
    std::cerr<<"STATS"<<std::endl;
}

void Comparator::finish_read(bam1_t *curAl){
    int ret_val = sam_write1(this->merged_al, merged_al_hdr, curAl);
}