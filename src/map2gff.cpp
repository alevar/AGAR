#include "map2gff.h"
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

Map2GFF::Map2GFF(const std::string& gffFP, const std::string& alFP){
    tlststream.open(gffFP.c_str(),std::ios::in);
    if (!tlststream.good()){
        std::cerr<<"FATAL: Couldn't open trascript data: "<<gffFP<<std::endl;
        exit(1);
    }
    std::ios::sync_with_stdio(false);
    
    std::cout<<"Reading the transcript data: "<<gffFP<<std::endl;
    std::string tline;
    while (std::getline(tlststream,tline)){
        if (tline.length()>4){
            GffTranscript* t=new GffTranscript(tline);
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

int Map2GFF::convert_cigar(int i,int cur_intron_len,int miss_length,GSeg *next_exon,int match_length,GVec<GSeg>& exon_list,
                           int &num_cigars,int read_start,bam1_t* curAl,std::string &cigar_str,int cigars[MAX_CIGARS]){
    int cur_pos=read_start;
    for (int c=0;c<curAl->core.n_cigar;++c){
        uint32_t *cigar_full=bam_get_cigar(curAl);
        int opcode=bam_cigar_op(cigar_full[c]);
        int length=bam_cigar_oplen(cigar_full[c]);
        
        if (opcode==BAM_CINS){
            cigars[num_cigars]=opcode|(length<<BAM_CIGAR_SHIFT);
            cigar_str+=std::to_string(bam_cigar_oplen(cigars[num_cigars]));
            cigar_str+=bam_cigar_opchr(opcode);
            ++num_cigars;
        }
        if (opcode==BAM_CSOFT_CLIP){
            cigars[num_cigars]=opcode|(length<<BAM_CIGAR_SHIFT);
            cigar_str+=std::to_string(bam_cigar_oplen(cigars[num_cigars]));
            cigar_str+=bam_cigar_opchr(opcode);
            ++num_cigars;
        }
        if (opcode != BAM_CMATCH && opcode != BAM_CDEL){
            continue;
        }
        int remaining_length=length;
        for (;i<exon_list.Count();++i){
            GSeg& cur_exon=exon_list[i];
            if (cur_pos>=(int)cur_exon.start && cur_pos+remaining_length-1<=(int)cur_exon.end){
                cigars[num_cigars]=opcode | (remaining_length <<BAM_CIGAR_SHIFT);
                cigar_str+=std::to_string(bam_cigar_oplen(cigars[num_cigars]));
                cigar_str+=bam_cigar_opchr(opcode);
                ++num_cigars;
                cur_pos+=remaining_length;
                break;
            }
            else if (cur_pos >= (int)cur_exon.start && cur_pos+remaining_length-1>(int)cur_exon.end){
                match_length=(int)cur_exon.end-cur_pos+1;
                if (match_length>0){
                    cigars[num_cigars]=opcode | (match_length << BAM_CIGAR_SHIFT);
                    cigar_str+=std::to_string(bam_cigar_oplen(cigars[num_cigars]));
                    cigar_str+=bam_cigar_opchr(opcode);
                    ++num_cigars;
                }
                if (i+1>=exon_list.Count()){
                    // std::cout<<"hola"<<std::endl;
                    return 0;
                }
                else{
                    next_exon=&(exon_list[i+1]);
                }
                miss_length=next_exon->start-cur_exon.end-1;
                cur_intron_len+=miss_length;
                
                cigars[num_cigars]=BAM_CREF_SKIP | (miss_length <<BAM_CIGAR_SHIFT);
                cigar_str+=std::to_string(bam_cigar_oplen(cigars[num_cigars]));
                cigar_str+=bam_cigar_opchr(bam_cigar_op(cigars[num_cigars]));
                ++num_cigars;

                cur_pos+=match_length+miss_length;
                remaining_length-=match_length;
                assert(cur_pos == (int)next_exon->start);
            }
        }
    }
    return 1;
}

void Map2GFF::convert_coords(const std::string& outFP, const std::string& genome_header){
    genome_al=hts_open(genome_header.c_str(),"r");
    genome_al_hdr=sam_hdr_read(genome_al);

    samFile *outSAM=sam_open(outFP.c_str(),"wb");
    bam_hdr_t *outSAM_header=bam_hdr_init();
    
    outSAM_header=bam_hdr_dup(genome_al_hdr);
    int res=sam_hdr_write(outSAM,outSAM_header);
    for (int i=0;i<genome_al_hdr->n_targets;++i){
        ref_to_id[genome_al_hdr->target_name[i]]=i;
    }
    bam_hdr_destroy(outSAM_header);

    GffTranscript *p_trans=NULL;
    GffTranscript *mate_p_trans=NULL;

    bool success;

    std::unordered_map<std::string,MateRead*> curReadGroup_tmp;
    std::unordered_map<std::string,MatePair*> curReadGroup_paired;
    std::unordered_map<std::string,bam1_t*> curReadGroup;
    std::string curReadName="";

    bool start=true;

    while(sam_read1(al,al_hdr,curAl)>0){
        // std::cout<<bam_get_qname(curAl)<<std::endl;
        std::string newReadName=bam_get_qname(curAl);
        // std::cout<<newReadName<<std::endl;
        if (newReadName.compare(curReadName)==0 || start){
            // if (newReadName.compare("read21068/1304")==0){
            //     std::cout<<"HELLO READ"<<std::endl;
            // }
            curReadName=newReadName;
            start=false;
            success=true;
            size_t read_start=0;
            size_t mate_read_start=0;

            const char *target_name=al_hdr->target_name[curAl->core.tid];
            int trans_idx=atoi(target_name);
            p_trans=tidx_to_t[target_name];

            GVec<GSeg>& exon_list=p_trans->exons;

            GSeg *next_exon=NULL;
            int cigars[MAX_CIGARS];

            int cur_pos,match_length,miss_length,cur_intron_len=0,i=0,num_cigars=0;

            bool ret_val=Map2GFF::get_read_start(exon_list,curAl->core.pos,read_start,i);
            if (!ret_val){
                std::cerr<<"SOMETHING WRONG WITH GETTING GENOMIC READ START"<<std::endl;
            }

            std::string cigar_str="";
            int cigar_ret_val=Map2GFF::convert_cigar(i,cur_intron_len,miss_length,next_exon,match_length,exon_list,num_cigars,read_start,curAl,cigar_str,cigars);
            if (cigar_ret_val){
                bool paired=curAl->core.flag           &    1;
                bool pairedly_aligned=curAl->core.flag &    2;
                // evaluate as a pair if a singleton is encountered - forget about it
                int mate_i=0;
                const char *mate_target_name=al_hdr->target_name[curAl->core.mtid];
                int mate_trans_idx=atoi(mate_target_name);
                mate_p_trans=tidx_to_t[mate_target_name];
                bool mate_ret_val=Map2GFF::get_read_start(exon_list,curAl->core.mpos,mate_read_start,mate_i);
                // if (!mate_ret_val){
                //     std::cerr<<"SOMETHING WRONG WITH GETTING GENOMIC MATE READ START"<<std::endl;
                //     // need aditional cases here
                // }

                int origRefID=curAl->core.tid, origPos=curAl->core.pos, origMPos=curAl->core.mpos;
                curAl->core.tid=ref_to_id[p_trans->refID];
                curAl->core.pos=read_start-1;
                if (mate_ret_val){
                    curAl->core.mtid=ref_to_id[mate_p_trans->refID];
                    curAl->core.mpos=mate_read_start-1;
                }
                else{
                    curAl->core.mtid=ref_to_id[p_trans->refID];
                    curAl->core.mpos=read_start-1;
                }

                int old_n_cigar=curAl->core.n_cigar;
                if (num_cigars != old_n_cigar){
                    int data_len=curAl->l_data+4*(num_cigars-old_n_cigar);
                    int m_data=std::max(data_len,(int)curAl->m_data);
                    kroundup32(m_data);

                    uint8_t* data = (uint8_t*)calloc(m_data,1);

                    int copy1_len = (uint8_t*)bam_get_cigar(curAl) - curAl->data;
                    memcpy(data, curAl->data, copy1_len);

                    int copy2_len = num_cigars * 4;
                    memcpy(data + copy1_len, cigars, copy2_len);

                    int copy3_len = curAl->l_data - copy1_len - (old_n_cigar * 4);
                    memcpy(data + copy1_len + copy2_len, bam_get_seq(curAl), copy3_len);

                    curAl->core.n_cigar = num_cigars;

                    free(curAl->data);
                    curAl->data = data;
                    curAl->l_data = data_len;
                    curAl->m_data = m_data;
                }
                
                char strand=p_trans->strand;
                uint8_t* ptr=bam_aux_get(curAl,"XS");
                if(ptr){
                    bam_aux_del(curAl,ptr);
                }
                if (strand=='-'){
                    bam_aux_append(curAl,"XS",'A',1,(const unsigned char*)"-");
                }
                if (strand=='+'){
                    bam_aux_append(curAl,"XS",'A',1,(const unsigned char*)"+");
                }

                // add optional tag with transcript ID
                uint8_t* ptr_op=bam_aux_get(curAl,"OP");
                if(ptr_op){
                    bam_aux_del(curAl,ptr_op);
                }
                bam_aux_append(curAl,"OP",'i',4,(uint8_t*)&origRefID);

                if (mate_ret_val && paired && pairedly_aligned){
                    if (curAl->core.flag&64){
                        std::string tmpKey_pre=std::to_string(origRefID)+"$"+std::to_string(origPos)+"_"+std::to_string(origMPos);
                        if (curReadGroup_tmp.find(tmpKey_pre)==curReadGroup_tmp.end()){
                            curReadGroup_tmp[tmpKey_pre]=new MateRead(bam_dup1(curAl),read_start,cigar_str);
                        }
                        else{
                            MateRead *mr=curReadGroup_tmp[tmpKey_pre];
                            if (mr->al->core.flag&128){
                                std::string tmpKey=p_trans->refID+"$"+std::to_string(read_start)+"_"+cigar_str+"&"+std::to_string(mr->pos)+"_"+mr->cigar;
                                if (curReadGroup_paired.find(tmpKey)==curReadGroup_paired.end()){
                                    curReadGroup_paired[tmpKey]=new MatePair(bam_dup1(curAl),bam_dup1(mr->al));
                                }
                                delete curReadGroup_tmp[tmpKey_pre];
                                curReadGroup_tmp.erase(tmpKey_pre);
                            }
                            else{
                                std::cout<<"BIG MISTAKE AND DON'T KNOW WHAT TO DO"<<std::endl;
                            }
                        }
                    }
                    else if (curAl->core.flag&128){
                        std::string tmpKey_pre=std::to_string(origRefID)+"$"+std::to_string(origMPos)+"_"+std::to_string(origPos);
                        if (curReadGroup_tmp.find(tmpKey_pre)==curReadGroup_tmp.end()){
                            curReadGroup_tmp[tmpKey_pre]=new MateRead(bam_dup1(curAl),read_start,cigar_str);
                        }
                        else{
                            // now write data where need be
                            MateRead *mr=curReadGroup_tmp[tmpKey_pre];
                            if (mr->al->core.flag&64){
                                std::string tmpKey=p_trans->refID+"$"+std::to_string(mr->pos)+"_"+mr->cigar+"&"+std::to_string(read_start)+"_"+cigar_str;
                                if (curReadGroup_paired.find(tmpKey)==curReadGroup_paired.end()){
                                    curReadGroup_paired[tmpKey]=new MatePair(bam_dup1(mr->al),bam_dup1(curAl));
                                }
                                delete curReadGroup_tmp[tmpKey_pre];
                                curReadGroup_tmp.erase(tmpKey_pre);
                            }
                            else{
                                std::cout<<"BIG MISTAKE #2 AND DON'T KNOW WHAT TO DO"<<std::endl;
                            }
                        }
                    }
                    else{
                        std::cout<<"mate_error"<<std::endl;
                        exit(-1);
                    }
                }
                else{
                    std::string tmpKey=p_trans->refID+"$"+std::to_string(read_start)+"_"+cigar_str;
                    curReadGroup[tmpKey]=bam_dup1(curAl);
                }
            }
        }
        else{
            if (!curReadGroup_paired.empty()){
                if (curReadGroup_tmp.empty()){ // Make sure nothing is left in the tmp group
                    int nh=curReadGroup_paired.size();
                    for (std::pair<std::string,MatePair*> kv: curReadGroup_paired){
                        // Process first mate
                        uint8_t* ptr_nh_1=bam_aux_get(kv.second->firstMate,"NH");
                        if(ptr_nh_1){
                            bam_aux_del(kv.second->firstMate,ptr_nh_1);
                        }
                        bam_aux_append(kv.second->firstMate,"NH",'i',4,(uint8_t*)&nh);

                        kv.second->firstMate->core.bin = hts_reg2bin(kv.second->firstMate->core.pos, kv.second->firstMate->core.pos + bam_cigar2rlen(kv.second->firstMate->core.n_cigar, bam_get_cigar(kv.second->firstMate)), 14, 5);
                        int ret_sam=sam_write1(outSAM,genome_al_hdr,kv.second->firstMate);

                        // Now process second mate
                        uint8_t* ptr_nh_2=bam_aux_get(kv.second->secondMate,"NH");
                        if(ptr_nh_2){
                            bam_aux_del(kv.second->secondMate,ptr_nh_2);
                        }
                        bam_aux_append(kv.second->secondMate,"NH",'i',4,(uint8_t*)&nh);

                        kv.second->secondMate->core.bin = hts_reg2bin(kv.second->secondMate->core.pos, kv.second->secondMate->core.pos + bam_cigar2rlen(kv.second->secondMate->core.n_cigar, bam_get_cigar(kv.second->secondMate)), 14, 5);
                        ret_sam=sam_write1(outSAM,genome_al_hdr,kv.second->secondMate);

                        delete kv.second;
                    }
                }
                else{
                    std::cout<<"MISTAKE. READS LEFT IN TMP"<<std::endl;
                }
            }
            else if (!curReadGroup.empty()){
                int nh=curReadGroup.size();
                for (std::pair<std::string,bam1_t*> kv: curReadGroup){
                    uint8_t* ptr_nh=bam_aux_get(kv.second,"NH");
                    if(ptr_nh){
                        bam_aux_del(kv.second,ptr_nh);
                    }
                    bam_aux_append(kv.second,"NH",'i',4,(uint8_t*)&nh);

                    kv.second->core.bin = hts_reg2bin(kv.second->core.pos, kv.second->core.pos + bam_cigar2rlen(kv.second->core.n_cigar, bam_get_cigar(kv.second)), 14, 5);
                    int ret_sam=sam_write1(outSAM,genome_al_hdr,kv.second);
                    delete kv.second;
                }
            }
            else{
                std::cerr<<"wrong!!!"<<std::endl;
                exit(-1);
            }

            // prepare for the next round
            curReadGroup.clear();
            curReadGroup_tmp.clear();
            curReadGroup_paired.clear();

            curReadName=newReadName;

            success=true;
            size_t read_start=0;
            size_t mate_read_start=0;

            const char *target_name=al_hdr->target_name[curAl->core.tid];
            int trans_idx=atoi(target_name);
            p_trans=tidx_to_t[target_name];

            GVec<GSeg>& exon_list=p_trans->exons;

            GSeg *next_exon=NULL;
            int cigars[MAX_CIGARS];

            int cur_pos,match_length,miss_length,cur_intron_len=0,i=0,num_cigars=0;

            bool ret_val=Map2GFF::get_read_start(exon_list,curAl->core.pos,read_start,i);
            if (!ret_val){
                std::cerr<<"SOMETHING WRONG WITH GETTING GENOMIC READ START"<<std::endl;
            }

            std::string cigar_str="";
            int cigar_ret_val=Map2GFF::convert_cigar(i,cur_intron_len,miss_length,next_exon,match_length,exon_list,num_cigars,read_start,curAl,cigar_str,cigars);
            if (cigar_ret_val){
                bool paired=curAl->core.flag           &    1;
                bool pairedly_aligned=curAl->core.flag &    2;
                // evaluate as a pair if a singleton is encountered - forget about it
                int mate_i=0;
                const char *mate_target_name=al_hdr->target_name[curAl->core.mtid];
                int mate_trans_idx=atoi(mate_target_name);
                mate_p_trans=tidx_to_t[mate_target_name];
                bool mate_ret_val=Map2GFF::get_read_start(exon_list,curAl->core.mpos,mate_read_start,mate_i);
                // if (!mate_ret_val){
                //     std::cerr<<"SOMETHING WRONG WITH GETTING GENOMIC MATE READ START"<<std::endl;
                //     // need aditional cases here
                // }

                int origRefID=curAl->core.tid, origPos=curAl->core.pos, origMPos=curAl->core.mpos;
                curAl->core.tid=ref_to_id[p_trans->refID];
                curAl->core.pos=read_start-1;
                if (ret_val){
                    curAl->core.mtid=ref_to_id[mate_p_trans->refID];
                    curAl->core.mpos=mate_read_start-1;
                }
                else{
                    curAl->core.mtid=ref_to_id[p_trans->refID];
                    curAl->core.mpos=read_start-1;
                }

                int old_n_cigar=curAl->core.n_cigar;
                if (num_cigars != old_n_cigar){
                    int data_len=curAl->l_data+4*(num_cigars-old_n_cigar);
                    int m_data=std::max(data_len,(int)curAl->m_data);
                    kroundup32(m_data);

                    uint8_t* data = (uint8_t*)calloc(m_data,1);

                    int copy1_len = (uint8_t*)bam_get_cigar(curAl) - curAl->data;
                    memcpy(data, curAl->data, copy1_len);

                    int copy2_len = num_cigars * 4;
                    memcpy(data + copy1_len, cigars, copy2_len);

                    int copy3_len = curAl->l_data - copy1_len - (old_n_cigar * 4);
                    memcpy(data + copy1_len + copy2_len, bam_get_seq(curAl), copy3_len);

                    curAl->core.n_cigar = num_cigars;

                    free(curAl->data);
                    curAl->data = data;
                    curAl->l_data = data_len;
                    curAl->m_data = m_data;
                }
                
                char strand=p_trans->strand;
                uint8_t* ptr=bam_aux_get(curAl,"XS");
                if(ptr){
                    bam_aux_del(curAl,ptr);
                }
                if (strand=='-'){
                    bam_aux_append(curAl,"XS",'A',1,(const unsigned char*)"-");
                }
                if (strand=='+'){
                    bam_aux_append(curAl,"XS",'A',1,(const unsigned char*)"+");
                }

                // add optional tag with transcript ID
                uint8_t* ptr_op=bam_aux_get(curAl,"OP");
                if(ptr_op){
                    bam_aux_del(curAl,ptr_op);
                }
                bam_aux_append(curAl,"OP",'i',4,(uint8_t*)&origRefID);

                if (mate_ret_val && paired && pairedly_aligned){
                    if (curAl->core.flag&64){
                        std::string tmpKey_pre=std::to_string(origRefID)+"$"+std::to_string(origPos)+"_"+std::to_string(origMPos);
                        curReadGroup_tmp[tmpKey_pre]=new MateRead(bam_dup1(curAl),read_start,cigar_str);;
                    }
                    else if (curAl->core.flag&128){
                        std::string tmpKey_pre=std::to_string(origRefID)+"$"+std::to_string(origMPos)+"_"+std::to_string(origPos);
                        curReadGroup_tmp[tmpKey_pre]=new MateRead(bam_dup1(curAl),read_start,cigar_str);;
                    }
                    else{
                        std::cout<<"mate_error"<<std::endl;
                        exit(-1);
                    }
                }
                else{
                    std::string tmpKey=p_trans->refID+"$"+std::to_string(read_start)+"_"+cigar_str;
                    curReadGroup[tmpKey]=bam_dup1(curAl);
                }
            }
        }
    }

    sam_close(outSAM);
}

bool Map2GFF::get_read_start(GVec<GSeg>& exon_list,size_t gff_start,size_t& genome_start, int& exon_idx){
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

Map2GFF::~Map2GFF(){
    bam_destroy1(curAl);
    sam_close(al);
    bam_hdr_destroy(al_hdr);
    bam_hdr_destroy(genome_al_hdr);
}
