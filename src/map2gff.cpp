#include "map2gff.h"
#include "tokenize.h"

void tline_parserr(const std::string& tline, std::string add="") {
	std::cerr << "Error at parsing .tlst line " << add << ":"
			<< std::endl << '\t' << tline << std::endl;
	exit(1);
}

coord_range make_coord_range(int strand, int lower, int upper) {
            if(upper < lower) {
                return std::make_tuple(strand,upper,lower);
            }
            return std::make_tuple(strand,lower,upper);
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

Map2GFF::Map2GFF(const std::string& tlstFP, const std::string& alFP, const std::string& multiFP, const std::string& glstFP){
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
            GffTranscript* t=new GffTranscript(tline);
            transcripts.Add(t);
            tidx_to_t[std::to_string(t->numID)]=t;
        }
    }
    tlststream.close();
    std::cout<<"Loaded Transcript Data"<<std::endl;

    // reading in the multiFP file
    multistream.open(multiFP.c_str(),std::ios::in);
    if(!multistream.good()){
        std::cerr<<"FATAL: Couldn't open multimapper data: "<<multiFP<<std::endl;
        exit(1);
    }

    std::ios::sync_with_stdio(false);
    std::cout<<"Reading the Multimapper data: "<<multiFP<<std::endl;
    
    // GffTranscript t=transcripts[0]// single gfftranscript record used to access the chromosome map
    std::string mline;

    std::pair< std::map<
                    std::vector< std::pair< int,int> >,
                    std::vector<const std::vector<std::pair<int,int> >*>
                >::iterator,bool> exists_cur_coord;

    std::pair<std::unordered_map<std::string,int>::iterator,bool> exists_ri;
    std::pair<std::unordered_map<int,std::string>::iterator,bool> exists_ir;

    std::vector<std::pair<int,int> > cur_coords1,cur_coords2;
    std::stringstream ss(""),sub_ss("");
    std::string pretab,posttab,sub;
    int chrid,strand,start,end;

    int max_chrID=0;
    int count=0;
    while (std::getline(multistream,mline)){
        count++;
        if(count%10000==0){
            std::cout<<count<<std::endl;
        }
        ss.str(mline);
        ss.clear();

        std::getline(ss,pretab,'\t');
        std::getline(ss,posttab,'\t');

        ss.str(pretab);
        ss.clear();
        std::getline(ss,pretab,':');
        exists_ri=this->ref_to_id_mult.insert(std::make_pair(pretab,max_chrID));
        if(exists_ri.second){
            exists_ir=this->id_to_ref_mult.insert(std::make_pair(max_chrID,pretab));
            if(!exists_ir.second){
                std::cerr<<"ERROR. chromosome IDs messed up"<<std::endl;

            }
            max_chrID++;
        }
        std::getline(ss,pretab,'@');

        // unordered
        strand=int(pretab[0]);
        cur_coords1.push_back(std::make_pair(exists_ri.first->second,strand));
        // 3. get the coordinates
        while(std::getline(ss,pretab,',')){
            sub_ss.str(pretab);
            sub_ss.clear();
            // get start coordinate
            std::getline(sub_ss,sub,'-');
            start=atoi(sub.c_str());
            // get end coordinate
            std::getline(sub_ss,sub,'-');
            end=atoi(sub.c_str());
            cur_coords1.push_back(std::make_pair(start,end));
        }
        //now for the matched kmer coordinates
        ss.str(posttab);
        ss.clear();
        std::getline(ss,posttab,':');
        exists_ri=this->ref_to_id_mult.insert(std::make_pair(posttab,max_chrID));
        if(exists_ri.second){
            exists_ir=this->id_to_ref_mult.insert(std::make_pair(max_chrID,posttab));
            if(!exists_ir.second){
                std::cerr<<"ERROR. chromosome IDs messed up"<<std::endl;

            }
            max_chrID++;
        }
        std::getline(ss,posttab,'@');
        strand=int(posttab[0]);
        cur_coords2.push_back(std::make_pair(exists_ri.first->second,strand));
        // 3. get the coordinates
        while(std::getline(ss,posttab,',')){
            sub_ss.str(posttab);
            sub_ss.clear();
            // get start coordinate
            std::getline(sub_ss,sub,'-');
            start=atoi(sub.c_str());
            // get end coordinate
            std::getline(sub_ss,sub,'-');
            end=atoi(sub.c_str());
            cur_coords2.push_back(std::make_pair(start,end));
        }

        std::vector<std::pair<int,int>>* new_coords=new std::vector<std::pair<int,int>>(cur_coords2);
        exists_cur_coord=this->multimappers.insert(std::pair<std::vector<std::pair<int,int>>,std::vector<const std::vector<std::pair<int,int> >* > >(cur_coords1,{}));
        exists_cur_coord.first->second.push_back(new_coords);
        cur_coords1.clear();
        cur_coords2.clear();
    }
    multistream.close();
    std::cout<<"Loaded Multimapper Data"<<std::endl;

    // reading in the gene coordinate file
    glststream.open(glstFP.c_str(),std::ios::in);
    if(!glststream.good()){
        std::cerr<<"FATAL: Couldn't open gene coordinate data: "<<glstFP<<std::endl;
        exit(1);
    }

    std::ios::sync_with_stdio(false);
    std::cout<<"Reading the Gene Coordinate data: "<<glstFP<<std::endl;
    
    // GffTranscript t=transcripts[0]// single gfftranscript record used to access the chromosome map
    std::string gline;
    count=0;
    std::string geneID,gstrand,gstart,gend;
    while (std::getline(glststream,gline)){
        ss.str(gline);
        ss.clear();

        std::getline(ss,geneID,'\t');
        std::getline(ss,gstrand,'\t');
        std::getline(ss,gstart,'\t');
        std::getline(ss,gend,'\t');
        this->gene_coords.insert(std::make_pair(make_coord_range(atoi(gstrand.c_str()),atoi(gstart.c_str()),atoi(gend.c_str())),count));
        count++;
    }
    glststream.close();
    std::cout<<"Loaded Gene Coordinate Data"<<std::endl;

    // let's access the chromosome map through the transcript data as
    // p_trans->names.getSeqName()
    
    // now initialize the sam file
    al=hts_open(alFP.c_str(),"r");
    al_hdr=sam_hdr_read(al); // read the alignment header
    curAl=bam_init1(); // initialize the alignment
}

int Map2GFF::convert_cigar(int i,int cur_intron_len,int miss_length,GSeg *next_exon,int match_length,GVec<GSeg>& exon_list,
                           int &num_cigars,int read_start,bam1_t* curAl,std::string &cigar_str,int cigars[MAX_CIGARS],std::vector<std::pair<int,int>> &coords){
    
    coords.push_back(std::make_pair(read_start,0)); // add read start to the list of coordinates. end of the segment will have to be changed later

    int cur_total_pos=read_start; // same as cur_pos but includes the soft clipping bases
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
                cigar_str+=std::to_string(bam_cigar_oplen(cigars[num_cigars]));
                cigar_str+=bam_cigar_opchr(opcode);
                ++num_cigars;
                cur_pos+=remaining_length;
                cur_total_pos+=remaining_length;
                coords.back().second=cur_total_pos;
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
                    coords.back().second=cur_total_pos;
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
                // std::cout<<"added: "<<bam_cigar_opchr(bam_cigar_op(cigars[num_cigars]))<<std::endl;
                ++num_cigars;

                cur_pos+=match_length+miss_length;
                cur_total_pos+=match_length+miss_length;

                coords.back().second=cur_exon.end; // add end of the current exon
                coords.push_back(std::make_pair(next_exon->start,0)); // add start of the new exon

                remaining_length-=match_length;
                assert(cur_pos == (int)next_exon->start);
            }
        }
    }
    coords.back().second=cur_total_pos; //add end of the read
    return 1;
}

void Map2GFF::add_multimapper_pair(const std::vector<std::pair<int,int>> *cor1, const std::vector<std::pair<int,int>> *cor2, bam1_t *al1, bam1_t *al2){
    std::pair<std::vector<std::pair<int,int>>,std::vector<std::pair<int,int>>> coor_key=std::make_pair(*cor1,*cor2);
    this->exists_curReadGroup_paired=curReadGroup_paired.insert(std::make_pair(coor_key,new MatePair(bam_dup1(al1),bam_dup1(al2))));
    if(exists_curReadGroup_paired.second){ // if the key did not previously exist - we can modify the alignment record and insert it into the map
        // first change the read start, chromosome - no need to modify the reverse/non-reverse in flag right now, since multimappers don't take that into account
        this->curReadGroup_paired[coor_key]->firstMate->core.pos=cor1->operator[](1).first; // set read start
        this->curReadGroup_paired[coor_key]->secondMate->core.pos=cor2->operator[](1).first; // set read start

        this->curReadGroup_paired[coor_key]->firstMate->core.tid=this->ref_to_id[this->id_to_ref_mult[cor1->operator[](0).first]]; // set sequence id
        this->curReadGroup_paired[coor_key]->secondMate->core.tid=this->ref_to_id[this->id_to_ref_mult[cor2->operator[](0).first]]; // set sequence id
        // second need to change cigar string

        int mult_num_cigars=0; // defines the current position in the cigar string to which we are adding
        int mult_cigars[MAX_CIGARS];
        uint32_t *cigar_full=bam_get_cigar(al1);
        int opcode=bam_cigar_op(cigar_full[0]);
        int oplength=0;
        if(opcode==BAM_CSOFT_CLIP){ // first see if there is soft clipping at the start of the read
            oplength=bam_cigar_oplen(cigar_full[0]);
            mult_cigars[mult_num_cigars]=opcode|(oplength<<BAM_CIGAR_SHIFT);
            mult_num_cigars++;
        }
        //next see if there is soft clipping at the end of the read
        opcode=bam_cigar_op(cigar_full[al1->core.n_cigar]);
        int last_oplength=0;
        if(opcode==BAM_CSOFT_CLIP){
            last_oplength=bam_cigar_oplen(cigar_full[al1->core.n_cigar]);
            // this soft clipping should be added at the very end
        }
        // now modify the end of the multimapper coordinates to accomodate for the potential soft clipping at the end
        std::pair<int,int> cur_pair;
        int last_pos=cor1->at(cor1->size()-1).second;
        int num_unused_seg=0;
        for(int curPairPos=cor1->size()-1;curPairPos>0;curPairPos--){ //move back through exons, and record the position of where the end is
            cur_pair=cor1->operator[](curPairPos);
            int cur_seg_length=cur_pair.second-cur_pair.first;
            if(cur_seg_length>last_oplength){
                last_pos=cor1->at(curPairPos).second-last_oplength;
                break;
            }
            else{
                last_oplength-=(cur_pair.second-cur_pair.first);
                num_unused_seg++;
            }
        }

        // the rest between the coordinates, should be filled in with the values from the multimapping coordinate vector
        // int bases_seen=oplength; // keep track of the number of bases of a multimapper that have been encountered so far
        bool firstmatch=true;
        for(int i=1;i<cor1->size()-num_unused_seg;i++){
            cur_pair=cor1->operator[](i);
            if(oplength+cur_pair.first>=cur_pair.second){
                oplength=oplength-(cur_pair.second-cur_pair.first);
                continue;
            }
            else{
                if(!firstmatch){
                    mult_cigars[mult_num_cigars]=BAM_CREF_SKIP|((cor1->operator[](i-1).second-cur_pair.first)<<BAM_CIGAR_SHIFT);
                    mult_num_cigars++;
                    // bases_seen+=(cor->operator[](i-1).second-cur_pair.first);
                }
                int cur_pair_len=((cur_pair.second-cur_pair.first)-oplength);
                if(i==cor1->size()-1){ // test if we've reached the last seq in coordinate vector and account for the possible soft clipping
                    mult_cigars[mult_num_cigars]=BAM_CMATCH|(cur_pair_len<<BAM_CIGAR_SHIFT);
                    // bases_seen+=cur_pair_len;
                    oplength=0; // reset since it has been used
                    mult_num_cigars++;
                    firstmatch=false;
                }
                mult_cigars[mult_num_cigars]=BAM_CMATCH|(cur_pair_len<<BAM_CIGAR_SHIFT);
                // bases_seen+=cur_pair_len;
                oplength=0; // reset since it has been used
                mult_num_cigars++;
                firstmatch=false; // next time we get here - we need to add an exon
            }
        }
        if(last_oplength>0){
            mult_cigars[mult_num_cigars]=BAM_CSOFT_CLIP|(last_oplength<<BAM_CIGAR_SHIFT);
            mult_num_cigars++;
        }

        int data_len=this->curReadGroup_paired[coor_key]->firstMate->l_data+4*(mult_num_cigars-al1->core.n_cigar);
        int m_data=std::max(data_len,(int)this->curReadGroup_paired[coor_key]->firstMate->m_data);
        kroundup32(m_data);

        uint8_t* data = (uint8_t*)calloc(m_data,1);

        int copy1_len = (uint8_t*)bam_get_cigar(this->curReadGroup_paired[coor_key]->firstMate) - this->curReadGroup_paired[coor_key]->firstMate->data;
        memcpy(data, this->curReadGroup_paired[coor_key]->firstMate->data, copy1_len);

        int copy2_len = mult_num_cigars * 4;
        memcpy(data + copy1_len, mult_cigars, copy2_len);

        int copy3_len = this->curReadGroup_paired[coor_key]->firstMate->l_data - copy1_len - (al1->core.n_cigar * 4);
        memcpy(data + copy1_len + copy2_len, bam_get_seq(this->curReadGroup_paired[coor_key]->firstMate), copy3_len);

        this->curReadGroup_paired[coor_key]->firstMate->core.n_cigar = mult_num_cigars;

        free(this->curReadGroup_paired[coor_key]->firstMate->data);
        this->curReadGroup_paired[coor_key]->firstMate->data = data;
        this->curReadGroup_paired[coor_key]->firstMate->l_data = data_len;
        this->curReadGroup_paired[coor_key]->firstMate->m_data = m_data;
        memset(mult_cigars,0,sizeof(mult_cigars));
        //========================================
        // Now repeat for mate 2
        //========================================

        mult_num_cigars=0; // defines the current position in the cigar string to which we are adding
        cigar_full=bam_get_cigar(al2);
        opcode=bam_cigar_op(cigar_full[0]);
        oplength=0;
        if(opcode==BAM_CSOFT_CLIP){ // first see if there is soft clipping at the start of the read
            oplength=bam_cigar_oplen(cigar_full[0]);
            mult_cigars[mult_num_cigars]=opcode|(oplength<<BAM_CIGAR_SHIFT);
            mult_num_cigars++;
        }
        //next see if there is soft clipping at the end of the read
        opcode=bam_cigar_op(cigar_full[al2->core.n_cigar]);
        last_oplength=0;
        if(opcode==BAM_CSOFT_CLIP){
            last_oplength=bam_cigar_oplen(cigar_full[al2->core.n_cigar]);
            // this soft clipping should be added at the very end
        }
        // now modify the end of the multimapper coordinates to accomodate for the potential soft clipping at the end
        last_pos=cor2->at(cor2->size()-1).second;
        num_unused_seg=0;
        for(int curPairPos=cor2->size()-1;curPairPos>0;curPairPos--){ //move back through exons, and record the position of where the end is
            cur_pair=cor2->operator[](curPairPos);
            int cur_seg_length=cur_pair.second-cur_pair.first;
            if(cur_seg_length>last_oplength){
                last_pos=cor2->at(curPairPos).second-last_oplength;
                break;
            }
            else{
                last_oplength-=(cur_pair.second-cur_pair.first);
                num_unused_seg++;
            }
        }

        // the rest between the coordinates, should be filled in with the values from the multimapping coordinate vector
        // int bases_seen=oplength; // keep track of the number of bases of a multimapper that have been encountered so far
        firstmatch=true;
        for(int i=1;i<cor2->size()-num_unused_seg;i++){
            cur_pair=cor2->operator[](i);
            if(oplength+cur_pair.first>=cur_pair.second){
                oplength=oplength-(cur_pair.second-cur_pair.first);
                continue;
            }
            else{
                if(!firstmatch){
                    mult_cigars[mult_num_cigars]=BAM_CREF_SKIP|((cor2->operator[](i-1).second-cur_pair.first)<<BAM_CIGAR_SHIFT);
                    mult_num_cigars++;
                    // bases_seen+=(cor->operator[](i-1).second-cur_pair.first);
                }
                int cur_pair_len=((cur_pair.second-cur_pair.first)-oplength);
                if(i==cor2->size()-1){ // test if we've reached the last seq in coordinate vector and account for the possible soft clipping
                    mult_cigars[mult_num_cigars]=BAM_CMATCH|(cur_pair_len<<BAM_CIGAR_SHIFT);
                    // bases_seen+=cur_pair_len;
                    oplength=0; // reset since it has been used
                    mult_num_cigars++;
                    firstmatch=false;
                }
                mult_cigars[mult_num_cigars]=BAM_CMATCH|(cur_pair_len<<BAM_CIGAR_SHIFT);
                // bases_seen+=cur_pair_len;
                oplength=0; // reset since it has been used
                mult_num_cigars++;
                firstmatch=false; // next time we get here - we need to add an exon
            }
        }
        if(last_oplength>0){
            mult_cigars[mult_num_cigars]=BAM_CSOFT_CLIP|(last_oplength<<BAM_CIGAR_SHIFT);
            mult_num_cigars++;
        }

        data_len=this->curReadGroup_paired[coor_key]->secondMate->l_data+4*(mult_num_cigars-al2->core.n_cigar);
        m_data=std::max(data_len,(int)this->curReadGroup_paired[coor_key]->secondMate->m_data);
        kroundup32(m_data);

        data = (uint8_t*)calloc(m_data,1);

        copy1_len = (uint8_t*)bam_get_cigar(this->curReadGroup_paired[coor_key]->secondMate) - this->curReadGroup_paired[coor_key]->secondMate->data;
        memcpy(data, this->curReadGroup_paired[coor_key]->secondMate->data, copy1_len);

        copy2_len = mult_num_cigars * 4;
        memcpy(data + copy1_len, mult_cigars, copy2_len);

        copy3_len = this->curReadGroup_paired[coor_key]->secondMate->l_data - copy1_len - (al2->core.n_cigar * 4);
        memcpy(data + copy1_len + copy2_len, bam_get_seq(this->curReadGroup_paired[coor_key]->secondMate), copy3_len);

        this->curReadGroup_paired[coor_key]->secondMate->core.n_cigar = mult_num_cigars;

        free(this->curReadGroup_paired[coor_key]->secondMate->data);
        this->curReadGroup_paired[coor_key]->secondMate->data = data;
        this->curReadGroup_paired[coor_key]->secondMate->l_data = data_len;
        this->curReadGroup_paired[coor_key]->secondMate->m_data = m_data;
    }
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

    bool start=true;

    while(sam_read1(al,al_hdr,curAl)>0){
        std::cout<<bam_get_qname(curAl)<<std::endl;
        std::string newReadName=bam_get_qname(curAl);
        // std::cout<<newReadName<<std::endl;
        if (newReadName.compare(curReadName)==0 || start){
            // if (newReadName.compare("read21068/1304")==0){
            //     std::cout<<"HELLO READ"<<std::endl;
            // }
            curReadName=newReadName;
            start=false;
            size_t read_start=0;
            size_t mate_read_start=0;

            const char *target_name=al_hdr->target_name[curAl->core.tid];
            p_trans=tidx_to_t[target_name];

            GVec<GSeg>& exon_list=p_trans->exons;

            GSeg *next_exon=NULL;
            int cigars[MAX_CIGARS];

            int cur_pos,match_length,miss_length,cur_intron_len=0,i=0,num_cigars=0;

            bool ret_val=Map2GFF::get_read_start(exon_list,curAl->core.pos,read_start,i);
            if (!ret_val){
                std::cerr<<"SOMETHING WRONG WITH GETTING GENOMIC READ START"<<std::endl;
            }
            // std::cout<<curAl->core.pos<<"\t"<<read_start<<std::endl;

            std::string cigar_str="";
            std::vector<std::pair<int,int>> al_coords={}; //coordinate vector of current alignment for lookup in multimappers
            al_coords.push_back(std::make_pair(ref_to_id_mult[p_trans->refID],p_trans->strand));
            int cigar_ret_val=Map2GFF::convert_cigar(i,cur_intron_len,miss_length,next_exon,match_length,exon_list,num_cigars,read_start,curAl,cigar_str,cigars,al_coords);
            if (cigar_ret_val){
                bool paired=curAl->core.flag           &    1;
                bool pairedly_aligned=curAl->core.flag &    2;
                // evaluate as a pair if a singleton is encountered - forget about it
                int mate_i=0;
                const char *mate_target_name=al_hdr->target_name[curAl->core.mtid];
                mate_p_trans=tidx_to_t[mate_target_name];
                bool mate_ret_val=Map2GFF::get_read_start(exon_list,curAl->core.mpos,mate_read_start,mate_i);
                // std::cout<<curAl->core.mpos<<"\t"<<mate_read_start<<std::endl;
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
                            curReadGroup_tmp[tmpKey_pre]=new MateRead(bam_dup1(curAl),read_start,cigar_str,al_coords);
                        }
                        else{
                            MateRead *mr=curReadGroup_tmp[tmpKey_pre];
                            if (mr->al->core.flag&128){
                                // std::string tmpKey=p_trans->refID+"$"+std::to_string(read_start)+"_"+cigar_str+"&"+std::to_string(mr->pos)+"_"+mr->cigar;
                                std::pair<std::vector<std::pair<int,int>>, std::vector<std::pair<int,int>>> tmpKey=std::make_pair(al_coords,mr->coords);
                                if (curReadGroup_paired.find(tmpKey)==curReadGroup_paired.end()){
                                    curReadGroup_paired[tmpKey]=new MatePair(bam_dup1(curAl),bam_dup1(mr->al));
                                    // this is where it now needs to search for the multimappers
                                    exists_al_coord_mate1=this->multimappers.find(al_coords);
                                    if(this->exists_al_coord_mate1!=this->multimappers.end()){ // found multimappers for one mate
                                        exists_al_coord_mate2=this->multimappers.find(mr->coords);
                                        if(this->exists_al_coord_mate2!=this->multimappers.end()){ // found multimappers for the second mate as well
                                            // now need to identify matches that belong to the same gene
                                            for(auto cor1:this->exists_al_coord_mate1->second){
                                                for(auto cor2:this->exists_al_coord_mate2->second){
                                                    this->exists_geneID_mate1=gene_coords.find(make_coord_range(cor1->operator[](0).second,cor1->operator[](1).first,cor1->operator[](1).first));
                                                    this->exists_geneID_mate2=gene_coords.find(make_coord_range(cor2->operator[](0).second,cor2->operator[](1).first,cor2->operator[](1).first));
                                                    if(this->exists_geneID_mate1!=this->gene_coords.end() &&
                                                       this->exists_geneID_mate2!=this->gene_coords.end() &&
                                                       this->exists_geneID_mate1->second==this->exists_geneID_mate2->second){ // both identified, and geneIDs match
                                                        // can now safely add to the multimappers
                                                        // std::cout<<"ADDING PAIR MULTI #1"<<std::endl;
                                                        add_multimapper_pair(cor1,cor2,curAl,mr->al);
                                                    }
                                                }                                            
                                            }
                                        }
                                    }
                                }
                                delete curReadGroup_tmp[tmpKey_pre];
                                curReadGroup_tmp.erase(tmpKey_pre);
                                // i think this is where we can evaluate the multimappers for a pair
                                // since this is the point when we know that two mates of a pair have been identified
                                // and we proceed to add them to the final map
                            }
                            else{
                                std::cerr<<"BIG MISTAKE AND DON'T KNOW WHAT TO DO"<<std::endl;
                            }
                        }
                    }
                    else if (curAl->core.flag&128){
                        std::string tmpKey_pre=std::to_string(origRefID)+"$"+std::to_string(origMPos)+"_"+std::to_string(origPos);
                        if (curReadGroup_tmp.find(tmpKey_pre)==curReadGroup_tmp.end()){
                            curReadGroup_tmp[tmpKey_pre]=new MateRead(bam_dup1(curAl),read_start,cigar_str,al_coords);
                        }
                        else{
                            // now write data where need be
                            MateRead *mr=curReadGroup_tmp[tmpKey_pre];
                            if (mr->al->core.flag&64){
                                // std::string tmpKey=p_trans->refID+"$"+std::to_string(mr->pos)+"_"+mr->cigar+"&"+std::to_string(read_start)+"_"+cigar_str;
                                std::pair<std::vector<std::pair<int,int>>, std::vector<std::pair<int,int>>> tmpKey=std::make_pair(mr->coords,al_coords);
                                if (curReadGroup_paired.find(tmpKey)==curReadGroup_paired.end()){
                                    curReadGroup_paired[tmpKey]=new MatePair(bam_dup1(mr->al),bam_dup1(curAl));
                                    // this is where it now needs to search for the multimappers
                                    exists_al_coord_mate2=this->multimappers.find(mr->coords);
                                    if(this->exists_al_coord_mate2!=this->multimappers.end()){ // found multimappers for one mate
                                        exists_al_coord_mate1=this->multimappers.find(al_coords);
                                        if(this->exists_al_coord_mate1!=this->multimappers.end()){ // found multimappers for the second mate as well
                                            // now need to identify matches that belong to the same gene
                                            for(auto cor2:this->exists_al_coord_mate2->second){
                                                for(auto cor1:this->exists_al_coord_mate1->second){
                                                    this->exists_geneID_mate2=gene_coords.find(make_coord_range(cor2->operator[](0).second,cor2->operator[](1).first,cor2->operator[](1).first));
                                                    this->exists_geneID_mate1=gene_coords.find(make_coord_range(cor1->operator[](0).second,cor1->operator[](1).first,cor1->operator[](1).first));
                                                    if(this->exists_geneID_mate2!=this->gene_coords.end() &&
                                                       this->exists_geneID_mate1!=this->gene_coords.end() &&
                                                       this->exists_geneID_mate2->second==this->exists_geneID_mate1->second){ // both identified, and geneIDs match
                                                        // can now safely add to the multimappers
                                                        // std::cout<<"ADDING PAIR MULTI #2"<<std::endl;
                                                        add_multimapper_pair(cor2,cor1,mr->al,curAl);
                                                    }
                                                }                                            
                                            }
                                        }
                                    }
                                }
                                delete curReadGroup_tmp[tmpKey_pre];
                                curReadGroup_tmp.erase(tmpKey_pre);
                                // i think this is where we can evaluate the multimappers for a pair
                                // since this is the point when we know that two mates of a pair have been identified
                                // and we proceed to add them to the final map
                            }
                            else{
                                std::cout<<"BIG MISTAKE #2 AND DON'T KNOW WHAT TO DO"<<std::endl;
                            }
                        }
                    }
                    else{
                        std::cerr<<"mate_error"<<std::endl;
                        exit(-1);
                    }
                }
                else{
                    curReadGroup[al_coords]=bam_dup1(curAl);
                    // now add the multimappers to the dictionary
                    //first check if any exist in the map
                    this->exists_al_coord=this->multimappers.find(al_coords);
                    if(this->exists_al_coord!=this->multimappers.end()){ // found multimappers and need to add them to the current groups
                        // first need to somehow evaluate whether both mates of the pair have been identified
                        for(auto cor:this->exists_al_coord->second){ // add each of the multimappers to the group
                            exists_curReadGroup=curReadGroup.insert(std::make_pair(*cor,bam_dup1(curAl)));
                            if(exists_curReadGroup.second){ // if the key did not previously exist - we can modify the alignment record and insert it into the map
                                // first change the read start, chromosome - no need to modify the reverse/non-reverse in flag right now, since multimappers don't take that into account
                                curReadGroup[*cor]->core.pos=cor->operator[](1).first; //read start
                                curReadGroup[*cor]->core.tid=ref_to_id[id_to_ref_mult[cor->operator[](0).first]];
                                // second need to change cigar string

                                int mult_num_cigars=0; // defines the current position in the cigar string to which we are adding
                                int mult_cigars[MAX_CIGARS];
                                uint32_t *cigar_full=bam_get_cigar(curAl);
                                int opcode=bam_cigar_op(cigar_full[0]);
                                int oplength=0;
                                if(opcode==BAM_CSOFT_CLIP){ // first see if there is soft clipping at the start of the read
                                    oplength=bam_cigar_oplen(cigar_full[0]);
                                    mult_cigars[mult_num_cigars]=opcode|(oplength<<BAM_CIGAR_SHIFT);
                                    mult_num_cigars++;
                                }
                                //next see if there is soft clipping at the end of the read
                                opcode=bam_cigar_op(cigar_full[curAl->core.n_cigar]);
                                int last_oplength=0;
                                if(opcode==BAM_CSOFT_CLIP){
                                    last_oplength=bam_cigar_oplen(cigar_full[curAl->core.n_cigar]);
                                    // this soft clipping should be added at the very end
                                }
                                // now modify the end of the multimapper coordinates to accomodate for the potential soft clipping at the end
                                std::pair<int,int> cur_pair;
                                int last_pos=cor->at(cor->size()-1).second;
                                int num_unused_seg=0;
                                for(int curPairPos=cor->size()-1;curPairPos>0;curPairPos--){ //move back through exons, and record the position of where the end is
                                    cur_pair=cor->operator[](curPairPos);
                                    int cur_seg_length=cur_pair.second-cur_pair.first;
                                    if(cur_seg_length>last_oplength){
                                        // cor->at(curPairPos).second-=last_oplength;
                                        last_pos=cor->at(curPairPos).second-last_oplength;
                                        break;
                                    }
                                    else{
                                        last_oplength-=(cur_pair.second-cur_pair.first);
                                        num_unused_seg++;
                                        // cor->pop_back();
                                    }
                                }

                                // the rest between the coordinates, should be filled in with the values from the multimapping coordinate vector
                                // int bases_seen=oplength; // keep track of the number of bases of a multimapper that have been encountered so far
                                bool firstmatch=true;
                                for(int i=1;i<cor->size()-num_unused_seg;i++){
                                    cur_pair=cor->operator[](i);
                                    if(oplength+cur_pair.first>=cur_pair.second){
                                        oplength=oplength-(cur_pair.second-cur_pair.first);
                                        continue;
                                    }
                                    else{
                                        if(!firstmatch){
                                            mult_cigars[mult_num_cigars]=BAM_CREF_SKIP|((cor->operator[](i-1).second-cur_pair.first)<<BAM_CIGAR_SHIFT);
                                            mult_num_cigars++;
                                            // bases_seen+=(cor->operator[](i-1).second-cur_pair.first);
                                        }
                                        int cur_pair_len=((cur_pair.second-cur_pair.first)-oplength);
                                        if(i==cor->size()-1){ // test if we've reached the last seq in coordinate vector and account for the possible soft clipping
                                            mult_cigars[mult_num_cigars]=BAM_CMATCH|(cur_pair_len<<BAM_CIGAR_SHIFT);
                                            // bases_seen+=cur_pair_len;
                                            oplength=0; // reset since it has been used
                                            mult_num_cigars++;
                                            firstmatch=false;
                                        }
                                        mult_cigars[mult_num_cigars]=BAM_CMATCH|(cur_pair_len<<BAM_CIGAR_SHIFT);
                                        // bases_seen+=cur_pair_len;
                                        oplength=0; // reset since it has been used
                                        mult_num_cigars++;
                                        firstmatch=false; // next time we get here - we need to add an exon
                                    }
                                }
                                if(last_oplength>0){
                                    mult_cigars[mult_num_cigars]=BAM_CSOFT_CLIP|(last_oplength<<BAM_CIGAR_SHIFT);
                                    mult_num_cigars++;
                                }

                                int data_len=curReadGroup[*cor]->l_data+4*(mult_num_cigars-curAl->core.n_cigar);
                                int m_data=std::max(data_len,(int)curReadGroup[*cor]->m_data);
                                kroundup32(m_data);

                                uint8_t* data = (uint8_t*)calloc(m_data,1);

                                int copy1_len = (uint8_t*)bam_get_cigar(curReadGroup[*cor]) - curReadGroup[*cor]->data;
                                memcpy(data, curReadGroup[*cor]->data, copy1_len);

                                int copy2_len = mult_num_cigars * 4;
                                memcpy(data + copy1_len, mult_cigars, copy2_len);

                                int copy3_len = curReadGroup[*cor]->l_data - copy1_len - (curAl->core.n_cigar * 4);
                                memcpy(data + copy1_len + copy2_len, bam_get_seq(curReadGroup[*cor]), copy3_len);

                                curReadGroup[*cor]->core.n_cigar = mult_num_cigars;

                                free(curReadGroup[*cor]->data);
                                curReadGroup[*cor]->data = data;
                                curReadGroup[*cor]->l_data = data_len;
                                curReadGroup[*cor]->m_data = m_data;
                            }
                            // std::cout<<this->id_to_ref_mult[cor->operator[](0).first]<<" "<<cor->operator[](0).second<<std::endl;
                        }

                        // std::cout<<"hello"<<std::endl;
                        // int cpt=0; //counter of the total number of bases
                        // bool cptf=false; //has the first value passed?
                        // for(auto v: al_coords){
                        //     if (cptf){
                        //         for(int ygh=v.first;ygh<v.second+1;ygh++){
                        //             cpt++;
                        //         }
                        //     }
                        //     else{
                        //         cptf=true;
                        //     }
                        // }
                        // std::cout<<"\n================ "<<p_trans->gffID<<std::endl;
                        // std::cout<<cpt<<" xxx ";
                        // std::cout<<read_start<<" : "<<cigar_str<<std::endl;
                        // for(auto v: al_coords){
                        //     std::cout<<v.first<<":"<<v.second<<" ; ";
                        // }
                        // std::cout<<std::endl;
                    }
                }
            }
        }
        else{
            if (!curReadGroup_paired.empty()){
                if (curReadGroup_tmp.empty()){ // Make sure nothing is left in the tmp group
                    // // first process to append multimappers
                    // for (std::pair<std::string,MatePair*> kv: curReadGroup_paired){
                    //     // get multimappers of both mates
                    //     exists_mate1=this->multimappers.find(kv.second->firstMate.coords);
                    //     exists_mate2=this->multimappers.find(kv.second->firstMate.coords);
                    // }
                    int nh=curReadGroup_paired.size();
                    bool prim=true;
                    for (std::pair<std::pair<std::vector<std::pair<int,int>>, std::vector<std::pair<int,int>>>,MatePair*> kv: curReadGroup_paired){
                        // Process first mate
                        uint8_t* ptr_nh_1=bam_aux_get(kv.second->firstMate,"NH");
                        if(ptr_nh_1){
                            bam_aux_del(kv.second->firstMate,ptr_nh_1);
                        }
                        bam_aux_append(kv.second->firstMate,"NH",'i',4,(uint8_t*)&nh);

                        kv.second->firstMate->core.bin = hts_reg2bin(kv.second->firstMate->core.pos, kv.second->firstMate->core.pos + bam_cigar2rlen(kv.second->firstMate->core.n_cigar, bam_get_cigar(kv.second->firstMate)), 14, 5);
                        
                        // set first alignment as primary and the rest to secondary
                        if (prim){
                            // if (nh>=2){
                            //     std::cout<<"FIRST "<<nh<<" "<<bam_get_qname(kv.second->firstMate)<<" "<<bam_flag2str(kv.second->firstMate->core.flag)<<" new: ";
                            // }
                            kv.second->firstMate->core.flag &= ~BAM_FSECONDARY;
                        }
                        else{
                            // if (nh>=2){
                            //     std::cout<<"SECOND "<<bam_get_qname(kv.second->firstMate)<<" "<<bam_flag2str(kv.second->firstMate->core.flag)<<" new: ";
                            // }
                            kv.second->firstMate->core.flag |= BAM_FSECONDARY;
                        }
                        int ret_sam=sam_write1(outSAM,genome_al_hdr,kv.second->firstMate);
                        // if (nh>=2){
                        //     std::cout<<bam_flag2str(kv.second->firstMate->core.flag)<<std::endl;
                        // }

                        // Now process second mate
                        uint8_t* ptr_nh_2=bam_aux_get(kv.second->secondMate,"NH");
                        if(ptr_nh_2){
                            bam_aux_del(kv.second->secondMate,ptr_nh_2);
                        }
                        bam_aux_append(kv.second->secondMate,"NH",'i',4,(uint8_t*)&nh);

                        kv.second->secondMate->core.bin = hts_reg2bin(kv.second->secondMate->core.pos, kv.second->secondMate->core.pos + bam_cigar2rlen(kv.second->secondMate->core.n_cigar, bam_get_cigar(kv.second->secondMate)), 14, 5);
                        
                        // set first alignment as primary and the rest to secondary
                        if (prim){
                            kv.second->secondMate->core.flag &= ~BAM_FSECONDARY;
                        }
                        else{
                            kv.second->secondMate->core.flag |= BAM_FSECONDARY;
                        }
                        ret_sam=sam_write1(outSAM,genome_al_hdr,kv.second->secondMate);

                        delete kv.second;
                        prim=false;
                    }
                }
                else{
                    std::cerr<<"MISTAKE. READS LEFT IN TMP"<<std::endl;
                }
            }
            else if (!curReadGroup.empty()){
                int nh=curReadGroup.size();
                bool prim=true;
                for (std::pair<std::vector<std::pair<int,int>>,bam1_t*> kv: curReadGroup){
                    uint8_t* ptr_nh=bam_aux_get(kv.second,"NH");
                    if(ptr_nh){
                        bam_aux_del(kv.second,ptr_nh);
                    }
                    bam_aux_append(kv.second,"NH",'i',4,(uint8_t*)&nh);

                    kv.second->core.bin = hts_reg2bin(kv.second->core.pos, kv.second->core.pos + bam_cigar2rlen(kv.second->core.n_cigar, bam_get_cigar(kv.second)), 14, 5);
                    if (prim){
                        kv.second->core.flag &= ~BAM_FSECONDARY;
                    }
                    else{
                        kv.second->core.flag |= BAM_FSECONDARY;
                    }
                    int ret_sam=sam_write1(outSAM,genome_al_hdr,kv.second);
                    delete kv.second;
                    prim=false;
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

            size_t read_start=0;
            size_t mate_read_start=0;

            const char *target_name=al_hdr->target_name[curAl->core.tid];
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
            std::vector<std::pair<int,int>> al_coords={}; //coordinate vector of current alignment for lookup in multimappers
            al_coords.push_back(std::make_pair(ref_to_id_mult[p_trans->refID],p_trans->strand));
            int cigar_ret_val=Map2GFF::convert_cigar(i,cur_intron_len,miss_length,next_exon,match_length,exon_list,num_cigars,read_start,curAl,cigar_str,cigars,al_coords);

            // test al_coords
            // if (al_coords.size()>2){
            //     int cpt=0; //counter of the total number of bases
            //     bool cptf=false; //has the first value passed?
            //     for(auto v: al_coords){
            //         if (cptf){
            //             for(int ygh=v.first;ygh<v.second+1;ygh++){
            //                 cpt++;
            //             }
            //         }
            //         else{
            //             cptf=true;
            //         }
            //     }
            //     std::string si="I";
            //     std::string sd="D";
            //     if((cpt!=77) && (cigar_str.find(si) == std::string::npos) && (cigar_str.find(sd) == std::string::npos)){
            //         std::cout<<"\n================ "<<p_trans->gffID<<std::endl;
            //         std::cout<<cpt<<" xxx ";
            //         std::cout<<read_start<<" : "<<cigar_str<<std::endl;
            //         for(auto v: al_coords){
            //             std::cout<<v.first<<":"<<v.second<<" ; ";
            //         }
            //     }
            // }

            if (cigar_ret_val){
                bool paired=curAl->core.flag           &    1;
                bool pairedly_aligned=curAl->core.flag &    2;
                // evaluate as a pair if a singleton is encountered - forget about it
                int mate_i=0;
                const char *mate_target_name=al_hdr->target_name[curAl->core.mtid];
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
                        curReadGroup_tmp[tmpKey_pre]=new MateRead(bam_dup1(curAl),read_start,cigar_str,al_coords);
                        // now add the multimappers to the dictionary
                        //first check if any exist in the map
                        // this->exists_al_coord=this->multimappers.find(al_coords);
                        // if(this->exists_al_coord!=this->multimappers.end()){ // found multimappers and need to add them to the current groups
                        //     // first need to somehow evaluate whether both mates of the pair have been identified
                        // }
                    }
                    else if (curAl->core.flag&128){
                        std::string tmpKey_pre=std::to_string(origRefID)+"$"+std::to_string(origMPos)+"_"+std::to_string(origPos);
                        curReadGroup_tmp[tmpKey_pre]=new MateRead(bam_dup1(curAl),read_start,cigar_str,al_coords);
                        // now add the multimappers to the dictionary
                        //first check if any exist in the map
                        // this->exists_al_coord=this->multimappers.find(al_coords);
                        // if(this->exists_al_coord!=this->multimappers.end()){ // found multimappers and need to add them to the current groups
                        //     // first need to somehow evaluate whether both mates of the pair have been identified

                        //     // when searching for the mate - can do the following
                        //     // extract the transcript id of the original multimapper and of the mate multimapper
                        //     // if both exist then need to figure out if the distance is good enough to be considered a multimapper
                        //     // we can see if they start position of the mate is within the coordinates of the first multimapper
                        //     // if so, then we can check the distance between them, to be smaller than n (max intron length)
                        // }
                    }
                    else{
                        std::cerr<<"mate_error"<<std::endl;
                        exit(-1);
                    }
                }
                else{
                    // std::string tmpKey=p_trans->refID+"$"+std::to_string(read_start)+"_"+cigar_str;
                    curReadGroup[al_coords]=bam_dup1(curAl);
                    // now add the multimappers to the dictionary
                    //first check if any exist in the map
                    this->exists_al_coord=this->multimappers.find(al_coords);
                    if(this->exists_al_coord!=this->multimappers.end()){ // found multimappers and need to add them to the current groups
                        // first need to somehow evaluate whether both mates of the pair have been identified
                        for(auto cor:this->exists_al_coord->second){ // add each of the multimappers to the group
                            exists_curReadGroup=curReadGroup.insert(std::make_pair(*cor,bam_dup1(curAl)));
                            if(exists_curReadGroup.second){ // if the key did not previously exist - we can modify the alignment record and insert it into the map
                                // first change the read start, chromosome - no need to modify the reverse/non-reverse in flag right now, since multimappers don't take that into account
                                curReadGroup[*cor]->core.pos=cor->operator[](1).first; //read start
                                curReadGroup[*cor]->core.tid=ref_to_id[id_to_ref_mult[cor->operator[](0).first]];
                                // second need to change cigar string

                                int mult_num_cigars=0; // defines the current position in the cigar string to which we are adding
                                int mult_cigars[MAX_CIGARS];
                                uint32_t *cigar_full=bam_get_cigar(curAl);
                                int opcode=bam_cigar_op(cigar_full[0]);
                                int oplength=0;
                                if(opcode==BAM_CSOFT_CLIP){ // first see if there is soft clipping at the start of the read
                                    oplength=bam_cigar_oplen(cigar_full[0]);
                                    mult_cigars[mult_num_cigars]=opcode|(oplength<<BAM_CIGAR_SHIFT);
                                    mult_num_cigars++;
                                }
                                //next see if there is soft clipping at the end of the read
                                opcode=bam_cigar_op(cigar_full[curAl->core.n_cigar]);
                                int last_oplength=0;
                                if(opcode==BAM_CSOFT_CLIP){
                                    last_oplength=bam_cigar_oplen(cigar_full[curAl->core.n_cigar]);
                                    // this soft clipping should be added at the very end
                                }
                                // now modify the end of the multimapper coordinates to accomodate for the potential soft clipping at the end
                                std::pair<int,int> cur_pair;
                                int last_pos=cor->at(cor->size()-1).second;
                                int num_unused_seg=0;
                                for(int curPairPos=cor->size()-1;curPairPos>0;curPairPos--){ //move back through exons, and record the position of where the end is
                                    cur_pair=cor->operator[](curPairPos);
                                    int cur_seg_length=cur_pair.second-cur_pair.first;
                                    if(cur_seg_length>last_oplength){
                                        // cor->at(curPairPos).second-=last_oplength;
                                        last_pos=cor->at(curPairPos).second-last_oplength;
                                        break;
                                    }
                                    else{
                                        last_oplength-=(cur_pair.second-cur_pair.first);
                                        num_unused_seg++;
                                        // cor->pop_back();
                                    }
                                }

                                // the rest between the coordinates, should be filled in with the values from the multimapping coordinate vector
                                // int bases_seen=oplength; // keep track of the number of bases of a multimapper that have been encountered so far
                                bool firstmatch=true;
                                for(int i=1;i<cor->size()-num_unused_seg;i++){
                                    cur_pair=cor->operator[](i);
                                    if(oplength+cur_pair.first>=cur_pair.second){
                                        oplength=oplength-(cur_pair.second-cur_pair.first);
                                        continue;
                                    }
                                    else{
                                        if(!firstmatch){
                                            mult_cigars[mult_num_cigars]=BAM_CREF_SKIP|((cor->operator[](i-1).second-cur_pair.first)<<BAM_CIGAR_SHIFT);
                                            mult_num_cigars++;
                                            // bases_seen+=(cor->operator[](i-1).second-cur_pair.first);
                                        }
                                        int cur_pair_len=((cur_pair.second-cur_pair.first)-oplength);
                                        if(i==cor->size()-1){ // test if we've reached the last seq in coordinate vector and account for the possible soft clipping
                                            mult_cigars[mult_num_cigars]=BAM_CMATCH|(cur_pair_len<<BAM_CIGAR_SHIFT);
                                            // bases_seen+=cur_pair_len;
                                            oplength=0; // reset since it has been used
                                            mult_num_cigars++;
                                            firstmatch=false;
                                        }
                                        mult_cigars[mult_num_cigars]=BAM_CMATCH|(cur_pair_len<<BAM_CIGAR_SHIFT);
                                        // bases_seen+=cur_pair_len;
                                        oplength=0; // reset since it has been used
                                        mult_num_cigars++;
                                        firstmatch=false; // next time we get here - we need to add an exon
                                    }
                                }
                                if(last_oplength>0){
                                    mult_cigars[mult_num_cigars]=BAM_CSOFT_CLIP|(last_oplength<<BAM_CIGAR_SHIFT);
                                    mult_num_cigars++;
                                }

                                int data_len=curReadGroup[*cor]->l_data+4*(mult_num_cigars-curAl->core.n_cigar);
                                int m_data=std::max(data_len,(int)curReadGroup[*cor]->m_data);
                                kroundup32(m_data);

                                uint8_t* data = (uint8_t*)calloc(m_data,1);

                                int copy1_len = (uint8_t*)bam_get_cigar(curReadGroup[*cor]) - curReadGroup[*cor]->data;
                                memcpy(data, curReadGroup[*cor]->data, copy1_len);

                                int copy2_len = mult_num_cigars * 4;
                                memcpy(data + copy1_len, mult_cigars, copy2_len);

                                int copy3_len = curReadGroup[*cor]->l_data - copy1_len - (curAl->core.n_cigar * 4);
                                memcpy(data + copy1_len + copy2_len, bam_get_seq(curReadGroup[*cor]), copy3_len);

                                curReadGroup[*cor]->core.n_cigar = mult_num_cigars;

                                free(curReadGroup[*cor]->data);
                                curReadGroup[*cor]->data = data;
                                curReadGroup[*cor]->l_data = data_len;
                                curReadGroup[*cor]->m_data = m_data;
                            }
                            // std::cout<<this->id_to_ref_mult[cor->operator[](0).first]<<" "<<cor->operator[](0).second<<std::endl;
                        }

                        // std::cout<<"hello"<<std::endl;
                        // int cpt=0; //counter of the total number of bases
                        // bool cptf=false; //has the first value passed?
                        // for(auto v: al_coords){
                        //     if (cptf){
                        //         for(int ygh=v.first;ygh<v.second+1;ygh++){
                        //             cpt++;
                        //         }
                        //     }
                        //     else{
                        //         cptf=true;
                        //     }
                        // }
                        // std::cout<<"\n================ "<<p_trans->gffID<<std::endl;
                        // std::cout<<cpt<<" xxx ";
                        // std::cout<<read_start<<" : "<<cigar_str<<std::endl;
                        // for(auto v: al_coords){
                        //     std::cout<<v.first<<":"<<v.second<<" ; ";
                        // }
                        // std::cout<<std::endl;
                    }
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