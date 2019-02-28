//
//  gtfToFasta.cpp
//  TopHat
//
//  Created by Harold Pimentel on 10/26/11.
//
#include "GTFToFasta.h"

std::string GTFToFasta::get_exonic_sequence(GffObj &p_trans,
                                FastaRecord &rec, std::string& coords, int kmer_length, std::ofstream& multimap)
{
    GList<GffExon>& exon_list = p_trans.exons;

    std::string exon_seq("");
    size_t length;
    coords.clear();
    std::stringstream ss;
    // base algorithm for multimappers as of now
    // for each transcript
    // 1. get start-end coordinates of a kmer
    // 2. insert into kmer_coords
    // 3. check return value: if second element is true - did not exist before; if false - was already present, and we should skip it
    std::vector<std::pair<int,int>> cur_coords; // the coordinates encountered from the previous exons
    int cur_len=0; // current length left from the previous exon
    int cur_pos=0;
    std::pair<std::set<std::vector<std::pair<int,int>>>::iterator,bool> exists;
    std::pair<std::unordered_map<std::string,std::vector<std::vector<std::pair<int,int>>>>::iterator,bool> exists2;
    std::string sub_seq("");
    cur_coords.push_back(std::make_pair(p_trans.gseq_id,(int)p_trans.strand));
    for (int i = 0; i < exon_list.Count(); ++i) {
        GffExon& cur_exon = *(exon_list.Get(i));
        length = cur_exon.end - cur_exon.start + 1;

        // get coordinates into the map
        if (length>1){ // sanity check for 0 and 1 baes exons
            for(int j=0;j<length;j++){ // iterate over all kmers in the given exon
                if ((length-j)+cur_len<kmer_length){ // not enough coordinates - need to look at the next exon
                    cur_coords.push_back(std::pair<int,int>(cur_exon.start+j,cur_exon.end));
                    cur_len+=(cur_exon.end-(cur_exon.start+j)); // save the length that has been seen thus far
                    break;
                }
                else{ // otherwise we have all the bases we need and can evaluate their uniqueness
                    if (cur_len!=0){ // some information is left from previous exons
                        for (int g=kmer_length-cur_len;g<kmer_length;g++){ // build new sequences using past coordinates
                            sub_seq="";
                            cur_coords.push_back(std::make_pair(cur_exon.start+j,cur_exon.start+j+g));
                            exists=this->kmer_coords.insert(cur_coords);
                            if (exists.second){
                                for (int d=1;d<cur_coords.size();d++){
                                    sub_seq+=rec.seq_.substr(cur_coords[d].first-1,cur_coords[d].second-cur_coords[d].first);
                                }
                                if(sub_seq.length()>kmer_length){
                                    std::cerr<<"1: "<<sub_seq.length()<<" "<<p_trans.getID()<<std::endl;
                                    for(int y=1;y<cur_coords.size();y++){
                                        std::cerr<<cur_coords[y].first<<"-"<<cur_coords[y].second<<";";
                                    }
                                    std::cerr<<std::endl;
                                }
                                std::transform(sub_seq.begin(), sub_seq.end(), sub_seq.begin(), ::toupper);
                                exists2=this->kmers.insert(std::pair<std::string,std::vector<std::vector<std::pair<int,int>>>>(sub_seq,{}));
                                exists2.first->second.push_back(cur_coords);
                                if(!exists2.second){
                                    // if exists2 is false, meaning the element already existed
                                    // we can then proceed to write the coordinates to the output file

                                    // for each previously added element
                                    // output a combination with every other element
                                    int u1=exists2.first->second.size()-1;
                                    for (int u2=0;u2<u1;u2++){
                                        multimap<<p_trans.names->gseqs.getName(exists2.first->second[u1][0].first)<<":"<<(char)exists2.first->second[u1][0].second<<"@";
                                        for (int k1=1;k1<exists2.first->second[u1].size();k1++){
                                            if (k1!=exists2.first->second[u1].size()-1){
                                                multimap<<exists2.first->second[u1][k1].first<<"-"<<exists2.first->second[u1][k1].second<<",";
                                            }
                                            else{
                                                multimap<<exists2.first->second[u1][k1].first<<"-"<<exists2.first->second[u1][k1].second;
                                            }
                                        }
                                        multimap<<"\t"<<p_trans.names->gseqs.getName(exists2.first->second[u2][0].first)<<":"<<(char)exists2.first->second[u2][0].second<<"@";
                                        for (int k2=1;k2<exists2.first->second[u2].size();k2++){
                                            if (k2!=exists2.first->second[u2].size()-1){
                                                multimap<<exists2.first->second[u2][k2].first<<"-"<<exists2.first->second[u2][k2].second<<",";
                                            }
                                            else{
                                                multimap<<exists2.first->second[u2][k2].first<<"-"<<exists2.first->second[u2][k2].second;
                                            }
                                        }
                                        multimap<<std::endl;
                                    }
                                }
                            }
                            cur_len-=1;
                            cur_coords.pop_back();
                            cur_coords[1].first+=1;
                            if (cur_coords[1].first==cur_coords[1].second){
                                cur_coords.erase(cur_coords.begin()+1); // delete the first element
                            }
                        }
                        // add new coordinates first
                    }
                    // need to resume from the current index
                    if ((length-j)+cur_len<kmer_length){ // not enough coordinates - need to look at the next exon
                        cur_coords.push_back(std::pair<int,int>(cur_exon.start+j,cur_exon.end));
                        cur_len+=(cur_exon.end-(cur_exon.start+j));
                        break;
                    }
                    else{
                        cur_coords.push_back(std::make_pair(cur_exon.start+j,cur_exon.start+j+kmer_length));
                        exists=this->kmer_coords.insert(cur_coords);
                        if (exists.second){ // was successfully inserted
                            // get sequence
                            sub_seq=rec.seq_.substr(cur_exon.start+j-1,kmer_length);
                            if(sub_seq.length()>kmer_length){
                                std::cerr<<"2: "<<sub_seq.length()<<" "<<p_trans.getID()<<" ";
                                for(int y=1;y<cur_coords.size();y++){
                                    std::cerr<<cur_coords[y].first<<"-"<<cur_coords[y].second<<";";
                                }
                                std::cerr<<std::endl;
                            }
                            std::transform(sub_seq.begin(), sub_seq.end(), sub_seq.begin(), ::toupper);
                            exists2=this->kmers.insert(std::pair<std::string,std::vector<std::vector<std::pair<int,int>>>>(sub_seq,{}));
                            exists2.first->second.push_back(cur_coords);
                            if(!exists2.second){
                                // if exists2 is false, meaning the element already existed
                                // we can then proceed to write the coordinates to the output file

                                // for each previously added element
                                // output a combination with every other element
                                int u1=exists2.first->second.size()-1; // need to only look at the current cur_coords, since otherwise we output duplicates
                                for (int u2=0;u2<u1;u2++){
                                    multimap<<p_trans.names->gseqs.getName(exists2.first->second[u1][0].first)<<":"<<(char)exists2.first->second[u1][0].second<<"@";
                                    for (int k1=1;k1<exists2.first->second[u1].size();k1++){
                                        if (k1!=exists2.first->second[u1].size()-1){
                                            multimap<<exists2.first->second[u1][k1].first<<"-"<<exists2.first->second[u1][k1].second<<",";
                                        }
                                        else{
                                            multimap<<exists2.first->second[u1][k1].first<<"-"<<exists2.first->second[u1][k1].second;
                                        }
                                    }
                                    multimap<<"\t"<<p_trans.names->gseqs.getName(exists2.first->second[u2][0].first)<<":"<<(char)exists2.first->second[u2][0].second<<"@";
                                    for (int k2=1;k2<exists2.first->second[u2].size();k2++){
                                        if (k2!=exists2.first->second[u2].size()-1){
                                            multimap<<exists2.first->second[u2][k2].first<<"-"<<exists2.first->second[u2][k2].second<<",";
                                        }
                                        else{
                                            multimap<<exists2.first->second[u2][k2].first<<"-"<<exists2.first->second[u2][k2].second;
                                        }
                                    }
                                    multimap<<std::endl;
                                }
                            }
                        }
                    }
                    // since we went into this conditional, that means no previously encountered exons are left
                    // and we can reset some things
                    cur_coords.erase(cur_coords.begin()+1,cur_coords.end());
                    cur_len=0;
                }
            }
        }
        exon_seq += rec.seq_.substr(cur_exon.start - 1, length);
        ss << ',' << cur_exon.start << '-' << cur_exon.end;
    }
    coords = ss.str().substr(1);
    return exon_seq;
}


GTFToFasta::GTFToFasta(std::string gtf_fname, std::string genome_fname)
{
    gtf_fname_ = gtf_fname;
    gtf_fhandle_ = fopen(gtf_fname_.c_str(), "r");
    if (gtf_fhandle_ == NULL)
    {
        std::cerr << "FATAL: Couldn't open annotation: " << gtf_fname_
        << std::endl;
        exit(1);
    }
    std::cout << "Reading the annotation file: " << gtf_fname_ << std::endl;
    gtfReader_.init(gtf_fhandle_, true); //load recognizable transcript features only
    gtfReader_.readAll();

    genome_fname_ = genome_fname;

    // Make a map from the GffObj
    transcript_map();
}

GTFToFasta::~GTFToFasta()
{
    ContigTransMap::iterator it;

    for (it = contigTransMap_.begin(); it != contigTransMap_.end(); ++it) {
        delete it->second;
    }

}

void GTFToFasta::make_transcriptome(std::string out_fname, int kmer_length)
{
    std::vector<int> *p_contig_vec;

    FastaReader fastaReader(genome_fname_);
    FastaWriter fastaWriter(out_fname);
    std::string tlst_fname(out_fname);
    tlst_fname.append(".tlst");
    std::ofstream tlst(tlst_fname.c_str());
    std::string multimap_fname(out_fname);
    multimap_fname.append(".multi");
    std::ofstream multimap(multimap_fname.c_str());
    FastaRecord cur_contig;

    std::map<std::string,std::tuple<int,int,int>> geneMap; // stores minimum and maximum gene coordinates of transcripts in a given gene
    std::pair< std::map<
                std::string,
                std::tuple<int,int,int>
            >::iterator,bool> exists_cur_gene;
    while (fastaReader.good()) {

        fastaReader.next(cur_contig);

        // If this contig isn't in the map, then there are no transcripts
        // associated with it. Skip it.
        if (contigTransMap_.find(cur_contig.id_) ==
            contigTransMap_.end())
        {
            continue;
        }

        p_contig_vec = contigTransMap_[cur_contig.id_];

        FastaRecord out_rec;
        for (size_t i = 0; i < p_contig_vec->size(); ++i) {
            int trans_idx = (*p_contig_vec)[i];
            GffObj *p_trans = gtfReader_.gflst.Get(trans_idx);
            //if (p_trans->isDiscarded() || p_trans->exons.Count()==0) continue;
            std::string coordstr;
            out_rec.seq_ = get_exonic_sequence(*p_trans, cur_contig, coordstr, kmer_length,multimap);
            if (out_rec.seq_.empty()) continue;
            std::stringstream ss;
            ss << trans_idx;
            out_rec.id_ = ss.str();
            //ss.str(std::string()); //clear ss
            out_rec.desc_=p_trans->getID();
            out_rec.desc_.push_back(' ');
            //ss << p_trans->getID() << ' ' << p_trans->getGSeqName() << p_trans->strand << '\t' << coordstr ;
            out_rec.desc_.append(cur_contig.id_);
            out_rec.desc_.push_back(p_trans->strand);
            out_rec.desc_.push_back(' ');
            out_rec.desc_.append(coordstr); //list of exon coordinates
            tlst << out_rec.id_ << ' ' << out_rec.desc_ << std::endl;
            //out_rec.desc_ = "";
            //out_rec.desc_ = ss.str();
            //out_rec.seq_ = exon_seq;
            fastaWriter.write(out_rec);

            // populate the genemap
            int nst=p_trans->start;
            int nen=p_trans->end;
            char* geneID=p_trans->getGeneID();
            if(geneID==NULL){
                geneID=p_trans->getGeneName();
                if(geneID==NULL){
                    std::cout<<"wrong geneID"<<std::endl;
                    continue;
                }
            }
            exists_cur_gene=geneMap.insert(std::make_pair(std::string(geneID),std::make_tuple(p_trans->strand,nst,nen)));
            if(!exists_cur_gene.second){ // the key did exist - update start and min
                int &st=std::get<1>(exists_cur_gene.first->second);
                int &en=std::get<2>(exists_cur_gene.first->second);
                if(nen>en){ //update end
                    en=nen;
                }
                if(nst<st){ //update start
                    st=nst;
                }
            }
        }
    }
    tlst.close();
    multimap.close();

    // write genes to file
    std::string gene_fname(out_fname);
    gene_fname.append(".glst");
    std::ofstream genefp(gene_fname.c_str());

    std::map<std::string,std::tuple<int,int,int>>::iterator it=geneMap.begin();
    while(it!=geneMap.end()){
        genefp<<it->first<<"\t"<<std::get<0>(it->second)<<"\t"<<std::get<1>(it->second)<<"\t"<<std::get<2>(it->second)<<std::endl;
        // std::cout<<it->first<<"\t"<<std::get<0>(it->second)<<"\t"<<std::get<1>(it->second)<<"\t"<<std::get<2>(it->second)<<std::endl;
        it++;
    }
    genefp.close();
}

void GTFToFasta::transcript_map()
{
    GffObj *p_gffObj;
    const char *p_contig_name;
    std::vector<int> *p_contig_vec;

    for (int i = 0; i < gtfReader_.gflst.Count(); ++i) 
    {
        p_gffObj = gtfReader_.gflst.Get(i);
        if (p_gffObj->isDiscarded() || p_gffObj->exons.Count()==0)
           continue;

        p_contig_name = p_gffObj->getRefName();
        std::string contig_name(p_contig_name);

        // Check if the current contig exists in the map
        // If it doesn't, add it
        if (contigTransMap_.find(contig_name) == contigTransMap_.end())
        {
            p_contig_vec = new std::vector<int>;
            contigTransMap_[contig_name] = p_contig_vec;
        }
        else
        {
            p_contig_vec = contigTransMap_[contig_name];
        }

        p_contig_vec->push_back(i);
    }
}

void GTFToFasta::print_mapping()
{
    std::ofstream out_file("out.names");
    GffObj *p_gffObj;

    for (int i = 0; i < gtfReader_.gflst.Count(); ++i) 
    {
        p_gffObj = gtfReader_.gflst.Get(i);
        if (p_gffObj->isDiscarded() || p_gffObj->exons.Count()==0) continue;
        out_file << i << "\t" << p_gffObj->getID() << std::endl;
    }

    out_file.close();
}

void gtf2fasta_print_usage() 
{
    std::cerr << "Usage: gtf_to_fasta kmer_length transcripts.gtf genome.fa out_file" << std::endl;
}

enum Opt {GFF_FP   = 'a',
          REF_FA     = 'r',
          OUT_FA    = 'o',
          KMER_LEN     = 'k'};

int main(int argc, char *argv[])
{
    ArgParse args("Map2GFF");
    args.add_string(Opt::GFF_FP,"gff","","");
    args.add_string(Opt::REF_FA,"ref","","");
    args.add_string(Opt::OUT_FA,"output","","");
    args.add_int(Opt::KMER_LEN,"kmer",76,"");
    
    args.parse_args(argc,argv);
    
    int kmer_length=args.get_int(Opt::KMER_LEN);
    std::string gtf_fname(args.get_string(Opt::GFF_FP));
    std::string genome_fname(args.get_string(Opt::REF_FA));
    std::string out_fname(args.get_string(Opt::OUT_FA));

    GTFToFasta gtfToFasta(gtf_fname, genome_fname);
    gtfToFasta.make_transcriptome(out_fname, kmer_length);
    //gtfToFasta.print_mapping();
    
    return 0;
}
