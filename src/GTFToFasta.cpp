#include <utility>

//
//  gtfToFasta.cpp
//  TopHat
//
//  Created by Harold Pimentel on 10/26/11.
//
#include "GTFToFasta.h"

std::string GTFToFasta::get_exonic_sequence(GffObj &p_trans,FastaRecord &rec, std::string& coords){
    GList<GffExon>& exon_list = p_trans.exons;

    std::string exon_seq;
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
    std::pair<std::map<std::vector<std::pair<int,int> >, std::pair<std::string,int> >::iterator,bool> exists;
    std::pair<std::unordered_map<std::string,std::vector<std::vector<std::pair<int,int>>>>::iterator,bool> exists2;
    std::string sub_seq;
    cur_coords.emplace_back(std::make_pair(p_trans.gseq_id,(int)p_trans.strand));
    bool exhaustedExon=false;
    for (int i = 0; i < exon_list.Count(); ++i) {
        GffExon& cur_exon = *(exon_list.Get(i));
        length = (cur_exon.end+1) - cur_exon.start;

        // get coordinates into the map
        if (length>1 and this->multi){ // sanity check for 0 and 1 base exons // TODO: need to refactor the multimapping index
            for(int j=0;j<length;j++){ // iterate over all kmers in the given exon
                if ((length-j)+cur_len<this->kmerlen){ // not enough coordinates - need to look at the next exon
                    cur_coords.emplace_back(std::pair<int,int>(cur_exon.start+j,cur_exon.end));
                    cur_len+=((cur_exon.end+1)-(cur_exon.start+j)); // save the length that has been seen thus far
                    break;
                }
                else{ // otherwise we have all the bases we need and can evaluate their uniqueness
                    if (cur_len!=0){ // some information is left from previous exons
                        for (int g=this->kmerlen-cur_len;g<this->kmerlen;g++){ // build new sequences using past coordinates
                            sub_seq="";
                            cur_coords.emplace_back(std::make_pair(cur_exon.start+j,cur_exon.start+j+g-1));
                            exists=this->kmer_coords.insert(std::make_pair(cur_coords,std::make_pair(p_trans.getID(),1)));
                            if (exists.second){ // if this set of genomic coordinates was not previously observed, we can proceed to evaluate
                                for (int d=1;d<cur_coords.size();d++){
                                    sub_seq+=rec.seq_.substr(cur_coords[d].first-1,(cur_coords[d].second+1)-cur_coords[d].first);
                                }
                                if(sub_seq.length()!=this->kmerlen){
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
                                        this->multimap<< GffObj::names->gseqs.getName(exists2.first->second[u1][0].first)<<":"<<(char)exists2.first->second[u1][0].second<<"@";
                                        for (int k1=1;k1<exists2.first->second[u1].size();k1++){
                                            if (k1!=exists2.first->second[u1].size()-1){
                                                this->multimap<<exists2.first->second[u1][k1].first<<"-"<<exists2.first->second[u1][k1].second<<",";
                                            }
                                            else{
                                                this->multimap<<exists2.first->second[u1][k1].first<<"-"<<exists2.first->second[u1][k1].second;
                                            }
                                        }
                                        this->multimap<<"\t"<< GffObj::names->gseqs.getName(exists2.first->second[u2][0].first)<<":"<<(char)exists2.first->second[u2][0].second<<"@";
                                        for (int k2=1;k2<exists2.first->second[u2].size();k2++){
                                            if (k2!=exists2.first->second[u2].size()-1){
                                                this->multimap<<exists2.first->second[u2][k2].first<<"-"<<exists2.first->second[u2][k2].second<<",";
                                            }
                                            else{
                                                this->multimap<<exists2.first->second[u2][k2].first<<"-"<<exists2.first->second[u2][k2].second;
                                            }
                                        }
                                        this->multimap<<std::endl;
                                    }
                                }
                            }
                            else{ // this is where we can perform evaluation of the number of unique kmers per transcript
                                // means that the genomic coordinates already existed before, and a counter needs to be updated
//                                std::cout<<"inc"<<std::endl;
                                exists.first->second.second++;
                            }
                            cur_len-=1;
                            cur_coords.pop_back();
                            cur_coords[1].first+=1;
                            if (cur_coords[1].first-1==cur_coords[1].second){
                                cur_coords.erase(cur_coords.begin()+1); // delete the first element
                            }
                            if((length - j) + cur_len < this->kmerlen){ // check again, since we are decreasing the cur_len with each iteration, thus it is possible that the length may not be sufficient
                                cur_coords.emplace_back(std::pair<int,int>(cur_exon.start+j,cur_exon.end));
                                cur_len += ((cur_exon.end+1) - (cur_exon.start + j)); // save the length that has been seen thus far
                                exhaustedExon=true; //set flag to exit from both for loops
                                break;
                            }
                        }
                        if(exhaustedExon){
                            exhaustedExon=false;
                            break;
                        }
                        // add new coordinates first
                    }
                    // need to resume from the current index
                    if ((length-j)+cur_len<this->kmerlen){ // not enough coordinates - need to look at the next exon
                        cur_coords.emplace_back(std::pair<int,int>(cur_exon.start+j,cur_exon.end));
                        cur_len+=((cur_exon.end+1)-(cur_exon.start+j));
                        break;
                    }
                    else{
                        cur_coords.emplace_back(std::make_pair(cur_exon.start+j,cur_exon.start+j+this->kmerlen));
                        exists=this->kmer_coords.insert(std::make_pair(cur_coords,std::make_pair(p_trans.getID(),1)));
                        if (exists.second){ // was successfully inserted
                            // get sequence
                            sub_seq=rec.seq_.substr(cur_exon.start+j-1,this->kmerlen);
                            if(sub_seq.length()>this->kmerlen){
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
                                    this->multimap<< GffObj::names->gseqs.getName(exists2.first->second[u1][0].first)<<":"<<(char)exists2.first->second[u1][0].second<<"@";
                                    for (int k1=1;k1<exists2.first->second[u1].size();k1++){
                                        if (k1!=exists2.first->second[u1].size()-1){
                                            this->multimap<<exists2.first->second[u1][k1].first<<"-"<<exists2.first->second[u1][k1].second<<",";
                                        }
                                        else{
                                            this->multimap<<exists2.first->second[u1][k1].first<<"-"<<exists2.first->second[u1][k1].second;
                                        }
                                    }
                                    this->multimap<<"\t"<< GffObj::names->gseqs.getName(exists2.first->second[u2][0].first)<<":"<<(char)exists2.first->second[u2][0].second<<"@";
                                    for (int k2=1;k2<exists2.first->second[u2].size();k2++){
                                        if (k2!=exists2.first->second[u2].size()-1){
                                            this->multimap<<exists2.first->second[u2][k2].first<<"-"<<exists2.first->second[u2][k2].second<<",";
                                        }
                                        else{
                                            this->multimap<<exists2.first->second[u2][k2].first<<"-"<<exists2.first->second[u2][k2].second;
                                        }
                                    }
                                    this->multimap<<std::endl;
                                }
                            }
                        }
                        else{
                            exists.first->second.second++;
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

GTFToFasta::GTFToFasta(std::string gtf_fname, std::string genome_fname,const std::string& out_fname, int kmerlen,bool multi)
{
    gtf_fname_ = gtf_fname;
    gtf_fhandle_ = fopen(gtf_fname_.c_str(), "r");
    if (gtf_fhandle_ == nullptr)
    {
        std::cerr << "FATAL: Couldn't open annotation: " << gtf_fname_
        << std::endl;
        exit(1);
    }
    std::cout << "Reading the annotation file: " << gtf_fname_ << std::endl;
    gtfReader_.init(gtf_fhandle_, true); //load recognizable transcript features only
    gtfReader_.readAll();
    std::cout << "loaded the annotation"<<std::endl;

    std::string multimap_fname(out_fname);
    multimap_fname.append(".multi");
    this->multimap = std::ofstream(multimap_fname.c_str());

    std::string tlst_fname(out_fname);
    tlst_fname.append(".tlst");
    this->tlst = std::ofstream(tlst_fname.c_str());

    std::string unique_fname(out_fname);
    unique_fname.append(".unq");
    this->uniquefp = std::ofstream(unique_fname.c_str());

    std::string gene_fname(out_fname);
    gene_fname.append(".glst");
    this->genefp = std::ofstream(gene_fname.c_str());

    std::string info_fname(out_fname);
    info_fname.append(".info");
    this->infofp = std::ofstream(info_fname.c_str());

    genome_fname_ = std::move(genome_fname);

    this->multi =  multi;
    this->kmerlen = kmerlen;
    this->out_fname = out_fname;

    // Make a map from the GffObj
    transcript_map();
}

GTFToFasta::~GTFToFasta(){
    ContigTransMap::iterator it;
    for (it = contigTransMap_.begin(); it != contigTransMap_.end(); ++it) {
        delete it->second;
    }
    this->infofp.close();
    this->tlst.close();
    this->multimap.close();
    this->uniquefp.close();
    this->genefp.close();
    this->out_file.close();
}

void GTFToFasta::make_transcriptome()
{
    std::vector<int> *p_contig_vec;

    FastaReader fastaReader(genome_fname_);
    FastaWriter fastaWriter(this->out_fname);
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
            contigTransMap_.end()){
            continue;
        }

        p_contig_vec = contigTransMap_[cur_contig.id_];

        FastaRecord out_rec;
        for (int trans_idx : *p_contig_vec) {
            GffObj *p_trans = gtfReader_.gflst.Get(trans_idx);
            std::string coordstr;
            out_rec.seq_ = get_exonic_sequence(*p_trans, cur_contig, coordstr);
            if (out_rec.seq_.empty()) continue;
            std::stringstream ss;
            ss << trans_idx;
            out_rec.id_ = ss.str();
            if(std::stoi(out_rec.id_)>this->topTransID){ // check if current id is greater than the highest previously observed
                topTransID = std::stoi(out_rec.id_);
            }
            out_rec.desc_=p_trans->getID();
            out_rec.desc_.push_back(' ');
            out_rec.desc_.append(cur_contig.id_);
            out_rec.desc_.push_back(p_trans->strand);
            out_rec.desc_.push_back(' ');
            out_rec.desc_.append(coordstr); //list of exon coordinates
            this->tlst << out_rec.id_ << ' ' << out_rec.desc_ << std::endl;
            fastaWriter.write(out_rec);

            // populate the genemap
            int nst=p_trans->start;
            int nen=p_trans->end;
            char* geneID=p_trans->getGeneID();
            if(geneID==nullptr){
                geneID=p_trans->getGeneName();
                if(geneID==nullptr){
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

    // now time to write the unique kmers for transcripts
    std::unordered_map<std::string,int> uniq_cnt; // counts of unique kmers per transcript
    std::pair<std::unordered_map<std::string,int>::iterator,bool> ex_ucnt; // exists or not

    auto it_unq = kmer_coords.begin();
    while(it_unq!=kmer_coords.end()){
        if(it_unq->second.second == 1){ // if only one instance was observed - write out
            ex_ucnt = uniq_cnt.insert(std::make_pair(it_unq->second.first,1));
            if(!ex_ucnt.second){ // did not exist
                ex_ucnt.first->second++;
            }
        }
        it_unq++;
    }

    auto it_unq_cnt = uniq_cnt.begin();
    while(it_unq_cnt!=uniq_cnt.end()){
        this->uniquefp<<it_unq_cnt->first<<","<<it_unq_cnt->second<<std::endl;
        it_unq_cnt++;
    }

    // write genes to file
    auto it=geneMap.begin();
    while(it!=geneMap.end()){
        this->genefp<<it->first<<"\t"<<std::get<0>(it->second)<<"\t"<<std::get<1>(it->second)<<"\t"<<std::get<2>(it->second)<<std::endl;
        it++;
    }

    // write the information about the index now
    this->infofp<<this->topTransID<<std::endl;
}

void GTFToFasta::transcript_map()
{
    GffObj *p_gffObj;
    const char *p_contig_name;
    std::vector<int> *p_contig_vec;

    for (int i = 0; i < gtfReader_.gflst.Count(); ++i){
        p_gffObj = gtfReader_.gflst.Get(i);
        if (p_gffObj->isDiscarded() || p_gffObj->exons.Count()==0){
           continue;
        }

        p_contig_name = p_gffObj->getRefName();
        std::string contig_name(p_contig_name);

        // Check if the current contig exists in the map
        // If it doesn't, add it
        if (contigTransMap_.find(contig_name) == contigTransMap_.end()){
            p_contig_vec = new std::vector<int>;
            contigTransMap_[contig_name] = p_contig_vec;
        }
        else{
            p_contig_vec = contigTransMap_[contig_name];
        }
        p_contig_vec->push_back(i);
    }
}

void GTFToFasta::print_mapping(){
    std::ofstream out_file("out.names");
    GffObj *p_gffObj;

    for (int i = 0; i < gtfReader_.gflst.Count(); ++i){
        p_gffObj = gtfReader_.gflst.Get(i);
        if (p_gffObj->isDiscarded() || p_gffObj->exons.Count()==0) continue;
        out_file << i << "\t" << p_gffObj->getID() << std::endl;
    }
}

void gtf2fasta_print_usage(){
    std::cerr << "Usage: gtf_to_fasta kmer_length transcripts.gtf genome.fa out_file" << std::endl;
}

enum Opt {GFF_FP   = 'a',
          REF_FA   = 'r',
          OUT_FA   = 'o',
          KMER_LEN = 'k',
          UNIQ     = 'u',
          MULTI    = 'm'};

int main(int argc, char *argv[])
{
    ArgParse args("Map2GFF");
    args.add_string(Opt::GFF_FP,"gff","","path to the annotation of the genome is GFF/GTF format",true);
    args.add_string(Opt::REF_FA,"ref","","path to the reference genome in the FASTA format",true);
    args.add_string(Opt::OUT_FA,"output","","base name for the output files",true);
    args.add_int(Opt::KMER_LEN,"kmer",76,"kmer length to use for building the index",false);
    args.add_flag(Opt::UNIQ,"uniq","get a separate output with uniq kmers per each transcript",false);
    args.add_flag(Opt::MULTI,"multi","identify all multimappers present in the input transcriptome",false);

    if(strcmp(argv[1],"--help")==0){
        std::cerr<<args.get_help()<<std::endl;
        exit(1);
    }
    
    args.parse_args(argc,argv);
    
    int kmer_length=args.get_int(Opt::KMER_LEN);
    std::string gtf_fname(args.get_string(Opt::GFF_FP));
    std::string genome_fname(args.get_string(Opt::REF_FA));
    std::string out_fname(args.get_string(Opt::OUT_FA));
    bool multi=args.get_flag(Opt::MULTI); // TODO: re-enable the multimapper mode

    GTFToFasta gtfToFasta(gtf_fname, genome_fname,out_fname, kmer_length,multi);
    gtfToFasta.make_transcriptome();
    
    return 0;
}


/* * can think of three different ways to reassign multimappers
   * 1. on the gene/locus level
        * compute average coverage across the entire locus
        * coverage varies when multiple isoforms and may incorrectly reflect the coverage at the region of interest
   * 2. on the transcript level
        * compute average number of reads for each transcript individually
        * problem with overlapping transcripts - how can we handle them
        * 1. build an index to lookup how many isoforms cover a given range of positions
            * count the coverage as the number of reads divided by the number of overlaps for each read
            * alternatively, we could potentially precompute this index for all transcripts and use it to normalize the read counts for that transcript
                * in which case we do not need to perform expensive lookus during runtime and the index stays small
        * 2. preliminary transcriptome quantification
   * 3. on the window level
        * select a window around a multimapper and compute the average coverage
   * */


// Important!!!:: do not need a two-pass algorithm. We can achieve the same easier in a single pass
// Similar to the PID algorithm which for each read will change the target (coverage) and at the next multimapper would compensate given an error