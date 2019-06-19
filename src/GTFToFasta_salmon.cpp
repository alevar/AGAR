#include <utility>

//
//  gtfToFasta.cpp
//  TopHat
//
//  Created by Harold Pimentel on 10/26/11.
//
#include "GTFToFasta_salmon.h"

std::string GTFToFasta::get_exonic_sequence(GffObj &p_trans,FastaRecord &rec, std::string& coords){
    GList<GffExon>& exon_list = p_trans.exons;

    std::string exon_seq;
    size_t length;
    coords.clear();
    std::stringstream ss;

    for (int i = 0; i < exon_list.Count(); ++i) {
        GffExon& cur_exon = *(exon_list.Get(i));
        length = (cur_exon.end+1) - cur_exon.start;

        exon_seq += rec.seq_.substr(cur_exon.start - 1, length);
        ss << ',' << cur_exon.start << '_' << cur_exon.end;
    }

    // get coordinates into the map
    if (length>this->kmerlen and this->multi){ // sanity check for 0 and 1 base exons
        this->found_gene = this->geneMap.find(std::string(p_trans.getGeneID()));
        if(this->found_gene != this->geneMap.end()){
            this->mmap.add_sequence(exon_seq,p_trans,std::get<0>(this->found_gene->second));
        }
        else{
            std::cerr<<"something went wrong with gene ID assignment"<<std::endl;
            std::cerr<<"looking up: "<<p_trans.getGeneID()<<"\t"<<std::string(p_trans.getGeneID())<<std::endl;
            exit(1);
        }
    }

    coords = ss.str().substr(1);
    return exon_seq;
}

GTFToFasta::GTFToFasta(std::string gtf_fname, std::string genome_fname,const std::string& out_fname, int kmerlen,bool multi){
    gtf_fname_ = gtf_fname;
    gtf_fhandle_ = fopen(gtf_fname_.c_str(), "r");
    if (gtf_fhandle_ == nullptr)
    {
        std::cerr << "FATAL: Couldn't open annotation: " << gtf_fname_
        << std::endl;
        exit(1);
    }
    std::cout << "Reading the annotation file: " << gtf_fname_ << std::endl;
    gtfReader_.init(gtf_fhandle_, true); // load recognizable transcript features only
    gtfReader_.readAll();
    std::cout << "loaded the annotation"<<std::endl;

    std::string tlst_fname(out_fname);
    tlst_fname.append(".tlst");
    this->tlst = new std::ofstream(tlst_fname.c_str());

    std::string gene_fname(out_fname);
    gene_fname.append(".glst");
    this->genefp = new std::ofstream(gene_fname.c_str());

    std::string info_fname(out_fname);
    info_fname.append(".info");
    this->infofp = new std::ofstream(info_fname.c_str());

    this->trans_fastaname = out_fname;
    this->trans_fastaname.append(".fasta");

    genome_fname_ = std::move(genome_fname);

    this->multi =  multi;
    this->kmerlen = kmerlen;
    this->mmap.set_kmerlen(kmerlen);
    this->out_fname = out_fname;

    // Make a map from the GffObj
    transcript_map();
}

GTFToFasta::~GTFToFasta(){
    ContigTransMap::iterator it;
    for (it = contigTransMap_.begin(); it != contigTransMap_.end(); ++it) {
        delete it->second;
    }
    this->infofp->close();
    this->tlst->close();
    this->genefp->close();
}

void GTFToFasta::add_to_geneMap(GffObj &p_trans){
    // populate the genemap
    int nst=p_trans.start;
    int nen=p_trans.end;
    char* geneID=p_trans.getGeneID();
    if(geneID==nullptr){
        geneID=p_trans.getGeneName();
        if(geneID==nullptr){
            std::cout<<"wrong geneID"<<std::endl;
            exit(1);
        }
    }
    exists_cur_gene=this->geneMap.insert(std::make_pair(std::string(geneID),std::make_tuple(this->curGeneID,p_trans.strand,nst,nen)));
    if(!exists_cur_gene.second){ // the key did exist - update start and min
        int &st=std::get<2>(exists_cur_gene.first->second);
        int &en=std::get<3>(exists_cur_gene.first->second);
        if(nen>en){ //update end
            en=nen;
        }
        if(nst<st){ //update start
            st=nst;
        }
    }
    else{ // saw a new gene, so need to create a new identifier
        this->curGeneID++;
    }
}

void GTFToFasta::make_transcriptome(){
    std::vector<int> *p_contig_vec;

    FastaReader fastaReader(genome_fname_);
    FastaWriter fastaWriter(this->trans_fastaname);
    FastaRecord cur_contig;

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

            // add the gene record to the geneMap first
            this->add_to_geneMap(*p_trans);

            std::string coordstr;
            out_rec.seq_ = get_exonic_sequence(*p_trans, cur_contig, coordstr);
            if (out_rec.seq_.empty()) continue;
            std::stringstream ss;
            ss << trans_idx;
            out_rec.id_ = ss.str();
            if(std::stoi(out_rec.id_)>this->topTransID){ // check if current id is greater than the highest previously observed
                topTransID = std::stoi(out_rec.id_);
            }
            out_rec.desc_="";
//            out_rec.desc_.append(p_trans->getID());
//            out_rec.desc_.push_back(':');
            this->found_gene = this->geneMap.find(p_trans->getGeneID());
            if(this->found_gene == this->geneMap.end()) { // gene not found
                std::cout << "an error in GeneID ocurred" << std::endl;
                exit(-1);
            }
            out_rec.desc_.append(std::to_string(std::get<0>(this->found_gene->second)));
            out_rec.desc_.push_back('@');
            out_rec.desc_.append(std::to_string(p_trans->gseq_id));
            out_rec.desc_.push_back(p_trans->strand);
            out_rec.desc_.append(coordstr); //list of exon coordinates
            *this->tlst << out_rec.id_ << '\t' << out_rec.desc_ << std::endl;
            fastaWriter.write(out_rec);
        }
    }

    // now time to write the unique kmers for transcripts
    std::unordered_map<std::string,int> uniq_cnt; // counts of unique kmers per transcript
    std::pair<std::unordered_map<std::string,int>::iterator,bool> ex_ucnt; // exists or not

    // write genes to file
    int max_locid = 0; // find what the maximum locus ID is in the current data to be saved into the info file
    auto it=this->geneMap.begin();
    while(it!=this->geneMap.end()){
        max_locid = (std::get<0>(it->second)>max_locid) ? std::get<0>(it->second) : max_locid;
        *this->genefp<<std::get<0>(it->second)<<"@"<<char(std::get<1>(it->second))<<":"<<std::get<2>(it->second)<<" "<<std::get<3>(it->second)<<std::endl;
        it++;
    }

    // write the information about the index now
    *this->infofp<<this->topTransID<<std::endl;
    *this->infofp<<max_locid<<std::endl;
}

void GTFToFasta::transcript_map(){
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

void GTFToFasta::print_mmap(){
    this->mmap.print();
}

void GTFToFasta::save(bool multi, bool uniq){
    if(multi){
        std::string multimap_fname(this->out_fname);
        multimap_fname.append(".multi");
        this->mmap.save_multimappers(multimap_fname);
    }
    if(uniq){
        std::string uniq_fname(this->out_fname);
        uniq_fname.append(".unq");
        this->mmap.save_unique(uniq_fname);
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
    bool multi=args.get_flag(Opt::MULTI);

    GTFToFasta gtfToFasta(gtf_fname, genome_fname,out_fname, kmer_length,multi);
    gtfToFasta.make_transcriptome();
    gtfToFasta.save(Opt::MULTI,Opt::UNIQ);
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