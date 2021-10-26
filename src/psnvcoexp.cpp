
/*
Copyright (c) 2018-2020 Gavin W. Wilson

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
#include "psnvcoexp.hpp"

using namespace sctx;

argagg::parser ProgCoExp::parser() const {
    argagg::parser argparser {{
        { "txidx", {"-i", "--index"},
          "Transcript Index", 1},
        { "barcodes", {"-b", "--barcodes"},
          "Barcode count prefix", 1},
        { "output", {"-o", "--output"},
          "Output prefix", 1},
        { "snvs", {"-s", "--snvs"},
          "Tab separated list of strand specific SNVs must have the following columns chrom, pos, ref, alt, strand", 1},
        { "cellranger", {"-c", "--cellranger"},
          "Indicates the merged bam file is from cell ranger", 0},
        { "tags", {"-t", "--tags"},
          "If this bam file is collapsed write tags that can be used to look for poorly mapping barcode UMI gene combinations", 0},
        { "library", {"-l", "--library"},
          "libary type (V2)", 1},
        { "help", {"-h", "--help"},
          "shows this help message", 1},
      }};
    return argparser;
}


void ProgCoExp::load() {
    /*
    iprefix_ = args_["txidx"].as<std::string>();
    isnvs_ = args_["snvs"].as<std::string>();
    bcin_ = args_["barcodes"].as<std::string>();
    cellranger_ = args_["cellranger"];
    //tags_ = args_["tags"];
    lib_ = args_["library"].as<std::string>("V2");
    outp_ = args_["output"].as<std::string>();
    if(args_.pos.size() != 1){
        throw std::runtime_error("Missing the output option");
    }
    bamin_ = args_.as<std::string>(0);
    */
}


void ProgCoExp::parse_snvs_(){
    /*
    std::unordered_map<std::string, int> tidmap;
    for(auto & r : idx_.refs()){
        tidmap[r.name] = r.tid;
    }
    snvs_.resize(tidmap.size());

    map_idx_ = 0;

    FileWrapper fb(isnvs_);
    {
        ParserTokens toks;
        std::vector<unsigned int> tmp;
        unsigned int chromi, posi, refi, alti, strandi;
        chromi = posi = refi = alti = strandi = std::numeric_limits<unsigned int>::max();
        fb.tokenize_line(toks);
        for(size_t i = 0; i < toks.size(); i++){
            if(toks[i] == "chrom")       chromi = i;
            else if(toks[i] == "pos")    posi = i;
            else if(toks[i] == "ref")    refi = i;
            else if(toks[i] == "alt")    alti = i;
            else if(toks[i] == "strand") strandi = i;
        }

        bool failed = false;
        if(chromi == std::numeric_limits<unsigned int>::max()){
            failed = true;
            std::cout << "Could not find a chrom column\n";
        }
        if(posi == std::numeric_limits<unsigned int>::max()){
            failed = true;
            std::cout << "Could not find a pos column\n";
        }
        if(refi == std::numeric_limits<unsigned int>::max()){
            failed = true;
            std::cout << "Could not find a ref column\n";
        }
        if(alti == std::numeric_limits<unsigned int>::max()){
            failed = true;
            std::cout << "Could not find a alt column\n";
        }
        if(strandi == std::numeric_limits<unsigned int>::max()){
            failed = true;
            std::cout << "Could not find a strand column\n";
        }
        if(failed) exit(EXIT_FAILURE);
        while(fb.tokenize_line(toks) >= 0){
            int tid = tidmap[toks[chromi]];
            unsigned int pos = std::stoi(toks[posi]);
            char ref = toks[refi][0];
            char alt = toks[alti][0];
            char strand = toks[strandi][0];
            if((size_t)tid >= snvs_.size()){
                snvs_.resize(tid + 1);
            }
            snvs_[tid].push_back(SNV(tid, pos, ref, alt, strand, snv_count_));
            snv_count_++;
        }
        size_t tid = 0;
        for(auto & snvs : snvs_){
            if(snvs.empty()) continue;
            //std::cout << " " << tid << " " << idx_.ref(tid).name << " snvs = " << snvs.size() << "\n";
            std::sort(snvs.begin(), snvs.end());
            tid++;
        }
    }
    map_idx_ = snv_count_;
    */
}
