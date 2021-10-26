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
#include "pdebug.hpp"
#include "reader.hpp"
#include "bam_genes.hpp"

using namespace sctx;
using namespace gwsc;

argagg::parser ProgDebug::parser() const {
    argagg::parser argparser {{
        { "txidx", {"-i", "--index"},
          "Transcript Index", 1},
        { "library", {"-l", "--library"},
          "libary type (V2)", 1},
        { "passed", {"-p", "--passed"},
         "Only process reads from the barcodes in this file", 1},
        { "help", {"-h", "--help"},

          "shows this help message", 1},
      }};
    return argparser;
}


void ProgDebug::load() {
    iprefix_ = args_["txidx"].as<std::string>();
    passed_ = args_["passed"].as<std::string>();
    bamin_ = args_.as<std::string>(0);
    lib_type_ = args_["library"].as<std::string>("V2");
}

void ProgDebug::read_passed_(){
    FileWrapper in(passed_);
    std::string line;
    in.get_line(line);
    size_t index = 0;
    while(in.get_line(line) > -1){
        barcodes_.push_back(line);
        bchash_[line] = index++;
    }
    tout << "Read " << barcodes_.size() << " passed barcodes\n";
}

template <typename T, typename P>
int ProgDebug::run_wrap_(){
    read_passed_();
    tout << "Loading the transcriptome index\n";

    typename BamGeneReaderFiltered<T, BamReader, P>::bcfilter_func fp = std::bind(&ProgDebug::filter_func, *this, std::placeholders::_1, std::placeholders::_2);
    BamGeneReaderFiltered<T, BamReader, P> br(fp); //br(cellranger_);

    br.index.load(iprefix_);
    br.index.build_splice_site_index();
    br.set_bam(bamin_);


    BamBuffer * rbuffer = new BamBuffer();
    unsigned int tot = 0;
    while((tot = br.read_genes(*rbuffer, 500)) > 0){
        std::cout << "Read " << tot << "\n";
        BamBuffer::rpair range;
        unsigned int prgt = 0;
        int ptid = -1;
        while(rbuffer->get_next(range)){
            unsigned int rgt = 0, lft = std::numeric_limits<unsigned int>::max();
            int tid = -1;
            for(auto it = range.first; it != range.second; it++){
                BamDetail &d = **it;
                if(d.b->core.tid != ptid){
                    prgt = 0;
                    ptid = -1;
                }
                if(tid == -1) tid = d.b->core.tid;
                else if(tid != d.b->core.tid){
                    std::cout << "Two different tid in the same range tid = " << tid << " other = " << d.b->core.tid << "\n";
                }
                lft = std::min((unsigned int)d.b->core.pos, lft);
                rgt = std::max((unsigned int)bam_endpos(d.b), rgt);
            }
            if(lft < prgt){
                std::cout << "Overlapping range " << lft << " - " << rgt << " overlaps previous rgt of " << prgt << "\n";
            }
            prgt = rgt;
        }
    }

    return EXIT_SUCCESS;
}

int ProgDebug::run() {
    if(lib_type_ == "V2"){
        return run_wrap_<Reader10X_V2, BamScSNVProcessor>();
    }else if(lib_type_ == "V3"){
        return run_wrap_<Reader10X_V3, BamScSNVProcessor>();
    }else if(lib_type_ == "V2_5P"){
        return run_wrap_<Reader10X_V2_5P, BamScSNVProcessor>();
    }else if(lib_type_ == "V3_5P"){
        return run_wrap_<Reader10X_V3_5P, BamScSNVProcessor>();
    }
    return EXIT_SUCCESS;
}
