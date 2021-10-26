#pragma once
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

#include "pbase.hpp"
#include <exception>
#include "sparsepp/sparsepp/spp.h"
namespace sctx{

class ProgDebug : public gwsc::ProgBase{
    public:
        using BarcodeHash = spp::sparse_hash_map<std::string, unsigned int>;

        argagg::parser parser() const;
        std::string usage() const {
            return "sctx snvcoexp -i index_prefix -s snvs.tsv -b barcode_counts.txt.gz -o out_prefix bam_in";
        }

        int run();
        void load();

    private:
        template <typename T, typename P>
        int run_wrap_();

        void read_passed_();

        unsigned int filter_func(const std::string & barcode, unsigned int fno){
            (void)fno;
            auto it = bchash_.find(barcode);
            if(it != bchash_.end()){
                return it->second;
            }
            return std::numeric_limits<unsigned int>::max();
        }



        std::vector<std::string>      barcodes_;
        BarcodeHash                   bchash_;
        std::string iprefix_;
        std::string passed_;
        std::string lib_;
        std::string bamin_;
        std::string lib_type_;

};

}
